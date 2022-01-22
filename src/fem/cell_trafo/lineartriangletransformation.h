// Copyright (C) 2011-2017 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for
// more details.
//
// You should have received a copy of the European Union Public Licence (EUPL)
// v1.2 along with HiFlow3.  If not, see
// <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

#ifndef __FEM_LINEAR_TRIANGLE_TRANSFORMATION_H_
#define __FEM_LINEAR_TRIANGLE_TRANSFORMATION_H_

#include "fem/cell_trafo/cell_transformation.h"
#include <cmath>
#include <iomanip>
#include <iostream>

namespace hiflow {
namespace doffem {

///
/// \class LinearTriangleTransformation lineartriangletransformation.h
/// \brief Linear transformation mapping from reference to physical cell for a
/// Triangle \author Michael Schick<br>Martin Baumann
///

template < class DataType, int DIM >
class LinearTriangleTransformation : public CellTransformation< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  explicit LinearTriangleTransformation(ConstRefCellPtr<DataType, DIM> ref_cell);

  bool differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const;
   
  bool inverse(Coord co_phy, Coord &co_ref) const;
  bool inverse_2Dto3D(Vec< 3, DataType > co_phy, Coord &co_ref) const;

  void reinit(const std::vector<DataType> &coord_vtx);
  
  DataType x(const Coord &coord_ref) const;
  DataType x_x(const Coord &coord_ref) const;
  DataType x_y(const Coord &coord_ref) const;

  DataType y(const Coord &coord_ref) const;
  DataType y_x(const Coord &coord_ref) const;
  DataType y_y(const Coord &coord_ref) const;
  
  DataType z(const Coord &coord_ref) const;
  DataType z_x(const Coord &coord_ref) const;
  DataType z_y(const Coord &coord_ref) const;

private:
  Mat<DIM, DIM, DataType> A_;
  Mat<DIM, DIM, DataType> Ainv_;
};

template < class DataType, int DIM >
LinearTriangleTransformation< DataType, DIM >::LinearTriangleTransformation(ConstRefCellPtr<DataType, DIM> ref_cell)
    : CellTransformation< DataType, DIM >(ref_cell) 
{
  this->order_ = 1;
  this->fixed_ref_cell_type_ = REF_CELL_TRI_STD;
  this->name_ = "Triangle";
}

template < class DataType, int DIM >
bool LinearTriangleTransformation< DataType, DIM >::differs_by_translation_from (ConstCellTrafoPtr<DataType, DIM> rhs) const
{
  assert(this->ref_cell_);
  
  if (this->ref_cell_->type() != rhs->get_ref_cell()->type())
  {
    return false;
  }
  
  std::vector< Coord > ref_coords = this->ref_cell_->get_coords();
  assert (ref_coords.size() == 3);
  
  Coord my_p01 = this->transform(ref_coords[0]) - this->transform(ref_coords[1]); 
  Coord rhs_p01 = rhs->transform(ref_coords[0]) - rhs->transform(ref_coords[1]); 
  if (my_p01 != rhs_p01)
  {
    return false;
  }
  Coord my_p02 = this->transform(ref_coords[0]) - this->transform(ref_coords[2]); 
  Coord rhs_p02 = rhs->transform(ref_coords[0]) - rhs->transform(ref_coords[2]); 
  if (my_p02 != rhs_p02)
  {
    return false;
  }
  return true;
}

template < class DataType, int DIM >
void LinearTriangleTransformation< DataType, DIM >::reinit(const std::vector<DataType> &coord_vtx) {
  assert (DIM == 2);
  
  assert (coord_vtx.size() == DIM * 3);
  
  this->coord_vtx_.clear();
  this->coord_vtx_.resize(3);
  
  for (int i=0; i<3; ++i)
  {
    for (int d=0; d<DIM; ++d)
    {
      this->coord_vtx_[i][d] = coord_vtx[this->ij2ind(i,d)];
    }
  }
  
  this->A_(0,0) = this->coord_vtx_[1][0] - this->coord_vtx_[0][0];
  this->A_(1,0) = this->coord_vtx_[1][1] - this->coord_vtx_[0][1];
  this->A_(0,1) = this->coord_vtx_[2][0] - this->coord_vtx_[0][0];
  this->A_(1,1) = this->coord_vtx_[2][1] - this->coord_vtx_[0][1];

#ifndef NDEBUG   
  DataType determ = this->A_(0,0) * this->A_(1,1) - this->A_(0,1) * this->A_(1,0);
  assert (std::abs(determ) > 1e-8);
#endif

  inv(this->A_, this->Ainv_);
}

template < class DataType, int DIM >
bool LinearTriangleTransformation< DataType, DIM >::inverse(Coord co_phy, Coord &co_ref) const {
  assert (DIM == 2);
  Coord rhs = co_phy - this->coord_vtx_[0];
  this->Ainv_.VectorMult(rhs, co_ref);
  
  const DataType eps = this->ref_cell_->eps();
  if (std::abs(co_ref[0]) < eps) 
  {
    co_ref[0] = 0.;
  } 
  else if (std::abs(co_ref[0] - 1.) < eps) 
  {
    co_ref[0] = 1.;
  }
    
  if (std::abs(co_ref[1]) < eps) 
  {
    co_ref[1] = 0.;
  } 
  else if (std::abs(co_ref[1] - 1.) < eps) 
  {
    co_ref[1] = 1.;
  }
  
  /*
  DataType a11 = this->coord_vtx_[1][0] -
                 this->coord_vtx_[0][0];
  DataType a12 = this->coord_vtx_[2][0] -
                 this->coord_vtx_[0][0];
  DataType a21 = this->coord_vtx_[1][1] -
                 this->coord_vtx_[0][1];
  DataType a22 = this->coord_vtx_[2][1] -
                 this->coord_vtx_[0][1];

  DataType det = a11 * a22 - a21 * a12;

  assert(det != 0.0);

  co_ref[0] = (1.0 / det) * (a22 * (co_phy[0] - this->coord_vtx_[0][0]) -
                         a12 * (co_phy[1] - this->coord_vtx_[0][1]));
  co_ref[1] = (1.0 / det) * (-a21 * (co_phy[0] - this->coord_vtx_[0][0]) +
                         a11 * (co_phy[1] - this->coord_vtx_[0][1]));
  */
  return this->contains_reference_point(co_ref);
}

template < class DataType, int DIM >
DataType LinearTriangleTransformation< DataType, DIM >::x(const Coord &coord_ref) const {
  assert (DIM == 2);
  return (this->coord_vtx_[1][0] - this->coord_vtx_[0][0]) * coord_ref[0] 
        +(this->coord_vtx_[2][0] - this->coord_vtx_[0][0]) * coord_ref[1] 
        + this->coord_vtx_[0][0];
}

template < class DataType, int DIM >
DataType LinearTriangleTransformation< DataType, DIM >::x_x(const Coord &coord_ref) const {
  assert (DIM == 2);
  return (this->coord_vtx_[1][0] - this->coord_vtx_[0][0]);
}

template < class DataType, int DIM >
DataType LinearTriangleTransformation< DataType, DIM >::x_y(const Coord &coord_ref) const {
  assert (DIM == 2);
  return (this->coord_vtx_[2][0] - this->coord_vtx_[0][0]);
}

template < class DataType, int DIM >
DataType LinearTriangleTransformation< DataType, DIM >::y(const Coord &coord_ref) const {
  assert (DIM == 2);
  return (this->coord_vtx_[1][1] - this->coord_vtx_[0][1]) * coord_ref[0] 
        +(this->coord_vtx_[2][1] - this->coord_vtx_[0][1]) * coord_ref[1] 
        + this->coord_vtx_[0][1];
}

template < class DataType, int DIM >
DataType LinearTriangleTransformation< DataType, DIM >::y_x(const Coord &coord_ref) const {
  assert (DIM == 2);
  return (this->coord_vtx_[1][1] - this->coord_vtx_[0][1]);
}

template < class DataType, int DIM >
DataType LinearTriangleTransformation< DataType, DIM >::y_y(const Coord &coord_ref) const {
  assert (DIM == 2);
  return (this->coord_vtx_[2][1] - this->coord_vtx_[0][1]);
}

// begin preliminary test functions when the dimension of the physical point is > 2

template < class DataType, int DIM >
DataType LinearTriangleTransformation< DataType, DIM >::z(const Coord &coord_ref) const {
  
  return (this->coord_vtx_[1][2] - this->coord_vtx_[0][2]) * coord_ref[0] 
        +(this->coord_vtx_[2][2] - this->coord_vtx_[0][2]) * coord_ref[1] 
        + this->coord_vtx_[0][2];
}

template < class DataType, int DIM >
DataType LinearTriangleTransformation< DataType, DIM >::z_x(const Coord &coord_ref) const {
  assert (DIM == 2);
  return (this->coord_vtx_[1][2] - this->coord_vtx_[0][2]);
}

template < class DataType, int DIM >
DataType LinearTriangleTransformation< DataType, DIM >::z_y(const Coord &coord_ref) const {
  assert (DIM == 2);
  return (this->coord_vtx_[2][2] - this->coord_vtx_[0][2]);
}

//to invert the matrix, we have to solve a 3 x 2 linear system. Since one row depends linearly on the other two rows,
//we indentify two linearly independent rows and solve the remaining 2x2 system.
template < class DataType, int DIM >
bool LinearTriangleTransformation< DataType, DIM >::inverse_2Dto3D(Vec< 3, DataType > co_phy, Coord &co_ref) const {
  
  DataType a11 = this->coord_vtx_[1][0] -
                 this->coord_vtx_[0][0];
  DataType a12 = this->coord_vtx_[2][0] -
                 this->coord_vtx_[0][0];
  DataType a21 = this->coord_vtx_[1][1] -
                 this->coord_vtx_[0][1];
  DataType a22 = this->coord_vtx_[2][1] -
                 this->coord_vtx_[0][1];
  DataType a31 = this->coord_vtx_[1][2] -
                 this->coord_vtx_[0][2];
  DataType a32 = this->coord_vtx_[2][2] -
                 this->coord_vtx_[0][2];
 

  DataType det1 = a11 * a22 - a21 * a12;
  DataType det2 = a11 * a32 - a31 * a12;
  DataType det3 = a21 * a32 - a31 * a22;
  
  if (det1 != 0.0) 
  {
    co_ref[0] = (1.0 / det1) * (a22 * (co_phy[0] - this->coord_vtx_[0][0]) -
                           a12 * (co_phy[1] - this->coord_vtx_[0][1]));
    co_ref[1] = (1.0 / det1) * (-a21 * (co_phy[0] - this->coord_vtx_[0][0]) +
                           a11 * (co_phy[1] - this->coord_vtx_[0][1]));
  }
  else if (det2 != 0) 
  {
    co_ref[0] = (1.0 / det2) * (a32 * (co_phy[0] - this->coord_vtx_[0][0]) -
                           a12 * (co_phy[2] - this->coord_vtx_[0][2]));
    co_ref[1] = (1.0 / det2) * (-a31 * (co_phy[0] - this->coord_vtx_[0][0]) +
                           a11 * (co_phy[2] - this->coord_vtx_[0][2]));
  }
  else if (det3 != 0) 
  {
    co_ref[0] = (1.0 / det3) * (a32 * (co_phy[1] - this->coord_vtx_[0][1]) -
                           a22 * (co_phy[2] - this->coord_vtx_[0][2]));
    co_ref[1] = (1.0 / det3) * (-a31 * (co_phy[1] - this->coord_vtx_[0][1]) +
                           a21 * (co_phy[2] - this->coord_vtx_[0][2]));
  }
   
  else 
  {
    return false; //triangle is degenerate
  }
  
  return this->contains_reference_point(co_ref);
}





} // namespace doffem
} // namespace hiflow

#endif
