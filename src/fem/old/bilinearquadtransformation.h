// Copyright (C) 2011-2021 Vincent Heuveline
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

#ifndef __FEM_BILINEAR_QUAD_TRANSFORMATION_H_
#define __FEM_BILINEAR_QUAD_TRANSFORMATION_H_

#include "fem/cell_transformation.h"

namespace hiflow {
namespace doffem {

///
/// \class BiLinearQuadTransformation bilinearquadtransformation.h
/// \brief Bilinear transformation mapping from reference to physical cell for a
/// Quadrilateral \author Michael Schick<br>Martin Baumann
///

template < class DataType >
class BiLinearQuadTransformation : public CellTransformation< DataType > {
public:
  typedef std::vector< DataType > Coord;

  explicit BiLinearQuadTransformation(int gdim);

  std::vector< std::vector< DataType > > get_reference_coordinates() const;

  bool inverse(DataType x_phy, DataType &x_ref) const;
  bool inverse(DataType x_phy, DataType y_phy, DataType &x_ref,
               DataType &y_ref) const;
  bool inverse(DataType x_phy, DataType y_phy, DataType z_phy, DataType &x_ref,
               DataType &y_ref, DataType &z_ref) const;

  DataType x(const Coord &coord_ref) const;
  DataType x_x(const Coord &coord_ref) const;
  DataType x_y(const Coord &coord_ref) const;
  DataType x_xx(const Coord &coord_ref) const;
  DataType x_xy(const Coord &coord_ref) const;
  DataType x_yy(const Coord &coord_ref) const;

  DataType y(const Coord &coord_ref) const;
  DataType y_x(const Coord &coord_ref) const;
  DataType y_y(const Coord &coord_ref) const;
  DataType y_xx(const Coord &coord_ref) const;
  DataType y_xy(const Coord &coord_ref) const;
  DataType y_yy(const Coord &coord_ref) const;

  bool contains_reference_point(const Coord &coord_ref) const;

protected:
  bool inverse_by_decomposition(DataType x_phy, DataType y_phy, DataType &x_ref,
                                DataType &y_ref) const;

  int8_t triangle_decomp_ind_[2][3];
};

} // namespace doffem
} // namespace hiflow

#endif
