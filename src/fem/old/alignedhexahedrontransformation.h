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

#ifndef __FEM_ALIGNED_HEXAHEDRON_TRANSFORMATION_H_
#define __FEM_ALIGNED_HEXAHEDRON_TRANSFORMATION_H_

#include "fem/cell_transformation.h"
#include "fem/trilinearhexahedrontransformation.h"
#include <cassert>

#define SAVE_MODE

namespace hiflow {
namespace doffem {

///
/// \class AlignedHexahedronTransformation alignedhexahedrontransformation.h
/// \brief Trilinear transformation mapping from reference to physical cell for
/// a Hexahedron which is assumed to be axis aligned \author Philipp Gerstner
///

template < class DataType >
class AlignedHexahedronTransformation : public CellTransformation< DataType > {
public:
  typedef std::vector< DataType > Coord;

  explicit AlignedHexahedronTransformation(int gdim);

  void reinit(const Coord &coord_vtx);

  std::vector< std::vector< DataType > > get_reference_coordinates() const;

  bool inverse(DataType x_phy, DataType &x_ref) const;
  bool inverse(DataType x_phy, DataType y_phy, DataType &x_ref,
               DataType &y_ref) const;
  bool inverse(DataType x_phy, DataType y_phy, DataType z_phy, DataType &x_ref,
               DataType &y_ref, DataType &z_ref) const;

  DataType x(const Coord &coord_ref) const;
  DataType x_x(const Coord &coord_ref) const;
  DataType x_y(const Coord &coord_ref) const;
  DataType x_z(const Coord &coord_ref) const;
  DataType x_xx(const Coord &coord_ref) const;
  DataType x_xy(const Coord &coord_ref) const;
  DataType x_xz(const Coord &coord_ref) const;
  DataType x_yy(const Coord &coord_ref) const;
  DataType x_yz(const Coord &coord_ref) const;
  DataType x_zz(const Coord &coord_ref) const;

  DataType y(const Coord &coord_ref) const;
  DataType y_x(const Coord &coord_ref) const;
  DataType y_y(const Coord &coord_ref) const;
  DataType y_z(const Coord &coord_ref) const;
  DataType y_xx(const Coord &coord_ref) const;
  DataType y_xy(const Coord &coord_ref) const;
  DataType y_xz(const Coord &coord_ref) const;
  DataType y_yy(const Coord &coord_ref) const;
  DataType y_yz(const Coord &coord_ref) const;
  DataType y_zz(const Coord &coord_ref) const;

  DataType z(const Coord &coord_ref) const;
  DataType z_x(const Coord &coord_ref) const;
  DataType z_y(const Coord &coord_ref) const;
  DataType z_z(const Coord &coord_ref) const;
  DataType z_xx(const Coord &coord_ref) const;
  DataType z_xy(const Coord &coord_ref) const;
  DataType z_xz(const Coord &coord_ref) const;
  DataType z_yy(const Coord &coord_ref) const;
  DataType z_yz(const Coord &coord_ref) const;
  DataType z_zz(const Coord &coord_ref) const;

  bool contains_reference_point(const Coord &coord_ref) const;

protected:
  TriLinearHexahedronTransformation< DataType > *T3_;

  DataType dx_;
  DataType dy_;
  DataType dz_;

  int origin_;
  DataType x0_;
  DataType y0_;
  DataType z0_;
};

} // namespace doffem
} // namespace hiflow

#endif
