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

#ifndef __QUADRATURE_QGAUSSLINE_H_
#define __QUADRATURE_QGAUSSLINE_H_

#include "quadrature/quadraturetype.h"

namespace hiflow {

///
/// \class QuadratureGaussLine qgaussline.h
/// \brief Derived class of all gaussian quadratures on lines
///        which is templated by a DataType (e.g. double)
/// \author Michael Schick
///

template < class DataType >
class QuadratureGaussLine : public QuadratureType< DataType > {
public:
  /// Constructor which sets up all necessary information about the
  /// Gaussian quadrature on a line for all implemented sizes
  QuadratureGaussLine();

  virtual QuadratureGaussLine< DataType > *clone() const;

  using QuadratureType< DataType >::x_;
  using QuadratureType< DataType >::y_;
  using QuadratureType< DataType >::z_;
  using QuadratureType< DataType >::w_;
  using QuadratureType< DataType >::index_field_;
  using QuadratureType< DataType >::order_size_map_;
};

} // namespace hiflow

#endif
