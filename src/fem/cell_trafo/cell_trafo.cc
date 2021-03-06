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

#include "fem/cell_trafo/cell_transformation.h"
#include "fem/cell_trafo/linearlinetransformation.h"
#include "fem/cell_trafo/lineartriangletransformation.h"
#include "fem/cell_trafo/lineartetrahedrontransformation.h"
#include "fem/cell_trafo/linearpyramidtransformation.h"
#include "fem/cell_trafo/bilinearquadtransformation.h"
#include "fem/cell_trafo/trilinearhexahedrontransformation.h"
#include "fem/cell_trafo/alignedquadtransformation.h"
#include "fem/cell_trafo/alignedhexahedrontransformation.h"

namespace hiflow {
namespace doffem {

template class CellTransformation <double, 3>;
template class LinearLineTransformation< double, 1 >;
template class LinearTriangleTransformation< double, 2 >;
template class LinearTetrahedronTransformation< double, 3 >;
template class LinearPyramidTransformation< double, 3 >;
template class BiLinearQuadTransformation< double, 2 >;
template class TriLinearHexahedronTransformation< double, 3 >;
template class AlignedQuadTransformation< double, 2 >;
template class AlignedHexahedronTransformation< double, 3 >;

} // namespace doffem
} // namespace hiflow
