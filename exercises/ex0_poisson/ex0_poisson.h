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

// System includes.
#include <mpi.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "hiflow.h"

// TODO exercise: uncomment the following line for preparing sub exercise d)
// #define SUBEX_D

// All names are imported for simplicity.
using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

// Shorten some datatypes with typedefs.
typedef LADescriptorCoupledD LAD;
typedef LAD::DataType DataType;
typedef LAD::VectorType VectorType;
typedef LAD::MatrixType MatrixType;

// DIM of the problem.
const int DIM = 2;

typedef Vec<DIM, DataType> Coord;

// Rank of the master process.
const int MASTER_RANK = 0;

// Parameters M, N and O in solution. These decide the period in x-, y- and
// z-direction respectively.
const int M = 1;
const int N = 1;
const int O = 1;

// Functor used to impose u(x) = c on the boundary.
struct DirichletConstant {
  DirichletConstant(DataType c) : c_(c) {}

  void evaluate(const mesh::Entity &face,
                const Vec<DIM, DataType> &coords_on_face,
                std::vector<DataType> &vals) const {
    // return array with Dirichlet values for dof:s on boundary face
    vals = std::vector<DataType>(1, c_);
  }

  size_t nb_comp() const { return 1; }

  size_t nb_func() const { return 1; }

  size_t iv2ind(size_t j, size_t v) const { return v; }

  DataType c_;
};

// Functor used for the local assembly of the stiffness matrix and load vector.
template <class ExactSol>
class LocalPoissonAssembler : private AssemblyAssistant<DIM, DataType> {
 public:
  // compute local matrix
  // [in]  element:    contains information about current cell
  // [in]  quadrature: quadrature rule to be used for approximating the
  // integrals [out] lm: contribution of the current cell to the global system
  // matrix

  void operator()(const Element<DataType, DIM> &element,
                  const Quadrature<DataType> &quadrature, LocalMatrix &lm) {
    const bool need_basis_hessians = false;

    // AssemblyAssistant sets up the local FE basis functions for the current // cell
    AssemblyAssistant<DIM, DataType>::initialize_for_element(element, quadrature, need_basis_hessians);
    const size_t num_dof = this->num_dofs_total();  // number of degrees of freedom on current cell
    const int num_q = this->num_quadrature_points();  // number of quadrature points

    for (int q = 0; q < num_q; ++q) {  // loop over quadrature points quadrature weight
      const DataType wq = w(q);
      const DataType dJ = std::abs(this->detJ(q)); // volume element of cell transformation

      for (int i = 0; i < num_dof; ++i) { // loop over test DOFs <-> test function v
        for (int j = 0; j < num_dof; ++j) { // loop over trrial DOFs <-> trial function u
          lm(i, j) += wq * dot(this->grad_Phi(j, q, 0), this->grad_Phi(i, q, 0)) * dJ;
        }
      }

    }
  }

  // compute local right hand side vector
  // [in]  element:    contains information about current cell
  // [in]  quadrature: quadrature rule to be used for approximating the
  // integrals [out] lv: contribution of the current cell to the global system
  // right hand side
  void operator()(const Element<DataType, DIM> &element,
                  const Quadrature<DataType> &quadrature, LocalVector &lv) {
    const bool need_basis_hessians = false;
    // AssemblyAssistant sets up the local FE basis functions for the current cell
    AssemblyAssistant<DIM, DataType>::initialize_for_element(element, quadrature, need_basis_hessians);
    const size_t num_dof = this->num_dofs_total(); // number of degrees of freedom on current cell
    const int num_q = this->num_quadrature_points(); // number of quadrature points
    
    for (int q = 0; q < num_q; ++q){ // loop over quadrature points
      const DataType wq = w(q); // quadrature weight
      const DataType dJ = std::abs(this->detJ(q)); // volume element of cell transformation
      
      for (int i = 0; i < num_dof; ++i) { // loop over test DOFs <-> test function v
        // x(q): physical coordinates of quadrature point q
        lv[i] += wq * f(this->x(q)) * this->Phi(i, q, 0) * dJ;
      }

    }
  }

  DataType f(Vec<DIM, DataType> pt) {
    ExactSol sol;
    DataType rhs_sol;
    rhs_sol = 4. * M_PI * M_PI * (M * M + N * N) * sol(pt);
    return rhs_sol;
  }
};


struct ExactSol {
  size_t nb_comp() const { return 1; }

  size_t nb_func() const { return 1; }

  size_t weight_size() const { return nb_func() * nb_comp(); }

  size_t iv2ind(size_t j, size_t v) const { return v; }

  // wrapper needed to make ExactSol compatible with FE interpolation
  void evaluate(const mesh::Entity &cell, const Vec<DIM, DataType> &pt,
                std::vector<DataType> &vals) const {
    vals.clear();
    vals.resize(1, 0.);
    vals[0] = this->operator()(pt);
  }

  // evaluate analytical function at given point
  DataType operator()(const Vec<DIM, DataType> &pt) const {
    const DataType x = pt[0];
    const DataType y = (DIM > 1) ? pt[1] : 0;
    const DataType z = (DIM > 2) ? pt[2] : 0;
    const DataType pi = M_PI;
    DataType solution;

    switch (DIM) {
      case 2: {
        solution = 10.0 * std::sin(2. * M * pi * x) * std::sin(2. * N * pi * y);
        break;
      }
      case 3: {
        solution = 10.0 * std::sin(2. * M * pi * x) *
                   std::sin(2. * N * pi * y) * std::sin(2. * O * pi * z);
        break;
      }

      default:
        assert(0);
    }
    return solution;
  }

  // evaluate gradient of analytical function at given point
  Vec<DIM, DataType> eval_grad(const Vec<DIM, DataType> &pt) const {
    Vec<DIM, DataType> grad;
    const DataType x = pt[0];
    const DataType y = (DIM > 1) ? pt[1] : 0;
    const DataType z = (DIM > 2) ? pt[2] : 0;
    const DataType pi = M_PI;

    switch (DIM) {
      case 2: {
        grad[0] = 20. * M * pi * std::cos(2. * M * pi * x) *
                  std::sin(2. * N * pi * y);
        grad[1] = 20. * N * pi * std::sin(2. * M * pi * x) *
                  std::cos(2. * N * pi * y);
        break;
      }
      case 3: {
        grad[0] = 20. * M * pi * std::cos(2. * M * pi * x) *
                  std::sin(2. * N * pi * y) * std::sin(2. * O * pi * z);
        grad[1] = 20. * N * pi * std::sin(2. * M * pi * x) *
                  std::cos(2. * N * pi * y) * std::sin(2. * O * pi * z);
        grad[2] = 20. * O * pi * std::sin(2. * M * pi * x) *
                  std::sin(2. * N * pi * y) * std::cos(2. * O * pi * z);
        break;
      }
      default:
        assert(0);
    }
    return grad;
  }
};
