// Functor to evaluate the exact solution u of the Poisson problem
// with Dirichlet BC, and its gradient \grad u.


// struct ExactSol {
//   size_t nb_comp() const { return 1; }

//   size_t nb_func() const { return 1; }

//   size_t weight_size() const { return nb_func() * nb_comp(); }

//   size_t iv2ind(size_t j, size_t v) const { return v; }

//   // wrapper needed to make ExactSol compatible with FE interpolation
//   void evaluate(const mesh::Entity &cell, const Vec<DIM, DataType> &pt,
//                 std::vector<DataType> &vals) const {
//     vals.clear();
//     vals.resize(1, 0.);
//     vals[0] = this->operator()(pt);
//   }

//   // evaluate analytical function at given point
//   DataType operator()(const Vec<DIM, DataType> &pt) const {
//     const DataType x = pt[0];
//     const DataType y = (DIM > 1) ? pt[1] : 0;
//     const DataType z = (DIM > 2) ? pt[2] : 0;
//     const DataType pi = M_PI;
//     DataType solution;

//     switch (DIM) {
//       case 2: {
//         solution = 10.0 * std::sin(2. * M * pi * x) * std::sin(2. * N * pi * y);
//         break;
//       }
//       case 3: {
//         solution = 10.0 * std::sin(2. * M * pi * x) *
//                    std::sin(2. * N * pi * y) * std::sin(2. * O * pi * z);
//         break;
//       }

//       default:
//         assert(0);
//     }
//     return solution;
//   }

//   // evaluate gradient of analytical function at given point
//   Vec<DIM, DataType> eval_grad(const Vec<DIM, DataType> &pt) const {
//     Vec<DIM, DataType> grad;
//     const DataType x = pt[0];
//     const DataType y = (DIM > 1) ? pt[1] : 0;
//     const DataType z = (DIM > 2) ? pt[2] : 0;
//     const DataType pi = M_PI;

//     switch (DIM) {
//       case 2: {
//         grad[0] = 20. * M * pi * std::cos(2. * M * pi * x) *
//                   std::sin(2. * N * pi * y);
//         grad[1] = 20. * N * pi * std::sin(2. * M * pi * x) *
//                   std::cos(2. * N * pi * y);
//         break;
//       }
//       case 3: {
//         grad[0] = 20. * M * pi * std::cos(2. * M * pi * x) *
//                   std::sin(2. * N * pi * y) * std::sin(2. * O * pi * z);
//         grad[1] = 20. * N * pi * std::sin(2. * M * pi * x) *
//                   std::cos(2. * N * pi * y) * std::sin(2. * O * pi * z);
//         grad[2] = 20. * O * pi * std::sin(2. * M * pi * x) *
//                   std::sin(2. * N * pi * y) * std::cos(2. * O * pi * z);
//         break;
//       }
//       default:
//         assert(0);
//     }
//     return grad;
//   }
// };






// Functor used for the local evaluation of the square of the L2-norm of the
// error on each element.

// template <class ExactSol>
// class L2ErrorIntegrator : private AssemblyAssistant<DIM, DataType> {
//  public:
//   L2ErrorIntegrator(VectorType &pp_sol) : pp_sol_(pp_sol) {}

//   void operator()(const Element<DataType, DIM> &element,
//                   const Quadrature<DataType> &quadrature, DataType &value) {
//     const bool need_basis_hessians = false;

//     // AssemblyAssistant sets up the local FE basis functions for the current
//     // cell
//     AssemblyAssistant<DIM, DataType>::initialize_for_element(
//         element, quadrature, need_basis_hessians);

//     // Evaluate the computed solution at all quadrature points.
//     evaluate_fe_function(pp_sol_, 0, approx_sol_);

//     const int num_q = num_quadrature_points();
//     for (int q = 0; q < num_q; ++q) {
//       const DataType wq = w(q);
//       const DataType dJ = std::abs(detJ(q));
//       const DataType delta = sol_(x(q)) - approx_sol_[q];

//       value += wq * delta * delta * dJ;
//     }
//   }

//  private:
//   // coefficients of the computed solution
//   const VectorType &pp_sol_;

//   // functor to evaluate exact solution
//   ExactSol sol_;

//   // vector with values of computed solution evaluated at each quadrature point
//   FunctionValues<DataType> approx_sol_;
// };

// // Functor used for the local evaluation of the square of the H1-norm of the
// // error on each element.

// template <class ExactSol>
// class H1ErrorIntegrator : private AssemblyAssistant<DIM, DataType> {
//  public:
//   H1ErrorIntegrator(VectorType &pp_sol) : pp_sol_(pp_sol) {}

//   void operator()(const Element<DataType, DIM> &element,
//                   const Quadrature<DataType> &quadrature, DataType &value) {
//     const bool need_basis_hessians = false;

//     // AssemblyAssistant sets up the local FE basis functions for the current
//     // cell
//     AssemblyAssistant<DIM, DataType>::initialize_for_element(
//         element, quadrature, need_basis_hessians);

//     // Evaluate the gradient of the computed solution at all quadrature points.
//     evaluate_fe_function_gradients(pp_sol_, 0, approx_grad_u_);

//     const int num_q = num_quadrature_points();
//     for (int q = 0; q < num_q; ++q) {
//       const DataType wq = w(q);
//       const DataType dJ = std::abs(detJ(q));
//       const Vec<DIM, DataType> grad_u = sol_.eval_grad(x(q));

//       value += wq * dot(grad_u - approx_grad_u_[q], grad_u - approx_grad_u_[q]) * dJ;
//     }
//   }

//  private:
//   // coefficients of the computed solution
//   const VectorType &pp_sol_;

//   // functor to evaluate exact solution
//   ExactSol sol_;

//   // gradient of computed solution evaluated at each quadrature point
//   FunctionValues<Vec<DIM, DataType> > approx_grad_u_;
// };
