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

#ifndef HIFLOW_LINEARSOLVER_HYPRE_LINEAR_SOLVER_H_
#define HIFLOW_LINEARSOLVER_HYPRE_LINEAR_SOLVER_H_

#include "common/log.h"
#include "common/timer.h"
#include "linear_solver/hypre_preconditioner.h"
#include "linear_solver/linear_solver.h"
#include <cstdlib>

namespace hiflow {
namespace la {
/// @author Simon Gawlok

/// @brief Base class for wrappers to solvers in Hypre library
///
/// A linear solver is in particular a preconditioner.

template < class LAD >
class HypreLinearSolver : public HyprePreconditioner< LAD >,
                          public LinearSolver< LAD > {
public:
  typedef typename LAD::MatrixType OperatorType;
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::DataType DataType;

  HypreLinearSolver() : LinearSolver< LAD >(), hypre_precond_(nullptr) {
    this->SetInitialized(false);
    this->comm_ = MPI_COMM_NULL;
    this->is_critical_hypre_solver_ = false;
    this->SetReuse(false);
  };

  virtual ~HypreLinearSolver(){};

  /// Clear allocated data

  virtual void Clear() {
    if (this->comm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&this->comm_);
    }

    LinearSolver< LAD >::Clear();
    this->SetInitialized(false);
    this->SetModifiedParam(false);
    this->is_critical_hypre_solver_ = false;
  }

  /// Initialize convergence control
  /// @param[in] maxits Maximum number of iterations
  /// @param[in] abstol Absolute tolerance for residual
  /// @param[in] reltol Relative tolerance for residual
  /// @param[in] divtol Divergence tolerance for residual (not used by Hypre
  /// solvers so far)

  virtual void InitControl(int maxits, DataType abstol, DataType reltol,
                           DataType divtol = 1e6) {
    this->maxits_ = maxits;
    this->abstol_ = abstol;
    this->reltol_ = reltol;
    this->divtol_ = divtol;
  }

  /// Sets up the operator, e.g. the system matrix.
  /// If implemented, it must be invoked before
  /// @c ApplyPreconditioner is called.
  /// For instance it could compute an ILU decomposition of a GlobalMatrix.
  /// Sets up the operator for the linear solver.
  /// @param op linear operator

  virtual void SetupOperator(OperatorType &op) {
    this->op_ = &op;
    this->SetModifiedOperator(true);

    if (this->hypre_precond_ != nullptr) {
      if (this->hypre_precond_->GetUseSolverOperator()) {
        this->hypre_precond_->SetupOperator(op);
        this->hypre_precond_->SetUseSolverOperator(true);
      }
    }
  }

  /// Applies the preconditioner which is possibly an inexact linear solver.
  /// @param b right hand side vector
  /// @param x solution vector
  /// @return status if preconditioning succeeded
  virtual LinearSolverState Solve(const VectorType &b, VectorType *x) {
    assert(x->is_initialized());
    assert(b.is_initialized());

    if (!this->GetReuse() || !this->initialized_) {
      this->Init();
    }

    if (this->hypre_precond_ != nullptr) {
      if (!this->GetReuse() || !this->hypre_precond_->GetReuse() ||
          !this->GetState() || !this->hypre_precond_->GetState()) {
        this->Build(&b, x);
      }
    } else {
      if (!this->GetReuse() || !this->GetState()) {
        this->Build(&b, x);
      }
    }

    if (this->print_level_ > 2) {
      LOG_INFO(this->name_, "Solve");
    }

    Timer timer;
    timer.start();
    LinearSolverState state = this->SolveImpl(b, x);
    timer.stop();

    this->num_solve_++;
    this->time_solve_ += timer.get_duration();
    this->acc_iter_ += this->iter();

    if (this->info_ != nullptr) {
      this->info_->add("iter", this->iter_);
      this->info_->add("time", timer.get_duration());
    }

    return state;
  }

  /// Sets up a Hypre preconditioner for the linear solver.
  /// @param precond preconditioner

  virtual void SetupPreconditioner(HyprePreconditioner< LAD > &precond) {
    this->hypre_precond_ = &precond;
  }

  /// Return preconditioner

  virtual HyprePreconditioner< LAD > *GetPreconditioner() {
    return this->hypre_precond_;
  }

protected:
  /// Preconditioner
  HyprePreconditioner< LAD > *hypre_precond_;
};

} // namespace la
} // namespace hiflow

#endif // HIFLOW_LINEARSOLVER_HYPRELINEAR_SOLVER_H_
