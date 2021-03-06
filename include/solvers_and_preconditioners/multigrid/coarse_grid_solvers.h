/*
 * mg_coarse_grid_solvers.h
 *
 *  Created on: Nov 23, 2016
 *      Author: fehn
 */

#ifndef INCLUDE_SOLVERS_AND_PRECONDITIONERS_MGCOARSEGRIDSOLVERS_H_
#define INCLUDE_SOLVERS_AND_PRECONDITIONERS_MGCOARSEGRIDSOLVERS_H_

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/multigrid/mg_base.h>

#include "../preconditioner/block_jacobi_preconditioner.h"
#include "../preconditioner/jacobi_preconditioner.h"
#include "../preconditioner/preconditioner_base.h"

#include "../solvers/iterative_solvers_dealii_wrapper.h"
#include "../solvers/solver_data.h"

#include "../../functionalities/set_zero_mean_value.h"
#include "../preconditioner/preconditioner_amg.h"

enum class KrylovSolverType
{
  CG,
  GMRES
};

template<typename Operator>
class MGCoarseKrylov
  : public MGCoarseGridBase<LinearAlgebra::distributed::Vector<typename Operator::value_type>>
{
public:
  typedef double TrilinosNumber;

  typedef typename Operator::value_type MultigridNumber;

  typedef LinearAlgebra::distributed::Vector<MultigridNumber> VectorType;

  typedef LinearAlgebra::distributed::Vector<TrilinosNumber> VectorTypeTrilinos;

  struct AdditionalData
  {
    /**
     * Constructor.
     */
    AdditionalData()
      : solver_type(KrylovSolverType::CG),
        solver_data(SolverData(1e4, 1.e-12, 1.e-3, 100)),
        operator_is_singular(false),
        preconditioner(MultigridCoarseGridPreconditioner::None),
        amg_data(AMGData())
    {
    }

    // Type of Krylov solver
    KrylovSolverType solver_type;

    // Solver data
    SolverData solver_data;

    // in case of singular operators (with constant vectors forming the nullspace) the rhs vector
    // has to be projected onto the space of vectors with zero mean prior to solving the coarse
    // grid problem
    bool operator_is_singular;

    // Preconditioner
    MultigridCoarseGridPreconditioner preconditioner;

    // Configuration of AMG settings
    AMGData amg_data;
  };

  MGCoarseKrylov(Operator const & matrix, AdditionalData const & additional_data)
    : coarse_matrix(matrix), additional_data(additional_data)
  {
    if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::PointJacobi)
    {
      preconditioner.reset(new JacobiPreconditioner<Operator>(coarse_matrix));
      std::shared_ptr<JacobiPreconditioner<Operator>> jacobi =
        std::dynamic_pointer_cast<JacobiPreconditioner<Operator>>(preconditioner);
      AssertDimension(jacobi->get_size_of_diagonal(), coarse_matrix.m());
    }
    else if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::BlockJacobi)
    {
      preconditioner.reset(new BlockJacobiPreconditioner<Operator>(coarse_matrix));
    }
    else if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::AMG)
    {
      preconditioner_trilinos.reset(
        new PreconditionerAMG<Operator, TrilinosNumber>(matrix, additional_data.amg_data));
    }
    else
    {
      AssertThrow(
        additional_data.preconditioner == MultigridCoarseGridPreconditioner::None ||
          additional_data.preconditioner == MultigridCoarseGridPreconditioner::PointJacobi ||
          additional_data.preconditioner == MultigridCoarseGridPreconditioner::BlockJacobi ||
          additional_data.preconditioner == MultigridCoarseGridPreconditioner::AMG,
        ExcMessage("Specified preconditioner for PCG coarse grid solver not implemented."));
    }
  }

  virtual ~MGCoarseKrylov()
  {
  }

  void
  update()
  {
    if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::None)
    {
      // do nothing
    }
    else if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::PointJacobi ||
            additional_data.preconditioner == MultigridCoarseGridPreconditioner::BlockJacobi)
    {
      preconditioner->update();
    }
    else if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::AMG)
    {
      preconditioner_trilinos->update();
    }
    else
    {
      AssertThrow(false, ExcMessage("Not implemented."));
    }
  }

  virtual void
  operator()(unsigned int const, VectorType & dst, VectorType const & src) const
  {
    VectorType r(src);
    if(additional_data.operator_is_singular)
      set_zero_mean_value(r);

    if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::AMG)
    {
#ifdef DEAL_II_WITH_TRILINOS
      // create temporal vectors of type TrilinosNumber (double)
      VectorTypeTrilinos dst_tri;
      dst_tri.reinit(dst, false);
      VectorTypeTrilinos src_tri;
      src_tri.reinit(r, true);
      src_tri.copy_locally_owned_data_from(r);

      ReductionControl solver_control(additional_data.solver_data.max_iter,
                                      additional_data.solver_data.abs_tol,
                                      additional_data.solver_data.rel_tol);

      if(additional_data.solver_type == KrylovSolverType::CG)
      {
        SolverCG<VectorTypeTrilinos>                                 solver(solver_control);
        std::shared_ptr<PreconditionerAMG<Operator, TrilinosNumber>> coarse_operator =
          std::dynamic_pointer_cast<PreconditionerAMG<Operator, TrilinosNumber>>(
            preconditioner_trilinos);
        solver.solve(coarse_operator->system_matrix, dst_tri, src_tri, *preconditioner_trilinos);
      }
      else if(additional_data.solver_type == KrylovSolverType::GMRES)
      {
        typename SolverGMRES<VectorTypeTrilinos>::AdditionalData gmres_data;
        gmres_data.max_n_tmp_vectors     = additional_data.solver_data.max_krylov_size;
        gmres_data.right_preconditioning = true;

        SolverGMRES<VectorTypeTrilinos> solver(solver_control, gmres_data);
        std::shared_ptr<PreconditionerAMG<Operator, TrilinosNumber>> coarse_operator =
          std::dynamic_pointer_cast<PreconditionerAMG<Operator, TrilinosNumber>>(
            preconditioner_trilinos);
        solver.solve(coarse_operator->system_matrix, dst_tri, src_tri, *preconditioner_trilinos);
      }
      else
      {
        AssertThrow(false, ExcMessage("Not implemented."));
      }

      // convert TrilinosNumber (double) -> MultigridNumber (float)
      dst.copy_locally_owned_data_from(dst_tri);
#endif
    }
    else
    {
      std::shared_ptr<IterativeSolverBase<VectorType>> solver;

      if(additional_data.solver_type == KrylovSolverType::CG)
      {
        CGSolverData solver_data;
        solver_data.max_iter             = additional_data.solver_data.max_iter;
        solver_data.solver_tolerance_abs = additional_data.solver_data.abs_tol;
        solver_data.solver_tolerance_rel = additional_data.solver_data.rel_tol;

        if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::None)
        {
          solver_data.use_preconditioner = false;
        }
        else if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::PointJacobi ||
                additional_data.preconditioner == MultigridCoarseGridPreconditioner::BlockJacobi)
        {
          solver_data.use_preconditioner = true;
        }
        else
        {
          AssertThrow(false, ExcMessage("Not implemented."));
        }

        solver.reset(new CGSolver<Operator, PreconditionerBase<MultigridNumber>, VectorType>(
          coarse_matrix, *preconditioner, solver_data));
      }
      else if(additional_data.solver_type == KrylovSolverType::GMRES)
      {
        GMRESSolverData solver_data;

        solver_data.max_iter             = additional_data.solver_data.max_iter;
        solver_data.solver_tolerance_abs = additional_data.solver_data.abs_tol;
        solver_data.solver_tolerance_rel = additional_data.solver_data.rel_tol;
        solver_data.max_n_tmp_vectors    = additional_data.solver_data.max_krylov_size;

        if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::None)
        {
          solver_data.use_preconditioner = false;
        }
        else if(additional_data.preconditioner == MultigridCoarseGridPreconditioner::PointJacobi ||
                additional_data.preconditioner == MultigridCoarseGridPreconditioner::BlockJacobi)
        {
          solver_data.use_preconditioner = true;
        }
        else
        {
          AssertThrow(false, ExcMessage("Not implemented."));
        }

        solver.reset(new GMRESSolver<Operator, PreconditionerBase<MultigridNumber>, VectorType>(
          coarse_matrix, *preconditioner, solver_data));
      }
      else
      {
        AssertThrow(false, ExcMessage("Not implemented."));
      }

      // Note that the preconditioner has already been updated
      solver->solve(dst, src, false);
    }
  }

private:
  const Operator & coarse_matrix;

  std::shared_ptr<PreconditionerBase<MultigridNumber>> preconditioner;

  // we need a separate object here because Trilinos needs double precision
  std::shared_ptr<PreconditionerBase<TrilinosNumber>> preconditioner_trilinos;

  AdditionalData additional_data;
};


template<typename Vector, typename InverseOperator>
class MGCoarseChebyshev : public MGCoarseGridBase<Vector>
{
public:
  MGCoarseChebyshev(std::shared_ptr<InverseOperator const> inverse) : inverse_operator(inverse)
  {
  }

  virtual ~MGCoarseChebyshev()
  {
  }

  virtual void
  operator()(const unsigned int level, Vector & dst, const Vector & src) const
  {
    AssertThrow(inverse_operator.get() != 0,
                ExcMessage("MGCoarseChebyshev: inverse_operator is not initialized."));

    AssertThrow(level == 0, ExcNotImplemented());

    inverse_operator->vmult(dst, src);
  }

  std::shared_ptr<InverseOperator const> inverse_operator;
};

template<typename Operator>
class MGCoarseAMG
  : public MGCoarseGridBase<LinearAlgebra::distributed::Vector<typename Operator::value_type>>
{
private:
  typedef double TrilinosNumber;

  typedef LinearAlgebra::distributed::Vector<TrilinosNumber> VectorTypeTrilinos;

  typedef LinearAlgebra::distributed::Vector<typename Operator::value_type> VectorTypeMultigrid;

public:
  MGCoarseAMG(Operator const & op, AMGData data = AMGData())
  {
    amg_preconditioner.reset(new PreconditionerAMG<Operator, TrilinosNumber>(op, data));
  }

  void
  update()
  {
    amg_preconditioner->update();
  }

  void
  operator()(unsigned int const /*level*/,
             VectorTypeMultigrid &       dst,
             VectorTypeMultigrid const & src) const
  {
    // create temporal vectors of type VectorTypeTrilinos (double)
    VectorTypeTrilinos dst_trilinos;
    dst_trilinos.reinit(dst, false);
    VectorTypeTrilinos src_trilinos;
    src_trilinos.reinit(src, true);

    // convert: VectorTypeMultigrid -> VectorTypeTrilinos
    src_trilinos.copy_locally_owned_data_from(src);

    // use Trilinos to perform AMG
    amg_preconditioner->vmult(dst_trilinos, src_trilinos);

    // convert: VectorTypeTrilinos -> VectorTypeMultigrid
    dst.copy_locally_owned_data_from(dst_trilinos);
  }

private:
  std::shared_ptr<PreconditionerAMG<Operator, TrilinosNumber>> amg_preconditioner;
};


#endif /* INCLUDE_SOLVERS_AND_PRECONDITIONERS_MGCOARSEGRIDSOLVERS_H_ */
