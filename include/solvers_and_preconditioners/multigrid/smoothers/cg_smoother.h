/*
 * cg_smoother.h
 *
 *  Created on: Mar 24, 2017
 *      Author: fehn
 */

#ifndef INCLUDE_SOLVERS_AND_PRECONDITIONERS_CG_SMOOTHER_H_
#define INCLUDE_SOLVERS_AND_PRECONDITIONERS_CG_SMOOTHER_H_

// deal.II
#include <deal.II/lac/solver_cg.h>

#include "smoother_base.h"

#include "../../preconditioner/block_jacobi_preconditioner.h"
#include "../../preconditioner/jacobi_preconditioner.h"

#include "../multigrid_input_parameters.h"

template<typename Operator, typename VectorType>
class CGSmoother : public SmootherBase<VectorType>
{
public:
  CGSmoother() : underlying_operator(nullptr), preconditioner(nullptr)
  {
  }

  ~CGSmoother()
  {
    delete preconditioner;
    preconditioner = nullptr;
  }

  CGSmoother(CGSmoother const &) = delete;

  CGSmoother &
  operator=(CGSmoother const &) = delete;

  struct AdditionalData
  {
    /**
     * Constructor.
     */
    AdditionalData() : preconditioner(PreconditionerSmoother::None), number_of_iterations(5)
    {
    }

    // preconditioner
    PreconditionerSmoother preconditioner;

    // number of CG iterations per smoothing step
    unsigned int number_of_iterations;
  };

  void
  initialize(Operator & operator_in, AdditionalData const & additional_data_in)
  {
    underlying_operator = &operator_in;
    data                = additional_data_in;

    if(data.preconditioner == PreconditionerSmoother::PointJacobi)
    {
      preconditioner = new JacobiPreconditioner<Operator>(*underlying_operator);
    }
    else if(data.preconditioner == PreconditionerSmoother::BlockJacobi)
    {
      preconditioner = new BlockJacobiPreconditioner<Operator>(*underlying_operator);
    }
    else
    {
      AssertThrow(data.preconditioner == PreconditionerSmoother::None,
                  ExcMessage("Specified preconditioner not implemented for CG smoother"));
    }
  }

  void
  update()
  {
    if(preconditioner != nullptr)
      preconditioner->update();
  }

  // same as step(), but sets dst-vector to zero
  void
  vmult(VectorType & dst, VectorType const & src) const
  {
    dst = 0.0;
    step(dst, src);
  }

  void
  step(VectorType & dst, VectorType const & src) const
  {
    IterationNumberControl control(data.number_of_iterations, 1.e-20, 1.e-10);

    SolverCG<VectorType> solver(control);

    if(preconditioner != nullptr)
      solver.solve(*underlying_operator, dst, src, *preconditioner);
    else
      solver.solve(*underlying_operator, dst, src, PreconditionIdentity());
  }

private:
  Operator *     underlying_operator;
  AdditionalData data;

  PreconditionerBase<typename Operator::value_type> * preconditioner;
};


#endif /* INCLUDE_SOLVERS_AND_PRECONDITIONERS_CG_SMOOTHER_H_ */
