/*
 * postprocessor.h
 *
 *  Created on: Oct 12, 2016
 *      Author: fehn
 */

#ifndef INCLUDE_CONVECTION_DIFFUSION_POSTPROCESSOR_H_
#define INCLUDE_CONVECTION_DIFFUSION_POSTPROCESSOR_H_

// deal.II
#include <deal.II/lac/la_parallel_vector.h>

#include "convection_diffusion/user_interface/analytical_solution.h"
#include "output_generator.h"
#include "postprocessor/error_calculation.h"
#include "postprocessor/output_data.h"
#include "postprocessor_base.h"

namespace ConvDiff
{
template<int dim>
struct PostProcessorData
{
  PostProcessorData()
  {
  }

  OutputDataBase            output_data;
  ErrorCalculationData<dim> error_data;
};

template<int dim, typename Number>
class PostProcessor : public PostProcessorBase<dim, Number>
{
private:
  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

public:
  PostProcessor(PostProcessorData<dim> const & pp_data_in);

  void
  setup(DoFHandler<dim> const & dof_handler, Mapping<dim> const & mapping);

  void
  do_postprocessing(VectorType const & solution,
                    double const       time             = 0.0,
                    int const          time_step_number = -1);

private:
  PostProcessorData<dim> pp_data;

  OutputGenerator<dim, Number> output_generator;
  ErrorCalculator<dim, Number> error_calculator;
};

} // namespace ConvDiff


#endif /* INCLUDE_CONVECTION_DIFFUSION_POSTPROCESSOR_H_ */
