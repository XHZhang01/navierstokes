/*
 * postprocessor.cpp
 *
 *  Created on: May 13, 2019
 *      Author: fehn
 */

#include "postprocessor.h"

namespace ConvDiff
{
template<int dim, typename Number>
PostProcessor<dim, Number>::PostProcessor(PostProcessorData<dim> const & pp_data_in)
  : pp_data(pp_data_in)
{
}

template<int dim, typename Number>
void
PostProcessor<dim, Number>::setup(DoFHandler<dim> const & dof_handler, Mapping<dim> const & mapping)
{
  error_calculator.setup(dof_handler, mapping, pp_data.error_data);

  output_generator.setup(dof_handler, mapping, pp_data.output_data);
}

template<int dim, typename Number>
void
PostProcessor<dim, Number>::do_postprocessing(VectorType const & solution,
                                              double const       time,
                                              int const          time_step_number)
{
  error_calculator.evaluate(solution, time, time_step_number);

  output_generator.evaluate(solution, time, time_step_number);
}

template class PostProcessor<2, float>;
template class PostProcessor<3, float>;

template class PostProcessor<2, double>;
template class PostProcessor<3, double>;

} // namespace ConvDiff
