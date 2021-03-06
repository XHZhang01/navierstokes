/*
 * rhs_operator.cpp
 *
 *  Created on: Dec 3, 2018
 *      Author: fehn
 */

#include "rhs_operator.h"

namespace ConvDiff
{
template<int dim, typename Number>
RHSOperator<dim, Number>::RHSOperator() : matrix_free(nullptr), time(0.0)
{
}

template<int dim, typename Number>
void
RHSOperator<dim, Number>::reinit(MatrixFree<dim, Number> const & matrix_free_in,
                                 RHSOperatorData<dim> const &    data_in)
{
  this->matrix_free = &matrix_free_in;
  this->data        = data_in;

  kernel.reinit(data.kernel_data);
}

template<int dim, typename Number>
void
RHSOperator<dim, Number>::evaluate(VectorType & dst, double const evaluation_time) const
{
  dst = 0;
  evaluate_add(dst, evaluation_time);
}

template<int dim, typename Number>
void
RHSOperator<dim, Number>::evaluate_add(VectorType & dst, double const evaluation_time) const
{
  this->time = evaluation_time;

  VectorType src;
  matrix_free->cell_loop(&This::cell_loop, this, dst, src);
}

template<int dim, typename Number>
void
RHSOperator<dim, Number>::do_cell_integral(IntegratorCell & integrator) const
{
  for(unsigned int q = 0; q < integrator.n_q_points; ++q)
  {
    integrator.submit_value(kernel.get_volume_flux(integrator, q, time), q);
  }

  integrator.integrate(true, false);
}

template<int dim, typename Number>
void
RHSOperator<dim, Number>::cell_loop(MatrixFree<dim, Number> const & matrix_free,
                                    VectorType &                    dst,
                                    VectorType const &              src,
                                    Range const &                   cell_range) const
{
  (void)src;

  IntegratorCell integrator(matrix_free, data.dof_index, data.quad_index);

  for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
  {
    integrator.reinit(cell);

    do_cell_integral(integrator);

    integrator.distribute_local_to_global(dst);
  }
}

template class RHSOperator<2, float>;
template class RHSOperator<2, double>;

template class RHSOperator<3, float>;
template class RHSOperator<3, double>;

} // namespace ConvDiff
