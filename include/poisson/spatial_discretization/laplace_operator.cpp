#include "laplace_operator.h"

#include "weak_boundary_conditions.h"

namespace Poisson
{
template<int dim, typename Number>
void
LaplaceOperator<dim, Number>::reinit(MatrixFree<dim, Number> const &   matrix_free,
                                     AffineConstraints<double> const & constraint_matrix,
                                     LaplaceOperatorData<dim> const &  data)
{
  Base::reinit(matrix_free, constraint_matrix, data);

  kernel.reinit(matrix_free, data.kernel_data, data.dof_index);

  this->integrator_flags = kernel.get_integrator_flags();
}

template<int dim, typename Number>
void
LaplaceOperator<dim, Number>::reinit_face(unsigned int const face) const
{
  Base::reinit_face(face);

  kernel.reinit_face(*this->integrator_m, *this->integrator_p);
}

template<int dim, typename Number>
void
LaplaceOperator<dim, Number>::reinit_boundary_face(unsigned int const face) const
{
  Base::reinit_boundary_face(face);

  kernel.reinit_boundary_face(*this->integrator_m);
}

template<int dim, typename Number>
void
LaplaceOperator<dim, Number>::reinit_face_cell_based(unsigned int const       cell,
                                                     unsigned int const       face,
                                                     types::boundary_id const boundary_id) const
{
  Base::reinit_face_cell_based(cell, face, boundary_id);

  kernel.reinit_face_cell_based(boundary_id, *this->integrator_m, *this->integrator_p);
}

template<int dim, typename Number>
void
LaplaceOperator<dim, Number>::do_cell_integral(IntegratorCell & integrator) const
{
  for(unsigned int q = 0; q < integrator.n_q_points; ++q)
  {
    integrator.submit_gradient(kernel.get_volume_flux(integrator, q), q);
  }
}

template<int dim, typename Number>
void
LaplaceOperator<dim, Number>::do_face_integral(IntegratorFace & integrator_m,
                                               IntegratorFace & integrator_p) const
{
  for(unsigned int q = 0; q < integrator_m.n_q_points; ++q)
  {
    scalar value_m = integrator_m.get_value(q);
    scalar value_p = integrator_p.get_value(q);

    scalar gradient_flux = kernel.calculate_gradient_flux(value_m, value_p);

    scalar normal_gradient_m = integrator_m.get_normal_derivative(q);
    scalar normal_gradient_p = integrator_p.get_normal_derivative(q);

    scalar value_flux =
      kernel.calculate_value_flux(normal_gradient_m, normal_gradient_p, value_m, value_p);

    integrator_m.submit_normal_derivative(gradient_flux, q);
    integrator_p.submit_normal_derivative(gradient_flux, q);

    integrator_m.submit_value(-value_flux, q);
    integrator_p.submit_value(value_flux, q); // + sign since n??? = -n???
  }
}

template<int dim, typename Number>
void
LaplaceOperator<dim, Number>::do_face_int_integral(IntegratorFace & integrator_m,
                                                   IntegratorFace & integrator_p) const
{
  (void)integrator_p;

  for(unsigned int q = 0; q < integrator_m.n_q_points; ++q)
  {
    // set exterior value to zero
    scalar value_m = integrator_m.get_value(q);
    scalar value_p = make_vectorized_array<Number>(0.0);

    scalar gradient_flux = kernel.calculate_gradient_flux(value_m, value_p);

    // set exterior value to zero
    scalar normal_gradient_m = integrator_m.get_normal_derivative(q);
    scalar normal_gradient_p = make_vectorized_array<Number>(0.0);

    scalar value_flux =
      kernel.calculate_value_flux(normal_gradient_m, normal_gradient_p, value_m, value_p);

    integrator_m.submit_normal_derivative(gradient_flux, q);
    integrator_m.submit_value(-value_flux, q);
  }
}

template<int dim, typename Number>
void
LaplaceOperator<dim, Number>::do_face_ext_integral(IntegratorFace & integrator_m,
                                                   IntegratorFace & integrator_p) const
{
  (void)integrator_m;

  for(unsigned int q = 0; q < integrator_p.n_q_points; ++q)
  {
    // set value_m to zero
    scalar value_p = integrator_p.get_value(q);
    scalar value_m = make_vectorized_array<Number>(0.0);

    scalar gradient_flux = kernel.calculate_gradient_flux(value_p, value_m);

    // set gradient_m to zero
    scalar normal_gradient_m = make_vectorized_array<Number>(0.0);
    // minus sign to get the correct normal vector n??? = -n???
    scalar normal_gradient_p = -integrator_p.get_normal_derivative(q);

    scalar value_flux =
      kernel.calculate_value_flux(normal_gradient_p, normal_gradient_m, value_p, value_m);

    integrator_p.submit_normal_derivative(-gradient_flux, q); // opposite sign since n??? = -n???
    integrator_p.submit_value(-value_flux, q);
  }
}

template<int dim, typename Number>
void
LaplaceOperator<dim, Number>::do_boundary_integral(IntegratorFace &           integrator_m,
                                                   OperatorType const &       operator_type,
                                                   types::boundary_id const & boundary_id) const
{
  BoundaryType boundary_type = this->data.bc->get_boundary_type(boundary_id);

  for(unsigned int q = 0; q < integrator_m.n_q_points; ++q)
  {
    scalar value_m = calculate_interior_value(q, integrator_m, operator_type);
    scalar value_p = calculate_exterior_value(value_m,
                                              q,
                                              integrator_m,
                                              operator_type,
                                              boundary_type,
                                              boundary_id,
                                              this->data.bc,
                                              this->time);

    scalar gradient_flux = kernel.calculate_gradient_flux(value_m, value_p);

    scalar normal_gradient_m = calculate_interior_normal_gradient(q, integrator_m, operator_type);
    scalar normal_gradient_p = calculate_exterior_normal_gradient(normal_gradient_m,
                                                                  q,
                                                                  integrator_m,
                                                                  operator_type,
                                                                  boundary_type,
                                                                  boundary_id,
                                                                  this->data.bc,
                                                                  this->time);

    scalar value_flux =
      kernel.calculate_value_flux(normal_gradient_m, normal_gradient_p, value_m, value_p);

    integrator_m.submit_normal_derivative(gradient_flux, q);
    integrator_m.submit_value(-value_flux, q);
  }
}

template<int dim, typename Number>
void
LaplaceOperator<dim, Number>::do_verify_boundary_conditions(
  types::boundary_id const             boundary_id,
  LaplaceOperatorData<dim> const &     data,
  std::set<types::boundary_id> const & periodic_boundary_ids) const
{
  unsigned int counter = 0;
  if(data.bc->dirichlet_bc.find(boundary_id) != data.bc->dirichlet_bc.end())
    counter++;

  if(data.bc->neumann_bc.find(boundary_id) != data.bc->neumann_bc.end())
    counter++;

  if(periodic_boundary_ids.find(boundary_id) != periodic_boundary_ids.end())
    counter++;

  AssertThrow(counter == 1, ExcMessage("Boundary face with non-unique boundary type found."));
}

template class LaplaceOperator<2, float>;
template class LaplaceOperator<2, double>;

template class LaplaceOperator<3, float>;
template class LaplaceOperator<3, double>;

} // namespace Poisson
