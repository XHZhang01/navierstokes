/*
 * comp_navier_stokes_operators.h
 *
 *  Created on: 2018
 *      Author: fehn
 */

#ifndef INCLUDE_COMPRESSIBLE_NAVIER_STOKES_SPATIAL_DISCRETIZATION_COMP_NAVIER_STOKES_OPERATORS_H_
#define INCLUDE_COMPRESSIBLE_NAVIER_STOKES_SPATIAL_DISCRETIZATION_COMP_NAVIER_STOKES_OPERATORS_H_

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/fe_evaluation_notemplate.h>

#include "../include/functionalities/evaluate_functions.h"
#include "operators/interior_penalty_parameter.h"

#include "../../compressible_navier_stokes/user_interface/boundary_descriptor.h"
#include "../../compressible_navier_stokes/user_interface/input_parameters.h"


#include <iostream>

namespace CompNS
{
template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
    VectorizedArray<Number>
    calculate_pressure(Tensor<1, dim, VectorizedArray<Number>> const & rho_u,
                       Tensor<1, dim, VectorizedArray<Number>> const & u,
                       VectorizedArray<Number> const &                 rho_E,
                       Number const &                                  gamma)
{
  return (gamma - 1.0) * (rho_E - 0.5 * scalar_product(rho_u, u));
}

template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
  VectorizedArray<Number>
  calculate_pressure(VectorizedArray<Number> const &                 rho,
                     Tensor<1, dim, VectorizedArray<Number>> const & u,
                     VectorizedArray<Number> const &                 E,
                     Number const &                                  gamma)
{
  return (gamma - 1.0) * rho * (E - 0.5 * scalar_product(u, u));
}

template<typename Number>
inline DEAL_II_ALWAYS_INLINE //
  VectorizedArray<Number>
  calculate_temperature(VectorizedArray<Number> const & p,
                        VectorizedArray<Number> const & rho,
                        Number const &                  R)
{
  return p / (rho * R);
}

template<int dim, typename Number>
inline VectorizedArray<Number>
calculate_energy(VectorizedArray<Number> const &                 T,
                 Tensor<1, dim, VectorizedArray<Number>> const & u,
                 Number const &                                  c_v)
{
  return c_v * T + 0.5 * scalar_product(u, u);
}

template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
  Tensor<1, dim, VectorizedArray<Number>>
  calculate_grad_E(VectorizedArray<Number> const &                 rho_inverse,
                   VectorizedArray<Number> const &                 rho_E,
                   Tensor<1, dim, VectorizedArray<Number>> const & grad_rho,
                   Tensor<1, dim, VectorizedArray<Number>> const & grad_rho_E)
{
  VectorizedArray<Number> E = rho_inverse * rho_E;

  return rho_inverse * (grad_rho_E - E * grad_rho);
}

template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
  Tensor<2, dim, VectorizedArray<Number>>
  calculate_grad_u(VectorizedArray<Number> const &                 rho_inverse,
                   Tensor<1, dim, VectorizedArray<Number>> const & rho_u,
                   Tensor<1, dim, VectorizedArray<Number>> const & grad_rho,
                   Tensor<2, dim, VectorizedArray<Number>> const & grad_rho_u)
{
  Tensor<2, dim, VectorizedArray<Number>> out;
  for(unsigned int d = 0; d < dim; ++d)
  {
    VectorizedArray<Number> ud = rho_inverse * rho_u[d];
    for(unsigned int e = 0; e < dim; ++e)
      out[d][e] = rho_inverse * (grad_rho_u[d][e] - ud * grad_rho[e]);
  }
  return out;
}

template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim, VectorizedArray<Number>>
    calculate_grad_T(Tensor<1, dim, VectorizedArray<Number>> const & grad_E,
                     Tensor<1, dim, VectorizedArray<Number>> const & u,
                     Tensor<2, dim, VectorizedArray<Number>> const & grad_u,
                     Number const &                                  gamma,
                     Number const &                                  R)
{
  return (gamma - 1.0) / R * (grad_E - u * grad_u);
}

template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
    Tensor<2, dim, VectorizedArray<Number>>
    calculate_stress_tensor(Tensor<2, dim, VectorizedArray<Number>> const & grad_u,
                            Number const &                                  viscosity)
{
  Tensor<2, dim, VectorizedArray<Number>> out;
  const VectorizedArray<Number>           divu = (2. / 3.) * trace(grad_u);
  for(unsigned int d = 0; d < dim; ++d)
  {
    for(unsigned int e = 0; e < dim; ++e)
      out[d][e] = viscosity * (grad_u[d][e] + grad_u[e][d]);
    out[d][d] -= viscosity * divu;
  }
  return out;
}

/*
 * Calculate exterior state "+" for a scalar quantity depending on interior state "-" and boundary
 * conditions.
 */
template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
  VectorizedArray<Number>
  calculate_exterior_value(
    VectorizedArray<Number> const &                          value_m,
    BoundaryType const &                                     boundary_type,
    std::shared_ptr<CompNS::BoundaryDescriptor<dim>> const & boundary_descriptor,
    types::boundary_id const &                               boundary_id,
    Point<dim, VectorizedArray<Number>> const &              q_point,
    Number const &                                           time)
{
  VectorizedArray<Number> value_p = make_vectorized_array<Number>(0.);

  if(boundary_type == BoundaryType::Dirichlet)
  {
    typename std::map<types::boundary_id, std::shared_ptr<Function<dim>>>::iterator it =
      boundary_descriptor->dirichlet_bc.find(boundary_id);
    VectorizedArray<Number> g = evaluate_scalar_function(it->second, q_point, time);

    value_p = -value_m + make_vectorized_array<Number>(2.0) * g;
  }
  else if(boundary_type == BoundaryType::Neumann)
  {
    value_p = value_m;
  }
  else
  {
    AssertThrow(false, ExcMessage("Not implemented."));
  }

  return value_p;
}

/*
 * Calculate exterior state "+" for a scalar quantity depending on interior state "-" and boundary
 * conditions.
 */
template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim, VectorizedArray<Number>>
    calculate_exterior_value(
      Tensor<1, dim, VectorizedArray<Number>> const &          value_m,
      BoundaryType const                                       boundary_type,
      std::shared_ptr<CompNS::BoundaryDescriptor<dim>> const & boundary_descriptor,
      types::boundary_id const &                               boundary_id,
      Point<dim, VectorizedArray<Number>> const &              q_point,
      Number const &                                           time)
{
  Tensor<1, dim, VectorizedArray<Number>> value_p;

  if(boundary_type == BoundaryType::Dirichlet)
  {
    typename std::map<types::boundary_id, std::shared_ptr<Function<dim>>>::iterator it =
      boundary_descriptor->dirichlet_bc.find(boundary_id);
    Tensor<1, dim, VectorizedArray<Number>> g =
      evaluate_vectorial_function(it->second, q_point, time);

    value_p = -value_m + make_vectorized_array<Number>(2.0) * g;
  }
  else if(boundary_type == BoundaryType::Neumann)
  {
    value_p = value_m;
  }
  else
  {
    AssertThrow(false, ExcMessage("Not implemented."));
  }

  return value_p;
}

/*
 * Calculate exterior state of normal stresses depending on interior data and boundary conditions.
 */
template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim, VectorizedArray<Number>>
    calculate_exterior_normal_grad(
      Tensor<1, dim, VectorizedArray<Number>> const &          tau_M_normal,
      BoundaryType const &                                     boundary_type,
      std::shared_ptr<CompNS::BoundaryDescriptor<dim>> const & boundary_descriptor,
      types::boundary_id const &                               boundary_id,
      Point<dim, VectorizedArray<Number>> const &              q_point,
      Number const &                                           time)
{
  Tensor<1, dim, VectorizedArray<Number>> tau_P_normal;

  if(boundary_type == BoundaryType::Dirichlet)
  {
    // do nothing
    tau_P_normal = tau_M_normal;
  }
  else if(boundary_type == BoundaryType::Neumann)
  {
    typename std::map<types::boundary_id, std::shared_ptr<Function<dim>>>::iterator it =
      boundary_descriptor->neumann_bc.find(boundary_id);
    Tensor<1, dim, VectorizedArray<Number>> h =
      evaluate_vectorial_function(it->second, q_point, time);

    tau_P_normal = -tau_M_normal + make_vectorized_array<Number>(2.0) * h;
  }
  else
  {
    AssertThrow(false, ExcMessage("Not implemented."));
  }

  return tau_P_normal;
}

/*
 * Calculate exterior temperature gradient depending on interior data and boundary conditions.
 */
template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
  VectorizedArray<Number>
  calculate_exterior_normal_grad(
    VectorizedArray<Number> const &                          grad_T_M_normal,
    BoundaryType const &                                     boundary_type,
    std::shared_ptr<CompNS::BoundaryDescriptor<dim>> const & boundary_descriptor,
    types::boundary_id const &                               boundary_id,
    Point<dim, VectorizedArray<Number>> const &              q_point,
    Number const &                                           time)
{
  VectorizedArray<Number> grad_T_P_normal = make_vectorized_array<Number>(0.0);

  if(boundary_type == BoundaryType::Dirichlet)
  {
    // do nothing
    grad_T_P_normal = grad_T_M_normal;
  }
  else if(boundary_type == BoundaryType::Neumann)
  {
    typename std::map<types::boundary_id, std::shared_ptr<Function<dim>>>::iterator it =
      boundary_descriptor->neumann_bc.find(boundary_id);
    VectorizedArray<Number> h = evaluate_scalar_function(it->second, q_point, time);

    grad_T_P_normal = -grad_T_M_normal + make_vectorized_array<Number>(2.0) * h;
  }
  else
  {
    AssertThrow(false, ExcMessage("Not implemented."));
  }

  return grad_T_P_normal;
}

/*
 * This function calculates the Lax-Friedrichs flux for the momentum equation
 */
template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim, VectorizedArray<Number>>
    calculate_flux(Tensor<2, dim, VectorizedArray<Number>> const & momentum_flux_M,
                   Tensor<2, dim, VectorizedArray<Number>> const & momentum_flux_P,
                   Tensor<1, dim, VectorizedArray<Number>> const & rho_u_M,
                   Tensor<1, dim, VectorizedArray<Number>> const & rho_u_P,
                   VectorizedArray<Number> const &                 lambda,
                   Tensor<1, dim, VectorizedArray<Number>> const & normal)
{
  Tensor<1, dim, VectorizedArray<Number>> out;
  for(unsigned int d = 0; d < dim; ++d)
  {
    VectorizedArray<Number> sum = VectorizedArray<Number>();
    for(unsigned int e = 0; e < dim; ++e)
      sum += (momentum_flux_M[d][e] + momentum_flux_P[d][e]) * normal[e];
    out[d] = 0.5 * (sum + lambda * (rho_u_M[d] - rho_u_P[d]));
  }
  return out;
}

/*
 * This function calculates the Lax-Friedrichs flux for scalar quantities (density/energy)
 */
template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
    VectorizedArray<Number>
    calculate_flux(Tensor<1, dim, VectorizedArray<Number>> const & flux_M,
                   Tensor<1, dim, VectorizedArray<Number>> const & flux_P,
                   VectorizedArray<Number> const &                 value_M,
                   VectorizedArray<Number> const &                 value_P,
                   VectorizedArray<Number> const &                 lambda,
                   Tensor<1, dim, VectorizedArray<Number>> const & normal)
{
  Tensor<1, dim, VectorizedArray<Number>> average_flux = 0.5 * (flux_M + flux_P);

  return average_flux * normal + 0.5 * lambda * (value_M - value_P);
}

/*
 * Calculation of lambda for Lax-Friedrichs flux according to Hesthaven:
 *   lambda = max( |u_M| + sqrt(|gamma * p_M / rho_M|) , |u_P| + sqrt(|gamma * p_P / rho_P|) )
 */
template<int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE //
  VectorizedArray<Number>
  calculate_lambda(VectorizedArray<Number> const &                 rho_m,
                   VectorizedArray<Number> const &                 rho_p,
                   Tensor<1, dim, VectorizedArray<Number>> const & u_m,
                   Tensor<1, dim, VectorizedArray<Number>> const & u_p,
                   VectorizedArray<Number> const &                 p_m,
                   VectorizedArray<Number> const &                 p_p,
                   Number const &                                  gamma)
{
  VectorizedArray<Number> lambda_m = u_m.norm() + std::sqrt(std::abs(gamma * p_m / rho_m));
  VectorizedArray<Number> lambda_p = u_p.norm() + std::sqrt(std::abs(gamma * p_p / rho_p));

  return std::max(lambda_m, lambda_p);
}

template<int dim>
struct BodyForceOperatorData
{
  BodyForceOperatorData() : dof_index(0), quad_index(0)
  {
  }

  unsigned int dof_index;
  unsigned int quad_index;

  std::shared_ptr<Function<dim>> rhs_rho;
  std::shared_ptr<Function<dim>> rhs_u;
  std::shared_ptr<Function<dim>> rhs_E;
};

template<int dim, typename Number>
class BodyForceOperator
{
public:
  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

  typedef BodyForceOperator<dim, Number> This;

  typedef CellIntegrator<dim, 1, Number>   CellIntegratorScalar;
  typedef CellIntegrator<dim, dim, Number> CellIntegratorVector;

  typedef VectorizedArray<Number>                 scalar;
  typedef Tensor<1, dim, VectorizedArray<Number>> vector;
  typedef Tensor<2, dim, VectorizedArray<Number>> tensor;
  typedef Point<dim, VectorizedArray<Number>>     point;

  BodyForceOperator() : matrix_free(nullptr), eval_time(0.0)
  {
  }

  void
  initialize(MatrixFree<dim, Number> const &    matrix_free_in,
             BodyForceOperatorData<dim> const & data_in)
  {
    this->matrix_free = &matrix_free_in;
    this->data        = data_in;
  }

  void
  evaluate(VectorType & dst, VectorType const & src, double const evaluation_time) const
  {
    dst = 0;
    evaluate_add(dst, src, evaluation_time);
  }

  void
  evaluate_add(VectorType & dst, VectorType const & src, double const evaluation_time) const
  {
    this->eval_time = evaluation_time;

    matrix_free->cell_loop(&This::cell_loop, this, dst, src);
  }

  inline DEAL_II_ALWAYS_INLINE //
    std::tuple<scalar, vector, scalar>
    get_volume_flux(CellIntegratorScalar & density,
                    CellIntegratorVector & momentum,
                    unsigned int const     q) const
  {
    point  q_points = density.quadrature_point(q);
    scalar rho      = density.get_value(q);
    vector u        = momentum.get_value(q) / rho;

    scalar rhs_density  = evaluate_scalar_function(data.rhs_rho, q_points, eval_time);
    vector rhs_momentum = evaluate_vectorial_function(data.rhs_u, q_points, eval_time);
    scalar rhs_energy   = evaluate_scalar_function(data.rhs_E, q_points, eval_time);

    return std::make_tuple(rhs_density, rhs_momentum, rhs_momentum * u + rhs_energy);
  }

private:
  void
  cell_loop(MatrixFree<dim, Number> const &               matrix_free,
            VectorType &                                  dst,
            VectorType const &                            src,
            std::pair<unsigned int, unsigned int> const & cell_range) const
  {
    CellIntegratorScalar density(matrix_free, data.dof_index, data.quad_index, 0);
    CellIntegratorVector momentum(matrix_free, data.dof_index, data.quad_index, 1);
    CellIntegratorScalar energy(matrix_free, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      density.reinit(cell);
      density.gather_evaluate(src, true, false);

      momentum.reinit(cell);
      momentum.gather_evaluate(src, true, false);

      energy.reinit(cell);

      for(unsigned int q = 0; q < density.n_q_points; ++q)
      {
        std::tuple<scalar, vector, scalar> flux = get_volume_flux(density, momentum, q);

        density.submit_value(std::get<0>(flux), q);
        momentum.submit_value(std::get<1>(flux), q);
        energy.submit_value(std::get<2>(flux), q);
      }

      density.integrate_scatter(true, false, dst);
      momentum.integrate_scatter(true, false, dst);
      energy.integrate_scatter(true, false, dst);
    }
  }

  MatrixFree<dim, Number> const * matrix_free;

  BodyForceOperatorData<dim> data;

  double mutable eval_time;
};

struct MassMatrixOperatorData
{
  MassMatrixOperatorData() : dof_index(0), quad_index(0)
  {
  }

  unsigned int dof_index;
  unsigned int quad_index;
};

template<int dim, typename Number>
class MassMatrixOperator
{
public:
  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

  typedef MassMatrixOperator<dim, Number> This;

  typedef CellIntegrator<dim, 1, Number>   CellIntegratorScalar;
  typedef CellIntegrator<dim, dim, Number> CellIntegratorVector;

  MassMatrixOperator() : matrix_free(nullptr)
  {
  }

  void
  initialize(MatrixFree<dim, Number> const & matrix_free_in, MassMatrixOperatorData const & data_in)
  {
    this->matrix_free = &matrix_free_in;
    this->data        = data_in;
  }

  // apply matrix vector multiplication
  void
  apply(VectorType & dst, VectorType const & src) const
  {
    dst = 0;
    apply_add(dst, src);
  }

  void
  apply_add(VectorType & dst, VectorType const & src) const
  {
    matrix_free->cell_loop(&This::cell_loop, this, dst, src);
  }

private:
  void
  cell_loop(MatrixFree<dim, Number> const &               matrix_free,
            VectorType &                                  dst,
            VectorType const &                            src,
            std::pair<unsigned int, unsigned int> const & cell_range) const
  {
    CellIntegratorScalar density(matrix_free, data.dof_index, data.quad_index, 0);
    CellIntegratorVector momentum(matrix_free, data.dof_index, data.quad_index, 1);
    CellIntegratorScalar energy(matrix_free, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      density.reinit(cell);
      density.gather_evaluate(src, true, false);

      momentum.reinit(cell);
      momentum.gather_evaluate(src, true, false);

      energy.reinit(cell);
      energy.gather_evaluate(src, true, false);

      for(unsigned int q = 0; q < density.n_q_points; ++q)
      {
        density.submit_value(density.get_value(q), q);
        momentum.submit_value(momentum.get_value(q), q);
        energy.submit_value(energy.get_value(q), q);
      }

      density.integrate_scatter(true, false, dst);
      momentum.integrate_scatter(true, false, dst);
      energy.integrate_scatter(true, false, dst);
    }
  }

  MatrixFree<dim, Number> const * matrix_free;
  MassMatrixOperatorData          data;
};

template<int dim>
struct ConvectiveOperatorData
{
  ConvectiveOperatorData()
    : dof_index(0), quad_index(0), heat_capacity_ratio(1.4), specific_gas_constant(287.0)
  {
  }

  unsigned int dof_index;
  unsigned int quad_index;

  std::shared_ptr<CompNS::BoundaryDescriptor<dim>>       bc_rho;
  std::shared_ptr<CompNS::BoundaryDescriptor<dim>>       bc_u;
  std::shared_ptr<CompNS::BoundaryDescriptor<dim>>       bc_p;
  std::shared_ptr<CompNS::BoundaryDescriptorEnergy<dim>> bc_E;

  double heat_capacity_ratio;
  double specific_gas_constant;
};

template<int dim, typename Number>
class ConvectiveOperator
{
public:
  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

  typedef ConvectiveOperator<dim, Number> This;

  typedef CellIntegrator<dim, 1, Number>   CellIntegratorScalar;
  typedef FaceIntegrator<dim, 1, Number>   FaceIntegratorScalar;
  typedef CellIntegrator<dim, dim, Number> CellIntegratorVector;
  typedef FaceIntegrator<dim, dim, Number> FaceIntegratorVector;

  typedef VectorizedArray<Number>                 scalar;
  typedef Tensor<1, dim, VectorizedArray<Number>> vector;
  typedef Tensor<2, dim, VectorizedArray<Number>> tensor;
  typedef Point<dim, VectorizedArray<Number>>     point;

  ConvectiveOperator() : matrix_free(nullptr)
  {
  }

  void
  initialize(MatrixFree<dim, Number> const &     matrix_free_in,
             ConvectiveOperatorData<dim> const & data_in)
  {
    this->matrix_free = &matrix_free_in;
    this->data        = data_in;

    gamma = data.heat_capacity_ratio;
    R     = data.specific_gas_constant;
    c_v   = R / (gamma - 1.0);
  }

  void
  evaluate(VectorType & dst, VectorType const & src, Number const evaluation_time) const
  {
    dst = 0;
    evaluate_add(dst, src, evaluation_time);
  }

  void
  evaluate_add(VectorType & dst, VectorType const & src, Number const evaluation_time) const
  {
    this->eval_time = evaluation_time;

    matrix_free->loop(
      &This::cell_loop, &This::face_loop, &This::boundary_face_loop, this, dst, src);
  }

  void
  set_evaluation_time(double const & evaluation_time) const
  {
    eval_time = evaluation_time;
  }

  inline DEAL_II_ALWAYS_INLINE //
    std::tuple<vector, tensor, vector>
    get_volume_flux(CellIntegratorScalar & density,
                    CellIntegratorVector & momentum,
                    CellIntegratorScalar & energy,
                    unsigned int const     q) const
  {
    scalar rho_inv = 1.0 / density.get_value(q);
    vector rho_u   = momentum.get_value(q);
    scalar rho_E   = energy.get_value(q);
    vector u       = rho_inv * rho_u;
    scalar p       = calculate_pressure(rho_u, u, rho_E, gamma);

    tensor momentum_flux = outer_product(rho_u, u);
    for(unsigned int d = 0; d < dim; ++d)
      momentum_flux[d][d] += p;

    vector energy_flux = (rho_E + p) * u;

    return std::make_tuple(rho_u, momentum_flux, energy_flux);
  }

  inline DEAL_II_ALWAYS_INLINE //
    std::tuple<scalar, vector, scalar>
    get_flux(FaceIntegratorScalar & density_m,
             FaceIntegratorScalar & density_p,
             FaceIntegratorVector & momentum_m,
             FaceIntegratorVector & momentum_p,
             FaceIntegratorScalar & energy_m,
             FaceIntegratorScalar & energy_p,
             unsigned int const     q) const
  {
    vector normal = momentum_m.get_normal_vector(q);

    // get values
    scalar rho_M   = density_m.get_value(q);
    scalar rho_P   = density_p.get_value(q);
    vector rho_u_M = momentum_m.get_value(q);
    vector rho_u_P = momentum_p.get_value(q);
    vector u_M     = rho_u_M / rho_M;
    vector u_P     = rho_u_P / rho_P;
    scalar rho_E_M = energy_m.get_value(q);
    scalar rho_E_P = energy_p.get_value(q);

    // calculate pressure
    scalar p_M = calculate_pressure(rho_u_M, u_M, rho_E_M, gamma);
    scalar p_P = calculate_pressure(rho_u_P, u_P, rho_E_P, gamma);

    // calculate lambda
    scalar lambda = calculate_lambda(rho_M, rho_P, u_M, u_P, p_M, p_P, gamma);

    // flux density
    scalar flux_density = calculate_flux(rho_u_M, rho_u_P, rho_M, rho_P, lambda, normal);

    // flux momentum
    tensor momentum_flux_M = outer_product(rho_u_M, u_M);
    for(unsigned int d = 0; d < dim; ++d)
      momentum_flux_M[d][d] += p_M;

    tensor momentum_flux_P = outer_product(rho_u_P, u_P);
    for(unsigned int d = 0; d < dim; ++d)
      momentum_flux_P[d][d] += p_P;

    vector flux_momentum =
      calculate_flux(momentum_flux_M, momentum_flux_P, rho_u_M, rho_u_P, lambda, normal);

    // flux energy
    vector energy_flux_M = (rho_E_M + p_M) * u_M;
    vector energy_flux_P = (rho_E_P + p_P) * u_P;

    scalar flux_energy =
      calculate_flux(energy_flux_M, energy_flux_P, rho_E_M, rho_E_P, lambda, normal);

    return std::make_tuple(flux_density, flux_momentum, flux_energy);
  }

  inline DEAL_II_ALWAYS_INLINE //
    std::tuple<scalar, vector, scalar>
    get_flux_boundary(FaceIntegratorScalar &         density,
                      FaceIntegratorVector &         momentum,
                      FaceIntegratorScalar &         energy,
                      BoundaryType const &           boundary_type_density,
                      BoundaryType const &           boundary_type_velocity,
                      BoundaryType const &           boundary_type_pressure,
                      BoundaryType const &           boundary_type_energy,
                      EnergyBoundaryVariable const & boundary_variable,
                      types::boundary_id const &     boundary_id,
                      unsigned int const             q) const
  {
    vector normal = momentum.get_normal_vector(q);

    // element e???
    scalar rho_M   = density.get_value(q);
    vector rho_u_M = momentum.get_value(q);
    vector u_M     = rho_u_M / rho_M;
    scalar rho_E_M = energy.get_value(q);
    scalar E_M     = rho_E_M / rho_M;
    scalar p_M     = calculate_pressure(rho_M, u_M, E_M, gamma);

    // element e???

    // calculate rho_P
    scalar rho_P = calculate_exterior_value<dim, Number>(rho_M,
                                                         boundary_type_density,
                                                         data.bc_rho,
                                                         boundary_id,
                                                         density.quadrature_point(q),
                                                         this->eval_time);

    // calculate u_P
    vector u_P = calculate_exterior_value<dim, Number>(u_M,
                                                       boundary_type_velocity,
                                                       data.bc_u,
                                                       boundary_id,
                                                       momentum.quadrature_point(q),
                                                       this->eval_time);

    vector rho_u_P = rho_P * u_P;

    // calculate p_P
    scalar p_P = calculate_exterior_value<dim, Number>(p_M,
                                                       boundary_type_pressure,
                                                       data.bc_p,
                                                       boundary_id,
                                                       density.quadrature_point(q),
                                                       this->eval_time);

    // calculate E_P
    scalar E_P = make_vectorized_array<Number>(0.0);
    if(boundary_variable == EnergyBoundaryVariable::Energy)
    {
      E_P = calculate_exterior_value<dim, Number>(E_M,
                                                  boundary_type_energy,
                                                  data.bc_E,
                                                  boundary_id,
                                                  energy.quadrature_point(q),
                                                  this->eval_time);
    }
    else if(boundary_variable == EnergyBoundaryVariable::Temperature)
    {
      scalar T_M = calculate_temperature(p_M, rho_M, R);
      scalar T_P = calculate_exterior_value<dim, Number>(T_M,
                                                         boundary_type_energy,
                                                         data.bc_E,
                                                         boundary_id,
                                                         energy.quadrature_point(q),
                                                         this->eval_time);

      E_P = calculate_energy(T_P, u_P, c_v);
    }
    scalar rho_E_P = rho_P * E_P;

    // calculate lambda
    scalar lambda = calculate_lambda(rho_M, rho_P, u_M, u_P, p_M, p_P, gamma);

    // flux density
    scalar flux_density = calculate_flux(rho_u_M, rho_u_P, rho_M, rho_P, lambda, normal);

    // flux momentum
    tensor momentum_flux_M = outer_product(rho_u_M, u_M);
    for(unsigned int d = 0; d < dim; ++d)
      momentum_flux_M[d][d] += p_M;

    tensor momentum_flux_P = outer_product(rho_u_P, u_P);
    for(unsigned int d = 0; d < dim; ++d)
      momentum_flux_P[d][d] += p_P;

    vector flux_momentum =
      calculate_flux(momentum_flux_M, momentum_flux_P, rho_u_M, rho_u_P, lambda, normal);

    // flux energy
    vector energy_flux_M = (rho_E_M + p_M) * u_M;
    vector energy_flux_P = (rho_E_P + p_P) * u_P;
    scalar flux_energy =
      calculate_flux(energy_flux_M, energy_flux_P, rho_E_M, rho_E_P, lambda, normal);

    return std::make_tuple(flux_density, flux_momentum, flux_energy);
  }

private:
  void
  cell_loop(MatrixFree<dim, Number> const &               matrix_free,
            VectorType &                                  dst,
            VectorType const &                            src,
            std::pair<unsigned int, unsigned int> const & cell_range) const
  {
    CellIntegratorScalar density(matrix_free, data.dof_index, data.quad_index, 0);
    CellIntegratorVector momentum(matrix_free, data.dof_index, data.quad_index, 1);
    CellIntegratorScalar energy(matrix_free, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      density.reinit(cell);
      density.gather_evaluate(src, true, false);

      momentum.reinit(cell);
      momentum.gather_evaluate(src, true, false);

      energy.reinit(cell);
      energy.gather_evaluate(src, true, false);

      for(unsigned int q = 0; q < momentum.n_q_points; ++q)
      {
        std::tuple<vector, tensor, vector> flux = get_volume_flux(density, momentum, energy, q);

        density.submit_gradient(-std::get<0>(flux), q);
        momentum.submit_gradient(-std::get<1>(flux), q);
        energy.submit_gradient(-std::get<2>(flux), q);
      }

      density.integrate_scatter(false, true, dst);
      momentum.integrate_scatter(false, true, dst);
      energy.integrate_scatter(false, true, dst);
    }
  }

  void
  face_loop(MatrixFree<dim, Number> const &               matrix_free,
            VectorType &                                  dst,
            VectorType const &                            src,
            std::pair<unsigned int, unsigned int> const & face_range) const
  {
    FaceIntegratorScalar density_m(matrix_free, true, data.dof_index, data.quad_index, 0);
    FaceIntegratorScalar density_p(matrix_free, false, data.dof_index, data.quad_index, 0);

    FaceIntegratorVector momentum_m(matrix_free, true, data.dof_index, data.quad_index, 1);
    FaceIntegratorVector momentum_p(matrix_free, false, data.dof_index, data.quad_index, 1);

    FaceIntegratorScalar energy_m(matrix_free, true, data.dof_index, data.quad_index, 1 + dim);
    FaceIntegratorScalar energy_p(matrix_free, false, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int face = face_range.first; face < face_range.second; face++)
    {
      // density
      density_m.reinit(face);
      density_m.gather_evaluate(src, true, false);

      density_p.reinit(face);
      density_p.gather_evaluate(src, true, false);

      // velocity
      momentum_m.reinit(face);
      momentum_m.gather_evaluate(src, true, false);

      momentum_p.reinit(face);
      momentum_p.gather_evaluate(src, true, false);

      // energy
      energy_m.reinit(face);
      energy_m.gather_evaluate(src, true, false);

      energy_p.reinit(face);
      energy_p.gather_evaluate(src, true, false);

      for(unsigned int q = 0; q < momentum_m.n_q_points; ++q)
      {
        std::tuple<scalar, vector, scalar> flux =
          get_flux(density_m, density_p, momentum_m, momentum_p, energy_m, energy_p, q);

        density_m.submit_value(std::get<0>(flux), q);
        // - sign since n??? = -n???
        density_p.submit_value(-std::get<0>(flux), q);

        momentum_m.submit_value(std::get<1>(flux), q);
        // - sign since n??? = -n???
        momentum_p.submit_value(-std::get<1>(flux), q);

        energy_m.submit_value(std::get<2>(flux), q);
        // - sign since n??? = -n???
        energy_p.submit_value(-std::get<2>(flux), q);
      }

      density_m.integrate_scatter(true, false, dst);
      density_p.integrate_scatter(true, false, dst);

      momentum_m.integrate_scatter(true, false, dst);
      momentum_p.integrate_scatter(true, false, dst);

      energy_m.integrate_scatter(true, false, dst);
      energy_p.integrate_scatter(true, false, dst);
    }
  }

  void
  boundary_face_loop(MatrixFree<dim, Number> const &               matrix_free,
                     VectorType &                                  dst,
                     VectorType const &                            src,
                     std::pair<unsigned int, unsigned int> const & face_range) const
  {
    FaceIntegratorScalar density(matrix_free, true, data.dof_index, data.quad_index, 0);
    FaceIntegratorVector momentum(matrix_free, true, data.dof_index, data.quad_index, 1);
    FaceIntegratorScalar energy(matrix_free, true, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int face = face_range.first; face < face_range.second; face++)
    {
      density.reinit(face);
      density.gather_evaluate(src, true, false);

      momentum.reinit(face);
      momentum.gather_evaluate(src, true, false);

      energy.reinit(face);
      energy.gather_evaluate(src, true, false);

      types::boundary_id boundary_id = matrix_free.get_boundary_id(face);

      BoundaryType boundary_type_density  = data.bc_rho->get_boundary_type(boundary_id);
      BoundaryType boundary_type_velocity = data.bc_u->get_boundary_type(boundary_id);
      BoundaryType boundary_type_pressure = data.bc_p->get_boundary_type(boundary_id);
      BoundaryType boundary_type_energy   = data.bc_E->get_boundary_type(boundary_id);

      EnergyBoundaryVariable boundary_variable = data.bc_E->get_boundary_variable(boundary_id);

      for(unsigned int q = 0; q < density.n_q_points; ++q)
      {
        std::tuple<scalar, vector, scalar> flux = get_flux_boundary(density,
                                                                    momentum,
                                                                    energy,
                                                                    boundary_type_density,
                                                                    boundary_type_velocity,
                                                                    boundary_type_pressure,
                                                                    boundary_type_energy,
                                                                    boundary_variable,
                                                                    boundary_id,
                                                                    q);

        density.submit_value(std::get<0>(flux), q);
        momentum.submit_value(std::get<1>(flux), q);
        energy.submit_value(std::get<2>(flux), q);
      }

      density.integrate_scatter(true, false, dst);
      momentum.integrate_scatter(true, false, dst);
      energy.integrate_scatter(true, false, dst);
    }
  }

  MatrixFree<dim, Number> const * matrix_free;

  ConvectiveOperatorData<dim> data;

  // heat capacity ratio
  Number gamma;

  // specific gas constant
  Number R;

  // specific heat at constant volume
  Number c_v;

  mutable Number eval_time;
};


template<int dim>
struct ViscousOperatorData
{
  ViscousOperatorData()
    : dof_index(0),
      quad_index(0),
      degree(1),
      IP_factor(1.0),
      dynamic_viscosity(1.0),
      reference_density(1.0),
      thermal_conductivity(0.0262),
      heat_capacity_ratio(1.4),
      specific_gas_constant(287.058)
  {
  }

  unsigned int dof_index;
  unsigned int quad_index;

  unsigned int degree;
  double       IP_factor;

  std::shared_ptr<CompNS::BoundaryDescriptor<dim>>       bc_rho;
  std::shared_ptr<CompNS::BoundaryDescriptor<dim>>       bc_u;
  std::shared_ptr<CompNS::BoundaryDescriptorEnergy<dim>> bc_E;

  double dynamic_viscosity;
  double reference_density;
  double thermal_conductivity;
  double heat_capacity_ratio;
  double specific_gas_constant;
};

template<int dim, typename Number>
class ViscousOperator
{
public:
  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

  typedef ViscousOperator<dim, Number> This;

  typedef CellIntegrator<dim, 1, Number>   CellIntegratorScalar;
  typedef FaceIntegrator<dim, 1, Number>   FaceIntegratorScalar;
  typedef CellIntegrator<dim, dim, Number> CellIntegratorVector;
  typedef FaceIntegrator<dim, dim, Number> FaceIntegratorVector;

  typedef VectorizedArray<Number>                 scalar;
  typedef Tensor<1, dim, VectorizedArray<Number>> vector;
  typedef Tensor<2, dim, VectorizedArray<Number>> tensor;
  typedef Point<dim, VectorizedArray<Number>>     point;

  ViscousOperator() : matrix_free(nullptr)
  {
  }

  void
  initialize(Mapping<dim> const &             mapping,
             MatrixFree<dim, Number> const &  matrix_free_in,
             ViscousOperatorData<dim> const & data_in)
  {
    this->matrix_free = &matrix_free_in;
    this->data        = data_in;

    gamma  = data.heat_capacity_ratio;
    R      = data.specific_gas_constant;
    c_v    = R / (gamma - 1.0);
    mu     = data.dynamic_viscosity;
    nu     = mu / data.reference_density;
    lambda = data.thermal_conductivity;

    IP::calculate_penalty_parameter<dim, Number>(
      array_penalty_parameter, *matrix_free, mapping, data.degree, data.dof_index);
  }

  void
  evaluate(VectorType & dst, VectorType const & src, Number const evaluation_time) const
  {
    dst = 0;
    evaluate_add(dst, src, evaluation_time);
  }

  void
  evaluate_add(VectorType & dst, VectorType const & src, Number const evaluation_time) const
  {
    this->eval_time = evaluation_time;

    matrix_free->loop(
      &This::cell_loop, &This::face_loop, &This::boundary_face_loop, this, dst, src);
  }

  void
  set_evaluation_time(double const & evaluation_time) const
  {
    eval_time = evaluation_time;
  }

  inline DEAL_II_ALWAYS_INLINE //
    scalar
    get_penalty_parameter(FaceIntegratorScalar & fe_eval_m, FaceIntegratorScalar & fe_eval_p) const
  {
    scalar tau = std::max(fe_eval_m.read_cell_data(array_penalty_parameter),
                          fe_eval_p.read_cell_data(array_penalty_parameter)) *
                 IP::get_penalty_factor<Number>(data.degree, data.IP_factor) * nu;

    return tau;
  }

  inline DEAL_II_ALWAYS_INLINE //
    scalar
    get_penalty_parameter(FaceIntegratorScalar & fe_eval) const
  {
    scalar tau = fe_eval.read_cell_data(array_penalty_parameter) *
                 IP::get_penalty_factor<Number>(data.degree, data.IP_factor) * nu;

    return tau;
  }

  inline DEAL_II_ALWAYS_INLINE //
    std::tuple<vector, tensor, vector>
    get_volume_flux(CellIntegratorScalar & density,
                    CellIntegratorVector & momentum,
                    CellIntegratorScalar & energy,
                    unsigned int const     q) const
  {
    scalar rho_inv  = 1.0 / density.get_value(q);
    vector grad_rho = density.get_gradient(q);

    vector rho_u      = momentum.get_value(q);
    vector u          = rho_inv * rho_u;
    tensor grad_rho_u = momentum.get_gradient(q);

    scalar rho_E      = energy.get_value(q);
    vector grad_rho_E = energy.get_gradient(q);

    // calculate flux momentum
    tensor grad_u = calculate_grad_u(rho_inv, rho_u, grad_rho, grad_rho_u);
    tensor tau    = calculate_stress_tensor(grad_u, mu);

    // calculate flux energy
    vector grad_E      = calculate_grad_E(rho_inv, rho_E, grad_rho, grad_rho_E);
    vector grad_T      = calculate_grad_T(grad_E, u, grad_u, gamma, R);
    vector energy_flux = tau * u + lambda * grad_T;

    return std::make_tuple(vector() /* dummy */, tau, energy_flux);
  }

  inline DEAL_II_ALWAYS_INLINE //
    std::tuple<scalar, vector, scalar>
    get_gradient_flux(FaceIntegratorScalar & density_m,
                      FaceIntegratorScalar & density_p,
                      FaceIntegratorVector & momentum_m,
                      FaceIntegratorVector & momentum_p,
                      FaceIntegratorScalar & energy_m,
                      FaceIntegratorScalar & energy_p,
                      scalar const &         tau_IP,
                      unsigned int const     q) const
  {
    vector normal = momentum_m.get_normal_vector(q);

    // density
    scalar rho_M      = density_m.get_value(q);
    vector grad_rho_M = density_m.get_gradient(q);

    scalar rho_P      = density_p.get_value(q);
    vector grad_rho_P = density_p.get_gradient(q);

    // velocity
    vector rho_u_M      = momentum_m.get_value(q);
    tensor grad_rho_u_M = momentum_m.get_gradient(q);

    vector rho_u_P      = momentum_p.get_value(q);
    tensor grad_rho_u_P = momentum_p.get_gradient(q);

    // energy
    scalar rho_E_M      = energy_m.get_value(q);
    vector grad_rho_E_M = energy_m.get_gradient(q);

    scalar rho_E_P      = energy_p.get_value(q);
    vector grad_rho_E_P = energy_p.get_gradient(q);

    // flux density
    scalar jump_density          = rho_M - rho_P;
    scalar gradient_flux_density = -tau_IP * jump_density;

    // flux momentum
    scalar rho_inv_M = 1.0 / rho_M;
    tensor grad_u_M  = calculate_grad_u(rho_inv_M, rho_u_M, grad_rho_M, grad_rho_u_M);
    tensor tau_M     = calculate_stress_tensor(grad_u_M, mu);

    scalar rho_inv_P = 1.0 / rho_P;
    tensor grad_u_P  = calculate_grad_u(rho_inv_P, rho_u_P, grad_rho_P, grad_rho_u_P);
    tensor tau_P     = calculate_stress_tensor(grad_u_P, mu);

    vector jump_momentum          = rho_u_M - rho_u_P;
    vector gradient_flux_momentum = 0.5 * (tau_M + tau_P) * normal - tau_IP * jump_momentum;

    // flux energy
    vector u_M      = rho_inv_M * rho_u_M;
    vector grad_E_M = calculate_grad_E(rho_inv_M, rho_E_M, grad_rho_M, grad_rho_E_M);
    vector grad_T_M = calculate_grad_T(grad_E_M, u_M, grad_u_M, gamma, R);

    vector u_P      = rho_inv_P * rho_u_P;
    vector grad_E_P = calculate_grad_E(rho_inv_P, rho_E_P, grad_rho_P, grad_rho_E_P);
    vector grad_T_P = calculate_grad_T(grad_E_P, u_P, grad_u_P, gamma, R);

    vector flux_energy_average = 0.5 * (tau_M * u_M + tau_P * u_P + lambda * (grad_T_M + grad_T_P));

    scalar jump_energy          = rho_E_M - rho_E_P;
    scalar gradient_flux_energy = flux_energy_average * normal - tau_IP * jump_energy;

    return std::make_tuple(gradient_flux_density, gradient_flux_momentum, gradient_flux_energy);
  }


  inline DEAL_II_ALWAYS_INLINE //
    std::tuple<scalar, vector, scalar>
    get_gradient_flux_boundary(FaceIntegratorScalar &         density,
                               FaceIntegratorVector &         momentum,
                               FaceIntegratorScalar &         energy,
                               scalar const &                 tau_IP,
                               BoundaryType const &           boundary_type_density,
                               BoundaryType const &           boundary_type_velocity,
                               BoundaryType const &           boundary_type_energy,
                               EnergyBoundaryVariable const & boundary_variable,
                               types::boundary_id const &     boundary_id,
                               unsigned int const             q) const
  {
    vector normal = momentum.get_normal_vector(q);

    // density
    scalar rho_M      = density.get_value(q);
    vector grad_rho_M = density.get_gradient(q);

    scalar rho_P = calculate_exterior_value<dim, Number>(rho_M,
                                                         boundary_type_density,
                                                         data.bc_rho,
                                                         boundary_id,
                                                         density.quadrature_point(q),
                                                         this->eval_time);

    scalar jump_density          = rho_M - rho_P;
    scalar gradient_flux_density = -tau_IP * jump_density;

    // velocity
    vector rho_u_M      = momentum.get_value(q);
    tensor grad_rho_u_M = momentum.get_gradient(q);

    scalar rho_inv_M = 1.0 / rho_M;
    vector u_M       = rho_inv_M * rho_u_M;

    vector u_P = calculate_exterior_value<dim, Number>(u_M,
                                                       boundary_type_velocity,
                                                       data.bc_u,
                                                       boundary_id,
                                                       momentum.quadrature_point(q),
                                                       this->eval_time);

    vector rho_u_P = rho_P * u_P;

    tensor grad_u_M = calculate_grad_u(rho_inv_M, rho_u_M, grad_rho_M, grad_rho_u_M);
    tensor tau_M    = calculate_stress_tensor(grad_u_M, mu);

    vector tau_P_normal = calculate_exterior_normal_grad(tau_M * normal,
                                                         boundary_type_velocity,
                                                         data.bc_u,
                                                         boundary_id,
                                                         momentum.quadrature_point(q),
                                                         this->eval_time);

    vector jump_momentum          = rho_u_M - rho_u_P;
    vector gradient_flux_momentum = 0.5 * (tau_M * normal + tau_P_normal) - tau_IP * jump_momentum;

    // energy
    scalar rho_E_M      = energy.get_value(q);
    vector grad_rho_E_M = energy.get_gradient(q);

    scalar E_M = rho_inv_M * rho_E_M;
    scalar E_P = make_vectorized_array<Number>(0.0);
    if(boundary_variable == EnergyBoundaryVariable::Energy)
    {
      E_P = calculate_exterior_value<dim, Number>(E_M,
                                                  boundary_type_energy,
                                                  data.bc_E,
                                                  boundary_id,
                                                  energy.quadrature_point(q),
                                                  this->eval_time);
    }
    else if(boundary_variable == EnergyBoundaryVariable::Temperature)
    {
      scalar p_M = calculate_pressure(rho_M, u_M, E_M, gamma);
      scalar T_M = calculate_temperature(p_M, rho_M, R);
      scalar T_P = calculate_exterior_value<dim, Number>(T_M,
                                                         boundary_type_energy,
                                                         data.bc_E,
                                                         boundary_id,
                                                         energy.quadrature_point(q),
                                                         this->eval_time);

      E_P = calculate_energy(T_P, u_P, c_v);
    }

    scalar rho_E_P = rho_P * E_P;

    vector grad_E_M = calculate_grad_E(rho_inv_M, rho_E_M, grad_rho_M, grad_rho_E_M);
    vector grad_T_M = calculate_grad_T(grad_E_M, u_M, grad_u_M, gamma, R);

    scalar grad_T_M_normal = grad_T_M * normal;
    scalar grad_T_P_normal = calculate_exterior_normal_grad<dim, Number>(grad_T_M_normal,
                                                                         boundary_type_energy,
                                                                         data.bc_E,
                                                                         boundary_id,
                                                                         energy.quadrature_point(q),
                                                                         this->eval_time);

    scalar jump_energy          = rho_E_M - rho_E_P;
    scalar gradient_flux_energy = 0.5 * (u_M * tau_M * normal + u_P * tau_P_normal +
                                         lambda * (grad_T_M * normal + grad_T_P_normal)) -
                                  tau_IP * jump_energy;

    return std::make_tuple(gradient_flux_density, gradient_flux_momentum, gradient_flux_energy);
  }

  inline DEAL_II_ALWAYS_INLINE //
    std::tuple<vector /*dummy_M*/,
               tensor /*value_flux_momentum_M*/,
               vector /*value_flux_energy_M*/,
               vector /*dummy_P*/,
               tensor /*value_flux_momentum_P*/,
               vector /*value_flux_energy_P*/>
    get_value_flux(FaceIntegratorScalar & density_m,
                   FaceIntegratorScalar & density_p,
                   FaceIntegratorVector & momentum_m,
                   FaceIntegratorVector & momentum_p,
                   FaceIntegratorScalar & energy_m,
                   FaceIntegratorScalar & energy_p,
                   unsigned int const     q) const
  {
    vector normal = momentum_m.get_normal_vector(q);

    // density
    scalar rho_M = density_m.get_value(q);
    scalar rho_P = density_p.get_value(q);

    // velocity
    vector rho_u_M = momentum_m.get_value(q);
    vector rho_u_P = momentum_p.get_value(q);

    // energy
    scalar rho_E_M = energy_m.get_value(q);
    scalar rho_E_P = energy_p.get_value(q);

    vector jump_rho   = (rho_M - rho_P) * normal;
    tensor jump_rho_u = outer_product(rho_u_M - rho_u_P, normal);
    vector jump_rho_E = (rho_E_M - rho_E_P) * normal;

    scalar rho_inv_M = 1.0 / rho_M;
    scalar rho_inv_P = 1.0 / rho_P;

    vector u_M = rho_inv_M * rho_u_M;
    vector u_P = rho_inv_P * rho_u_P;

    // value flux momentum
    tensor grad_u_using_jumps_M = calculate_grad_u(rho_inv_M,
                                                   rho_u_M,
                                                   jump_rho /*instead of grad_rho*/,
                                                   jump_rho_u /*instead of grad_rho_u*/);

    tensor tau_using_jumps_M     = calculate_stress_tensor(grad_u_using_jumps_M, mu);
    tensor value_flux_momentum_M = -0.5 * tau_using_jumps_M;

    tensor grad_u_using_jumps_P = calculate_grad_u(rho_inv_P,
                                                   rho_u_P,
                                                   jump_rho /*instead of grad_rho*/,
                                                   jump_rho_u /*instead of grad_rho_u*/);

    tensor tau_using_jumps_P     = calculate_stress_tensor(grad_u_using_jumps_P, mu);
    tensor value_flux_momentum_P = -0.5 * tau_using_jumps_P;

    // value flux energy
    vector grad_E_using_jumps_M = calculate_grad_E(rho_inv_M,
                                                   rho_E_M,
                                                   jump_rho /*instead of grad_rho*/,
                                                   jump_rho_E /*instead of grad_rho_E*/);

    vector grad_T_using_jumps_M =
      calculate_grad_T(grad_E_using_jumps_M, u_M, grad_u_using_jumps_M, gamma, R);
    vector value_flux_energy_M = -0.5 * (tau_using_jumps_M * u_M + lambda * grad_T_using_jumps_M);

    vector grad_E_using_jumps_P = calculate_grad_E(rho_inv_P,
                                                   rho_E_P,
                                                   jump_rho /*instead of grad_rho*/,
                                                   jump_rho_E /*instead of grad_rho_E*/);

    vector grad_T_using_jumps_P =
      calculate_grad_T(grad_E_using_jumps_P, u_P, grad_u_using_jumps_P, gamma, R);
    vector value_flux_energy_P = -0.5 * (tau_using_jumps_P * u_P + lambda * grad_T_using_jumps_P);

    return std::make_tuple(vector() /*dummy*/,
                           value_flux_momentum_M,
                           value_flux_energy_M,
                           vector() /*dummy*/,
                           value_flux_momentum_P,
                           value_flux_energy_P);
  }

  inline DEAL_II_ALWAYS_INLINE //
    std::tuple<vector /*dummy_M*/, tensor /*value_flux_momentum_M*/, vector /*value_flux_energy_M*/>
    get_value_flux_boundary(FaceIntegratorScalar &         density,
                            FaceIntegratorVector &         momentum,
                            FaceIntegratorScalar &         energy,
                            BoundaryType const &           boundary_type_density,
                            BoundaryType const &           boundary_type_velocity,
                            BoundaryType const &           boundary_type_energy,
                            EnergyBoundaryVariable const & boundary_variable,
                            types::boundary_id const &     boundary_id,
                            unsigned int const             q) const
  {
    vector normal = momentum.get_normal_vector(q);

    // density
    scalar rho_M = density.get_value(q);
    scalar rho_P = calculate_exterior_value<dim, Number>(rho_M,
                                                         boundary_type_density,
                                                         data.bc_rho,
                                                         boundary_id,
                                                         density.quadrature_point(q),
                                                         this->eval_time);

    scalar rho_inv_M = 1.0 / rho_M;

    // velocity
    vector rho_u_M = momentum.get_value(q);
    vector u_M     = rho_inv_M * rho_u_M;

    vector u_P = calculate_exterior_value<dim, Number>(u_M,
                                                       boundary_type_velocity,
                                                       data.bc_u,
                                                       boundary_id,
                                                       momentum.quadrature_point(q),
                                                       this->eval_time);

    vector rho_u_P = rho_P * u_P;

    // energy
    scalar rho_E_M = energy.get_value(q);
    scalar E_M     = rho_inv_M * rho_E_M;

    scalar E_P = make_vectorized_array<Number>(0.0);
    if(boundary_variable == EnergyBoundaryVariable::Energy)
    {
      E_P = calculate_exterior_value<dim, Number>(E_M,
                                                  boundary_type_energy,
                                                  data.bc_E,
                                                  boundary_id,
                                                  energy.quadrature_point(q),
                                                  this->eval_time);
    }
    else if(boundary_variable == EnergyBoundaryVariable::Temperature)
    {
      scalar p_M = calculate_pressure(rho_M, u_M, E_M, gamma);
      scalar T_M = calculate_temperature(p_M, rho_M, R);
      scalar T_P = calculate_exterior_value<dim, Number>(T_M,
                                                         boundary_type_energy,
                                                         data.bc_E,
                                                         boundary_id,
                                                         energy.quadrature_point(q),
                                                         this->eval_time);

      E_P = calculate_energy(T_P, u_P, c_v);
    }
    scalar rho_E_P = rho_P * E_P;

    vector jump_rho   = (rho_M - rho_P) * normal;
    tensor jump_rho_u = outer_product(rho_u_M - rho_u_P, normal);
    vector jump_rho_E = (rho_E_M - rho_E_P) * normal;

    // value flux momentum
    tensor grad_u_using_jumps_M = calculate_grad_u(rho_inv_M,
                                                   rho_u_M,
                                                   jump_rho /*instead of grad_rho*/,
                                                   jump_rho_u /*instead of grad_rho_u*/);

    tensor tau_using_jumps_M     = calculate_stress_tensor(grad_u_using_jumps_M, mu);
    tensor value_flux_momentum_M = -0.5 * tau_using_jumps_M;

    // value flux energy
    vector grad_E_using_jumps_M = calculate_grad_E(rho_inv_M,
                                                   rho_E_M,
                                                   jump_rho /*instead of grad_rho*/,
                                                   jump_rho_E /*instead of grad_rho_E*/);

    vector grad_T_using_jumps_M =
      calculate_grad_T(grad_E_using_jumps_M, u_M, grad_u_using_jumps_M, gamma, R);
    vector value_flux_energy_M = -0.5 * (tau_using_jumps_M * u_M + lambda * grad_T_using_jumps_M);

    return std::make_tuple(vector() /*dummy*/, value_flux_momentum_M, value_flux_energy_M);
  }

private:
  void
  cell_loop(MatrixFree<dim, Number> const &               matrix_free,
            VectorType &                                  dst,
            VectorType const &                            src,
            std::pair<unsigned int, unsigned int> const & cell_range) const
  {
    CellIntegratorScalar density(matrix_free, data.dof_index, data.quad_index, 0);
    CellIntegratorVector momentum(matrix_free, data.dof_index, data.quad_index, 1);
    CellIntegratorScalar energy(matrix_free, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      density.reinit(cell);
      density.gather_evaluate(src, true, true);

      momentum.reinit(cell);
      momentum.gather_evaluate(src, true, true);

      energy.reinit(cell);
      energy.gather_evaluate(src, true, true);

      for(unsigned int q = 0; q < momentum.n_q_points; ++q)
      {
        std::tuple<vector, tensor, vector> flux = get_volume_flux(density, momentum, energy, q);

        momentum.submit_gradient(std::get<1>(flux), q);
        energy.submit_gradient(std::get<2>(flux), q);
      }

      momentum.integrate_scatter(false, true, dst);
      energy.integrate_scatter(false, true, dst);
    }
  }

  void
  face_loop(MatrixFree<dim, Number> const &               matrix_free,
            VectorType &                                  dst,
            VectorType const &                            src,
            std::pair<unsigned int, unsigned int> const & face_range) const
  {
    FaceIntegratorScalar density_m(matrix_free, true, data.dof_index, data.quad_index, 0);
    FaceIntegratorScalar density_p(matrix_free, false, data.dof_index, data.quad_index, 0);

    FaceIntegratorVector momentum_m(matrix_free, true, data.dof_index, data.quad_index, 1);
    FaceIntegratorVector momentum_p(matrix_free, false, data.dof_index, data.quad_index, 1);

    FaceIntegratorScalar energy_m(matrix_free, true, data.dof_index, data.quad_index, 1 + dim);
    FaceIntegratorScalar energy_p(matrix_free, false, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int face = face_range.first; face < face_range.second; face++)
    {
      // density
      density_m.reinit(face);
      density_m.gather_evaluate(src, true, true);

      density_p.reinit(face);
      density_p.gather_evaluate(src, true, true);

      // momentum
      momentum_m.reinit(face);
      momentum_m.gather_evaluate(src, true, true);

      momentum_p.reinit(face);
      momentum_p.gather_evaluate(src, true, true);

      // energy
      energy_m.reinit(face);
      energy_m.gather_evaluate(src, true, true);

      energy_p.reinit(face);
      energy_p.gather_evaluate(src, true, true);

      scalar tau_IP = get_penalty_parameter(density_m, density_p);

      for(unsigned int q = 0; q < density_m.n_q_points; ++q)
      {
        std::tuple<scalar, vector, scalar> gradient_flux = get_gradient_flux(
          density_m, density_p, momentum_m, momentum_p, energy_m, energy_p, tau_IP, q);

        std::tuple<vector, tensor, vector, vector, tensor, vector> value_flux =
          get_value_flux(density_m, density_p, momentum_m, momentum_p, energy_m, energy_p, q);

        density_m.submit_value(-std::get<0>(gradient_flux), q);
        // + sign since n??? = -n???
        density_p.submit_value(std::get<0>(gradient_flux), q);

        momentum_m.submit_gradient(std::get<1>(value_flux), q);
        // note that value_flux_momentum is not conservative
        momentum_p.submit_gradient(std::get<4>(value_flux), q);

        momentum_m.submit_value(-std::get<1>(gradient_flux), q);
        // + sign since n??? = -n???
        momentum_p.submit_value(std::get<1>(gradient_flux), q);

        energy_m.submit_gradient(std::get<2>(value_flux), q);
        // note that value_flux_energy is not conservative
        energy_p.submit_gradient(std::get<5>(value_flux), q);

        energy_m.submit_value(-std::get<2>(gradient_flux), q);
        // + sign since n??? = -n???
        energy_p.submit_value(std::get<2>(gradient_flux), q);
      }

      density_m.integrate_scatter(true, false, dst);
      density_p.integrate_scatter(true, false, dst);

      momentum_m.integrate_scatter(true, true, dst);
      momentum_p.integrate_scatter(true, true, dst);

      energy_m.integrate_scatter(true, true, dst);
      energy_p.integrate_scatter(true, true, dst);
    }
  }

  void
  boundary_face_loop(MatrixFree<dim, Number> const &               matrix_free,
                     VectorType &                                  dst,
                     VectorType const &                            src,
                     std::pair<unsigned int, unsigned int> const & face_range) const
  {
    FaceIntegratorScalar density(matrix_free, true, data.dof_index, data.quad_index, 0);
    FaceIntegratorVector momentum(matrix_free, true, data.dof_index, data.quad_index, 1);
    FaceIntegratorScalar energy(matrix_free, true, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int face = face_range.first; face < face_range.second; face++)
    {
      density.reinit(face);
      density.gather_evaluate(src, true, true);

      momentum.reinit(face);
      momentum.gather_evaluate(src, true, true);

      energy.reinit(face);
      energy.gather_evaluate(src, true, true);

      scalar tau_IP = get_penalty_parameter(density);

      types::boundary_id boundary_id = matrix_free.get_boundary_id(face);

      BoundaryType boundary_type_density  = data.bc_rho->get_boundary_type(boundary_id);
      BoundaryType boundary_type_velocity = data.bc_u->get_boundary_type(boundary_id);
      BoundaryType boundary_type_energy   = data.bc_E->get_boundary_type(boundary_id);

      EnergyBoundaryVariable boundary_variable = data.bc_E->get_boundary_variable(boundary_id);

      for(unsigned int q = 0; q < density.n_q_points; ++q)
      {
        std::tuple<scalar, vector, scalar> gradient_flux =
          get_gradient_flux_boundary(density,
                                     momentum,
                                     energy,
                                     tau_IP,
                                     boundary_type_density,
                                     boundary_type_velocity,
                                     boundary_type_energy,
                                     boundary_variable,
                                     boundary_id,
                                     q);

        std::tuple<vector, tensor, vector> value_flux =
          get_value_flux_boundary(density,
                                  momentum,
                                  energy,
                                  boundary_type_density,
                                  boundary_type_velocity,
                                  boundary_type_energy,
                                  boundary_variable,
                                  boundary_id,
                                  q);

        density.submit_value(-std::get<0>(gradient_flux), q);

        momentum.submit_gradient(std::get<1>(value_flux), q);
        momentum.submit_value(-std::get<1>(gradient_flux), q);

        energy.submit_gradient(std::get<2>(value_flux), q);
        energy.submit_value(-std::get<2>(gradient_flux), q);
      }

      density.integrate_scatter(true, false, dst);
      momentum.integrate_scatter(true, true, dst);
      energy.integrate_scatter(true, true, dst);
    }
  }

  MatrixFree<dim, Number> const * matrix_free;

  ViscousOperatorData<dim> data;

  // heat capacity ratio
  Number gamma;

  // specific gas constant
  Number R;

  // specific heat at constant volume
  Number c_v;

  // dynamic viscosity
  Number mu;

  // kinematic viscosity
  Number nu;

  // thermal conductivity
  Number lambda;

  AlignedVector<VectorizedArray<Number>> array_penalty_parameter;

  mutable Number eval_time;
};

template<int dim>
struct CombinedOperatorData
{
  CombinedOperatorData() : dof_index(0), quad_index(0)
  {
  }

  unsigned int dof_index;
  unsigned int quad_index;

  std::shared_ptr<CompNS::BoundaryDescriptor<dim>>       bc_rho;
  std::shared_ptr<CompNS::BoundaryDescriptor<dim>>       bc_u;
  std::shared_ptr<CompNS::BoundaryDescriptor<dim>>       bc_p;
  std::shared_ptr<CompNS::BoundaryDescriptorEnergy<dim>> bc_E;
};

template<int dim, typename Number>
class CombinedOperator
{
public:
  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

  typedef ConvectiveOperator<dim, Number> ConvectiveOp;
  typedef ViscousOperator<dim, Number>    ViscousOp;
  typedef CombinedOperator<dim, Number>   This;

  typedef CellIntegrator<dim, 1, Number>   CellIntegratorScalar;
  typedef FaceIntegrator<dim, 1, Number>   FaceIntegratorScalar;
  typedef CellIntegrator<dim, dim, Number> CellIntegratorVector;
  typedef FaceIntegrator<dim, dim, Number> FaceIntegratorVector;

  typedef VectorizedArray<Number>                 scalar;
  typedef Tensor<1, dim, VectorizedArray<Number>> vector;
  typedef Tensor<2, dim, VectorizedArray<Number>> tensor;
  typedef Point<dim, VectorizedArray<Number>>     point;

  CombinedOperator() : matrix_free(nullptr), convective_operator(nullptr), viscous_operator(nullptr)
  {
  }

  void
  initialize(MatrixFree<dim, Number> const &   matrix_free_in,
             CombinedOperatorData<dim> const & data_in,
             ConvectiveOp const &              convective_operator_in,
             ViscousOp const &                 viscous_operator_in)
  {
    this->matrix_free = &matrix_free_in;
    this->data        = data_in;

    this->convective_operator = &convective_operator_in;
    this->viscous_operator    = &viscous_operator_in;
  }

  void
  evaluate(VectorType & dst, VectorType const & src, Number const evaluation_time) const
  {
    dst = 0;
    evaluate_add(dst, src, evaluation_time);
  }

  void
  evaluate_add(VectorType & dst, VectorType const & src, Number const evaluation_time) const
  {
    convective_operator->set_evaluation_time(evaluation_time);
    viscous_operator->set_evaluation_time(evaluation_time);

    matrix_free->loop(
      &This::cell_loop, &This::face_loop, &This::boundary_face_loop, this, dst, src);

    // perform cell integrals only for performance measurements
    //    matrix_free->cell_loop(&This::cell_loop, this, dst, src);
  }

private:
  void
  cell_loop(MatrixFree<dim, Number> const &               matrix_free,
            VectorType &                                  dst,
            VectorType const &                            src,
            std::pair<unsigned int, unsigned int> const & cell_range) const
  {
    CellIntegratorScalar density(matrix_free, data.dof_index, data.quad_index, 0);
    CellIntegratorVector momentum(matrix_free, data.dof_index, data.quad_index, 1);
    CellIntegratorScalar energy(matrix_free, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      density.reinit(cell);
      density.gather_evaluate(src, true, true);

      momentum.reinit(cell);
      momentum.gather_evaluate(src, true, true);

      energy.reinit(cell);
      energy.gather_evaluate(src, true, true);

      for(unsigned int q = 0; q < momentum.n_q_points; ++q)
      {
        std::tuple<vector, tensor, vector> conv_flux =
          convective_operator->get_volume_flux(density, momentum, energy, q);

        std::tuple<vector, tensor, vector> visc_flux =
          viscous_operator->get_volume_flux(density, momentum, energy, q);

        density.submit_gradient(-std::get<0>(conv_flux), q);
        momentum.submit_gradient(-std::get<1>(conv_flux) + std::get<1>(visc_flux), q);
        energy.submit_gradient(-std::get<2>(conv_flux) + std::get<2>(visc_flux), q);
      }

      density.integrate_scatter(false, true, dst);
      momentum.integrate_scatter(false, true, dst);
      energy.integrate_scatter(false, true, dst);
    }
  }

  void
  face_loop(MatrixFree<dim, Number> const &               matrix_free,
            VectorType &                                  dst,
            VectorType const &                            src,
            std::pair<unsigned int, unsigned int> const & face_range) const
  {
    FaceIntegratorScalar density_m(matrix_free, true, data.dof_index, data.quad_index, 0);
    FaceIntegratorScalar density_p(matrix_free, false, data.dof_index, data.quad_index, 0);
    FaceIntegratorVector momentum_m(matrix_free, true, data.dof_index, data.quad_index, 1);
    FaceIntegratorVector momentum_p(matrix_free, false, data.dof_index, data.quad_index, 1);
    FaceIntegratorScalar energy_m(matrix_free, true, data.dof_index, data.quad_index, 1 + dim);
    FaceIntegratorScalar energy_p(matrix_free, false, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int face = face_range.first; face < face_range.second; face++)
    {
      // density
      density_m.reinit(face);
      density_m.gather_evaluate(src, true, true);

      density_p.reinit(face);
      density_p.gather_evaluate(src, true, true);

      // momentum
      momentum_m.reinit(face);
      momentum_m.gather_evaluate(src, true, true);

      momentum_p.reinit(face);
      momentum_p.gather_evaluate(src, true, true);

      // energy
      energy_m.reinit(face);
      energy_m.gather_evaluate(src, true, true);

      energy_p.reinit(face);
      energy_p.gather_evaluate(src, true, true);

      scalar tau_IP = viscous_operator->get_penalty_parameter(density_m, density_p);

      for(unsigned int q = 0; q < density_m.n_q_points; ++q)
      {
        std::tuple<scalar, vector, scalar> conv_flux = convective_operator->get_flux(
          density_m, density_p, momentum_m, momentum_p, energy_m, energy_p, q);

        std::tuple<scalar, vector, scalar> visc_grad_flux = viscous_operator->get_gradient_flux(
          density_m, density_p, momentum_m, momentum_p, energy_m, energy_p, tau_IP, q);

        std::tuple<vector, tensor, vector, vector, tensor, vector> visc_value_flux =
          viscous_operator->get_value_flux(
            density_m, density_p, momentum_m, momentum_p, energy_m, energy_p, q);

        density_m.submit_value(std::get<0>(conv_flux) - std::get<0>(visc_grad_flux), q);
        // - sign since n??? = -n???
        density_p.submit_value(-std::get<0>(conv_flux) + std::get<0>(visc_grad_flux), q);

        momentum_m.submit_value(std::get<1>(conv_flux) - std::get<1>(visc_grad_flux), q);
        // - sign since n??? = -n???
        momentum_p.submit_value(-std::get<1>(conv_flux) + std::get<1>(visc_grad_flux), q);

        momentum_m.submit_gradient(std::get<1>(visc_value_flux), q);
        // note that value_flux_momentum is not conservative
        momentum_p.submit_gradient(std::get<4>(visc_value_flux), q);

        energy_m.submit_value(std::get<2>(conv_flux) - std::get<2>(visc_grad_flux), q);
        // - sign since n??? = -n???
        energy_p.submit_value(-std::get<2>(conv_flux) + std::get<2>(visc_grad_flux), q);

        energy_m.submit_gradient(std::get<2>(visc_value_flux), q);
        // note that value_flux_energy is not conservative
        energy_p.submit_gradient(std::get<5>(visc_value_flux), q);
      }

      density_m.integrate_scatter(true, false, dst);
      density_p.integrate_scatter(true, false, dst);

      momentum_m.integrate_scatter(true, true, dst);
      momentum_p.integrate_scatter(true, true, dst);

      energy_m.integrate_scatter(true, true, dst);
      energy_p.integrate_scatter(true, true, dst);
    }
  }

  void
  boundary_face_loop(MatrixFree<dim, Number> const &               matrix_free,
                     VectorType &                                  dst,
                     VectorType const &                            src,
                     std::pair<unsigned int, unsigned int> const & face_range) const
  {
    FaceIntegratorScalar density(matrix_free, true, data.dof_index, data.quad_index, 0);
    FaceIntegratorVector momentum(matrix_free, true, data.dof_index, data.quad_index, 1);
    FaceIntegratorScalar energy(matrix_free, true, data.dof_index, data.quad_index, 1 + dim);

    for(unsigned int face = face_range.first; face < face_range.second; face++)
    {
      density.reinit(face);
      density.gather_evaluate(src, true, true);

      momentum.reinit(face);
      momentum.gather_evaluate(src, true, true);

      energy.reinit(face);
      energy.gather_evaluate(src, true, true);

      scalar tau_IP = viscous_operator->get_penalty_parameter(density);

      types::boundary_id boundary_id = matrix_free.get_boundary_id(face);

      BoundaryType boundary_type_density  = data.bc_rho->get_boundary_type(boundary_id);
      BoundaryType boundary_type_velocity = data.bc_u->get_boundary_type(boundary_id);
      BoundaryType boundary_type_pressure = data.bc_p->get_boundary_type(boundary_id);
      BoundaryType boundary_type_energy   = data.bc_E->get_boundary_type(boundary_id);

      EnergyBoundaryVariable boundary_variable = data.bc_E->get_boundary_variable(boundary_id);

      for(unsigned int q = 0; q < density.n_q_points; ++q)
      {
        std::tuple<scalar, vector, scalar> conv_flux =
          convective_operator->get_flux_boundary(density,
                                                 momentum,
                                                 energy,
                                                 boundary_type_density,
                                                 boundary_type_velocity,
                                                 boundary_type_pressure,
                                                 boundary_type_energy,
                                                 boundary_variable,
                                                 boundary_id,
                                                 q);

        std::tuple<scalar, vector, scalar> visc_grad_flux =
          viscous_operator->get_gradient_flux_boundary(density,
                                                       momentum,
                                                       energy,
                                                       tau_IP,
                                                       boundary_type_density,
                                                       boundary_type_velocity,
                                                       boundary_type_energy,
                                                       boundary_variable,
                                                       boundary_id,
                                                       q);

        std::tuple<vector, tensor, vector> visc_value_flux =
          viscous_operator->get_value_flux_boundary(density,
                                                    momentum,
                                                    energy,
                                                    boundary_type_density,
                                                    boundary_type_velocity,
                                                    boundary_type_energy,
                                                    boundary_variable,
                                                    boundary_id,
                                                    q);

        density.submit_value(std::get<0>(conv_flux) - std::get<0>(visc_grad_flux), q);

        momentum.submit_value(std::get<1>(conv_flux) - std::get<1>(visc_grad_flux), q);
        momentum.submit_gradient(std::get<1>(visc_value_flux), q);

        energy.submit_value(std::get<2>(conv_flux) - std::get<2>(visc_grad_flux), q);
        energy.submit_gradient(std::get<2>(visc_value_flux), q);
      }

      density.integrate_scatter(true, false, dst);
      momentum.integrate_scatter(true, true, dst);
      energy.integrate_scatter(true, true, dst);
    }
  }

  MatrixFree<dim, Number> const * matrix_free;

  CombinedOperatorData<dim> data;

  ConvectiveOperator<dim, Number> const * convective_operator;
  ViscousOperator<dim, Number> const *    viscous_operator;
};

} // namespace CompNS

#endif /* INCLUDE_COMPRESSIBLE_NAVIER_STOKES_SPATIAL_DISCRETIZATION_COMP_NAVIER_STOKES_OPERATORS_H_ \
        */
