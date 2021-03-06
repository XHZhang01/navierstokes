/*
 * gradient_operator.h
 *
 *  Created on: Nov 5, 2018
 *      Author: fehn
 */

#ifndef INCLUDE_INCOMPRESSIBLE_NAVIER_STOKES_SPATIAL_DISCRETIZATION_OPERATORS_GRADIENT_OPERATOR_H_
#define INCLUDE_INCOMPRESSIBLE_NAVIER_STOKES_SPATIAL_DISCRETIZATION_OPERATORS_GRADIENT_OPERATOR_H_

#include <deal.II/matrix_free/fe_evaluation_notemplate.h>

#include "../../../operators/mapping_flags.h"
#include "../../user_interface/input_parameters.h"
#include "weak_boundary_conditions.h"

using namespace dealii;

namespace IncNS
{
namespace Operators
{
template<int dim, typename Number>
class GradientKernel
{
private:
  typedef CellIntegrator<dim, 1, Number> CellIntegratorP;

  typedef VectorizedArray<Number>                 scalar;
  typedef Tensor<1, dim, VectorizedArray<Number>> vector;

public:
  static MappingFlags
  get_mapping_flags()
  {
    MappingFlags flags;

    flags.cells          = update_JxW_values | update_gradients;
    flags.inner_faces    = update_JxW_values | update_normal_vectors;
    flags.boundary_faces = update_JxW_values | update_quadrature_points | update_normal_vectors;

    return flags;
  }

  /*
   *  This function implements the central flux as numerical flux function.
   */
  inline DEAL_II_ALWAYS_INLINE //
    scalar
    calculate_flux(scalar const & value_m, scalar const & value_p) const
  {
    return 0.5 * (value_m + value_p);
  }

  /*
   * Volume flux, i.e., the term occurring in the volume integral for
   * weak formulation (performing integration-by-parts)
   */
  inline DEAL_II_ALWAYS_INLINE //
    scalar
    get_volume_flux_weak(CellIntegratorP & pressure, unsigned int const q) const
  {
    // minus sign due to integration by parts
    return -pressure.get_value(q);
  }

  /*
   * Volume flux, i.e., the term occurring in the volume integral for
   * strong formulation (no integration-by-parts, or integration-by-parts performed twice)
   */
  inline DEAL_II_ALWAYS_INLINE //
    vector
    get_volume_flux_strong(CellIntegratorP & pressure, unsigned int const q) const
  {
    return pressure.get_gradient(q);
  }
};
} // namespace Operators

template<int dim>
struct GradientOperatorData
{
  GradientOperatorData()
    : dof_index_velocity(0),
      dof_index_pressure(1),
      quad_index(0),
      integration_by_parts(true),
      use_boundary_data(true)
  {
  }

  unsigned int dof_index_velocity;
  unsigned int dof_index_pressure;

  unsigned int quad_index;

  bool integration_by_parts;
  bool use_boundary_data;

  std::shared_ptr<BoundaryDescriptorP<dim>> bc;
};

template<int dim, typename Number>
class GradientOperator
{
public:
  typedef GradientOperator<dim, Number> This;

  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

  typedef VectorizedArray<Number>                 scalar;
  typedef Tensor<1, dim, VectorizedArray<Number>> vector;

  typedef std::pair<unsigned int, unsigned int> Range;

  typedef CellIntegrator<dim, dim, Number> CellIntegratorU;
  typedef CellIntegrator<dim, 1, Number>   CellIntegratorP;

  typedef FaceIntegrator<dim, dim, Number> FaceIntegratorU;
  typedef FaceIntegrator<dim, 1, Number>   FaceIntegratorP;

  GradientOperator();

  void
  reinit(MatrixFree<dim, Number> const & matrix_free_in, GradientOperatorData<dim> const & data_in);

  void
  set_scaling_factor_pressure(double const & scaling_factor);

  GradientOperatorData<dim> const &
  get_operator_data() const;

  // homogeneous operator
  void
  apply(VectorType & dst, const VectorType & src) const;

  void
  apply_add(VectorType & dst, const VectorType & src) const;

  // inhomogeneous operator
  void
  rhs(VectorType & dst, Number const evaluation_time) const;

  void
  rhs_add(VectorType & dst, Number const evaluation_time) const;

  // full operator, i.e., homogeneous and inhomogeneous contributions
  void
  evaluate(VectorType & dst, VectorType const & src, Number const evaluation_time) const;

  void
  evaluate_add(VectorType & dst, VectorType const & src, Number const evaluation_time) const;

private:
  void
  do_cell_integral_weak(CellIntegratorP & pressure, CellIntegratorU & velocity) const;

  void
  do_cell_integral_strong(CellIntegratorP & pressure, CellIntegratorU & velocity) const;

  void
  do_face_integral(FaceIntegratorP & pressure_m,
                   FaceIntegratorP & pressure_p,
                   FaceIntegratorU & velocity_m,
                   FaceIntegratorU & velocity_p) const;

  void
  do_boundary_integral(FaceIntegratorP &          pressure,
                       FaceIntegratorU &          velocity,
                       OperatorType const &       operator_type,
                       types::boundary_id const & boundary_id) const;

  void
  cell_loop(MatrixFree<dim, Number> const & matrix_free,
            VectorType &                    dst,
            VectorType const &              src,
            Range const &                   cell_range) const;

  void
  face_loop(MatrixFree<dim, Number> const & matrix_free,
            VectorType &                    dst,
            VectorType const &              src,
            Range const &                   face_range) const;

  void
  boundary_face_loop_hom_operator(MatrixFree<dim, Number> const & matrix_free,
                                  VectorType &                    dst,
                                  VectorType const &              src,
                                  Range const &                   face_range) const;

  void
  boundary_face_loop_full_operator(MatrixFree<dim, Number> const & matrix_free,
                                   VectorType &                    dst,
                                   VectorType const &              src,
                                   Range const &                   face_range) const;

  void
  cell_loop_inhom_operator(MatrixFree<dim, Number> const & matrix_free,
                           VectorType &                    dst,
                           VectorType const &              src,
                           Range const &                   cell_range) const;

  void
  face_loop_inhom_operator(MatrixFree<dim, Number> const & matrix_free,
                           VectorType &                    dst,
                           VectorType const &              src,
                           Range const &                   face_range) const;

  void
  boundary_face_loop_inhom_operator(MatrixFree<dim, Number> const & matrix_free,
                                    VectorType &                    dst,
                                    VectorType const &              src,
                                    Range const &                   face_range) const;

  MatrixFree<dim, Number> const * matrix_free;

  GradientOperatorData<dim> data;

  mutable double time;

  // if the continuity equation of the incompressible Navier-Stokes
  // equations is scaled by a constant factor, the system of equations
  // is solved for a modified pressure p^* = 1/scaling_factor * p. Hence,
  // when applying the gradient operator to this modified pressure we have
  // to make sure that we also apply the correct boundary conditions for p^*,
  // i.e., g_p^* = 1/scaling_factor * g_p
  double inverse_scaling_factor_pressure;

  Operators::GradientKernel<dim, Number> kernel;
};

} // namespace IncNS



#endif /* INCLUDE_INCOMPRESSIBLE_NAVIER_STOKES_SPATIAL_DISCRETIZATION_OPERATORS_GRADIENT_OPERATOR_H_ \
        */
