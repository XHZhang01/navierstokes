/*
 * convection_diffusion_operator_merged.h
 *
 *  Created on: Jun 6, 2019
 *      Author: fehn
 */

#ifndef INCLUDE_CONVECTION_DIFFUSION_SPATIAL_DISCRETIZATION_OPERATORS_CONVECTION_DIFFUSION_OPERATOR_MERGED_H_
#define INCLUDE_CONVECTION_DIFFUSION_SPATIAL_DISCRETIZATION_OPERATORS_CONVECTION_DIFFUSION_OPERATOR_MERGED_H_

#include "convective_operator.h"
#include "diffusive_operator.h"
#include "mass_operator.h"

namespace ConvDiff
{
template<int dim>
struct ConvectionDiffusionOperatorMergedData : public OperatorBaseData
{
  ConvectionDiffusionOperatorMergedData()
    : OperatorBaseData(0 /* dof_index */, 0 /* quad_index */),
      unsteady_problem(false),
      convective_problem(false),
      diffusive_problem(false)
  {
  }

  bool unsteady_problem;
  bool convective_problem;
  bool diffusive_problem;

  Operators::ConvectiveKernelData<dim> convective_kernel_data;
  Operators::DiffusiveKernelData       diffusive_kernel_data;

  std::shared_ptr<ConvDiff::BoundaryDescriptor<dim>> bc;
};

template<int dim, typename Number>
class ConvectionDiffusionOperatorMerged
  : public OperatorBase<dim, Number, ConvectionDiffusionOperatorMergedData<dim>>
{
private:
  typedef OperatorBase<dim, Number, ConvectionDiffusionOperatorMergedData<dim>> Base;

  typedef typename Base::IntegratorCell IntegratorCell;
  typedef typename Base::IntegratorFace IntegratorFace;

  typedef typename Base::VectorType VectorType;

  typedef VectorizedArray<Number>                 scalar;
  typedef Tensor<1, dim, VectorizedArray<Number>> vector;

public:
  void
  reinit(MatrixFree<dim, Number> const &                    matrix_free,
         AffineConstraints<double> const &                  constraint_matrix,
         ConvectionDiffusionOperatorMergedData<dim> const & operator_data) const;

  LinearAlgebra::distributed::Vector<Number> const &
  get_velocity() const;

  void
  set_velocity_copy(VectorType const & velocity) const;

  void
  set_velocity_ptr(VectorType const & velocity) const;

  void
  set_scaling_factor_mass_matrix(Number const & number) const;

private:
  void
  reinit_cell(unsigned int const cell) const;

  void
  reinit_face(unsigned int const face) const;

  void
  reinit_boundary_face(unsigned int const face) const;

  void
  reinit_face_cell_based(unsigned int const       cell,
                         unsigned int const       face,
                         types::boundary_id const boundary_id) const;

  void
  do_cell_integral(IntegratorCell & integrator) const;

  void
  do_face_integral(IntegratorFace & integrator_m, IntegratorFace & integrator_p) const;

  void
  do_face_int_integral(IntegratorFace & integrator_m, IntegratorFace & integrator_p) const;

  void
  do_face_ext_integral(IntegratorFace & integrator_m, IntegratorFace & integrator_p) const;

  void
  do_boundary_integral(IntegratorFace &           integrator_m,
                       OperatorType const &       operator_type,
                       types::boundary_id const & boundary_id) const;

  // TODO can be removed later once matrix-free evaluation allows accessing neighboring data for
  // cell-based face loops
  void
  do_face_int_integral_cell_based(IntegratorFace & integrator_m,
                                  IntegratorFace & integrator_p) const;

  void
  do_verify_boundary_conditions(types::boundary_id const                           boundary_id,
                                ConvectionDiffusionOperatorMergedData<dim> const & operator_data,
                                std::set<types::boundary_id> const & periodic_boundary_ids) const;

  Operators::MassMatrixKernel<dim, Number> mass_kernel;
  Operators::ConvectiveKernel<dim, Number> convective_kernel;
  Operators::DiffusiveKernel<dim, Number>  diffusive_kernel;
};

} // namespace ConvDiff


#endif /* INCLUDE_CONVECTION_DIFFUSION_SPATIAL_DISCRETIZATION_OPERATORS_CONVECTION_DIFFUSION_OPERATOR_MERGED_H_ \
        */