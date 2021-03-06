#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include "mg_transfer_mf_c.h"

template<int dim, typename Number, typename VectorType, int components>
MGTransferMFC<dim, Number, VectorType, components>::MGTransferMFC(
  const MF &                        data_dg,
  const MF &                        data_cg,
  const AffineConstraints<double> & cm_dg,
  const AffineConstraints<double> & cm_cg,
  const unsigned int                level,
  const unsigned int                fe_degree,
  const unsigned int                dof_handler_index)
  : fe_degree(fe_degree)
{
  std::vector<const DoFHandler<dim> *> dofhandlers = {&data_cg.get_dof_handler(dof_handler_index),
                                                      &data_dg.get_dof_handler(dof_handler_index)};

  std::vector<const AffineConstraints<double> *> constraint_matrices = {&cm_cg, &cm_dg};
  QGauss<1>                                      quadrature(1);

  typename MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.level_mg_handler = level;
  data_composite.reinit(dofhandlers, constraint_matrices, quadrature, additional_data);
}

template<int dim, typename Number, typename VectorType, int components>
MGTransferMFC<dim, Number, VectorType, components>::~MGTransferMFC()
{
}

template<int dim, typename Number, typename VectorType, int components>
template<int degree>
void
MGTransferMFC<dim, Number, VectorType, components>::do_interpolate(VectorType &       dst,
                                                                   const VectorType & src) const
{
  FEEvaluation<dim, degree, 1, components, Number> fe_eval_cg(data_composite, 0);
  FEEvaluation<dim, degree, 1, components, Number> fe_eval_dg(data_composite, 1);

  VectorType vec__dg;
  data_composite.initialize_dof_vector(vec__dg, 1);
  vec__dg.copy_locally_owned_data_from(src);

  for(unsigned int cell = 0; cell < data_composite.n_macro_cells(); ++cell)
  {
    fe_eval_cg.reinit(cell);
    fe_eval_dg.reinit(cell);

    fe_eval_dg.read_dof_values(vec__dg);

    for(unsigned int i = 0; i < fe_eval_cg.static_dofs_per_cell; i++)
      fe_eval_cg.begin_dof_values()[i] = fe_eval_dg.begin_dof_values()[i];

    fe_eval_cg.set_dof_values(dst);
  }
}

template<int dim, typename Number, typename VectorType, int components>
template<int degree>
void
MGTransferMFC<dim, Number, VectorType, components>::do_restrict_and_add(
  VectorType &       dst,
  const VectorType & src) const
{
  FEEvaluation<dim, degree, 1, components, Number> fe_eval_cg(data_composite, 0);
  FEEvaluation<dim, degree, 1, components, Number> fe_eval_dg(data_composite, 1);

  VectorType vec__dg;
  data_composite.initialize_dof_vector(vec__dg, 1);
  vec__dg.copy_locally_owned_data_from(src);

  for(unsigned int cell = 0; cell < data_composite.n_macro_cells(); ++cell)
  {
    fe_eval_cg.reinit(cell);
    fe_eval_dg.reinit(cell);

    fe_eval_dg.read_dof_values(vec__dg);

    for(unsigned int i = 0; i < fe_eval_cg.static_dofs_per_cell; i++)
      fe_eval_cg.begin_dof_values()[i] = fe_eval_dg.begin_dof_values()[i];

    fe_eval_cg.distribute_local_to_global(dst);
  }

  dst.compress(VectorOperation::add);
}

template<int dim, typename Number, typename VectorType, int components>
template<int degree>
void
MGTransferMFC<dim, Number, VectorType, components>::do_prolongate(VectorType &       dst,
                                                                  const VectorType & src) const
{
  src.update_ghost_values();

  FEEvaluation<dim, degree, 1, components, Number> fe_eval_cg(data_composite, 0);
  FEEvaluation<dim, degree, 1, components, Number> fe_eval_dg(data_composite, 1);

  VectorType vec__dg;
  data_composite.initialize_dof_vector(vec__dg, 1);

  for(unsigned int cell = 0; cell < data_composite.n_macro_cells(); ++cell)
  {
    fe_eval_cg.reinit(cell);
    fe_eval_dg.reinit(cell);

    fe_eval_cg.read_dof_values(src);

    for(unsigned int i = 0; i < fe_eval_cg.static_dofs_per_cell; i++)
      fe_eval_dg.begin_dof_values()[i] = fe_eval_cg.begin_dof_values()[i];

    fe_eval_dg.distribute_local_to_global(vec__dg);
  }
  dst.copy_locally_owned_data_from(vec__dg);
}

template<int dim, typename Number, typename VectorType, int components>
void
MGTransferMFC<dim, Number, VectorType, components>::interpolate(const unsigned int level,
                                                                VectorType &       dst,
                                                                const VectorType & src) const
{
  (void)level;

  switch(this->fe_degree)
  {
      // clang-format off
    case  1: do_interpolate< 1>(dst, src); break;
    case  2: do_interpolate< 2>(dst, src); break;
    case  3: do_interpolate< 3>(dst, src); break;
    case  4: do_interpolate< 4>(dst, src); break;
    case  5: do_interpolate< 5>(dst, src); break;
    case  6: do_interpolate< 6>(dst, src); break;
    case  7: do_interpolate< 7>(dst, src); break;
    case  8: do_interpolate< 8>(dst, src); break;
    case  9: do_interpolate< 9>(dst, src); break;
    case 10: do_interpolate<10>(dst, src); break;
    case 11: do_interpolate<11>(dst, src); break;
    case 12: do_interpolate<12>(dst, src); break;
    case 13: do_interpolate<13>(dst, src); break;
    case 14: do_interpolate<14>(dst, src); break;
    case 15: do_interpolate<15>(dst, src); break;
    default:
      AssertThrow(false, ExcMessage("MGTransferMFC::interpolate not implemented for this degree!"));
      // clang-format on
  }
}

template<int dim, typename Number, typename VectorType, int components>
void
MGTransferMFC<dim, Number, VectorType, components>::restrict_and_add(const unsigned int /*level*/,
                                                                     VectorType &       dst,
                                                                     const VectorType & src) const
{
  switch(this->fe_degree)
  {
      // clang-format off
    case  1: do_restrict_and_add< 1>(dst, src); break;
    case  2: do_restrict_and_add< 2>(dst, src); break;
    case  3: do_restrict_and_add< 3>(dst, src); break;
    case  4: do_restrict_and_add< 4>(dst, src); break;
    case  5: do_restrict_and_add< 5>(dst, src); break;
    case  6: do_restrict_and_add< 6>(dst, src); break;
    case  7: do_restrict_and_add< 7>(dst, src); break;
    case  8: do_restrict_and_add< 8>(dst, src); break;
    case  9: do_restrict_and_add< 9>(dst, src); break;
    case 10: do_restrict_and_add<10>(dst, src); break;
    case 11: do_restrict_and_add<11>(dst, src); break;
    case 12: do_restrict_and_add<12>(dst, src); break;
    case 13: do_restrict_and_add<13>(dst, src); break;
    case 14: do_restrict_and_add<14>(dst, src); break;
    case 15: do_restrict_and_add<15>(dst, src); break;
    default:
      AssertThrow(false, ExcMessage("MGTransferMFC::restrict_and_add not implemented for this degree!"));
      // clang-format on
  }
}

template<int dim, typename Number, typename VectorType, int components>
void
MGTransferMFC<dim, Number, VectorType, components>::prolongate(const unsigned int /*level*/,
                                                               VectorType &       dst,
                                                               const VectorType & src) const
{
  switch(this->fe_degree)
  {
      // clang-format off
    case  1: do_prolongate< 1>(dst, src); break;
    case  2: do_prolongate< 2>(dst, src); break;
    case  3: do_prolongate< 3>(dst, src); break;
    case  4: do_prolongate< 4>(dst, src); break;
    case  5: do_prolongate< 5>(dst, src); break;
    case  6: do_prolongate< 6>(dst, src); break;
    case  7: do_prolongate< 7>(dst, src); break;
    case  8: do_prolongate< 8>(dst, src); break;
    case  9: do_prolongate< 9>(dst, src); break;
    case 10: do_prolongate<10>(dst, src); break;
    case 11: do_prolongate<11>(dst, src); break;
    case 12: do_prolongate<12>(dst, src); break;
    case 13: do_prolongate<13>(dst, src); break;
    case 14: do_prolongate<14>(dst, src); break;
    case 15: do_prolongate<15>(dst, src); break;
    default:
      AssertThrow(false, ExcMessage("MGTransferMFC::prolongate not implemented for this degree!"));
      // clang-format on
  }
}

typedef dealii::LinearAlgebra::distributed::Vector<float>  VectorTypeFloat;
typedef dealii::LinearAlgebra::distributed::Vector<double> VectorTypeDouble;

template class MGTransferMFC<2, float, VectorTypeFloat, 1>;
template class MGTransferMFC<2, float, VectorTypeFloat, 2>;

template class MGTransferMFC<3, float, VectorTypeFloat, 1>;
template class MGTransferMFC<3, float, VectorTypeFloat, 3>;

template class MGTransferMFC<2, double, VectorTypeDouble, 1>;
template class MGTransferMFC<2, double, VectorTypeDouble, 2>;

template class MGTransferMFC<3, double, VectorTypeDouble, 1>;
template class MGTransferMFC<3, double, VectorTypeDouble, 3>;
