/*
 * resolve_templates_double_2d.cpp
 *
 *  Created on: May 2, 2019
 *      Author: fehn
 */

#include "which_degrees.h"

#include <deal.II/matrix_free/evaluation_template_factory.templates.h>

template struct dealii::internal::FEEvaluationFactory<2,1,double>;
template struct dealii::internal::FEEvaluationFactory<2,2,double>;

template struct dealii::internal::FEFaceEvaluationFactory<2,1,double>;
template struct dealii::internal::FEFaceEvaluationFactory<2,2,double>;

// inverse mass
template struct dealii::internal::CellwiseInverseMassFactory<2,1,double>;
template struct dealii::internal::CellwiseInverseMassFactory<2,2,double>;
template struct dealii::internal::CellwiseInverseMassFactory<2,2+2,double>;

