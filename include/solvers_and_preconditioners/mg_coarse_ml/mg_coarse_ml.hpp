#include "mg_coarse_ml.h"

#include <navier-constants.h>

#include "../../operators/matrix_operator_base_new.h"

#if DIM_2
template class MGCoarseML<MatrixOperatorBaseNew<2, float>, float>;
template class MGCoarseML<MatrixOperatorBaseNew<2, float>, double>;
#endif

#if DIM_3
template class MGCoarseML<MatrixOperatorBaseNew<3, float>, float>;
template class MGCoarseML<MatrixOperatorBaseNew<3, float>, double>;
#endif
