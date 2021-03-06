#include "laplace_operator.h"

#include <navierstokes/config.h>

#if DIM_2 && DEGREE_1
template class LaplaceOperator<2, 1, double>;
#endif

#if DIM_2 && DEGREE_2
template class LaplaceOperator<2, 2, double>;
#endif

#if DIM_2 && DEGREE_3
template class LaplaceOperator<2, 3, double>;
#endif

#if DIM_2 && DEGREE_4
template class LaplaceOperator<2, 4, double>;
#endif

#if DIM_2 && DEGREE_5
template class LaplaceOperator<2, 5, double>;
#endif

#if DIM_2 && DEGREE_6
template class LaplaceOperator<2, 6, double>;
#endif

#if DIM_2 && DEGREE_7
template class LaplaceOperator<2, 7, double>;
#endif

#if DIM_2 && DEGREE_8
template class LaplaceOperator<2, 8, double>;
#endif

#if DIM_2 && DEGREE_9
template class LaplaceOperator<2, 9, double>;
#endif

#if DIM_3 && DEGREE_1
template class LaplaceOperator<3, 1, double>;
#endif

#if DIM_3 && DEGREE_2
template class LaplaceOperator<3, 2, double>;
#endif