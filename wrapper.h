#ifndef WRAPPER_H
#define WRAPPER_H

#include "biqbin.h"

#ifdef __cplusplus
# define EXTERN_C extern "C"
#else
# define EXTERN_C
#endif

EXTERN_C double wrapped_heuristic(Problem *P0, Problem *P, BabNode *node, int *x, int num);


#endif