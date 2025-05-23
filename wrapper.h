#ifndef WRAPPER_H
#define WRAPPER_H

#ifdef __cplusplus
# define EXTERN_C extern "C"
#else
# define EXTERN_C
#endif

EXTERN_C double wrapped_heuristic(Problem *P0, Problem *P, BabNode *node, int *x, int num);
EXTERN_C int wrapped_read_data(const char *instance);
EXTERN_C void clean_python_references(void);
EXTERN_C void copy_solution(void);

#endif