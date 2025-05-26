import numpy as np

from biqbin import default_heuristic

def heuristic(L0, L, xfixed, sol_X, x):

    #print(xfixed.sum())
    #print(L.shape)
    #print(L)
    #print(L0)
    #print(x)

    result = default_heuristic(L0, L, xfixed, sol_X, x)

    #print(x)

    # exit(1)

    return result