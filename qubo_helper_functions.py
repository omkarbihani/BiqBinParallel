import numpy as np
import scipy as sp
import joblib


def qubo_to_maxcut(q: np.ndarray):
    q_sym = 1/2*(q.T + q)
    Qe_plus_c = -np.array([(np.sum(q_sym, 1))])
    np.fill_diagonal(q_sym, 0)

    return np.block([
        [np.zeros((1, 1)), Qe_plus_c],
        [Qe_plus_c.T,     q_sym]
    ])


def maxcut_to_L(Adj):
    num_vertices = Adj.shape[0]

    # Initialize Laplacian matrix L
    L = np.zeros((num_vertices, num_vertices))

    # Compute Adje = Adj * e (sum of each row)
    Adje = Adj.sum(axis=1)

    # Compute Diag(Adje) - diagonal matrix with Adje values
    tmp = np.diag(Adje)

    # Construct Laplacian matrix: L = tmp - Adj, excluding the last row and column
    L[:-1, :-1] = tmp[:-1, :-1] - Adj[:-1, :-1]

    # Compute vector parts and constant part
    sum_row = L[:-1, :-1].sum(axis=1)
    sum_total = sum_row.sum()

    # Fill vector parts and constant term
    L[:-1, -1] = sum_row
    L[-1, :-1] = sum_row
    L[-1, -1] = sum_total

    return L


def maxcut_solution_to_qubo_solution(_x: np.ndarray, n_mc: int):
    _x_mc = np.array(_x, dtype=int)-1
    x_mc_sol = np.ones(n_mc)
    x_mc_sol[_x_mc] = -1
    x_mc_sol *= -x_mc_sol[0]
    x_mc_sol
    y = 1/2*(x_mc_sol+1)[1:]
    return y


def read_maxcut_input(filename):
    with open(filename, 'r') as f:
        # Read number of vertices and edges
        num_vertices, num_edges = map(int, f.readline().split())

        # Allocate adjacency matrix as a contiguous array (row-major)
        adj_matrix = np.zeros(
            (num_vertices, num_vertices), dtype=np.float64)

        # Read edges
        for _ in range(num_edges):
            i, j, weight = f.readline().split()
            i, j = int(i) - 1, int(j) - 1  # Convert to zero-based indexing
            weight = float(weight)

            adj_matrix[i, j] = weight
            adj_matrix[j, i] = weight  # Since the graph is undirected

        # Create the MaxCutInputData struct
        return adj_matrix
