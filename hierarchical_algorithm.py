import numpy as np
import sympy as sp
from matrix import rref_null, row_delete, basis_matrix, null_identity
from graph_preprocess import dag_processing

def path_homology(edgelist: list[tuple[str, str]], calculate_basis: bool = True) -> tuple[int, list[list[str]]]:
    """
    This function computes the ll-dimensional path homology of the graph from the edge list, where ll is
    the maximum length of the longest path in the graph. And a basis will be calculated if specified.

    Args:
        edgelist (list[tuple[str, str]]): A list of tuples where each tuple represents a directed edge in the graph.
        calculate_basis (bool):           A flag indicating whether to calculate the basis of the null space. Default is True.

    Returns:
        tuple[int, int, list[list[str]]], or tuple[int, int, str]:
            ll (int):                     The order of the path homology group to be computed.
            sum_dim (int):                   The sum of the highest order of the path homology group across all subgraphs.
            basis (list[list[str]] or None): A list of basis vectors if `calculate_basis` is True; otherwise, `None`.
    """
    subgraph_dict, N, ll, num_subgraphs, graph_list = dag_processing(edgelist)
    sum_dim = 0
    basis = [] if calculate_basis else None
    if not subgraph_dict:
        return ll, num_subgraphs, basis
    assert len(N) == num_subgraphs
    for idx_subgraphs in range(num_subgraphs):
        if any(num == 1 for num in N[idx_subgraphs]):
            continue
        dim = N[idx_subgraphs][0] - 1
        V = np.vstack((-np.ones((1, dim), dtype = np.int8), np.identity(dim, dtype = np.int8)))    
        partition = {node: 1 for node in subgraph_dict[idx_subgraphs][0]}

        if calculate_basis:
            keys = list(partition.keys())
            symbols = {key: sp.symbols(key) for key in keys}
            basis_iteration = sp.Matrix([symbols[key]- symbols[keys[0]] for key in keys[1:]]).T
        for layer_idx in range(1, ll + 1):  
            col_list, P, Status = [], {}, False
            V_rows = V.shape[0]
            for idx, node in enumerate(subgraph_dict[idx_subgraphs][layer_idx]):
                predecessors = list(graph_list[idx_subgraphs].predecessors(node))
                sum_pred = sum(partition.get(node, 0) for node in predecessors)
                
                if sum_pred == V_rows:
                    P[node] = dim
                    col_list.append(-np.identity(dim, dtype = np.int8))
                    if not Status:
                        Status, index = True, idx
                else:
                    Reduced_matrix = row_delete(V, partition, predecessors)
                    A_x = rref_null(Reduced_matrix)

                    if A_x.size == 0:
                        P[node] = 0
                        continue
                
                    dim_A_x = A_x.shape[1]
                    P[node] = dim_A_x
                    if dim_A_x == dim:
                        col_list.append(-np.identity(dim, dtype = np.int8))
                        if not Status:
                            Status, index = True, idx
                    else:
                        col_list.append(-A_x)

            if not col_list:
                dim = 0
                break
            
            A = np.hstack(col_list)
            if Status:
                V, dim = null_identity(P, index, A, dim)                 
            else:
                V = rref_null(A)
                dim = V.shape[1]
            if dim == 0:
                break
            if calculate_basis:
                basis_iteration = basis_iteration * basis_matrix(col_list, P, V)
            partition = P
            
        sum_dim += dim
        if calculate_basis and dim != 0:
            basis.extend(basis_iteration.tolist())
    
    return ll, sum_dim, basis