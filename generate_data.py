# generate matrices x from the adjency matrix of a graph

import networkx as nx
import numpy as np
import pickle

# different datasets require specific parameters that might require calibration.

# dataset_id = 'soc-hamsterster_v2'
# dataset_id='erdos_renyi'
dataset_id = 'email-Eu-core'

if dataset_id == 'soc-hamsterster_v2':

    dataset_name='soc-hamsterster_v2.csv'
    m = 5000         # number of influence samples
    p_ic = 0.035     # probability  independent cascade model (ICM)
    N = 10          # number of ICM realizations

if dataset_id == 'erdos_renyi':

    dataset_name = 'erdos_renyi.csv'
    m = 2000  # number of influence samples
    p_ic = 0.03 # probability  independent cascade model (ICM)
    N = 50 # number of ICM realizations

if dataset_id == 'email-Eu-core':

    dataset_name = 'email-Eu-core.csv'
    m = 2000  # number of influence samples
    p_ic = 0.0155  # probability  independent cascade model (ICM)
    N = 50  # number of ICM realizations



def generate_live_arc_graph(adj_matrix, p_ic,N):
    """
        Returns an influence cascade graph of N nodes
        using probability p_ic using the adjency matrix of the
        original graph.
       """

    list_live_arcs=[]

    for i in range(N):

        UP = np.triu(adj_matrix, k=1)
        boolean_mask = np.random.choice([True, False], len(np.where(UP == 1)[0]), p=[1 - p_ic, p_ic])
        UP[np.where(UP == 1)[0][boolean_mask], np.where(UP == 1)[1][boolean_mask]] = 0
        cascade = UP + UP.T
        G = nx.from_numpy_matrix(cascade)
        G= G.to_undirected()

        list_live_arcs.append(G)

    return list_live_arcs

def built_matrix_X(m, live_arc_graph):
    """
       Returns the matrix X of dimension mxn
       with m number of nodes sampled from the influence cascade and n number of nodes.
       """

    indices = list(live_arc_graph.nodes())
    n=len(indices)
    matrix_X = np.zeros((m, n))

    uniform_nodes = np.random.choice(indices, size=m, replace=True)

    for i in range(m):

        matrix_X[i, :] = np.array(
            [nx.has_path(live_arc_graph, j, uniform_nodes[i]) for j in range(n)])

    return (matrix_X)



def set_matrices(adj_matrix_init,p_ic,m, N ):
    """
    This method saves a list of N matrices X, i.e. N independent realizations
    of the influence cascade model.
    """
    list_live_arcs=generate_live_arc_graph(adj_matrix_init,p_ic,N)
    output_matrix= [built_matrix_X(m,list_live_arcs[i]) for i in range(N)]

    file_name='./matrices_x.pkl'

    with open(file_name, 'wb') as f:
        pickle.dump(output_matrix, f)

    return 0



if __name__ == "__main__":

    if dataset_id=='soc-hamsterster_v2':

        adj_matrix_init=np.genfromtxt("soc-hamsterster_v2.csv", delimiter=",")
        set_matrices(adj_matrix_init,p_ic,m, N)

    if dataset_id=='erdos_renyi':

        adj_matrix_init = nx.to_numpy_array(nx.erdos_renyi_graph(n=200, p=0.15, seed=100, directed=False))
        set_matrices(adj_matrix_init, p_ic, m, N)

    if dataset_id == 'email-Eu-core':
        adj_matrix_init = np.genfromtxt("email-Eu-core.csv", delimiter=",")
        set_matrices(adj_matrix_init, p_ic, m, N)
