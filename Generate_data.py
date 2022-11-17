import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
import pickle
from scipy.stats import bernoulli

# Parameters of the network

n = 1005                    # Number of nodes
m = 100                     # Number of nodes samples to built matrix X for each influence cascade
k = 4                       # Number of the seed set
p_ic = 0.03                 # Probability for Independent Cascade Model
s_bulk=40                   # Number of Influence cascade model
# epsilon = 1e-1



G_base = nx.read_edgelist("email-Eu-core.txt.gz", nodetype=int, data=(("Type", str),))


def live_edges_saving(n, p, G):
    #    Function to generate Influence cascades.
    #   Input parameter are n: number of nodes, p: probability of the ICM, the network G.

    L = len(G.edges)

    tabular = np.zeros((L, 1))

    Lista = bernoulli.rvs(p=p, size=L)
    tabular = pd.DataFrame(tabular)
    tabular[0] = G.edges()

    sampled_edges = list(tabular[[i == 1 for i in Lista]][0])

    # print(sampled_edges)
    ini_nodes = set(G.nodes())
    sG = G.edge_subgraph(sampled_edges).copy()
    s_nodes = set(sG.nodes())
    diff_nodes = list(ini_nodes.difference(s_nodes))

    for i in diff_nodes:
        sG.add_node(i)

    # path='Live_edges\live_edges'+ str(step)+'_'+str(sample_size)+'.pkl'
    # pickle.dump(sG, open(path, 'wb'))

    return sG


def save_graph_in_list(m, G_base, p_ic):
    #   This function generates a list where each element is a influence cascade network
    #   Input as the parameters defined globally.

    test_keys = list(G_base.nodes())
    test_values = list(range(len(test_keys)))
    dict_to_map = {test_keys[i]: test_values[i] for i in range(len(test_keys))}

    H = nx.relabel_nodes(G_base, dict_to_map)
    H = H.to_undirected()

    n = len(H.nodes())
    n_e = len(H.edges())

    m_live_edges = [live_edges_saving(n_e, p_ic, H) for i in range(m)]

    # with open('live_edges.pkl', 'wb') as f:
    #     pickle.dump(m_live_edges, f)

    return m_live_edges


# save_graph_in_list(m, G_base, p_ic)


def load_subgraphs():
    # Function to load pickle files with the list of Influence cascade networks

    with open('live_edges.pkl', 'rb') as f:
        m_live_edges = pickle.load(f)

    return m_live_edges


def built_matrix_X(m, n, m_live_edges):
    # This function generates the matrix X of dimension mxn
    # Where m and n are the global parameters

    matrix_X = np.zeros((m, n))
    indices = list(m_live_edges[0].nodes())

    uniform_nodes = np.random.choice(indices, size=m, replace=True)

    print(uniform_nodes)

    for i in range(m):
        matrix_X[i, :] = np.array(
            [nx.has_path(m_live_edges[i], j, uniform_nodes[i]) for j in range(n)])

    return matrix_X


def bulk_matrices(s_bulk, m_live_edges, m, n, p_ic):
    #   This function generates and save a set of X matrices.
    #   The number of matrices are defined by bulk_matrices parameter

    bulk_matrices = []

    for i in range(s_bulk):
        G_float = copy.deepcopy(G_base)
        m_live_edges = save_graph_in_list(m, G_float, p_ic)
        bulk_matrices.append(built_matrix_X(m, n, m_live_edges))

    with open('matrices_X.pkl', 'wb') as f:
        pickle.dump(bulk_matrices, f)

    return bulk_matrices
