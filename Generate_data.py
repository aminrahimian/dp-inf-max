# Auxiliary script to generate matrices X from a networkx graph

import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
import pickle
from scipy.stats import bernoulli


# Settings

m = 250  # Number of nodes samples to built matrix X for each influence cascade
p_ic = 0.05  # Probability for Independent Cascade Model
N=1 # Number of Influence cascade model

G_base = nx.read_edgelist("soc-hamsterster_v2.edgelist.txt")
G_base=G_base.to_undirected()
previous_n_labels=G_base.nodes()
new_labes=list(range(len(previous_n_labels)))
dictionary_mapping = dict(zip(previous_n_labels, new_labes))
G_base = nx.relabel_nodes(G_base, dictionary_mapping)


def generate_live_arc_graph(adj_matrix, p_ic):
    UP = np.triu(adj_matrix, k=1)
    boolean_mask = np.random.choice([True, False], len(np.where(UP == 1)[0]), p=[1 - p_ic, p_ic])
    UP[np.where(UP == 1)[0][boolean_mask], np.where(UP == 1)[1][boolean_mask]] = 0
    cascade = UP + UP.T
    G = nx.from_numpy_matrix(cascade)

    return G

def list_live_arc_graph(adj_matrix, p_ic,N):
    #   This function generates a list where each element is an influence cascade network
    #   Input as the parameters defined globally.
    m_live_edges = [generate_live_arc_graph(adj_matrix, p_ic) for _ in range(N)]

    return m_live_edges


def load_subgraphs():
    # Function to load pickle files with the list of influence cascade networks

    with open('live_edges.pkl', 'rb') as f:
        m_live_edges = pickle.load(f)

    return m_live_edges

def built_matrix_X(m, live_edges):

    # This function generates the matrix X of dimension mxn
    # Where m and n are the global parameters

    indices = list(live_edges.nodes())
    n=len(indices)
    matrix_X = np.zeros((m, n))

    uniform_nodes = np.random.choice(indices, size=m, replace=True)

    for i in range(m):
        matrix_X[i, :] = np.array(
            [nx.has_path(live_edges, j, uniform_nodes[i]) for j in range(n)])

    return matrix_X

def generate_random_order(G_base):

    order_nodes=list(G_base.nodes())
    np.random.shuffle(order_nodes)

    file_name = 'order_nodes.pkl'

    with open(file_name, 'wb') as f:
        pickle.dump(order_nodes, f)

    return order_nodes

def load_order_nodes():

    with open('order_nodes.pkl', 'rb') as f:
        order_nodes = pickle.load(f)

    return order_nodes

def set_matrices( m, N, p_ic):
    #   This function generates and save a set of X matrices.
    #   The number of matrices are defined by N parameter
    #   Return a list with the matrices X and save it as pkl file

    G_float = copy.deepcopy(G_base)
    m_live_edges = save_graph_in_list(G_float, p_ic,N)
    output_matrix= [built_matrix_X(m,m_live_edges[i]) for i in range(N)]

    file_name='./matrices_X.pkl'

    with open(file_name, 'wb') as f:
        pickle.dump(output_matrix, f)

    return 0

def load_matrix_X(iter):
    # Function to load pickle file with already generate X matrices

    file_name = './Matricesx/matrices_X_' + str(iter) + '.pkl'

    with open(file_name, 'rb') as f:
        matrix_X = pickle.load(f)

    return matrix_X

def load_matrices_X():
    # Function to load pickle file with already generate X matrices

    file_name = 'Matrix_X.pkl'

    with open(file_name, 'rb') as f:
        matrix_X = pickle.load(f)

    return matrix_X


if __name__ == "__main__":

    set_matrices(m, N, p_ic)


