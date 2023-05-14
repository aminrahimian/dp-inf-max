# Auxiliary script to generate matrices X from a networkx graph

import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
import pickle
from scipy.stats import bernoulli


# Settings

m = 10  # Number of nodes samples to built matrix X for each influence cascade
p_ic = 0.05  # Probability for Independent Cascade Model
N=5 # Number of Influence cascade model

G_base = nx.read_edgelist("soc-hamsterster_v2.edgelist.txt")
G_base=G_base.to_undirected()

adj_matrix_init=nx.to_numpy_array(G_base)

def generate_live_arc_graph(adj_matrix, p_ic,N):

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

def load_subgraphs():
    # Function to load pickle files with the list of influence cascade networks

    with open('live_edges.pkl', 'rb') as f:
        m_live_edges = pickle.load(f)

    return m_live_edges

def built_matrix_X(m, live_arc_graph):

    # This function generates the matrix X of dimension mxn
    # Where m and n are the global parameters

    indices = list(live_arc_graph.nodes())
    n=len(indices)
    matrix_X = np.zeros((m, n))

    uniform_nodes = np.random.choice(indices, size=m, replace=True)

    for i in range(m):
        matrix_X[i, :] = np.array(
            [nx.has_path(live_arc_graph, j, uniform_nodes[i]) for j in range(n)])

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

def set_matrices(adj_matrix_init,p_ic,m, N ):
    #   This function generates and save a set of X matrices.
    #   The number of matrices are defined by N parameter
    #   Return a list with the matrices X and save it as pkl file
    list_live_arcs=generate_live_arc_graph(adj_matrix_init,p_ic,N)


    output_matrix= [built_matrix_X(m,list_live_arcs[i]) for i in range(N)]

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


    set_matrices(adj_matrix_init,p_ic,m, N)


