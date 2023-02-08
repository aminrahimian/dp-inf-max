import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
import pickle
from scipy.stats import bernoulli
# import matplotlib.pyplot as plt
import multiprocessing
import time

# Parameters of the network
# Parameters of College Network


n = 2426                    # Number of nosde
m = 250                     # Number of nodes samples to built matrix X for each influence cascade
k = 4                       # Number of the seed set
p_ic = 0.05                 # Probability for Independent Cascade Model
s_bulk= 80                 # Number of Influence cascade model
# epsilon = 1e-1


# # Parameters of the network
# # Parameters of Tweeter Network
# # #
# n = 3225                    # Number of nodes
# m = 250                     # Number of nodes samples to built matrix X for each influence cascade
# k = 8                       # Number of the seed set
# p_ic = 0.3                # Probability for Independent Cascade Model
# s_bulk=10                   # Number of Influence cascade model
# # # epsilon = 1e-1
# #


    #
# G_base = nx.read_edgelist("Contact-diaries-network_data_2013.csv.gz",nodetype=int, data=(("Type", str),))
# G_base=G_base.to_undirected()




G_base = nx.read_edgelist("soc-hamsterster_v2.edgelist.txt")
G_base=G_base.to_undirected()



# G_base = nx.read_edgelist("rt_occupy_v2.txt")
# G_base=G_base.to_undirected()
previous_n_labels=G_base.nodes()
new_labes=list(range(len(previous_n_labels)))
dictionary_mapping = dict(zip(previous_n_labels, new_labes))
G_base = nx.relabel_nodes(G_base, dictionary_mapping)




# previous_n_labels=G_base.nodes()
# new_labes=list(range(len(previous_n_labels)))
# dictionary_mapping = dict(zip(previous_n_labels, new_labes))
# G_base = nx.relabel_nodes(G_base, dictionary_mapping)
#


# deg_dis=[G_base.degree(i) for i in range(1,1899)]



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
    #   This function generates a list where each element is an influence cascade network
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

    for i in range(m):
        matrix_X[i, :] = np.array(
            [nx.has_path(m_live_edges[i], j, uniform_nodes[i]) for j in range(n)])

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





def bulk_matrices( m, n, p_ic,iter):
    #   This function generates and save a set of X matrices.
    #   The number of matrices are defined by bulk_matrices parameter
    #   Return a list with the matrices X and save it as pkl file


    G_float = copy.deepcopy(G_base)
    m_live_edges = save_graph_in_list(m, G_float, p_ic)
    output_matrix= (built_matrix_X(m, n, m_live_edges))


    file_name='./Matricesx/matrices_X_' +str(iter) + '.pkl'

    with open(file_name, 'wb') as f:
        pickle.dump(output_matrix, f)

    return output_matrix




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

    processes=[]

    for t in range(s_bulk):

        p=multiprocessing.Process(target=bulk_matrices, args=[ m, n, p_ic,t])
        p.start()
        processes.append(p)


    for process in processes:

        process.join()

    # p1=multiprocessing.Process(target=bulk_matrices, args=[ m, n, p_ic,0])
    # p2=multiprocessing.Process(target=bulk_matrices, args=[ m, n, p_ic,1])
    #
    # p1.start()
    # p2.start()
    #
    # p1.join()
    # p2.join()

    Matrix_X=[]

    for i in range(s_bulk):

        Matrix_X.append(load_matrix_X(i))


    file_name = 'Matrix_X' + '.pkl'

    with open(file_name, 'wb') as f:

        pickle.dump(Matrix_X,f)



    with open('Matrix_X.pkl', 'rb') as f:
        matrix_X = pickle.load(f)

    # bulk_matrices(s_bulk, m, n, p_ic,1)

