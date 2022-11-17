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

G_base = nx.read_edgelist("email-Eu-core.txt.gz",nodetype=int, data=(("Type", str),))

def live_edges_saving(n, p, G):

    #    Function to generate Influence cascades.
    #   Input parameter are n: number of nodes, p: probability of the ICM, the network G.

    L=len(G.edges)
   
    tabular=np.zeros((L,1))
    

    Lista=bernoulli.rvs(p=p, size=L)
    tabular=pd.DataFrame(tabular)
    tabular[0]=G.edges()
    
    sampled_edges=list(tabular[[i==1 for i in Lista]][0])
    
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

def bulk_matrices(s_bulk,m_live_edges,m,n,p_ic):

    #   This function generates and save a set of X matrices.
    #   The number of matrices are defined by bulk_matrices parameter

    bulk_matrices=[]
    
    for i in range(s_bulk):
    
        G_float=copy.deepcopy(G_base)
        m_live_edges = save_graph_in_list(m, G_float, p_ic)
        bulk_matrices.append(built_matrix_X(m, n, m_live_edges))

    with open('matrices_X.pkl', 'wb') as f:
        pickle.dump(bulk_matrices, f)

    return bulk_matrices


def load_matrix_X():
    # Function to load pickle file with already generate X matrices
    with open('matrices_X.pkl', 'rb') as f:
        matrix_X = pickle.load(f)

    return matrix_X




def i_x_s(S, x):

    if S == []:

        return 0

    else:

        sub_x = x[:, S]
        temp = np.sum(sub_x, axis=1)

        return sum(temp >= 1)


matrices_X = load_matrix_X()


def m_zero(k,epsilon):

    #   Function for expected value spread when no information is used
    
    vals=[]
    
    for phi in range(40):
        
        
        indices = list(G_base.nodes())
        greedy_seeds= list(np.random.choice(indices, size=k, replace=False))
        vals.append(i_x_s(greedy_seeds, matrices_X[phi])*(n/m))
       
       
     
    name_file = "k" +str(k) + "alg3_"+str(0) + "epsilon_"+str(eps)+".csv"
    
    np.savetxt(name_file, vals, delimiter=",")


    name_file = "k" +str(k) + "alg4_"+str(0) + "epsilon_"+str(eps)+".csv"
    
    np.savetxt(name_file, vals, delimiter=",")




seed_set = []

A = set(G_base.nodes())


def exponential_probabilities(Ix, epsilon, n, m):

    #   Function that return exponential value of Ix

    return np.exp(((epsilon*n)/(m*2))*Ix)




#-----------------------------------------------------------------

# List that defines the number of rows of X used.
m_values = [10,20,30,40,50,80]

# List of the epsilon values
epsilon_values = [0.01, 1]



for eps in epsilon_values:
    
    m_zero(k,eps)

    # Iteration over the size of the submatrix of X of size m.

    for h in m_values:

        Expected_spread = []

        c_matrix_C = copy.deepcopy(matrices_X)

        # Iteration over the number of influence cascades

        for w in range(40):
            
            sub_matrix_X = c_matrix_C[w][0:h, :]
            print("Realization " + str(w) + "size of m " + str(h))

            # Times that algortihm runs for influence cascade, each value f m, and epsilon.
            
            for chi in range(5):
                
                seed_set = []

                # Iteration over the size of the seed set
    
                for i in range(k):
    
                    Q = copy.deepcopy(A)
                    B = set(seed_set)
                    S_union_v = copy.deepcopy(seed_set)
                    seed_set_s = copy.deepcopy(seed_set)
    
                    iter_set = Q.difference(B)
                    list_S_minus_seed = list(iter_set)
                    table = np.zeros((len(iter_set), 3))


                    # Iteration over the nodes that have not been selected yet
    
                    for j in range(len(iter_set)):
                        
                    
                        S_union_v.append(list_S_minus_seed[j])
    
                        # print("seeding " +str(i_x_s(S_union_v, sub_matrix_X)))
                        # print("which j value " +str(j))
                        table[j, 1] = i_x_s(S_union_v, sub_matrix_X) - \
                            i_x_s(seed_set_s, sub_matrix_X)
                        table[j, 0] = list_S_minus_seed[j]
    
                        S_union_v = copy.deepcopy(seed_set)
    
                    # Transform from I_x_s values to probabilities

                    probabilities = np.array([exponential_probabilities(
                        table[f, 1], eps, n, h) for f in range(len(iter_set))])
                    table[:, 2] = probabilities/np.sum(probabilities)

                    #Sampling a node with the probablities found above

                    seed_set.append(int(np.random.choice(
                        list(table[:, 0]), 1, p=list(table[:, 2]))[0]))
    
                Expected_spread.append(i_x_s(seed_set, c_matrix_C[w])*(n/m))

        
        name_file = "k" +str(k) + "alg3_"+str(h) + "epsilon_"+str(eps)+".csv"
        
        np.savetxt(name_file, Expected_spread, delimiter=",")



        
        


