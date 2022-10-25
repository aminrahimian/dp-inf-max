import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
import pickle
from scipy.stats import bernoulli

# Parameters


n = 1005
m = 100
k = 4
p_ic = 0.03
s_bulk=40
# epsilon = 1e-1


G_base = nx.read_edgelist("email-Eu-core.txt.gz",nodetype=int, data=(("Type", str),))



# n is the number of edges not nodes
def live_edges_saving(n, p, G):

    
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

    with open('live_edges.pkl', 'rb') as f:
        m_live_edges = pickle.load(f)

    return m_live_edges


def built_matrix_X(m, n, m_live_edges):

    matrix_X = np.zeros((m, n))
    indices = list(m_live_edges[0].nodes())

    uniform_nodes = np.random.choice(indices, size=m, replace=True)
    
    print(uniform_nodes)

    for i in range(m):

        matrix_X[i, :] = np.array(
            [nx.has_path(m_live_edges[i], j, uniform_nodes[i]) for j in range(n)])

    return matrix_X




def bulk_matrices():

    bulk_matrices=[]
    
    for i in range(s_bulk):
    
        G_float=copy.deepcopy(G_base)
        m_live_edges = save_graph_in_list(m, G_float, p_ic)
        bulk_matrices.append(built_matrix_X(m, n, m_live_edges))

    return bulk_matrices




def save_matrix_X():

    with open('matrices_X.pkl', 'wb') as f:
        pickle.dump(bulk_matrices, f)


def load_matrix_X():

    with open('matrices_X.pkl', 'rb') as f:
        matrix_X = pickle.load(f)

    return matrix_X


# save_matrix_X()
# bulk_matrices()




def i_x_s(S, x):

    if S == []:

        return 0

    else:

        sub_x = x[:, S]
        temp = np.sum(sub_x, axis=1)

        return sum(temp >= 1)


matrices_X = load_matrix_X()



def m_zero(k,epsilon):
    
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

    return np.exp(((epsilon*n)/(m*2))*Ix)


# chose a matrix according to the size m that we want to plot


m_values = [10,20,30,40,50,75]

epsilon_values = [0.01, 1]


for eps in epsilon_values:
    
    m_zero(k,eps)

    for h in m_values:

        Expected_spread = []

        c_matrix_C = copy.deepcopy(matrices_X)

     
        for w in range(40):
            
            sub_matrix_X = c_matrix_C[w][0:h, :]
            print("Realization " + str(w) + "size of m " + str(h))
            
            for chi in range(5):
                
                seed_set = []
                
    
                for i in range(k):
    
                    Q = copy.deepcopy(A)
                    B = set(seed_set)
                    S_union_v = copy.deepcopy(seed_set)
                    seed_set_s = copy.deepcopy(seed_set)
    
                    iter_set = Q.difference(B)
                    list_S_minus_seed = list(iter_set)
                    table = np.zeros((len(iter_set), 3))
    
                    for j in range(len(iter_set)):
                        
                    
                        S_union_v.append(list_S_minus_seed[j])
    
                        # print("seeding " +str(i_x_s(S_union_v, sub_matrix_X)))
                        # print("which j value " +str(j))
                        table[j, 1] = i_x_s(S_union_v, sub_matrix_X) - \
                            i_x_s(seed_set_s, sub_matrix_X)
                        table[j, 0] = list_S_minus_seed[j]
    
                        S_union_v = copy.deepcopy(seed_set)
    
                    # transform to probability
    
                    probabilities = np.array([exponential_probabilities(
                        table[f, 1], eps, n, h) for f in range(len(iter_set))])
                    table[:, 2] = probabilities/np.sum(probabilities)
    
                    seed_set.append(int(np.random.choice(
                        list(table[:, 0]), 1, p=list(table[:, 2]))[0]))
    
                Expected_spread.append(i_x_s(seed_set, c_matrix_C[w])*(n/m))

        
        name_file = "k" +str(k) + "alg3_"+str(h) + "epsilon_"+str(eps)+".csv"
        
        np.savetxt(name_file, Expected_spread, delimiter=",")



        
        


