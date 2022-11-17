import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
import pickle
from scipy.stats import bernoulli
import Generate_data

# Parameters of the network

n = 1005                    # Number of nodes
m = 100                     # Number of nodes samples to built matrix X for each influence cascade
k = 4                       # Number of the seed set
p_ic = 0.03                 # Probability for Independent Cascade Model
s_bulk=40                   # Number of Influence cascade model
# epsilon = 1e-1

G_base = nx.read_edgelist("email-Eu-core.txt.gz",nodetype=int, data=(("Type", str),))

def i_x_s(S, x):

    if S == []:

        return 0

    else:

        sub_x = x[:, S]
        temp = np.sum(sub_x, axis=1)

        return sum(temp >= 1)


matrices_X = Generate_data.load_matrix_X()


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



        
        


