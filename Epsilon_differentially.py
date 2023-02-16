import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
import pickle
from scipy.stats import bernoulli
import Generate_data
import multiprocessing
from itertools import product
from multiprocessing import Pool
import sys


# Parameters of the network
# Parameters of College Network



n = 2426  # Number of nosde
m = 1000  # Number of nodes samples to built matrix X for each influence cascade
p_ic = 0.05  # Probability for Independent Cascade Model
s_bulk = 80  # Number of Influence cascade model
# epsilon = 1
penalty = 0   # Penalty value for the regularization method
runs_alg = 5
algo_arg = 3
number_CPU = 320

# n = 3225                    # Number of nodes
# m = 250                     # Number of nodes samples to built matrix X for each influence cascade
# k = 4                       # Number of the seed set
# p_ic = 0.1                # Probability for Independent Cascade Model
# s_bulk=10 
#                   # Number of Influence cascade model
# # # epsilon = 1e-1
# #


# G_base = nx.read_edgelist("Contact-diaries-network_data_2013.csv.gz",nodetype=int, data=(("Type", str),))
# G_base=G_base.to_undirected()

G_base = nx.read_edgelist("soc-hamsterster_v2.edgelist.txt")
G_base=G_base.to_undirected()
previous_n_labels=G_base.nodes()
new_labes=list(range(len(previous_n_labels)))
dictionary_mapping = dict(zip(previous_n_labels, new_labes))
G_base = nx.relabel_nodes(G_base, dictionary_mapping)


def i_x_s(S, x):

    if S == []:

        return 0

    else:

        sub_x = x[:, S]
        temp = np.sum(sub_x, axis=1)

        return sum(temp >= 1)


def m_zero(k):

    #   Function for expected value spread when no information is used
    
    vals=[]

    epsilon_list=[0.01, 0.1, 0.8]

    for phi in range(s_bulk):
        
        
        indices = list(G_base.nodes())
        greedy_seeds= list(np.random.choice(indices, size=k, replace=False))
        vals.append(i_x_s(greedy_seeds, matrices_X[phi])*(n/m))
       
       
    for eps in epsilon_list:

        name_file = "cnk" +str(k) + "alg3_"+str(0) + "epsilon_"+str(eps)+".csv"

        np.savetxt(name_file, vals, delimiter=",")


        name_file = "cnk" +str(k) + "alg4_"+str(0) + "epsilon_"+str(eps)+".csv"

        np.savetxt(name_file, vals, delimiter=",")

        name_file = "cnk" + str(k) + "algo5_" + str(0) + "epsilon_" + str(eps) + ".csv"

        np.savetxt(name_file, vals, delimiter=",")

        name_file = "cnk" + str(k) + "algo6_" + str(0) + "epsilon_" + str(eps) + ".csv"

        np.savetxt(name_file, vals, delimiter=",")

seed_set = []

A = set(G_base.nodes())


def exponential_probabilities(Ix, epsilon, n, m):

    #   Function that return exponential value of Ix

    return np.exp(((epsilon*n)/(m*2))*Ix)



def exp_diff_spread(k,eps,h,iter,matrices_X):

    Expected_spread = []

    c_matrix_C = copy.deepcopy(matrices_X)

    sub_matrix_X = c_matrix_C[iter][0:h, :]
    print("Realization " + str(iter) + "size of m " + str(h) + "_random_" +str(eps))

    # Times that algortihm runs for influence cascade, each value f m, and epsilon.

    for _ in range(runs_alg):

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
            table[:, 2] = probabilities / np.sum(probabilities)

            # Sampling a node with the probablities found above

            seed_set.append(int(np.random.choice(
                list(table[:, 0]), 1, p=list(table[:, 2]))[0]))

        Expected_spread.append(i_x_s(seed_set, c_matrix_C[iter]) * (n / m))

    name_file = "./alg3/cnk" + str(k) + "alg3_" + str(h) + "epsilon_" + str(eps) + "iter_" + str(iter) + ".csv"

    with open(name_file, 'wb') as f:
        np.savetxt(f, Expected_spread, fmt='%-7.8f', delimiter=',')



#-----------------------------------------------------------------

if __name__ == "__main__":

    matrices_X = Generate_data.load_matrices_X()

    list_k=[8,12]
    epsilon_values = [0.01, 0.1, 0.8]
    m_values = [25,50, 100, 200, 400, 600, 800]
    iter_val=list(range(s_bulk))

    arguments=list(product(list_k, epsilon_values,m_values,iter_val))

    arguments=[list(i) for i in arguments]

    new_args=[]
    for i in range(len(arguments)):
        temp=arguments[i]
        temp.append(matrices_X)
        new_args.append(temp)


    processes=[]


    with multiprocessing.Pool(processes=number_CPU) as pool:
        pool.starmap(exp_diff_spread, new_args)

    pool.close()
    pool.join()





