# script to model exponential mechanism and randomized response algorithms

import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
import pickle
from scipy.stats import bernoulli
import multiprocessing
from multiprocessing import Pool
import sys
from math import comb
from itertools import product


# Parameters

# m m_values                    # samples from influence cascades
# iter_arc_live = 40           # number of influence cascade models
# epsilon = 1                  # privacy budget
# runs_alg = 5                 # times to run dp-algorithm fixing parameters
# number_CPU = 4               # multiprocessing parameters
# list_k                       # list size of seed set

# iter_matrix=list_matrices_x[0]
# iter_arc_live=0
# m=20
# k=4
# epsilon=0.1

def expect_spread_exp_mechanism(iter_matrix,m,k,epsilon):
    """
       Returns I_x(S) using exponential mechanism given the m influence samples,
       total number of seeds k, and privacy budget epsilon.
       """

    iter_matrix_m=iter_matrix[0:m,:]
    n_nodes=iter_matrix.shape[1]
    seed_set=np.array([],dtype=int)

    for i in range(k):

        node_list=np.arange(0,n_nodes, 1, dtype=int)
        candidates_node_list=node_list[np.logical_not(np.isin(node_list, seed_set))]

        boolean_mask = np.where(np.sum(iter_matrix_m[:, seed_set], axis=1) > 0, False, True)
        tile_boolean_mask = np.tile(boolean_mask, (n_nodes - len(seed_set), 1)).T
        boolean_x_diff_s = np.array(iter_matrix_m[:, candidates_node_list], dtype=bool)

        weights=np.array([np.exp((epsilon*m*i)/(2*n_nodes)) for i in np.sum(np.logical_and(boolean_x_diff_s, tile_boolean_mask), axis=0)])
        weights=weights/np.sum(weights)
        node_nu=int(np.random.choice(candidates_node_list,1,p=weights)[0])
        seed_set=np.append(seed_set,node_nu)

    print(str(np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0))*(n_nodes/iter_matrix.shape[0])))
    return  np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0))*(n_nodes/iter_matrix.shape[0])


def m_zero(k):

    #   Function for expected value spread when no information is used
    
    vals=[]

    epsilon_list=[0.01, 0.1, 0.8]
    #
    # for phi in range(s_bulk):
    #
    #
    #     indices = list(G_base.nodes())
    #     greedy_seeds= list(np.random.choice(indices, size=k, replace=False))
    #     vals.append(i_x_s(greedy_seeds, matrices_x[phi])*(n/m))
    #
    #
    # for eps in epsilon_list:
    #
    #     name_file = "cnk" +str(k) + "alg3_"+str(0) + "epsilon_"+str(eps)+".csv"
    #
    #     np.savetxt(name_file, vals, delimiter=",")
    #
    #
    #     name_file = "cnk" +str(k) + "alg4_"+str(0) + "epsilon_"+str(eps)+".csv"
    #
    #     np.savetxt(name_file, vals, delimiter=",")
    #
    #     name_file = "cnk" + str(k) + "algo5_" + str(0) + "epsilon_" + str(eps) + ".csv"
    #
    #     np.savetxt(name_file, vals, delimiter=",")
    #
    #     name_file = "cnk" + str(k) + "algo6_" + str(0) + "epsilon_" + str(eps) + ".csv"
    #
    #     np.savetxt(name_file, vals, delimiter=",")

def generate_x_tilde(matrix_x,epsilon):

    matrix_x_tilde=copy.deepcopy(matrix_x)
    p_flip=1/(np.exp(epsilon)+1)

    boolean_mask= np.random.choice([True, False],matrix_x_tilde.size,
                                   p=[p_flip, 1-p_flip]).reshape((matrix_x_tilde.shape[0],
                                                                  matrix_x_tilde.shape[1]))

    matrix_x_tilde[np.logical_and(np.array(matrix_x_tilde==1),boolean_mask)]=0
    matrix_x_tilde[np.logical_and(np.array(matrix_x_tilde==0),boolean_mask)]=1

    return matrix_x_tilde

def fill_matrix_c(list_k,epsilon_values):

    keys_c=list(product(list(range(int(max(list_k)))),epsilon_values))
    dict_matrices_C=dict.fromkeys(keys_c)

    for k,eps in keys_c:

        l=k+1
        C = np.zeros((l, l))
        rho = 1 / (np.exp(eps) + 1)

        for a,b in product(range(l), range(l)):

            start = max([0, (b - a)])
            end = min([(l - a), b])

            C[a, b] = np.sum(
                [comb(b, i) * comb((l - b), (a - b + i)) * (rho ** (a - b + 2 * i)) * (1 - rho) ** (l - a + b - 2 * i)
                 for i in np.arange(start, end + 1)])

        dict_matrices_C[(k,eps)]=np.linalg.inv(C)

    return dict_matrices_C








