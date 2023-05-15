# script to model exponential mechanism and randomized response algorithms

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


# models settings

m = 20  # Number of nodes samples to built matrix X for each influence cascade
s_bulk = 40  # Number of Influence cascade model
# epsilon = 1
runs_alg = 5
number_CPU = 4


def load_matrices_X():
    # Function to load pickle file with already generate X matrices

    file_name = 'matrices_X.pkl'

    with open(file_name, 'rb') as f:
        matrix_X = pickle.load(f)

    return matrix_X

list_matrices_X=load_matrices_X()

iter_arc_live=0
m=20
k=4
epsilon=0.1


def expect_spread_exp_mechanism(list_matrices_X,iter_arc_live,m,k,epsilon):
    """
       Returns I_x(S) using exponential mechanism given the m influence samples,
       total number of seeds k, and privacy budget epsilon.
       """

    iter_matrix=list_matrices_X[iter_arc_live][0:m,:]
    ref_matrix=list_matrices_X[iter_arc_live]
    n_nodes=iter_matrix.shape[1]
    seed_set=np.array([],dtype=int)

    for i in range(k):

        node_list=np.arange(0,n_nodes, 1, dtype=int)
        candidates_node_list=node_list[np.logical_not(np.isin(node_list, seed_set))]

        boolean_mask = np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, False, True)
        tile_boolean_mask = np.tile(boolean_mask, (n_nodes - len(seed_set), 1)).T
        boolean_x_diff_s = np.array(iter_matrix[:, candidates_node_list], dtype=bool)

        weights=np.array([np.exp((epsilon*m*i)/(2*n_nodes)) for i in np.sum(np.logical_and(boolean_x_diff_s, tile_boolean_mask), axis=0)])
        weights=weights/np.sum(weights)
        node_nu=int(np.random.choice(candidates_node_list,1,p=weights)[0])
        seed_set=np.append(seed_set,node_nu)


    return  np.sum(np.where(np.sum(ref_matrix[:, seed_set], axis=1) > 0, 1, 0))*(n_nodes/ref_matrix.shape[0])


t=[]


for i in range(40):
    t.append(expect_spread_exp_mechanism(list_matrices_X,2,80,4,0.1))

print(np.mean(t))

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




