# Script to model exponential mechanism and randomized response algorithms

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

m = 1000  # Number of nodes samples to built matrix X for each influence cascade
p_ic = 0.05  # Probability for Independent Cascade Model
s_bulk = 80  # Number of Influence cascade model
# epsilon = 1
runs_alg = 5
number_CPU = 4


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





