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



n = 2426  # Number of nosde
m = 1000  # Number of nodes samples to built matrix X for each influence cascade
p_ic = 0.05  # Probability for Independent Cascade Model
s_bulk = 80  # Number of Influence cascade model
# epsilon = 1
penalty = 0   # Penalty value for the regularization method
runs_alg = 5
algo_arg = 3
number_CPU = 320



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





