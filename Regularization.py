import Randomize_version
import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
from scipy.stats import bernoulli
from scipy.special import comb
import pickle
import Generate_data
from itertools import product
import multiprocessing
from multiprocessing import Pool
import sys
import os

# Parameters of the network

# Parameters of College Network

n = 2426  # Number of nosde
m = 250  # Number of nodes samples to built matrix X for each influence cascade
p_ic = 0.05  # Probability for Independent Cascade Model
s_bulk = 80  # Number of Influence cascade model
runs_alg = 5
number_CPU = 8

if __name__ == "__main__":

    save_data_x_tilde=False
    matrices_X = Generate_data.load_matrices_X()
    order_nodes=Generate_data.load_order_nodes()

    epsilon_values = [0.8]
    list_k = [4]
    m_values = [50]

    list_s_bulk = list(range(s_bulk))

    algo_arg_list = [601, 602, 603, 604, 605, 6075, 61, 6125,615, 6175, 62, 62_5, 63, 64, 65, 675,610, 620, 650]

    dict_penalties={601: 0.1, 602: 0.2,603:0.3, 604:0.4, 605:0.5, 6075:0.75, 61: 1,
                    6125:1.25,615:1.5, 6175:1.75, 62: 2, 62_5:2.5, 63: 3, 64:4, 65:5, 675:7.5,610:10, 620:20, 650:50}


    # m_zero(runs_alg, s_bulk, algo_arg, epsilon_values, list_k)

    file_name = "Matrices_X_tilde.pkl"

    with open(file_name, 'rb') as f:
        dict_matrices_x_tilde = pickle.load(f)


    dict_set_of_matrices = {}

    for algo_arg in algo_arg_list:

        for eps in epsilon_values:

            for k in list_k:

                rho = (np.exp(eps) + 1) ** -1
                dict_set_of_matrices[(k,eps,algo_arg)]=Randomize_version.set_of_matrices(k, rho, algo_arg,dict_penalties[algo_arg])



    # sys.exit(1)

    arguments = list(product(list_k, epsilon_values, m_values,list_s_bulk, algo_arg_list))

    arguments = [list(i) for i in arguments]

    for i in arguments:
        i.append(dict_penalties[i[-1]])


    # generate_matrix_x_tilde(1, 0.01, m, n, matrices_X)


    new_args = []
    for i in range(len(arguments)):
        temp = arguments[i]
        temp.append(dict_set_of_matrices)
        temp.append(dict_matrices_x_tilde)
        temp.append(matrices_X)
        temp.append(order_nodes)
        new_args.append(temp)


    with multiprocessing.Pool(processes=number_CPU) as pool:
        pool.starmap(Randomize_version.randomize_response, new_args)

    pool.close()
    pool.join()





