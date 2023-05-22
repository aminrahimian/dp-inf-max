# script to model exponential mechanism and randomized response (all variations) algorithms

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




def expect_spread_exp_mechanism(iter_matrix,m,k,epsilon):
    """
       Returns I_x(S) using exponential mechanism given the m influence samples,
       total number of seeds k, and privacy budget epsilon.
       """


    if m == 0:
        n_nodes = iter_matrix.shape[1]
        seed_set = random.choices(np.arange(iter_matrix.shape[1]), k=k)

        return np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0])

    else:

        # print("Iteration m:" + str(m)+ ",k: "+ str(k)+ " ,e: "+str(epsilon))
        iter_matrix_m=iter_matrix[0:m,:]
        n_nodes=iter_matrix.shape[1]
        seed_set=np.array([],dtype=int)

        for i in range(k):

            node_list=np.arange(0,n_nodes, 1, dtype=int)
            candidates_node_list=node_list[np.logical_not(np.isin(node_list, seed_set))]

            boolean_mask = np.where(np.sum(iter_matrix_m[:, seed_set], axis=1) > 0, False, True)
            tile_boolean_mask = np.tile(boolean_mask, (n_nodes - len(seed_set), 1)).T
            boolean_x_diff_s = np.array(iter_matrix_m[:, candidates_node_list], dtype=bool)


            weights=np.array([round(np.exp((epsilon*m*i)/(2*n_nodes)), 2) for i in np.sum(np.logical_and(boolean_x_diff_s, tile_boolean_mask), axis=0)])
            weights = np.nan_to_num(weights, posinf=np.inf)
            print(" weights infsda " + str((weights)))
            weights=weights/np.sum(weights)
            # weights=np.nan_to_num(weights)
            print("sum weights " + str(np.sum(weights)))
            print("weights after " + str((weights)))
            node_nu=int(np.random.choice(candidates_node_list,1,p=weights)[0])
            seed_set=np.append(seed_set,node_nu)

        # print(str(np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0))*(n_nodes/iter_matrix.shape[0])))
        return  np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0))*(n_nodes/iter_matrix.shape[0])





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
    dict_matrices_c=dict.fromkeys(keys_c)

    for k,eps in keys_c:

        l=k+1
        C = np.zeros((l, l))
        rho = 1 / (np.exp(eps) + 1)

        for a,b in product(range(l), range(l)):

            a+=1
            b+=1
            start = max([0, (b - a)])
            end = min([(l - a), b])

            C[(a-1), (b-1)] = np.sum(
                [comb(b, i) * comb((l - b), (a - b + i)) * (rho ** (a - b + 2 * i)) * (1 - rho) ** (l - a + b - 2 * i)
                 for i in np.arange(start, end + 1)])

        dict_matrices_c[(k,eps)]=np.linalg.inv(C)

    return dict_matrices_c

def expect_spread_randomized_resp(iter_matrix, m, k, epsilon,dict_matrices_c,runs_alg):

    # print("This is matrix x " +str(iter_matrix[0:m,:]))

    mu_ixs = []

    for t in range(runs_alg):

        # print("runnnung " +str(t))

        if m == 0:

            n_nodes = iter_matrix.shape[1]
            seed_set = random.choices(np.arange(iter_matrix.shape[1]), k=k)

            # print(str(np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (
            #         n_nodes / iter_matrix.shape[0])))

            mu_ixs.append(
                np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0]))
            # return np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0])

        else:

            iter_matrix_m = generate_x_tilde(iter_matrix[0:m, :], epsilon)

            print("Coincidencia" + str(np.sum(iter_matrix_m!=iter_matrix[0:m,:])))


            n_nodes = iter_matrix.shape[1]
            seed_set = np.array([], dtype=int)
            # seed_set = np.array([10], dtype=int)

            for i in range(k):

                node_list = np.arange(0, n_nodes, 1, dtype=int)
                candidates_node_list = node_list[np.logical_not(np.isin(node_list, seed_set))]

                partial_s = np.sum(iter_matrix_m[:, seed_set], axis=1)
                rep_partial_s = np.broadcast_to(partial_s, ((n_nodes - len(seed_set)), m)).T
                addition = np.add(rep_partial_s, iter_matrix_m[:, candidates_node_list])

                vector_f_tilde = np.zeros(((len(seed_set) + 1), (n_nodes - len(seed_set))))

                for f in np.arange(len(seed_set) + 1):
                    vector_f_tilde[f, :] = np.sum(np.where(addition == f, True, False), axis=0) / m

                matrix_C = dict_matrices_c[(i, epsilon)]
                vector_f = matrix_C @ vector_f_tilde

                indice = np.max(vector_f[0, :]) == vector_f[0, :]
                node_nu = random.choices(candidates_node_list[indice], k=1)[0]
                seed_set = np.append(seed_set, node_nu)

                # print("seed va por: " +str(node_nu))

            # print(str(np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0])))
            mu_ixs.append(
                np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0]))

            # return np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0])

    # print(" for  m: " + str(m) + ' k: ' + str(k) + 'eps: ' + str(epsilon) + 'having: ' + str(mu_ixs))
    return (mu_ixs)


def expect_spread_randomized_resp_wpp(iter_matrix, m, k, epsilon, dict_matrices_c, runs_alg):
    mu_ixs = []

    for t in range(runs_alg):

        # print("runnnung " +str(t))

        if m == 0:

            n_nodes = iter_matrix.shape[1]
            seed_set = random.choices(np.arange(iter_matrix.shape[1]), k=k)

            # print(str(np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (
            #         n_nodes / iter_matrix.shape[0])))

            mu_ixs.append(
                np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0]))
            # return np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0])

        else:

            iter_matrix_m = generate_x_tilde(iter_matrix[0:m, :], epsilon)
            n_nodes = iter_matrix.shape[1]
            seed_set = np.array([], dtype=int)
            # seed_set = np.array([10], dtype=int)

            for i in range(k):

                node_list = np.arange(0, n_nodes, 1, dtype=int)
                candidates_node_list = node_list[np.logical_not(np.isin(node_list, seed_set))]

                partial_s = np.sum(iter_matrix_m[:, seed_set], axis=1)
                rep_partial_s = np.broadcast_to(partial_s, ((n_nodes - len(seed_set)), m)).T
                addition = np.add(rep_partial_s, iter_matrix_m[:, candidates_node_list])

                vector_f_tilde = np.zeros(((len(seed_set) + 1), (n_nodes - len(seed_set))))

                for f in np.arange(len(seed_set) + 1):
                    vector_f_tilde[f, :] = np.sum(np.where(addition == f, True, False), axis=0) / m

                matrix_C = dict_matrices_c[(i, epsilon)]
                vector_f = vector_f_tilde

                indice = np.max(vector_f[0, :]) == vector_f[0, :]
                node_nu = random.choices(candidates_node_list[indice], k=1)[0]
                seed_set = np.append(seed_set, node_nu)

                # print("seed va por: " +str(node_nu))

            # print(str(np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0])))
            mu_ixs.append(
                np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0]))

            # return np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0])

    # print(" for  m: " + str(m) + ' k: ' + str(k) + 'eps: ' + str(epsilon) + 'having: ' + str(mu_ixs))
    return (mu_ixs)


def expect_spread_greedy_algorithm(iter_matrix, k):
    """
         Returns I_x(S) using exponential mechanism given the m influence samples,
         total number of seeds k, and privacy budget epsilon.
         """

    # print("Iteration m:" + str(m)+ ",k: "+ str(k)+ " ,e: "+str(epsilon))
    iter_matrix_m = iter_matrix
    n_nodes = iter_matrix.shape[1]
    seed_set = np.array([], dtype=int)

    for i in range(k):
        node_list = np.arange(0, n_nodes, 1, dtype=int)
        candidates_node_list = node_list[np.logical_not(np.isin(node_list, seed_set))]

        boolean_mask = np.where(np.sum(iter_matrix_m[:, seed_set], axis=1) > 0, False, True)
        tile_boolean_mask = np.tile(boolean_mask, (n_nodes - len(seed_set), 1)).T
        boolean_x_diff_s = np.array(iter_matrix_m[:, candidates_node_list], dtype=bool)

        partial_sum = np.sum(np.logical_and(boolean_x_diff_s, tile_boolean_mask), axis=0)

        indice = np.max(partial_sum) == partial_sum
        node_nu = random.choices(candidates_node_list[indice], k=1)[0]
        seed_set = np.append(seed_set, node_nu)

    # print(str(np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0))*(n_nodes/iter_matrix.shape[0])))
    return np.sum(np.where(np.sum(iter_matrix[:, seed_set], axis=1) > 0, 1, 0)) * (n_nodes / iter_matrix.shape[0])
