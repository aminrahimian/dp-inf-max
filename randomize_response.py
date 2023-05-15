import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
from scipy.stats import bernoulli
from scipy.special import comb
import pickle
import generate_data
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
# epsilon = 1
penalty = 0   # Penalty value for the regularization method
runs_alg = 5
algo_arg = 4
number_CPU = 8
data_path_directory="/Users/carloshurtado/Documents"


# G_base=nx.read_edgelist("email-Eu-core.txt.gz", nodetype=int, data=(("Type", str),))
# G_base = nx.read_edgelist("Contact-diaries-network_data_2013.csv.gz",nodetype=int, data=(("Type", str),))
# G_base=G_base.to_undirected()

G_base = nx.read_edgelist("soc-hamsterster_v2.edgelist.txt")
G_base = G_base.to_undirected()
previous_n_labels = G_base.nodes()
new_labes = list(range(len(previous_n_labels)))
dictionary_mapping = dict(zip(previous_n_labels, new_labes))
G_base = nx.relabel_nodes(G_base, dictionary_mapping)


def i_x_s(S,x): 
    
    if S==[]:
    
        return 0
    
    else:
        
        sub_x= x[:,S] 
        temp=np.sum(sub_x, axis=1)
        
        return sum(temp>=1)
    

def generate_matrix_x_tilde(position, eps,m,n,list_matrices_x):
    
    #   Function to generate the matrix x_tilde from matrix X of dimensions mxn
    #   Rho is the flipping probability
    matrix_X=list_matrices_x[position]

    rho = (np.exp(eps) + 1) ** -1

    matrix_x_tilde_output=np.zeros((m,n))
    
    for i in range(m):
        
        row_i= matrix_X[i,:]
        
        flipped_row= 1- row_i
        
        for j in range(n):
            
            if bernoulli.rvs(p=rho, size=1)[0]==1:
                
                matrix_x_tilde_output[i,j]=flipped_row[j]
                
            else:
                
                matrix_x_tilde_output[i,j]=matrix_X[i,j]

        
    file_name='./Matrices_x_tilde/matrix_X_tilde' +str(position) + 'eps_' +str(eps)+'.pkl'

    with open(file_name, 'wb') as f:
        pickle.dump(matrix_x_tilde_output, f)

    
    return matrix_x_tilde_output


def combinatorics(a,b,l_curv,rho, l):

    # Function used in sum_combinatorics(a,b,l,rho) define value of the entry a,b in the matrix C
    
    return comb(b,l_curv)* comb((l-b), (a-b+l_curv))*(rho**(a-b + 2*l_curv))* ((1-rho)**(l-a+b -2*l_curv))

def sum_combinatorics(a,b,l,rho):

    # Function to define value of the entry a,b in the matrix C
    
    lower= max(0, b-a)
    upper= min((l-a),b)
    
    
    cont=0
    
    for i in range(lower, upper+1):

        cont+=combinatorics(a,b,i,rho,l)

    return cont

def fill_matrix_c(matrix_C, dim_matrix_c, rho,l):

    #   Function to built matrix C_{rho,l}
    
    for i in range(dim_matrix_c):
        
        for j in range(dim_matrix_c):
            
            matrix_C[i,j]=sum_combinatorics(i,j,l, rho)

    return matrix_C


def f_tilde_a(seed_set,m,a,matrix_x_tilde):

    # Function to calculate f_tilde defined in the line 1 of algorithm #5

    cont=0
    
    for i in range(m):
        
        x_t_tilde=np.sum(matrix_x_tilde[i,seed_set])
        
        if x_t_tilde==a:
            
            cont+=1 
        
        else:
            
            pass
        
    
    return cont*(1/m)


def set_of_matrices(k, rho, algo_arg,penalty):

    # Function that returns all the inverse matrices X in function of the seed set k
    if algo_arg==4:

        set_of_matrices = []

        for o in range(k + 1):

            dim_matrix_c = o + 1
            matrix_C = np.zeros((dim_matrix_c, dim_matrix_c))
            matrix_C_up = fill_matrix_c(matrix_C, dim_matrix_c, rho, o)
            # matrix_C_up=np.append(matrix_C_up,[ dim_matrix_c*[1]], axis=0)
            psd_inverse = np.linalg.inv(matrix_C_up)
            set_of_matrices.append(psd_inverse)

        return set_of_matrices


    elif algo_arg==5:

        set_of_matrices = []

        for o in range(k + 1):

            dim_matrix_c = o + 1
            set_of_matrices.append(np.eye(dim_matrix_c))

        return set_of_matrices

    else:

        set_of_matrices = []

        for o in range(k + 1):
            dim_matrix_c = o + 1
            matrix_C = np.zeros((dim_matrix_c, dim_matrix_c))
            matrix_C_up = fill_matrix_c(matrix_C, dim_matrix_c, rho, o)
            # matrix_C_up=np.append(matrix_C_up,[ dim_matrix_c*[1]], axis=0)
            matrix_penalty=penalty*np.eye(dim_matrix_c)
            psd_inverse =np.matmul(np.linalg.inv((matrix_C_up.T)@matrix_C_up + (matrix_penalty.T)@matrix_penalty),(matrix_C_up.T))
            set_of_matrices.append(psd_inverse)

        return set_of_matrices


def binary_converter(l, list_nodes):

    # Function for test purposes

    if l in list_nodes:
        
        return 1
    else:
        
        return 0


def f_tilded(seed_set, matrix_x_tilde,m):

    # Function to calculate f_tilde defined in line 1 algorithm 5

    if seed_set==[]:
        
        return np.array([1])

    else:

        f_tilde=[]
        dim_f_tilde=len(seed_set)+1
        sub_x = matrix_x_tilde[:, seed_set]
        vector_one=np.ones((1,len(seed_set)))
        counts= np.matmul(sub_x, vector_one.T)
        
        for j in range(dim_f_tilde):

            f_tilde.append(sum(counts==j)[0]*(1/m))

        return np.array(f_tilde)


# set_of_matrices_in=set_of_matrices(k,rho)

def j_value(seed_set,m,matrix_x_tilde,set_of_matrices_in, algo_arg, penalty):

    #   Function to calculate the value of f_0 using the seed set and matrix_x_tilde
    if not algo_arg==6:

        f_tilde=f_tilded(seed_set, matrix_x_tilde,m)
        f_0=np.matmul(set_of_matrices_in[len(seed_set)], f_tilde)[0]

    else:

        f_tilde = f_tilded(seed_set, matrix_x_tilde, m)
        f_0 = np.matmul(set_of_matrices_in[len(seed_set)], f_tilde)[0]


    return n*(1-f_0)

def my_indices(lst, item):

    # Auxiliary function to randomize selection process when max function is tight.

   return [i for i, x in enumerate(lst) if x == item]



def randomize_response(k,eps,h,position,algo_arg,penalty,dict_set_of_matrices,dict_matrices_x_tilde,matrices_X, order_nodes):

    rho = (np.exp(eps) + 1) ** -1
    set_of_matrices_array= dict_set_of_matrices[(k, eps, algo_arg)]

    Expected_spread = []

    sample_rows = np.random.choice(list(range(m)), size=h, replace=False)
    matrix_x_tilde = dict_matrices_x_tilde[eps][position][sample_rows, :]

    print("iteration  :" + str(position) + " m value" + str(h) + " eps " + str(eps) + "alg " + str(algo_arg))

    for _ in range(runs_alg):

        seed_set = []

        for w in range(k):

            # print("seed  :" +str(w))

            set_s = copy.deepcopy(seed_set)

            A = set(set_s)
            B = set(G_base.nodes())

            iter_set = B.difference(A)

            set_s_union_v = copy.deepcopy(seed_set)
            list_iter_set = list(iter_set)

            j_lists = []

            # Iteration over the nodes that have not been selected yet

            for i in range(len(list_iter_set)):
                set_s_union_v.append(i)

                j_lists.append(j_value(set_s_union_v, h, matrix_x_tilde, set_of_matrices_array, algo_arg, penalty))

                set_s_union_v = copy.deepcopy(seed_set)

            candidates = my_indices(j_lists, max(j_lists))

            sub_order = [order_nodes.index(i) for i in candidates]

            seed_set.append(order_nodes[min(sub_order)])



        Expected_spread.append(i_x_s(seed_set, matrices_X[position]) * (n / m))

    if ((algo_arg==4) or (algo_arg==3)):

        newpath = data_path_directory+ "/alg" +str(algo_arg)

        if not os.path.exists(newpath):
            os.makedirs(newpath)

        name_file = data_path_directory + "/alg" +str(algo_arg)+ "/cnk" + str(k) + "alg" +str(algo_arg)+"_" + str(h) + "epsilon_" + str(eps) + "iter_" + str(position) + ".csv"

        # np.savetxt(name_file, Expected_spread, delimiter=",")

        with open(name_file, 'wb') as f:
            np.savetxt(f, Expected_spread, fmt='%-7.8f', delimiter=',')

        return(1)


    else:

        newpath =  data_path_directory+ "/alg6"

        if not os.path.exists(newpath):
            os.makedirs(newpath)

        name_file =  data_path_directory+ "/alg6/cnk" + str(k) + "alg" + str(
            algo_arg) + "_" + str(h) + "epsilon_" + str(eps) + "iter_" + str(position) + ".csv"

        # np.savetxt(name_file, Expected_spread, delimiter=",")

        with open(name_file, 'wb') as f:
            np.savetxt(f, Expected_spread, fmt='%-7.8f', delimiter=',')

        return (1)


def m_zero(runs_alg,s_bulk,algo_arg_list,epsilon_values,list_k):

    for k in list_k:

        for iter in range(s_bulk):

            vals=[]
            for _ in range(10):

                indices = list(G_base.nodes())

                if k==4:
                    greedy_seeds = list(np.random.choice(indices, size=k, replace=False))
                    vals.append(i_x_s(greedy_seeds, matrices_X[iter]) * (n / m))

                elif k==8:

                    greedy_seeds = list(np.random.choice(indices, size=5, replace=False))
                    vals.append(i_x_s(greedy_seeds, matrices_X[iter]) * (n / m))

                elif k==12:

                    greedy_seeds = list(np.random.choice(indices, size=6, replace=False))
                    vals.append(i_x_s(greedy_seeds, matrices_X[iter]) * (n / m))

            for eps in epsilon_values:

                for algo_arg in algo_arg_list:

                    name_file = data_path_directory+"/alg" + str(algo_arg) + "/cnk" + str(k) + "alg" + str(
                        algo_arg) + "_" + str(0) + "epsilon_" + str(eps) + "iter_" + str(iter) + ".csv"

                    print("Here: " + name_file)

                    with open(name_file, 'wb') as f:
                        np.savetxt(f, vals, fmt='%-7.8f', delimiter=',')



def load_matrices_tilde(position, eps):

    file_name ='./Matrices_x_tilde/matrix_X_tilde' +str(position) + 'eps_' +str(eps)+'.pkl'

    with open(file_name, 'rb') as f:
        matrix_X = pickle.load(f)

    return matrix_X

def save_list_matrices_x_tilde(list_matrices_x_tilde):

    file_name = 'Matrices_X_tilde.pkl'

    with open(file_name, 'wb') as f:
        pickle.dump(list_matrices_x_tilde, f)


if __name__ == "__main__":

    save_data_x_tilde=False
    matrices_X = Generate_data.load_matrices_X()
    order_nodes=Generate_data.load_order_nodes()

    epsilon_values = [0.01,0.1,0.8]
    list_k = [4,8,12]
    m_values = [10,30,50,80,100,150,200]
    list_s_bulk = list(range(s_bulk))

    # algo_arg_list=[3,4,5,61,62,67,65,69]

    # m_zero(runs_alg, s_bulk, algo_arg_list, epsilon_values, list_k)

    # sys.exit(1)

    if save_data_x_tilde:

        list_s_bulk=list(range(s_bulk))
        arguments = list(product(list_s_bulk, epsilon_values))
        arguments = [list(i) for i in arguments]

        # generate_matrix_x_tilde(1, 0.01, m, n, matrices_X)


        new_args=[]
        for i in range(len(arguments)):
            temp=arguments[i]
            temp.append(m)
            temp.append(n)
            temp.append(matrices_X)
            new_args.append(temp)

        with multiprocessing.Pool(processes=number_CPU) as pool:
            pool.starmap(generate_matrix_x_tilde, new_args)

        pool.close()
        pool.join()

        dict_matrices_x_tilde={}

        for eps in epsilon_values:

            temp_list_matrix=[]

            for i in list_s_bulk:

                temp_list_matrix.append(load_matrices_tilde(i, eps))

            dict_matrices_x_tilde[eps]=temp_list_matrix


        save_list_matrices_x_tilde(dict_matrices_x_tilde)

    sys.exit(1)

    file_name = "Matrices_X_tilde.pkl"

    with open(file_name, 'rb') as f:
        dict_matrices_x_tilde = pickle.load(f)


    dict_set_of_matrices = {}

    for eps in epsilon_values:

        for k in list_k:

            rho = (np.exp(eps) + 1) ** -1
            dict_set_of_matrices[(k,eps,algo_arg)]=set_of_matrices(k, rho, algo_arg,penalty)


    arguments = list(product(list_k, epsilon_values, m_values,list_s_bulk))

    arguments = [list(i) for i in arguments]

    # generate_matrix_x_tilde(1, 0.01, m, n, matrices_X)

    processes = []

    new_args = []
    for i in range(len(arguments)):
        temp = arguments[i]
        temp.append(dict_set_of_matrices)
        temp.append(dict_matrices_x_tilde)
        temp.append(penalty)
        temp.append(matrices_X)
        temp.append(algo_arg)
        temp.append(order_nodes)
        new_args.append(temp)


    with multiprocessing.Pool(processes=number_CPU) as pool:
        pool.starmap(randomize_response, new_args)

    pool.close()
    pool.join()


