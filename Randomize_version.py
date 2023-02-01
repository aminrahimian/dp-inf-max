import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
from scipy.stats import bernoulli
from scipy.special import comb
import pickle
import Generate_data


# Parameters of the network

# Parameters of College Network

n = 2426  # Number of nosde
m = 250  # Number of nodes samples to built matrix X for each influence cascade
k = 4  # Number of the seed set
p_ic = 0.05  # Probability for Independent Cascade Model
s_bulk = 2  # Number of Influence cascade model
# epsilon = 1
penalty = 0   # Penalty value for the regularization method


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
    

def generate_matrix_x_tilde(matrix_X, rho,m,n):
    
    #   Function to generate the matrix x_tilde from matrix X of dimensions mxn
    #   Rho is the flipping probability

    matrix_x_tilde=np.zeros((m,n))
    
    for i in range(m):
        
        row_i= matrix_X[i,:]
        
        flipped_row= 1- row_i
        
        for j in range(n):
            
            if bernoulli.rvs(p=rho, size=1)[0]==1:
                
                matrix_x_tilde[i,j]=flipped_row[j]
                
            else:
                
                matrix_x_tilde[i,j]=matrix_X[i,j]
        
    
    
    return matrix_x_tilde
    


# rho=(np.exp(epsilon)+1)**-1

# matrix_x_tilde=generate_matrix_x_tilde(matrix_X, rho,m,n)


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


def set_of_matrices(k, rho, algo_arg):

    # Function that returns all the inverse matrices X in function of the seed set k
    if not algo_arg=="5":

        set_of_matrices = []

        for o in range(k + 1):

            dim_matrix_c = o + 1
            matrix_C = np.zeros((dim_matrix_c, dim_matrix_c))
            matrix_C_up = fill_matrix_c(matrix_C, dim_matrix_c, rho, o)
            # matrix_C_up=np.append(matrix_C_up,[ dim_matrix_c*[1]], axis=0)
            psd_inverse = np.linalg.inv(matrix_C_up)
            set_of_matrices.append(psd_inverse)

        return set_of_matrices


    else:

        set_of_matrices = []

        for o in range(k + 1):

            dim_matrix_c = o + 1
            set_of_matrices.append(np.eye(dim_matrix_c))

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
    if not algo_arg=="6":

        f_tilde=f_tilded(seed_set, matrix_x_tilde,m)
        f_0=np.matmul(set_of_matrices_in[len(seed_set)], f_tilde)[0]

    else:

        f_tilde = f_tilded(seed_set, matrix_x_tilde, m)
        f_0 = np.matmul(set_of_matrices_in[len(seed_set)], f_tilde)[0]


    return n*(1-f_0)



def my_indices(lst, item):

    # Auxiliary function to randomize selection process when max function is tight.

   return [i for i, x in enumerate(lst) if x == item]


if __name__ == "__main__":

    matrices_X = Generate_data.load_matrices_X()

    algo_arg="4"

    epsilon_values = [0.01, 1]

    dict_matrices_x_tilde = {}

    for eps in epsilon_values:
        rho = (np.exp(eps) + 1) ** -1
        dict_matrices_x_tilde[eps] = [generate_matrix_x_tilde(i, rho, m, n) for i in matrices_X]

    list_nodes = list(G_base.nodes())

    # List that defines the number of rows of X used.

    m_values = [10, 30, 50, 80, 100, 150, 200]

    for eps in epsilon_values:

        rho = (np.exp(eps) + 1) ** -1

        set_of_matrices_array = set_of_matrices(k, rho, algo_arg, penalty)

        # Iteration over the size of the submatrix of X of size m.

        for h in m_values:

            Expected_spread = []

            # Iteration over the number of influence cascades
            # no_ind_cascades
            for z in range(s_bulk):

                matrix_x_tilde = dict_matrices_x_tilde[eps][z][0:h, :]
                print("iteration  :" + str(z) + " m value" + str(h) + " eps " + str(eps))

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

                        j_lists.append(j_value(set_s_union_v, h, matrix_x_tilde, set_of_matrices_array,algo_arg, penalty))

                        set_s_union_v = copy.deepcopy(seed_set)

                    candidates = my_indices(j_lists, max(j_lists))

                    seed_set.append(random.choice(candidates))

                Expected_spread.append(i_x_s(seed_set, matrices_X[z]) * (n / m))


            ######
            Expected_spread

            name_file =  "cnk" + str(k) + "algo"+ algo_arg  + "_"+str(h)+"epsilon_" + str(eps) + ".csv"

            np.savetxt(name_file, Expected_spread, delimiter=",")




