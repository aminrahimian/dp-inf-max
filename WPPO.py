import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
from scipy.stats import bernoulli
from scipy.special import comb
import pickle
import gurobipy as gp
from gurobipy import GRB


#Parameters

n = 1005
m = 100
k = 4
p_ic = 0.05
s_bulk=40
# epsilon = 1e-1

G_base=nx.read_edgelist("email-Eu-core.txt.gz", nodetype=int, data=(("Type", str),))


def load_subgraphs():
    
    with open('live_edges.pkl', 'rb') as f:
        m_live_edges = pickle.load(f)
        
    return m_live_edges



# m_live_edges=load_subgraphs()


def load_matrix_X():
    
    
    with open('matrices_X.pkl', 'rb') as f:
        matrix_X = pickle.load(f)
        
    return matrix_X


matrices_X=load_matrix_X()



def i_x_s(S,x): 
    
    if S==[]:
    
        return 0
    
    else:
        
        sub_x= x[:,S] 
        temp=np.sum(sub_x, axis=1)
        
        return sum(temp>=1)
    



def generate_matrix_x_tilde(matrix_X, rho,m,n):
    
    
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
    
    
    return comb(b,l_curv)* comb((l-b), (a-b+l_curv))*(rho**(a-b + 2*l_curv))* ((1-rho)**(l-a+b -2*l_curv))



def sum_combinatorics(a,b,l,rho):
    
    
    lower= max(0, b-a)
    upper= min((l-a),b)
    
    
    cont=0
    
    for i in range(lower, upper+1):
        
        
        cont+=combinatorics(a,b,i,rho,l)
    

    return cont



def fill_matrix_c(matrix_C, dim_matrix_c, rho,l):
    
    for i in range(dim_matrix_c):
        
        for j in range(dim_matrix_c):
            
            matrix_C[i,j]=sum_combinatorics(i,j,l, rho)
            
            
    return matrix_C


def f_tilde_a(seed_set,m,a,matrix_x_tilde):
    
    cont=0
    
    for i in range(m):
        
        x_t_tilde=np.sum(matrix_x_tilde[i,seed_set])
        
        if x_t_tilde==a:
            
            cont+=1 
        
        else:
            
            pass
        
    
    return cont*(1/m)



def set_of_matrices(k,rho):
    
    
    set_of_matrices=[]
    
    
    for o in range(k+1):
        
        dim_matrix_c=o+1
        
        matrix_C=np.zeros((dim_matrix_c,dim_matrix_c))
        
        matrix_C_up=fill_matrix_c(matrix_C, dim_matrix_c, rho,o)
        
        matrix_C_up=np.append(matrix_C_up,[ dim_matrix_c*[1]], axis=0)

        psd_inverse=np.linalg.inv(np.matmul(matrix_C_up.T, matrix_C_up))
        
        set_of_matrices.append(np.matmul(psd_inverse,matrix_C_up.T))
        
        

    return set_of_matrices



def set_of_matrices_without_inversion(k,rho):
    
    
    set_of_matrices=[]
    
    
    for o in range(k+1):
        
        dim_matrix_c=o+1
        
        matrix_C=np.zeros((dim_matrix_c,dim_matrix_c))
        
        matrix_C_up=fill_matrix_c(matrix_C, dim_matrix_c, rho,o)
        
        # matrix_C_up=np.append(matrix_C_up, [dim_matrix_c*[1]], axis=0)

        
        set_of_matrices.append(np.linalg.inv(matrix_C_up))
    
    
    return set_of_matrices






def binary_converter(l, list_nodes):
    
    
    if l in list_nodes:
        
        return 1
    
    else:
        
        return 0


def f_tilded(seed_set, matrix_x_tilde,m):
    
    
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

def j_value(seed_set,m,matrix_x_tilde,set_of_matrices_in):
    
    
    f_tilde=f_tilded(seed_set, matrix_x_tilde,m)  
    
    
    
    f_0=f_tilde[0]
    

    return n*(1-f_0)




def j0_value(seed_set,m,matrix_x_tilde,set_of_matrices_in):
    
    
    f_tilde=f_tilded(seed_set, matrix_x_tilde,m)  
    
    # f_0=np.matmul(set_of_matrices_in[len(seed_set)], f_tilde)[0]
    
    # if f_0<=0:
        
    #     f_0=0
    
    return f_tilde[0]



def my_indices(lst, item):
   return [i for i, x in enumerate(lst) if x == item]



list_nodes=list(G_base.nodes())


m_values = [5,10,20,30,40,60,80]

epsilon_values = [1,10]

for eps in epsilon_values:
    
    rho=(np.exp(eps)+1)**-1
    
    set_of_matrices_array=set_of_matrices_without_inversion(k,rho)

    
    for h in m_values:
        
        
        Expected_spread=[]
    
        ## be sure you load bulk matrices
        
        c_matrix_C = copy.deepcopy(matrices_X)
        
        ## Iteration over influence samples 
        
        
        for z in range(40):
            
                
            sub_matrix_X = c_matrix_C[z][0:h, :]
            print("iteration  :" +str(z) + " m value" +str(h) +" eps " +str(eps))
            
            
            ## Iteration over runnings of the algorithm
           
            for chi in range(5):
                
               
                
                seed_set=[]
                
                matrix_x_tilde=generate_matrix_x_tilde(sub_matrix_X, rho,h,n)
                
                
                for w in range(k):
                    
                    # print("seed  :" +str(w))
                     
                    set_s= copy.deepcopy(seed_set)
                    
                    A=set(set_s)
                    B=set(G_base.nodes())
                    
                    iter_set=B.difference(A)
                    
                    set_s_union_v=copy.deepcopy(seed_set)
                    list_iter_set=list(iter_set)
                
                
                    j_lists=[]
                    
                    # j_o_vals=[]
                    
                    for i in range(len(list_iter_set)):
                        
                        
                        set_s_union_v.append(i)
                        
                        # print("Value of J_ms" +str(j_value(set_s_union_v,h,matrix_x_tilde,set_of_matrices_array)))
                        j_lists.append(j_value(set_s_union_v,h,matrix_x_tilde,set_of_matrices_array))
                        # j_o_vals.append(j0_value(set_s_union_v,h,matrix_x_tilde,set_of_matrices_array))
                        
                        
                        set_s_union_v=copy.deepcopy(seed_set)
                     
                      
                    
                    candidates=my_indices(j_lists, max(j_lists))    
                    
                    seed_set.append(random.choice(candidates))
                    
                    
              
                    
                Expected_spread.append(i_x_s(seed_set, c_matrix_C[z])*(n/m))
            
        
            
        Expected_spread
        
        name_file="k" +str(k) +"alg6_"+str(h) + "epsilon_"+str(eps)+".csv"
        
        np.savetxt(name_file, Expected_spread, delimiter=",")
    
        

