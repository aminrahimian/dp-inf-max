import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
from scipy.stats import bernoulli
from scipy.special import comb
import pickle
import Randomize_version
import generate_data

#Parameters


n = 1005
m = 100
k = 4
p_ic = 0.05
no_ind_cascades=40
algo_iter=5
# epsilon = 1e-1

G_base=nx.read_edgelist("email-Eu-core.txt.gz", nodetype=int, data=(("Type", str),))

matrices_X=Generate_data.load_matrix_X()

list_nodes=list(G_base.nodes())


m_values = [5,10,20,30,40,60,80]

epsilon_values = [1,10]

for eps in epsilon_values:
    
    rho=(np.exp(eps)+1)**-1
    
    set_of_matrices_array=Randomize_version.set_of_matrices(k,rho)

    
    for h in m_values:
        
        
        Expected_spread=[]
    
        ## be sure you load bulk matrices
        
        c_matrix_C = copy.deepcopy(matrices_X)
        
        ## Iteration over influence samples 
        
        
        for z in range(no_ind_cascades):
            

            sub_matrix_X = c_matrix_C[z][0:h, :]
            print("iteration  :" +str(z) + " m value" +str(h) +" eps " +str(eps))
            
            
            ## Iteration over runnings of the algorithm
           
            for chi in range(algo_iter):
                

                
                seed_set=[]
                
                matrix_x_tilde=Randomize_version.generate_matrix_x_tilde(sub_matrix_X, rho,h,n)
                
                
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
                        j_lists.append(Randomize_version.j0_value(set_s_union_v,h,matrix_x_tilde,set_of_matrices_array))
                        # j_o_vals.append(j0_value(set_s_union_v,h,matrix_x_tilde,set_of_matrices_array))
                        

                        set_s_union_v=copy.deepcopy(seed_set)
                     
                      
                    
                    candidates=Randomize_version.my_indices(j_lists, max(j_lists))
                    
                    seed_set.append(random.choice(candidates))
                    
                    
              
                    
                Expected_spread.append(Randomize_version.i_x_s(seed_set, c_matrix_C[z])*(n/m))
            
        
            
        Expected_spread
        
        name_file="k" +str(k) +"alg6_"+str(h) + "epsilon_"+str(eps)+".csv"
        
        np.savetxt(name_file, Expected_spread, delimiter=",")
    
        


