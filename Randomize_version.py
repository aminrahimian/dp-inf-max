import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
from scipy.stats import bernoulli
from math import comb
import pickle


#Parameters


      
n=4039
m=1000
k=10
p_ic=0.1
epsilon=1e-1





G_base=nx.read_edgelist("facebook_combined.txt.gz", nodetype=int, data=(("Type", str),))
G_base.edges()


def save_graph_in_list(m,G_base):
    
    

    test_keys = list(G_base.nodes())
    test_values = list(range(len(test_keys)))
    dict_to_map = {test_keys[i]: test_values[i] for i in range(len(test_keys))}
      
    H = nx.relabel_nodes(G_base, dict_to_map)
    H=H.to_undirected()
    
    
    
    n=len(H.nodes())
    n_e=len(H.edges())
    
    
    
     ### n is the number of edges not nodes
    def live_edges_saving(n,p,G):
        
        k=int(n*p)
        sampled_edges = random.sample(G.edges, k)
        ini_nodes=set(G.nodes())
        sG=G.edge_subgraph(sampled_edges).copy()
        s_nodes=set(sG.nodes())
        diff_nodes=list(ini_nodes.difference(s_nodes))
        
        for i in diff_nodes:
            sG.add_node(i)
    
        
        # path='Live_edges\live_edges'+ str(step)+'_'+str(sample_size)+'.pkl'
        # pickle.dump(sG, open(path, 'wb'))
        
        return sG
    
    m_live_edges=[live_edges_saving(n_e,p_ic,H)  for i in range(m)]
    
    
    with open('live_edges.pkl', 'wb') as f:
        pickle.dump(m_live_edges, f)








def load_subgraphs():
    
    
    with open('live_edges.pkl', 'rb') as f:
        m_live_edges = pickle.load(f)
        
    return m_live_edges




def built_matrix_X(m,n,m_live_edges):
    
    
    matrix_X=np.zeros((m,n))
    indices=list(m_live_edges[0].nodes())
    
    uniform_nodes=np.random.choice(indices, size=m, replace=True)
    
    
    for i in range(m):
    
        matrix_X[i,:]=np.array([nx.has_path(m_live_edges[i], j, uniform_nodes[i]) for j in range(n)])
    

    return matrix_X




m_live_edges=load_subgraphs()

matrix_X=built_matrix_X(m,n,m_live_edges)





seed_set=[]




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
    
    
    


rho=(np.exp(epsilon)+1)**-1

matrix_x_tilde=generate_matrix_x_tilde(matrix_X, rho,m,n)





def combinatorics(a,b,l,rho):
    
    
    return comb(b,l)* comb((l-b), (a-b+l))* (rho**(a-b + 2*l))* (rho **(l-a+b -2*l))




def sum_combinatorics(a,b,S,rho):
    
    
    lower= max(0, b-a)
    upper= min((len(S) - a),b)
    
    
    cont=0
    
    for i in range(lower, upper+1):
        
        
        cont+=combinatorics(a,b,i,rho)
    
    

    return cont



def fill_matrix_c(matrix_C, dim_matrix_c, rho,seed_set):
    
    for i in range(dim_matrix_c):
        
        for j in range(dim_matrix_c):
            
            matrix_C[i,j]=sum_combinatorics(i,j,seed_set, rho)
            
            
            
            


def f_tilde_a(seed_set,m,a):
    
    cont=0
    
    for i in range(m):
        
        x_t_tilde=np.sum(matrix_x_tilde[i,:][seed_set])
        
        if x_t_tilde==a:
            
            cont+=1 
        
        else:
            
            pass
        
    
    return cont*(1/m)







def j_value(seed_set,n,m):
    
    dim_matrix_c=len(seed_set)+1
    
    matrix_C=np.zeros((dim_matrix_c,dim_matrix_c))
    
    
    fill_matrix_c(matrix_C, dim_matrix_c, rho,seed_set)
    
    f_tilde=[]
    

    
    
    for i in range(dim_matrix_c):
        
        f_tilde.append(f_tilde_a(seed_set,m,i))

    
    
    
    f_tilde=np.array(f_tilde)

    return n*(1-np.matmul(np.linalg.inv(matrix_C),f_tilde)[0])




Expected_spread=[]



for z in range(10):

    
    seed_set=[]
    
    
    
    for w in range(k):
        
        set_s= copy.deepcopy(seed_set)
        
        A=set(set_s)
        B=set(m_live_edges[0].nodes())
        
        iter_set=B.difference(A)
        
        set_s_union_v=copy.deepcopy(seed_set)
        list_iter_set=list(iter_set)
    
    
        j_lists=[]
        
        for i in range(len(list_iter_set)):
            
            set_s_union_v.append(i)
            
        
            j_lists.append(j_value(set_s_union_v,n,m) - j_value(seed_set,n,m))
            
         
          
         
        seed_set.append(list_iter_set[np.argmax(j_lists)])
        
        
  
        
    Expected_spread.append(i_x_s(seed_set, matrix_X)*(n/m))





Expected_spread

























