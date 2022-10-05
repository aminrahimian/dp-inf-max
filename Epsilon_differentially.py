import networkx as nx
import numpy as np
import pandas as pd
import random
import copy
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



def i_x_s(S,x): 
    
    if S==[]:
    
        return 0
    
    else:
        
        sub_x= x[:,S] 
        temp=np.sum(sub_x, axis=1)
        
        return sum(temp>=1)
    


seed_set=[]


A=set(m_live_edges[0].nodes())


def exponential_probabilities(Ix, epsilon,n,m):
    
    return np.exp(((epsilon*n)/(m*2))*Ix)


Expected_spread=[]



###chose a matrix according to the size m that we want to plot

# sub_matrix_X= matrix_X[0:10,:]



for w in range(5):
    
    seed_set=[]
        
    for i in range(k):
        
        
        Q=copy.deepcopy(A)
        B=set(seed_set)
        S_union_v = copy.deepcopy(seed_set) 
        seed_set_s=copy.deepcopy(seed_set)
        
      
        iter_set= Q.difference(B)
        list_S_minus_seed=list(iter_set)
        table=np.zeros((len(iter_set),3))
        
        
        for j in range(len(iter_set)):
            
            S_union_v.append(list_S_minus_seed[j])
            table[j,1]=i_x_s(S_union_v, matrix_X) - i_x_s(seed_set_s, matrix_X) 
            table[j,0]=list_S_minus_seed[j]
                 
            S_union_v=copy.deepcopy(seed_set)
            
          
            
        #transform to probability 
        
        
        probabilities=np.array([exponential_probabilities(table[k,1], epsilon,n,m) for k in range(len(iter_set))])
        table[:,2]=probabilities/np.sum(probabilities)
        
    
    
        seed_set.append(int(np.random.choice(list(table[:,0]), 1, p=list(table[:,2]))[0]))
        
        
        
        
    Expected_spread.append(i_x_s(seed_set, matrix_X)*(n/m))





Expected_spread








































