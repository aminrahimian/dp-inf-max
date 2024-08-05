# dp-inf-max

This project provides python scripts to simulate the following differentially private algorithms for influence maximization problem: exponential mechanism and randomized response. R scripts plot the results of the simulations.

## getting started

### software 

+ RStudio 2022.07.0+548 
+ Python 3.10

*Depending on the size of the dataset and parameters simulations, user might require access to high performance computing (HPC)*

## simulations execution:

### data preprocessing 

User should transform network dataset information to run the simulations:

+ **generate_data.py:** It receives the adjency matrix of the network data set as csv file and generates N independent matrices X, where each matrix have m (number of influence samples obtained from an independent network cascades realization) rows and n columns (number of nodes of the network). User can modify the following parameters:

- m : number of influence samples.
- p_ic : probability  independent cascade model (ICM).
- N : number of ICM realizations.

**predetermined datasets:**

+ email-Eu-core dataset : [source.](https://snap.stanford.edu/data/email-Eu-core.html)
+ soc-hamsterster_v2 dataset : [source.](https://networkrepository.com/soc-hamsterster.php)
+ Erdos-Renyi: synthetic dataset using Erdos-Renyi model with parameters ($n=200, p=0.15$)
+ MSM network: MSM(Men who have Sex with Men) network using temporal exponential random graph models (TERGMs) based on data from the egocentrically sampled survey ARTnet, conducted in the United States from 2017 to 2019, with 4,909 participants and 16,198 sexual partnerships. [source.](https://github.com/EpiModel/ARTnet, https://github.com/EpiModel/NetAnalysis-SF-ATL)

### execution:

+ **computing_mean_spread_size.py:** generates a dataset for each algorithm (setting to TRUE) given a list of parameters specified by the user (m ,k, epsilon):  

  + do_computation_exp_mech : simulation using exponential mechanism algorithm.
  + do_computation_randomized : simulation using randomized response algorithm.
  + do_computation_randomized_without_post_pr : simulation using randomized response algorithm skiping line 2  and using $\tilde{f}_0$ in line 3 of algorihm 5.
  + do_computation_greedy_algortihm: simulation using a non-private  algorithm. 
  + do_computation_greedy_algortihm_reference: simulation using a non-private  algorithm when $m \rightarrow \infty$.


+ **output structure**: 

| m  | k   | epsilon| Ixs_mu | Ixs_sd | Ixs_n |
| --- | ---- | --- | ---- | --- | ---- |
| 100 | 4 | 0.1| 120 | 20 | 10 |
| . | . | .| . | . | .|


-**computing_mean_spread_l2_regularization.py:** generates a dataset to evaluate  L2-regularization sweeping over different penalty values for a fixed m and epsilon. 

+ **output structure**: 

| k | penalty  | Ixs_mu | Ixs_sd | Ixs_n |
| --- | ---- | --- | ---- | --- |
| 4 | 0.1 | 120 | 20 | 10 |
| . | . | . | . | .|


+ **auxiliary script:**

  + **dp-models.py:** contains all dp algorithms and functions to generate all additional parameters needed in the simulations (e.g., matrix C used in algortihm 5)


### Visualization

+ **figure 1**: R script to generate a pdf file for figure 1 using outputs for computing_mean_spread_size.py setting *dataset_id = 'email-Eu-core'*. 

+ **figure 2**: R script to generate a pdf file for figure 2 using outputs for computing_mean_spread_size.py setting *dataset_id = 'erdos_renyi'*.
  
+ **figure 3**: R script to generate a pdf file for figure 3 using outputs for computing_mean_spread_size.py setting *dataset_id = 'MSM_network'*.



