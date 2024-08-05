#compute mean expected spread size using l2-regularization for a given m and epsilonfrom dp_models import *#dataset_id='soc-hamsterster_v2'dataset_id = 'MSM_network'parameters_randomized_version = {    'matrix_iter_val': list(range(50)),    'm_values': [500],    'list_k': [10,20],    'epsilon_values': [5],    'penalty_values':[0],    'runs_alg': 10,    'save_computation': True,    'number_CPU': 4,    'output_file_name': 'randomized_response_l2.csv'}def load_matrices_x():    # Function to load pickle file with already generate X matrices    file_name = 'matrices_x.pkl'    with open(file_name, 'rb') as f:        matrix_X = pickle.load(f)    return matrix_Xdef dump_expected_spread_rr_l2_regularization(list_matrices_x,matrix_iter_val, m_values,list_k,                                  epsilon_values,penalty_values,runs_alg,save_computation,number_CPU,                                  output_file_name):    arguments = product(matrix_iter_val, m_values, list_k, epsilon_values,penalty_values)    arguments = [list(i) for i in arguments]    new_args = []    dict_matrices_c =fill_matrix_c_regularization(list_k, epsilon_values[0], penalty_values)    copy_arguments = copy.deepcopy(arguments)    for i in range(len(arguments)):        temp = arguments[i]        # print("This is " + str((temp)))        temp[0] = list_matrices_x[int(temp[0])]        temp.append(dict_matrices_c)        temp.append(runs_alg)        new_args.append(temp)    # new_args = new_args * runs_alg    table_parameters = np.array(copy_arguments)    with multiprocessing.Pool(processes=number_CPU) as pool:        result = pool.starmap(compute_exp_spread_randomized_l2_regularization, new_args)    pool.close()    pool.join()    print("structure_Data " + str(np.array(result)))    hp = []    for t in range(runs_alg):        hp.append(np.concatenate((table_parameters,                                  (np.array(result)[:, t]).reshape(table_parameters.shape[0], 1)),                                 axis=1))    # print("hp :  " + str(np.array(hp).reshape((runs_alg * table_parameters.shape[0]), table_parameters.shape[1] + 1)))    if save_computation:        df = pd.DataFrame(np.array(hp).reshape((runs_alg * table_parameters.shape[0]), table_parameters.shape[1] + 1),                          columns=['arc_live', 'm', 'k', 'epsilon', 'penalty','ixs'])        df.groupby(['m', 'k', 'penalty']).agg(            ixs_mu=('ixs', 'mean'),            ixs_sd=('ixs', 'std'),            ixs_n=('ixs', 'count')).to_csv(output_file_name)        print(str(df))if __name__ == "__main__":    list_matrices_x = load_matrices_x()    dump_expected_spread_rr_l2_regularization(list_matrices_x, **parameters_randomized_version)