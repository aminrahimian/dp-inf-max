# expected spread size for different  DP algorithms

from dp_models import *

do_computation_exp_mech=True
do_computation_randomized=True
do_generate_matrix_x_randomized=False
do_computation_randomized_without_post_pr=True
do_computation_greedy_algortihm=False


parameters_greedy_alg={
            'matrix_iter_val': list(range(50)),
            'list_k' : [4,8,12],
            'save_computation':True,
            'number_CPU':4,
            'output_file_name':'greedy_alg.csv'}

parameters_exp_mech={
            'matrix_iter_val': list(range(50)),
            'm_values' : [0,5,10,20,30,40,60,80],
            'list_k' : [4,8,12],
            'epsilon_values' : [0.1,1,10],
            'runs_alg':10,
            'save_computation':True,
            'number_CPU':4,
            'output_file_name':'exp_mech.csv'}


parameters_randomized_version={
            'matrix_iter_val': list(range(50)),
            'm_values' :  [0,5,10,20,30,40,60,80],
            'list_k' : [4,8,12],
            'epsilon_values':[0.1,1,10],
            'runs_alg':10,
            'save_computation':True,
            'number_CPU':4,
            'output_file_name':'randomized_version.csv'}



def load_matrices_x():

    # Function to load pickle file with already generate X matrices
    file_name = 'matrices_x.pkl'

    with open(file_name, 'rb') as f:
        matrix_X = pickle.load(f)

    return matrix_X


def dump_expected_spread_exp_mech(list_matrices_x,matrix_iter_val, m_values,list_k,
                                  epsilon_values,runs_alg,save_computation,number_CPU
                                  ,output_file_name):

    arguments=product(matrix_iter_val, m_values,list_k, epsilon_values)
    arguments=[list(i) for i in arguments]
    new_args=[]

    copy_arguments=copy.deepcopy(arguments)

    for i in range(len(arguments)):
        temp=arguments[i]
        temp[0]=list_matrices_x[int(temp[0])]
        new_args.append(temp)

    new_args=new_args*runs_alg

    table_parameters=np.array(copy_arguments*runs_alg)

    with multiprocessing.Pool(processes=number_CPU) as pool:
        result=pool.starmap(expect_spread_exp_mechanism, new_args)

    pool.close()
    pool.join()

    if save_computation:

        df=pd.DataFrame(np.concatenate((table_parameters, np.array(result).reshape(table_parameters.shape[0],1)), axis=1),
                             columns=['arc_live','m','k','epsilon','ixs'])

        df.groupby(['m', 'k', 'epsilon'], as_index=False).agg(
            ixs_mu=('ixs', 'mean'),
            ixs_sd=('ixs', 'std'),
            ixs_n=('ixs', 'count')).to_csv(output_file_name, index=False)

def dump_expected_spread_randomized_response(list_matrices_x,matrix_iter_val, m_values,list_k,
                                  epsilon_values,runs_alg,save_computation,number_CPU,
                                  output_file_name):

    arguments = product(matrix_iter_val, m_values, list_k, epsilon_values)
    arguments = [list(i) for i in arguments]
    new_args = []

    dict_matrices_c = fill_matrix_c(list_k, epsilon_values)
    copy_arguments = copy.deepcopy(arguments)

    for i in range(len(arguments)):
        temp = arguments[i]
        # print("This is " + str((temp)))
        temp[0] = list_matrices_x[int(temp[0])]
        temp.append(dict_matrices_c)
        temp.append(runs_alg)
        new_args.append(temp)

    # new_args = new_args * runs_alg

    table_parameters = np.array(copy_arguments)

    with multiprocessing.Pool(processes=number_CPU) as pool:
        result = pool.starmap(expect_spread_randomized_resp, new_args)

    pool.close()
    pool.join()

    print("structure_Data " + str(np.array(result)))

    hp = []

    for t in range(runs_alg):
        hp.append(np.concatenate((table_parameters,
                                  (np.array(result)[:, t]).reshape(table_parameters.shape[0], 1)),
                                 axis=1))

    # print("hp :  " + str(np.array(hp).reshape((runs_alg * table_parameters.shape[0]), table_parameters.shape[1] + 1)))

    if save_computation:
        df = pd.DataFrame(np.array(hp).reshape((runs_alg * table_parameters.shape[0]), table_parameters.shape[1] + 1),
                          columns=['arc_live', 'm', 'k', 'epsilon', 'ixs'])

        df.groupby(['m', 'k', 'epsilon'], as_index=False).agg(
            ixs_mu=('ixs', 'mean'),
            ixs_sd=('ixs', 'std'),
            ixs_n=('ixs', 'count')).to_csv(output_file_name, index=False)

        print(str(df))


def dump_expected_spread_randomized_response_wpp(list_matrices_x,matrix_iter_val, m_values,list_k,
                                  epsilon_values,runs_alg,save_computation,number_CPU,
                                  output_file_name):

    arguments = product(matrix_iter_val, m_values, list_k, epsilon_values)
    arguments = [list(i) for i in arguments]
    new_args = []

    dict_matrices_c = fill_matrix_c(list_k, epsilon_values)
    copy_arguments = copy.deepcopy(arguments)

    for i in range(len(arguments)):
        temp = arguments[i]
        print("This is " + str((temp)))
        temp[0] = list_matrices_x[int(temp[0])]
        temp.append(dict_matrices_c)
        temp.append(runs_alg)
        new_args.append(temp)

    # new_args = new_args * runs_alg

    table_parameters = np.array(copy_arguments)

    with multiprocessing.Pool(processes=number_CPU) as pool:
        result = pool.starmap(expect_spread_randomized_resp_wpp, new_args)

    pool.close()
    pool.join()

    # print("structure_Data " + str(np.array(result)))

    hp = []

    for t in range(runs_alg):
        hp.append(np.concatenate((table_parameters,
                                  (np.array(result)[:, t]).reshape(table_parameters.shape[0], 1)),
                                 axis=1))

    # print("hp :  " + str(np.array(hp).reshape((runs_alg * table_parameters.shape[0]), table_parameters.shape[1] + 1)))

    if save_computation:
        df = pd.DataFrame(np.array(hp).reshape((runs_alg * table_parameters.shape[0]), table_parameters.shape[1] + 1),
                          columns=['arc_live', 'm', 'k', 'epsilon', 'ixs'])

        df.groupby(['m', 'k', 'epsilon'], as_index=False).agg(
            ixs_mu=('ixs', 'mean'),
            ixs_sd=('ixs', 'std'),
            ixs_n=('ixs', 'count')).to_csv(output_file_name, index=False)


def dump_expected_spread_greedy_alg(list_matrices_x,matrix_iter_val,
                                    list_k,save_computation,number_CPU
                                  ,output_file_name):

    arguments=product(matrix_iter_val,list_k)
    arguments=[list(i) for i in arguments]
    new_args=[]

    copy_arguments=copy.deepcopy(arguments)

    for i in range(len(arguments)):
        temp=arguments[i]
        temp[0]=list_matrices_x[int(temp[0])]
        new_args.append(temp)

    new_args=new_args

    table_parameters=np.array(copy_arguments)

    with multiprocessing.Pool(processes=number_CPU) as pool:
        result=pool.starmap(expect_spread_greedy_algorithm, new_args)

    pool.close()
    pool.join()

    if save_computation:

        df=pd.DataFrame(np.concatenate((table_parameters, np.array(result).reshape(table_parameters.shape[0],1)), axis=1),
                             columns=['arc_live','k','ixs'])

        df.groupby([ 'k', ], as_index=False).agg(
            ixs_mu=('ixs', 'mean'),
            ixs_sd=('ixs', 'std'),
            ixs_n=('ixs', 'count')).to_csv(output_file_name, index=False)


if __name__ == "__main__":

    list_matrices_x = load_matrices_x()

    if do_computation_exp_mech:

        dump_expected_spread_exp_mech(list_matrices_x,**parameters_exp_mech)

    if do_computation_randomized:

        dump_expected_spread_randomized_response(list_matrices_x, **parameters_randomized_version)

    if do_computation_randomized_without_post_pr:
        parameters_randomized_version['output_file_name']='randomized_version_wpp.csv'
        dump_expected_spread_randomized_response_wpp(list_matrices_x, **parameters_randomized_version)

    if do_computation_greedy_algortihm:
        dump_expected_spread_greedy_alg(list_matrices_x, **parameters_greedy_alg)

