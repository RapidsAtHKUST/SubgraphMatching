import sys
import os
import glob
from subprocess import Popen, PIPE


def generate_args(binary, *params):
    arguments = [binary]
    arguments.extend(list(params))
    return arguments


def execute_binary(args):
    process = Popen(' '.join(args), shell=True, stdout=PIPE, stderr=PIPE)
    (std_output, std_error) = process.communicate()
    process.wait()
    rc = process.returncode

    return rc, std_output, std_error


def check_correctness(binary_path, data_graph_path, query_folder_path, expected_results):
    # find all query graphs.
    query_graph_path_list = glob.glob('{0}/*.graph'.format(query_folder_path))

    # check the correctness of Filtering.
    filter_type_list = ['LDF', 'NLF', 'GQL', 'TSO', 'CFL', 'DPiso']

    for filter_type in filter_type_list:
        for query_graph_path in query_graph_path_list:
            execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', filter_type,
                                       '-order', 'GQL', '-engine', 'LFTJ')

            (rc, std_output, std_error) = execute_binary(execution_args)
            query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
            expected_embedding_num = expected_results[query_graph_name]

            if rc == 0:
                embedding_num = 0
                std_output_list = std_output.split('\n')
                for line in std_output_list:
                    if '#Embeddings' in line:
                        embedding_num = int(line.split(':')[1].strip())
                        break
                if embedding_num != expected_embedding_num:
                    print 'Filter type {0} {1} is wrong. Expected {2}, Output {3}'.format(filter_type, query_graph_name, expected_embedding_num,
                                                                      embedding_num)
                    exit(-1)
            else:
                print 'Filter type {0} {1} error.'.format(filter_type, query_graph_name)
                exit(-1)

        print 'Filter type {0} passes the correctness check.'.format(filter_type)

    print 'All filter types pass the correctness check.'

    # check the correctness of ordering.
    order_type_list = ['QSI', 'GQL', 'TSO', 'CFL', 'DPiso', 'RI', 'VF2PP']
    for order_type in order_type_list:
        for query_graph_path in query_graph_path_list:
            execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', 'LDF',
                                       '-order', order_type, '-engine', 'LFTJ')

            (rc, std_output, std_error) = execute_binary(execution_args)
            query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
            expected_embedding_num = expected_results[query_graph_name]

            if rc == 0:
                embedding_num = 0
                std_output_list = std_output.split('\n')
                for line in std_output_list:
                    if '#Embeddings' in line:
                        embedding_num = int(line.split(':')[1].strip())
                        break
                if embedding_num != expected_embedding_num:
                    print 'Order type {0} {1} is wrong. Expected {2}, Output {3}'.format(order_type, query_graph_name,
                            expected_embedding_num, embedding_num)
                    exit(-1)
            else:
                print 'Order type {0} {1} error.'.format(order_type, query_graph_name)
                exit(-1)

        print 'Order type {0} passes the correctness check.'.format(order_type)

    print 'All order types pass the correctness check.'

    # check the correctness of LFTJ engine
    for query_graph_path in query_graph_path_list:
        execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', 'LDF',
                                       '-order', 'GQL', '-engine', 'LFTJ')

        (rc, std_output, std_error) = execute_binary(execution_args)
        query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
        expected_embedding_num = expected_results[query_graph_name]

        if rc == 0:
            embedding_num = 0
            std_output_list = std_output.split('\n')
            for line in std_output_list:
                if '#Embeddings' in line:
                    embedding_num = int(line.split(':')[1].strip())
                    break
            if embedding_num != expected_embedding_num:
                print 'LFTJ engine {0} is wrong. Expected {1}, Output {2}'.format(query_graph_name,
                        expected_embedding_num, embedding_num)
                exit(-1)
        else:
            print 'LFTJ engine {0} error.'.format(query_graph_name)
            exit(-1)
    print 'LFTJ engine pass the correctness check.'

    # check the correctness of EXPLORE engine
    for query_graph_path in query_graph_path_list:
        execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', 'LDF',
                                       '-order', 'GQL', '-engine', 'EXPLORE')

        (rc, std_output, std_error) = execute_binary(execution_args)
        query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
        expected_embedding_num = expected_results[query_graph_name]

        if rc == 0:
            embedding_num = 0
            std_output_list = std_output.split('\n')
            for line in std_output_list:
                if '#Embeddings' in line:
                    embedding_num = int(line.split(':')[1].strip())
                    break
            if embedding_num != expected_embedding_num:
                print 'EXPLORE engine {0} is wrong. Expected {1}, Output {2}'.format(query_graph_name,
                        expected_embedding_num, embedding_num)
                exit(-1)
        else:
            print 'EXPLORE engine {0} error.'.format(query_graph_name)
            exit(-1)
    print 'EXPLORE engine pass the correctness check.'

    # check the correctness of GraphQL engine
    for query_graph_path in query_graph_path_list:
        execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', 'GQL',
                                       '-order', 'GQL', '-engine', 'GQL')

        (rc, std_output, std_error) = execute_binary(execution_args)
        query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
        expected_embedding_num = expected_results[query_graph_name]

        if rc == 0:
            embedding_num = 0
            std_output_list = std_output.split('\n')
            for line in std_output_list:
                if '#Embeddings' in line:
                    embedding_num = int(line.split(':')[1].strip())
                    break
            if embedding_num != expected_embedding_num:
                print 'GQL engine {0} is wrong. Expected {1}, Output {2}'.format(query_graph_name,
                        expected_embedding_num, embedding_num)
                exit(-1)
        else:
            print 'GQL engine {0} error.'.format(query_graph_name)
            exit(-1)
    print 'GQL engine pass the correctness check.'

    # check the correctness of QSI engine
    for query_graph_path in query_graph_path_list:
        execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', 'LDF',
                                       '-order', 'QSI', '-engine', 'QSI')

        (rc, std_output, std_error) = execute_binary(execution_args)
        query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
        expected_embedding_num = expected_results[query_graph_name]

        if rc == 0:
            embedding_num = 0
            std_output_list = std_output.split('\n')
            for line in std_output_list:
                if '#Embeddings' in line:
                    embedding_num = int(line.split(':')[1].strip())
                    break
            if embedding_num != expected_embedding_num:
                print 'QSI engine {0} is wrong. Expected {1}, Output {2}'.format(query_graph_name,
                        expected_embedding_num, embedding_num)
                exit(-1)
        else:
            print 'QSI engine {0} error.'.format(query_graph_name)
            exit(-1)
    print 'QSI engine pass the correctness check.'

    # check the correctness of DPiso engine.
    for query_graph_path in query_graph_path_list:
        execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', 'DPiso',
                                       '-order', 'DPiso', '-engine', 'DPiso')

        (rc, std_output, std_error) = execute_binary(execution_args)
        query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
        expected_embedding_num = expected_results[query_graph_name]

        if rc == 0:
            embedding_num = 0
            std_output_list = std_output.split('\n')
            for line in std_output_list:
                if '#Embeddings' in line:
                    embedding_num = int(line.split(':')[1].strip())
                    break
            if embedding_num != expected_embedding_num:
                print 'DPiso engine {0} is wrong. Expected {1}, Output {2}'.format(query_graph_name,
                        expected_embedding_num, embedding_num)
                exit(-1)
        else:
            print 'DPiso engine {0} error.'.format(query_graph_name)
            exit(-1)
    print 'DPiso engine pass the correctness check.'

    # check the correctness of CECI engine.
    for query_graph_path in query_graph_path_list:
        execution_args = generate_args(binary_path, '-d', data_graph_path, '-q', query_graph_path, '-filter', 'CECI',
                                       '-order', 'CECI', '-engine', 'CECI')

        (rc, std_output, std_error) = execute_binary(execution_args)
        query_graph_name = os.path.splitext(os.path.basename(query_graph_path))[0]
        expected_embedding_num = expected_results[query_graph_name]

        if rc == 0:
            embedding_num = 0
            std_output_list = std_output.split('\n')
            for line in std_output_list:
                if '#Embeddings' in line:
                    embedding_num = int(line.split(':')[1].strip())
                    break
            if embedding_num != expected_embedding_num:
                print 'CECI engine {0} is wrong. Expected {1}, Output {2}'.format(query_graph_name,
                        expected_embedding_num, embedding_num)
                exit(-1)
        else:
            print 'CECI engine {0} error.'.format(query_graph_name)
            exit(-1)
    print 'CECI engine pass the correctness check.'


if __name__ == '__main__':
    input_binary_path = sys.argv[1]
    if not os.path.isfile(input_binary_path):
        print 'The binary {0} does not exist.'.format(input_binary_path)
        exit(-1)

    # load expected results.
    input_expected_results = {}
    input_expected_results_file = 'expected_output.res'
    with open(input_expected_results_file, 'r') as f:
        for line in f:
            if line:
                result_item = line.split(':')
                input_expected_results[result_item[0].strip()] = int(result_item[1].strip())

    dir_path = os.path.dirname(os.path.realpath(__file__))
    input_data_graph_path = '{0}/data_graph/HPRD.graph'.format(dir_path)
    input_query_graph_folder_path = '{0}/query_graph/'.format(dir_path)

    check_correctness(input_binary_path, input_data_graph_path, input_query_graph_folder_path, input_expected_results)

