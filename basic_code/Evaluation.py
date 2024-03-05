import MachineLearning as ml
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import scipy as sp
import numpy as np
import pandas as pd
import random
from decimal import Decimal

def right_round(num,keep_n):
    if isinstance(num,float):
        num = str(num)
    return Decimal(num).quantize((Decimal('0.' + '0'*keep_n)))

def performance_estimate(predict_value_list, real_value_list):
    PCC = sp.stats.pearsonr(predict_value_list, real_value_list)[0]
    mse1 = mean_squared_error(predict_value_list, real_value_list)
    mse2 = pow(mse1, 0.5)
    RMSE = mse2
    R2 = r2_score(real_value_list, predict_value_list)
    performance_estimator = {'PCC': PCC, 'RMSE': RMSE, 'R^2': R2}
    return performance_estimator


###### ensemble functions
default_params = {'n_estimators': 4000, 'max_depth': 7, 'min_samples_split': 2,
                  'learning_rate': 0.01, 'loss': 'squared_error', 'max_features': 'sqrt', 'subsample': 0.7}


def comparison_models(path, character_table_list, label_csv, label_name, index_list_names, num_crossvalidation,
                      project_name, mode='random', machine_learning=ml.gradient_boosting, params=default_params, need_machine_learning=True):
    print(project_name)
    dataset_file = pd.read_csv(path + label_csv)
    label_list = dataset_file[label_name].tolist()
    num_elements = len(label_list)
    num_index_lists = len(index_list_names)
    ### Generate cross validation index list
    if mode == 'random':
        index_list_names_set = []
        for index in range(num_index_lists):
            random_num = random.sample(range(0, num_elements), num_elements)
            crossvalidation_list = []
            for i in range(num_crossvalidation - 1):
                crossvalidation_list.append(
                    random_num[
                    int(i * num_elements / num_crossvalidation):int((i + 1) * num_elements / num_crossvalidation)])
            crossvalidation_list.append(
                random_num[int((num_crossvalidation - 1) * num_elements / num_crossvalidation):])
            ### Record the cross validation index list
            ml.Listtotxt(crossvalidation_list,
                         path + 'IndexSet' + str(num_crossvalidation) + 'crossvalidation_' + project_name + str(
                             index) + '.txt')
            index_list_names_set.append(
                'IndexSet' + str(num_crossvalidation) + 'crossvalidation_' + project_name + str(index) + '.txt')
        cross_validation_list_names = index_list_names_set
    else:
        cross_validation_list_names = index_list_names

    ### proceed machine learning
    if need_machine_learning:
        for character_table in character_table_list:
            for i in range(len(cross_validation_list_names)):
                CV_list_name = cross_validation_list_names[i]
                file_name = project_name + character_table + str(i)
                ml.CrossValidationMachineLearning(path, 'character_table' + character_table + '.txt', label_csv, label_name,
                                                  num_crossvalidation, name=file_name, validation_name=CV_list_name,
                                                  params=params, model=machine_learning)

    ### Evaluation
    data_list = []
    data_std_list = []
    column = ['', 'PCC', 'RMSE', 'R^2']
    for character_table in character_table_list:
        PCC_list = []
        RMSE_list = []
        R2_list = []
        for i in range(len(cross_validation_list_names)):
            index_list = ml.TXTtoList(
                path + 'IndexSet' + str(num_crossvalidation) + 'crossvalidation_' + project_name + str(i) + '.txt')
            file_name = project_name + character_table + str(i)
            predict = np.loadtxt(path + 'Final_prediction_' + file_name + '.txt')
            PCC = 0
            RMSE = 0
            R2 = 0
            for indexes in index_list:
                part_predict_value = []
                part_real_value = []
                for index in indexes:
                    part_predict_value.append(predict[int(index)])
                    part_real_value.append(label_list[int(index)])
                result = performance_estimate(part_predict_value, part_real_value)
                PCC = PCC + result['PCC']
                RMSE = RMSE + result['RMSE']
                R2 = R2 + result['R^2']
            PCC_list.append(PCC / num_crossvalidation)
            RMSE_list.append(RMSE / num_crossvalidation)
            R2_list.append(R2 / num_crossvalidation)
        print(character_table + 'mean PCC=', np.mean(PCC_list))
        print(character_table + 'PCC std=', np.std(PCC_list))
        print(character_table + 'mean RMSE=', np.mean(RMSE_list))
        print(character_table + 'RMSE std=', np.std(RMSE_list))
        print(character_table + 'mean R^2=', np.mean(R2_list))
        print(character_table + 'R^2 std=', np.std(R2_list))

        row_elements = []
        row_elements.append(character_table)
        result_string = str(right_round(np.mean(PCC_list), 3))
        row_elements.append(result_string)
        result_string = str(right_round(np.mean(RMSE_list), 3))
        row_elements.append(result_string)
        result_string = str(right_round(np.mean(R2_list), 3))
        row_elements.append(result_string)

        data_list.append(row_elements)

        row_elements2 = []
        row_elements2.append(character_table)
        result_string = str(right_round(np.std(PCC_list), 3))
        row_elements2.append(result_string)
        result_string = str(right_round(np.std(RMSE_list), 3))
        row_elements2.append(result_string)
        result_string = str(right_round(np.std(R2_list), 3))
        row_elements2.append(result_string)
        
        data_std_list.append(row_elements2)

    df = pd.DataFrame(columns=column, data=data_list)
    df.to_csv(path + project_name + '_Performance(mean).csv')
    df = pd.DataFrame(columns=column, data=data_std_list)
    df.to_csv(path + project_name + '_Performance(std).csv')

def indexes_to_aug_indexes(aug_index_list):
    num_elements = len(aug_index_list)
    num_indexes = int(aug_index_list[-1]) + 1
    aug_index_set = []
    for i in range(num_indexes):
        aug_index_set.append([])

    for i in range(num_elements):
        aug_index_set[int(aug_index_list[i])].append(i)
    return aug_index_set

def comparison_models_aug(path, character_table_list, label_csv, label_name, index_list_names, num_crossvalidation,
                          project_name, mode='random', machine_learning=ml.gradient_boosting, params=default_params, need_machine_learning=True):
    print(project_name)
    dataset_file = pd.read_csv(path + label_csv)
    label_list = dataset_file[label_name].tolist()
    num_elements = len(label_list)
    num_index_lists = len(index_list_names)
    index_list = dataset_file['index'].tolist()
    num_index = index_list[-1] + 1

    aug_index_set = indexes_to_aug_indexes(index_list)

    ### Generate cross validation index list
    if mode == 'random':
        index_list_names_set = []
        for index in range(num_index_lists):
            random_num = random.sample(range(0, num_index), num_index)
            crossvalidation_list = []
            for i in range(num_crossvalidation - 1):
                crossvalidation_list.append(
                    random_num[int(i * num_index / num_crossvalidation):int((i + 1) * num_index / num_crossvalidation)])
            crossvalidation_list.append(random_num[int((num_crossvalidation - 1) * num_index / num_crossvalidation):])
            ### Record the cross validation index list
            ml.Listtotxt(crossvalidation_list,
                         path + 'IndexSet' + str(num_crossvalidation) + 'crossvalidation_' + project_name + str(
                             index) + '.txt')
            index_list_names_set.append(
                'IndexSet' + str(num_crossvalidation) + 'crossvalidation_' + project_name + str(index) + '.txt')
        cross_validation_list_names = index_list_names_set
    else:
        cross_validation_list_names = index_list_names

    ### proceed machine learning
    if need_machine_learning:
        for character_table in character_table_list:
            for i in range(len(cross_validation_list_names)):
                CV_list_name = cross_validation_list_names[i]
                file_name = project_name + character_table + str(i)
                ml.CrossValidationMachineLearning_aug(path, 'character_table_' + project_name + character_table + '.txt', label_csv, label_name, num_crossvalidation, name=file_name,validation_name=CV_list_name, params=params, model=machine_learning)

    ### Evaluations
    value_type = ['predict value of the original data','predict value of the mean of all augmented data','predict value of the median of all augmented data']
    column = ['', 'PCC', 'RMSE', 'R^2']
    data_list = []
    data_std_list = []
    for character_table in character_table_list:
        PCC_list = [[],[],[]]
        RMSE_list = [[],[],[]]
        R2_list = [[],[],[]]
        for i in range(len(cross_validation_list_names)):
            part_index_list = ml.TXTtoList(path + 'IndexSet' + str(
                num_crossvalidation) + 'crossvalidation_' + project_name + str(i) + '.txt')
            file_name = project_name + character_table + str(i)
            predict = np.loadtxt(path + 'Final_prediction_' + file_name + '.txt')
            PCC1 = 0
            PCC2 = 0
            PCC3 = 0
            RMSE1 = 0
            RMSE2 = 0
            RMSE3 = 0
            R21 = 0
            R22 = 0
            R23 = 0
            for indexes in part_index_list:
                part_predict_value1 = []
                part_predict_value2 = []
                part_predict_value3 = []
                part_real_value = []
                for index in indexes:
                    aug_indexes = aug_index_set[int(index)]
                    value1 = predict[aug_indexes[0]]
                    values = []
                    for aug in aug_indexes:
                        values.append(predict[aug])
                    value2 = np.mean(values)
                    value3 = np.median(values)

                    part_predict_value1.append(value1)
                    part_predict_value2.append(value2)
                    part_predict_value3.append(value3)
                    part_real_value.append(label_list[aug_indexes[0]])
                result1 = performance_estimate(part_predict_value1, part_real_value)
                result2 = performance_estimate(part_predict_value2, part_real_value)
                result3 = performance_estimate(part_predict_value3, part_real_value)
                PCC1 = PCC1 + result1['PCC']
                PCC2 = PCC2 + result2['PCC']
                PCC3 = PCC3 + result3['PCC']
                RMSE1 = RMSE1 + result1['RMSE']
                RMSE2 = RMSE2 + result2['RMSE']
                RMSE3 = RMSE3 + result3['RMSE']
                R21 = R21 + result1['R^2']
                R22 = R22 + result2['R^2']
                R23 = R23 + result3['R^2']
            PCC_list[0].append(PCC1 / num_crossvalidation)
            PCC_list[1].append(PCC2 / num_crossvalidation)
            PCC_list[2].append(PCC3 / num_crossvalidation)
            RMSE_list[0].append(RMSE1 / num_crossvalidation)
            RMSE_list[1].append(RMSE2 / num_crossvalidation)
            RMSE_list[2].append(RMSE3 / num_crossvalidation)
            R2_list[0].append(R21 / num_crossvalidation)
            R2_list[1].append(R22 / num_crossvalidation)
            R2_list[2].append(R23 / num_crossvalidation)
        for i in range(3):
            print('Use '+value_type[i]+' as final predict value for the original data')
            print(character_table + 'mean PCC=', np.mean(PCC_list[i]))
            print(character_table + 'PCC std=', np.std(PCC_list[i]))
            print(character_table + 'mean RMSE=', np.mean(RMSE_list[i]))
            print(character_table + 'RMSE std=', np.std(RMSE_list[i]))
            print(character_table + 'mean R^2=', np.mean(R2_list[i]))
            print(character_table + 'R^2 std=', np.std(R2_list[i]))

        row_elements = []
        row_elements.append(character_table)
        result_string = ''
        for j in range(2):
            result_string = result_string + str(right_round(np.mean(PCC_list[j]),3)) + '/'
        result_string = result_string + str(right_round(np.mean(PCC_list[2]),3))
        row_elements.append(result_string)
        result_string = ''
        for j in range(2):
            result_string = result_string + str(right_round(np.mean(RMSE_list[j]),2)) + '/'
        result_string = result_string + str(right_round(np.mean(RMSE_list[2]),2))
        row_elements.append(result_string)
        result_string = ''
        for j in range(2):
            result_string = result_string + str(right_round(np.mean(R2_list[j]),2)) + '/'
        result_string = result_string + str(right_round(np.mean(R2_list[2]),2))
        row_elements.append(result_string)
        data_list.append(row_elements)

        row_elements2 = []
        row_elements2.append(character_table)
        result_string = ''
        for j in range(2):
            result_string = result_string + str(right_round(np.std(PCC_list[j]), 3)) + '/'
        result_string = result_string + str(right_round(np.std(PCC_list[2]), 3))
        row_elements2.append(result_string)
        result_string = ''
        for j in range(2):
            result_string = result_string + str(right_round(np.std(RMSE_list[j]), 3)) + '/'
        result_string = result_string + str(right_round(np.std(RMSE_list[2]), 3))
        row_elements2.append(result_string)
        result_string = ''
        for j in range(2):
            result_string = result_string + str(right_round(np.std(R2_list[j]), 3)) + '/'
        result_string = result_string + str(right_round(np.std(R2_list[2]), 3))
        row_elements2.append(result_string)
        data_std_list.append(row_elements2)

    df = pd.DataFrame(columns=column,data= data_list)
    df.to_csv(path + project_name + '_Performance(mean).csv')
    df = pd.DataFrame(columns=column, data=data_std_list)
    df.to_csv(path + project_name + '_Performance(std).csv')

