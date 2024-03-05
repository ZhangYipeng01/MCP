# Author: Zhang Yipeng
# Email: yipeng001@e.ntu.edu.sg
# =================================================
import gudhi as gd
import numpy as np
import math
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


####### Basic file-dealing functions
def TXTtoList(filepath, return_type='str'):
    fp = open(filepath, 'r', encoding='utf-8')
    content = fp.read()
    fp.close()
    rowlist = content.splitlines()
    recordlist = [row.split() for row in rowlist if row != 'END' if row.strip()]
    if return_type == 'str':
        return recordlist
    if return_type == 'int':
        int_recordlist = [[int(element) for element in inner_list] for inner_list in recordlist]
        return int_recordlist
    if return_type == 'float':
        float_recordlist = [[float(element) for element in inner_list] for inner_list in recordlist]
        return float_recordlist


def Listtotxt(List, filepath):
    file = open(filepath, 'w')
    for row in List:
        line = ' '.join(str(element) for element in row) + '\n'
        file.write(line)
    file.close()

def DistributionHist(path,data,name,num_bins='auto',qt=[1,99]):
    sns.histplot(data, kde=True, color='steelblue', bins=num_bins, edgecolor='black')
    mean_value = np.mean(data)
    q1 = np.percentile(data, qt[0])
    q2 = np.percentile(data, qt[1])
    plt.axvline(x=mean_value, color='green', linestyle='--', label='Mean')
    plt.axvline(x=q1, color='blue', linestyle='--', label=str(qt[0])+' Percentile')
    plt.axvline(x=q2, color='purple', linestyle='--', label=str(qt[1])+' Percentile')
    plt.title('Distribution of ' + name)
    plt.legend()
    plt.savefig(path + 'Distribution of ' + name + '.jpg')
    plt.close()

######## TDA basic functions
def AtomCombinationCondition(atom_type, comb):
    if comb[0] == 'all':
        return True
    else:
        if comb[0] == 'no':
            return atom_type not in comb
        else:
            return atom_type in comb


def SeperateAtomCombination(filepath, savepath, filename, atom_combination):
    List_reader = TXTtoList(filepath + filename + '.txt')
    for comb in atom_combination:
        CombinationCoords = []
        for row in List_reader:
            if AtomCombinationCondition(row[3], comb):
                CombinationCoords.append([float(row[0]), float(row[1]), float(row[2])])
        comb_string = ''
        for atom in comb:
            comb_string = comb_string + atom
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        Listtotxt(CombinationCoords, savepath + filename + '_' + comb_string + '.txt')


def RipsFiltration(coords, dim=2,
                   max_range=8):  # coords is the coordinates like [[0.5,3,2],[1.3,4,4],...]; dim is the max homology dimension
    rips_complex = gd.RipsComplex(points=coords, max_edge_length=max_range)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=dim)
    result = []
    for filtration_value in simplex_tree.get_filtration():
        result.append(filtration_value)
    return result


def RipsBarcode(coords, dim=2, max_range=8):
    rips_complex = gd.RipsComplex(points=coords, max_edge_length=max_range)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=dim + 1)
    barcodes = simplex_tree.persistence()
    betti = []
    for i in range(dim + 1):
        betti.append([])
    for bar in barcodes:
        betti[bar[0]].append([bar[1][0], bar[1][1]])
    return betti


def AlphaFiltration(coords):
    alpha_complex = gd.AlphaComplex(points=coords)
    simplex_tree = alpha_complex.create_simplex_tree()
    result = []
    for filtration_value in simplex_tree.get_filtration():
        result.append(filtration_value)
    return result


def AlphaBarcode(coords, dim=2):
    alpha_complex = gd.AlphaComplex(points=coords)
    simplex_tree = alpha_complex.create_simplex_tree()
    barcodes = simplex_tree.persistence()
    betti = []
    for i in range(dim + 1):
        betti.append([])
    for bar in barcodes:
        if bar[0] <= dim:
            betti[bar[0]].append([bar[1][0], bar[1][1]])
    return betti


def SimplicialFiltrationtoBarcode(filtration, dim=2):
    st = gd.SimplexTree()
    for simplex in filtration:
        st.insert(simplex[0], simplex[1])
    barcodes = st.persistence()
    betti = []
    for i in range(dim + 1):
        betti.append([])
    for bar in barcodes:
        if bar[0] <= dim:
            betti[bar[0]].append([bar[1][0], bar[1][1]])
    return betti


def ListStatistics(input_list):
    result = []
    ## total number
    total_num = len(input_list)
    # result.append(total_num)
    ## max
    if total_num == 0:
        result.append(0)
    else:
        result.append(np.max(input_list))
    ## min
    if total_num == 0:
        result.append(0)
    else:
        result.append(np.min(input_list))
    ## mean
    if total_num == 0:
        result.append(0)
    else:
        result.append(np.mean(input_list))
    ## sum
    if total_num == 0:
        result.append(0)
    else:
        result.append(np.sum(input_list))

    ## std
    if total_num <= 1:
        result.append(0)
    else:
        result.append(np.std(input_list))

    ##num_elememt
    result.append(total_num)

    return result


def BarcodetoStatistics(barcode, bar_range=[0, 8], style='Normal', function=ListStatistics):
    birth = []
    death = []
    persistence = []
    for bar in barcode:
        if bar[0] >= bar_range[0] and bar[1] <= bar_range[1]:
            birth.append(bar[0])
            death.append(bar[1])
            persistence.append(bar[1] - bar[0])
        feature_vector = ['sorry, no such style']
    if style == 'Normal':
        feature_vector = function(birth) + function(death) + function(persistence)
    if style == 'Endpoints':
        feature_vector = function(birth+death)
    return feature_vector


def BarcodetoStatistics_bins(barcode, num_bins, bar_range=[0, 8], style='Normal', function=ListStatistics):
    birth = []
    death = []
    persistence = []
    for i in range(num_bins):
        birth.append([])
        death.append([])
        persistence.append([])

    for bar in barcode:
        if bar[0] >= bar_range[0] and bar[1] < bar_range[1]:
            birth[math.floor(num_bins * (bar[0] - bar_range[0]) / (bar_range[1] - bar_range[0]))].append(bar[0])
            death[math.floor(num_bins * (bar[1] - bar_range[0]) / (bar_range[1] - bar_range[0]))].append(bar[1])
            persistence[math.floor(num_bins * (bar[1] - bar[0]) / (bar_range[1] - bar_range[0]))].append(bar[1] - bar[0])
    feature_vector = ['No such style']
    if style == 'Normal':
        feature_vector = []
        for i in range(num_bins):
            feature_vector = feature_vector + function(birth[i])
            feature_vector = feature_vector + function(death[i])
            feature_vector = feature_vector + function(persistence[i])

    if style == 'Endpoints':
        feature_vector = []
        for i in range(num_bins):
            feature_vector = feature_vector + function(birth[i]+death[i])

    return feature_vector

def Barcodetobetticurve(barcode,barcode_range,num_grid):
    betticurve = np.zeros(num_grid)
    for bar in barcode:
        for i in range(num_grid):
            if (bar[1]>barcode_range[0]+(barcode_range[1]-barcode_range[0])*i/num_grid) and (bar[0]<barcode_range[0]+(barcode_range[1]-barcode_range[0])*(i+1)/num_grid):
                betticurve[i] += 1
    return betticurve.tolist()


######## Ensemble functions: Dealing with all files in a file folder
def AtomCombination(filepath, savepath, atom_combination):
    item_list = os.listdir(filepath)
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    for item in item_list:
        item_path = os.path.join(filepath, item)
        if os.path.isfile(item_path) and item.endswith(".txt"):
            file_name = os.path.splitext(item)[0]
            SeperateAtomCombination(filepath, savepath, file_name, atom_combination)

def PreEstimate_datasize(path_data,save_path,project_name,num_bins='auto',qt=[1,99],key_word = ''):
    name_list = []
    num_coords_list = []
    item_list = os.listdir(path_data)
    for item in item_list:
        item_path = os.path.join(path_data, item)
        if item.endswith(".txt") and key_word in item:
            num_coords_list.append(len(TXTtoList(item_path)))
            name_list.append(item_path)
    DistributionHist(save_path,num_coords_list,'number of elements('+project_name+')',num_bins=num_bins,qt=qt)
    min_num = np.percentile(num_coords_list,0)
    mean_num = np.percentile(num_coords_list,50)
    max_num = np.percentile(num_coords_list,100)
    index_min = np.searchsorted(num_coords_list, min_num)
    index_mean = np.searchsorted(num_coords_list, mean_num)
    index_max = np.searchsorted(num_coords_list, max_num)

    return [min_num,mean_num,max_num,name_list[index_min],name_list[index_mean],name_list[index_max-1]]

def PreEstimate_filtrations(path_filtrations,save_path,project_name, num_bins='auto',qt=[1,99],max_homology_dim=1):
    lower_bound = []
    upper_bound = []
    for i in range(max_homology_dim+1):
        lower_bound.append([])
        upper_bound.append([])

    item_list = os.listdir(path_filtrations)
    for item in item_list:
        item_path = os.path.join(path_filtrations, item)
        if item.endswith(".txt"):
            raw_filtration = TXTtoList(item_path)
            filtration = []
            for row in raw_filtration:
                vertices_set = list(map(int, row[:-1]))
                filtration_value = math.sqrt(float(row[-1]))
                filtration.append([vertices_set, filtration_value])
            barcodes = SimplicialFiltrationtoBarcode(filtration, dim=max_homology_dim+1)

            for i in range(max_homology_dim+1):
                birth = []
                death = []
                if len(barcodes[i])>0:
                    for bar in barcodes[i]:
                        if bar[1]<=100:
                            birth.append(bar[0])
                            death.append(bar[1])
                    if len(birth)>0:
                        lower_bound[i].append(np.percentile(birth,0))
                        upper_bound[i].append(np.percentile(death,1))

    result = []
    for i in range(max_homology_dim+1):
        DistributionHist(save_path, lower_bound[i], 'min_birthtime(' + project_name + ',h' + str(i) + ')', num_bins=num_bins, qt=qt)
        DistributionHist(save_path, upper_bound[i], 'max_deathtime(' + project_name + ',h' + str(i) + ')', num_bins=num_bins, qt=qt)
        result.append([np.percentile(lower_bound[i],qt[0]),np.percentile(upper_bound[i],qt[1])])
    return result

def FiltrationstoCharactertable(filepath, savepath, csvfile, project_name, name_function, bar_range=[0, 8], style='Normal', function=ListStatistics):
    csv_file = pd.read_csv(csvfile)
    character_table = []
    for idx, row in csv_file.iterrows():
        current_index = row['index']
        feature_vector = []
        for filename in name_function(current_index):
            raw_filtration = TXTtoList(filepath + filename)
            filtration = []
            for row in raw_filtration:
                raw_vertices_set = list(map(int, row[:-1]))
                vertices_set = [abs(x) for x in raw_vertices_set]
                filtration_value = math.sqrt(float(row[-1]))
                filtration.append([vertices_set, filtration_value])
            barcode = SimplicialFiltrationtoBarcode(filtration)
            feature_vector = feature_vector + BarcodetoStatistics(barcode[0], bar_range=bar_range, style=style,
                                                                  function=function)
            feature_vector = feature_vector + BarcodetoStatistics(barcode[1], bar_range=bar_range, style=style,
                                                                  function=function)
        character_table.append(feature_vector)
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
    np.savetxt(savepath + 'character_table_' + project_name + '.txt', character_table)
    
def FiltrationstoCharactertable_op(filepath, savepath, csvfile, project_name, name_function, bar_range=[-8, 0], style='Normal', function=ListStatistics):
    csv_file = pd.read_csv(csvfile)
    character_table = []
    for idx, row in csv_file.iterrows():
        current_index = row['index']
        feature_vector = []
        for filename in name_function(current_index):
            raw_filtration = TXTtoList(filepath + filename)
            filtration = []
            for row in raw_filtration:
                raw_vertices_set = list(map(int, row[:-1]))
                vertices_set = [abs(x) for x in raw_vertices_set]
                filtration_value = -float(row[-1])
                filtration.append([vertices_set, filtration_value])
            barcode = SimplicialFiltrationtoBarcode(filtration)
            feature_vector = feature_vector + BarcodetoStatistics(barcode[0], bar_range=bar_range, style=style,
                                                                  function=function)
            feature_vector = feature_vector + BarcodetoStatistics(barcode[1], bar_range=bar_range, style=style,
                                                                  function=function)
        character_table.append(feature_vector)
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
    np.savetxt(savepath + 'character_table_' + project_name + '.txt', character_table)

def FiltrationstoCharactertable_aug(filepath, savepath, csvfile, project_name, name_function, bar_range=[0, 8], style='Normal', function=ListStatistics):
    csv_file = pd.read_csv(csvfile)
    last_index = -1
    num_aug = 0
    character_table = []
    for idx, row in csv_file.iterrows():
        current_index = row['index']
        if current_index != last_index:
            num_aug = 0
            last_index = current_index
        else:
            num_aug = num_aug + 1
        feature_vector = []
        for filename in name_function(current_index, num_aug):
            raw_filtration = TXTtoList(filepath + filename)
            filtration = []
            for row in raw_filtration:
                raw_vertices_set = list(map(int, row[:-1]))
                vertices_set = [abs(x) for x in raw_vertices_set]
                filtration_value = math.sqrt(float(row[-1]))
                filtration.append([vertices_set, filtration_value])
            barcode = SimplicialFiltrationtoBarcode(filtration)
            feature_vector = feature_vector + BarcodetoStatistics(barcode[0], bar_range=bar_range, style=style,
                                                                  function=function)
            feature_vector = feature_vector + BarcodetoStatistics(barcode[1], bar_range=bar_range, style=style,
                                                                  function=function)
        character_table.append(feature_vector)
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
    np.savetxt(savepath + 'character_table_' + project_name + '.txt', character_table)

def FiltrationstoCharactertable_aug_op(filepath, savepath, csvfile, project_name, name_function, bar_range=[-8, 0], style='Normal', function=ListStatistics):
    csv_file = pd.read_csv(csvfile)
    last_index = -1
    num_aug = 0
    character_table = []
    for idx, row in csv_file.iterrows():
        current_index = row['index']
        if current_index != last_index:
            num_aug = 0
            last_index = current_index
        else:
            num_aug = num_aug + 1
        feature_vector = []
        for filename in name_function(current_index, num_aug):
            raw_filtration = TXTtoList(filepath + filename)
            filtration = []
            for row in raw_filtration:
                raw_vertices_set = list(map(int, row[:-1]))
                vertices_set = [abs(x) for x in raw_vertices_set]
                filtration_value = -float(row[-1])
                filtration.append([vertices_set, filtration_value])
            barcode = SimplicialFiltrationtoBarcode(filtration)
            feature_vector = feature_vector + BarcodetoStatistics(barcode[0], bar_range=bar_range, style=style,
                                                                  function=function)
            feature_vector = feature_vector + BarcodetoStatistics(barcode[1], bar_range=bar_range, style=style,
                                                                  function=function)
        character_table.append(feature_vector)
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
    np.savetxt(savepath + 'character_table_' + project_name + '.txt', character_table)

def FiltrationstoCharactertable_bins(filepath, savepath, csvfile, project_name, name_function, num_bins, bar_range_set=[[0, 8],[0,8]], style='Normal',need_pre_estimate=False, function=ListStatistics):
    if need_pre_estimate:
        bar_range_set = PreEstimate_filtrations(filepath,savepath,project_name)
        print(project_name + ' bar_range_set=', bar_range_set)
    csv_file = pd.read_csv(csvfile)
    last_index = -1
    num_aug = 0
    character_table = []
    for idx, row in csv_file.iterrows():
        current_index = row['index']
        if current_index != last_index:
            num_aug = 0
            last_index = current_index
        else:
            num_aug = num_aug + 1
        feature_vector = []
        for filename in name_function(current_index, num_aug):
            raw_filtration = TXTtoList(filepath + filename)
            filtration = []
            for row in raw_filtration:
                vertices_set = list(map(int, row[:-1]))
                filtration_value = math.sqrt(float(row[-1]))
                filtration.append([vertices_set, filtration_value])
            barcode = SimplicialFiltrationtoBarcode(filtration)
            feature_vector = feature_vector + BarcodetoStatistics_bins(barcode[0], num_bins, barcode_range=bar_range_set[0],
                                                                       style=style, function=function)
            feature_vector = feature_vector + BarcodetoStatistics_bins(barcode[1], num_bins, barcode_range=bar_range_set[1],
                                                                       style=style, function=function)
        character_table.append(feature_vector)
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
    np.savetxt(savepath + 'character_table('+project_name+').txt', character_table)

def FiltrationstoBetticurvetable(filepath, savepath, csvfile, project_name, name_function, num_bins, bar_range_set=[[0,10],[0,10]], need_pre_estimate=False):
    if need_pre_estimate:
        bar_range_set = PreEstimate_filtrations(filepath,savepath,project_name)
        print(project_name + ' bar_range_set=', bar_range_set)
    csv_file = pd.read_csv(csvfile)
    last_index = -1
    num_aug = 0
    character_table = []
    for idx, row in csv_file.iterrows():
        current_index = row['index']
        if current_index != last_index:
            num_aug = 0
            last_index = current_index
        else:
            num_aug = num_aug + 1
        feature_vector = []
        for filename in name_function(current_index, num_aug):
            raw_filtration = TXTtoList(filepath + filename)
            filtration = []
            for row in raw_filtration:
                vertices_set = list(map(int, row[:-1]))
                filtration_value = math.sqrt(float(row[-1]))
                filtration.append([vertices_set, filtration_value])
            barcode = SimplicialFiltrationtoBarcode(filtration)
            feature_vector = feature_vector + Barcodetobetticurve(barcode[0], barcode_range=bar_range_set[0], num_grid = num_bins)
            feature_vector = feature_vector + Barcodetobetticurve(barcode[0], barcode_range=bar_range_set[0], num_grid = num_bins)
        character_table.append(feature_vector)
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
    np.savetxt(savepath + 'betticurve_table('+project_name+').txt', character_table)