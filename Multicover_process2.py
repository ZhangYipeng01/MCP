import os
import sys
sys.path.append('./basic_code/')
import TDAbasic as tda
import MultiCover as mc
import shutil
import yaml

finetune_config = yaml.load(open("config.yaml", "r"), Loader=yaml.FullLoader)

project_name = finetune_config['dataset_name']
save_path_pre = finetune_config['save_path_pre'] + project_name

### seperate the fslices into filtrations with different order
filepath = save_path_pre + '/filtrations/'
savepath = save_path_pre + '/filtrations/seperate_filtration/'
max_order = 4
mc.SeperateFslices_folder(filepath, savepath, max_order)

### transfer the filtration into betti curve table

error_indexes = []
filepath = save_path_pre + '/filtrations/seperate_filtration/'
savepath = './result/' + project_name + '_character/'
csv_file = './data/' + project_name + '.csv'
os.makedirs(savepath, mode = 0o777, exist_ok=True)
def name_function(index,aug_index):
    name_list = []
    for order in range(max_order):
        for comb in ['all','C','CN','CNO','CO','NO','noCH','noH']:
            if index not in error_indexes:
                name_list.append(str(index) +'_aug' + str(aug_index) + '_' + comb + '_fslices(k=' + str(order+1) + ').txt')
    return name_list

print(tda.PreEstimate_datasize(path_data = save_path_pre + '/augmented/', save_path = savepath, project_name = project_name))
tda.FiltrationstoCharactertable(filepath, savepath, csv_file, project_name, name_function, bar_range=[0, 8], style='Normal')
#tda.FiltrationstoBetticurvetable(filepath=filepath,savepath=savepath,csvfile=csv_file,project_name=project_name,name_function=name_function, num_bins = 60, bar_range_set = [[0,6],[0,6]])