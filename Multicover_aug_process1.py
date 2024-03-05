import os
import sys
sys.path.append('./basic_code/')
import TDAbasic as tda
import SMILEprocess as smile
import yaml

finetune_config = yaml.load(open("config.yaml", "r"), Loader=yaml.FullLoader)
print(finetune_config)

project_name = finetune_config['dataset_name']
### Create necessary file folder
save_path_pre = finetune_config['save_path_pre'] + project_name
os.makedirs(save_path_pre + '/augmented', mode=0o755, exist_ok=True)
os.makedirs(save_path_pre + '/coordinates', mode=0o755, exist_ok=True)
os.makedirs(save_path_pre + '/filtrations', mode=0o755, exist_ok=True)
os.makedirs(save_path_pre + '/filtrations/seperate_filtration', mode=0o755, exist_ok=True)


### SMILE string to coordinates with atom type (default is to get the augmentation file first)
pre = "./data/"
filename = pre + project_name
savepath = save_path_pre + '/augmented/'
label_list = ['value']
error_list = smile.SMILEsCsvtoCoordsSet_aug(filename,savepath,label_list,style = 'pending process')
print('error list=', error_list)

### seperate different coordinates set based on atom combination
filepath = save_path_pre +'/augmented/'
savepath = save_path_pre +'/coordinates/'
atom_combination = [['all'],['C'],['C','N'],['C','N','O'],['C','O'],['N','O'],['no','C','H'],['no','H']]
tda.AtomCombination(filepath=filepath,savepath=savepath,atom_combination=atom_combination)

### Compute Delaunay slice for each coordinates set
coords_path = save_path_pre +'/coordinates/'
filtrations_path = save_path_pre +'/filtrations/'
item_list = os.listdir(coords_path)
for item in item_list:
    item_path = os.path.join(coords_path, item)
    if os.path.isfile(item_path) and item.endswith(".txt"):
        file_name = os.path.splitext(item)[0]
        os.system(
            "cd " + finetune_config['rhomboidtiling_path'] + " && ./main "+ item_path +" "+ filtrations_path + file_name +"_fslices.txt 3 4 fslices")

