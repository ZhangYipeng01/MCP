import os
####### Basic file-dealing functions
def TXTtoList(filepath,return_type='str'):
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

def Listtotxt(List,filepath):
    file = open(filepath, 'w')
    for row in List:
        line = ' '.join(str(element) for element in row) + '\n'
        file.write(line)
    file.close()

######## Multi-cover functions
def Combnumber(k,n): #Calculate the combination number

    numerator_index = 1
    denominator_index = n
    numerator = 1
    denominator = 1
    for i in range(k):
        numerator = numerator * numerator_index
        denominator = denominator * denominator_index
        numerator_index = numerator_index +1
        denominator_index = denominator_index -1
    return int(denominator/numerator)

def GiveIndex(s): #give an index to s='[x_1,x_2,...,x_k]'
    s = s.strip('[]')
    l = s.split(',')

    sum = 1
    for i in range(len(l)):
        sum = sum + Combnumber(i+1,int(l[i]))
    return sum

def SeperateFslices(filepath,savepath,filename,max_order):

    fp = open(filepath + filename + '.txt', 'r', encoding='utf-8')
    content = fp.read()
    fp.close()
    rowlist = content.splitlines()

    GivePtsIndex = {"null": "null"}
    index = 1  # max index = 4
    slice_set = []
    for i in range(max_order):
        slice_set.append([])

    for i in range(len(rowlist)):
        if rowlist[i] == 'Slice ' + str(index) + ':':
            index = index + 1
        else:
            row = []
            row_word = rowlist[i].strip('()')
            row_word = row_word.replace(' ', '')
            row_mid = row_word.split('],')
            if len(row_mid) == 2:
                GivePtsIndex[row_mid[0][1:-1] + ']'] = GiveIndex(row_mid[0][1:-1])
                row.append(GivePtsIndex[row_mid[0][1:-1] + ']'])
            if (len(row_mid) > 2):
                row.append(GivePtsIndex[row_mid[0][1:len(row_mid[0])] + ']'])
                for j in range(1, len(row_mid) - 2):
                    row.append(GivePtsIndex[row_mid[j] + ']'])
            if (len(row_mid) > 2):
                row.append(GivePtsIndex[row_mid[len(row_mid) - 2][:-1] + ']'])
            row.append(row_mid[len(row_mid) - 1])
            slice_set[index - 2].append(row)

    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
    for i in range(max_order):
        Listtotxt(slice_set[i], savepath+filename+'(k='+str(i+1)+').txt')

def SeperateFslices_folder(filepath,savepath,max_order):
    item_list = os.listdir(filepath)
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
    for item in item_list:
        item_path = os.path.join(filepath, item)
        if os.path.isfile(item_path) and item.endswith(".txt"):
            file_name = os.path.splitext(item)[0]
            SeperateFslices(filepath,savepath,file_name,max_order)

def bifi_select(bifi_list,fix_param,fix_param_value):
    seperate_list = []
    if fix_param == 'r':
        fix_param_index = 3
        fil_index = 2
    if fix_param == 'k':
        fix_param_index = 2
        fil_index = 3
    boundary = {'null':'null'}

    processed_list = [
        [int(item) if i != 3 else item for i, item in enumerate(sublist)]
        for sublist in bifi_list
    ]

    for complex in processed_list:
        if complex[fix_param_index]<=fix_param_value*fix_param_value:
            if complex[1]== 0:
                seperate_list.append([complex[0],complex[fil_index]])
            if complex[1] == 1:
                seperate_list.append(complex[4:]+[complex[fil_index]])
                boundary[complex[0]] = complex[4:]
            if complex[1] == 2:
                seperate_list.append([complex[0],complex[fil_index]])
                for edge in complex[4:]:
                    for vertex in boundary[edge]:
                        seperate_list.append([vertex,complex[0],complex[fil_index]])
                    seperate_list.append(boundary[edge]+[complex[0],complex[fil_index]])

    return seperate_list

def SeperateBifi(filepath,savepath,filename,fix_param,max_fix_param_value):
    input_list = TXTtoList(filepath + filename + '.txt', 'float')
    for i in range(max_fix_param_value):
        Listtotxt(bifi_select(input_list, fix_param, i+1), savepath+filename+'('+fix_param+'='+str(i+1)+').txt')

def SeperateBifi_folder(filepath,savepath,fix_param,max_fix_param_value):
    item_list = os.listdir(filepath)
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
    for item in item_list:
        item_path = os.path.join(filepath, item)
        if os.path.isfile(item_path) and item.endswith(".txt"):
            file_name = os.path.splitext(item)[0]
            SeperateBifi(filepath,savepath,file_name,fix_param,max_fix_param_value)