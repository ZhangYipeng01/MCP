# Author: Zhang Yipeng
# Email: yipeng001@e.ntu.edu.sg
#=================================================
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem,rdForceFieldHelpers
from sentence_transformers import SentenceTransformer
import os
import re
import copy

###### Basic file-dealing functions
def Listtotxt(List,filepath):
    file = open(filepath, 'w')
    for row in List:
        line = ' '.join(str(element) for element in row) + '\n'
        file.write(line)
    file.close()

# The below code is taken from: https://www.programcreek.com/python/example/123334/rdkit.Chem.CombineMols Example4
def connect_frags(mol1, mol2, atom1, atom2):
    combined = Chem.CombineMols(mol1, mol2)
    emol = Chem.EditableMol(combined)
    neighbor1_idx = atom1.GetNeighbors()[0].GetIdx()
    neighbor2_idx = atom2.GetNeighbors()[0].GetIdx()
    atom1_idx = atom1.GetIdx()
    atom2_idx = atom2.GetIdx()
    bond_order = atom2.GetBonds()[0].GetBondType()
    emol.AddBond(neighbor1_idx,
        neighbor2_idx + mol1.GetNumAtoms(),
        order=bond_order)
    emol.RemoveAtom(atom2_idx + mol1.GetNumAtoms())
    emol.RemoveAtom(atom1_idx)
    mol = emol.GetMol()
    return mol

def connect_frags(mol1, mol2, atom1, atom2):
    combined = Chem.CombineMols(mol1, mol2)
    emol = Chem.EditableMol(combined)
    neighbor1_idx = atom1.GetNeighbors()[0].GetIdx()
    neighbor2_idx = atom2.GetNeighbors()[0].GetIdx()
    atom1_idx = atom1.GetIdx()
    atom2_idx = atom2.GetIdx()
    bond_order = atom2.GetBonds()[0].GetBondType()
    emol.AddBond(neighbor1_idx,
        neighbor2_idx + mol1.GetNumAtoms(),
        order=bond_order)
    emol.RemoveAtom(atom2_idx + mol1.GetNumAtoms())
    emol.RemoveAtom(atom1_idx)
    mol = emol.GetMol()
    return mol

def generate_augmented_smiles(input_smiles):
    canon_smiles = Chem.CanonSmiles(input_smiles)
    monomer = Chem.MolFromSmiles(canon_smiles)
    dum_atms_monomer = [atom for atom in monomer.GetAtoms() if atom.GetSymbol() == '*']
    if len(dum_atms_monomer)!=2:
        print(input_smiles)
        return 'Error'
    dimer = connect_frags(monomer, monomer, dum_atms_monomer[1], dum_atms_monomer[0])
    dum_atm_idx_dimer = [idx for idx,atom in enumerate(dimer.GetAtoms()) if atom.GetSymbol() == '*']
    sp_bw_dummy_atoms = Chem.rdmolops.GetShortestPath(dimer, dum_atm_idx_dimer[0], dum_atm_idx_dimer[1])[1:-1] #atoms in shortest path between dummy atoms
    sp_bonds = [dimer.GetBondBetweenAtoms(first, second).GetIdx() for first, second in zip(sp_bw_dummy_atoms, sp_bw_dummy_atoms[1:])] #get indices of bonds in shortest path between dummy atoms
    sp_nonring_bonds = [idx for idx in sp_bonds if not dimer.GetBonds()[idx].IsInRing()] #get indices of bonds in SP that do not belong to a ring
    del sp_nonring_bonds[(len(sp_nonring_bonds)-1)//2] #delete bond that joins the two monomers together
    bond_pair_list = np.reshape(sp_nonring_bonds, (2, -1)) #get pairs of bonds to be fragmented in order generate alternate repeating units
    new_ru_list = [] #initialize list of alternate units
    for bond_pair in bond_pair_list.T:
        bond_pair = [val.item() for val in list(bond_pair)] #convert numpy data type to python data type (for rdkit compatibility)
        mol_f = Chem.FragmentOnBonds(dimer, bond_pair, addDummies = True, dummyLabels=[(0, 0), (0, 0)]) #fragment dimer at selected bond pairs
        mol_f1 = list(Chem.rdmolops.GetMolFrags(mol_f, asMols = True))
        new_ru_list.append(mol_f1[1]) #middle fragment corresponds to the alternate repeating unit
    augmented_reps = [monomer] + new_ru_list
    augmented_reps = [Chem.MolToSmiles(each) for each in augmented_reps]
    augmented_reps = [Chem.CanonSmiles(each) for each in augmented_reps]
    return sorted(set(augmented_reps), key=augmented_reps.index)

class DataAugmentation:
    # The below code is taken from: https://github.com/ChangwenXu98/TransPolymer/blob/master/dataset.py
    def __init__(self, aug_indicator=10):
        super(DataAugmentation, self).__init__()
        self.aug_indicator = aug_indicator

    def rotate_atoms(self, li, x):
        return (li[x % len(li):] + li[:x % len(li)])

    """Generate SMILES"""
    def generate_smiles(self, smiles):
        smiles_list = []
        try:
            mol = Chem.MolFromSmiles(smiles)
        except:
            mol = None
        if mol != None:
            n_atoms = mol.GetNumAtoms()
            n_atoms_list = [nat for nat in range(n_atoms)]
            if n_atoms != 0:
                for iatoms in range(n_atoms):
                    n_atoms_list_tmp = self.rotate_atoms(n_atoms_list, iatoms)  # rotate atoms' index
                    nmol = Chem.RenumberAtoms(mol, n_atoms_list_tmp)  # renumber atoms in mol
                    try:
                        smiles = Chem.MolToSmiles(nmol,
                                                  isomericSmiles=True,  # keep isomerism
                                                  kekuleSmiles=False,  # kekulize or not
                                                  rootedAtAtom=-1,  # default
                                                  canonical=False,  # canonicalize or not
                                                  allBondsExplicit=False,  #
                                                  allHsExplicit=False)  #
                    except:
                        smiles = 'None'
                    smiles_list.append(smiles)
            else:
                smiles = 'None'
                smiles_list.append(smiles)
        else:
            try:
                smiles = Chem.MolToSmiles(mol,
                                          isomericSmiles=True,  # keep isomerism
                                          kekuleSmiles=False,  # kekulize or not
                                          rootedAtAtom=-1,  # default
                                          canonical=False,  # canonicalize or not
                                          allBondsExplicit=False,  #
                                          allHsExplicit=False)  #
            except:
                smiles = 'None'
            smiles_list.append(smiles)
        smiles_array = pd.DataFrame(smiles_list).drop_duplicates().values
        # """
        if self.aug_indicator is not None:
            smiles_aug = smiles_array[1:, :]
            np.random.shuffle(smiles_aug)
            smiles_array = np.vstack((smiles_array[0, :], smiles_aug[:self.aug_indicator-1, :]))
        regular_smiles_array = []
        for smiles in smiles_array:
            regular_smiles_array.append(smiles[0])
        regular_smiles_array = np.array(regular_smiles_array)
        return regular_smiles_array

    def smiles_augmentation(self, df, style = 'stable'):
        column_list = df.columns
        data_aug = np.zeros((1, df.shape[1]))
        for i in range(df.shape[0]):
            smiles = df.iloc[i, 0]
            prop = df.iloc[i, 1:]
            if style == 'stable':
                smiles_array = np.array(generate_augmented_smiles(smiles))
            if style == 'extreme':
                smiles_array = self.generate_smiles(smiles)
            smiles_array = smiles_array[:, np.newaxis]
            if 'None' not in smiles_array:
                prop = np.tile(prop, (len(smiles_array), 1))
                data_new = np.hstack((smiles_array, prop))
                data_aug = np.vstack((data_aug, data_new))
        df_aug = pd.DataFrame(data_aug[1:, :], columns=column_list)
        return df_aug

def cap_with_next_atom(smiles_string):
    try:
        smiles_string = Chem.CanonSmiles(smiles_string)
        molecule = Chem.MolFromSmiles(smiles_string)
        new_atoms = []
        for atom in molecule.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol == "*":
                new_atoms.append(''.join([nbr.GetSymbol() for nbr in atom.GetNeighbors()]))
        new_atoms = ["[Si]" if (each=="Si" or each=="si") else each for each in new_atoms] #because parsing Silicon atoms is tricky, see https://github.com/PatWalters/rd_filters/blob/master/Notes.txt
        capped_molecule = Chem.ReplaceSubstructs(molecule, Chem.MolFromSmiles('*'), Chem.MolFromSmiles(new_atoms[1]), replaceAll = False)[0]
        capped_molecule = Chem.ReplaceSubstructs(capped_molecule, Chem.MolFromSmiles('*'), Chem.MolFromSmiles(new_atoms[0]), replaceAll = False)[0]
        return Chem.CanonSmiles(Chem.MolToSmiles(capped_molecule))
    except:
        print("CappingFailed for", smiles_string)
        return "CappingFailed"

## adapted from Kevin

def SMILESextension(smiles,n_mer):
    if re.search(r"\[\*\].*\[\*\]", smiles):
        updated_smiles = smiles
    else:
        updated_smiles = re.sub(r"(?<!\[\*)\*(?!\*\])", "[*]", smiles)

    smiles = updated_smiles

    # Initialize variables
    highest_num = 0

    # Find the highest number in the SMILES string
    for char in smiles:
        try:
            if int(char) > highest_num:
                highest_num = int(char)
        except:
            pass

    # Start with the original SMILES string in the list
    smiles_all = [smiles, ]

    # Generate n-mer versions of the SMILES string
    for n in range(n_mer - 1):
        smiles_n = copy.deepcopy(smiles)
        if highest_num > 0:
            for j in range(0, highest_num):
                if j + 1 + (n + 1) * highest_num < 10:
                    smiles_n = re.sub(f'{j + 1}', f'{j + 1 + (n + 1) * highest_num}', smiles_n)
                else:
                    pattern = re.compile(f"(?<!%[0-9]|.%){j + 1}")  # Avoid matching digits part of larger numbers
                    smiles_n = pattern.sub(fr'%{j + 1 + (n + 1) * highest_num}', smiles_n)
        smiles_n = re.sub('\[[0-9]H\]|\[\%[0-9][0-9]H\]', r'[2H]', smiles_n)
        smiles_all.append(smiles_n)

    # Reverse SMILES strings for concatenation
    rsmiles_all = [a[::-1] for a in smiles_all]

    # Concatenate reversed SMILES strings
    for m in range(n_mer - 1):
        rsmiles_all[0] = rsmiles_all[0].replace(r']*[', rsmiles_all[m + 1][:-3], 1)

    # Return the dimerized SMILES string
    multimer_smiles_smiles = rsmiles_all[0][::-1]

    return multimer_smiles_smiles

# Input csv file should use "smiles" to name the column of SMILES string
def PolymerSMILESextension(filename, rept_times):
    df = pd.read_csv(filename + ".csv")
    df['processed_strings'] = df['smiles'].apply(lambda x: SMILESextension(x,rept_times))
    df.to_csv(filename + "_ext_rept"+str(rept_times)+".csv", index=False, encoding='utf-8')

#  Input is the name of a Csv file contain 3 kinds of columns: "smiles", "index" and a list of name of the label. Output is an augmented csv file.    
def PolymerSMILEAugmentation(filename,label_list,aug_indicator = 10,style = 'stable'):
    data = pd.read_csv(filename + ".csv")
    DataAug = DataAugmentation(aug_indicator)
    data_Aug = DataAug.smiles_augmentation(data.loc[:, ["smiles", "index"] + label_list],style)
    smile_strings = data_Aug["smiles"].tolist()
    processed_string = []
    for string in smile_strings:
        result = cap_with_next_atom(string)
        if result != "CappingFailed":
            processed_string.append(result)
    data_Aug['processed_strings'] = processed_string
    if style == 'stable':
        data_Aug.to_csv(filename + "_aug.csv",index=False, encoding='utf-8')
    else:
        data_Aug.to_csv(filename + "_aug(aug_indicator="+str(aug_indicator)+",style="+style+").csv", index=False, encoding='utf-8')       

def is_valid_smiles(smiles_string):
    try:
        molecule = Chem.MolFromSmiles(smiles_string)
        if molecule is not None and molecule.GetNumAtoms() > 0:
            return True
        else:
            return False
    except Exception as e:
        # Handle exceptions (e.g., invalid SMILES syntax)
        return False

def SMILEtoDenseFingerprint(smile_string):
    canon_smile_string = Chem.CanonSmiles(smile_string)
    polyBERT = SentenceTransformer('kuelumbus/polyBERT')
    embeddings = polyBERT.encode(canon_smile_string)
    return embeddings

def SMILEtoMol(smile_string):
    # Replace asterisks with placeholder atoms (e.g., bromine)
    smile_string = smile_string.replace('*', 'Br')
    if is_valid_smiles(smile_string) == False:
        return 'unsuccessful'
    max_tries = 10
    num_tries = 0
    coords = []  
    for i in range(max_tries):
        num_tries += 1
        try:
            mol = Chem.MolFromSmiles(smile_string)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=10, useRandomCoords=True)
            rdForceFieldHelpers.MMFFOptimizeMolecule(mol, maxIters=2147483647)
            atoms = mol.GetAtoms()
            for atom in atoms:
                x, y, z = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                symbol = atom.GetSymbol()
                # Handle placeholder atoms (e.g., remove or mark them)
                #if symbol == 'Br':
                    #symbol = '*'  # Replace back with asterisk or handle as needed
                coords.append([x, y, z, symbol])

        except BaseException as be:
            if num_tries < max_tries:
                continue
            else:
                return 'unsuccessful'
        else:
            return coords

##### ensemble functions
def SMILEsCsvtoDenseFingerprint(filename,savepath,label_list,style='pending process',name='',aug_indicator = 10,aug_style = 'stable'):
    style_key = 'smiles'
    if style == 'pending process':
        PolymerSMILEAugmentation(filename,label_list,aug_indicator=aug_indicator,style=aug_style)
        style_key = 'processed_strings'
        if aug_style == 'stable':
            filename = filename + '_aug'
        else:
            filename = filename + '_aug(aug_indicator='+str(aug_indicator)+',style='+aug_style+')'
    if style == 'processed':
        style_key = 'processed_strings'
    
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
        
    polymers_aug_df = pd.read_csv(filename + '.csv')
    
    character_table = []
    for idx, row in polymers_aug_df.iterrows():
        smile_string = row[style_key]
        dense_feature = SMILEtoDenseFingerprint(smile_string)
        character_table.append(dense_feature)
    np.savetxt(savepath + 'character_table'+name+'(DensityFingerprint).txt',character_table)
    
def SMILEsCsvtoCoordsSet(filename, savepath, savename='', style_key = 'smiles'):
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
        
    polymers_aug_df = pd.read_csv(filename + '.csv')
    error_information = []
    for idx, row in polymers_aug_df.iterrows():
        current_index = row['index']
        smile_string = row[style_key]
        coords = SMILEtoMol(smile_string)
        if coords == 'unsuccessful':
            error_information.append([idx,smile_string])
        else:
            Listtotxt(coords, savepath + str(
                current_index)+ savename + '.txt')
    return error_information
    
def SMILEsCsvtoCoordsSet_aug(filename,savepath,label_list,style='pending process',aug_indicator = 10,aug_style = 'stable'):
    style_key = 'smiles'
    if style == 'pending process':
        PolymerSMILEAugmentation(filename,label_list,aug_indicator=aug_indicator,style=aug_style)
        style_key = 'processed_strings'
        if aug_style == 'stable':
            filename = filename + '_aug'
        else:
            filename = filename + '_aug(aug_indicator='+str(aug_indicator)+',style='+aug_style+')'
    if style == 'processed':
        style_key = 'processed_strings'
    
    if not os.path.exists(savepath):
        os.makedirs(savepath, mode=0o755, exist_ok=True)
        
    polymers_aug_df = pd.read_csv(filename + '.csv')
    last_index = -1
    num_aug = 0
    error_information = []
    for idx, row in polymers_aug_df.iterrows():
        current_index = row['index']
        if current_index != last_index:
            num_aug = 0
            last_index = current_index
        else:
            num_aug = num_aug + 1

        smile_string = row[style_key]
        coords = SMILEtoMol(smile_string)
        if coords == 'unsuccessful':
            error_information.append([idx,smile_string])
        else:
            Listtotxt(coords, savepath + str(
                current_index) + '_aug' + str(num_aug) + '.txt')
    return error_information