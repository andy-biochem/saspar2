# This file is part of SASPAR2 software package.
# Copyright 2020 Andrei M. Semenov
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#  
#     http://www.apache.org/licenses/LICENSE-2.0
#  
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import sys
import warnings
from pathlib import Path
import periodictable as pt
import numpy

import configuration as config



def load_amino_files_dict(filepath):
    amino_dict = {}
    amino_dir = os.listdir(filepath)
    for fname in amino_dir:
        fd = open(filepath+'/'+fname)
        lines = fd.readlines()
        fd.close()
        c_di ={}
        for line in lines:
            arr = line.strip().split(' ')
            c_di[arr[0].strip()] = arr[1].strip()
        amino_dict[fname[0:3]] = c_di
    return amino_dict

def get_group_dict(fname):
    fd = open(fname)
    lines = fd.readlines()
    fd.close()
    res = {}
    for line in lines:
        arr = line.split('\t')
        res[arr[0]] = [arr[2].strip(),int(arr[1])]
    return res
def load_pdb_data(lines):
    aminoacid_data = []
    for line in lines:
        linea = line.strip()
        if linea.startswith('ATOM') or linea.startswith('HETATM'):
            aminoacid_data.append((linea[17:20],linea[13:16].strip(),int('0'),float(linea[30:38]),float(linea[38:46]),float(linea[46:54])))
    return aminoacid_data


def get_protein_data(protein_data,amino_dict):
    protein_data_output = []
    for atom_i in protein_data:
        if atom_i[0] in amino_dict:
            if atom_i[1].strip().upper() in amino_dict[atom_i[0]]:
                gr_name = amino_dict[atom_i[0]][atom_i[1].strip().upper()]
                protein_data_output.append((atom_i[2], atom_i[3], atom_i[4], atom_i[5], gr_name))
            else:
                print('WARNING! unknown atom '+atom_i[1]+' in residue '+atom_i[0])
        else:
            print('WARNING!!! Unknown residue : '+atom_i[0])
    return protein_data_output

def data_to_file(data, filename):
    root_path = Path(config.OUTPUT_PATH)
    os.chdir(root_path)

    with open(filename, "w") as file:
        numpy.savetxt(file, data, fmt='%5s', delimiter=' ', newline='\n')
