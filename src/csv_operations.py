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
from subprocess import Popen
from controllers import protein_reader as reader
import configuration as config
import pycuda.gpuarray as gpuarray
import numpy
import time
import math
import sys

def csv_read(fname):
    fd = open(fname)
    line = fd.readline()
    group_list = []
    gpu_frames_x = []
    gpu_frames_y = []
    gpu_frames_z = []
    factors = []
    i = int(line.strip().split(';')[0])
    f_frame_0 = []
    f_frame_1 = []
    f_frame_2 = []
    factor_arr= []
    while not line.startswith('END'):
        crr_arr = line.strip().split(';')
        try:
            frame_num = int(crr_arr[0])
        except:
            continue
        if int(crr_arr[0]) != i:
            n_frame_0 = numpy.asarray(f_frame_0)
            n_frame_1 = numpy.asarray(f_frame_1)
            n_frame_2 = numpy.asarray(f_frame_2)
            gpu_frame_0 = gpuarray.to_gpu(n_frame_0)
            gpu_frame_1 = gpuarray.to_gpu(n_frame_1)
            gpu_frame_2 = gpuarray.to_gpu(n_frame_2)
            gpu_frames_x.append(gpu_frame_0)
            gpu_frames_y.append(gpu_frame_1)
            gpu_frames_z.append(gpu_frame_2)
            factors.append(factor_arr)
            f_frame_0 = []
            f_frame_1 = []
            f_frame_2 = []
            factor_arr= []
            i = int(crr_arr[0])
        f_frame_0.append(float(crr_arr[1]))
        f_frame_1.append(float(crr_arr[2]))
        f_frame_2.append(float(crr_arr[3]))
        group_name = crr_arr[4].strip()
        factor_arr.append(group_name)
        if not (group_name in group_list):
            group_list.append(group_name)
        line = fd.readline()
    fd.close()
    n_frame_0 = numpy.asarray(f_frame_0)
    n_frame_1 = numpy.asarray(f_frame_1)
    n_frame_2 = numpy.asarray(f_frame_2)
    gpu_frame_0 = gpuarray.to_gpu(n_frame_0)
    gpu_frame_1 = gpuarray.to_gpu(n_frame_1)
    gpu_frame_2 = gpuarray.to_gpu(n_frame_2)
    gpu_frames_x.append(gpu_frame_0)
    gpu_frames_y.append(gpu_frame_1)
    gpu_frames_z.append(gpu_frame_2)
    factors.append(factor_arr)
    return [gpu_frames_x,gpu_frames_y,gpu_frames_z,factors,group_list]
def csv_write(fname,proteins):
    fd  = open(fname,'w+')
    for i in range(len(proteins)):
        for j in range(len(proteins[i])):
            line = str(i)+';'+str(proteins[i][j][1])+';'+str(proteins[i][j][2])+';'+str(proteins[i][j][3])+';'+str(proteins[i][j][4])+';\n'
            fd.write(line)
    fd.write('END \n')
    fd.close()
def csv_append(fname,protein,num):
    fd  = open(fname,'a+')
    for j in range(len(protein)):
        line = str(num)+';'+str(protein[j][1])+';'+str(protein[j][2])+';'+str(protein[j][3])+';'+str(protein[j][4])+';\n'
        fd.write(line)
    fd.close()
    return num+1
def csv_final(fname):
    fd  = open(fname,'a+')
    fd.write('END \n')
    fd.close()

def pdb2csv(thrds):
    print("Files loading....")
    time_file_load = time.time()
    amino_dict = reader.load_amino_files_dict(config.AMINOACID_PATH)
    if thrds == 2:
        pr = Popen(['python3',config.ROOT_PATH+'/prot2csv.py'])
    else:
        proteins = []
        for file in os.listdir(config.AMINO_PATH):
            if file.endswith(".pdb"):
                protein_pdb = reader.load_pdb_data(file, config.AMINO_PATH + file)
                protein = reader.get_protein_data(protein_pdb,amino_dict)
                proteins.append(protein)

        csv_write(config.ROOT_PATH+'/proteins.csv',proteins)
    waters = []
    for file in os.listdir(config.WATER_PATH):
        print('Processing...'+file)
        if file.endswith(".pdb"):
            water_pdb = reader.load_pdb_data(file, config.WATER_PATH + file)
            water = reader.get_protein_data(water_pdb,amino_dict)
            waters.append(water)
    csv_write(config.ROOT_PATH+'/waters.csv',waters)
    print('Time spent on the files loading (s) = ' + str((math.ceil((time.time() - time_file_load) * 10) * 0.1)))
    if thrds == 2:
        while pr.poll() == None:
            time.sleep(3)
    print('Time spent on the protein and water loading (s) = ' + str((math.ceil((time.time() - time_file_load) * 10) * 0.1)))




