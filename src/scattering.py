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
from controllers import protein_reader as reader
from controllers import intensity_controller as controller
import configuration as config
import time
import math
import sys
import csv_operations as csv_op
import numpy
import pycuda.autoinit
import pycuda.tools as cutools
import pycuda.driver as cuda
import pycuda.gpuarray as gpuarray
import pycuda.cumath as gpumath

def intensity(vect_aa):
    print("CSV loading....")
    time_file_load = time.time()

    proteins = csv_op.csv_read(config.ROOT_PATH+'/proteins.csv')
    waters = csv_op.csv_read(config.ROOT_PATH+'/waters.csv')
    amino_dict = reader.load_amino_files_dict(config.AMINOACID_PATH)

    print('Time spent on the CSV loading (s) = ' + str((math.ceil((time.time() - time_file_load) * 10) * 0.1)))


    for vect in vect_aa.strip().split(';'):
        print("Start intensity compute for q = "+vect)
        time_intensity_compute = time.time()
        controller.intensities_compute(proteins, waters, amino_dict, float(vect))
        print('Time spent on the intensity computing (s) = ' + str((math.ceil((time.time() - time_intensity_compute) * 10) * 0.1)))

intensity(sys.argv[1])



