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

import time

import math
import numpy
import pycuda.tools as cutools
import pycuda.driver as cuda
import pycuda.gpuarray as gpuarray
import pycuda.cumath as gpumath
import periodictable as pt
import configuration as config
from controllers import protein_reader as reader
from controllers import vector_controller as vectors


def sin_c(x):
    if x == 0:
        return 1.0
    return math.sin(x) / x

def factors_to_gpu(frames,atom_factors):
    gpu_factors = []
    for frame in frames:
        f_factor = []
        for atom in frame:
            f_factor.append(atom_factors[atom][1])
        n_factor = numpy.asarray(f_factor)
        gpu_factor = gpuarray.to_gpu(n_factor)
        gpu_factors.append(gpu_factor)
    return gpu_factors

def amplitude_compute_gpu_version_4(gpu_frame_x,gpu_frame_y,gpu_frame_z,gpu_factor,vector):
    gpu_frame_0 = gpu_frame_x
    gpu_frame_1 = gpu_frame_y
    gpu_frame_2 = gpu_frame_z
    gpu_dot_0 = gpu_frame_0.__mul__(vector[0])
    gpu_dot_1 = gpu_frame_1.__mul__(vector[1])
    gpu_dot_2 = gpu_frame_2.__mul__(vector[2])
    gpu_dot = gpu_dot_0.__add__(gpu_dot_1)
    gpu_dot_0.gpudata.free()
    gpu_dot_1.gpudata.free()
    gpu_dot_res = gpu_dot.__add__(gpu_dot_2)
    gpu_dot_2.gpudata.free()
    gpu_dot.gpudata.free()
    gpu_sin = gpumath.sin(gpu_dot_res)
    gpu_cos = gpumath.cos(gpu_dot_res)
    gpu_imag = gpu_sin.__mul__(gpu_factor)
    gpu_sin.gpudata.free()
    gpu_real = gpu_cos.__mul__(gpu_factor)
    gpu_cos.gpudata.free()
    g_real = gpuarray.sum(gpu_real)
    gpu_real.gpudata.free()
    g_imag = gpuarray.sum(gpu_imag)
    gpu_imag.gpudata.free()
    return g_real.get(), g_imag.get()

def atom_factor_compute(data, module, groups_array):
    res = {} 
    for item in data:
        if item.strip() in groups_array:
            el_num = int(groups_array[item][1])
            factor = 0
            for atom in groups_array[item][0].split(':'):
                factor += pt.elements.isotope(atom).xray.f0(module)
        else:
            if len(item) >1:
                atom = item[0]+item[1].lower()
            else:
                atom = item
            el_num = pt.elements.isotope(atom).number
            factor = pt.elements.isotope(atom).xray.f0(module)
        res[item] = [el_num,factor]
    return res

def get_cell_size(fname):
    fd  = open(fname)
    x = float(fd.readline())
    y = float(fd.readline())
    z = float(fd.readline())
    fd.close()
    return [x,y,z]
def intensities_compute(protein_frames, buffer_frames, data, v_m):
    uniform_v = vectors.generate_vectors(1, config.N_Z, config.N_FI)
    atom_factor_one = {}
    cell_size = get_cell_size(config.ROOT_PATH +'/dimensions.dat')
    group_dict = reader.get_group_dict(config.GROUP_DESCR)
    atom_factor = atom_factor_compute(protein_frames[4], v_m/10, group_dict)
    intensity_one_module(protein_frames, buffer_frames, atom_factor, uniform_v, v_m,cell_size)


def intensity_one_module(protein_frames, buffer_frames, atom_factors, uniforn_vectors, vector_module,cell_size):
    t1 = time.time()
    intensity = 0
    intensity_one = 0
    intensity_two = 0
    f_len_a = len(protein_frames[3])
    f_len_b = len(buffer_frames[3])
    num_wat_mols = 0
    el_num_in_wat = 0
    for frame in buffer_frames[3]:
        for item in frame:
            el_num_in_wat += atom_factors[item][0]
    el_num_in_wat = el_num_in_wat /f_len_b
    water_frames_x = buffer_frames[0]
    water_frames_y = buffer_frames[1]
    water_frames_z = buffer_frames[2]
    water_factors = factors_to_gpu(buffer_frames[3],atom_factors)
    prot_frames_x = protein_frames[0] 
    prot_frames_y = protein_frames[1]
    prot_frames_z = protein_frames[2]
    prot_factors = factors_to_gpu(protein_frames[3],atom_factors)
    for uniform_vector in uniforn_vectors:
        a_mean_mod_2 = 0
        b_mean_mod_2 = 0
        vector = (
            uniform_vector[0] * vector_module, uniform_vector[1] * vector_module, uniform_vector[2] * vector_module)
        par_amplitude = el_num_in_wat * sin_c(0.5 * cell_size[0] * vector[0]) * sin_c(0.5 * cell_size[1] * vector[1]) * sin_c(0.5 * cell_size[2] * vector[2])
        for frame_num in range(len(protein_frames[3])):
            a_real, a_imag = amplitude_compute_gpu_version_4(prot_frames_x[frame_num],prot_frames_y[frame_num],prot_frames_z[frame_num],prot_factors[frame_num],vector)
            a_mean_mod_2 += math.pow(a_real - par_amplitude, 2) + math.pow(a_imag, 2)
        for frame_num in range(len(buffer_frames[3])):
            b_real, b_imag = amplitude_compute_gpu_version_4(water_frames_x[frame_num],water_frames_y[frame_num],water_frames_z[frame_num],water_factors[frame_num],vector)
            b_mean_mod_2 += math.pow(b_real - par_amplitude, 2) + math.pow(b_imag, 2)
        a_mean_mod_2 = a_mean_mod_2 / f_len_a
        b_mean_mod_2 = b_mean_mod_2 / f_len_b
        intensity += a_mean_mod_2 - b_mean_mod_2
        intensity_one += a_mean_mod_2
        intensity_two += b_mean_mod_2
    print('for q =' + str(vector_module) + ' intensity = ' +
          str(intensity) + 'Time spent on the intensity computing(s) = ' +
          str((math.ceil((time.time() - t1) * 10) * 0.1)))
    intensity = math.log10(intensity / (config.N_Z * config.N_FI))
    reader.data_to_file((vector_module, intensity), "POINT_Q=" +
                        str(math.floor((vector_module + 0.0001) * 1000.0) / 1000.0) + "_NZ=" +
                        str(config.N_Z) + "_NFI=" + str(config.N_FI) + "_P@Wfr = " + str(f_len_a) + "_CWfr = " +
                        str(f_len_b) + "_indent = " + str(config.INDENT) + ".rslt")

     