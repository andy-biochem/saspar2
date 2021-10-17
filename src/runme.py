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

from subprocess import Popen
import os
import configuration as config
from time import sleep
import csv_operations as csv_op
import trimmer
def result_glue():
    res_arr = []
    for fname in os.listdir(config.OUTPUT_PATH):
        if fname.endswith('.rslt'):
            fd = open(config.OUTPUT_PATH+'/'+fname)
            q = float(fd.readline())
            r = float(fd.readline())
            res_arr.append([q,r])
    res_arr.sort()
    out_fn = os.path.split(config.ROOT_PATH)[1]
    fd = open(config.OUTPUT_PATH+'/'+out_fn+'.csv','w+')
    for item in res_arr:
        fd.write(str(item[0])+';'+str(item[1])+';\n')
    fd.close()

v_m = config.VECTOR_START
proces_arr = []
print(1)
if not (os.path.exists(config.ROOT_PATH+'/proteins.csv') and os.path.exists(config.ROOT_PATH+'/waters.csv')):
    trimmer.trim()
for i in range(0,config.VECTOR_COUNT):
    v_m += config.VECTOR_DIF
    pr = Popen(['python3',config.ROOT_PATH+'/scattering.py',str(v_m)])
    proces_arr.append(pr)
    if len(proces_arr) < config.MAX_PROC:
        continue
    while len(proces_arr) == config.MAX_PROC:
        for pr in proces_arr[:]:
            if not pr.poll() == None:
                proces_arr.remove(pr)
        sleep(3)
while not len(proces_arr) == 0:
    for pr in proces_arr[:]:
        if not pr.poll() == None:
            proces_arr.remove(pr)
    sleep(3)
result_glue()
