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

import os,gzip
import bigpdb
import configuration as config
def water_split():
	fname = config.NEW_WATER_PATH
	print(fname)
	if fname.endswith('gz'):
		prot_file = gzip.open(fname,mode='rt')
	else:
		prot_file = open(fname,mode='rt')
	line = prot_file.readline()
	model_num = 0
	current_mdl = []
	while not line.strip() == '':
		if line.startswith('MODEL'):
			print (line.strip())
			if not model_num  == 0:
				current_mdl = bigpdb.rm_element(current_mdl,'H')
				lim_pos = bigpdb.get_lim_pos(current_mdl,line_type = 'ALL')
				shift_arr = [(lim_pos[1]+lim_pos[0])/2,(lim_pos[3]+lim_pos[2])/2,(lim_pos[5]+lim_pos[4])/2]
				current_mdl = bigpdb.shift_pdb(current_mdl,shift_arr)
				bigpdb.write_file(current_mdl,config.WATER_PATH+'/water_model_'+str(model_num)+'.pdb.gz')
			current_mdl = [line]
			model_num = int(line.strip().split(' ')[-1])
		else:
			current_mdl.append(line)
		line = prot_file.readline()
	current_mdl = bigpdb.rm_element(current_mdl,'H')
	lim_pos = bigpdb.get_lim_pos(current_mdl,line_type = 'ALL')
	shift_arr = [(lim_pos[1]+lim_pos[0])/2,(lim_pos[3]+lim_pos[2])/2,(lim_pos[5]+lim_pos[4])/2]
	current_mdl = bigpdb.shift_pdb(current_mdl,shift_arr)
	bigpdb.write_file(current_mdl,config.WATER_PATH+'/water_model_'+str(model_num)+'.pdb.gz')
	prot_file.close()
water_split()