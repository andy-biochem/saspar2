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
import bigpdb
import configuration as config
for fname in os.listdir(config.NEW_WATER_PATH):
	print('Shifting...'+fname)
	current_mdl = bigpdb.read_file(config.WATER_PATH+fname)
	current_mdl = bigpdb.rm_element(current_mdl,'H')
	lim_pos = bigpdb.get_lim_pos(current_mdl,line_type='ALL')
	shift_arr = [(lim_pos[1]+lim_pos[0])/2,(lim_pos[3]+lim_pos[2])/2,(lim_pos[5]+lim_pos[4])/2]
	current_mdl = bigpdb.shift_pdb(current_mdl,shift_arr)
	bigpdb.write_file(current_mdl,config.WATER_PATH+fname+'.gz')