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
from controllers import protein_reader as reader
import csv_operations as csv_op
def min_max(arr):
	min_x = arr[0]
	max_x = arr[0]
	for item in arr:
		if item < min_x:
			min_x = item
		if item > max_x:
			max_x  = item
	return min_x,max_x
def border_check(pdb_arr,coord_arr):
	atom_coord_x = []
	atom_coord_y = []
	atom_coord_z = []
	for loc_line in pdb_arr:
		if (loc_line.startswith('HETATM') or loc_line.startswith('ATOM')):
			loc_arr = bigpdb.to_arr(loc_line)
			atom_coord_x.append(loc_arr[0])
			atom_coord_y.append(loc_arr[1])
			atom_coord_z.append(loc_arr[2])
	min_x,max_x = min_max(atom_coord_x)
	min_y,max_y = min_max(atom_coord_y)
	min_z,max_z = min_max(atom_coord_z)
	if coord_arr[0] > min_x and coord_arr[1] < max_x and coord_arr[2] > min_y and coord_arr[3] < max_y and coord_arr[4] > min_z and coord_arr[5] < max_z:
		return True
	else:
		return False
def trim():
	atom_coord_x = []
	atom_coord_y = []
	atom_coord_z = []
	non_prot_groups = ['SOL','NA','CL','A','L']
	if os.path.exists(config.ROOT_PATH +'/dimensions.dat'):
		fd = open(config.ROOT_PATH +'/dimensions.dat')
		lines = fd.readlines()
		fd.close()
		dist_x = float(lines[0])
		dist_y = float(lines[1])
		dist_z = float(lines[2])
		border_min_x = float(lines[3])
		border_max_x = float(lines[4])
		border_min_y = float(lines[5])
		border_max_y = float(lines[6])
		border_min_z = float(lines[7])
		border_max_z = float(lines[8])
	else:
		if config.MULIT_MODEL_FILE.upper() == 'YES':
			if os.path.isdir(config.AMINO_PATH):
				fname =config.AMINO_PATH+os.listdir(config.AMINO_PATH)[0]
			else:
				fname = config.AMINO_PATH
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
						lim_pos = bigpdb.get_lim_pos(current_mdl,bad_res = non_prot_groups)
						shift_arr = [(lim_pos[1]+lim_pos[0])/2,(lim_pos[3]+lim_pos[2])/2,(lim_pos[5]+lim_pos[4])/2]
						current_mdl = bigpdb.shift_pdb(current_mdl,shift_arr)
						for loc_line in current_mdl:
							if loc_line.startswith('ATOM'):
								loc_arr = bigpdb.to_arr(loc_line)
								if bigpdb.get_residue(loc_line) not in non_prot_groups:
									atom_coord_x.append(loc_arr[0])
									atom_coord_y.append(loc_arr[1])
									atom_coord_z.append(loc_arr[2])
						bigpdb.write_file(current_mdl,config.TEMP_PATH+'/model_'+str(model_num)+'.pdb.gz')
					current_mdl = [line]
					model_num = int(line.strip().split(' ')[-1])
				else:
					current_mdl.append(line)
				line = prot_file.readline()
			current_mdl = bigpdb.rm_element(current_mdl,'H')
			lim_pos = bigpdb.get_lim_pos(current_mdl,bad_res = non_prot_groups)
			shift_arr = [(lim_pos[1]+lim_pos[0])/2,(lim_pos[3]+lim_pos[2])/2,(lim_pos[5]+lim_pos[4])/2]
			current_mdl = bigpdb.shift_pdb(current_mdl,shift_arr)
			for loc_line in current_mdl:
				if loc_line.startswith('ATOM'):
					loc_arr = bigpdb.to_arr(loc_line)
					if bigpdb.get_residue(loc_line) not in non_prot_groups:
						atom_coord_x.append(loc_arr[0])
						atom_coord_y.append(loc_arr[1])
						atom_coord_z.append(loc_arr[2])
			bigpdb.write_file(current_mdl,config.TEMP_PATH+'/model_'+str(model_num)+'.pdb.gz')
			prot_file.close()
		else:
			for fname in os.listdir(config.AMINO_PATH):
				print('Shifting...'+fname)
				current_mdl = bigpdb.read_file(config.AMINO_PATH+fname)
				current_mdl = bigpdb.rm_element(current_mdl,'H')
				lim_pos = bigpdb.get_lim_pos(current_mdl,bad_res = non_prot_groups)
				shift_arr = [(lim_pos[1]+lim_pos[0])/2,(lim_pos[3]+lim_pos[2])/2,(lim_pos[5]+lim_pos[4])/2]
				current_mdl = bigpdb.shift_pdb(current_mdl,shift_arr)
				for loc_line in current_mdl:
					if loc_line.startswith('ATOM'):
						loc_arr = bigpdb.to_arr(loc_line)
						if bigpdb.get_residue(loc_line) not in non_prot_groups:
							atom_coord_x.append(loc_arr[0])
							atom_coord_y.append(loc_arr[1])
							atom_coord_z.append(loc_arr[2])
				bigpdb.write_file(current_mdl,config.TEMP_PATH+fname+'.gz')
		min_x,max_x = min_max(atom_coord_x)
		min_y,max_y = min_max(atom_coord_y)
		min_z,max_z = min_max(atom_coord_z)
		border_min_x = min_x - config.INDENT
		border_min_y = min_y - config.INDENT
		border_min_z = min_z - config.INDENT
		border_max_x = max_x + config.INDENT
		border_max_y = max_y + config.INDENT
		border_max_z = max_z + config.INDENT
		dist_x = border_max_x - border_min_x
		dist_y = border_max_y - border_min_y
		dist_z = border_max_z - border_min_z
		fd = open(config.ROOT_PATH +'/dimensions.dat','w+')
		fd.write("%8.3f" % dist_x)
		fd.write('\n')
		fd.write("%8.3f" % dist_y)
		fd.write('\n')
		fd.write("%8.3f" % dist_z)
		fd.write('\n')
		fd.write("%8.3f" % border_min_x)
		fd.write('\n')
		fd.write("%8.3f" % border_max_x)
		fd.write('\n')
		fd.write("%8.3f" % border_min_y)
		fd.write('\n')
		fd.write("%8.3f" % border_max_y)
		fd.write('\n')
		fd.write("%8.3f" % border_min_z)
		fd.write('\n')
		fd.write("%8.3f" % border_max_z)
		fd.write('\n')
		fd.close()
	amino_dict = reader.load_amino_files_dict(config.AMINOACID_PATH)
	if not os.path.exists(config.ROOT_PATH+'/proteins.csv'):
		proteins = []
		pdb_lst = open('pdb_lst.txt','w+')
		flist = os.listdir(config.TEMP_PATH)
		flist.sort()
		flist = flist[::config.USE_EVERY_FRAME]
		num = 1
		for file in flist:
			if file.endswith(".pdb.gz"):
				now_mod = bigpdb.read_file(config.TEMP_PATH + file)
				if border_check(now_mod,[border_min_x - config.INDENT, border_max_x + config.INDENT, border_min_y - config.INDENT, border_max_y + config.INDENT, border_min_z - config.INDENT, border_max_z + config.INDENT]):
					pdb_lst.write(file+'...is OK\n')
					print(file+'...is OK')
					new_mod = bigpdb.extract_box(now_mod,[border_min_x, border_max_x  , border_min_y, border_max_y  , border_min_z, border_max_z  ])
					protein_pdb = reader.load_pdb_data(new_mod)
					protein = reader.get_protein_data(protein_pdb,amino_dict)
					num = csv_op.csv_append(config.ROOT_PATH+'/proteins.csv',protein,num)
				else:
					pdb_lst.write(file+'...is out of border\n')
					print(file+'...is out of border')
			
		pdb_lst.close()
		csv_op.csv_final(config.ROOT_PATH+'/proteins.csv')
	if not os.path.exists(config.ROOT_PATH+'/waters.csv'):
		flist = os.listdir(config.WATER_PATH)
		flist.sort()
		num = 1
		for file in flist:
			if file.endswith(".pdb") or file.endswith(".pdb.gz"):
				now_mod = bigpdb.read_file(config.WATER_PATH+file)
				if border_check(now_mod,[border_min_x - config.INDENT, border_max_x + config.INDENT, border_min_y - config.INDENT, border_max_y + config.INDENT, border_min_z - config.INDENT, border_max_z + config.INDENT]):
					print(file+'...is OK')
					new_mod = bigpdb.extract_box(now_mod,[border_min_x, border_max_x  , border_min_y, border_max_y  , border_min_z, border_max_z  ])
					water_pdb = reader.load_pdb_data(new_mod)
					water = reader.get_protein_data(water_pdb,amino_dict)
					num = csv_op.csv_append(config.ROOT_PATH+'/waters.csv',water,num)
				else:
					print(file+'...is out of border')
		csv_op.csv_final(config.ROOT_PATH+'/waters.csv')