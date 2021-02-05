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

import math,gzip

# read_file - reading pdb file
# input: file name
# output: cortage of strings
def read_file(file_name):
	if file_name.endswith('.gz'):
		fd = gzip.open(file_name,mode='rt')
	else:
		fd = open(file_name)
	lines = fd.readlines()
	fd.close() 
	return lines

# mod_num - counting models in pdb file 
# input: cortage of strings of pdb file;
# output: number of model in pdb file
def mod_num(pdb_lines):
	models = 0
	for line in pdb_lines:
		if line.startswith('MODEL'):
			models +=1
	return models
	
# mod_split - splitting multi-model pdb file into many single-model
# input: 1. cortage of strings of pdb file; 2. prefix of output files
# output: cortage of output files names
def mod_split(pdb_lines,prefix,first_num = 0):
	f_names = []
	header = []
	models = first_num
	now_write = 0
	for line in pdb_lines:
		if line.startswith('MODEL'):
			break
		header.append(line)
	for line in pdb_lines:
		if line.startswith('MODEL'):
			models +=1
			f_name = prefix+'model_'+str(models)+'.pdb'
			f_names.append(f_name)
			fd = open(f_name,'w+')
			for hline in header:
				fd.write(hline)
			fd.write('MODEL        1\n')
			now_write = 1
			continue
		if now_write == 1:
			fd.write(line)
		if line.startswith('ENDMDL'):
			now_write = 0
			fd.close()
			continue 
	return f_names

# to_arr - convert single pdb line to array
# input: pdb line
# output: four item cortage (x,y,z,element name)
def to_arr(pdb_line):
	coord_x = float(pdb_line[30:38])
	coord_y = float(pdb_line[38:46])
	coord_z = float(pdb_line[46:54])
	el_sym = pdb_line[76:78].strip()
	el_sym = el_sym.upper()
	return [coord_x,coord_y,coord_z,el_sym]

# rm_records - remove records of determined type from cortage of pdb lines;
# input: 1. cortage of strings of pdb file; 2. Type of record
# output: cortage of strings of pdb file without records of this type;
def rm_records(pdb_lines,rec_type):
	bad_lines = []
	for i in range(len(pdb_lines)):
		line = pdb_lines[i]
		if line.startswith(rec_type):
			bad_lines.append(i)
	bad_lines.sort(reverse=True)
	for i in bad_lines:
		del pdb_lines[i]
	return pdb_lines

# rm_element - remove ATOM or HETATM records of determined element from cortage of pdb lines;
# input: 1. cortage of strings of pdb file; 2. Element symbol
# output:  cortage of strings of pdb file without records with this element;
def rm_element(pdb_lines,el_sym):
	at = []
	bad_lines = []
	for i in range(len(pdb_lines)):
		line = pdb_lines[i]
		if line.startswith('ATOM') or line.startswith('HETATM'):
			atom_arr = to_arr(line)
			if not atom_arr[-1] in at:
				at.append(atom_arr[-1])
			if atom_arr[-1] == el_sym:
				bad_lines.append(i)
	bad_lines.sort(reverse=True)
	for i in bad_lines:
		del pdb_lines[i]
	return pdb_lines

# get_lim_pos - get minimum and maximum of coordinates
# input: cortage of strings of pdb file; record_type (ATOM or HETATM)
# output:  cortage of coordinates
def get_lim_pos(pdb_lines,line_type ='ATOM',bad_res=[]):
	prot_x = []
	prot_y = []
	prot_z = []
	for line in pdb_lines:
		if ((line.startswith('HETATM') and line_type == 'ALL') or line.startswith('ATOM')) and (not (get_residue(line) in bad_res)):
			atom_arr = to_arr(line)
			prot_x.append(atom_arr[0])
			prot_y.append(atom_arr[1])
			prot_z.append(atom_arr[2])
	return [min(prot_x),max(prot_x),min(prot_y),max(prot_y),min(prot_z),max(prot_z)]

def shift_pdb(pdb_lines,shift_arr):
	result_lines = []
	for line in pdb_lines:
		if (line.startswith('HETATM') or line.startswith('ATOM') ):
			atom_arr = to_arr(line)
			new_coord_arr = [atom_arr[0]-shift_arr[0],atom_arr[1]-shift_arr[1],atom_arr[2]-shift_arr[2]]
			result_lines.append(coord_upd(line,new_coord_arr))
		else:
			result_lines.append(line)
	return result_lines

def extract_box(pdb_lines, coords,group_id='' ):
	prot_min_x = coords[0]
	prot_max_x = coords[1]
	prot_min_y = coords[2]
	prot_max_y = coords[3]
	prot_min_z = coords[4]
	prot_max_z = coords[5]
	
	new_pdb = []
	for line in pdb_lines:
		if (line.startswith('HETATM') or line.startswith('ATOM') )and group_id in line:
			atom_arr = to_arr(line)
			if ((atom_arr[0] > prot_min_x) and (atom_arr[0] <prot_max_x) and (atom_arr[1] >prot_min_y) and (atom_arr[1]<prot_max_y) and (atom_arr[2] >prot_min_z) and (atom_arr[2] <prot_max_z)):
				new_pdb.append(line)
		else:
			new_pdb.append(line)
	info_lines = ["REMARK 250 FOLLOWING LABELS CREATED BY BIGPDB LIBRARY\n"]
	info_lines.append("REMARK 250 MIN_X="+str(prot_min_x)+" MAX_X="+str(prot_max_x)+" DIFFERENCE_X="+str(prot_max_x-prot_min_x)+"\n")
	info_lines.append("REMARK 250 MIN_Y="+str(prot_min_y)+" MAX_Y="+str(prot_max_y)+" DIFFERENCE_Y="+str(prot_max_y-prot_min_y)+"\n")
	info_lines.append("REMARK 250 MIN_Z="+str(prot_min_z)+" MAX_Z="+str(prot_max_z)+" DIFFERENCE_Z="+str(prot_max_z-prot_min_z)+"\n")
	info_lines.append("REMARK 250 TOTAL VOLUME="+str((prot_max_x-prot_min_x)*(prot_max_y-prot_min_y)*(prot_max_z-prot_min_z))+"\n")
	for i in range(len(new_pdb)):
		if new_pdb[i].startswith('TITLE'):
			new_pdb = new_pdb[:i]+info_lines+new_pdb[i:]
			break
	return new_pdb

def cut_sphere(pdb_lines, radius ,group_id='' ):
	bad_lines = []
	for i in range(len(pdb_lines)):
		line = pdb_lines[i]
		if (line.startswith('HETATM') or line.startswith('ATOM') )and group_id in line:
			atom_arr = to_arr(line)
			atom_radius = math.sqrt((atom_arr[1] * atom_arr[1]) +(atom_arr[2] * atom_arr[2]) +(atom_arr[0] * atom_arr[0]) )
			if (atom_radius>radius):
				bad_lines.append(i)
	bad_lines.sort(reverse=True)
	for i in bad_lines:
		del pdb_lines[i]
	info_lines = ["REMARK 250 FOLLOWING LABELS CREATED BY BIGPDB LIBRARY\n"]
	info_lines.append("REMARK 250 RADIUS="+str(radius)+"TOTAL VOLUME="+str((4*math.pi*radius*radius*radius)/3) +"\n")
	for i in range(len(pdb_lines)):
		if pdb_lines[i].startswith('TITLE'):
			pdb_lines = pdb_lines[:i]+info_lines+pdb_lines[i:]
			break
	return pdb_lines

# cut_offset_box - cut hetheroatoms farther from peptide atom than offset value
# input: 1. cortage of strings of pdb file; 2. offset value in angstroms; group id; 
# output:  cortage of strings of pdb file ;
def cut_offset_box(pdb_lines, offset,group_id=''):
	prot_min_max =  get_lim_pos(pdb_lines)
	prot_min_x = prot_min_max[0] - offset
	prot_max_x = prot_min_max[1] + offset
	prot_min_y = prot_min_max[2] - offset
	prot_max_y = prot_min_max[3] + offset
	prot_min_z = prot_min_max[4] - offset
	prot_max_z = prot_min_max[5] + offset
	
	bad_lines = []
	for i in range(len(pdb_lines)):
		line = pdb_lines[i]
		if line.startswith('HETATM') and group_id in line:
			atom_arr = to_arr(line)
			if not ((atom_arr[0] > prot_min_x) and (atom_arr[0] <prot_max_x) and (atom_arr[1] >prot_min_y) and (atom_arr[1]<prot_max_y) and (atom_arr[2] >prot_min_z) and (atom_arr[2] <prot_max_z)):
				bad_lines.append(i)
	bad_lines.sort(reverse=True)
	for i in bad_lines:
		del pdb_lines[i]
	info_lines = ["REMARK 250 FOLLOWING LABELS CREATED BY BIGPDB LIBRARY\n"]
	info_lines.append("REMARK 250 MIN_X="+str(prot_min_x)+" MAX_X="+str(prot_max_x)+" DIFFERENCE_X="+str(prot_max_x-prot_min_x)+"\n")
	info_lines.append("REMARK 250 MIN_Y="+str(prot_min_y)+" MAX_Y="+str(prot_max_y)+" DIFFERENCE_Y="+str(prot_max_y-prot_min_y)+"\n")
	info_lines.append("REMARK 250 MIN_Z="+str(prot_min_z)+" MAX_Z="+str(prot_max_z)+" DIFFERENCE_Z="+str(prot_max_z-prot_min_z)+"\n")
	info_lines.append("REMARK 250 TOTAL VOLUME="+str((prot_max_x-prot_min_x)*(prot_max_y-prot_min_y)*(prot_max_z-prot_min_z))+" OFFSET="+str(offset)+"\n")
	for i in range(len(pdb_lines)):
		if pdb_lines[i].startswith('TITLE'):
			pdb_lines = pdb_lines[:i]+info_lines+pdb_lines[i:]
			break
	return pdb_lines

# atom_num_upd - correct atom number
# input: pdb line, new atom number
# output: pdb line with new atom number
def atom_num_upd(line,atom_cnt):
	atom_num = str(atom_cnt)
	if atom_cnt <10:
		atom_num = ' '+atom_num
	if atom_cnt <100:
		atom_num = ' '+atom_num
	if atom_cnt <1000:
		atom_num = ' '+atom_num
	if atom_cnt <10000:
		atom_num = ' '+atom_num
	return line[0:6]+atom_num+line[11:]

# res_id_upd - correct residue id
# input: pdb line, new residue id
# output: pdb line with new residue id
def res_id_upd(line,atom_cnt):
	atom_num = str(atom_cnt)
	if atom_cnt <10:
		atom_num = ' '+atom_num
	if atom_cnt <100:
		atom_num = ' '+atom_num
	if atom_cnt <1000:
		atom_num = ' '+atom_num
	return line[0:22]+atom_num+line[26:]

# upd_num - update atoms numeration and residue ids
# input:  cortage of strings of pdb file; 
# output:  cortage of strings of pdb file with correct numeration;
def upd_num(pdb_lines):
	new_pdb_lines =[]
	atom_cnt = 0
	hetatm_cnt = 0

	for line in pdb_lines:
		if line.startswith('ATOM'):
			atom_cnt +=1
			line = atom_num_upd(line,atom_cnt)
			new_pdb_lines.append(line)
		elif line.startswith('HETATM'):
			atom_cnt +=1
			arr = to_arr(line)
			hetatm_cnt +=1
			line = atom_num_upd(line,atom_cnt)	
			line = res_id_upd(line,hetatm_cnt)
			new_pdb_lines.append(line)
		else:
			new_pdb_lines.append(line)
	return new_pdb_lines

# write_file - writing pdb file
# input: 1. pdb lines cortage; 2. file name
# output: none
def write_file(pdb_lines,file_name):
	if file_name.endswith('.gz'):
		fd = gzip.open(file_name,mode='wt+')
	else:
		fd =  open(file_name,'w+')
	for line in pdb_lines:
		fd.write(line)
	fd.close()

def coord_upd(line,coords):
	x = "%8.3f" % coords[0]
	y = "%8.3f" % coords[1]
	z = "%8.3f" % coords[2]
	line = line[:30]+x+y+z+line[54:]
	return line
def mean(a):
	return sum(a)/len(a)
	
def pdb_idx(fname,type='prot',ou='line'):
	pdb_str = read_file(fname)
	arr_x = []
	arr_y = []
	arr_z = []
	for line in pdb_str:
		if (type == 'prot' and line.startswith('ATOM')) or (type == 'het' and line.startswith('HETATM')):
			arr = to_arr(line)
			arr_x.append(arr[0])
			arr_y.append(arr[1])
			arr_z.append(arr[2])
	arr_x.sort()
	arr_y.sort()
	arr_z.sort()
	if ou =='line':
		line = fname+';'+str(arr_x[0])+';'+str(mean(arr_x))+';'+str(arr_x[-1])+';'+str(arr_y[0])+';'+str(mean(arr_y))+';'+';'+str(arr_y[-1])+';'+str(arr_z[0])+';'+str(mean(arr_z))+';'+';'+str(arr_z[-1])+';'+str((arr_x[-1]-arr_x[0])*(arr_y[-1]-arr_y[0])*(arr_z[-1]-arr_z[0]))+';\n'
		return line
	elif ou == 'arr':
		arr = [arr_x[0],arr_x[-1],arr_y[0],arr_y[-1],arr_z[0],arr_z[-1]]
		return arr

def find_farest(pdb_lines,distan=0):
	atom_lns = []
	dist = 0
	at_a = []
	at_b = []
	for line in pdb_lines:
		if line.startswith('ATOM'):
			atom_lns.append(line)
	for i in range(len(atom_lns)):
		for j in range(i,len(atom_lns)):
			at_c = to_arr(atom_lns[i])
			at_d = to_arr(atom_lns[j])
			now_dist = math.sqrt((at_c[0]-at_d[0])*(at_c[0]-at_d[0])+(at_c[1]-at_d[1])*(at_c[1]-at_d[1])+(at_c[2]-at_d[2])*(at_c[2]-at_d[2]))
			if now_dist > dist:
				dist = now_dist
				at_a = at_c[:]
				at_b  = at_d[:]
	if distan == 0:
		ret_arr = [(at_a[0]-at_b[0])/dist,(at_a[1]-at_b[1])/dist]
		return ret_arr
	elif distan == 1:
		return dist
	
def find_farest_XY(pdb_lines):
	atom_lns = []
	dist = 0
	at_a = []
	at_b = []
	for line in pdb_lines:
		if line.startswith('ATOM'):
			atom_lns.append(line)
	for i in range(len(atom_lns)):
		for j in range(i,len(atom_lns)):
			at_c = to_arr(atom_lns[i])
			at_d = to_arr(atom_lns[j])
			now_dist = math.sqrt((at_c[0]-at_d[0])*(at_c[0]-at_d[0])+(at_c[1]-at_d[1])*(at_c[1]-at_d[1]))
			if now_dist > dist:
				dist = now_dist
				at_a = at_c[:]
				at_b  = at_d[:]
	ret_arr = [(at_a[0]-at_b[0])/dist,(at_a[1]-at_b[1])/dist,(at_a[2]-at_b[2])/dist]
	return ret_arr

def get_residue(pdb_line):
	return pdb_line[17:20].strip()
