# This file is part of SASPAR2 software package.
# Copyright 2020 Andrei M. Semenov
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#  
#	 http://www.apache.org/licenses/LICENSE-2.0
#  
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import configuration as config
import os,sys
import tkinter as tk
from tkinter import filedialog

def question(question,answer):
	print(question+'Enter for ['+answer+']')
	s = input()
	if s == '':
		return answer
	else:
		return s
def question_path(question,answer,obj = 'dir'):
	print(question+'Enter for ['+answer+'], E for open dialog')
	s = input()
	if s == '':
		return answer
	elif s.upper() == 'E':
		if obj == 'file':
			file_path = filedialog.askopenfilename()
			return file_path
		elif obj == 'dir':
			file_path = filedialog.askdirectory()
			return file_path
	else:
		return s
def question_yn(question,answer):
	print(question+'Enter for ['+answer+']')
	s = input()
	while  s.upper() not in ['YES','NO','']:
		s = input()
	if s == '':
		return answer
	else:
		return s.upper()
def interactivity():
	config.MULIT_MODEL_FILE = question_yn('Protein file is multi model? ',config.MULIT_MODEL_FILE)
	config.ROOT_PATH = question_path('Current working path: ',os.getcwd())
	config.OUTPUT_PATH = question_path('Output files path: ',config.ROOT_PATH + "/outputs")
	config.PDBS_PATH = question_path('Input files path:',config.ROOT_PATH + "/pdbs/")
	if config.MULIT_MODEL_FILE == 'YES':
		obj = 'file'
	else:
		obj = 'dir'
	config.WATER_PATH = question_path('Path to water PDBs: ',config.PDBS_PATH + "waters/",obj)
	config.AMINO_PATH = question_path('Path to protein PDBs: ',config.PDBS_PATH + "aminos/",obj)
	config.TEMP_PATH = question_path('Temporary files path: ',config.PDBS_PATH + "temp/")
	config.USE_EVERY_FRAME = int(question('Use every frame:',str(config.USE_EVERY_FRAME)))
	config.VECTOR_START = float(question('Start value for q: ',str(config.VECTOR_START)))
	config.VECTOR_COUNT = int(question('Number of q points: ',str(config.VECTOR_COUNT)))
	config.VECTOR_DIF = float(question('Difference between nearest q: ',str(config.VECTOR_DIF)))
	config.MAX_PROC = int(question('Number of q calculated parallelly: ',str(config.MAX_PROC)))
	config.N_Z = int(question('Number of vectors for Z: ',str(config.N_Z)))
	config.N_FI = int(question('Number of vectors for Phi: ',str(config.N_FI)))
	config.INDENT = int(question('Value of indentation: ',str(config.INDENT)))
	if question_yn('Water is preprocessed? ','NO') == 'NO':
		if question_yn('Water file is multi model? ',config.MULIT_MODEL_FILE) == 'YES':
			config.NEW_WATER_PATH = question_path('Path to new water file','','file')
			water_cmd = 'python3 water_split.py'
		else:
			config.NEW_WATER_PATH = question_path('Path to new water file','','dir')
			water_cmd = 'python3 water_shift.py'
	else:
		water_cmd = 'NULL'
	fd = open('configuration.py','w+')
	fd.write('ROOT_PATH = \''+config.ROOT_PATH+'\'\n')
	fd.write('AMINOACID_PATH = \''+ config.ROOT_PATH + '/aminoacids/atoms\'\n')
	fd.write('GROUP_DESCR = \''+ config.ROOT_PATH + '/aminoacids/data/atomic_groups.txt\'\n')
	fd.write('OUTPUT_PATH = \''+config.OUTPUT_PATH+'\'\n')
	fd.write('PDBS_PATH = \''+config.PDBS_PATH+'\'\n')
	fd.write('WATER_PATH = \''+config.WATER_PATH+'\'\n')
	fd.write('NEW_WATER_PATH = \''+config.NEW_WATER_PATH+'\'\n')
	fd.write('AMINO_PATH = \''+config.AMINO_PATH+'\'\n')
	fd.write('TEMP_PATH = \''+config.TEMP_PATH+'\'\n')
	fd.write('USE_EVERY_FRAME = '+str(config.USE_EVERY_FRAME)+'\n')
	fd.write('VECTOR_START = '+str(config.VECTOR_START)+'\n')
	fd.write('VECTOR_COUNT= '+str(config.VECTOR_COUNT)+'\n')
	fd.write('VECTOR_DIF = '+str(config.VECTOR_DIF)+'\n')
	fd.write('MAX_PROC = '+str(config.MAX_PROC)+'\n')
	fd.write('N_Z = '+str(config.N_Z)+'\n')
	fd.write('N_FI = '+str(config.N_FI)+'\n')
	fd.write('MULIT_MODEL_FILE = \''+config.MULIT_MODEL_FILE+'\'\n')
	fd.write('INDENT = '+str(config.INDENT)+'\n')
	fd.close()
	if not water_cmd == 'NULL':
		os.system(water_cmd)
root = tk.Tk()
root.withdraw()
if not sys.argv[-1] == '-a':
	interactivity()
try:
	import pycuda
except:
	print('PyCUDA is not installed. Try to install it by typing "pip3 install pycuda" in command prompt and press Enter')
	exit()
os.system('python3 runme.py')
