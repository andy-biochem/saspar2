User manual for SASPAR2 software
1. Introduction
The SASPAR2 program is designed for the calculation of small-angle X-ray scattering curves of globular proteins using the molecular dynamics data (the calculation method for molecular dynamics is described in the Appendix of this document). The program was written by Andrei Semenov in 2020 (the author is grateful to Konstantin Mazokha for the assistance with developing the algorithm). If you are using this program, please reference this work in the following citation form <publication pending>.

2. System requirements
This program requires a graphics card supporting NVidia CUDA technology. The author has tested the work of this program using Nvidia GTX 1060 and NVidia RTX 3070, and it is most likely that the program will work with other graphics cards that are supported by PyCuda library. The author has achieved stable operation of the SASPAR2 program only for Linux. For correct work of the program, Python 3 version, as well as the installed libraries PyCuda, periodictable and Tkinter, are required.

3. Installation
3.1 Download CUDA package from the website https://developer.nvidia.com/cuda-zone and install this package
3.2 Install Tkinter according to the software documentation of the Linux in use
3.3 To install the required libraries, execute a console command pip3 install numpy, then pip3 install periodictable, and then pip3 install pycuda
3.4 To create a local copy of SASPAR2, it is possible to create a local copy of the repository using the command git clone https://github.com/andy-biochem/saspar2

4. Usage
The program can be used in two modes –– in interactive mode or in batch mode. To run the program in the interactive mode, execute the command python runme-dialog.py in the working directory of the program. After that, a number of questions will appear, aimed to save the settings in the file configuration.py (the names of the corresponding parameters are given in brackets)
4.1 Paths to the root directory of the program, to the input directories for water and protein, to the output directory and to the files containing input  water frames (in case of using multi-model files): (ROOT_PATH, PDBS_PATH, WATER_PATH, AMINO_PATH, OUTPUT_PATH,  NEW_WATER_PATH, respectively). After preprocessing, the water frames are saved; to use the previously saved water, agree with the default path and answer “yes” to the question of whether the water is preprocessed. Pdb files with frames can be packed with Gzip archiver: in this case, the presence of the extension .pdb.gz is important 
4.2 Whether to use all the frames or only some of the frames in the calculation (USE_EVERY_FRAME) (it is possible to use every second frame, every third frame, etc., to speed up the program, at the expense of accuracy).
4.3 The start value of scattering vector q (VECTOR_START), the difference between the closest values of the scattering vector (VECTOR_DIF) and the number of its values (VECTOR_COUNT). All the values are in inverse angstroms. The first value on the curve will be VECTOR_START + VECTOR_DIF.
4.4 The number of simultaneously calculated q values (MAX_PROC). It is recommended to use one thread per each 1.5 GB of video memory.
4.5 Number of averages in reciprocal space (N_Z and N_FI). These values can be changed to find a balance between speed and accuracy, but 28 is a good choice in most cases.
4.6 The value of the indentation from the protein and the edge of the cell, in angstroms (INDENT). The best value is 7, as shown in the work Fedorov, B. A., Smirnov, A. V., Yaroshenko, V. V. & Porozov, Y. B. (2019). Biophysics 64, 38–48.
To run the program in batch mode, manually set the corresponding parameters in the configuration.py file and run the command python3 runme.py 
After the program ends, the file with the calculated curve in the csv format can be found in the outputs directory.

Appendix: Calculating molecular dynamics
Molecular dynamics frames can be calculated using the program Gromacs, which is available to download from the website Gromacs.org and can be installed according to the instructions at this website. The distribution package of SASPAR includes catalogues gmx_water and gmx_prot, which contain configuration and command files for Gromacs. To calculate water frames, execute commands in gmx_water catalogue: first, sh runme.sh, then sh runme2.sh, and then sh runme3.sh It is recommended to use AMBER99SB-ILDN force field, and TIP4P water model (see the explanation in Chen, P. & Hub J. S. (2014). Biophys. J. 107, 435–447). To calculate protein frames, save the pdb file that contains the protein structure in the catalogue gmx_prot and write the name of this protein to the file protein.name. Then execute the same commands. The molecular dynamics results will be saved to the file with .pdb.gz extension.