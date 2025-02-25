%RWF=/scratch/liuchw/Gau-ethanol/,50GB
%Nosave
%Chk=ethanol-opt_1.chk
%Mem=16GB
%Nproc=8
#P opt=(ModRedundant,maxcycles=400,LOOSE) MP2/6-31G* MaxDisk=50GB 

ethanol Gaussian OPT Calculation on bme-earth.bme.utexas.edu

0 1
 O    1.144800   -0.560000   -0.793100
 C    0.459900    0.246900    0.144500
 C   -1.016300   -0.115200    0.103100
 H    0.591200    1.311100   -0.101800
 H    0.858700    0.076900    1.155700
 H   -1.407400    0.055400   -0.895600
 H   -1.569000    0.490600    0.815200
 H   -1.142800   -1.165600    0.348500
 H    2.080900   -0.340100   -0.776400

9 1 2 4 F
4 2 3 6 F


