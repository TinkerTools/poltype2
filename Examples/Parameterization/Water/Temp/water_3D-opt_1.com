%RWF=/scratch/liuchw/Gau-water_3D/,12672GB
%Nosave
%Chk=water_3D-opt_1.chk
%Mem=209GB
%Nproc=51
#P opt=(ModRedundant,maxcycles=400,Loose) MP2/6-31G* MaxDisk=12672GB 

water_3D Gaussian OPT Calculation on node165.bme.utexas.edu

0 1
 O    0.005800    0.397800    0.000000
 H    0.760500   -0.210000    0.000000
 H   -0.766300   -0.187700    0.000000


