%RWF=/scratch/liuchw/Gau-water_3D/,12672GB
%Nosave
%Chk=water_3D-dma_temp.chk
%Mem=209GB
%Nproc=51
#P MP2/6-311G** Sp Density=MP2 MaxDisk=12672GB 

water_3D Gaussian SP Calculation on node165.bme.utexas.edu

0 1
 O    0.000979    0.066556    0.000000
 H    0.755430   -0.540910    0.000000
 H   -0.770887   -0.518624    0.000000



