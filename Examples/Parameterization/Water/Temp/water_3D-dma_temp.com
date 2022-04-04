%RWF=/scratch/bdw2292/Gau-water_3D/,317GB
%Nosave
%Chk=water_3D-dma_temp.chk
%Mem=49GB
%Nproc=12
#P MP2/6-311G** Sp Density=MP2 MaxDisk=317GB 

water_3D Gaussian SP Calculation on node74.bme.utexas.edu

0 1
 O    0.001578    0.066727    0.000000
 H    0.750498   -0.547545    0.000000
 H   -0.775539   -0.511461    0.000000



