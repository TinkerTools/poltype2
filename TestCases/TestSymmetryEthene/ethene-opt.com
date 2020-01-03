%RWF=/scratch/bdw2292/Gau-ethene/,20GB
%Nosave
%Chk=ethene-opt.chk
%Mem=5GB
%Nproc=4
#P opt=(maxcycle=400) MP2/6-31G* Guess=INDO MaxDisk=20GB

ethene Gaussian SP Calculation on bme-nova.bme.utexas.edu

0 1
 C   -0.667200    0.000000    0.000000
 C    0.667200    0.000000    0.000000
 H   -1.221300   -0.929000    0.070800
 H   -1.221200    0.929000   -0.070800
 H    1.221300    0.929000   -0.070800
 H    1.221300   -0.929000    0.070800


