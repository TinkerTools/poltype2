%RWF=/scratch/bdw2292/Gau-ammonia/,20GB
%Nosave
%Chk=ammonia-opt.chk
%Mem=5GB
%Nproc=4
#P opt=(maxcycle=400) MP2/6-31G* Guess=INDO MaxDisk=20GB

ammonia Gaussian SP Calculation on bme-nova.bme.utexas.edu

0 1
 N    1.075300   -0.018600    0.047100
 H    0.763700    0.618300    0.779000
 H    0.763700    0.390600   -0.832600
 H    2.092600    0.039200    0.039000


