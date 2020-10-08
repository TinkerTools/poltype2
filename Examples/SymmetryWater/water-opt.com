%RWF=/scratch/bdw2292/Gau-water/,100GB
%Nosave
%Chk=water-opt.chk
%Mem=700MB
%Nproc=1
#P opt=(maxcycle=3,Loose) MP2/6-31G* Guess=INDO MaxDisk=100GB

water Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
 O   -0.005700    0.385200   -0.000000
 H   -0.796100   -0.194700   -0.000000
 H    0.801800   -0.190500    0.000000

