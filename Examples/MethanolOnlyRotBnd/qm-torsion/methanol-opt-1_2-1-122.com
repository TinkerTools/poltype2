%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-122.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749435    0.106251    0.000985
 C   -0.666513   -0.014265   -0.002711
 H   -1.017157   -0.538724    0.898503
 H   -1.028522   -0.527481   -0.907282
 H   -1.121852    0.987055    0.003121
 H    1.096310   -0.321437   -0.798016

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

