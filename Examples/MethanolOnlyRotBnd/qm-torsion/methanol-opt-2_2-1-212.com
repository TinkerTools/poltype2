%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-212.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O   -0.758040   -0.120540    0.012410
 C    0.678000    0.010900    0.001610
 H    1.095370   -0.869680   -0.457380
 H    1.086180    0.165960    0.993440
 H    1.011010    0.867220   -0.580290
 H   -1.196250    0.735420   -0.064660

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

