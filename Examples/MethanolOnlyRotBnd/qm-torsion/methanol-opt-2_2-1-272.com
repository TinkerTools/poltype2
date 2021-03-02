%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-272.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.758060   -0.120530    0.012400
 C   -0.677990    0.010900    0.001610
 H   -1.086310    0.165820    0.993390
 H   -1.095260   -0.869680   -0.457430
 H   -1.011130    0.867180   -0.580250
 H    1.196140    0.735490   -0.064600

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

