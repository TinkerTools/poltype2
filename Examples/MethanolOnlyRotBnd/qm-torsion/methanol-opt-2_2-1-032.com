%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-032.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.756990   -0.121490   -0.000190
 C   -0.677790    0.022040    0.002450
 H   -1.021460    0.916930   -0.505220
 H   -1.088990    0.016300    1.005010
 H   -1.084620   -0.826220   -0.527460
 H    1.205960    0.732650    0.014430

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

