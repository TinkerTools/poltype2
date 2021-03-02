%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-122.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.757910   -0.122460    0.006570
 C   -0.679440    0.020910   -0.006380
 H   -1.103570   -0.517980   -0.841360
 H   -1.011350    1.052590   -0.035800
 H   -1.091220   -0.402810    0.901380
 H    1.219550    0.722480   -0.038480

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

