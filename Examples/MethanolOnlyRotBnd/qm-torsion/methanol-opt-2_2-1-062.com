%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-062.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O   -0.756310    0.120220    0.000000
 C    0.676190   -0.021500    0.000000
 H    1.052950   -0.534560   -0.878820
 H    1.052910   -0.534590    0.878820
 H    1.078630    0.978910    0.000010
 H   -1.191180   -0.742530   -0.000020

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

