%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-002.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O   -0.757910   -0.122460    0.006600
 C    0.679450    0.020910   -0.006390
 H    1.011370    1.052580   -0.035970
 H    1.103510   -0.518190   -0.841270
 H    1.091310   -0.402630    0.901410
 H   -1.219580    0.722460   -0.038590

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

