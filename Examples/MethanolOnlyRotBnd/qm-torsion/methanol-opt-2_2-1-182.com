%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-182.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O   -0.757340    0.120110   -0.005060
 C    0.676310   -0.012300   -0.005180
 H    1.090570    0.980440   -0.019730
 H    1.051440   -0.607980   -0.831180
 H    1.041450   -0.513470    0.888820
 H   -1.182610   -0.746140    0.033690

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

