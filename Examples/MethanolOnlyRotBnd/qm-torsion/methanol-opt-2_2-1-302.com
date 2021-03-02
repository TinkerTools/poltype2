%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-302.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.757340    0.120110   -0.005070
 C   -0.676310   -0.012300   -0.005190
 H   -1.051430   -0.608020   -0.831160
 H   -1.090610    0.980420   -0.019740
 H   -1.041460   -0.513420    0.888840
 H    1.182670   -0.746100    0.033730

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

