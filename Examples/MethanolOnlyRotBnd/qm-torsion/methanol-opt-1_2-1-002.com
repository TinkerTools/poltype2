%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-002.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749665    0.104560    0.003950
 C   -0.666402   -0.014690    0.002039
 H   -1.031841   -0.534705    0.901308
 H   -1.014891   -0.531603   -0.904427
 H   -1.120422    0.987207    0.002673
 H    1.093491   -0.329769    0.800657

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

