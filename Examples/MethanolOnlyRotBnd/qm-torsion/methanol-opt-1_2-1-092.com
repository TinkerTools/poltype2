%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-092.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749540    0.105070    0.001937
 C   -0.666253   -0.014562   -0.001859
 H   -1.019012   -0.538010    0.899336
 H   -1.024789   -0.528158   -0.907372
 H   -1.120848    0.987018    0.002766
 H    1.121062   -0.662059   -0.461909

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

