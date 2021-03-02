%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-212.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749119    0.107687    0.002221
 C   -0.666958   -0.010971   -0.000719
 H   -1.017446   -0.538752    0.898563
 H   -1.016639   -0.522474   -0.909937
 H   -1.132633    0.986662    0.008268
 H    0.987856    0.932149   -0.450497

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

