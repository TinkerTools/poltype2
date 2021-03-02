%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-152.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.757670   -0.121400    0.005830
 C   -0.677880    0.016540   -0.007710
 H   -1.096940   -0.863890   -0.465600
 H   -1.021520    0.912410   -0.513900
 H   -1.076920    0.088510    1.000180
 H    1.201310    0.734920   -0.021010

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

