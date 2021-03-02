%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-332.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O   -0.757700   -0.121400    0.005830
 C    0.677910    0.016550   -0.007720
 H    1.021600    0.912400   -0.513900
 H    1.096970   -0.863910   -0.465560
 H    1.077090    0.088530    1.000120
 H   -1.201490    0.734860   -0.021000

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

