%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-032.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749674    0.104222    0.002724
 C   -0.666190   -0.014885    0.001018
 H   -1.028020   -0.535656    0.901120
 H   -1.016762   -0.530819   -0.905387
 H   -1.119919    0.987068    0.003068
 H    1.119817   -0.666630    0.461458

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

