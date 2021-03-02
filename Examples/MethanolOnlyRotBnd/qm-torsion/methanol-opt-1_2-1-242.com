%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-1_2-1-242.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) HF/3-21G MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O    0.749045    0.107751    0.004030
 C   -0.667230   -0.010480    0.000088
 H   -1.018523   -0.539402    0.898509
 H   -1.014676   -0.521054   -0.910287
 H   -1.136109    0.985832    0.009241
 H    0.979992    1.050253    0.014019

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

