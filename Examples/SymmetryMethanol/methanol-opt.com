%RWF=/scratch/bdw2292/Gau-methanol/,100GB
%Nosave
%Chk=methanol-opt.chk
%Mem=700MB
%Nproc=1
#P opt=(maxcycle=3,Loose) MP2/6-31G* Guess=INDO MaxDisk=100GB

methanol Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
 O    0.708000   -0.000000    0.000000
 C   -0.708000   -0.000000    0.000000
 H   -1.073000   -0.769000    0.685000
 H   -1.073000   -0.195000   -1.011000
 H   -1.063000    0.979000    0.331000
 H    0.994000   -0.880000   -0.298000

