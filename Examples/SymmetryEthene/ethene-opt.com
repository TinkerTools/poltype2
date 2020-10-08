%RWF=/scratch/bdw2292/Gau-ethene/,100GB
%Nosave
%Chk=ethene-opt.chk
%Mem=700MB
%Nproc=1
#P opt=(maxcycle=3,Loose) MP2/6-31G* Guess=INDO MaxDisk=100GB

ethene Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
 C   -0.668000    0.000000   -0.000000
 C    0.668000    0.000000   -0.000000
 H   -1.228000   -0.927000    0.071000
 H   -1.228000    0.927000   -0.071000
 H    1.228000    0.927000   -0.071000
 H    1.228000   -0.927000    0.071000

