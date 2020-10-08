%RWF=/scratch/bdw2292/Gau-DiMethylPhosphate1H/,100GB
%Nosave
%Chk=DiMethylPhosphate1H-opt.chk
%Mem=20GB
%Nproc=4
#P opt=(maxcycle=3,Loose) MP2/6-31G* Guess=INDO MaxDisk=100GB

DiMethylPhosphate1H Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
 P    0.069000    0.472000    0.156000
 O    1.134000   -0.494000    0.866000
 O   -0.890000   -0.480000   -0.710000
 O    0.920000    1.235000   -0.966000
 O   -0.638000    1.429000    1.058000
 C    1.922000   -1.338000    0.053000
 C   -1.855000   -1.220000    0.010000
 H    1.289000   -1.973000   -0.572000
 H    2.529000   -1.978000    0.699000
 H    2.591000   -0.744000   -0.576000
 H   -2.578000   -0.543000    0.473000
 H   -1.371000   -1.829000    0.778000
 H   -2.382000   -1.879000   -0.685000
 H    0.518000    2.104000   -1.147000

