%RWF=/scratch/bdw2292/Gau-MonoMethylPhosphate1H/,150GB
%Nosave
%Chk=MonoMethylPhosphate1H-opt.chk
%Mem=100GB
%Nproc=4
#P opt=(maxcycle=3,Loose) wB97XD/6-31G* Guess=INDO MaxDisk=150GB

MonoMethylPhosphate1H Gaussian OPT Calculation on node37.bme.utexas.edu

-1 1
 H    0.447000    1.202000   -1.651000
 O    0.378000    1.432000   -0.708000
 P    0.435000   -0.024000   -0.018000
 O   -0.943000    0.027000    0.827000
 O    1.637000   -0.069000    0.876000
 O    0.269000   -0.939000   -1.202000
 C   -2.145000   -0.092000    0.103000
 H   -2.177000    0.624000   -0.723000
 H   -2.981000    0.116000    0.776000
 H   -2.254000   -1.110000   -0.282000

