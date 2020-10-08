%RWF=/scratch/bdw2292/Gau-DiMethylPhosphate0H/,100GB
%Nosave
%Chk=DiMethylPhosphate0H-opt.chk
%Mem=20GB
%Nproc=4
#P opt=(maxcycle=3,Loose) wB97XD/6-31G* Guess=INDO MaxDisk=100GB

DiMethylPhosphate0H Gaussian OPT Calculation on node37.bme.utexas.edu

-1 1
 P   -0.015000    0.453000   -0.040000
 O    0.978000   -0.557000    0.762000
 O   -0.986000   -0.654000   -0.735000
 O    0.773000    1.138000   -1.128000
 O   -0.819000    1.227000    0.974000
 C    1.947000   -1.223000   -0.017000
 C   -1.941000   -1.258000    0.108000
 H    1.490000   -1.699000   -0.889000
 H    2.415000   -1.998000    0.596000
 H    2.720000   -0.519000   -0.337000
 H   -2.730000   -0.542000    0.356000
 H   -1.475000   -1.633000    1.024000
 H   -2.391000   -2.101000   -0.422000

