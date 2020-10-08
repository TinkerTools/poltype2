%RWF=/scratch/bdw2292/Gau-MonoMethylPhosphate0H/,150GB
%Nosave
%Chk=MonoMethylPhosphate0H-opt.chk
%Mem=100GB
%Nproc=8
#P opt=(maxcycle=3,Loose) wB97XD/6-31G* Guess=INDO MaxDisk=150GB SCRF=(PCM)

MonoMethylPhosphate0H Gaussian OPT Calculation on node37.bme.utexas.edu

-2 1
 P    0.488000    0.001000   -0.050000
 O   -0.913000   -0.004000    0.808000
 O    1.564000   -0.005000    1.025000
 O    0.406000    1.286000   -0.865000
 O    0.406000   -1.274000   -0.879000
 C   -2.084000   -0.000000    0.027000
 H   -2.137000    0.899000   -0.593000
 H   -2.947000   -0.004000    0.697000
 H   -2.136000   -0.892000   -0.603000

