%RWF=/scratch/bdw2292/Gau-MethylDihydrogenPhosphate/,100GB
%Nosave
%Chk=MethylDihydrogenPhosphate-opt.chk
%Mem=20GB
%Nproc=4
#P opt=(maxcycle=3,Loose) MP2/6-31G* Guess=INDO MaxDisk=100GB

MethylDihydrogenPhosphate Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
 P    0.522000   -0.213000    0.025000
 O   -0.785000    0.698000   -0.086000
 O    0.390000   -1.333000   -1.103000
 O    0.388000   -1.033000    1.386000
 O    1.824000    0.508000   -0.061000
 C   -2.045000    0.064000   -0.012000
 H   -2.825000    0.825000   -0.101000
 H   -2.161000   -0.445000    0.949000
 H   -2.160000   -0.653000   -0.830000
 H    1.005000   -1.144000   -1.835000
 H    1.000000   -0.674000    2.054000

