%RWF=/scratch/bdw2292/Gau-methylamine/,20GB
%Nosave
%Chk=methylamine-opt.chk
%Mem=5GB
%Nproc=4
#P opt=(maxcycle=3,Loose) wB97XD/6-31G* Guess=INDO MaxDisk=20GB

methylamine Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
 N   -0.727000   -0.000000   -0.003000
 C    0.725000   -0.000000   -0.002000
 H    1.109000    0.092000    1.018000
 H    1.103000    0.837000   -0.595000
 H    1.103000   -0.929000   -0.436000
 H   -1.066000    0.857000    0.432000
 H   -1.066000   -0.766000    0.578000

