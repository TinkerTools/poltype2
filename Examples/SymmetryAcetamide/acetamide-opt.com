%RWF=/scratch/bdw2292/Gau-acetamide/,100GB
%Nosave
%Chk=acetamide-opt.chk
%Mem=700MB
%Nproc=1
#P opt=(maxcycle=3,Loose) MP2/6-31G* Guess=INDO MaxDisk=100GB

acetamide Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
 O    0.642000   -1.121000    0.044000
 N    0.750000    1.131000   -0.033000
 C   -1.434000    0.090000    0.003000
 C    0.062000   -0.044000   -0.001000
 H   -1.867000   -0.691000   -0.628000
 H   -1.801000   -0.016000    1.027000
 H   -1.745000    1.061000   -0.393000
 H    1.761000    1.082000   -0.004000
 H    0.305000    2.036000    0.006000

