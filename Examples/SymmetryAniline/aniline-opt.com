%RWF=/scratch/bdw2292/Gau-aniline/,100GB
%Nosave
%Chk=aniline-opt.chk
%Mem=700MB
%Nproc=1
#P opt=(maxcycle=3,Loose) MP2/6-31G* Guess=INDO MaxDisk=100GB

aniline Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
 N   -2.409000   -0.000000    0.179000
 C   -1.020000   -0.000000    0.005000
 C   -0.307000    1.201000    0.033000
 C   -0.307000   -1.202000    0.035000
 C    1.090000    1.204000   -0.011000
 C    1.090000   -1.205000   -0.008000
 C    1.789000   -0.000000   -0.040000
 H   -0.835000    2.149000    0.090000
 H   -0.834000   -2.149000    0.095000
 H    1.630000    2.147000   -0.006000
 H    1.630000   -2.147000   -0.001000
 H    2.875000   -0.000000   -0.066000
 H   -2.861000   -0.843000   -0.160000
 H   -2.861000    0.842000   -0.162000
