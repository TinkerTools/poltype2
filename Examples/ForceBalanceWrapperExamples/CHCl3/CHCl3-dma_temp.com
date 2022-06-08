%RWF=/scratch/bdw2292/Gau-CHCl3/,100GB
%Nosave
%Chk=CHCl3-dma_temp.chk
%Mem=15GB
%Nproc=4
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB 

CHCl3 Gaussian SP Calculation on node37.bme.utexas.edu

0 1
Cl   -0.838973    1.458857   -0.042477
Cl    1.682209   -0.001232   -0.068257
Cl   -0.842947   -1.454004   -0.093933
 C   -0.000871   -0.008250    0.466083
 H    0.000354   -0.027410    1.551839



