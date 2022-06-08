%RWF=/scratch/bdw2292/Gau-CH3Cl/,100GB
%Nosave
%Chk=CH3Cl-dma_temp.chk
%Mem=15GB
%Nproc=4
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB 

CH3Cl Gaussian SP Calculation on node37.bme.utexas.edu

0 1
Cl    0.552658    0.050919   -0.021174
 C   -1.215710   -0.112024    0.046574
 H   -1.547587   -0.709546   -0.800067
 H   -1.489498   -0.602069    0.978726
 H   -1.663412    0.878730    0.001487



