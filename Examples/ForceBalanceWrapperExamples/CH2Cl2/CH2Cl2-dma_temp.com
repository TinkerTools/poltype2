%RWF=/scratch/bdw2292/Gau-CH2Cl2/,100GB
%Nosave
%Chk=CH2Cl2-dma_temp.chk
%Mem=15GB
%Nproc=4
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB 

CH2Cl2 Gaussian SP Calculation on node37.bme.utexas.edu

0 1
Cl   -1.471689    0.195045   -0.026316
Cl    1.476423    0.159223    0.006006
 C   -0.010673   -0.795293    0.045352
 H   -0.027427   -1.357472    0.974726
 H   -0.009749   -1.465284   -0.810032



