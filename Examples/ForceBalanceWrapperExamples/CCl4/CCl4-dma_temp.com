%RWF=/scratch/bdw2292/Gau-CCl4/,100GB
%Nosave
%Chk=CCl4-dma_temp.chk
%Mem=15GB
%Nproc=4
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB 

CCl4 Gaussian SP Calculation on node37.bme.utexas.edu

0 1
Cl   -0.900181   -0.507459    1.435987
Cl   -0.713579    1.482673   -0.649813
Cl   -0.074096   -1.275691   -1.223782
Cl    1.687885    0.300450    0.437515
 C   -0.000087    0.000080    0.000272



