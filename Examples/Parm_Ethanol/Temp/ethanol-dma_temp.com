%RWF=/scratch/liuchw/Gau-ethanol/,50GB
%Nosave
%Chk=ethanol-dma_temp.chk
%Mem=16GB
%Nproc=8
#P MP2/6-311G** Sp Density=MP2 MaxDisk=50GB 

ethanol Gaussian SP Calculation on bme-earth.bme.utexas.edu

0 1
 O    0.891446   -0.422961   -0.602533
 C    0.221490    0.403002    0.348589
 C   -1.245550    0.036625    0.300863
 H    0.352841    1.467475    0.108242
 H    0.618999    0.241487    1.360557
 H   -1.640995    0.209345   -0.702207
 H   -1.817813    0.637061    1.013725
 H   -1.374439   -1.019518    0.546242
 H    1.835754   -0.194927   -0.581956



