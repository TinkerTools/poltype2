%RWF=/scratch/bdw2292/Gau-ammonia/,20GB
%Nosave
%Chk=ammonia-dma_temp.chk
%Mem=5GB
%Nproc=4
#P MP2/6-311G** Sp Density=MP2 MaxDisk=20GB

 ammonia Gaussian Optimization on bme-nova.bme.utexas.edu

0 1
N	
H	
H	
H	
	
