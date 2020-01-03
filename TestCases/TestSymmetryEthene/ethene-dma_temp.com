%RWF=/scratch/bdw2292/Gau-ethene/,20GB
%Nosave
%Chk=ethene-dma_temp.chk
%Mem=5GB
%Nproc=4
#P MP2/6-311G** Sp Density=MP2 MaxDisk=20GB

 ethene Gaussian Optimization on bme-nova.bme.utexas.edu

0 1
C	
C	
H	
H	
H	
H	
	
