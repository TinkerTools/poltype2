%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-esp.chk
%Mem=100GB
%Nproc=8
#P MP2/aug-cc-pVTZ Sp Density=MP2 SCF=Save MaxDisk=200GB 
 MaxDisk=200GB 

methanol Gaussian SP Calculation on node37.bme.utexas.edu

0 1
 O    0.749531    0.122589   -0.000002
 C   -0.664447   -0.019643    0.000001
 H   -1.034255   -0.542785    0.891370
 H   -1.034254   -0.542822   -0.891347
 H   -1.076406    0.989749   -0.000020
 H    1.135352   -0.766996    0.000010



