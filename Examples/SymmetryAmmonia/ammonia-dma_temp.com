%RWF=/scratch/bdw2292/Gau-ammonia/,100GB
%Nosave
%Chk=ammonia-dma_temp.chk
%Mem=700MB
%Nproc=1
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB

ammonia Gaussian SP Calculation on node37.bme.utexas.edu

0 1
 N   -0.000162    0.000373    0.068626
 H    0.918802   -0.203885   -0.314671
 H   -0.633383   -0.695552   -0.315878
 H   -0.283169    0.894261   -0.322965

$nbo bndidx $end

