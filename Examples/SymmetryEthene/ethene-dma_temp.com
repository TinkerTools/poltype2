%RWF=/scratch/bdw2292/Gau-ethene/,100GB
%Nosave
%Chk=ethene-dma_temp.chk
%Mem=700MB
%Nproc=1
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB

ethene Gaussian SP Calculation on node37.bme.utexas.edu

0 1
 C   -0.664246   -0.000000    0.000000
 C    0.664246    0.000000    0.000000
 H   -1.232024   -0.926223    0.000000
 H   -1.232024    0.926223   -0.000000
 H    1.232024    0.926223    0.000000
 H    1.232024   -0.926223   -0.000000

$nbo bndidx $end

