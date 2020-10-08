%RWF=/scratch/bdw2292/Gau-water/,100GB
%Nosave
%Chk=water-dma_temp.chk
%Mem=700MB
%Nproc=1
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB

water Gaussian SP Calculation on node37.bme.utexas.edu

0 1
 O   -0.000305    0.066251    0.000000
 H   -0.757757   -0.529210    0.000000
 H    0.762593   -0.522243    0.000000

$nbo bndidx $end

