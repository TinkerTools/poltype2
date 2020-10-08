%RWF=/scratch/bdw2292/Gau-methylamine/,100GB
%Nosave
%Chk=methylamine-dma_temp.chk
%Mem=700MB
%Nproc=1
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB

methylamine Gaussian SP Calculation on node37.bme.utexas.edu

0 1
 N   -0.722156   -0.004423   -0.048588
 C    0.737742   -0.000901   -0.010587
 H    1.186936    0.089475    0.992767
 H    1.109350    0.828761   -0.620083
 H    1.109460   -0.926021   -0.462115
 H   -1.078071    0.848689    0.372537
 H   -1.077952   -0.768718    0.518054

$nbo bndidx $end

