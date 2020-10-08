%RWF=/scratch/bdw2292/Gau-methylamine/,20GB
%Nosave
%Chk=methylamine-dma_temp.chk
%Mem=5GB
%Nproc=4
#P MP2/6-311G** Sp Density=MP2 MaxDisk=20GB

methylamine Gaussian SP Calculation on node37.bme.utexas.edu

0 1
 N    0.749598    0.000001   -0.124965
 C   -0.703842    0.000004    0.017478
 H   -1.078538   -0.000156    1.054655
 H   -1.112986   -0.881109   -0.486486
 H   -1.112943    0.881261   -0.486248
 H    1.140166   -0.812079    0.343959
 H    1.140166    0.812056    0.344005

$nbo bndidx $end

