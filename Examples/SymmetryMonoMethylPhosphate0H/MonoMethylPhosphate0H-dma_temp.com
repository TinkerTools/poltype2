%RWF=/scratch/bdw2292/Gau-MonoMethylPhosphate0H/,150GB
%Nosave
%Chk=MonoMethylPhosphate0H-dma_temp.chk
%Mem=100GB
%Nproc=8
#P MP2/6-311G** Sp Density=MP2 MaxDisk=150GB

MonoMethylPhosphate0H Gaussian SP Calculation on node37.bme.utexas.edu

-2 1
 P   -0.521330   -0.056132   -0.000372
 O    0.949739    0.803804    0.002597
 O   -1.517104    1.105985   -0.021027
 O   -0.496386   -0.871477    1.299332
 O   -0.475126   -0.900691   -1.280814
 C    2.127013    0.039371    0.000850
 H    2.185060   -0.619830    0.880044
 H    2.988589    0.719826    0.015292
 H    2.195233   -0.595211   -0.895561

$nbo bndidx $end

