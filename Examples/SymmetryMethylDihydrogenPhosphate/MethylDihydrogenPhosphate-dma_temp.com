%RWF=/scratch/bdw2292/Gau-MethylDihydrogenPhosphate/,100GB
%Nosave
%Chk=MethylDihydrogenPhosphate-dma_temp.chk
%Mem=20GB
%Nproc=4
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB

MethylDihydrogenPhosphate Gaussian SP Calculation on node37.bme.utexas.edu

0 1
 P    0.419543    0.035763   -0.003811
 O   -0.918183    0.900561   -0.109049
 O    0.243659   -1.121635   -1.107835
 O    0.241179   -0.827402    1.342439
 O    1.639149    0.847350   -0.100026
 C   -2.196321    0.267114   -0.034655
 H   -2.937952    1.060778   -0.130527
 H   -2.317562   -0.239791    0.927320
 H   -2.315476   -0.452611   -0.850007
 H    0.848570   -0.960707   -1.843517
 H    0.842713   -0.495243    2.021315

$nbo bndidx $end

