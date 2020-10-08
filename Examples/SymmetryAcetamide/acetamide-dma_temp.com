%RWF=/scratch/bdw2292/Gau-acetamide/,100GB
%Nosave
%Chk=acetamide-dma_temp.chk
%Mem=700MB
%Nproc=1
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB

acetamide Gaussian SP Calculation on node37.bme.utexas.edu

0 1
 O    0.581361   -1.174787    0.048614
 N    0.778573    1.077397   -0.061957
 C   -1.454601    0.103784    0.005639
 C    0.053374   -0.079144    0.001739
 H   -1.904533   -0.698689   -0.581327
 H   -1.819284    0.010348    1.033939
 H   -1.774571    1.071370   -0.390694
 H    1.783065    1.008511   -0.012460
 H    0.355112    1.990065   -0.047997

$nbo bndidx $end

