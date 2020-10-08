%RWF=/scratch/bdw2292/Gau-methanol/,100GB
%Nosave
%Chk=methanol-dma_temp.chk
%Mem=700MB
%Nproc=1
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB

methanol Gaussian SP Calculation on node37.bme.utexas.edu

0 1
 O    0.690256    0.021078    0.007121
 C   -0.722885    0.029156    0.009880
 H   -1.147667   -0.713337    0.701390
 H   -1.147701   -0.140836   -0.990406
 H   -1.030822    1.021944    0.345847
 H    0.978599   -0.849443   -0.287482

$nbo bndidx $end

