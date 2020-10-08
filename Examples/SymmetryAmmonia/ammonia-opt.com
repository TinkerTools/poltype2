%RWF=/scratch/bdw2292/Gau-ammonia/,100GB
%Nosave
%Chk=ammonia-opt.chk
%Mem=700MB
%Nproc=1
#P opt=(maxcycle=3,Loose) MP2/6-31G* Guess=INDO MaxDisk=100GB

ammonia Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
 N   -0.004700    0.010500    0.239900
 H    0.951900   -0.223400   -0.078500
 H   -0.652300   -0.746700   -0.079000
 H   -0.294900    0.959600   -0.082400

