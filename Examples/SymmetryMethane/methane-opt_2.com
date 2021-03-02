%RWF=/scratch/bdw2292/Gau-methane/,100GB
%Nosave
%Chk=methane-opt_2.chk
%Mem=700MB
%Nproc=1
#P opt=(maxcycle=3,Loose) MP2/6-31G* MaxDisk=100GB 

methane Gaussian OPT Calculation on g2-node38.bme.utexas.edu

0 1
 C   -0.000084    0.000082    0.000095
 H   -0.676518    0.853273   -0.085717
 H   -0.394038   -0.836372   -0.581791
 H    0.086307   -0.293963    1.048557
 H    0.985247    0.276089   -0.382179

