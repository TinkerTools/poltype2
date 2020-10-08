%RWF=/scratch/bdw2292/Gau-DiMethylPhosphate0H/,100GB
%Nosave
%Chk=DiMethylPhosphate0H-esp.chk
%Mem=20GB
%Nproc=4
#P MP2/aug-cc-pVTZ Sp Density=MP2 SCF=Save MaxDisk=100GB Pop=NBORead

DiMethylPhosphate0H Gaussian SP Calculation on node37.bme.utexas.edu

-1 1
 P   -0.000031   -0.612185    0.000044
 O   -0.949996    0.478239    0.834052
 O    0.949973    0.478261   -0.834088
 O   -0.831467   -1.287115   -1.036970
 O    0.831619   -1.286979    1.036967
 C   -1.909019    1.164295    0.072138
 C    1.908943    1.164323   -0.072171
 H   -1.451625    1.685437   -0.781417
 H   -2.384729    1.908984    0.724723
 H   -2.671880    0.481335   -0.322283
 H    2.672440    0.481620    0.321417
 H    1.451710    1.684672    0.781922
 H    2.383986    1.909764   -0.724508

$nbo bndidx $end

