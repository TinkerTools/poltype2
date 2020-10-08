%RWF=/scratch/bdw2292/Gau-MonoMethylPhosphate1H/,150GB
%Nosave
%Chk=MonoMethylPhosphate1H-esp.chk
%Mem=100GB
%Nproc=4
#P MP2/aug-cc-pVTZ Sp Density=MP2 SCF=Save MaxDisk=150GB Pop=NBORead

MonoMethylPhosphate1H Gaussian SP Calculation on node37.bme.utexas.edu

-1 1
 H   -0.661223    1.921661   -0.554905
 O   -0.547329    1.113913   -1.068151
 P   -0.495897   -0.076403    0.102728
 O    0.945526   -0.782244   -0.310171
 O   -1.539802   -1.106281   -0.156908
 O   -0.337719    0.666011    1.403674
 C    2.107585   -0.053092   -0.015782
 H    2.078270    0.954421   -0.456979
 H    2.958455   -0.596680   -0.448465
 H    2.252028    0.054017    1.066574

$nbo bndidx $end

