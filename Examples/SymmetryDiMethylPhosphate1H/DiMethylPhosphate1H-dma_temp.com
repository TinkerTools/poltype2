%RWF=/scratch/bdw2292/Gau-DiMethylPhosphate1H/,100GB
%Nosave
%Chk=DiMethylPhosphate1H-dma_temp.chk
%Mem=20GB
%Nproc=4
#P MP2/6-311G** Sp Density=MP2 MaxDisk=100GB

DiMethylPhosphate1H Gaussian SP Calculation on node37.bme.utexas.edu

0 1
 P   -0.036356    0.468360    0.100405
 O    1.023256   -0.503527    0.789603
 O   -0.956225   -0.505850   -0.785920
 O    0.829093    1.190542   -1.044275
 O   -0.744946    1.357526    1.028579
 C    1.829940   -1.372383   -0.011804
 C   -1.939129   -1.284835   -0.095384
 H    1.202337   -2.014080   -0.638039
 H    2.407097   -1.981124    0.684392
 H    2.503641   -0.787759   -0.643698
 H   -2.664279   -0.631278    0.396654
 H   -1.465105   -1.933673    0.650339
 H   -2.432494   -1.897018   -0.852116
 H    0.466928    2.069592   -1.216351

$nbo bndidx $end

