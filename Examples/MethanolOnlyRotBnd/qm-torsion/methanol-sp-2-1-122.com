%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-122.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O    0.754461   -0.115843    0.006262
 C   -0.660793    0.011942   -0.005852
 H   -1.096172   -0.535900   -0.854693
 H   -0.972931    1.067669   -0.032975
 H   -1.078205   -0.423790    0.913903
 H    1.143596    0.771949   -0.038964

