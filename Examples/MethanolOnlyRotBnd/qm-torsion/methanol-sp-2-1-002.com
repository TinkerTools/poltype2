%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-002.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O   -0.754462   -0.115837    0.006267
 C    0.660782    0.011932   -0.005848
 H    0.972939    1.067656   -0.033015
 H    1.096179   -0.535955   -0.854650
 H    1.078195   -0.423752    0.913940
 H   -1.143600    0.771952   -0.039008

