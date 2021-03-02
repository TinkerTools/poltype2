%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-182.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O   -0.757399    0.103718   -0.003935
 C    0.659608   -0.000253   -0.005427
 H    1.125453    0.996164   -0.021160
 H    1.013399   -0.603384   -0.855429
 H    0.999309   -0.503631    0.912356
 H   -1.127115   -0.792591    0.036653

