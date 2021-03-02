%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-302.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O    0.757338    0.103672   -0.003958
 C   -0.659638   -0.000208   -0.005421
 H   -1.013346   -0.603400   -0.855430
 H   -1.125472    0.996159   -0.021195
 H   -0.999248   -0.503568    0.912447
 H    1.127110   -0.792620    0.036697

