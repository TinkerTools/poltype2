%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-152.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O    0.756018   -0.109548    0.010713
 C   -0.660178    0.005806   -0.009713
 H   -1.118161   -0.879865   -0.474313
 H   -0.983345    0.920488   -0.530540
 H   -1.042083    0.068975    1.020087
 H    1.135817    0.783077   -0.013946

