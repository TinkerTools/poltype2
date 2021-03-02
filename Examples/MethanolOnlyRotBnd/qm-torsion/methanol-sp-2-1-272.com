%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-272.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O    0.756052   -0.108171    0.015999
 C   -0.660711    0.000672   -0.000457
 H   -1.056029    0.148223    1.015701
 H   -1.115714   -0.885571   -0.466936
 H   -0.969274    0.873396   -0.596565
 H    1.132565    0.783074   -0.059425

