%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-062.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O   -0.756918    0.104669   -0.000017
 C    0.659174   -0.009248   -0.000042
 H    1.016097   -0.526394   -0.903636
 H    1.016087   -0.526414    0.903666
 H    1.109857    0.994023    0.000039
 H   -1.131847   -0.790372   -0.000011

