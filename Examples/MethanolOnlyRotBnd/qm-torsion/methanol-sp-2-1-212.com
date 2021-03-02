%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-212.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O   -0.756022   -0.108141    0.016011
 C    0.660761    0.000677   -0.000469
 H    1.115720   -0.885646   -0.466822
 H    1.055991    0.148355    1.015716
 H    0.969220    0.873339   -0.596691
 H   -1.132542    0.783085   -0.059529

