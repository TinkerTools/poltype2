%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-332.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O   -0.756054   -0.109589    0.010753
 C    0.660090    0.005787   -0.009750
 H    0.983372    0.920448   -0.530538
 H    1.118235   -0.879830   -0.474315
 H    1.041949    0.068959    1.020089
 H   -1.135693    0.783138   -0.013881

