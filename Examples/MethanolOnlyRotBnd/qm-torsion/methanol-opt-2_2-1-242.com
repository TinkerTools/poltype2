%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-242.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O   -0.758650   -0.122140    0.000020
 C    0.679630    0.010990   -0.000010
 H    1.100710   -0.433480   -0.890240
 H    1.100770   -0.433210    0.890330
 H    1.001060    1.049260   -0.000120
 H   -1.211140    0.728610   -0.000080

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

