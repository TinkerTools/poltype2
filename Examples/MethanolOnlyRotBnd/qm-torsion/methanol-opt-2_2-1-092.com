%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-opt-2_2-1-092.chk
%Mem=100GB
%Nproc=8
#P opt=(ModRedundant,maxcycle=2,Loose) wB97XD/6-31G* MaxDisk=200GB 

methanol Rotatable Bond Optimization on node37.bme.utexas.edu

0 1
 O   -0.757010   -0.121480   -0.000170
 C    0.677820    0.022040    0.002450
 H    1.089140    0.016420    1.004960
 H    1.021510    0.916870   -0.505310
 H    1.084590   -0.826300   -0.527390
 H   -1.206060    0.732630    0.014380

3 2 1 6 F
4 2 1 6 F
5 2 1 6 F

