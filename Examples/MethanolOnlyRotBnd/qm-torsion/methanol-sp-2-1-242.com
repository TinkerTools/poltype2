%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-242.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O   -0.754564   -0.114807    0.000044
 C    0.661750    0.002891   -0.000013
 H    1.092104   -0.452082   -0.904486
 H    1.092134   -0.451847    0.904578
 H    0.958870    1.063207   -0.000148
 H   -1.138344    0.776514   -0.000090
