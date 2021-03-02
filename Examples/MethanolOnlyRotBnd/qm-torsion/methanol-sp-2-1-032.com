%RWF=/scratch/bdw2292/Gau-methanol/,200GB
%Nosave
%Chk=methanol-sp-2-1-032.chk
%Mem=100GB
%Nproc=8
#P wB97XD/6-311+G* SP SCF=(qc,maxcycle=800) Pop=NBORead MaxDisk=200GB 

methanol Rotatable Bond SP Calculation on node37.bme.utexas.edu

0 1
 O    0.755772   -0.110481    0.004885
 C   -0.659873    0.011512   -0.000339
 H   -0.984489    0.923933   -0.524130
 H   -1.060161   -0.002827    1.024407
 H   -1.101149   -0.842058   -0.535430
 H    1.138390    0.781132    0.022505

