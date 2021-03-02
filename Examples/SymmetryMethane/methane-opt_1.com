%RWF=/scratch/bdw2292/Gau-methane/,100GB
%Nosave
%Chk=methane-opt_1.chk
%Mem=700MB
%Nproc=1
#P opt=(maxcycle=3,Loose) wB97X-D/6-31G* MaxDisk=100GB 

methane Gaussian OPT Calculation on g2-node38.bme.utexas.edu

0 1
 C   -0.022100    0.003200    0.016500
 H   -0.669000    0.889400   -0.100900
 H   -0.377800   -0.857800   -0.588300
 H    0.096400   -0.315100    1.063800
 H    0.972500    0.280300   -0.391100

