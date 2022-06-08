%RWF=/scratch/bdw2292/Gau-CH3Cl/,100GB
%Nosave
%Chk=CH3Cl-opt_1.chk
%Mem=15GB
%Nproc=4
#P opt=(Cartesian,maxcycles=400) MP2/6-31G* MaxDisk=100GB 

CH3Cl Gaussian OPT Calculation on node37.bme.utexas.edu

0 1
Cl    1.654600    0.143200   -0.074400
 C   -0.122000    0.001700    0.018000
 H   -0.475400   -0.621700   -0.842400
 H   -0.438100   -0.524100    0.944600
 H   -0.619100    1.001000   -0.045800

