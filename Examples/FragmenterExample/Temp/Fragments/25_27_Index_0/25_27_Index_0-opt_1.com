%RWF=/scratch/liuchw/Gau-25_27_Index_0-farH_in17/,1267GB
%Nosave
%Chk=25_27_Index_0-opt_1.chk
%Mem=20GB
%Nproc=5
#P opt=(ModRedundant,maxcycles=400,Loose) MP2/6-31G* MaxDisk=1267GB 

25_27_Index_0 Gaussian OPT Calculation on node165.bme.utexas.edu

0 1
 N    2.967700   -1.225900    0.938300
 C    4.207900   -0.750700    0.463500
 C    5.417000   -0.912000    1.425500
 N    6.609300   -0.198000    0.919400
 O    4.378600   -0.195400   -0.626700
 H    5.638500   -1.984900    1.468800
 H    5.161100   -0.544900    2.425100
 H    2.126700   -1.148600    0.348700
 H    2.903100   -1.650200    1.874600
 H    7.171200    0.138500    1.714400
 H    7.175600   -0.839100    0.345700


