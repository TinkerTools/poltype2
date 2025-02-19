%RWF=/scratch/liuchw/Gau-35_38_Index_0-farH_in17/,1267GB
%Nosave
%Chk=35_38_Index_0-opt_1.chk
%Mem=20GB
%Nproc=5
#P opt=(ModRedundant,maxcycles=400,Loose) MP2/6-31G* MaxDisk=1267GB 

35_38_Index_0 Gaussian OPT Calculation on node165.bme.utexas.edu

0 1
 H  -11.047200    0.233600    0.388200
 C  -10.636100    0.611800    1.329500
 O   -9.386700    1.266900    1.102500
 C   -8.427800    0.459100    0.585700
 H  -10.517400   -0.198700    2.055500
 H  -11.333400    1.346400    1.741500
 O   -8.557700   -0.725300    0.321900
 H   -7.445300    0.911800    0.386500


