## Installation

* Download POLTYPE
Make sure to download the ZIP folder from the poltype2 branch, NOT the master branch

```https://github.com/pren/poltype/tree/poltype2 ```

* Or simply using the following git command
```shell
git clone git@github.com:pren/poltype.git --branch poltype2
```

* Install anaconda 

https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Miniconda installer for Linux.(We also use Anaconda)
In your terminal window, run:
Miniconda:
```shell
bash Miniconda3-latest-Linux-x86_64.sh
```
You can install system wide or in your own directory (e.g /home/username/miniconda3/bin). If you do the later, make sure conda is in your PATH.

(For Anaconda: bash Anaconda-latest-Linux-x86_64.sh . We installed system wide conda: /opt/software/anaconda/)
After you have conda in your PATH, run the following commands. These will be installed under "/home/username/miniconda3/"

```shell
conda create -n p4env -c psi4/label/dev psi4 
conda activate p4env
conda install scipy 
conda install matplotlib 
conda install -c openbabel openbabel=2.4.0
conda install git pip
conda install -c rdkit rdkit
conda install -c conda-forge mdanalysis
conda install -c conda-forge svgutils
conda install -c conda-forge cairosvg
pip install git+git://github.com/ParmEd/ParmEd.git 
```
                   
### Environment var for running POLTYPE:
 Create environments for Tinker,Gaussian 09,and GDMA. Put these in a file and source it before running POLTYPE

```shell      
export PATH="/opt/software/anaconda/bin:$PATH"
conda activate p4env
export g09root=/opt/g09gh/gaussian
source $g09root/g09/bsd/g09.profile
export GAUSS_SCRDIR=/scratch/your_username/
export GDMADIR=/opt/gdma-2.2/bin/
export PATH=/opt/gdma-2.2/bin/gdma:$PATH
export PSI_SCRATCH=/scratch/your_username/ 
export PATH=/home/liuchw/Softwares/tinkers/Tinker-latest/source/:$PATH
```
