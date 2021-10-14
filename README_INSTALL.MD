## Installation


### Required Programs
* Tinker, version >=8.9.4 (potential fitting errors previous versions)
* GDMA, version 2.3 and greater
* PSI4 (or Gaussian 09/16)
* scipy
* matplotlib
* openbabel, version 2.4
* rdkit
* svgutils
* cairosvg


### Tinker CPU Installation
```
git clone git@github.com:TinkerTools/Tinker.git tinker
cd tinker
cd fftw
./configure --prefix=/path_to_tinkerfolder/fftw/ --enable-openmp --enable-threads
make -j 8
make install
cd ..
cp ./make/Makefile ./source
cd source
```
* Open Makefile and comment out default operating system (MAC, unless you have mac)
* Uncomment the default linux installation lines
* Keep the RENAME line commented out

* Download GDMA-2.3 and greater
* Ren lab users can just add the following lines to your poltype .bashrc
```shell
export GDMADIR=/opt/gdma/gdma-2.3.3/bin/
export PATH=/opt/gdma/gdma-2.3.3/bin/:$PATH
```

* If you dont have a GitHub account, make one
* Click the watch button on https://github.com/pren/poltype for important updates


* Download POLTYPE

* Simply using the following git command
```shell
git clone https://github.com/TinkerTools/poltype2.git master
```
* Make sure you are on a node with an updated OS (as much as possible)


* Install anaconda 

https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Miniconda installer for Linux. (We also use Anaconda)
In your terminal window, run:
Miniconda:
```shell
bash Miniconda3-latest-Linux-x86_64.sh
```
You can install system wide or in your own directory (e.g /home/username/miniconda3/bin). If you do the later, make sure conda is in your PATH.

(For Anaconda: `bash Anaconda-latest-Linux-x86_64.sh`. We installed system wide conda: /opt/software/anaconda/)

You may need to run "conda init bash" and then source ~/.bashrc first.

After you have conda in your PATH, run the following commands. These will be installed under "/home/username/miniconda3/". You can use the provided `environment.yml` file to create a conda environment named "poltype" with `conda env create -f environment.yml`, or run the following commands:

```shell
conda create -n poltype python=3.6 psi4=1.3.2+ecbda83 -c psi4 --yes
conda activate poltype
conda install scipy --yes
conda install matplotlib --yes
conda install -c conda-forge openbabel=2.4.1 --yes
conda install git pip --yes
conda install -c conda-forge rdkit --yes
conda install -c conda-forge mdanalysis --yes
conda install -c conda-forge svgutils --yes
conda install -c conda-forge cairosvg --yes
```

### Environment var for running POLTYPE:
 Create environments for Tinker,Gaussian 09,and GDMA. Put these in a file and source it before running POLTYPE

* Note do not put source in your default .bashrc profile
* Note do not keep multiple sources for different purposes (unless necessary) in the same .bashrc source file, this can cause internal conflicts to arise

* Generic bashrc
```shell
#!/usr/bin/bash
#
##Python
##Set whichever conda you are using in next line
export PATH="path_to_conda/bin:$PATH"

#Poltype
conda activate poltype
# you may add alias for submitting - replace path_to_poltype
#alias run_poltype="nohup python path_to_poltype/poltype.py &"
#
#Gaussian, GDMA and scratch
#replace paths with your own
export g09root=path_to_g09/gaussian
source $g09root/g09/path_to_g09_profile/g09.profile
export GAUSS_SCRDIR=/scratch/$USER
export GDMADIR=/path_to_gdma/bin
export PATH=$GDMADIR:$PATH
export PSI_SCRATCH=/scratch/$USER

#Tinker
export PATH=path_to_tinker/source/:$PATH

```


* Ren lab users
```shell      
#Set whichever conda you are using in next line
export PATH="/opt/software/anaconda/bin:$PATH"
conda activate poltype
export g09root=/opt/g09gh/gaussian
source $g09root/g09/bsd/g09.profile
export GAUSS_SCRDIR=/scratch/your_username/
export GDMADIR=/opt/gdma/gdma-2.3.3/bin/
export PATH=/opt/gdma/gdma-2.3.3/bin/:$PATH
export PSI_SCRATCH=/scratch/your_username/ 
export PATH=/home/bdw2292/TheLatestTinker/bin/:$PATH
# you may add alias for submitting - replace path_to_poltype
#alias run_poltype="nohup python path_to_poltype/poltype.py &"
```


