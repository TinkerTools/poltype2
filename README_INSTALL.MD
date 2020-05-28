## Installation

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
```
make -j 8
make install


* Download GDMA-2.2

* Download POLTYPE
Make sure to download the ZIP folder from the poltype2 branch, NOT the master branch

```https://github.com/pren/poltype/tree/poltype2 ```

* Or simply using the following git command
```shell
git clone git@github.com:pren/poltype.git --branch poltype2
```

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

After you have conda in your PATH, run the following commands. These will be installed under "/home/username/miniconda3/". You can use the provided `environment.yml` file to create a conda environment named "poltype" with `conda env create -f environment.yml`, or run the following commands:

```shell
conda create -n poltype -c psi4/label/dev psi4 python=3.6
conda activate poltype
conda install scipy
conda install matplotlib
conda install -c conda-forge openbabel=2.4.1
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

* Now clone the poltype2 branch
```
git clone -b poltype2 git@github.com:pren/poltype.git poltype-poltype2
```

* If you dont have a GitHub account, make one
* Click the watch button on https://github.com/pren/poltype for important updates
