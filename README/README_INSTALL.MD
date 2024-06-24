## Installation


### Required Programs
* Tinker, version >=8.9.4
* GDMA, version 2.3 and greater
* PSI4 (or Gaussian 09/16)


### Required Tinker CPU Installation
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


## Optional Tinker GPU Installation
* **Not required for parameterization purpose.**
* Install tinker CPU first. 
```
git clone git@github.com:TinkerTools/tinker9.git
cd tinker9
cp -r TinkerCPU/* tinker/.
mkdir compile
cd compile
export CUDAHOME=/usr/local/cuda-10.2
export CUDACXX=$CUDAHOME/bin/nvcc
export FC=/usr/bin/gfortran
export CXX=/usr/bin/g++
export ACC=/home/liuchw/shared/nvidia/hpc_sdk/Linux_x86_64/21.1/compilers/bin/nvc++
export opt=release
export host=0
export prec=m
export compute_capability=80
export cuda_dir=$CUDAHOME
export CMAKEHOME=/home/liuchw/shared/cmake3.12/bin/
$CMAKEHOME/cmake ../
make -j
```


### Python Envioronment
* Simply using the following git command to download Poltype
```shell
git clone https://github.com/TinkerTools/poltype2.git master
```

* Install anaconda
* You can use the provided `environment.yml` file (under Environments folder) to create a conda environment named "poltype" with `conda env create -f environment.yml`

* Download GDMA-2.3 and greater
* Note do not put source in your default .bashrc profile
* Note do not keep multiple sources for different purposes (unless necessary) in the same .bashrc source file, this can cause internal conflicts to arise

### Environment var for running POLTYPE:
 Create environments for Tinker,Gaussian 09,and GDMA. Put these in a file and source it before running POLTYPE so that programs such as "analyze.x" "gdma" "g16" can be found. You can use "which" to verify.

* Example bashrc
```shell
OSVERSION=`cat /etc/system-release`
place="${OSVERSION//[!0-9]/}"
OSVERSION=${place::1}
VAL=`nvidia-smi &> /dev/null; echo $?`

if [ $VAL != 0 ]; then
  echo -e "\e[101mCUDA utility not installed on `hostname`\e[0m"
  declare -i VALUE=0
else

  for dir in /usr/local/*;
  do
      SUB='cuda'
      if [[ "$dir" == *"$SUB"* ]]; then
        VALUE=`echo $dir | tr -dc '0-9'`
      fi
       
  done

  
  if [ $VALUE == 112 ] ; then

       export PATH=/home/liuchw/Softwares/tinkers/Tinker9/2206/build_cuda11:$PATH
       export CUDA_DEVICE_ORDER=PCI_BUS_ID
       export CUDA_VISIBLE_DEVICES=0

  fi
fi

if [ $OSVERSION == 8 ] ; then
    export PATH=/home/bdw2292/NewestTinkerOS8/bin/:$PATH
     
else
    export PATH=/home/bdw2292/NewestTinker/bin/:$PATH

fi

export PATH=/opt/CCDC/CCDC2020/Discovery_2020/bin/:$PATH # for GOLD Docking
export myusername=`whoami`
conda activate amoebamdpoltype
export g09root=/opt/g09gh/gaussian
source $g09root/g09/bsd/g09.profile
export GAUSS_SCRDIR=/scratch/$myusername/
export GDMADIR=/opt/gdma/gdma-2.3.3/bin/
export PATH=/opt/gdma/gdma-2.3.3/bin/:$PATH
export PSI_SCRATCH=/scratch/$myusername/
```


### Explicit Conda Installation
* Only do this if .yml file doesnt work
```
conda create --name amoebamdpoltype -c psi4 psi4 --yes
conda install -c conda-forge geometric --yes
conda activate amoebamdpoltype
conda config --add channels conda-forge
conda install mdanalysis --yes
conda install -c conda-forge openbabel --yes
conda install -c rdkit rdkit --yes
pip install PyAstronomy
conda install -c conda-forge svgutils --yes
conda install -c conda-forge cairosvg --yes
conda install -c conda-forge psutil --yes
conda install -c conda-forge mdtraj --yes
conda install -c conda-forge pymbar --yes
conda install -c salilab modeller --yes
conda install -c anaconda openpyxl --yes
conda install -c conda-forge pdb2pqr --yes
conda install -c conda-forge pdbfixer --yes
```

### Explicit Conda Installation - XTB
```
conda create --name xtbenv
conda activate xtbenv
conda install -c conda-forge xtb --yes
```
* add the path to xtb binary in PATH


### Explicit Conda Installation - ANI
```
conda install -c conda-forge ase --yes
conda install -c conda-forge torchani --yes
conda install -c conda-forge openbabel --yes
```
* Install via ani.yml file
