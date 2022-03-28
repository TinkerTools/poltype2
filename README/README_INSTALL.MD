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


## Tinker GPU Installation
* Install tinker CPU first 
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

* Install anaconda https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
* You may need to run "conda init bash" and then source ~/.bashrc first.
* After you have conda in your PATH, run the following commands. These will be installed under "/home/username/miniconda3/". You can use the provided `environment.yml` file to create a conda environment named "poltype" with `conda env create -f environment.yml`, or run the following commands:

```shell
conda create -n amoebamdpoltype python=3.8 --yes
conda activate amoebamdpoltype
conda install scipy --yes
conda install matplotlib --yes
conda install -c conda-forge openbabel --yes
conda install -c conda-forge rdkit --yes
conda install -c conda-forge mdanalysis --yes
conda install -c conda-forge svgutils --yes
conda install -c conda-forge cairosvg --yes
conda install -c anaconda psutil --yes
conda install -c anaconda numpy --yes
conda install git pip --yes
pip install PyAstronomy
conda install -c conda-forge -c schrodinger pymol-bundle --yes
conda install -c anaconda scikit-learn --yes
conda install -c conda-forge mdtraj --yes
conda upgrade numpy --yes
conda install -c conda-forge forcebalance --yes
conda install -c conda-forge pymbar --yes
conda install -c anaconda packaging --yes
conda install -c conda-forge tqdm --yes
conda install -c psi4 psi4 --yes
```

### Environment var for running POLTYPE:
 Create environments for Tinker,Gaussian 09,and GDMA. Put these in a file and source it before running POLTYPE

* Download GDMA-2.3 and greater
* Note do not put source in your default .bashrc profile
* Note do not keep multiple sources for different purposes (unless necessary) in the same .bashrc source file, this can cause internal conflicts to arise

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
  

  if [ $VALUE == 102 ] ; then
       export TINKERGPUDIR=/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda10.2/
       function dynamic_gpu () {
           local TINKER9=/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda10.2
       
           $TINKER9/dynamic9.sh $@ 
       }
       export -f dynamic_gpu
       
       function bar_gpu () {
           local TINKER9=/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda10.2
       
           $TINKER9/bar9.sh $@  
       }
       export -f bar_gpu
       
       function analyze_gpu () {
           local TINKER9=/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda10.2
       
           $TINKER9/analyze9.sh $@
       }
       export -f analyze_gpu
       export GPUDYNAMICS=True # need this for ForceBalance to let program know to use dynamic_gpu
  elif [ $VALUE == 112 ] ; then
       export TINKERGPUDIR=/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda11.2
       function dynamic_gpu () {
           local TINKER9=/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda11.2
       
           $TINKER9/dynamic9.sh $@
       }
       export -f dynamic_gpu
       
       function bar_gpu () {
           local TINKER9=/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda11.2
       
           $TINKER9/bar9.sh $@
       }
       export -f bar_gpu
       
       function analyze_gpu () {
           local TINKER9=/home/liuchw/Softwares/tinkers/Tinker9-latest/build_cuda11.2
       
           $TINKER9/analyze9.sh $@
       }
       export -f analyze_gpu
       
       export GPUDYNAMICS=True
  fi

 


fi
if [ $OSVERSION == 8 ] ; then
    export PATH=/home/bdw2292/NewestTinkerOS8/bin/:$PATH
    export TINKERPATH=/home/bdw2292/NewestTinkerOS8/bin/ # this is for ForceBalance to recognize tinker path in different nodes
     
else
    export PATH=/home/bdw2292/NewestTinker/bin/:$PATH
    export TINKERPATH=/home/bdw2292/NewestTinker/bin/

fi

export myusername=`whoami`
conda activate amoebamdpoltype
export g09root=/opt/g09gh/gaussian
source $g09root/g09/bsd/g09.profile
export GAUSS_SCRDIR=/scratch/$myusername/
export GDMADIR=/opt/gdma/gdma-2.3.3/bin/
export PATH=/opt/gdma/gdma-2.3.3/bin/:$PATH
export PSI_SCRATCH=/scratch/$myusername/
```


### ForceBalance Wrapper Setup
* Install modified forcebalance for Tinker9 GPU dynamics
```
cd ~
git clone https://github.com/bdw2292/forcebalance.git
```

* Make sure to copy files from modified forcebalance to conda site packages (for example mine are located /home/bdw2292/miniconda3/envs/FBTest/lib/python3.8/site-packages/forcebalance/) 
```
conda install -c conda-forge forcebalance --yes
cp ~/forcebalance/src/*py /home/bdw2292/miniconda3/envs/amoebamdpoltype/lib/python3.6/site-packages/forcebalance/
cp ~/forcebalance/src/data/*py /home/bdw2292/miniconda3/envs/amoebamdpoltype/lib/python3.6/site-packages/forcebalance/data/

```


### ForceBalance Results Viewing

* Make sure pymol license file is located in /home/USERNAME/.pymol/license.lic

