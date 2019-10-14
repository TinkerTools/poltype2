# Poltype README
--------------

Poltype is an automatic parameterization tool for the AMOEBA polarizable force field.


https://biomol.bme.utexas.edu/tinkergpu/index.php/Poltype:Poltype

## INSTALLATION:

1)

https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Download the installer:

    Miniconda installer for Linux.

    Anaconda installer for Linux.

In your terminal window, run:

    Miniconda:

    ```bash Miniconda3-latest-Linux-x86_64.sh```

    Anaconda:

    ```bash Anaconda-latest-Linux-x86_64.sh```

2)
```
conda install -c psi4 psi4 
conda create -n p4env psi4 -c psi4 # creates enviormental enviorments for psi4
conda install numpy
conda install scipy
conda install matplotlib
conda install -c openbabel openbabel

conda activate p4env ``` 

3) 
Create enviorments for Tinker,Gaussian 09,and GDMA

```bash
export PATH=path_to_tinker_binary:$PATH
export g09root=path_to_gaussian_binary
export GAUSS_SCRDIR=path_to_scratch_directory
export GDMADIR=path_to_gdma_directory
export PATH=path_to_gdma_directory:$PATH
export PSI_SCRATCH=path_to_scratch_directory
```



=====================================================================
Poltype is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 (gpl-3.0.txt)
as published by the Free Software Foundation.

Poltype is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Poltype; if not, write to:
Free Software Foundation, Inc.
59 Temple Place, Suite 330
Boston, MA 02111-1307  USA

