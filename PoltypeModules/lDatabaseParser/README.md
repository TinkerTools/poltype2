# Automated Parameter Assignment for AMOEBA+ Model

## Introduction
This program (`lAssignAMOEBAplusPRM.py`) is designed to automatically assign the bonded and nonbonded parameters for the AMOEBA and AMOEBA+ models, based on chemical pattern matching. Here we will use phenol as an example (files can be found in example folder)

* !! If you don't have `pybel` and `numpy` in your python environment, please have them installed before you proceed.

### Bonded terms

It can automatically assign valence parameters to molecules based on current atomic types and parameter set. 

```shell
python lAssignAMOEBAplusPRM.py -potent bonded -xyz phenol.xyz -key phenol.key -konly NO
```
In addition, ranking tree and reversed searching algorithm can be used to assign the best bonded parameters if no exact match available.
* Format of the typing tree

  The whole structure of typing tree has been documented in file dat/typing_tree.log.
  The format for this file:
  
  `level   upper_index   general_index   type_index`
  
### Polarizability

Polarizability parameters have been derived for neutral molecules and common charged molecules. To assign this set of parameters, run the following command:

```shell
python lAssignAMOEBAplusPRM.py -potent polar -xyz phenol.xyz -key phenol.key
```

### Charge flux

* Note !!! Charge flux parameters have been re-optimized with bigger set of molecules. Details can be found in `CF_update_2022.docx` file.

Charge flux parameters can be assigned if `bond` and `angle` keywords exist in the `key` file. 

```shell
python lAssignAMOEBAplusPRM.py -potent CF -xyz phenol.xyz -key phenol.key
```

### Example

Go to the example folder and execute `sh run.sh`, you should be able to get the assigned parameter files including `phenol.key_boned`, `phenol.key_cf` and `phenol.key_polar`.

### Reference
* Non-bonded parameters (CP, CT, and VDW) (under development), 2022
* Atomic Polarizabilities for Interactive Dipole Induction Models. _J. Chem. Info. Model._ 2022, 62, 1, 79-87 [link](https://doi.org/10.1021/acs.jcim.1c01307)
* Valence Parameters for Organic molecules. _Journal of Computational Biophysics and Chemistry_ 0 0:0, 1-17 [link](https://doi.org/10.1142/S2737416521420047)
* CF parameters for organic molecules. _J. Chem. Phys._ 153, 064103 (2020); [link](https://doi.org/10.1063/5.0016376)
* CF model implementation & AMOEBA+ water. _J. Phys. Chem. Lett._ 2020, 11, 2, 419–426; [link](https://doi.org/10.1021/acs.jpclett.9b03489)
