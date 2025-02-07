import subprocess
import os

class PySCF_init_setup():

    """
    This class corresponds to the initial setup
    of a PySCF calculation.

    It reads a geometry and prepare the pyscf input file,
    torsion restraints and final output.
    
    !!!! WARNING !!!
    You need to make sure that pyscf version is 2.7.0 to use wb97x-d3 !!

    Inputs:
        -   poltype_obj: poltype class object (class obj)
        -   cure_dir: path for the current directory (string)
        -   COM_file: initial file name (string)
        -   mol: molecule object (babel or rdkit mol obj)
        -   skip_error: Skip or not the errors when checking for failures (bool)
        -   charge: charge of the molecule (None or int)
        

    """

    def __init__(self,poltype_obj,cure_dir,COM_file,mol,\
                 skip_error,charge):

        self.pol_obj = poltype_obj
        self.mol_obj = mol
        self.init_data = {'topdir': cure_dir, 'COM_file': COM_file,
                          'skip_error': skip_error,
                          }
        self.mol_data_init = {'charge': charge}


    def read_Gauss_inp_coord(self):

        """
        This function reads the intitial coordinate from a gaussian 
        input file. 
        This is one is used only when running the optimization of the
        parent ligand.
        The coordinate and atom list are saved in two lists.
    
        """
    
        # Initialize the atom list and coordinate
        self.mol_data_init['Atoms_list'] = []
        self.mol_data_init['Atoms_coord'] = []

        # Loop through the gaussian file and grab the coordinates
        with open(f"{self.init_data['topdir']}/{self.init_data['COM_file']}", 'r') as f:
            for lines in f:
                L = lines.strip('\n').split()
                if len(L)==4 and '#' not in L:
                    self.mol_data_init['Atoms_list'].append(L[0])
                    self.mol_data_init['Atoms_coord'].append([float(i) for i in L[1:]])

        return


    def read_Tinker_inp_coord(self):

        """
        This function reads the coordinate from a Tinker xyz file.
        This function is used when the torsion scan is being done.
        The coordinate and atom list are saved in two lists.

        """

        # Initialize the atom list and coordinate
        self.mol_data_init['Atoms_list'] = []
        self.mol_data_init['Atoms_coord'] = []

        # Loop through the tinker xyz file and grab the coordinates
        with open(f"{self.init_data['topdir']}/{self.init_data['COM_file']}", 'r') as f:
            next(f)
            next(f)
            for lines in f:
                L = lines.strip('\n').split()
                self.mol_data_init['Atoms_list'].append(L[0])
                self.mol_data_init['Atoms_coord'].append([float(i) for i in L[1:]])

        return


    def write_xyz(self,name):

        """
        This function writes a standard xyz file to be used as 
        input geometry for PySCF.
        The coordinates are taken from either a Gaussian .com 
        or tinker xyz file.

        Inputs:
            -   name: Name to give to the xyz file to be written (string)

        """

        # Define the PySCF input xyz file name
        self.PySCF_inp_xyz = f"{self.init_data['COM_file'].strip('.com')}_{name}.xyz"

        # Write the xyz file
        with open(f"{self.init_data['topdir']}/{self.PySCF_inp_xyz}", 'w') as f:
            f.write(f"{len(self.mol_data_init['Atoms_coord'])}\n\n")
            for A,C in zip(self.mol_data_init['Atoms_list'],self.mol_data_init['Atoms_coord']):
                f.write(f'{A:3s} {C[0]:9.6f} {C[1]:9.6f} {C[2]:9.6f}\n')

        return


    def write_geometric_tor_const(self,name,tor_const, tor_angle=None):

        """
        This function writes the torsion constraint file in the 
        GeomeTRIC file format. This file is being used by PySCF.

        Inputs:
            -   name: Name to give to the constraint file to be written (string)
            -   tor_const: the torsion atom index to be constrained (list of list)
            -   tor_angle: torsion angle to constrained (list of float)

        """

        # Define the PySCF constraint file name
        self.PySCF_inp_const = f"PySCF_{name}_constraint.txt"

        # Write the constraint file using the given torsion index
        # and tor angle, or by calculting it
        with open(f"{self.init_data['topdir']}/{self.PySCF_inp_const}", 'w') as f:
            f.write('$set\n')
            for i,tor in enumerate(tor_const):
                if tor_angle == None:
                    angle = self.mol_obj.GetTorsion(tor[0], tor[1], tor[2], tor[3]) % 360
                else:
                    angle = tor_angle[i]
                f.write(f'dihedral {tor[0]} {tor[1]} {tor[2]} {tor[3]} {angle:3.4f}\n')

        return


    def write_PySCF_input(self,is_opt=True,is_frag=False):

        """
        This function writes the PySCF input file.
        The input depends on: 
            - Is the molecule the parent ligand or a fragment
            - Is it an optimization or single point

        Inputs:
            -   is_opt: is it an optimization (bool)
            -   is_frag: is the molecule a fragment (bool)

        """

        # If the molecule us a fragment
        if is_frag:

            # Need to add pyscf in the input file name to make sure
            # that it is properly submitted per poltype.call_subsystem
            self.PySCF_inp_file = f"{self.init_data['COM_file']}_pyscf.py"

            # Remove pyscf in the log name for SP calculation 
            if is_opt:
                self.PySCF_out_file = f"{self.init_data['COM_file']}_pyscf.log"
            else:
                self.PySCF_out_file = f"{self.init_data['COM_file']}.log"
        else:
            self.PySCF_inp_file = f"{self.init_data['COM_file'].strip('.com')}_pyscf.py"
            self.PySCF_out_file = f"{self.init_data['COM_file'].strip('.com')}_pyscf.log" 

        # Write the PySCF header
        self.write_PySCF_header(is_opt)

        # If it is a geometry optimization
        # Write the PySCF params part
        if is_opt:
            self.write_PySCF_params()

        # Write the main core of the PySCF script
        self.write_PySCF_main(is_opt,is_frag)

        return


    def write_PySCF_header(self,is_opt):

        """
        This function writes the header in the PySCF input file.
        It only writes the library import

        Inputs:
            -   is_opt: is it an optimization (bool)
            
        """

        # Write the header part with all the import needed
        # If it not an optimization, no need to import optimize
        with open(f"{self.init_data['topdir']}/{self.PySCF_inp_file}", 'w') as f:
            f.write('import pyscf\n')
            f.write('from pyscf import gto, scf\n')
            f.write('from pyscf.dft import rks,uks\n')
            f.write('from pyscf import dft\n')
            if is_opt:
                f.write('from pyscf.geomopt.geometric_solver import optimize\n\n')

        return


    def write_PySCF_params(self):

        """
        This function writes the PySCF parameters needed
        to pass to GeomeTRIC. 
        It gives the contraints file and the convergence criteria

        """

        # Write the params to pass to GeomeTRIC
        with open(f"{self.init_data['topdir']}/{self.PySCF_inp_file}", 'a+') as f:
            f.write('params = {\n')
            f.write(f'"constraints": "{self.PySCF_inp_const}",\n')
            if self.pol_obj.optconvergence == 'LOOSE':
                f.write('            "convergence_energy": 1e-4,\n') 
                f.write('            "convergence_grms": 1.7e-3,\n') 
                f.write('            "convergence_gmax": 2.5e-3,\n') 
                f.write('            "convergence_drms": 6.7e-3,\n') 
                f.write('            "convergence_dmax": 1.0e-2,\n') 
            f.write('}\n\n')
 
        return 


    def write_PySCF_main(self,is_opt,is_frag):

        """
        This function writes the main part of the PySCF script.

        Inputs:
            -   is_opt: is it an optimization (bool)
            -   is_frag: is the molecule a fragment (bool) 

        """

        # Define the atom string == input geometry
        atm_s = f'atom="{self.PySCF_inp_xyz}"'

        # Define the basis string == basis used for QM
        # It is using the same basis as for Psi in gas phase
        if is_opt:
            if not is_frag:
                basis_s = f'basis = "{self.pol_obj.optbasisset}"'
            else:
                basis_s = f'basis = "{self.pol_obj.torspbasisset}"'
        else:
            if is_frag:
                basis_s = f'basis = "{self.pol_obj.torspbasisset}"'
            else:
                self.pol_obj.WriteToLog('PySCF is not set up to run PCM single point QM for parent ligand\n Talk to a developer\n\n')
                raise Exception('PySCF is not set up to run PCM single point QM for parent ligand\n Talk to a developer\n\n')
                
        # Define the spin string == multiplicity of the molecule (2S in PySCF)
        spin_s = f'spin = {int(self.mol_obj.GetTotalSpinMultiplicity())-1}'

        # Define the charge string == charge of the molecule 
        # Either given as input or determine with mol object
        if self.mol_data_init['charge'] == None:
            chg_s = f'charge = {self.mol_obj.GetTotalCharge()}'
        else:
            chg_s = f'charge = {self.mol_data_init["charge"]}'

        # Define the unit string for the geometry (always A so far)
        unit_s = 'unit = "A"'

        # Write the main part of the input file
        with open(f"{self.init_data['topdir']}/{self.PySCF_inp_file}", 'a+') as f:

            # Define the molecule and build it
            f.write(f'mol = pyscf.M({atm_s}, {basis_s}, {spin_s}, {chg_s}, {unit_s})\n')
            f.write('mol.build()\n')

            # For now use the same DFT functional for both fragment and parent
            # ligand
            if not is_frag:
                f.write(f'mf = pyscf.solvent.PCM(mol.UKS(xc = "{self.pol_obj.pyscf_opt_meth}").density_fit())\n')
            else:
                f.write(f'mf = pyscf.solvent.PCM(mol.UKS(xc = "{self.pol_obj.pyscf_opt_meth}").density_fit())\n')
            
            # Define the solvent part 
            f.write(f'mf.with_solvent.method = "{self.pol_obj.pyscf_sol_imp}"\n')
            f.write(f'mf.with_solvent.eps = {self.pol_obj.pyscf_sol_eps}\n\n')

            # Define the optimization part
            if is_opt:
                f.write(f'mol_eq = optimize(mf, maxsteps={self.pol_obj.optmaxcycle}, **params)\n')
                f.write('print(f"Final optimized structure:")\n')
                f.write('print(mol_eq.tostring())\n\n')
            
            # Define the single point part
            # "Normal Termination" in added in PySCF output to ensure that
            # poltype.CheckNormalTermination can correctly process good termination 
            else:
                f.write('mol_sp = mf.kernel()\n')
                f.write('print(f"Final Energy: {mol_sp}")\n')
                f.write('print("Normal Termination")\n')
        return

    
class PySCF_post_run():

    """
    This class corresponds to the PySCF post run analysis.
    For now, it only grabs the optimized geometry.
    The SP energy grab is done in torsionfit.py file

    Inputs:
        -   poltype_obj: poltype class object (class obj)
        -   cure_dir: path for the current directory (string)
        -   output_file: name of the pyscf outptu file (string)
        -   mol_obj: molecule object (babel or rdkit mol obj)

    """


    def __init__(self,poltype_obj,cure_dir,output_file,mol_obj):

        self.pol_obj = poltype_obj
        self.topdir = cure_dir
        self.out_file = output_file
        self.mol_obj = mol_obj


    def read_out(self):

        """
        This function reads the optimized structure in the 
        PySCF output file and write in a standard xyz format

        """

        # Use grep to get the optimized geometry
        val = subprocess.check_output(f"grep -A {self.mol_obj.NumAtoms()} 'Final optimized structure' {self.topdir}/{self.out_file}",shell=True)
        Geom = val.decode("utf-8").split()

        # Define the final optimized xyz file 
        self.final_opt_xyz = f"{self.out_file.split('.log')[0]}.xyz"

        # Write the optimized coordinate
        start = 3
        with open(f'{self.topdir}/{self.final_opt_xyz}', 'w') as file_1:
            file_1.write(f'{self.mol_obj.NumAtoms()}\n\n')
            for i in range(0,len(Geom[3:]),4):
                end = start + 4
                file_1.write(f'{" ".join([j for j in Geom[start:end]])}\n')
                start = end
        
        return 
