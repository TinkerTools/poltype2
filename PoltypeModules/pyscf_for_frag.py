# Import both PySCF setup class
from pyscf_setup import PySCF_init_setup
from pyscf_setup import PySCF_post_run
import os

def setup_frag_opt(poltype,newtorset,toranglelist,\
                    pre,inputmol,inputname,inputstruc):
    
    """
    This function setup the PySCF file to do
    a geometry optimization.

    Inputs:
        -   poltype: poltype class object (class obj)
        -   newtorset: list of torsion to restrain (list of list)
        -   toranglelist: list of torsion angle to restrain (list of float)
        -   pre: prefix name (string)
        -   inputmol: molecule object (babel or rdkit mol obj)
        -   inputname: name of the input file (string)
        -   inputstruc: name if the file with input geometry (string)

    Outputs:
        -   Opt_prep.PySCF_inp_file: Name of PySCF input file
        -   Opt_prep.PySCF_out_file: Name of PySCF output file
        -   cmd: command to run        

    """

    # Define the charge of the molecule 
    charge = inputmol.GetTotalCharge()

    # Instantiate the PySCF initial setup object
    Opt_prep = PySCF_init_setup(poltype,os.getcwd(),pre,inputmol,\
                               False,charge)

    # Define the input geometry file
    Opt_prep.PySCF_inp_xyz = inputstruc

    # Write the geometric constrain file
    Opt_prep.write_geometric_tor_const(pre,newtorset)

    # Write the PySCF input file: is_opt = True, is_frag = True
    Opt_prep.write_PySCF_input(True,True)

    # Define the command to run
    cmd = f'python {Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_inp_file} &> {Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}'

    return Opt_prep.PySCF_inp_file,Opt_prep.PySCF_out_file,cmd


def read_frag_opt(poltype,outputlog,mol):
   
    """
    This function creates the post run PySCf object and
    reads the optimized geometry
    
    Inputs:
        -   poltype: poltype class object (class obj)
        -   outputlog: Name of the output file (string)
        -   mol: molecule object (babel or rdkit mol obj)

    """

    # Instantiate the PySCF post run class 
    Opt_post = PySCF_post_run(poltype,os.getcwd(),outputlog,mol)

    # Read and create the optimized geometry file
    Opt_post.read_out()


def setup_frag_sp(poltype,inputmol,inputstruc,pre):

    """
    This function set up the fragment SP calculation

    Inputs:
        -   poltype: poltype class object (class obj)
        -   inputmol: molecule object (babel or rdkit mol obj)
        -   inputstruc: name if the file with input geometry (string)
        -   pre: prefix name (string)

    Outputs:
        -   Opt_prep.PySCF_inp_file: Name of PySCF input file
        -   Opt_prep.PySCF_out_file: Name of PySCF output file
        -   cmd: command to run        

    """

    # Define the charge of the molecule
    charge = inputmol.GetTotalCharge()

    # Instantiate the PySCF initial setup object 
    Opt_prep = PySCF_init_setup(poltype,os.getcwd(),pre,inputmol,\
                               False,charge)

    # Define the input geometry file
    Opt_prep.PySCF_inp_xyz = inputstruc

    # Write the PySCF input file: is_opt = False, is_frag = True
    Opt_prep.write_PySCF_input(False,True)

    # Define the command to run
    cmd = f'python {Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_inp_file} &> {Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}'

    return Opt_prep.PySCF_inp_file,Opt_prep.PySCF_out_file,cmd
 
