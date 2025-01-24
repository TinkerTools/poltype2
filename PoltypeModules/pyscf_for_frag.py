# Import both PySCF setup class
from pyscf_setup import PySCF_init_setup
from pyscf_setup import PySCF_post_run
import os

def setup_frag_opt(poltype,newtorset,toranglelist,\
                    pre,inputmol,inputname,inputstruc):

    
    charge = inputmol.GetTotalCharge()
    print('In setup_frag_opt')
    print(newtorset)
    print(toranglelist)
    print(pre)
    print(inputname)
    print(inputstruc)
    print(charge)

    Opt_prep = PySCF_init_setup(poltype,os.getcwd(),inputstruc,inputmol,\
                               False,charge)

    #Opt_prep.read_Tinker_inp_coord()

    #Opt_prep.write_xyz('init') 

    Opt_prep.PySCF_inp_xyz = Opt_prep.init_data['COM_file']

    Opt_prep.write_geometric_tor_const(pre,newtorset)

    Opt_prep.write_PySCF_input(True,True)
