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

    Opt_prep = PySCF_init_setup(poltype,os.getcwd(),pre,inputmol,\
                               False,charge)

    #Opt_prep.read_Tinker_inp_coord()

    #Opt_prep.write_xyz('init') 

    Opt_prep.PySCF_inp_xyz = inputstruc

    Opt_prep.write_geometric_tor_const(pre,newtorset)

    Opt_prep.write_PySCF_input(True,True)

    cmd = f'python {Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_inp_file} &> {Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}'

    return Opt_prep.PySCF_inp_file,Opt_prep.PySCF_out_file,cmd

    #cmd_to_run = {cmd: f'{Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}'}

    #if os.path.isfile(f'{Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}'):
    #    os.remove(f'{Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}')


    #if poltype.externalapi==None:
    #    #Temporary fix to work on CheckNormalTermination
    #    #if os.path.isfile(f'{Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}'):
    #    #    pass
    #    #else:
    #    finishedjobs,errorjobs=poltype.CallJobsSeriallyLocalHost(cmd_to_run,False)
    #    print(f'Optimization for Frag {Opt_prep.PySCF_inp_file}')
    #    print(finishedjobs,errorjobs)


def read_frag_opt(poltype,outputlog,mol):
    
    Opt_post = PySCF_post_run(poltype,os.getcwd(),outputlog,mol)

    Opt_post.read_out()


def setup_frag_sp(poltype,inputmol,inputstruc,pre):

    charge = inputmol.GetTotalCharge()
    print('In setup_frag_sp')
    print(pre)
    print(inputstruc)
    print(charge)

    Opt_prep = PySCF_init_setup(poltype,os.getcwd(),pre,inputmol,\
                               False,charge)


    Opt_prep.PySCF_inp_xyz = inputstruc
    Opt_prep.write_PySCF_input(False,True)


    cmd = f'python {Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_inp_file} &> {Opt_prep.init_data["topdir"]}/{Opt_prep.PySCF_out_file}'

    return Opt_prep.PySCF_inp_file,Opt_prep.PySCF_out_file,cmd
 
