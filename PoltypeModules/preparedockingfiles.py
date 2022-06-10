import sys
import getopt
import re
import os

'''
conda install -c bioconda mgltools --yes
'''






def receptorprep_usage():
    print("Usage: prepare_receptor4.py -r filename")
    print()
    print("    Description of command...")
    print("         -r   receptor_filename ")
    print("        supported file types include pdb,mol2,pdbq,pdbqs,pdbqt, possibly pqr,cif")
    print("    Optional parameters:")
    print("        [-v]  verbose output (default is minimal output)")
    print("        [-o pdbqt_filename]  (default is 'molecule_name.pdbqt')")
    print("        [-A]  type(s) of repairs to make: ")
    print("             'bonds_hydrogens': build bonds and add hydrogens ")
    print("             'bonds': build a single bond from each atom with no bonds to its closest neighbor") 
    print("             'hydrogens': add hydrogens")
    print("             'checkhydrogens': add hydrogens only if there are none already")
    print("             'None': do not make any repairs ")
    print("             (default is 'checkhydrogens')")
    print("        [-C]  preserve all input charges ie do not add new charges ")
    print("             (default is addition of gasteiger charges)")
    print("        [-p]  preserve input charges on specific atom types, eg -p Zn -p Fe")
    print("        [-U]  cleanup type:")
    print("             'nphs': merge charges and remove non-polar hydrogens")
    print("             'lps': merge charges and remove lone pairs")
    print("             'waters': remove water residues")
    print("             'nonstdres': remove chains composed entirely of residues of")
    print("                      types other than the standard 20 amino acids")
    print("             'deleteAltB': remove XX@B atoms and rename XX@A atoms->XX")
    print("             (default is 'nphs_lps_waters_nonstdres') ")
    print("        [-e]  delete every nonstd residue from any chain")
    print("              'True': any residue whose name is not in this list:")
    print("                      ['CYS','ILE','SER','VAL','GLN','LYS','ASN', ")
    print("                      'PRO','THR','PHE','ALA','HIS','GLY','ASP', ")
    print("                      'LEU', 'ARG', 'TRP', 'GLU', 'TYR','MET', ")
    print("                      'HID', 'HSP', 'HIE', 'HIP', 'CYX', 'CSS']")
    print("              will be deleted from any chain. ")
    print("              NB: there are no  nucleic acid residue names at all ")
    print("              in the list and no metals. ")
    print("             (default is False which means not to do this)")
    print("        [-M]  interactive ")
    print("             (default is 'automatic': outputfile is written with no further user input)")


def PrepareReceptorPDBQT(cmdstr):
    from MolKit import Read
    import MolKit.molecule
    import MolKit.protein
    from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation
    
    
    
    
    
    
    # process command arguments
    try:
        opt_list, args = getopt.getopt(cmdstr.split()[2:], 'r:vo:A:Cp:U:eM:')
    
    except (getopt.GetoptError, msg):
        print('prepare_receptor4.py: %s' %msg)
        receptorprep_usage()
        sys.exit(2)
    
    # initialize required parameters
    #-s: receptor
    receptor_filename =  None
    
    # optional parameters
    verbose = None
    #-A: repairs to make: add bonds and/or hydrogens or checkhydrogens
    repairs = ''
    #-C default: add gasteiger charges 
    charges_to_add = 'gasteiger'
    #-p preserve charges on specific atom types
    preserve_charge_types=None
    #-U: cleanup by merging nphs_lps, nphs, lps, waters, nonstdres
    cleanup  = "nphs_lps_waters_nonstdres"
    #-o outputfilename
    outputfilename = None
    #-m mode 
    mode = 'automatic'
    #-e delete every nonstd residue from each chain
    delete_single_nonstd_residues = None
    
    #'r:vo:A:Cp:U:eMh'
    for o, a in opt_list:
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print('set receptor_filename to ', a)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('set outputfilename to ', a)
        if o in ('-A', '--A'):
            repairs = a
            if verbose: print('set repairs to ', a)
        if o in ('-C', '--C'):
            charges_to_add = None
            if verbose: print('do not add charges')
        if o in ('-p', '--p'):
            if not preserve_charge_types:
                preserve_charge_types = a
            else:
                preserve_charge_types = preserve_charge_types + ','+ a
            if verbose: print('preserve initial charges on ', preserve_charge_types)
        if o in ('-U', '--U'):
            cleanup  = a
            if verbose: print('set cleanup to ', a)
        if o in ('-e', '--e'):
            delete_single_nonstd_residues  = True
            if verbose: print('set delete_single_nonstd_residues to True')
        if o in ('-M', '--M'):
            mode = a
            if verbose: print('set mode to ', a)
        if o in ('-h', '--'):
            receptorprep_usage()
            sys.exit()
    
    
    if not receptor_filename:
        print('prepare_receptor4: receptor filename must be specified.')
        receptorprep_usage()
        sys.exit()
    
    #what about nucleic acids???
    
    mols = Read(receptor_filename)
    if verbose: print('read ', receptor_filename)
    mol = mols[0]
    preserved = {}
    if charges_to_add is not None and preserve_charge_types is not None:
        preserved_types = preserve_charge_types.split(',') 
        if verbose: print("preserved_types=", preserved_types)
        for t in preserved_types:
            if verbose: print('preserving charges on type->', t)
            if not len(t): continue
            ats = mol.allAtoms.get(lambda x: x.autodock_element==t)
            if verbose: print("preserving charges on ", ats.name)
            for a in ats:
                if a.chargeSet is not None:
                    preserved[a] = [a.chargeSet, a.charge]
    
    if len(mols)>1:
        if verbose: print("more than one molecule in file")
        #use the molecule with the most atoms
        ctr = 1
        for m in mols[1:]:
            ctr += 1
            if len(m.allAtoms)>len(mol.allAtoms):
                mol = m
                if verbose: print("mol set to ", ctr, "th molecule with", len(mol.allAtoms), "atoms")
    mol.buildBondsByDistance()
    
    if verbose:
        print("setting up RPO with mode=", mode,"and outputfilename= ", outputfilename)
        print("charges_to_add=", charges_to_add)
        print("delete_single_nonstd_residues=", delete_single_nonstd_residues)
    
    RPO = AD4ReceptorPreparation(mol, mode, repairs, charges_to_add, 
                        cleanup, outputfilename=outputfilename,
                        preserved=preserved, 
                        delete_single_nonstd_residues=delete_single_nonstd_residues)    
    
    if charges_to_add is not None:
        #restore any previous charges
        for atom, chargeList in preserved.items():
            atom._charges[chargeList[0]] = chargeList[1]
            atom.chargeSet = chargeList[0]
    
    
    # To execute this command type:
    # prepare_receptor4.py -r pdb_file -o outputfilename -A checkhydrogens 




def gridprep_usage():
    print("Usage: prepare_gpf4.py -l pdbqt_file -r pdbqt_file ")
    print("     -l ligand_filename")
    print("     -r receptor_filename")
    print()
    print("Optional parameters:")
    print("    [-i reference_gpf_filename]")
    print("    [-o output_gpf_filename]")
    print("    [-x flexres_filename]")
    print("    [-p parameter=newvalue. For example: -p ligand_types='HD,Br,A,C,OA' ]")
    print("    [-d directory of ligands to use to set types]")
    print("    [-y boolean to center grids on center of ligand]")
    print("    [-n boolean to NOT size_box_to_include_ligand]")
    print("    [-v]")
    print()
    print("Prepare a grid parameter file (GPF) for AutoDock4.")
    print()
    print("   The GPF will by default be <receptor>.gpf. This")
    print("may be overridden using the -o flag.")


def ligandprep_usage():
    ("Print helpful, accurate usage statement to stdout.")
    print("Usage: prepare_ligand4.py -l filename")
    print()
    print("    Description of command...")
    print("         -l     ligand_filename (.pdb or .mol2 or .pdbq format)")
    print("    Optional parameters:")
    print("        [-v]    verbose output")
    print("        [-o pdbqt_filename] (default output filename is ligand_filename_stem + .pdbqt)")
    print("        [-d]    dictionary to write types list and number of active torsions ")
    
    print("        [-A]    type(s) of repairs to make:\n\t\t bonds_hydrogens, bonds, hydrogens (default is to do no repairs)")
    print("        [-C]    do not add charges (default is to add gasteiger charges)")
    print("        [-p]    preserve input charges on atom type, eg -p Zn")
    print("               (default is not to preserve charges on any specific atom type)")
    print("        [-U]    cleanup type:\n\t\t nphs_lps, nphs, lps, '' (default is 'nphs_lps') ")
    print("        [-B]    type(s) of bonds to allow to rotate ")
    print("               (default sets 'backbone' rotatable and 'amide' + 'guanidinium' non-rotatable)")
    print("        [-R]    index for root")
    print("        [-F]    check for and use largest non-bonded fragment (default is not to do this)")
    print("        [-M]    interactive (default is automatic output)")
    print("        [-I]    string of bonds to inactivate composed of ")
    print("                   of zero-based atom indices eg 5_13_2_10  ")
    print("                   will inactivate atoms[5]-atoms[13] bond ")
    print("                               and atoms[2]-atoms[10] bond ")
    print("                      (default is not to inactivate any specific bonds)")
    print("        [-Z]    inactivate all active torsions     ")
    print("                      (default is leave all rotatable active except amide and guanidinium)")
    print("        [-g]    attach all nonbonded fragments ")
    print("                      (default is not to do this)")

def prmfile_usage():
    print("Usage: prepare_dpf4.py -l pdbqt_file -r pdbqt_file")
    print("    -l ligand_filename")
    print("    -r receptor_filename")
    print()
    print("Optional parameters:")
    print("    [-o output dpf_filename]")
    print("    [-i template dpf_filename]")
    print("    [-x flexres_filename]")
    print("    [-p parameter_name=new_value]")
    print("    [-k list of parameters to write]")
    print("    [-e write epdb dpf ]")
    print("    [-v] verbose output")
    print("    [-L] use local search parameters")
    print("    [-S] use simulated annealing search parameters")
    print("    [-s] seed population using ligand's present conformation")
    print()
    print("Prepare a docking parameter file (DPF) for AutoDock4.")
    print()
    print("   The DPF will by default be <ligand>_<receptor>.dpf. This")
    print("may be overridden using the -o flag.")


def PrepareParameterFile(cmdstr):

    import string
    import os.path
    from MolKit import Read
    from AutoDockTools.DockingParameters import DockingParameters, DockingParameter4FileMaker, genetic_algorithm_list, genetic_algorithm_local_search_list4, local_search_list4,simulated_annealing_list4
                    
    
    
     
    
        
        
    
    try:
        opt_list, args = getopt.getopt(cmdstr.split()[2:], 'sLShvl:r:i:o:x:p:k:e')
    except (getopt.GetoptError, msg):
        print('prepare_dpf4.py: %s' % msg)
        prmfile_usage()
        sys.exit(2)
    
    receptor_filename = ligand_filename = None
    dpf_filename = None
    template_filename = None
    flexres_filename = None
    parameters = []
    parameter_list = genetic_algorithm_local_search_list4
    pop_seed = False
    verbose = None
    epdb_output = False
    for o, a in opt_list:
        if verbose: print("o=", o, ' a=', a)
        if o in ('-v', '--v'):
            verbose = 1
            if verbose: print('verbose output')
        if o in ('-l', '--l'):   #ligand filename
            ligand_filename = a
            if verbose: print('ligand_filename =', ligand_filename)
        if o in ('-r', '--r'):   #receptor filename
            receptor_filename = a
            if verbose: print('receptor_filename =', receptor_filename)
        if o in ('-x', '--x'):   #flexres_filename 
            flexres_filename = a
            if verbose: print('flexres_filename =', flexres_filename)
        if o in ('-i', '--i'):   #input reference
            template_filename = a
            if verbose: print('template_filename =', template_filename)
        if o in ('-o', '--o'):   #output filename
            dpf_filename = a
            if verbose: print('output dpf_filename =', dpf_filename)
        if o in ('-p', '--p'):   #parameter
            parameters.append(a)
            if verbose: print('parameters =', parameters)
        if o in ('-e', '--e'):
            epdb_output = True
            if verbose: print('output epdb file')
            parameter_list = epdb_list4_2
        if o in ('-k', '--k'):   #parameter_list_to_write
            parameter_list = a
            if verbose: print('parameter_list =', parameter_list)
        if o in ('-L', '--L'):   #parameter_list_to_write
            local_search = 1
            parameter_list = local_search_list4
            if verbose: print('parameter_list =', parameter_list)
        if o in ('-S', '--S'):   #parameter_list_to_write
            parameter_list = simulated_annealing_list4
            if verbose: print('parameter_list =', parameter_list)
        if o in ('-h', '--'):
            prmfile_usage()
            sys.exit()
        if o in ('-s'):
            pop_seed = True
    
    
    if (not receptor_filename) or (not ligand_filename):
        print("prepare_dpf4.py: ligand and receptor filenames")
        print("                    must be specified.")
        prmfile_usage()
        sys.exit()
    
    
    #9/2011: fixing local_search bugs:
    # specifically: 
    # 1. quaternion0 0 0 0 0  
    # 2. dihe0 0 0 0 0 0 <one per rotatable bond>
    # 3. about == tran0 
    # 4. remove tstep  qstep and dstep
    # 5. remove ls_search_freq
    local_search = parameter_list==local_search_list4
    dm = DockingParameter4FileMaker(verbose=verbose)
    if template_filename is not None:  #setup values by reading dpf
        dm.dpo.read(template_filename)
    dm.set_ligand(ligand_filename)
    dm.set_receptor(receptor_filename)
    if flexres_filename is not None:
        flexmol = Read(flexres_filename)[0]
        flexres_types = flexmol.allAtoms.autodock_element
        lig_types = dm.dpo['ligand_types']['value'].split()
        all_types = lig_types
        for t in flexres_types:
            if t not in all_types: 
                all_types.append(t)
        all_types_string = all_types[0]
        if len(all_types)>1:
            for t in all_types[1:]:
                all_types_string = all_types_string + " " + t
                if verbose: print("adding ", t, " to all_types->", all_types_string)
        dm.dpo['ligand_types']['value'] = all_types_string 
        dm.dpo['flexres']['value'] = flexres_filename
        dm.dpo['flexres_flag']['value'] = True
    #dm.set_docking_parameters( ga_num_evals=1750000,ga_pop_size=150, ga_run=20, rmstol=2.0)
    kw = {}    
    for p in parameters:
        key,newvalue = string.split(p, '=')
        #detect string reps of lists: eg "[1.,1.,1.]"
        if newvalue[0]=='[':
            nv = []
            for item in newvalue[1:-1].split(','):
                nv.append(float(item))
            #print "nv=", nv
            newvalue = nv
        if key=='epdb_flag':
            print("setting epdb_flag to", newvalue)
            kw['epdb_flag'] = 1
        elif key=='set_psw1':
            print("setting psw1_flag to", newvalue)
            kw['set_psw1'] = 1
            kw['set_sw1'] = 0
        elif key=='set_sw1':
            print("setting set_sw1 to", newvalue)
            kw['set_sw1'] = 1
            kw['set_psw1'] = 0
        elif key=='include_1_4_interactions_flag':
            kw['include_1_4_interactions'] = 1
        elif 'flag' in key:
            if newvalue in ['1','0']:
                newvalue = int(newvalue)
            if newvalue =='False':
                newvalue = False
            if newvalue =='True':
                newvalue = True
        elif local_search and 'about' in key:
            kw['about'] = newvalue
            kw['tran0'] = newvalue     
        else:         
            kw[key] = newvalue
        apply(dm.set_docking_parameters, (), kw)
        if key not in parameter_list:
            #special hack for output_pop_file
            if key=='output_pop_file':
                parameter_list.insert(parameter_list.index('set_ga'), key)
            else:
                parameter_list.append(key) 
    dm.write_dpf(dpf_filename, parameter_list, pop_seed)
   



def PrepareLigandPDBQT(cmdstr):
    import os 
    from MolKit import Read
    from AutoDockTools.MoleculePreparation import AD4LigandPreparation
    
    
    
    
    
    
    # process command arguments
    try:
        opt_list, args = getopt.getopt(cmdstr.split()[2:], 'l:vo:d:A:Cp:U:B:R:MFI:Zgh')
    except (getopt.GetoptError, msg):
        print('prepare_ligand4.py: %s' %msg)
        ligandprep_usage()
        sys.exit(2)
    
    # initialize required parameters
    #-l: ligand
    ligand_filename =  None
    # optional parameters
    verbose = None
    add_bonds = False
    #-A: repairs to make: add bonds and/or hydrogens
    repairs = ""
    #-C  default: add gasteiger charges 
    charges_to_add = 'gasteiger'
    #-p preserve charges on specific atom types
    preserve_charge_types=''
    #-U: cleanup by merging nphs_lps, nphs, lps
    cleanup  = "nphs_lps"
    #-B named rotatable bond type(s) to allow to rotate
    #allowed_bonds = ""
    allowed_bonds = "backbone"
    #-r  root
    root = 'auto'
    #-o outputfilename
    outputfilename = None
    #-F check_for_fragments
    check_for_fragments = False
    #-I bonds_to_inactivate
    bonds_to_inactivate = ""
    #-Z inactivate_all_torsions
    inactivate_all_torsions = False
    #-g attach_nonbonded_fragments
    attach_nonbonded_fragments = False
    #-m mode 
    mode = 'automatic'
    #-d dictionary
    dict = None
    
    #'l:vo:d:A:CKU:B:R:MFI:Zg'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-l', '--l'):
            ligand_filename = a
            if verbose: print('set ligand_filename to ', a)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print('set verbose to ', True)
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose: print('set outputfilename to ', a)
        if o in ('-d', '--d'):
            dict = a
            if verbose: print('set dict to ', a)
        if o in ('-A', '--A'):
            repairs = a
            if verbose: print('set repairs to ', a)
        if o in ('-C', '--C'):
            charges_to_add = None
            if verbose: print('do not add charges')
        if o in ('-p', '--p'):
            preserve_charge_types+=a
            preserve_charge_types+=','
            if verbose: print('preserve initial charges on ', preserve_charge_types)
        if o in ('-U', '--U'):
            cleanup  = a
            if verbose: print('set cleanup to merge ', a)
        if o in ('-B', '--B'):
            allowed_bonds = a
            if verbose: print('allow ', a, 'bonds set to rotate')
        if o in ('-R', '--R'):
            root = a
            if verbose: print('set root to ', root)
        if o in ('-F', '--F'):
            check_for_fragments = True
            if verbose: print('set check_for_fragments to True')
        if o in ('-M', '--M'):
            mode = a
            if verbose: print('set mode to ', a)
        if o in ('-I', '--I'):
            bonds_to_inactivate = a
            if verbose: print('set bonds_to_inactivate to ', a)
        if o in ('-Z', '--Z'):
            inactivate_all_torsions = True
            if verbose: print('set inactivate_all_torsions to ', inactivate_all_torsions)
        if o in ('-g', '--g'):
            attach_nonbonded_fragments = True
            if verbose: print('set attach_nonbonded_fragments to ', attach_nonbonded_fragments)
        if o in ('-h', '--'):
            ligandprep_usage()
            sys.exit()
    
    
    if not  ligand_filename:
        print('prepare_ligand4: ligand filename must be specified.')
        ligandprep_usage()
        sys.exit()
    
    mols = Read(ligand_filename)
    if verbose: print('read ', ligand_filename)
    mol = mols[0]
    if len(mols)>1:
        if verbose: 
            print("more than one molecule in file")
        #use the one molecule with the most atoms
        ctr = 1
        for m in mols[1:]:
            ctr += 1
            if len(m.allAtoms)>len(mol.allAtoms):
                mol = m
                if verbose:
                    print("mol set to ", ctr, "th molecule with", len(mol.allAtoms), "atoms")
    coord_dict = {}
    for a in mol.allAtoms: coord_dict[a] = a.coords
    
    
    mol.buildBondsByDistance()
    if charges_to_add is not None:
        preserved = {}
        preserved_types = preserve_charge_types.split(',') 
        for t in preserved_types:
            if not len(t): continue
            ats = mol.allAtoms.get(lambda x: x.autodock_element==t)
            for a in ats:
                if a.chargeSet is not None:
                    preserved[a] = [a.chargeSet, a.charge]
    
    
    
    if verbose:
        print("setting up LPO with mode=", mode,"and outputfilename= ", outputfilename)
        print("and check_for_fragments=", check_for_fragments)
        print("and bonds_to_inactivate=", bonds_to_inactivate)
    LPO = AD4LigandPreparation(mol, mode, repairs, charges_to_add, 
                            cleanup, allowed_bonds, root, 
                            outputfilename=outputfilename,
                            dict=dict, check_for_fragments=check_for_fragments,
                            bonds_to_inactivate=bonds_to_inactivate, 
                            inactivate_all_torsions=inactivate_all_torsions,
                            attach_nonbonded_fragments=attach_nonbonded_fragments)
    #do something about atoms with too many bonds (?)
    #FIX THIS: could be peptide ligand (???)
    #          ??use isPeptide to decide chargeSet??
    if charges_to_add is not None:
        #restore any previous charges
        for atom, chargeList in preserved.items():
            atom._charges[chargeList[0]] = chargeList[1]
            atom.chargeSet = chargeList[0]
    if verbose: print("returning ", mol.returnCode)
    bad_list = []
    for a in mol.allAtoms:
        if a.coords!=coord_dict[a]: bad_list.append(a)
    if len(bad_list):
        print(len(bad_list), ' atom coordinates changed!')
        for a in bad_list:
            print(a.name, ":", coord_dict[a], ' -> ', a.coords)
    else:
        if verbose: print("No change in atomic coordinates")
    if mol.returnCode!=0: 
        sys.stderr.write(mol.returnMsg+"\n")


def PrepareGridParameterFile(cmdstr):
    import string
    import os.path
    import glob
    from MolKit import Read
    from AutoDockTools.GridParameters import GridParameters, grid_parameter_list4
    from AutoDockTools.GridParameters import GridParameter4FileMaker
    from AutoDockTools.atomTypeTools import AutoDock4_AtomTyper
    
    
    
    
        
    
    try:
        opt_list, args = getopt.getopt(cmdstr.split()[2:], 'vl:r:i:x:o:p:d:yn')
    except (getopt.GetoptError, msg):
        print('prepare_gpf4.py: %s' % msg)
        gridprep_usage()
        sys.exit(2)
    
    receptor_filename = ligand_filename = None
    list_filename = gpf_filename = gpf_filename = None
    output_gpf_filename = None
    flexres_filename = None
    directory = None
    parameters = []
    verbose = None
    center_on_ligand = False
    size_box_to_include_ligand = True
    for o, a in opt_list:
        if o in ('-v', '--v'):
            verbose = 1
        if o in ('-l', '--l'):
            ligand_filename = a
            if verbose: print('ligand_filename=', ligand_filename)
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print('receptor_filename=', receptor_filename)
        if o in ('-i', '--i'):
            gpf_filename = a
            if verbose: print('reference_gpf_filename=', gpf_filename)
        if o in ('-x', '--x'):
            flexres_filename = a
            if verbose: print('flexres_filename=', flexres_filename)
        if o in ('-o', '--o'):
            output_gpf_filename = a
            if verbose: print('output_gpf_filename=', output_gpf_filename)
        if o in ('-p', '--p'):
            parameters.append(a)
            if verbose: print('parameters=', parameters)
        if o in ('-d', '--d'):
            directory = a
            if verbose: print('directory=', directory)
        if o in ('-y', '--y'):
            center_on_ligand = True
            if verbose: print('set center_on_ligand to ', center_on_ligand)
        if o in ('-n', '--n'):
            size_box_to_include_ligand = False
            if verbose: print('set size_box_to_include_ligand to ', size_box_to_include_ligand)
        if o in ('-h', '--'):
            gridprep_usage()
            sys.exit()
    
    
    if (not receptor_filename) or (ligand_filename is None and directory is None):
        print("prepare_gpf4.py: ligand and receptor filenames")
        print("                    must be specified.")
        gridprep_usage()
    gpfm = GridParameter4FileMaker(size_box_to_include_ligand=size_box_to_include_ligand,verbose=verbose)
    if gpf_filename is not None:
        gpfm.read_reference(gpf_filename)
    if ligand_filename is not None:
        gpfm.set_ligand(ligand_filename)
    gpfm.set_receptor(receptor_filename)
    if directory is not None:
        gpfm.set_types_from_directory(directory)
    if flexres_filename is not None:
        flexmol = Read(flexres_filename)[0]
        flexres_types = flexmol.allAtoms.autodock_element
        lig_types = gpfm.gpo['ligand_types']['value'].split()
        all_types = lig_types
        for t in flexres_types:
            if t not in all_types: 
                all_types.append(t)
        all_types_string = all_types[0]
        if len(all_types)>1:
            for t in all_types[1:]:
                all_types_string = all_types_string + " " + t
        gpfm.gpo['ligand_types']['value'] = all_types_string 
    for p in parameters:
        key,newvalue = string.split(p, '=')
        kw = {key:newvalue}
        apply(gpfm.set_grid_parameters, (), kw)
    #gpfm.set_grid_parameters(spacing=1.0)
    if center_on_ligand is True:
        gpfm.gpo['gridcenterAuto']['value'] = 0
        cenx,ceny,cenz = gpfm.ligand.getCenter()
        gpfm.gpo['gridcenter']['value'] = "%.3f %.3f %.3f" %(cenx,ceny,cenz)
    gpfm.write_gpf(output_gpf_filename)


def PrepareDockingFiles(ligandinputfilename,receptorinputfilename,dockgridcenter,dockgridsize,spacing):

    cmdstr='python prepare_receptor4.py -r '+receptorinputfilename
    receptorname=receptorinputfilename.replace('.pdb','.pdbqt')
    PrepareReceptorPDBQT(cmdstr)
    cmdstr='python prepare_ligand4.py '+'-l '+ligandinputfilename 
    ligandname=ligandinputfilename.replace('.pdb','.pdbqt')
    PrepareLigandPDBQT(cmdstr) 
    npts=[]
    for size in dockgridsize:
        points=int(size/spacing)
        npts.append(str(points))
    nptsstr=','.join(npts)
    dockgridcenter=[str(i) for i in dockgridcenter]
    dockgridcenterstr=','.join(dockgridcenter)
    cmdstr='python prepare_gpf4.py '+'-l '+ligandname+' -r '+receptorname+' -p '+'gridcenter='+dockgridcenterstr+' spacing='+str(spacing)+' -p '+'npts='+nptsstr
    PrepareGridParameterFile(cmdstr)
    cmdstr='python prepare_dpf4.py -l '+ligandname+' -r '+receptorname
    PrepareParameterFile(cmdstr)




opt_list, args = getopt.getopt(sys.argv[1:], 'dgc:dgz:s:ln:rn:',['dockgridcenter=','dockgridsize=','spacing=','ligandname=','receptorname='])
for o, a in opt_list:
    if o in ('-dgc','--dockgridcenter'):
        dockgridcenter=a.split(',')
        dockgridcenter=[float(i) for i in dockgridcenter]
    if o in ('--dockgridsize'):
        dockgridsize=a.split(',')
        dockgridsize=[float(i) for i in dockgridsize]
    if o in ('--spacing'):
        spacing=float(a)
    if o in ('--ligandname'):
        ligandfilename=a
    if o in ('--receptorname'):
        receptorfilename=a



PrepareDockingFiles(ligandfilename,receptorfilename,dockgridcenter,dockgridsize,spacing)
