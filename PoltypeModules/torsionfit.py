import torsiongenerator as torgen
import symmetry as symm
import optimization as opt
import electrostaticpotential as esp
import fragmenter as frag
import os
import sys
import shutil
import re
import numpy
from scipy import optimize
from scipy.optimize import fmin
import matplotlib
import pylab as plt
import time
from itertools import product,combinations



def insert_torprmdict_angle(poltype,angle, angledict):
    """
    Intent: Increase the count of this angle by one
    """
    anglekey = int(angle)
    if anglekey in angledict:
        angledict[anglekey] += 1
    else:
        angledict[anglekey] = 1

def tor_func_term (poltype,parms, x, nfold, C, angle, offset):
    """
    Intent: Returns energy vs dihedral angle profile of the torsion
    Input: 
        parms: torsion parameter estimate
        x: angle list in radians 
        nfold: fold #
        C: number of times this torsion exists
        angle: current dihedral angle
        offset: an offset for the current fold if it exists
    Output:
        energy profile
    Referenced By: fitfunc
    Description: -
    """
    return C*(parms/2.0)*(1+numpy.cos(nfold*(x+angle)+offset))


def fitfunc (poltype,parms, x,torset, torprmdict, debug = False):
    """
    Intent: Gives energy due to torsion around one rotatable bond 
    (torsion energy vs. dihedral angle) given a set of parameters
    This function is used to make the callable function 'errfunc' used by leastsq
    Input:
        parms: current parameter estimate
        x: angle list in radians 
        torprmdict: contains information about the torsions (like nfolds)
    Output:
        tor_energy: energy due to torsion at various dihedral angles (found using 'parms')
                    (for one rotatable bond)
    Referenced By: fit_rot_bond_tors 
    Description: Loops over torsions about this rotatable bond, then further loops about nfolds,
    making multiple calls to 'tor_func_term', summing its results in 'tor_energy'
    """
    tor_energy_array = [ 0.0 ] * len(x)
    offset = 0
    for j in range(len(x)):
        angtup=x[j]
        tor_energy=0
        for i in range(len(angtup)):
            ang=angtup[i]
            tor=torset[i]
            a,b,c,d=tor[0:4]
            clskey = torgen.get_class_key(poltype,a,b,c,d)
            torprm=torprmdict[clskey]
            # for each torsion about the same rotatable bond (clskey)
            for nfold in torprm['prmdict']:
                # for each nfold for this torsion (# of parameters)
                for clsangle, clscnt in torprm['phasedict'].items():
                    # current dihedral angles and how many torsions are this angle 
                    # current parameter for this 'fold'
                    prm = torprm['prmdict'][nfold]
                    # not called by 'eval'
                    if parms is not 'eval':
                        # get prm from parms array
                        prm = parms[torprm['prmdict'][nfold]]
                    if 'firstphaseprmindex' in torprm.keys() and nfold==1 and parms!='eval':
                        prmindex=torprm['firstphaseprmindex']
                        offset=parms[prmindex]
                    elif 'firstphaseprmindex' in torprm.keys() and nfold==1 and parms=='eval':
                        offset=torprm['firstphaseprmindex']
                    else:
                        offset=poltype.foldoffsetlist[nfold-1]
                    tor_energy += tor_func_term (poltype,prm, ang, nfold, clscnt, torgen.rads(poltype,clsangle),torgen.rads(poltype,offset))
        tor_energy_array[j]=tor_energy

    if parms is 'eval' and 'offset' in torprm:
        offset = torprm['offset']
    if parms is not 'eval':
        offset = parms[-1]
    if type(offset)==list:
        if len(tor_energy_array)==len(offset):
            tor_energy_array += offset
    return numpy.array(tor_energy_array)

def compute_qm_tor_energy(poltype,torset,mol,flatphaselist):
    """
    Intent: Store the QM Energies (vs. Dihedral Angle) found in 'gen_torsion' in a list
    Input:
        a: atom 1 in the torsion of interest
        b: atom 2 in the torsion of interest
        c: atom 3 in the torsion of interest
        d: atom 4 in the torsion of interest
        startangle: current or initial dihedral angle
        phase_list: list of phase offsets. default is 0-360 in intervals of 30
    Output:
        list(rows[1]): QM energies
        list(rows[0]): Dihedral angles
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms
    Description: Read in the *-sp-*.log files created in 'gen_torsion', find the energy
    values, and store them in a list.
    """

    energy_list = []
    angle_list = []
    WBOarray=[]
    for phaseangles in flatphaselist:
        prefix='%s-sp-'%(poltype.molecprefix)
        postfix='.log' 
        minstrctfname,angles=torgen.GenerateFilename(poltype,torset,phaseangles,prefix,postfix,mol)
        if not os.path.exists(minstrctfname): # if optimization failed then SP file will not exist
            tor_energy=None
            WBOvalues=None
        else:
            if poltype.use_gaus:
                WBOmatrix=frag.GrabWBOMatrixGaussian(poltype,minstrctfname,poltype.mol)
            else:
                WBOmatrix=frag.GrabWBOMatrixPsi4(poltype,minstrctfname,poltype.mol)

            WBOvalues=[]
            for tor in torset:
                a,b,c,d=tor[:]
                WBOvalue=WBOmatrix[b-1,c-1]
                WBOvalues.append(WBOvalue)
                
            tmpfh = open(minstrctfname, 'r')
            tor_energy = None
            if not poltype.use_gaus:
                mengi=esp.GrabFinalPsi4Energy(poltype,minstrctfname)
                if mengi==None:
                    tor_energy=None
                else:
                    tor_energy = float(mengi) * poltype.Hartree2kcal_mol
            else:
                for line in tmpfh:
                    if poltype.torspmethod=='MP2':
                        m = re.search(r'EUMP2 =\s+(\-*\d+\.\d+D\+\d+)',line)
                        if not m is None:
                            mengi = m.group(1).replace('D+', 'E+')
                            tor_energy = float(mengi) * poltype.Hartree2kcal_mol
                    else:
                        if 'SCF Done:' in line:
                            linesplit=line.split()
                            result=float(linesplit[4])
                            tor_energy = result* poltype.Hartree2kcal_mol


            tmpfh.close()

        WBOarray.append(WBOvalues)
        energy_list.append(tor_energy)
        angle_list.append(angles)
    rows = zip(*[angle_list, energy_list])
    rows=sorted(rows)
    rows0=list([i[0] for i in rows])
    rows1=list([i[1] for i in rows])
    return rows1,rows0,WBOarray

def compute_mm_tor_energy(poltype,mol,torset,designatexyz,flatphaselist,keyfile = None):
    """
    Intent: Use tinker analyze to find the Pre-fit MM Energy vs. Dihedral Angle profile
    Input:
        a: atom 1 in the torsion of interest
        b: atom 2 in the torsion of interest
        c: atom 3 in the torsion of interest
        d: atom 4 in the torsion of interest
        startangle: current or initial dihedral angle
        phase_list: list of phase offsets. default is 0-360 in intervals of 30
        keyfile: keyfile for tinker analyze
    Output:
        list(rows[1]): MM energies
        list(rows[0]): Dihedral angles
        list(rows[2]): Energy just due to torsion 
    Referenced By: get_qmmm_rot_bond_energy, ooooooot_bond_parms
    Description:
    1. For each phase offset
        a. Restrain the dihedral angle at (startangle + phaseangle)
        b. Run tinker analyze (for this new restraint)
        c. Read in and store the energy
    """
    torse_list = []
    angle_list = []
    energy_list=[] 
    for phaseangles in flatphaselist:
        prefix='%s-opt-' % (poltype.molecprefix)
        postfix='%s.xyz' % (designatexyz) 
        torxyzfname,angles=torgen.GenerateFilename(poltype,torset,phaseangles,prefix,postfix,mol)
        newtorxyzfname=torxyzfname.replace('.xyz','.xyz_2')
        toralzfname = os.path.splitext(torxyzfname)[0] + '.alz'
        tot_energy,tor_energy=GrabTinkerEnergy(poltype,toralzfname)

        energy_list.append(tot_energy)
        torse_list.append(tor_energy)
        angle_list.append(angles)
    if None in energy_list:
        string='Cannot analyze XYZ file for torsion %s'%(str(torset))
        errstr = [string, energy_list,angle_list]

    rows = zip(*[angle_list, energy_list, torse_list])
    rows=sorted(rows)
    rows0=list([i[0] for i in rows])
    rows1=list([i[1] for i in rows])
    rows2=list([i[2] for i in rows])

    return rows1,rows0,rows2


def GrabTinkerEnergy(poltype,toralzfname):
    tot_energy = None
    tor_energy = None
    if os.path.isfile(toralzfname):  
        tmpfh = open(toralzfname, 'r')
        for line in tmpfh:
            m = re.search(r'Potential Energy :\s+(\-*\d+\.\d+)',line)
            if not m is None:
                tot_energy = float(m.group(1))
            m = re.search(r'Torsional Angle\s+(\-*\d+\.\d+)',line)
            if not m is None:
                tor_energy = float(m.group(1))
        tmpfh.close()
    return tot_energy,tor_energy


def find_del_list(poltype,mme_list,current_ang_list):
    """
    Intent: Run through 'mme_list' and remove None objects;
            remove the corresponding angle as well
    Input:
        mme_list: list of MM energies (vs. angle)
        current_ang_list: list of angles
    Output: 
        del_ang_list: List of angles to remove since the MM energy was not found for that angle
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms 
    Description: -
    """
    del_ang_list = []
    for listidx in range(0,len(mme_list)):
        mm_eng = mme_list[listidx]
        if mm_eng is None:
            del_ang_list.append(current_ang_list[listidx])
    return del_ang_list

def sum_xy_list(poltype,x1,y1,x2,y2):
    for xx in x1:
        if xx in x2:
            idx1 = x1.index(xx)
            idx2 = x2.index(xx)
            y2[idx2] = y1[idx1]

def del_tor_from_fit(poltype,dellist, torprmdict,initialprms):
    for keylist in dellist:
        key=keylist[0]
        indicestodelete=[]
        for clskey in torprmdict.keys():
            if clskey==key:
                if 'firstphaseprmindex' in torprmdict[key].keys():
                    prmindex=torprmdict[key]['firstphaseprmindex']
                    indicestodelete.append(prmindex) 
   
                for fold,index in torprmdict[clskey]['prmdict'].items():
                    indicestodelete.append(index)
    
    minidx=min(indicestodelete)
    newinitialprms=[]
    for i in range(len(initialprms)):
        value=initialprms[i]
        if i not in indicestodelete:
            newinitialprms.append(value)
    newinitialprms.append(0) # vertical shift

    del torprmdict[key]
    for clskey in torprmdict.keys():
        for fold,index in torprmdict[clskey]['prmdict'].items():
            if index>minidx:
                torprmdict[clskey]['prmdict'][fold]-=minidx
  
        
    return newinitialprms,torprmdict

def find_least_connected_torsion(poltype,torprmdict,toralreadyremovedlist):
    """
    Find least connected torsion (i.e. the two outer atoms have the highest summed class number)
    """

    least_connected_tor = None
    highest_clssum = 0
    keylist = torprmdict.keys()

    for chkclskey in keylist:
        a,b,c,d = chkclskey.split()
        cur_clssum = int(a) + int(d)
        if (least_connected_tor is None or cur_clssum > highest_clssum) and chkclskey not in toralreadyremovedlist:
            least_connected_tor = chkclskey
            highest_clssum = cur_clssum
    return least_connected_tor

def prune_mme_error(poltype,del_ang_list,*arr_list):
    """
    Intent: Delete the given ids (del_ang_list) in every list in *arr_list
    Input:
        del_ang_list: list of angles to delete
        *arr_list: list of lists; the corresponding elements of each list will be removed
    Output: -
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms
    Description: -
    """
    arr_list_idx = 0
    x_list = arr_list[0]
    for del_ang in del_ang_list:
        if del_ang in x_list:
            del_idx = x_list.index(del_ang)
            for a_list in arr_list:
                del a_list[del_idx]
    return arr_list

def prune_qme_error(poltype,del_ang_list,*arr_list):
    """
    Intent: Delete the given ids (del_ang_list) in every list in *arr_list
    Input:
        del_ang_list: list of angles to delete
        *arr_list: list of lists; the corresponding elements of each list will be removed
    Output: -
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms
    Description: -
    """
    arr_list_idx = 0
    x_list = arr_list[0]
    for del_ang in del_ang_list:
        if del_ang in x_list:
            del_idx = x_list.index(del_ang)
            for a_list in arr_list:
                del a_list[del_idx]
    return arr_list

def get_qmmm_rot_bond_energy(poltype,mol,tmpkey1basename):
    """
    Intent: Form dicts for each torsion in torlist, mapping the torsion class key ('clskey') to 
    an energy profile (dihedral angle vs. energy). 'cls_mm_engy_dict' maps 'clskey' to pre-fit MM 
    calculated energy profiles, 'cls_qm_engy_dict' maps 'clskey' to QM calculated energy profiles
    Input:
        mol: OBMol Structure
        anglist: phase list. default: 0 - 360 in increments of 30
        tmpkey1basename: key file name for tinker
    Output:
        cls_mm_engy_dict: given a class key, this will provide a list of mm energies (vs. angles)
        cls_qm_engy_dict: given a class key, this will provide a list of qm energies (vs. angles)
        cls_angle_dict: given a class key, this provides the angles that the energies
        above are based on
    Referenced By: process_rot_bond_tors
    Description:
    1. For each rotatable bond
        a. Find QM Energy profile
        b. Find MM Energy profile
        c. delete parts of the list where MM energy was not able to be found
        d. Add the profiles found above to 'cls_qm_engy_dict' and 'cls_mm_engy_dict' for this
        clskey
    2. If a clskey has already been looked at more than once, find the average energies of all
       the torsions having the same clskey
    """
    cls_mm_engy_dict = {}
    cls_qm_engy_dict = {}
    cls_angle_dict = {}
    clscount_dict = {}
    for torset in poltype.torlist:
        if torset not in poltype.nonaroringtorsets and len(torset)==2 and poltype.torfit1Drotonly==True and poltype.torfit2Drotonly==False:
           continue
        elif torset not in poltype.nonaroringtorsets and len(torset)==1 and poltype.torfit1Drotonly==True and poltype.torfit2Drotonly==False:
           pass
        elif torset not in poltype.nonaroringtorsets and len(torset)==2 and poltype.torfit1Drotonly==False and poltype.torfit2Drotonly==True:
           pass
        elif poltype.torfit1Drotonly==False and poltype.torfit2Drotonly==False: 
            pass
        else:
            continue
        


        classkeylist=[]
        flatphaselist=poltype.torsettophaselist[tuple(torset)]
        for i in range(len(torset)):
            tor=torset[i]
            a,b,c,d = tor[0:4]
            clskey = torgen.get_class_key(poltype,a,b,c,d)
            classkeylist.append(clskey)
        tup=tuple(classkeylist)
        if tup not in clscount_dict.keys():
            clscount_dict[tup] = 0
            cls_mm_engy_dict[tup] = [0]*len(flatphaselist)
            cls_qm_engy_dict[tup] = [0]*len(flatphaselist)
            cls_angle_dict[tup] = [0]*len(flatphaselist)

        clscount_dict[tup] += 1
        mme_list = []  # MM Energy before fitting to QM torsion energy
        qme_list = []  # QM torsion energy
        # find qm, then mm energies of the various torsion values found for 'tor'
        qme_list,qang_list,WBOarray = compute_qm_tor_energy(poltype,torset,mol,flatphaselist)
        mme_list,mang_list,tor_e_list = compute_mm_tor_energy(poltype,mol,torset,'_postQMOPTprefit',flatphaselist,tmpkey1basename)
        # delete members of the list where the energy was not able to be found 
        del_ang_list = find_del_list(poltype,mme_list,mang_list)
        (cls_angle_dict[tup],cls_mm_engy_dict[tup])=prune_mme_error(poltype,del_ang_list,cls_angle_dict[tup],cls_mm_engy_dict[tup])
        (mang_list,mme_list,qme_list,qang_list,tor_e_list)=prune_mme_error(poltype,del_ang_list,mang_list,mme_list,qme_list,qang_list,tor_e_list)
        del_ang_list = find_del_list(poltype,qme_list,qang_list)
        (cls_angle_dict[tup],cls_mm_engy_dict[tup])=prune_qme_error(poltype,del_ang_list,cls_angle_dict[tup],cls_mm_engy_dict[tup])
        (mang_list,mme_list,qme_list,qang_list,tor_e_list)=prune_qme_error(poltype,del_ang_list,mang_list,mme_list,qme_list,qang_list,tor_e_list)

        cls_qm_engy_dict[tup] = [ runsum+eng for runsum,eng in zip (cls_qm_engy_dict[tup], qme_list)]
        cls_mm_engy_dict[tup] = [ runsum+eng for runsum,eng in zip (cls_mm_engy_dict[tup], mme_list)]
        cls_angle_dict[tup] = mang_list

    # if multiple class keys, take the average
    for tup in clscount_dict:
        cnt = clscount_dict[tup]
        cls_mm_engy_dict[tup] = [eng/cnt for eng in cls_mm_engy_dict[tup]]
        cls_qm_engy_dict[tup] = [eng/cnt for eng in cls_qm_engy_dict[tup]]

    return cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict

def check_cooperative_tor_terms(poltype,clskey, nfold, angle, xvals, torprmdict):
    """ 
    Intent: Check if a cosine term for a torsion is cooperative with the corresponding cosine term of another torsion
    i.e. do they have equal or opposite phase shifts
    Input:
        clskey: class key for the torsion
        nfold: the 'fold' of the cosine term for the torsion
        angle: angle that the torsion is at
        xvals: angle list
        tormprmdict: dictionary containing information about the torsions about a rotatable bond
    Output:
        ck: class key of the torsion that has a cooperative term
    Referenced By: insert_torprmdict
    Description:
    """
    # for each class key
    b=clskey[1]
    c=clskey[2]
    for ck,tp in torprmdict.items():
    # for each angle for the class key
        newb=ck[1]
        newc=ck[2]
        if b==newb and c==newc:
            for phase in tp['phasedict']:
                # if this torsion has the same number of folds
                if nfold in tp['prmdict'] and \
                   not isinstance(tp['prmdict'][nfold], tuple):
                    outer_tor_e = tor_func_term(poltype,1.0,xvals,nfold,1.0,torgen.rads(poltype,phase),
                                      torgen.rads(poltype,poltype.foldoffsetlist[nfold-1]))
                    inner_tor_e = tor_func_term(poltype,1.0,xvals,nfold,1.0,torgen.rads(poltype,angle),
                                      torgen.rads(poltype,poltype.foldoffsetlist[nfold-1]))
                    # sum the profiles, then subtract the min from the list
                    sum_e = outer_tor_e + inner_tor_e
                    sum_e -= min(sum_e)
                    # diff the profiles, subtract the min from the list
                    diff_e = outer_tor_e - inner_tor_e
                    diff_e -= min(diff_e)
                    # if when summed, or diffed, the profile becomes close to constant
                    # i.e. the difference from the max and the min is < 1e-10
                    # then return the class key for this similar torsion
                    if max(sum_e) < 1e-10 or max(diff_e) < 1e-10:
                        return ck
    return None

def insert_torphasedict (poltype,mol, toraboutbnd, torprmdict, initangle,write_prm_dict, keyfilter = None):
    """
    Intent: Adds torsion to be fitted to torprmdict.
    Input:
        mol: An openbabel molecule structure
        toraboutbnd: A list containing the quadruplet of atoms to add
        write_prm_dict: A dict of torsion parameters that don't need fitting
        tormprmdict: dictionary containing information about the torsions about a rotatable bond
        initangle: The torsion angle of the torsion in optimized geometry
        keyfilter: A list of class numbers to allow adding to torprmdict
    Output:
        Modifies appends toraboutbnd to torprmdict or appends to
        write_prm_dict (with 0s for parameters).
    Referenced By: fit_rot_bond_tors
    Description: Adds the torsion 'toraboutbnd' to 'torprmdict'.
    If a torsion with this class key already exists in torprmdict, increase the count
    Adds the current dihedral angle of the torsion to the phasedict of torprmdict
    """
    # quadruplet
    a2,b2,c2,d2 = toraboutbnd
    # create atom structures
    obaa = mol.GetAtom(a2)
    obab = mol.GetAtom(b2)
    obac = mol.GetAtom(c2)
    obad = mol.GetAtom(d2)
    # create a key
    # because it is using symmetry classes instead of atom id's, tpdkey can repeat
    tpdkey = torgen.get_class_key(poltype,a2, b2, c2, d2)

    # if the key passes the keyfilter or if the keyfilter does exist
    # and, the end atoms are not hydrogens
    if (keyfilter is None or keyfilter == tpdkey): # MODIFY THIS LINE TO ALLOW HYDROGEN TORSION
        # current torsion value (normalized by initangle)
        torabangle = round(mol.GetTorsion(obaa,obab,obac,obad) -initangle) % 360
        if tpdkey in torprmdict:
            # increase the count for this tpdkey
            torprmdict[tpdkey]['count'] += 1
            insert_torprmdict_angle(poltype,torabangle, torprmdict[tpdkey]['phasedict'])
        else:
            # set up dict
            torprmdict[tpdkey] = {}
            torprmdict[tpdkey]['count'] = 1
            torprmdict[tpdkey]['phasedict'] = {}
            torprmdict[tpdkey]['prmdict'] = {}
            if len(torprmdict) == 1:
                torprmdict[tpdkey]['offset'] = 1.
            # alter count for this current angle (torabangle)
            insert_torprmdict_angle(poltype,torabangle,torprmdict[tpdkey]['phasedict'])

    # else, force constants are set to 0
    else:
        if tpdkey in poltype.classkeytoinitialprmguess.keys():
            prms=poltype.classkeytoinitialprmguess[tpdkey]
            write_prm_dict[tpdkey] = {1:prms[0], 2:prms[1], 3:prms[2]}
        else: 
            write_prm_dict[tpdkey] = {1:0., 2:0., 3:0.}


def insert_torprmdict(poltype,mol, torprmdict):
    """
    Intent: Initialize the prmdicts in torprmdict 
    Give each torsion intially 3 folds
    Input:
        mol: OBMol object
        tormprmdict: dictionary containing information about the torsions about a rotatable bond
    Output:
        torprmdict is modified
        prmidx: the number of parameters to be fitted for
    Referenced By: fit_rot_bond_tors
    Description: -
    """
    initialprms=[]
    tmpx = numpy.arange ( 0.0, 360.0, 10)
    prmidx = 0
    # for each cls key
    for (chkclskey, torprm) in torprmdict.items():
        # for each parameter in the energy equation for this torsion
        # nfoldlist = [1,2,3]
        if poltype.fitfirsttorsionfoldphase==True:
            torprm['firstphaseprmindex']=prmidx
            prmidx+=1
            initialprms.append(0)
        for nfold in poltype.nfoldlist:
            # init array
            test_tor_energy = numpy.zeros(len(tmpx))
            # for each dihedral angle about this rotatable bond and the number of time it occurs 
            # sum up test_tor_energy
            for (angle, scale) in torprm['phasedict'].items():
                # create a test energy list (0 - 360, 10) starting from this angle
                test_tor_energy += tor_func_term(poltype,1.0, tmpx, nfold, scale, torgen.rads(poltype,angle),torgen.rads(poltype,poltype.foldoffsetlist[nfold-1]))

                # check if another torsion has a similar profile 
                basetorkeys = check_cooperative_tor_terms(poltype,chkclskey, nfold, angle, tmpx, torprmdict)

            # normalize
            test_tor_energy -= min(test_tor_energy)
            
            if torprm['phasedict']:
                if basetorkeys is not None:
                    # if the energy profile for this torsion/fold is not so dissimilar from
                    # another torsion/fold, then take its parameter 
                    # in this case, prmidx is not increased by one
                    torprmdict[chkclskey]['prmdict'][nfold] = torprmdict[basetorkeys]['prmdict'][nfold]


                # will come here first
                elif max(test_tor_energy) > 1e-10:
                    torprm['prmdict'][nfold] = prmidx
                    prmidx += 1
                    if chkclskey in poltype.classkeytoinitialprmguess.keys():
                        prms=poltype.classkeytoinitialprmguess[chkclskey]
                        initialprms.append(prms[nfold-1])
                    else: 
                        initialprms.append(0)

        if not torprmdict[chkclskey]['prmdict']:
            torprmdict[chkclskey]['count'] = 0
            torprmdict[chkclskey]['phasedict'] = {}

        
    prmidx += 1
    initialprms.append(0)
    return prmidx,initialprms

def is_torprmdict_all_empty (poltype,torprmdict):
    """
    Intent: Determines if all torsion parameters for fitting have been eliminated.
    Input:
        torprmdict: a dict containing list of torsions to be fit
    Ouptput:
        results: True if empty, False if not
    Referenced By: fit_rot_bond_tors
    Description: -
    """
    result = True
    for clskey in torprmdict:
        if torprmdict[clskey]['prmdict']:
            return False
    return result


def GenerateBoundaries(poltype,max_amp,refine,initialprms,torprmdict):
    lowerbounds=[]
    upperbounds=[]
    if refine==False: 
        bounds=[-max_amp,max_amp]

    else:
        for chkclskey in torprmdict:

            prms=[]
            for nfold in torprmdict[chkclskey]['prmdict']:
                parm  = initialprms[torprmdict[chkclskey]['prmdict'][nfold]]
                prms.append(parm)
            allzero=True
            for prm in prms:
                value=numpy.abs(prm)
                if value!=0:
                    allzero=False
            if 'firstphaseprmindex' in torprmdict[chkclskey].keys():
                lowerbounds.append(0)
                upperbounds.append(2*numpy.pi)

            for prmidx in range(len(prms)):
                prm=prms[prmidx]
                value=numpy.abs(prm)
                scaledvalue=.3*value
                if allzero==False:
                    if value!=0:
                        lowerbounds.append(prm-scaledvalue)
                        upperbounds.append(prm+scaledvalue)
                    else:
                        lowerbounds.append(prm)
                        upperbounds.append(.2) # cant bound by 0,0
                else:
                    lowerbounds.append(-max_amp)
                    upperbounds.append(max_amp)

        lowerbounds.append(-max_amp) # vertical shift at end
        upperbounds.append(max_amp)

        bounds=[lowerbounds,upperbounds]
    return tuple(bounds)
       

def fit_rot_bond_tors(poltype,mol,cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,clskeyswithbadfits):
    """
    Intent: Uses scipy's optimize.leastsq function to find estimates for the torsion 
    parameters based on energy values found at various angles using qm and mm
    Each rotatable bond is fit for one at a time
    Input:
        mol: OBMol structure
        cls_mm_engy_dict: given a class key, this will provide a list of mm energies (vs. angles)
        cls_qm_engy_dict: given a class key, this will provide a list of qm energies (vs. angles)
        cls_angle_dict: given a class key, this provides the angles that the energies
        above are based on
    Output:
        write_prm_dict: map from class key to parameter information. 
                        Used to write out new key file
        fitfunc_dict: energy profile (using the new parameters) to be plotted
    Referenced By: process_rot_bond_tors
    Description:
    1. Initialize 'fitfunc_dict' and 'write_prm_dict'
    2. For each tor in torlist (essentially, for each rotatable bond) 
    (the fit is done for each rotatable bond one at a time):
        a. Initialize 'torprmdict'
            i. For each torsion about the current rotatable bond, 'torprmdict' maps the torsion
               class key to a set of dictionaries containing information about that torsion.
               The three dictionaries containing information are: 
               *count: # times a torsion with this class key exists 
               (there can be multiple torsions about the rotatable bond that have the same class
               key, for example two H-CC-O 's can exist)
               *prmdict: parameters for this torsion
                         The energy equation for the torsion is a sum of cosines
                         'nfolds' is the number of cosines in the sum 
                         each 'fold' or cosine has a coefficient
                         prmdict contains information on the number of folds and the parameters
                         for each fold
               *phasedict: the various dihedral angles that this torsion currently exists at,
                and how many torsions are at this angle
        b. Get the atoms involved in the main torsion around this rotatable bond
        c. Get the current dihedral angle and the class key
        d. Fill in the phasedict portion of 'torprmdict' by calling method 'insert_torphasedict'
           on each torsion about the current rotatable bond 
        e. Edit 'torprmdict' by calling 'insert_torprmdict'
           'prmidx' now equals: number of parameters to be fit
        f. 'angle_list', 'mm_energy_list', 'qm_energy_list' are all initialized for the current
            classkey
        g. normalize the qm_energy_list and mm_energy_list by subtracting all values by the min 
        h. 'tor_energy_list" is created by subtracting mme from qme and is written to a file
        i. 'nfolds' list (of lists) is found. Should at least initially be [[1,2,3],[1,2,3],...]
            nfolds are the number of force constants per torsion in question
        j. 'max_amp', max - min of tor_energy_list, is found 
        k. 'pzero' is initialized
        l. Remove parameters while # of parameters > # of data points
           This can happen if two torsions are very similar so their parameters can be combined
           into one set
        m. now run optimize.leastsq. Keep rerunning it until the parameters no longer have to be 
           'sanitized', meaning that none of the parameter estimates found by leastsq
           are greater than max_amp
        n. If all of the parameter estimates ended up being deleted, leastsq is rerun, 
           this time fitting for only the main torsion
        o. fill in 'torprmdict' with parameter estimates found by leastsq
        p. write out a plot of the fit
        q. write out the parameter estimates
    """

    fitfunc_dict = {}
    write_prm_dict = {}
    if len(poltype.torlist)==0:
        return write_prm_dict,fitfunc_dict

    # For each rotatable bond 
    for torset in poltype.torlist:
        if torset not in poltype.nonaroringtorsets and len(torset)==2 and poltype.torfit1Drotonly==True:
           continue
        torprmdict = {}
        mm_energy_list2 = [] # MM Energy after fitting
        classkeylist=[]
        refine=False
        flatphaselist=poltype.torsettophaselist[tuple(torset)]

        for i in range(len(torset)):
            tor=torset[i]
            # get the atoms in the main torsion about this rotatable bond
            a,b,c,d = tor[0:4]
            # current torsion value
            torang = mol.GetTorsion(a,b,c,d)
            # class key; ie symmetry classes key
            clskey = torgen.get_class_key(poltype,a,b,c,d)
            classkeylist.append(clskey)
            if clskey in clskeyswithbadfits:
                useweights=True
            else:
                useweights=False
            if clskey in poltype.classkeytoinitialprmguess.keys():
                refine=True
            rotbndkey = '%d %d' % (b, c)
            initangle = mol.GetTorsion(a,b,c,d)
            # Identify all torsion parameters involved with current rotatable bond.
            if rotbndkey in poltype.rotbndlist.keys():
                for toraboutbnd in poltype.rotbndlist[rotbndkey]:
                    # However, initangle is the current angle for 'tor' not for 'toraboutbnd'
                    insert_torphasedict(poltype,mol, toraboutbnd, torprmdict, initangle, write_prm_dict)
            else:
                insert_torphasedict(poltype,mol, tor, torprmdict, initangle, write_prm_dict)

        tup=tuple(classkeylist)
        prmidx,initialprms = insert_torprmdict(poltype,mol, torprmdict)
        #poltype.WriteToLog('number of parameters to fit for '+clskey+' are '+str(prmidx))
        # get all the lists for the current clskey
        angle_list = cls_angle_dict[tup]  # Torsion angle for each corresponding energy
        mm_energy_list = cls_mm_engy_dict[tup]  # MM Energy before fitting to QM torsion energy
        qm_energy_list = cls_qm_engy_dict[tup]  # QM torsion energy
        # 'normalize'
        poltype.WriteToLog('about to normalize') 
        qm_energy_list = [en - min(qm_energy_list) for en in qm_energy_list]
        mm_energy_list = [en - min(mm_energy_list) for en in mm_energy_list]
        if len(qm_energy_list)<round(prmidx*.5): # then might not be great fit any way, too many QM failed
            raise ValueError('Too many QM jobs have failed for '+str(tor)+' '+os.getcwd())
        # Parameterize each group of rotatable bond (identified by
        #  atoms restrained during restrained rotation.
        # tor_energy_list is set as qm - mm

        weightlist=numpy.exp(-numpy.add(numpy.array(qm_energy_list),15)/2.5)

        tor_energy_list = numpy.array([qme - mme for qme,mme in zip(qm_energy_list,mm_energy_list)])
        if useweights==True:
            tor_energy_list=numpy.multiply(tor_energy_list,weightlist)
        tor_energy_list_noweight = [qme - mme for qme,mme in zip(qm_energy_list,mm_energy_list)]
        '''
        txtfname = "%s-fit-" % (poltype.molecprefix) 
        for i in range(len(torset)):
            tor=torset[i]
            a,b,c,d = tor[0:4]
            txtfname+="%d-%d" % (b, c)
            txtfname+='_'
        txtfname=txtfname[:-1]
        txtfname+=".txt"         # create initial fit file, initially it seems to be 2d instead of 3d
        torgen.write_arr_to_file(poltype,txtfname,[tor_energy_list])
        print('about to fit')
        poltype.WriteToLog('about to fit')
        ''' 
        # max amplitude of function
        max_amp = max(tor_energy_list) - min(tor_energy_list)
        boundstup=GenerateBoundaries(poltype,max_amp,refine,initialprms,torprmdict)
        # Remove parameters while # of parameters > # data points
        toralreadyremovedlist=[]
        while prmidx > len(mm_energy_list):
            
            dellist= []
            least_conn_tor = find_least_connected_torsion(poltype,torprmdict,toralreadyremovedlist)
            for nfold in torprmdict[least_conn_tor]['prmdict']:
                dellist.append((least_conn_tor,nfold))
            initialprms,torprmdict = del_tor_from_fit(poltype,dellist,torprmdict,initialprms)
            prmidx=len(initialprms)
            if len(dellist)>1:
                poltype.WriteToLog('torsion cosine terms that are being removed due to having too many parameters to fit '+str(dellist))
                poltype.WriteToLog('number of parameters to fit for '+clskey+' are '+str(prmidx))
                toralreadyremovedlist.append(least_conn_tor)
            
        pzero = initialprms
        # run leastsq until all the parameter estimates are reasonable
        parm_sanitized = False
        bypassrmsd=False
        while not parm_sanitized:
            parm_sanitized = True
            dellist = []
            keylist = list(torprmdict.keys())
            keylist.reverse()
            
            # creating a new function, errfunc
            # p: parameters
            # x: angle list
            # torprmdict: torsion information 
            # y: tor_energy_list 
            errfunc = lambda p, x, z, torprmdict, y: fitfunc(poltype,p, x,z, torprmdict) - y

            # optimize.leastsq is run, found in the scipy library
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html
            # Inputs:
            # errfunc : the callable function that we supply, based on fitfunc
            # This is the function for which the squares are attempted to be minimized
            # pzero : starting parameter estimates
            # args=(rads(numpy.array(angle_list)), torprmdict, tor_energy_list) :
            # These are the arguments that errfunc needs aside from 'p' (which will be supplied by
            # leastsq)
            # full_output : set to True, leastsq returns all optional outputs
            # Outputs:
            # p1 : parameter estimates returned by optimize.leastsq
            # covx : matrix that can be used to find covariance of parameter estimates
            # idict : dictionary containing info about the run
            # msg : string giving cause of failure if one exists
            # ier : an int. if ier is 1,2,3 or 4, a solution was found, else no solution was found
            #p1,covx,idict,msg,ier = optimize.leastsq(errfunc, pzero, args=(rads(numpy.array(angle_list)),torprmdict, tor_energy_list), full_output = True)
            array=optimize.least_squares(errfunc, pzero, jac='2-point', bounds=boundstup, args=(torgen.rads(poltype,numpy.array(angle_list)),torset,torprmdict, tor_energy_list))
            p1=array['x']
            # Remove parameters found by least.sq that aren't reasonable; 
            # remove parameters found that are greater than max_amp
            if refine==False:
                for chkclskey in keylist:
                    for nfold in torprmdict[chkclskey]['prmdict']:
                        # if the param value is greater than max_amp, remove it
                        if abs(p1[torprmdict[chkclskey]['prmdict'][nfold]]) > max_amp:
                            dellist.append((chkclskey,nfold))
                            parm_sanitized = False
                            break
            foldtoparmslist={}
            parmtokey={}
            for chkclskey in keylist:
                for nfold in torprmdict[chkclskey]['prmdict']:
                    if nfold not in foldtoparmslist.keys():
                        foldtoparmslist[nfold]=[]
                    # if the param value is greater than max_amp, remove it
                    parm=p1[torprmdict[chkclskey]['prmdict'][nfold]]
                    foldtoparmslist[nfold].append(parm)
                    parmtokey[parm]=chkclskey
                    if abs(parm) > max_amp:
                        dellist.append((chkclskey,nfold))
                        parm_sanitized = False
                        break
            for fold,parmslist in foldtoparmslist.items():
                combs=list(combinations(parmslist,2))
                for comb in combs:
                    Sum=numpy.abs(comb[0]+comb[1])
                     
                    if Sum<.01: # tolerance for torsions cancelling each other
                       
                       keytodelete=parmtokey[comb[0]]
                       max_amp=20
                       boundstup=GenerateBoundaries(poltype,max_amp,refine,initialprms,torprmdict)
                       bypassrmsd=True
                       parm_sanitized = False
                       for nfold in torprmdict[keytodelete]['prmdict']:
                           dellist.append((keytodelete,nfold)) 

            dellist = list(set(dellist))
            if len(dellist)>0:
                poltype.WriteToLog('torsion cosine terms that are being removed due to unreasonable parameters '+str(dellist))
                poltype.WriteToLog('number of parameters to fit for '+clskey+' are '+str(prmidx))
                initialprms,torprmdict = del_tor_from_fit(poltype,dellist,torprmdict,initialprms)
                prmidx=len(initialprms)
            # new parameter array since prm size may have changed due to deletions
            pzero = initialprms
        # Attempts to insert main torsion type if all are removed
        # Rerur leastsq, this time fitting for the force constants of the main torsion
        if is_torprmdict_all_empty(poltype,torprmdict):
            for i in range(len(torset)):
                tor=torset[i]
                a,b,c,d = tor[0:4]
                rotbndkey = '%d %d' % (b, c)
                toraboutbnd = poltype.rotbndlist[rotbndkey][0]
                insert_torphasedict(poltype,mol, toraboutbnd, torprmdict,initangle, write_prm_dict,keyfilter = clskey)
            prmidx,initialprms = insert_torprmdict(mol, torprmdict)

            pzero = initialprms
            errfunc = lambda p, x, z, torprmdict, y: fitfunc(poltype,p, x, z, torprmdict) - y
            array=optimize.least_squares(errfunc, pzero, jac='2-point', bounds=boundstup, args=(torgen.rads(poltype,numpy.array(angle_list)),torset,torprmdict, tor_energy_list))
            p1=array['x']

        # fill in torprmdict with the parameter estimates
        for chkclskey in torprmdict:
            if 'firstphaseprmindex' in torprmdict[chkclskey].keys():
                prmindex=torprmdict[chkclskey]['firstphaseprmindex']
                prm=p1[prmindex]
                torprmdict[chkclskey]['firstphaseprmindex']=prm

            for nfold in torprmdict[chkclskey]['prmdict']:
                parm  = p1[torprmdict[chkclskey]['prmdict'][nfold]]
                torprmdict[chkclskey]['prmdict'][nfold] = parm
                if numpy.abs(parm)>20:
                    raise ValueError('parameter way to big'+str(parm))
            write_prm_dict[chkclskey] = torprmdict[chkclskey]['prmdict']
            # if not found, set as 0
            if write_prm_dict[chkclskey] == {}:
                if chkclskey in poltype.classkeytoinitialprmguess.keys():
                    prms=poltype.classkeytoinitialprmguess[chkclskey]
                    write_prm_dict[chclskey] = {1:prms[0], 2:prms[1], 3:prms[2]}
                else: 
                    write_prm_dict[chkclskey] = {1:0., 2:0., 3:0.}

            # Check if no arguments were fitted.
            if 'offset' in torprmdict[chkclskey]:
                if isinstance(p1,numpy.ndarray):
                    torprmdict[chkclskey]['offset'] = p1[-1]
                else:
                    torprmdict[chkclskey]['offset'] = p1
        dim=len(cls_angle_dict[tup][0])
        figfname = '%s-fit-' % (poltype.molecprefix)
        for i in range(len(torset)):
            tor=torset[i]
            a,b,c,d = tor[0:4]
            figfname+='%d-%d' % (b,c)
            figfname+='_'
        figfname=figfname[:-1]
        if useweights==True:
            wstring='_useweights'
        else:
            wstring=''
        figfname+=wstring+'.png'
        Sx = numpy.array(cls_angle_dict[tup])
        out=[]


        out.append(cls_angle_dict[tup])
        for clskey in classkeylist:
            fitfunc_dict[clskey] = fitfunc(poltype,'eval',torgen.rads(poltype,Sx),torset,torprmdict,debug=False)
            out.append(fitfunc_dict[clskey])
            def RMSD(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.subtract(fitfunc_dict[clskey],tor_energy_list),c))))
            result=fmin(RMSD,.5)
            minRMSD=RMSD(result[0])

        out.append(tor_energy_list)
        numprms=1 # offset parameter incldued with torsion force constant parameters
        for classkey in torprmdict:
            numprms+= len(torprmdict[classkey]['prmdict'].keys())
        string=' , '.join(list(torprmdict.keys()))
        datapts=len(mm_energy_list)
        print('Torsions being fit',string,'RMSD(QM-MM1)',minRMSD)
        poltype.WriteToLog('Torsions being fit '+string+' RMSD(QM-MM1)'+str(minRMSD))
        if dim==1:
            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111)
            l1, = ax.plot(Sx,fitfunc_dict[clskey],'r',label='Fit')
            if useweights==False:
                lab='QM-MM1'
            else:
                lab='Weighted_(QM-MM1)'
            l2, = ax.plot(Sx,tor_energy_list,'b',label=lab)
            if useweights==True:
                l3, = ax.plot(Sx,tor_energy_list_noweight,'black',label='QM-MM1')

                plt.legend(handles=[l1,l2,l3],loc=9, bbox_to_anchor=(0.5, -0.1), ncol=3)
            else:
                plt.legend(handles=[l1,l2],loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)

            ax.text(0.05, 1.1, 'Torsions Being Fit =%s'%(string), transform=ax.transAxes, fontsize=10,verticalalignment='top')
            ax.text(0, -0.1, 'FoldNum=%s NumPrms=%s DataPts=%s RMSD(fit,QM-MM1),Abs=%s'%(str(len(poltype.nfoldlist)),str(numprms),str(len(mm_energy_list)),round(minRMSD,2)), transform=ax.transAxes, fontsize=10,verticalalignment='bottom')
            fig.savefig(figfname)
        elif dim==2:
            tormat,qmmat,mmmat,idealanglematrix,actualanglematrix=FillInEnergyTensors(poltype,flatphaselist,cls_angle_dict[tup],tor_energy_list,qm_energy_list,mm_energy_list,torset)
            PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,tormat,'QM-MM1 Heatmap (kcal/mol)','QM-MM1_Heatmap.png',numprms,datapts,minRMSD,string)
            PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,qmmat,'QM Heatmap (kcal/mol)','QM_Heatmap.png',numprms,datapts,minRMSD,string)
            PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,mmmat,'MM1 Heatmap (kcal/mol)','MM1_Heatmap.png',numprms,datapts,minRMSD,string)
        #torgen.write_arr_to_file(poltype,txtfname,out)
    return write_prm_dict,fitfunc_dict

def PlotHeatmap(poltype,idealanglematrix,actualanglematrix,matrix,title,figfname,numprms,datapts,minRMSD,textstring=None):
    fig, ax = plt.subplots()
    im = ax.imshow(matrix)

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('', rotation=-90, va="bottom")
    yangles=idealanglematrix[:,0,0]
    yangles=[round(i,1) for i in yangles]
    xangles=idealanglematrix[0,:,1]
    xangles=[round(i,1) for i in xangles]
    # We want to show all ticks...
    ax.set_xticks(numpy.arange(len(xangles)))
    ax.set_yticks(numpy.arange(len(yangles)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(xangles)
    ax.set_yticklabels(yangles)
    ax.set_ylim(len(matrix),-0.5, -0.5) 
    # Loop over data dimensions and create text annotations.
    for i in range(len(yangles)):
        for j in range(len(xangles)):
            energyvalue=str(round(matrix[i,j],1))
            value=actualanglematrix[i,j,:]
            string=str(round(value[0],1))+','+str(round(value[1],1))
            text = ax.text(j, i, energyvalue,ha="center", va="center", color="w")
    
    ax.set_title(title)
    if textstring!=None:
        ax.text(0.05, 1.1, 'Torsions Being Fit =%s'%(textstring), transform=ax.transAxes, fontsize=10,verticalalignment='top')
        ax.text(0, -0.1, 'FoldNum=%s NumPrms=%s DataPts=%s RMSD(fit,QM-MM1),Abs=%s'%(str(len(poltype.nfoldlist)),str(numprms),str(datapts),round(minRMSD,2)), transform=ax.transAxes, fontsize=10,verticalalignment='bottom')
    else:
        ax.text(0.05, 1.1, 'RMSD(MM2,QM)=%s'%(round(minRMSD,2)), transform=ax.transAxes, fontsize=12,verticalalignment='top')



    fig.tight_layout()
    plt.show() 
    fig.savefig(figfname)



def FillInEnergyTensors(poltype,phaseanglearray,actualanglearray,tor_energy_list,qm_energy_list,mm_energy_list,torset):
    phasetensor=poltype.tensorphases[torset]
    shape=list(phasetensor.shape)
    
    shape=shape[:-1]
    tormat=numpy.empty(shape)
    qmmat=numpy.empty(shape)
    mmmat=numpy.empty(shape)
    sqrt=int(numpy.sqrt(shape))
    shape=tuple([sqrt,sqrt])
    idealanglematrix=poltype.idealangletensor[torset]
    actualanglematrix=numpy.empty(phasetensor.shape)
    for i in range(len(phaseanglearray)):
        phasetup=phaseanglearray[i]
        torenergy=tor_energy_list[i]
        qmenergy=qm_energy_list[i]
        mmenergy=mm_energy_list[i]
        angletup=actualanglearray[i]
        indexes=numpy.where(numpy.all(phasetensor == numpy.array(phasetup), axis=-1))
        tormat[indexes]=torenergy
        qmmat[indexes]=qmenergy
        mmmat[indexes]=mmenergy
        actualanglematrix[indexes]=angletup
    tormat=tormat.reshape(shape)
    qmmat=qmmat.reshape(shape)
    mmmat=mmmat.reshape(shape)
    shape=list(shape)
    shape.append(2)
    shape=tuple(shape)
    idealanglematrix=idealanglematrix.reshape(shape)
    actualanglematrix=actualanglematrix.reshape(shape)
    
    return tormat,qmmat,mmmat,idealanglematrix,actualanglematrix


def write_key_file(poltype,write_prm_dict,tmpkey1basename,tmpkey2basename):
    """
    Intent: Output the new key file based on parameters in write_prm_dict
    """
    tmpfh1 = open(tmpkey1basename, "r")
    tmpfh2 = open(tmpkey2basename, "w")
    fitline="# Fitted torsion" +"\n"
    for line in tmpfh1:
        if '#' in line:
            continue
        m = re.search(r'torsion',line)
        if m is None:
            tmpfh2.write(line)
        else:
            linarr = line.split()
            cl = linarr[1:5]
            clskey = ' '.join(cl) # Order is fine (read from *.prm file)
            torline = line
            torvals = [float(ele) for ele in linarr[5:24:3]]
#            if clskey in write_prm_dict and \
#               torvals == [0.]*len(torvals):
            if clskey in write_prm_dict:
                torline = ' torsion %7s %4s %4s %4s   ' % (cl[0],cl[1],cl[2],cl[3])
                for (nfold, prm) in write_prm_dict[clskey].items():
                    torline += ' %7.3f %.1f %d' % (prm,poltype.foldoffsetlist[nfold - 1], nfold)
                torline += '\n'
            tmpfh2.write(fitline)
            tmpfh2.write(torline)
    tmpfh1.close()
    tmpfh2.close()

def eval_rot_bond_parms(poltype,mol,fitfunc_dict,tmpkey1basename,tmpkey2basename,count,clskeyswithbadfits):
    """
    Intent: 
    For each torsion whose parameters were fit for:
        Using the new parameters, find the new MM Energy vs. Dihedral Angle Profile
        Output the MM Energy (Pre-fit), MM Energy (Post-fit), and QM Energy profiles as 
        plots in the file: *energy*.png
    Ideally the profiles of MM Energy (Post-fit) will be much closer to the QM Energy profiles
    than the MM Energy (Pre-fit) profiles were. Look at the *png post running poltype to confirm.
    Input:
        mol: OBMol structure
        anglelist: dihedral angles, default 0-360, increments of 30
        fitfunc_dict: energy profile
        tmpkey1basename: Old key file, with old torsion parameters
        tmpkey2basename: New key file, with new torsion parameters
    Output:
        *energy*.png:
    Referenced By: process_rot_bond_tors
    Description:
    1. For each torsion whose parameters have been fit for (for each tor in torlist):
        a. Get each energy profile (MM pre, MM post, QM)
        b. Plot the profiles
    """
    # for each main torsion
    for torset in poltype.torlist:
        if torset not in poltype.nonaroringtorsets and len(torset)==2 and poltype.torfit1Drotonly==True:
            continue


        flatphaselist=poltype.torsettophaselist[tuple(torset)]
        for tor in torset:
            a,b,c,d = tor[0:4]
            key=str(b)+' '+str(c)
            anginc=poltype.rotbndtoanginc[key]
            maxrange=poltype.rotbndtomaxrange[key]
            anglelist=range(0,maxrange,anginc)
            torang = mol.GetTorsion(a,b,c,d)
            atmnuma = mol.GetAtom(a).GetAtomicNum()
            atmnumd = mol.GetAtom(d).GetAtomicNum()

            # clskey
            clskey = torgen.get_class_key(poltype,a, b, c, d)
            if clskey not in clskeyswithbadfits and count>0:
                continue
            if clskey not in fitfunc_dict.keys():
                continue
        mm_energy_list = []
        mm_energy_list2 = []
        qm_energy_list = []

        # get the qm energy profile
        qm_energy_list,qang_list,WBOarray = compute_qm_tor_energy(poltype,torset,mol,flatphaselist)
        tmpkeyfname = 'tmp.key'
        shutil.copy(tmpkey1basename, tmpkeyfname)
        # get the original mm energy profile
        mm_energy_list,mang_list,tor_e_list = compute_mm_tor_energy(poltype,mol,torset,'_postQMOPTprefit',flatphaselist,tmpkeyfname)
        # get the new mm energy profile (uses new parameters to find energies)
        mm2_energy_list,m2ang_list,tor_e_list2 = compute_mm_tor_energy(poltype,mol,torset,'_postQMOPTpostfit',flatphaselist,tmpkey2basename)
        # remove angles for which energy was unable to be found
        del_ang_list = find_del_list(poltype,mm_energy_list,mang_list)
        (mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2,WBOarray)=prune_mme_error(poltype,del_ang_list,mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2,WBOarray)
        del_ang_list = find_del_list(poltype,qm_energy_list,qang_list)
        (mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2,WBOarray)=prune_qme_error(poltype,del_ang_list,mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2,WBOarray)
        del_ang_list = find_del_list(poltype,mm2_energy_list,m2ang_list)
        (mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2,WBOarray)=prune_qme_error(poltype,del_ang_list,mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2,WBOarray)
        

        # normalize profiles
        qm_energy_list = [en - min(qm_energy_list) for en in qm_energy_list]
        mm_energy_list = [en - min(mm_energy_list) for en in mm_energy_list]
        mm2_energy_list = [en - min(mm2_energy_list) for en in mm2_energy_list]
        # find the difference between the two energy due to torsion profiles 
        tordif_list = [e2-e1 for (e1,e2) in zip(tor_e_list,tor_e_list2)]
        # normalize
        tordif_list = [en - min(tordif_list) for en in tordif_list]
        # find the difference between the two mm energy profiles
        tordifmm_list = [e1+e2 for (e1,e2) in zip (tordif_list,mm_energy_list)]
        tordifmm_list = [en - min(tordifmm_list) for en in tordifmm_list]
        # TBC
        ff_list = [aa+bb for (aa,bb) in zip(mm_energy_list,fitfunc_dict[clskey])]
        weight=numpy.exp(-numpy.add(numpy.array(qm_energy_list),15)/2.5)

        if len(ff_list)==len(mm2_energy_list):
            def RMSDW(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.multiply(numpy.subtract(mm2_energy_list,qm_energy_list),weight),c))))

            resultW=fmin(RMSDW,.5)
            minRMSDW=RMSDW(resultW[0])
            def RMSD(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.subtract(mm2_energy_list,qm_energy_list),c))))
            result=fmin(RMSD,.5)
            minRMSD=RMSD(result[0])
        # output the profiles as plots

        figfname = "%s-energy-" % (poltype.molecprefix)
        for tor in torset:
            a,b,c,d = tor[0:4]
            figfname+="%d-%d" % (b,c)
            figfname+='_'
        figfname=figfname[:-1]
        if count>0:
            wstring='_usewights'
        else:
            wstring=''
        figfname+=wstring+'.png'
        dim=len(mang_list[0])
        datapts=len(mm_energy_list)
        numprms=None 
        print('Torset',torset,'RMSD(MM2,QM)',minRMSD)
        poltype.WriteToLog('Torset'+str(torset)+' RMSD(MM2,QM) '+str(minRMSD))
        if dim==1: 
            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111)
            # energy profiles: mm (pre-fit), mm (post-fit), qm
            line1, =ax.plot(mang_list,mm_energy_list,'g',label='MM1 (prefit)')
            line2, =ax.plot(m2ang_list,mm2_energy_list,'r',label='MM2 (postfit)')
            line3, =ax.plot(qang_list,qm_energy_list,'b',label='QM')
            if count>0:
                ax.text(0.05, 1.1, 'RMSD_Weighted(MM2,QM)=%s'%(round(minRMSDW,2)), transform=ax.transAxes, fontsize=12,verticalalignment='top')
            else:
                ax.text(0.05, 1.1, 'RMSD(MM2,QM)=%s'%(round(minRMSD,2)), transform=ax.transAxes, fontsize=12,verticalalignment='top')

            # mm + fit
            line4, =ax.plot(mang_list,ff_list,'md-',label='MM1+Fit')
            ax2=ax.twinx()
            # make a plot with different y-axis using second axis object
            line5, =ax2.plot(qang_list,WBOarray,'y',label='WBO')
            ax2.set_ylabel("WBO",color="blue",fontsize=14)
            ax.set_xlabel('Dihedral Angle')
            ax.set_ylabel('SP Energy (kcal/mol)')
            plt.legend(handles=[line1,line2,line3,line4,line5],loc=9, bbox_to_anchor=(0.5, -0.1), ncol=5)

            fig = plt.gcf()
            plt.show()
            fig.savefig(figfname)
        elif dim==2:
            mmmat,mm2mat,fmat,idealanglematrix,actualanglematrix=FillInEnergyTensors(poltype,flatphaselist,mang_list,mm_energy_list,mm2_energy_list,ff_list,torset)
            PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,mm2mat,'MM2 Heatmap (kcal/mol)','MM2_Heatmap.png',numprms,datapts,minRMSD)
            PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,fmat,'MM1+Fit Heatmap (kcal/mol)','MM1+Fit_Heatmap.png',numprms,datapts,minRMSD)



        txtfname = figfname.replace('.png','.txt')
 
        out=[]
        out.append(mang_list)
        out.append(mm_energy_list)
        out.append(mm2_energy_list)
        out.append(qm_energy_list)
        out.append(tordif_list)
        torgen.write_arr_to_file(poltype,txtfname,out)
        if count>0: # use weighted RMSD if already failed at fitting
            RMSD=minRMSDW
        else:
            RMSD=minRMSD
        if float(RMSD)>poltype.maxtorRMSPD:
            poltype.WriteToLog('RMSPD of QM and MM torsion profiles is high, RMSPD = '+ str(RMSD)+' Tolerance is '+str(poltype.maxtorRMSPD)+' kcal/mol ')
            clskeyswithbadfits.append(clskey)
            if poltype.suppresstorfiterr==False and count>0 and bypassrmsd==False:
                raise ValueError('RMSPD of QM and MM torsion profile is high, tried fitting to minima and failed, RMSPD = '+str(RMSD))
    return clskeyswithbadfits
                 


def gen_toromit_list(poltype):
    """
    Intent: if 'omittorsion2' is True, read in the *.toromit file to see which torsions 
    should not be scanned for
    Input: *.toromit is read in
    Output: 
        toromit_list: list of torsions that scanning should be omitted for
    Referenced By: main
    Description: Read in file, append information to toromit_list
    """
    toromitf = open(poltype.molecprefix+".toromit")
    for l in toromitf:
        poltype.toromit_list.append(sorttorsion([int(l.split()[0]), int(l.split()[1]), int(l.split()[2]), int(l.split()[3])]))
    toromitf.close()

def sorttorsion(poltype,keylist):
    """
    Intent: Sort the torsion key by sorting the two outer terms and then two inner terms
    (e.g. 4-2-1-3 -> 3-1-2-4)
    Input: 
        keylist: torsion key to be sorted
    Output:
        keylist is updated
    Referenced By: get_torlist, gen_toromit_list
    Description: -
    """
    if(keylist[1] > keylist[2] or (keylist[1] == keylist[2] and keylist[0] > keylist[3])):
        temp1 = keylist[1]
        keylist[1] = keylist[2]
        keylist[2] = temp1 
        temp2 = keylist[3]
        keylist[3] = keylist[0]
        keylist[0] = temp2
    return keylist

# Fit torsion parameters for rotatable bonds
def process_rot_bond_tors(poltype,mol):
    """
    Intent: Fit torsion parameters for torsions about rotatable bonds 
    Input:
        mol: OBMol structure
    Output:
        *.key_5 is written out, with updated torsion parameters
    Referenced By: main
    Description:
    1. Get the QM Energy vs. Dihedral Angle and (Initial/Pre-Fit) MM Energy vs. Dihedral Angle 
       profiles for each rotatable bond. 
       Store these profiles in 'cls_qm_engy_dict' and 'cls_mm_engy_dict'.
    2. Use these profiles to fit for the torsion parameters by calling 'fit_rot_bond_tors'
    3. Evaluate the new parameters output by the fitting fuction and output informational plots 
       by calling 'eval_rot_bond_parms'
    4. Write out the new keyfile (*.key_5) with these new torsion parameters
    """

    #create list from 0 - 360 in increments of 30
    #anglist = range(0,360,30)
    tordir = 'qm-torsion'
    tmpkey1basename = 'tinker.key'
    tmpkey2basename = 'tinker.key_2'
    tmpkey1fname = tordir + '/' + tmpkey1basename
    assert os.path.isdir(tordir), \
       "ERROR: Directory '%s' does not exist" % (tordir) +' '+os.getcwd()
    # copy *.key_4 to the directory qm-torsion
    shutil.copy(poltype.key4fname, tmpkey1fname)
    # change directory to qm-torsion
    os.chdir(tordir)

    # Group all rotatable bonds with the same classes and identify
    # the torsion parameters that need to be fitted.

    # For each rotatable bond, get torsion energy profile from QM
    # and MM (with no rotatable bond torsion parameters)
    # Get QM and MM (pre-fit) energy profiles for torsion parameters
    DecomposeTorsionTorsion(poltype,mol)
    cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict = get_qmmm_rot_bond_energy(poltype,mol,tmpkey1basename)
    # if the fit has not been done already
    clskeyswithbadfits=[]
    count=0
    while 1:
        # do the fit
        if count==2:
            break # dont redo fitting forever
        write_prm_dict,fitfunc_dict = fit_rot_bond_tors(poltype,mol,cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,clskeyswithbadfits)
        # write out new keyfile
        write_key_file(poltype,write_prm_dict,tmpkey1basename,tmpkey2basename)
        if count>0:
            rmstr='rm *post*post*'
            os.system(rmstr) # delete previous files for post fitting
        if len(poltype.torlist)!=0:
            PostfitMinAlz(poltype,tmpkey2basename,'')
        # evaluate the new parameters
        clskeyswithbadfits=eval_rot_bond_parms(poltype,mol,fitfunc_dict,tmpkey1basename,tmpkey2basename,count,clskeyswithbadfits)
        if len(clskeyswithbadfits)==0:
            break
        count+=1
   
    if poltype.torfit1Drotonly==True:
        poltype.torfit1Drotonly=False
        poltype.torfit2Drotonly=True 
        cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict = get_qmmm_rot_bond_energy(poltype,mol,tmpkey1basename)
        print('about to tor torspline',flush=True)
        poltype.WriteToLog('tor tor spline')
        PrepareTorTorSplineInput(poltype,cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,mol)

    shutil.copy(tmpkey2basename,'../' + poltype.key5fname)
    os.chdir('..')

def PostfitMinAlz(poltype,keybasename,keybasepath):
    for outputlog in poltype.optoutputtotorsioninfo.keys():
        term,error=poltype.CheckNormalTermination(outputlog)
        [torset,optmol,variabletorlist,phaseangles,cartxyzname,bondtopology]=poltype.optoutputtotorsioninfo[outputlog]
        if term==True:    
            if not poltype.use_gaus:
                cartxyz,torxyzfname=torgen.tinker_minimize_analyze_QM_Struct(poltype,torset,optmol,variabletorlist,phaseangles,cartxyzname,poltype.torsionrestraint,'_postQMOPTpostfit',keybasename,keybasepath,bondtopology)
            else:
                cartxyz,torxyzfname=torgen.tinker_minimize_analyze_QM_Struct(poltype,torset,optmol,variabletorlist,phaseangles,outputlog,poltype.torsionrestraint,'_postQMOPTpostfit',keybasename,keybasepath,bondtopology)

def DecomposeTorsionTorsion(poltype,optmol):
    torstoadd=[]
    for torset in poltype.torlist:
        if torset not in poltype.nonaroringtorsets and len(torset)==2:
            for toridx in range(len(torset)):
                tor=torset[toridx]
                torstoadd.append(tuple([tor]))
                poltype.torsettofilenametorset[tuple([tor])]=torset
                poltype.torsettotortorindex[tuple([tor])]=toridx
                              
    for torset in torstoadd:
        poltype.torlist.append(torset)
    poltype.torfit1Drotonly=True

def PrepareTorTorSplineInput(poltype,cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,mol):
    temp=open(poltype.key5fname,'a')
    for torset in poltype.torlist:
        if torset not in poltype.nonaroringtorsets and len(torset)==2:
            pass
        else:
            continue
        classkeylist=[]
        for i in range(len(torset)):
            tor=torset[i]
            # get the atoms in the main torsion about this rotatable bond
            a,b,c,d = tor[0:4]
            # current torsion value
            torang = mol.GetTorsion(a,b,c,d)
            # class key; ie symmetry classes key
            clskey = torgen.get_class_key(poltype,a,b,c,d)
            classkeylist.append(clskey)

        tup=tuple(classkeylist)
        angle_list = cls_angle_dict[tup]  # Torsion angle for each corresponding energy
        mm_energy_list = cls_mm_engy_dict[tup]  # MM Energy before fitting to QM torsion energy
        qm_energy_list = cls_qm_engy_dict[tup]  # QM torsion energy
        qm_energy_list = [en - min(qm_energy_list) for en in qm_energy_list]
        mm_energy_list = [en - min(mm_energy_list) for en in mm_energy_list]
        tor_energy_list = numpy.array([qme - mme for qme,mme in zip(qm_energy_list,mm_energy_list)])
        flatphaselist=poltype.torsettophaselist[tuple(torset)]
        tormat,qmmat,mmmat,idealanglematrix,actualanglematrix=FillInEnergyTensors(poltype,flatphaselist,cls_angle_dict[tup],tor_energy_list,qm_energy_list,mm_energy_list,torset)
        firstclskey=classkeylist[0]
        secondclskey=classkeylist[1]
        firstclskeysplit=firstclskey.split()
        secondclskeysplit=secondclskey.split()
        for cls in secondclskeysplit:
            if cls in firstclskeysplit:
                continue
            else:
                firstclskeysplit.append(cls)
        tortorclskey=' '.join(firstclskeysplit)
        print('torset',torset,flush=True)
        print('tormat',tormat)
        print('idealanglematrix',idealanglematrix,flush=True)
        print('actualanglematrix',actualanglematrix,flush=True)
        print('flatphaselist',flatphaselist,flush=True)
        firstanglerow=idealanglematrix[0,:]
        firstrow=tormat[0,:]
        print('firstanglerow',firstanglerow,flush=True)
        newarray=numpy.array([firstanglerow])
        print('idealanglematrix shape',idealanglematrix.shape,'newarray.shape',newarray.shape)
        idealanglematrix=numpy.vstack((idealanglematrix,newarray))
        tormat=numpy.vstack((tormat,numpy.array([firstrow])))
        firstanglecol=idealanglematrix[:,0]
        firstcol=tormat[:,0]
        print('firstanglecol',firstanglecol,flush=True)
        newarray=numpy.array([firstanglecol])
        print('idealanglematrix shape',idealanglematrix.shape,'newarray.shape',newarray.shape)
        idealanglematrix=numpy.hstack((idealanglematrix,newarray))
        tormat=numpy.hstack((tormat,numpy.array([firstcol])))
        print('idealanglematrix',idealanglematrix,flush=True)
        print('tormat',tormat,flush=True)
        rowpts=idealanglematrix.size[0]
        colpts=idealanglematrix.size[1]
        tortorline='tortor '+tortorclskey+' '+str(rowpts)+' '+str(colpts) +'\n'
        temp.write(tortorline)
        for i in range(len(tormat)):
            erow=tormat[i] 
            anglerow=idealanglematrix[i]
            for j in range(len(erow)):
                ecol=erow[j]
                anglecol=anglerow[j]
                line=str(anglecol[0])+' '+str(anglecol[1])+' '+str(ecol)+'\n'
                temp.write(line)
    temp.close()
