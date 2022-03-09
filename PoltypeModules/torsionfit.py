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
from scipy.interpolate import interp1d
import databaseparser as db
from copy import deepcopy
from rdkit import Chem
import databaseparser


def CheckGeometricRestraintEnergy(poltype,alzfile):
    temp=open(alzfile,'r')
    results=temp.readlines()
    temp.close()
    energy=0
    for line in results:
        if 'Geometric Restraints' in line:
            linesplit=line.split()
            energy=float(linesplit[2])
    return energy


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


def fitfunc (poltype,parms, x,torset, torprmdict,keyonlylist=None,nfoldonlylist=None,debug = False):
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
    for clskey,torprm in torprmdict.items():
        if keyonlylist!=None:
            if clskey not in keyonlylist:
                continue
        for j in range(len(x)):
            angtup=x[j]
            tor_energy=0
            for i in range(len(angtup)):
                ang=angtup[i]
                for nfold in torprm['prmdict']:
                    if nfoldonlylist!=None:
                        if nfold not in nfoldonlylist:
                            continue
                    # for each nfold for this torsion (# of parameters)
                    for clsangle, clscnt in torprm['phasedict'].items():
                        # current dihedral angles and how many torsions are this angle 
                        # current parameter for this 'fold'
                        prm = torprm['prmdict'][nfold]
                        # not called by 'eval'
                        evaluate=True
                        if parms is not 'eval':
                            # get prm from parms array
                            prm = parms[torprm['prmdict'][nfold]]
                            evaluate=False
                        if 'firstphaseprmindex' in torprm.keys() and nfold==1 and evaluate==False:
                            prmindex=torprm['firstphaseprmindex']
                            offset=parms[prmindex]
                        elif 'firstphaseprmindex' in torprm.keys() and nfold==1 and evaluate==True:
                            offset=torprm['firstphaseprmindex']
                        else:
                            offset=poltype.foldoffsetlist[nfold-1]
                        tor_energy += tor_func_term (poltype,prm, ang, nfold, clscnt, torgen.rads(poltype,clsangle),torgen.rads(poltype,offset))
            tor_energy_array[j]+=tor_energy

    if parms is 'eval' and 'offset' in torprm:
        offset = torprm['offset']
    if parms is not 'eval':
        offset = parms[-1]
    if type(offset)==list:
        if len(tor_energy_array)==len(offset):
            tor_energy_array += offset
    return numpy.array(tor_energy_array)


def CheckIfLogFileUsingGaussian(poltype,f):
    use_gaus=False
    temp=open(f,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Entering Gaussian System' in line:
            use_gaus=True
            break 
    return use_gaus




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
    phaseangle_list=[]
    for phaseangles in flatphaselist:
        prefix='%s-sp-'%(poltype.molecprefix)
        postfix='.log' 
        minstrctfname,angles=torgen.GenerateFilename(poltype,torset,phaseangles,prefix,postfix,mol)
        if not os.path.exists(minstrctfname): # if optimization failed then SP file will not exist
            
            tor_energy=None
            WBOvalues=None
        else:
            use_gaus=False
            use_gaus=CheckIfLogFileUsingGaussian(poltype,minstrctfname)
            if use_gaus:
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
            if not use_gaus:
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
        phaseangle_list.append(phaseangles)
    nonecount=energy_list.count(None)
    normalpts=len(energy_list)-nonecount
    if torset in poltype.torsionsettonumptsneeded.keys():
        prmnum=poltype.torsionsettonumptsneeded[torset]
        if normalpts<prmnum:
            raise ValueError('Too many missing QM SP energy values for torsion set = '+str(torset)+' , need '+str(prmnum)+' points') 
    rows = zip(*[angle_list, energy_list])
    energytophaseangle=dict(zip(energy_list,phaseangle_list))
    rows=sorted(rows)
    rows0=list([i[0] for i in rows])
    rows1=list([i[1] for i in rows])
    return rows1,rows0,WBOarray,energytophaseangle

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


def ReadAnglesFromOutputFile(poltype,torset,newtorxyzfname):
    cartxyz=torgen.ConvertTinktoXYZ(poltype,newtorxyzfname,newtorxyzfname.replace('.xyz','_cart.xyz'))
    angles=[]
    themol = opt.load_structfile(poltype,cartxyz)
    for tor in torset:
        a,b,c,d=tor[:]
        torang = themol.GetTorsion(a,b,c,d)
        if torang<0:
            torang+=360
        angles.append(torang)
   
    return angles



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
        geom=CheckGeometricRestraintEnergy(poltype,toralzfname)
        tor_energy=tor_energy-geom
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

def prune_mme_error(poltype,indicesremoved,*arr_list):
    """
    Intent: Delete the given ids (del_ang_list) in every list in *arr_list
    Input:
        del_ang_list: list of angles to delete
        *arr_list: list of lists; the corresponding elements of each list will be removed
    Output: -
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms
    Description: -
    """

    x_list = arr_list[0]
    temp=[]
    new_arr_list=deepcopy(arr_list)
    indicesremoved.sort(reverse=True)
    for del_ang_idx in indicesremoved:
        temp.append(del_ang_idx)
        for a_list_idx in range(len(arr_list)):
            a_list=arr_list[a_list_idx]
            del new_arr_list[a_list_idx][del_ang_idx]
    return new_arr_list,temp

    
def prune_qme_error(poltype,indicesremoved,*arr_list):
    """
    Intent: Delete the given ids (del_ang_list) in every list in *arr_list
    Input:
        del_ang_list: list of angles to delete
        *arr_list: list of lists; the corresponding elements of each list will be removed
    Output: -
    Referenced By: get_qmmm_rot_bond_energy, eval_rot_bond_parms
    Description: -
    """
    x_list = arr_list[0]
    temp=[]
    new_arr_list=deepcopy(arr_list)
    indicesremoved.sort(reverse=True)
    for del_ang_idx in indicesremoved:
        temp.append(del_ang_idx)
        for a_list_idx in range(len(arr_list)):
            a_list=arr_list[a_list_idx]
            del new_arr_list[a_list_idx][del_ang_idx]
    return new_arr_list,temp

    
def FindRemovedIndices(poltype,ang_list,del_ang_list,indicesremoved=None):
    if indicesremoved!=None:
        pass
    else:
       indicesremoved=[]
    for ang in del_ang_list:
        index=ang_list.index(ang)
        if index not in indicesremoved:
            indicesremoved.append(index)
    return indicesremoved




def get_qmmm_rot_bond_energy(poltype,mol,tmpkey1basename,fileprefix):
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
    cls_mm_engy_dict_unmodified = {}
    cls_qm_engy_dict = {}
    cls_qm_engy_dict_unmodified = {}
    cls_angle_dict = {}
    cls_angle_dict_unmodified = {}
    clscount_dict = {}
    classkeylisttoindicesalreadydicremoved={}
    classkeylisttoindicesremoved={}
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
        
        indicesremoved=[]
        indicesalreadyremoved=[]
        classkeylist=[]
        flatphaselist=poltype.torsettophaselist[tuple(torset)]
        qme_list,qang_list,WBOarray,energytophaseangle = compute_qm_tor_energy(poltype,torset,mol,flatphaselist)
        mme_list,mang_list,tor_e_list = compute_mm_tor_energy(poltype,mol,torset,fileprefix,flatphaselist,tmpkey1basename)
        if len(poltype.onlyfittorstogether)!=0:
            torset=tuple(poltype.onlyfittorstogether)

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
            cls_angle_dict_unmodified[tup] = [0]*len(flatphaselist)
            cls_mm_engy_dict_unmodified[tup] = [0]*len(flatphaselist)
            cls_qm_engy_dict_unmodified[tup] = [0]*len(flatphaselist)
        if tup not in classkeylisttoindicesalreadydicremoved.keys():
            classkeylisttoindicesalreadydicremoved[tup]=[]
        if tup not in classkeylisttoindicesremoved.keys():
            classkeylisttoindicesremoved[tup]=[]

        clscount_dict[tup] += 1
        # find qm, then mm energies of the various torsion values found for 'tor'
        originalmang_list=mang_list.copy()
        originalqang_list=qang_list.copy()
        originalmme_list=mme_list.copy()
       
        originalqme_list=qme_list.copy()
        originaltor_e_list=tor_e_list.copy()
        # delete members of the list where the energy was not able to be found 
        del_ang_list = find_del_list(poltype,mme_list,mang_list)

        classkeylisttoindicesremoved[tup]=FindRemovedIndices(poltype,originalmang_list,del_ang_list,indicesremoved=classkeylisttoindicesremoved[tup])

        array=prune_mme_error(poltype,classkeylisttoindicesremoved[tup],originalmang_list,originalmme_list,originalqme_list,originalqang_list,originaltor_e_list)
        (mang_list,mme_list,qme_list,qang_list,tor_e_list)=array[:-1][0]

        del_ang_list = find_del_list(poltype,qme_list,qang_list)

        classkeylisttoindicesremoved[tup]=FindRemovedIndices(poltype,originalqang_list,del_ang_list,indicesremoved=classkeylisttoindicesremoved[tup])

        array=prune_qme_error(poltype,classkeylisttoindicesremoved[tup],originalmang_list,originalmme_list,originalqme_list,originalqang_list,originaltor_e_list)
        (mang_list,mme_list,qme_list,qang_list,tor_e_list)=array[:-1][0]
        classkeylisttoindicesalreadydicremoved[tup]=classkeylisttoindicesremoved[tup]
        cls_qm_engy_dict[tup] = [ runsum+eng for runsum,eng in zip (cls_qm_engy_dict[tup], qme_list)]
        cls_mm_engy_dict[tup] = [ runsum+eng for runsum,eng in zip (cls_mm_engy_dict[tup], mme_list)]
        cls_mm_engy_dict_unmodified[tup]=originalmme_list
        cls_qm_engy_dict_unmodified[tup]=originalqme_list
        cls_angle_dict[tup] = mang_list
        cls_angle_dict_unmodified[tup]=originalmang_list
    # if multiple class keys, take the average
    for tup in clscount_dict:
        cnt = clscount_dict[tup]
        cls_qm_engy_dict_unmodified[tup]=ConvertNoneToZero(poltype,cls_qm_engy_dict_unmodified[tup])
        cls_mm_engy_dict_unmodified[tup]=ConvertNoneToZero(poltype,cls_mm_engy_dict_unmodified[tup])
        
        cls_mm_engy_dict[tup] = [eng/cnt for eng in cls_mm_engy_dict[tup]]
        cls_qm_engy_dict[tup] = [eng/cnt for eng in cls_qm_engy_dict[tup]]
        cls_mm_engy_dict_unmodified[tup] = [eng/cnt for eng in cls_mm_engy_dict_unmodified[tup]]
        cls_qm_engy_dict_unmodified[tup] = [eng/cnt for eng in cls_qm_engy_dict_unmodified[tup]]

    return cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,classkeylisttoindicesremoved,cls_angle_dict_unmodified,cls_mm_engy_dict_unmodified,cls_qm_engy_dict_unmodified,classkeylisttoindicesalreadydicremoved

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
    aatomicnum=obaa.GetAtomicNum()
    datomicnum=obad.GetAtomicNum()
    babelindices=[a2,b2,c2,d2]
    torsionsmissing=databaseparser.ReadTorsionList(poltype,os.path.join(os.path.abspath('..'),poltype.torsionsmissingfilename))
    classes=[poltype.idxtosymclass[i] for i in babelindices]
    if classes not in torsionsmissing and classes[::-1] not in torsionsmissing: # then probably H torsions transferred
        return
    if (keyfilter is None or keyfilter == tpdkey): 
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

def insert_torprmdict(poltype,mol, torprmdict,max_amp):
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
    prmidx = -1
    # for each cls key
    for (chkclskey, torprm) in torprmdict.items():
        # for each parameter in the energy equation for this torsion
        # nfoldlist = [1,2,3]
        if poltype.fitfirsttorsionfoldphase==True:
            prmidx+=1
            torprm['firstphaseprmindex']=prmidx
            initialprms.append(0)
        for nfold in poltype.nfoldlist:
            # init array
            test_tor_energy = numpy.zeros(len(tmpx))
            # for each dihedral angle about this rotatable bond and the number of time it occurs 

            # normalize
            if torprm['phasedict']:
                prmidx += 1
                torprm['prmdict'][nfold] = prmidx
                if chkclskey in poltype.classkeytoinitialprmguess.keys():
                    prms=poltype.classkeytoinitialprmguess[chkclskey]
                    initialprms.append(prms[nfold-1])
                else: 
                    initialprms.append(max_amp)

        if not torprmdict[chkclskey]['prmdict']:
            torprmdict[chkclskey]['count'] = 0
            torprmdict[chkclskey]['phasedict'] = {}

        
    prmidx += 1
    initialprms.append(0)
    return prmidx,initialprms


def GenerateBoundaries(poltype,max_amp,refine,initialprms,torprmdict):
    lowerbounds=[]
    upperbounds=[]
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
            upperbounds.append(360)

        for prmidx in range(len(prms)):
            prm=prms[prmidx]
            value=numpy.abs(prm)
            scaledvalue=.3*value
            if refine==True:
                if allzero==False:
                    if value!=0:
                        lowerbounds.append(prm-scaledvalue)
                        upperbounds.append(prm+scaledvalue)
                    else:
                        lowerbounds.append(prm)
                        upperbounds.append(.2) # cant bound by 0,0
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
       
def GenerateNewParmGuessAndBounds(poltype,parm):
    bounds=[]
    if parm<0:
        newparmguess=-1.5
    elif parm>0:
        newparmguess=1.5
    
    adjustment=.001 
    bounds=[newparmguess-adjustment,newparmguess+adjustment]

    return newparmguess,bounds


def ModifyInitialGuessAndBoundries(poltype,parameterinfo,pzero,boundstup):
    boundstup=list(boundstup)
    lowerbounds=boundstup[0]
    upperbounds=boundstup[1]
    for ls in parameterinfo:
        newparmguess,bounds,index=ls[:]
        pzero[index]=newparmguess
        lowerbound=bounds[0]
        upperbound=bounds[1]
        lowerbounds[index]=lowerbound
        upperbounds[index]=upperbound
    
    boundstup=tuple([lowerbounds,upperbounds])

    return pzero,boundstup


def fit_rot_bond_tors(poltype,mol,cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,clskeyswithbadfits,indicesremoveddic,cls_angle_dict_unmodified,cls_mm_engy_dict_unmodified,cls_qm_engy_dict_unmodified):
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
    torsettobypassrmsd={}
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
        if len(poltype.onlyfittorstogether)!=0:
            torset=tuple(poltype.onlyfittorstogether)
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
        angle_list = cls_angle_dict[tup]  # Torsion angle for each corresponding energy
        mm_energy_list = cls_mm_engy_dict[tup]  # MM Energy before fitting to QM torsion energy
        qm_energy_list = cls_qm_engy_dict[tup]  # QM torsion energy

        mm_energy_list_unmodified = cls_mm_engy_dict_unmodified[tup]  # MM Energy before fitting to QM torsion energy
        qm_energy_list_unmodified = cls_qm_engy_dict_unmodified[tup]  # QM torsion energy
        # 'normalize'
        qm_energy_list = [en - min(qm_energy_list) for en in qm_energy_list]
        for e in qm_energy_list:
            if e>50:
                raise ValueError('Energy is greater than 50 kcal/mol for '+str(torset)+' '+str(qm_energy_list)+' '+str(angle_list))
        mm_energy_list = [en - min(mm_energy_list) for en in mm_energy_list]

        weightlist=numpy.exp(-numpy.array(qm_energy_list)/poltype.boltzmantemp)

        tor_energy_list = [qme - mme for qme,mme in zip(qm_energy_list,mm_energy_list)]
        qm_energy_list_unmodified = [en - min(qm_energy_list_unmodified) for en in qm_energy_list_unmodified]
        mm_energy_list_unmodified = [en - min(mm_energy_list_unmodified) for en in mm_energy_list_unmodified]
        tor_energy_list_unmodified = [qme - mme for qme,mme in zip(qm_energy_list_unmodified,mm_energy_list_unmodified)]
       
        max_amp = max(tor_energy_list) - min(tor_energy_list)
        amplist=[20,max_amp]
        max_amp=min(amplist)
        prmidx,initialprms = insert_torprmdict(poltype,mol, torprmdict,max_amp)
        pzero = initialprms
        boundstup=GenerateBoundaries(poltype,max_amp,refine,initialprms,torprmdict)
        # run leastsq until all the parameter estimates are reasonable
        parm_sanitized = False
        bypassrmsd=False
        maxiter=5
        count=0
        while not parm_sanitized:
            if count>=maxiter:
                break
            parm_sanitized = True
            keylist = list(torprmdict.keys())
            
            # creating a new function, errfunc
            # p: parameters
            # x: angle list
            # torprmdict: torsion information 
            # y: tor_energy_list
            if useweights==True: 
                errfunc = lambda p, x, z, torprmdict, y: weightlist*(fitfunc(poltype,p, x,z, torprmdict) - y)

            else:
                errfunc = lambda p, x, z, torprmdict, y: fitfunc(poltype,p, x,z, torprmdict) - y
            array=optimize.least_squares(errfunc, pzero,jac='2-point', bounds=boundstup,args=(torgen.rads(poltype,numpy.array(angle_list)),torset,torprmdict, tor_energy_list))
            p1=array['x']
            pzero,boundstup,parm_sanitized=CheckFitParameters(poltype,pzero,boundstup,parm_sanitized,refine,keylist,torprmdict,p1,angle_list,torset,max_amp)
            count+=1
 
        

        torprmdict,write_prm_dict,classkeytofoldtophase=FillInDictionariesParameterEstimates(poltype,torprmdict,p1,write_prm_dict)
        torsettobypassrmsd[torset]=bypassrmsd

        GeneratePlots(poltype,cls_angle_dict,torset,useweights,classkeylist,fitfunc_dict,torprmdict,mm_energy_list,tor_energy_list,flatphaselist,qm_energy_list,tup,indicesremoveddic,cls_angle_dict_unmodified,qm_energy_list_unmodified,mm_energy_list_unmodified,clskeyswithbadfits,weightlist)
    return write_prm_dict,fitfunc_dict,torsettobypassrmsd,classkeytofoldtophase



def CheckFitParameters(poltype,pzero,boundstup,parm_sanitized,refine,keylist,torprmdict,p1,angle_list,torset,max_amp):
    foldtoparmslist={}
    parameterinfo=[]
    parmtokey={}
    parmtoindex={}
    if refine==False:
        indicesmodified=[]
        for chkclskey in keylist:
            for nfold in torprmdict[chkclskey]['prmdict']:
                if nfold not in foldtoparmslist.keys():
                    foldtoparmslist[nfold]=[]
                index=torprmdict[chkclskey]['prmdict'][nfold]
                parm=p1[index]
                parmtokey[parm]=chkclskey
                parmtoindex[parm]=index
                foldtoparmslist[nfold].append(parm)
                if abs(parm) > max_amp:
                    parm_sanitized = False
                    indicesmodified.append(index)
                    newparmguess,bounds=GenerateNewParmGuessAndBounds(poltype,parm) 
                    parameterinfo.append([newparmguess,bounds,index])
                    break
        for fold,parmslist in foldtoparmslist.items():
            combs=list(combinations(parmslist,2))
            for comb in combs:
                keylist=[parmtokey[i] for i in comb]
                energyarray=fitfunc(poltype,p1,torgen.rads(poltype,numpy.array(angle_list)),torset,torprmdict,keyonlylist=keylist,nfoldonlylist=[fold])  
                minenergy=min(energyarray)
                energyarray=[i-minenergy for i in energyarray]
                energy=max(energyarray)
                maxV=max(comb)
                ratio=numpy.abs(energy/maxV)
                if ratio<.1 and numpy.abs(maxV)>=1.5: 
                    parm_sanitized = False 
                    parm=comb[0] # just fix one 
                    index=parmtoindex[parm]
                    indicesmodified.append(index)
                    newparmguess,bounds=GenerateNewParmGuessAndBounds(poltype,parm) 
                    parameterinfo.append([newparmguess,bounds,index])

                if fold==1 or fold==3:
                    allbig=True
                    for prm in comb:
                        if numpy.abs(prm)<15:
                            allbig=False
                    if allbig==True:
                        parm=comb[0] # just fix one 
                        index=parmtoindex[parm]
                        if index not in indicesmodified:
                            parm_sanitized = False
                            indicesmodified.append(index)
                            newparmguess,bounds=GenerateNewParmGuessAndBounds(poltype,parm) 
                            parameterinfo.append([newparmguess,bounds,index])
    

        pzero,boundstup=ModifyInitialGuessAndBoundries(poltype,parameterinfo,pzero,boundstup)
    return pzero,boundstup,parm_sanitized



def GeneratePlots(poltype,cls_angle_dict,torset,useweights,classkeylist,fitfunc_dict,torprmdict,mm_energy_list,tor_energy_list,flatphaselist,qm_energy_list,tup,indicesremoveddic,cls_angle_dict_unmodified,qm_energy_list_unmodified,mm_energy_list_unmodified,clskeyswithbadfits,weightlist):
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
    for clskey in classkeylist:
        fitfunc_dict[clskey] = fitfunc(poltype,'eval',torgen.rads(poltype,Sx),torset,torprmdict,debug=False)
        if clskey in clskeyswithbadfits:
            def RMSD(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(weightlist*(numpy.subtract(fitfunc_dict[clskey],tor_energy_list)),c))))


        else:
            def RMSD(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.subtract(fitfunc_dict[clskey],tor_energy_list),c))))
        result=fmin(RMSD,.5)
        minRMSD=RMSD(result[0])

    numprms=1 # offset parameter incldued with torsion force constant parameters
    for classkey in torprmdict:
        numprms+= len(torprmdict[classkey]['prmdict'].keys())
    string=' , '.join(list(torprmdict.keys()))
    datapts=len(mm_energy_list)
    poltype.WriteToLog('Torsions being fit '+string+' RMSD(QM-MM1)'+str(minRMSD))
    if dim==1:
        first=numpy.array([Sx[i][0] for i in range(len(Sx))])
         
        xpoints=numpy.array([Sx[i][0] for i in range(len(Sx))])
        x_new = numpy.linspace(xpoints.min(),xpoints.max(),500)
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        l1, = ax.plot(Sx,fitfunc_dict[clskey],'ro',color='red',label='Fit')
        fitarray=numpy.array(fitfunc_dict[clskey])
        torarray=numpy.array(tor_energy_list)
        f = interp1d(xpoints,fitarray, kind='quadratic')
        y_smooth=f(x_new)
        ax.plot(x_new,y_smooth,color='red')
        l2, = ax.plot(Sx,tor_energy_list,'bo',color='blue',label='QM-MM1')
        f = interp1d(xpoints,torarray, kind='quadratic')
        y_smooth=f(x_new)
        ax.plot(x_new,y_smooth,color='blue')

        plt.legend(handles=[l1,l2],loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
        ax.set_xlabel('Dihedral Angle (degrees)')
        ax.set_ylabel('Energy (kcal/mol)')

        ax.text(0.05, 1.1, 'Torsions Being Fit =%s'%(string), transform=ax.transAxes, fontsize=10,verticalalignment='top')
        ax.text(0, -0.1, 'FoldNum=%s NumPrms=%s DataPts=%s RMSD(fit,QM-MM1),Abs=%s'%(str(len(poltype.nfoldlist)),str(numprms),str(len(mm_energy_list)),round(minRMSD,2)), transform=ax.transAxes, fontsize=10,verticalalignment='bottom')
        fig.savefig(figfname)
    elif dim==2:
        tor_energy_list_unmodified = numpy.array([qme - mme for qme,mme in zip(qm_energy_list_unmodified,mm_energy_list_unmodified)])

        tormat,qmmat,mmmat,idealanglematrix,actualanglematrix=FillInEnergyTensors(poltype,flatphaselist,cls_angle_dict[tup],tor_energy_list,qm_energy_list,mm_energy_list,torset,indicesremoveddic,cls_angle_dict_unmodified[tup],tor_energy_list_unmodified,qm_energy_list_unmodified,mm_energy_list_unmodified)
        PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,tormat,'QM-MM1 Heatmap (kcal/mol)','QM-MM1_Heatmap.png',numprms,datapts,minRMSD,string)
        PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,qmmat,'QM Heatmap (kcal/mol)','QM_Heatmap.png',numprms,datapts,minRMSD,string)
        PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,mmmat,'MM1 Heatmap (kcal/mol)','MM1_Heatmap.png',numprms,datapts,minRMSD,string)



def FillInDictionariesParameterEstimates(poltype,torprmdict,p1,write_prm_dict):
    classkeytofoldtophase={}
    for chkclskey in torprmdict:
        classkeytofoldtophase[chkclskey]=poltype.foldoffsetlist.copy()
        if 'firstphaseprmindex' in torprmdict[chkclskey].keys():
            prmindex=torprmdict[chkclskey]['firstphaseprmindex']
            prm=p1[prmindex]
            phaselist=classkeytofoldtophase[chkclskey] 
            phaselist[0]=prm   
            classkeytofoldtophase[chkclskey]=phaselist
            torprmdict[chkclskey]['firstphaseprmindex']=prm

        for nfold in torprmdict[chkclskey]['prmdict']:
            parm  = p1[torprmdict[chkclskey]['prmdict'][nfold]]
            torprmdict[chkclskey]['prmdict'][nfold] = parm
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
    return torprmdict,write_prm_dict,classkeytofoldtophase



def PlotHeatmap(poltype,idealanglematrix,actualanglematrix,matrix,title,figfname,numprms,datapts,minRMSD,textstring=None):
    fig, ax = plt.subplots()
    im = ax.imshow(matrix)

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('', rotation=-90, va="bottom")
    yangles=idealanglematrix[:,0,0]
    yangles=[round(i) for i in yangles]
    yangles=[str(i) for i in yangles]
    yangles=[i.split('.')[0] for i in yangles]
    xangles=idealanglematrix[0,:,1]
    xangles=[round(i) for i in xangles]
    xangles=[str(i) for i in xangles]
    xangles=[i.split('.')[0] for i in xangles]

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
            energyvalue=str(round(matrix[i,j]))
            text = ax.text(j, i, energyvalue,ha="center", va="center", color="w")
    
    ax.set_title(title)
    if textstring!=None and numprms!=None and datapts!=None and minRMSD!=None:
        ax.text(0.05, 1.1, 'Torsions Being Fit =%s'%(textstring), transform=ax.transAxes, fontsize=10,verticalalignment='top')
        ax.text(0, -0.1, 'FoldNum=%s NumPrms=%s DataPts=%s RMSD(fit,QM-MM1),Abs=%s'%(str(len(poltype.nfoldlist)),str(numprms),str(datapts),round(minRMSD,2)), transform=ax.transAxes, fontsize=10,verticalalignment='bottom')
    elif textstring==None and numprms==None and datapts==None and minRMSD!=None:

        ax.text(0.05, 1.1, 'RMSD(MM2,QM)=%s'%(round(minRMSD,2)), transform=ax.transAxes, fontsize=12,verticalalignment='top')



    fig.tight_layout()
    plt.show() 
    fig.savefig(figfname)




def FillInEnergyTensors(poltype,phaseanglearray,actualanglearray,tor_energy_list,qm_energy_list,mm_energy_list,torset,indicesremoved,actualanglearrayunmodified,tor_energy_list_unmodified,qm_energy_list_unmodified,mm_energy_list_unmodified):

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
        haveenergy=True
        if i in indicesremoved:
            haveenergy=False
        phasetup=phaseanglearray[i]
        angletup=actualanglearrayunmodified[i]
        indexes=numpy.where(numpy.all(phasetensor == numpy.array(phasetup), axis=-1))
        if haveenergy==True:
            torenergy=tor_energy_list_unmodified[i]
            qmenergy=qm_energy_list_unmodified[i]
            mmenergy=mm_energy_list_unmodified[i]
        else:
            torenergy=0
            qmenergy=0
            mmenergy=0
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


def write_key_file(poltype,write_prm_dict,tmpkey1basename,tmpkey2basename,classkeytofoldtophase):
    """
    Intent: Output the new key file based on parameters in write_prm_dict
    """
    tmpfh1 = open(tmpkey1basename, "r")
    tmpfh2 = open(tmpkey2basename, "w")
    fitline="# Fitted torsion" +"\n"
    for line in tmpfh1:
        m = re.search(r'torsion',line)
        if m is None or '#' in line or 'Missing' in line:
            tmpfh2.write(line)
        else:
            linarr = line.split()
            cl = linarr[1:5]
            revcl=cl[::-1]
            revclskey= ' '.join(revcl) 
            clskey = ' '.join(cl) # Order is fine (read from *.prm file)
            torline = line
            torvals = [float(ele) for ele in linarr[5:24:3]]
            if clskey in classkeytofoldtophase.keys():
               phaselist=classkeytofoldtophase[clskey]
            elif revclskey in classkeytofoldtophase.keys():
                phaselist=classkeytofoldtophase[revclskey]
            else:
                phaselist=poltype.foldoffsetlist.copy()
            if clskey in write_prm_dict:
                torline = ' torsion %7s %4s %4s %4s   ' % (cl[0],cl[1],cl[2],cl[3])
                for (nfold, prm) in write_prm_dict[clskey].items():
                    torline += ' %7.3f %.1f %d' % (prm,phaselist[nfold - 1], nfold)
                torline += '\n'
                tmpfh2.write(fitline)

            elif revclskey in write_prm_dict:
                torline = ' torsion %7s %4s %4s %4s   ' % (cl[0],cl[1],cl[2],cl[3])
                for (nfold, prm) in write_prm_dict[revclskey].items():
                    torline += ' %7.3f %.1f %d' % (prm,phaselist[nfold - 1], nfold)
                torline += '\n'
                tmpfh2.write(fitline)

            tmpfh2.write(torline)
    tmpfh1.close()
    tmpfh2.close()

def ConvertNoneToZero(poltype,listinput):
    indices=[]
    for i in range(len(listinput)):
        value=listinput[i]
        if value==None:
            indices.append(i)
    for i in range(len(listinput)):
        if i in indices:
            listinput[i]=0
        
    return listinput

def AssignZeros(poltype,array):
    newarray=[]
    allindices=[]
    for ls in array:
        indices=[i for i in range(len(ls)) if ls[i]==0 or i==None]
        for index in indices:
            if index not in allindices:
                allindices.append(index)
    for ls in array:
        for index in allindices:
            ls[index]=0
        newarray.append(ls)

    return newarray


def eval_rot_bond_parms(poltype,mol,fitfunc_dict,tmpkey1basename,tmpkey2basename,count,clskeyswithbadfits,torsettobypassrmsd,classkeylisttoindicesalreadydicremoved,classkeylisttoindicesremoved,write_prm_dict):
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
    classkeytofitresults={}
    # for each main torsion
    for torset in poltype.torlist:
        if torset not in poltype.nonaroringtorsets and len(torset)==2 and poltype.torfit1Drotonly==True:
            continue
        flatphaselist=poltype.torsettophaselist[tuple(torset)]

        tmpkeyfname = 'tmp.key'
        qm_energy_list,qang_list,WBOarray,energytophaseangle = compute_qm_tor_energy(poltype,torset,mol,flatphaselist)
        mm_energy_list,mang_list,tor_e_list = compute_mm_tor_energy(poltype,mol,torset,'_preQMOPTprefit',flatphaselist,tmpkeyfname)
        prepostmm_energy_list,prepostmang_list,preposttor_e_list = compute_mm_tor_energy(poltype,mol,torset,'_preQMOPTpostfit',flatphaselist,tmpkeyfname)
        mm2_energy_list,m2ang_list,tor_e_list2 = compute_mm_tor_energy(poltype,mol,torset,'_postQMOPTpostfit',flatphaselist,tmpkey2basename)

        if len(poltype.onlyfittorstogether)!=0:
            torset=tuple(poltype.onlyfittorstogether)
        bypassrmsd=torsettobypassrmsd[torset]
        classkeylist=[]
        for tor in torset:
            a,b,c,d = tor[0:4]
            key=str(b)+' '+str(c)
            torang = mol.GetTorsion(a,b,c,d)
            atmnuma = mol.GetAtom(a).GetAtomicNum()
            atmnumd = mol.GetAtom(d).GetAtomicNum()
            # clskey
            clskey = torgen.get_class_key(poltype,a, b, c, d)
            classkeylist.append(clskey)
        
        # get the qm energy profile
        shutil.copy(tmpkey1basename, tmpkeyfname)
        # get the original mm energy profile
        tup=tuple(classkeylist)
        if tup not in classkeylisttoindicesremoved.keys():
            classkeylisttoindicesremoved[tup]=[]
        # get the new mm energy profile (uses new parameters to find energies)
        originalmang_list=mang_list.copy()
        originalqang_list=qang_list.copy()
        originalmm2ang_list=m2ang_list.copy()
        originalmm2_energy_list=mm2_energy_list.copy()
        originalmm_energy_list=mm_energy_list.copy()
        originalprepostmm_energy_list=prepostmm_energy_list.copy()
        originalqm_energy_list=qm_energy_list.copy()
        arrays=AssignZeros(poltype,[originalmm_energy_list,originalqm_energy_list,originalmm2_energy_list,originalprepostmm_energy_list])
        originalmm_energy_list=arrays[0]
        originalqm_energy_list=arrays[1]
        originalmm2_energy_list=arrays[2]
        originalprepostmm_energy_list=arrays[3]
        originaltor_e_list=tor_e_list.copy()
        originaltor_e_list2=tor_e_list2.copy()
        originalWBOarray=WBOarray.copy()
        originalm2:ng_list=m2ang_list.copy()
        # remove angles for which energy was unable to be found
        del_ang_list = find_del_list(poltype,originalmm_energy_list,originalmang_list)
        classkeylisttoindicesremoved[tup]=FindRemovedIndices(poltype,originalmang_list,del_ang_list,indicesremoved=classkeylisttoindicesremoved[tup])
        array=prune_mme_error(poltype,classkeylisttoindicesremoved[tup],originalmang_list,originalmm_energy_list,originalmm2ang_list,originalmm2_energy_list,originalqm_energy_list,originalqang_list,originaltor_e_list,originaltor_e_list2,originalWBOarray)

        (mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2,WBOarray)=array[:-1][0]

        x=array[-1]
        del_ang_list = find_del_list(poltype,originalqm_energy_list,originalqang_list)
        classkeylisttoindicesremoved[tup]=FindRemovedIndices(poltype,originalqang_list,del_ang_list,indicesremoved=classkeylisttoindicesremoved[tup])
        array=prune_qme_error(poltype,classkeylisttoindicesremoved[tup],originalmang_list,originalmm_energy_list,originalmm2ang_list,originalmm2_energy_list,originalqm_energy_list,originalqang_list,originaltor_e_list,originaltor_e_list2,originalWBOarray)
        (mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2,WBOarray)=array[:-1][0]

        x=array[-1]
        del_ang_list = find_del_list(poltype,originalmm2_energy_list,originalmm2ang_list)
        classkeylisttoindicesremoved[tup]=FindRemovedIndices(poltype,originalmm2ang_list,del_ang_list,indicesremoved=classkeylisttoindicesremoved[tup])
        array=prune_qme_error(poltype,classkeylisttoindicesremoved[tup],originalmang_list,originalmm_energy_list,originalmm2ang_list,originalmm2_energy_list,originalqm_energy_list,originalqang_list,originaltor_e_list,originaltor_e_list2,originalWBOarray)
         
        (mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2,WBOarray)=array[:-1][0]

        x=array[-1]
        del_ang_list = find_del_list(poltype,originalWBOarray,originalqang_list)
        classkeylisttoindicesremoved[tup]=FindRemovedIndices(poltype,originalqang_list,del_ang_list,indicesremoved=classkeylisttoindicesremoved[tup])
        array=prune_qme_error(poltype,classkeylisttoindicesremoved[tup],originalmang_list,originalmm_energy_list,originalmm2ang_list,originalmm2_energy_list,originalqm_energy_list,originalqang_list,originaltor_e_list,originaltor_e_list2,originalWBOarray)
        (mang_list,mm_energy_list,m2ang_list,mm2_energy_list,qm_energy_list,qang_list,tor_e_list,tor_e_list2,WBOarray)=array[:-1][0]
        y=array[-1]

        classkeylisttoindicesremoved[tup].sort(reverse=True)
        for del_ang_idx in classkeylisttoindicesremoved[tup]:
            if del_ang_idx not in classkeylisttoindicesalreadydicremoved[tup]: 
                transformedidx=TransformDeleteIdx(del_ang_idx,classkeylisttoindicesalreadydicremoved[tup])
                fitfunc_dict[clskey]=numpy.delete(fitfunc_dict[clskey],transformedidx)
        array=numpy.array(originalmm_energy_list)
        originalfitfuncarray=numpy.zeros(array.shape)
        for i in range(len(originalmm_energy_list)):
           e=originalmm_energy_list[i]
           if e==0:
               originalfitfuncarray[i]=0
           else:
               originalfitfuncarray[i]=1
        # normalize profiles
        qm_energy_list = [en - min(qm_energy_list) for en in qm_energy_list]
        mm_energy_list = [en - min(mm_energy_list) for en in mm_energy_list]
        mm2_energy_list = [en - min(mm2_energy_list) for en in mm2_energy_list]
        prepostmm_energy_list = [en - min(prepostmm_energy_list) for en in prepostmm_energy_list]


        originalmm2_energy_list=ConvertNoneToZero(poltype,originalmm2_energy_list)
        originalmm_energy_list=ConvertNoneToZero(poltype,originalmm_energy_list)
        originalprepostmm_energy_list=ConvertNoneToZero(poltype,originalprepostmm_energy_list)

        originalmm2_energy_list= [en - min(originalmm2_energy_list) for en in originalmm2_energy_list]
        originalmm_energy_list= [en - min(originalmm_energy_list) for en in originalmm_energy_list]
        originalprepostmm_energy_list= [en - min(originalprepostmm_energy_list) for en in originalprepostmm_energy_list]

        originalff_list = [aa+bb for (aa,bb) in zip(originalmm_energy_list,originalfitfuncarray)]
        # find the difference between the two energy due to torsion profiles 
        tordif_list = [e2-e1 for (e1,e2) in zip(tor_e_list,tor_e_list2)]
        # normalize
        tordif_list = [en - min(tordif_list) for en in tordif_list]
        # find the difference between the two mm energy profiles
        tordifmm_list = [e1+e2 for (e1,e2) in zip (tordif_list,mm_energy_list)]
        tordifmm_list = [en - min(tordifmm_list) for en in tordifmm_list]
        # TBC
        weights=numpy.exp(-numpy.array(qm_energy_list)/poltype.boltzmantemp)
        weights = [float(i)/sum(weights) for i in weights] 
        ff_list = [aa+bb for (aa,bb) in zip(mm_energy_list,fitfunc_dict[clskey])]
        shifted_mm2_energy_list=numpy.add(1,mm2_energy_list)
        shifted_qm_energy_list=numpy.add(1,qm_energy_list)
        
        final_tor_energy_list=numpy.subtract(mm2_energy_list,qm_energy_list)
        final_relative_tor_energy_list=numpy.subtract(shifted_mm2_energy_list,shifted_qm_energy_list)
        if clskey in clskeyswithbadfits:
            final_tor_energy_list=numpy.multiply(final_tor_energy_list,weights)
            final_relative_tor_energy_list=numpy.multiply(final_relative_tor_energy_list,weights)
        if len(ff_list)==len(mm2_energy_list):
            def RMSD(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(final_tor_energy_list,c))))
            result=fmin(RMSD,.5)
            minRMSD=RMSD(result[0])

            def RMSDRel(c):
                return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.divide(final_relative_tor_energy_list,shifted_qm_energy_list),c))))
            resultRel=fmin(RMSDRel,.5)
            minRMSDRel=RMSDRel(resultRel[0])
        # output the profiles as plots

        figfname = "%s-energy-" % (poltype.molecprefix)
        for tor in torset:
            a,b,c,d = tor[0:4]
            figfname+="%d-%d" % (b,c)
            figfname+='_'
        figfname=figfname[:-1]
        if clskey in clskeyswithbadfits:
            wstring='_use_weights'
            ostring='Boltzmann Fit'
        else:
            wstring=''
            ostring=''
        figfname+=wstring+'.png'
        dim=len(mang_list[0])
        datapts=len(mm_energy_list)
        numprms=None 
        clskeysplit=clskey.split()
        middle=[clskeysplit[1],clskeysplit[2]]
        for testkey in write_prm_dict.keys():
            testkeysplit=testkey.split()
            themiddle=[testkeysplit[1],testkeysplit[2]]
            if middle==themiddle or middle==themiddle[::-1]:
                thestring='Torsion '+str(testkey)+' RMSD(MM2,QM) '+str(minRMSD)+' '+'RelativeRMSD(MM2,QM) '+str(minRMSDRel)+' '+ostring+'\n'
                poltype.WriteToLog(thestring)
                classkeytofitresults[testkey]=thestring
        if dim==1: 
            numpy.savetxt('qmarray_'+str(b)+'-'+str(c)+'.txt',qm_energy_list)
            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111)
            # energy profiles: mm (pre-fit), mm (post-fit), qm
            line1, =ax.plot(mang_list,mm_energy_list,'go',color='green',label='MM1 (prefit)')
            xpoints=numpy.array([mang_list[i][0] for i in range(len(mang_list))])
            x_new = numpy.linspace(xpoints.min(),xpoints.max(),500)
            f = interp1d(xpoints,numpy.array(mm_energy_list), kind='quadratic')
            y_smooth=f(x_new)
            ax.plot(x_new,y_smooth,color='green')

            line2, =ax.plot(m2ang_list,mm2_energy_list,'ro',color='red',label='MM2 (postfit)')
            xpoints=numpy.array([m2ang_list[i][0] for i in range(len(m2ang_list))])
            x_new = numpy.linspace(xpoints.min(),xpoints.max(),500)
            f = interp1d(xpoints,numpy.array(mm2_energy_list), kind='quadratic')
            y_smooth=f(x_new)
            ax.plot(x_new,y_smooth,color='red')

            line3, =ax.plot(qang_list,qm_energy_list,'bo',color='blue',label='QM')
            xpoints=numpy.array([qang_list[i][0] for i in range(len(qang_list))])
            x_new = numpy.linspace(xpoints.min(),xpoints.max(),500)
            f = interp1d(xpoints,numpy.array(qm_energy_list), kind='quadratic')
            y_smooth=f(x_new)
            ax.plot(x_new,y_smooth,color='blue')

            ax.text(0.05, 1.1, 'RMSD(MM2,QM)=%s , RMSDRel(MM2,QM)=%s'%(round(minRMSD,2),round(minRMSDRel,2)), transform=ax.transAxes, fontsize=12,verticalalignment='top')

            # mm + fit
            line4, =ax.plot(mang_list,ff_list,'mo',color='magenta',label='MM1+Fit')
            xpoints=numpy.array([mang_list[i][0] for i in range(len(mang_list))])
            x_new = numpy.linspace(xpoints.min(), xpoints.max(),500)
            f = interp1d(xpoints,numpy.array(ff_list), kind='quadratic')
            y_smooth=f(x_new)
            ax.plot(x_new,y_smooth,color='magenta')

            # prefit tinker XYZ structure but with postfit parameters
            line5, =ax.plot(prepostmang_list,prepostmm_energy_list,'o',color='black',label='MM1XYZMM2Prm')
            xpoints=numpy.array([prepostmang_list[i][0] for i in range(len(mang_list))])
            x_new = numpy.linspace(xpoints.min(), xpoints.max(),500)
            f = interp1d(xpoints,numpy.array(prepostmm_energy_list), kind='quadratic')
            y_smooth=f(x_new)
            ax.plot(x_new,y_smooth,color='black')


            ax2=ax.twinx()
            # make a plot with different y-axis using second axis object
            try:
                line6, =ax2.plot(qang_list,WBOarray,'yo',color='yellow',label='WBO')
                xpoints=numpy.array([qang_list[i][0] for i in range(len(qang_list))])
                x_new = numpy.linspace(xpoints.min(), xpoints.max(),500)
                ypoints=numpy.array([WBOarray[i][0] for i in range(len(WBOarray))])
                f = interp1d(xpoints,ypoints, kind='quadratic')
                y_smooth=f(x_new)
                ax2.plot(x_new,y_smooth,color='yellow')

                ax2.set_ylabel("WBO",color="blue",fontsize=14)
            except:
                pass
            ax.set_xlabel('Dihedral Angle (degrees)')
            ax.set_ylabel('Energy (kcal/mol)')
            plt.legend(handles=[line1,line2,line3,line4,line5,line6],loc=9, bbox_to_anchor=(0.5, -0.1), ncol=5)

            fig = plt.gcf()
            plt.show()
            fig.savefig(figfname)
        elif dim==2:
            mmmat,mm2mat,fmat,idealanglematrix,actualanglematrix=FillInEnergyTensors(poltype,flatphaselist,mang_list,mm_energy_list,mm2_energy_list,ff_list,torset,classkeylisttoindicesremoved[tup],originalmang_list,originalmm_energy_list,originalmm2_energy_list,originalff_list)
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
        if float(minRMSD)>poltype.maxtorRMSPD and float(minRMSDRel)>poltype.maxtorRMSPDRel:
            poltype.WriteToLog('Absolute or Relative RMSPD of QM and MM torsion profiles is high, RMSPD = '+ str(minRMSD)+' Tolerance is '+str(poltype.maxtorRMSPD)+' kcal/mol '+'RMSPDRel ='+str(minRMSDRel)+' tolerance is '+str(poltype.maxtorRMSPDRel))
            if poltype.suppresstorfiterr==False and count>0 and bypassrmsd==False and poltype.tordebugmode==False and clskey in clskeyswithbadfits:
                raise ValueError('RMSPD of QM and MM torsion profile is high, tried fitting to minima and failed, RMSPD = '+str(minRMSD)+','+str(minRMSDRel))
            
            
            clskeyswithbadfits.append(clskey)
    return clskeyswithbadfits,classkeytofitresults
                 

def TransformDeleteIdx(del_ang_idx,indicesalreadyremoved):
    indicesleftofindex=0
    for index in indicesalreadyremoved:
        if index<del_ang_idx:
            indicesleftofindex+=1
    transformedidx=del_ang_idx-indicesleftofindex 

    return transformedidx
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
    # copy *.key_5 to the directory qm-torsion
    shutil.copy(poltype.key6fname, tmpkey1fname)
    # change directory to qm-torsion
    os.chdir(tordir)

    # Group all rotatable bonds with the same classes and identify
    # the torsion parameters that need to be fitted.

    # For each rotatable bond, get torsion energy profile from QM
    # and MM (with no rotatable bond torsion parameters)
    # Get QM and MM (pre-fit) energy profiles for torsion parameters
    if poltype.tortor==True: 
        DecomposeTorsionTorsion(poltype,mol)
    cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,indicesremoveddic, cls_angle_dict_unmodified,cls_mm_engy_dict_unmodified,cls_qm_engy_dict_unmodified,classkeylisttoindicesalreadydicremoved= get_qmmm_rot_bond_energy(poltype,mol,tmpkey1basename,'_preQMOPTprefit')
    # if the fit has not been done already
    clskeyswithbadfits=[]

    count=0
    while 1:
        # do the fit
        if count==3:
            break # dont redo fitting forever
        write_prm_dict,fitfunc_dict,torsettobypassrmsd,classkeytofoldtophase = fit_rot_bond_tors(poltype,mol,cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,clskeyswithbadfits,indicesremoveddic,cls_angle_dict_unmodified,cls_mm_engy_dict_unmodified,cls_qm_engy_dict_unmodified)
        # write out new keyfile
        write_key_file(poltype,write_prm_dict,tmpkey1basename,tmpkey2basename,classkeytofoldtophase)

        RemoveFiles(poltype,'post',2)  
        if len(poltype.torlist)!=0:
            PostfitMinAlz(poltype,tmpkey2basename,'')
        # evaluate the new parameters
        clskeyswithbadfits,classkeytofitresults=eval_rot_bond_parms(poltype,mol,fitfunc_dict,tmpkey1basename,tmpkey2basename,count,clskeyswithbadfits,torsettobypassrmsd,classkeylisttoindicesalreadydicremoved,indicesremoveddic,write_prm_dict)
        if len(clskeyswithbadfits)==0:
            break
        count+=1
    WriteOutFitResults(poltype,tmpkey2basename,classkeytofitresults)
    if poltype.torfit1Drotonly==True and poltype.tortor==True:
        poltype.torfit1Drotonly=False
        poltype.torfit2Drotonly=True 
        cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,indicesremoveddic,cls_angle_dict_unmodified,cls_mm_engy_dict_unmodified,cls_qm_engy_dict_unmodified,classkeylisttoindicesalreadydicremoved= get_qmmm_rot_bond_energy(poltype,mol,tmpkey1basename,'_postQMOPTpostfit')
        PrepareTorTorSplineInput(poltype,cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,mol,tmpkey2basename,indicesremoveddic,cls_angle_dict_unmodified,cls_mm_engy_dict_unmodified,cls_qm_engy_dict_unmodified)

    shutil.copy(tmpkey2basename,'../' + poltype.key7fname)
    os.chdir('..')


def WriteOutFitResults(poltype,tmpkey2basename,classkeytofitresults):
    temp=open(tmpkey2basename,'r')
    results=temp.readlines()
    temp.close()
    tempname=tmpkey2basename.replace('.key','_TEMP.key')
    temp=open(tempname,'w')
    for line in results:
        linesplit=line.split()
        if 'torsion' in line and '#' not in line:
            fwd='%d %d %d %d' % (int(linesplit[1]), int(linesplit[2]), int(linesplit[3]), int(linesplit[4]))       
            fwdsplit=fwd.split() 
            revsplit=fwdsplit[::-1]
            rev='%d %d %d %d' % (int(revsplit[0]), int(revsplit[1]), int(revsplit[2]), int(revsplit[3]))
            for classkey,fitresults in classkeytofitresults.items():
                if classkey==fwd or classkey==rev:
                    newresults='# '+fitresults
                    temp.write(newresults)     

        temp.write(line)
    temp.close()
    os.remove(tmpkey2basename)
    os.rename(tempname,tmpkey2basename)


def RemoveFiles(poltype,string,occurance):
    files=os.listdir()
    for f in files:
        count = f.count(string)  
        if count==occurance:
            os.remove(f)


def PostfitMinAlz(poltype,keybasename,keybasepath):
    for outputlog in poltype.optoutputtotorsioninfo.keys():
        term,error=poltype.CheckNormalTermination(outputlog)
        [torset,optmol,variabletorlist,phaseangles,bondtopology,optoutputlog,initialxyz]=poltype.optoutputtotorsioninfo[outputlog]
        if term==True:   
            if not poltype.use_gaus:
                cartxyzname=optoutputlog.replace('.log','.xyz')
                cartxyz,torxyzfname=torgen.tinker_minimize_analyze_QM_Struct(poltype,torset,optmol,variabletorlist,phaseangles,cartxyzname,poltype.torsionrestraint,'_postQMOPTpostfit',keybasename,keybasepath,bondtopology)
                cartxyz,torxyzfname=torgen.tinker_minimize_analyze_QM_Struct(poltype,torset,optmol,variabletorlist,phaseangles,cartxyzname,poltype.torsionrestraint,'_preQMOPTpostfit',keybasename,keybasepath,bondtopology,tinkerxyz=initialxyz,minimize=False)

            else:
                cartxyz,torxyzfname=torgen.tinker_minimize_analyze_QM_Struct(poltype,torset,optmol,variabletorlist,phaseangles,outputlog,poltype.torsionrestraint,'_postQMOPTpostfit',keybasename,keybasepath,bondtopology)
                cartxyz,torxyzfname=torgen.tinker_minimize_analyze_QM_Struct(poltype,torset,optmol,variabletorlist,phaseangles,outputlog,poltype.torsionrestraint,'_preQMOPTpostfit',keybasename,keybasepath,bondtopology,tinkerxyz=initialxyz,minimize=False)


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

def PrepareTorTorSplineInput(poltype,cls_mm_engy_dict,cls_qm_engy_dict,cls_angle_dict,mol,tmpkey2basename,indicesremoved,cls_angle_dict_unmodified,cls_mm_engy_dict_unmodified,cls_qm_engy_dict_unmodified):
    temp=open(tmpkey2basename,'a')
    for torset in poltype.torlist:
        if torset not in poltype.nonaroringtorsets and len(torset)==2:
            pass
        else:
            continue
        classkeylist=[]
        torsions=[] 
        for i in range(len(torset)):
            tor=torset[i]
            # get the atoms in the main torsion about this rotatable bond
            a,b,c,d = tor[0:4]
            # current torsion value
            torang = mol.GetTorsion(a,b,c,d)
            # class key; ie symmetry classes key
            clskey = torgen.get_class_key(poltype,a,b,c,d)
            classkeylist.append(clskey)
            torsions.append([a,b,c,d])
        tup=tuple(classkeylist)
        angle_list = cls_angle_dict[tup]  # Torsion angle for each corresponding energy
        mm_energy_list = cls_mm_engy_dict[tup]  # MM Energy before fitting to QM torsion energy
        qm_energy_list = cls_qm_engy_dict[tup]  # QM torsion energy
        qm_energy_list = [en - min(qm_energy_list) for en in qm_energy_list]
        mm_energy_list = [en - min(mm_energy_list) for en in mm_energy_list]
        tor_energy_list = numpy.array([qme - mme for qme,mme in zip(qm_energy_list,mm_energy_list)])
        mm_energy_list_unmodified = cls_mm_engy_dict_unmodified[tup]  # MM Energy before fitting to QM torsion energy
        qm_energy_list_unmodified = cls_qm_engy_dict_unmodified[tup]  # QM torsion energy
        qm_energy_list_unmodified = [en - min(qm_energy_list_unmodified) for en in qm_energy_list_unmodified]
        mm_energy_list_unmodified = [en - min(mm_energy_list_unmodified) for en in mm_energy_list_unmodified]
        tor_energy_list_unmodified = numpy.array([qme - mme for qme,mme in zip(qm_energy_list_unmodified,mm_energy_list_unmodified)])

        flatphaselist=poltype.torsettophaselist[tuple(torset)]
        
        tormat,qmmat,mmmat,idealanglematrix,actualanglematrix=FillInEnergyTensors(poltype,flatphaselist,cls_angle_dict[tup],tor_energy_list,qm_energy_list,mm_energy_list,torset,indicesremoved,cls_angle_dict_unmodified[tup],tor_energy_list_unmodified,qm_energy_list_unmodified,mm_energy_list_unmodified)
        firsttor=torsions[0]
        secondtor=torsions[1]
        tortorclskey,atomidxs=GenerateTorTorClasskey(poltype,firsttor,secondtor,poltype.idxtosymclass,poltype.rdkitmol)
        firstanglerow=idealanglematrix[0,:]
        firstrow=tormat[0,:]
        firstqmrow=qmmat[0,:]
        firstmmrow=mmmat[0,:]
        N=list(idealanglematrix.shape)[0]
        b = numpy.zeros((N+1,N,2))
        b[:-1,:,:] = idealanglematrix
        b[-1,:,:]=firstanglerow    
        idealanglematrix=numpy.copy(b)

        b = numpy.zeros((N+1,N))
        b[:-1,:] = tormat
        b[-1,:]=firstrow
        tormat=numpy.copy(b)


        b = numpy.zeros((N+1,N))
        b[:-1,:] = qmmat
        b[-1,:]=firstqmrow
        qmmat=numpy.copy(b)

        b = numpy.zeros((N+1,N))
        b[:-1,:] = mmmat
        b[-1,:]=firstmmrow
        mmmat=numpy.copy(b)


        firstanglecol=idealanglematrix[:,0]
        firstcol=tormat[:,0]
        firstqmcol=qmmat[:,0]
        firstmmcol=mmmat[:,0]

        b = numpy.zeros((N+1,N+1,2))
        b[:,:-1,:] = idealanglematrix
        b[:,-1,:]=firstanglecol        
        idealanglematrix=numpy.copy(b)

        b = numpy.zeros((N+1,N+1))
        b[:,:-1] = tormat
        b[:,-1]=firstcol
        tormat=numpy.copy(b)


        b = numpy.zeros((N+1,N+1))
        b[:,:-1] = mmmat
        b[:,-1]=firstmmcol
        mmmat=numpy.copy(b)

        b = numpy.zeros((N+1,N+1))
        b[:,:-1] = qmmat
        b[:,-1]=firstqmcol
        qmmat=numpy.copy(b)




        eps=.00001
        tormat[numpy.abs(tormat) < eps] = 0
        qmmat[numpy.abs(qmmat) < eps] = 0
        mmmat[numpy.abs(mmmat) < eps] = 0

        eps=10**7
        tormat[numpy.abs(tormat) > eps] = 0
        qmmat[numpy.abs(qmmat) > eps] = 0
        mmmat[numpy.abs(mmmat) > eps] = 0

        tormat=numpy.nan_to_num(tormat)
        tormat=numpy.around(tormat, decimals=3)
        qmmat=numpy.nan_to_num(qmmat)
        qmmat=numpy.around(qmmat, decimals=3)
        mmmat=numpy.nan_to_num(mmmat)
        mmmat=numpy.around(mmmat, decimals=3)

        idealanglematrix=numpy.around(idealanglematrix, decimals=3)
        actualanglematrix=numpy.around(actualanglematrix, decimals=3)

        PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,tormat,'QM-MM2 Heatmap (kcal/mol)','QM-MM2_Heatmap.png',None,None,None)
        PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,qmmat,'QM Heatmap (kcal/mol)','QM_Heatmap.png',None,None,None)
        PlotHeatmap(poltype,idealanglematrix,actualanglematrix ,mmmat,'MM2 Heatmap (kcal/mol)','MM2_Heatmap.png',None,None,None)

        rowpts=N+1
        colpts=N+1
        tortorline='tortors '+tortorclskey+' '+str(rowpts)+' '+str(colpts) +'\n'
        temp.write(tortorline)
        for i in range(len(tormat)):
            erow=tormat[i] 
            anglerow=idealanglematrix[i]
            for j in range(len(erow)):
                ecol=erow[j]
                anglecol=anglerow[j]
                line=str(anglecol[0])+' '+str(anglecol[1])+' '+str(ecol)+'\n'
                temp.write(line)
    temp.flush()
    temp.close()

def GenerateTorTorClasskey(poltype,firsttor,secondtor,idxtosymclass,mol):
    rdkitfirsttor=[i-1 for i in firsttor]
    rdkitsecondtor=[i-1 for i in secondtor]
    paths=Chem.rdmolops.FindAllPathsOfLengthN(mol,5,useBonds=False,useHs=True)
    for path in paths:
        allin=True
        for index in path:
            if index not in rdkitfirsttor and index not in rdkitsecondtor:
                allin=False
        if allin==True:
            break
    path=[i+1 for i in path]
    tortoratomidxs=[path[0],path[1],path[2],path[3],path[4]]
    tortortypeidxs=[idxtosymclass[j] for j in tortoratomidxs]
    tortortypeidxs=[str(j) for j in tortortypeidxs]
    tortorclskey=' '.join(tortortypeidxs)
    return tortorclskey,tortoratomidxs


