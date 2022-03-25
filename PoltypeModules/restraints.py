import keyfilemodifications as keymods
import math
import numpy as np
import mdtraj as md
import sys
import os
import time
import scipy.optimize
from scipy.spatial import distance_matrix
from scipy.cluster.hierarchy import fcluster, leaders, linkage
from scipy.spatial.distance import pdist
from collections import Counter, defaultdict
import openbabel


def ComputeCOM(poltype,atomidxtovecdic,atomidxtomassdic):
    num=np.array([0,0,0],dtype='float64')
    temp=np.array([0,0,0],dtype='float64')
    totmass=0
    for atomidx in atomidxtovecdic.keys():
        mass=atomidxtomassdic[atomidx]
        vec=atomidxtovecdic[atomidx]*10
        num+=mass*np.array(vec)
        temp+=np.array(vec)
        totmass+=mass
    COM=num/totmass
    return COM



def AddHarmonicRestrainGroupTermsToKeyFile(poltype,keyfilename,teatherdist):
    group1string='group '+str(1)+' '
    group2string='group '+str(2)+' '
    for num in poltype.restrainatomgroup1:
        group1string+=str(num)+' '
    keymods.AddKeyWord(poltype,keyfilename,group1string+'\n')       
    for num in poltype.restrainatomgroup2:
        group2string+=str(num)+' '
    keymods.AddKeyWord(poltype,keyfilename,group2string+'\n')

    if poltype.flatbotrest==False:
        restrainstring='restrain-groups '+str(1)+' '+ str(2)+' '+str(poltype.distancerestraintconstant)+' '+str(teatherdist)+' '+str(teatherdist)
    else:
        restrainstring='restrain-groups '+str(1)+' '+ str(2)+' '+str(poltype.distancerestraintconstant)+' '+'0'+' '+str(teatherdist)
    keymods.AddKeyWord(poltype,keyfilename,restrainstring+'\n')
        
    
def AverageCOMGroups(poltype,filename):
    poltype.WriteToLog('Averaging COM groups from '+filename,prin=True)
    try:
        t = md.load_arc(filename)
    except IOError:
        t = md.load(filename)

    g1COMtog2COMdistarray=[]
    g1atomidxtovecdic={}
    g1atomidxtomassdic={} 
    g2atomidxtovecdic={}
    g2atomidxtomassdic={}
    for framecount in range(t.n_frames):
        for a in t.topology.atoms:
            zeroindex=a.index    
            atomidx=zeroindex+1
            atomsymb=a.element.symbol
            mass=poltype.elementsymtomass[atomsymb]
            vec=t.xyz[framecount,zeroindex,:]
            if atomidx in poltype.restrainatomgroup1:
                g1atomidxtovecdic[atomidx]=vec
                g1atomidxtomassdic[atomidx]=mass
            elif atomidx in poltype.restrainatomgroup2:
                g2atomidxtovecdic[atomidx]=vec
                g2atomidxtomassdic[atomidx]=mass

        if len(g1atomidxtovecdic.keys())==len(poltype.restrainatomgroup1) and len(g2atomidxtovecdic.keys())==len(poltype.restrainatomgroup2):
            g1COM=ComputeCOM(poltype,g1atomidxtovecdic,g1atomidxtomassdic)
            g2COM=ComputeCOM(poltype,g2atomidxtovecdic,g2atomidxtomassdic)
            dispx=g1COM[0]-g2COM[0]
            dispy=g1COM[1]-g2COM[1]
            dispz=g1COM[2]-g2COM[2]
            dist=np.sqrt(dispx**2+dispy**2+dispz**2)
            g1COMtog2COMdistarray.append(dist)
            atomidxtovecdic={}
            atomidxtomassdic={}
    COM=round(np.mean(g1COMtog2COMdistarray),2)
    if COM<1:
        x=2
    else:
        x=COM+1
    return x

def AddInitialGroupRestraints(poltype):
    if poltype.restrainatomgroup1!=None and poltype.restrainatomgroup2!=None:
        teatherdist=AverageCOMGroups(poltype,poltype.xyzfilename)
        poltype.restraintdistance=teatherdist
        AddHarmonicRestrainGroupTermsToKeyFile(poltype,poltype.configkeyfilename,teatherdist)


def GroupRestraintFreeEnergyFix(poltype):
    fo=poltype.distancerestraintconstant
    pi=np.pi
    gasconst=8.314/4184
    temp=poltype.equilibriatescheme[-1]
    kt = temp * gasconst
    stdcon = (1.0*10**27) / (6.02214 * 10**23)
    if poltype.flatbotrest==True:
        ri=0
        fi=1
    else:
        ri=poltype.restraintdistance
        fi=poltype.distancerestraintconstant

    ro=poltype.restraintdistance
    v1 = 2.0*pi*ri*(-2.0+np.exp(-ri**2*fi/kt))*kt/fi + np.sqrt(kt*(pi/fi)**3)*(2.0*fi*ri*ri+kt)*math.erf(ri*np.sqrt(fi/kt))
    v2 = (4.0*pi/3.0) * (ro**3-ri**3)
    v3 = np.sqrt(kt*(pi/fo)**3) * (2.0*fo*ro*ro+kt+4.0*ro*np.sqrt(kt*fo/pi))
    vol = v1 + v2 + v3
    dv1 = 2.0*pi*ri**3*np.exp(-ri**2*fi/kt)/temp + 2.0*pi*ri*(-2.0+np.exp(-ri**2*fi/kt))*kt/(fi*temp) + 0.5*np.sqrt((pi/fi)**3)*np.sqrt(kt)*(2.0*ri**2*fi+kt)*math.erf(ri*np.sqrt(fi/kt))/temp - pi*ri*np.exp(-ri**2*fi/kt)*(2.0*ri**2*fi+kt)/(fi*temp) + np.sqrt((kt*pi/fi)**3)*math.erf(ri*np.sqrt(fi/kt))/temp
    dv2 = 0.0
    dv3 = np.sqrt(kt*(pi/fo)**3)*fo*ro*ro/temp + 4.0*kt*(pi/fo)*ro/temp + 1.5*np.sqrt((kt*pi/fo)**3)/temp
    dvol = dv1 + dv2 + dv3
    dg = -kt * np.log(vol/stdcon)
    ds = -dg/temp + kt*dvol/vol
    dh = dg + temp*ds
    poltype.rescorrection=dg
    return dg,dh,ds


def ComputeIdealRestraints(poltype,fxyz):
    ligidx = [_-1 for _ in poltype.ligandindices]
    rest_rotbond = False
    no_orient = False
    find_rotation_rest(poltype,fxyz, ligidx, "CA C N", rest_rotbond=rest_rotbond, rotbond_only=no_orient)


def get_adjlist(top):
    adjlist = defaultdict(list)
    for bond in top.bonds:
        a1, a2 = bond
        if not (a1.name.startswith('H') or a2.name.startswith('H') ):
            adjlist[a1.index].append(a2.index)
            adjlist[a2.index].append(a1.index)
    return adjlist

def get_tab_branch(top):
    '''
    Given mdtraj Trajectory.Topology, 
    returns a dictionary for the number of connected heavy atoms
    '''
    i = 0
    nbr = np.zeros(top.n_atoms, dtype=np.int)
    for bond in top.bonds:
        i += 1
        a1, a2 = bond
        if not (a1.name.startswith('H') or a2.name.startswith('H') ):
            nbr[a1.index] += 1
            nbr[a2.index] += 1
    return nbr

def write_tinker_idx(idxs):
    idx_out = []
    rs = []
    for i0 in sorted(idxs):
        if len(rs) == 0 or rs[-1][1]+1 < i0:
            rs.append([i0, i0])
        else: 
            rs[-1][1] = i0
    for r in rs:
        if r[0] == r[1]:
            idx_out.append(r[0])
        elif r[0] == r[1] - 1:
            idx_out.append(r[0])
            idx_out.append(r[1])
        else:
            idx_out.append(-r[0])
            idx_out.append(r[1])
    return idx_out
        
def target_disp(weights, coord, alpha=0.1):
    '''target function for center of mass displacement and number of atoms in the group
    '''
    assert len(weights) == coord.shape[0]
    wts = np.array(weights).reshape((-1, 1))
    wts = np.maximum(wts, 0)
    assert sum(wts) > 0
    wts *= 1.0/np.mean(wts)

    loss = np.sum(np.abs(np.sum(wts*coord, axis=0)))
    loss += alpha*np.sum(np.abs(np.abs(wts*coord)))

    for m in np.arange(2, 20):
        mask0 = (wts > m)
        loss += alpha*np.sum((np.abs((wts*mask0))))

    return loss

def calc_coord(traj, idxs):
    vals = np.zeros(traj.n_frames)
    idxs = np.array(idxs)
    NM_IN_ANG = 10.0
    RAD_IN_DEG = 180.0/np.pi
    if len(idxs) == 2:
        vals = md.compute_distances(traj, idxs.reshape(1, -1))
        vals *= NM_IN_ANG
    elif len(idxs) == 3:
        vals = md.compute_angles(traj, idxs.reshape(1, -1))
        vals *= RAD_IN_DEG
    elif len(idxs) == 4:
        vals = md.compute_dihedrals(traj, idxs.reshape(1, -1))
        vals *= RAD_IN_DEG
    return vals[:, 0]

def calc_idx_ortho(traj, i1, i2, idx3, nbr=None, method='long'):
    '''
    find the index that gives the largest ortho vector

    nbr: list of nr of branched atoms
    '''
    DELTA_R = 0.2
    DELTA_COS = 0.3

    if nbr is None:
        nbr = defaultdict(int)
    t = traj
    a1 = t.xyz[0, [i1], :] - t.xyz[0, [i2], :]
    u1 = a1/np.linalg.norm(a1, axis=1)

    vec2 = t.xyz[0, idx3, :] - t.xyz[0, [i2], :]
    d2 = np.linalg.norm(vec2, axis=1)
    d2p = np.abs(np.sum(vec2 * u1,axis=1))
    d2o = np.sqrt(d2**2.0 - d2p**2.0)
    d2cos = d2p / np.maximum(1e-5, d2)
    if method == 'short':
        # large angle (~90 deg), short distance
        d2on = np.array(list(zip([nbr[_]<1 for _ in idx3], -d2cos//DELTA_COS, -d2o)) , dtype=[('n', np.int),('c', np.float), ('r', np.float)])
    else:
        # large ortho vector
        d2on = np.array(list(zip([nbr[_] for _ in idx3], d2o)), dtype=[('n', np.int), ('r', np.float)])
    i3 = idx3[np.argsort(d2on)[-1]]
    return i3

def write_xyz_from_md(traj, idx):
    outp = '%d\nExtracted from MD\n'%(len(idx))
    for n, i in enumerate(idx):
        atom = traj.topology.atom(i)
        xyz = traj.xyz[0, i, :]*10
        outp += '%5s %12.6f %12.6f %12.6f\n'%(atom.element.symbol, xyz[0], xyz[1], xyz[2])
    return outp


def get_idx_pocket(traj, idx1, idx2, rcutoff=0.5):
    R_CLUSTER = rcutoff
    t = traj
    nbr = get_tab_branch(t.topology)
    sorted_i1, sorted_i2, dist1, dist2 = calc_idx_interface(traj, idx1, idx2, rcutoff=rcutoff)

    if len(sorted_i1) == 0:
        error_exit('No interface atom found within %.3f nm'%rcutoff)

    prot_idx = sorted_i1
    lig_idx = sorted_i2
    prot_dist = pdist(t.xyz[0, prot_idx, :]) 
    prot_distmat = distance_matrix(t.xyz[0, prot_idx, :],t.xyz[0, prot_idx, :]) 
    Z = linkage(prot_dist, method='average')

    clst_idx = fcluster(Z, R_CLUSTER, criterion='distance')
    clst_size = Counter(clst_idx)
    clst_size_list = sorted(clst_size.items(), key = lambda t:t[1])
    iclstm = clst_size_list[-1][0]

    flag_max = (clst_idx == iclstm)
    iprot = prot_idx[flag_max][0]

    lig_dist = np.linalg.norm(t.xyz[0, lig_idx, :] - t.xyz[0, [iprot], :], axis=1)
    dist_nr2 = np.array(list(zip(lig_dist, [nbr[_]<=1 for _ in lig_idx])), dtype=[('r', np.float), ('n', np.int)])
    ord_nr2 = np.argsort(dist_nr2, order=('n', 'r'))
    ilig = sorted_i2[ord_nr2[0]]
    iprot2 = calc_idx_ortho(traj, ilig, iprot, prot_idx, nbr)
    iprot3 = calc_idx_ortho(traj, iprot, iprot2, prot_idx[(prot_idx != iprot)*(prot_idx != iprot2)], nbr)


    ilig2 = calc_idx_ortho(traj, iprot, ilig, lig_idx[lig_idx != ilig], nbr, method='short')
    ilig3 = calc_idx_ortho(traj, ilig, ilig2, lig_idx[(lig_idx != ilig)*(lig_idx != ilig2)], nbr, method='short')

    int_idxs = [[iprot, ilig]]
    int_idxs.extend([[iprot2, iprot, ilig], [iprot, ilig, ilig2]])
    int_idxs.extend([[iprot3, iprot2, iprot, ilig], [iprot2, iprot, ilig, ilig2], [iprot, ilig, ilig2, ilig3]])

    return int_idxs

def calc_idx_interface(traj, idx1, idx2, rcutoff=0.6):
    t = traj
    distmat = distance_matrix(t.xyz[0, idx1, :], t.xyz[0, idx2, :])
    dist2 = np.min(distmat, axis=0)
    dist1 = np.min(distmat, axis=1)

    ord1 = np.argsort(dist1)
    ord2 = np.argsort(dist2)

    n1 = np.sum(dist1 <= rcutoff)
    n2 = np.sum(dist2 <= rcutoff)

    return idx1[ord1[:n1]], idx2[ord2[:n2]], dist1[ord1[:n1]], dist2[ord2[:n2]]

def write_rest(int_idxs, r0s, k0s, fmt='tinker'):
    rest_name = ['', '', 'restrain-distance', 'restrain-angle', 'restrain-torsion']
    outp =[]
    for ridx, r0, k0 in zip(int_idxs, r0s, k0s):
        nat = len(ridx)
        if nat >= len(rest_name):
            continue
        if fmt == 'tinker':
            outp.append( '%s %s %.6f %.6f %.6f\n'%(rest_name[nat], ' '.join('%4d'%(_+1) for _ in ridx), k0, r0, r0))
        else:
            print("Format %s not supported"%fmt)
            return ''
    return outp



def find_rotation_rest(poltype,fxyz, ligidx0, atomnames='CA', rcutoff=0.6, alpha=1.0, rest_rotbond=False, rotbond_only=False):
    try:
        t = md.load_arc(fxyz)
    except IOError:
        t = md.load(fxyz)
    if rotbond_only:
        return


    nbr = get_tab_branch(t.topology)

    ligidx = np.array(sorted(list(set(ligidx0) - set(t.topology.select('name H')))))
    if len(ligidx) == 0:
        error_exit("No ligand heavy atoms found")

    protidx0 = t.topology.select('name %s'%(atomnames))
    protidx0 = np.array(sorted(set(protidx0) - set(ligidx)))

    if len(protidx0) == 0:
        error_exit("No ligand atoms within %.3f nm of protein %s atoms"%(rcutoff, atomnames))

    poltype.int_idxs= get_idx_pocket(t, protidx0, ligidx)
    r0s = [calc_coord(t, _)[-1] for _ in poltype.int_idxs] # last frame
    RAD_IN_DEG = 180/np.pi
    k0s = np.zeros_like(r0s) + poltype.anglerestraintconstant
    k0s[0] = poltype.distancerestraintconstant
    restraintstrings=write_rest(poltype.int_idxs, r0s, k0s)
    keymods.RemoveKeyWord(poltype,poltype.configkeyfilename,'restrain')
    for string in restraintstrings:
        keymods.AddKeyWord(poltype,poltype.configkeyfilename,string+'\n')       

    RT = 8.314 * 298 / 4184
    # https://doi.org/10.1021/jp0217839
    # Eq. (14)
    assert len(k0s) == 6
    dgrest = -RT*np.log(1662*8*np.pi**2.0*np.sqrt(np.prod(k0s)*RAD_IN_DEG**10.0)/(r0s[0]**2.0*np.sin(r0s[1]/RAD_IN_DEG)*np.sin(r0s[2]/RAD_IN_DEG)*(2*np.pi*RT)**3))
    poltype.rescorrection=dgrest


def find_grp_idx(poltype,fxyz, ligidx0, atomnames='CA', rcutoff=1.2, alpha=1.0, extend=None):
    try:
        t = md.load_arc(fxyz)
    except IOError:
        t = md.load(fxyz)
    ligidx = np.array(sorted(list(set(ligidx0) - set(t.topology.select('name H')))))
    if len(ligidx) == 0:
        error_exit("No ligand heavy atoms found")

    protidx0 = t.topology.select('name %s'%(atomnames))
    idx0 = list(protidx0)
    if extend is not None:
        for delta_i in extend:
            idx0 += [_+delta_i for _ in protidx0]
    #protidx0 = np.array(sorted(set(protidx0) - set(ligidx)))
    protidx0 = np.array(sorted(set(idx0) - set(ligidx)))

    distmat = distance_matrix(t.xyz[0, ligidx, :], t.xyz[0, protidx0, :])
    distpro = np.min(distmat, axis=0)
    protidx = protidx0[distpro <= rcutoff]
    if len(protidx) == 0:
        return

    #print(distmat.shape)
    #print(np.min(distmat))
    imindist = np.argmin(distmat)
    iminlig = ligidx[imindist // len(protidx0)]

    com_lig = np.mean(t.xyz[0, list(ligidx), :], axis=0)
    com_lig = t.xyz[0, [iminlig], :]


    wts0 = np.ones(len(protidx))
    xyz_prot = t.xyz[0, protidx, :]
    xyz1_prot = xyz_prot - com_lig.reshape((1, -1))
    res = scipy.optimize.minimize(target_disp, wts0, args=(xyz1_prot, alpha))
    wts1 = np.array(res.x)
    wts1 = np.maximum(0, wts1)
    wts1 *= 1.0/np.sum(wts1)
    wtm = np.mean(wts1[wts1 > 0])
    mask1 = wts1 > 0.4*wtm
    #print(wts1)
    #print(np.mean(wts1))
    #print(mask1)
    com_p1 = (np.mean(t.xyz[0, protidx[mask1], :], axis=0)).reshape((1, -1))

    xyz_lig = t.xyz[0, ligidx, :]
    dists = np.linalg.norm(xyz_lig - com_p1, axis=1)
    #print('DIST', dists)
    imin = np.argmin(dists)

    dcom = np.linalg.norm(com_lig - com_p1)
    dmin = np.linalg.norm(xyz_lig[imin] - com_p1)


    idx_tinker = write_tinker_idx([_+1 for _ in protidx[mask1]])
    #print('#', (' '.join(['%5d'%(_) for _ in protidx[mask1]])))
    outp = ''
    sgrp = ''
    for n in idx_tinker:
        sgrp += ' %5d'%n
        if len(sgrp) > 50 and n > 0:
            #print('group 1 %s'%sgrp)
            outp += ('group 1 %s\n'%sgrp)
            sgrp = ''
    array=[ligidx[imin] + 1]
    newarray=ConvertRangeSyntaxToNumberList(idx_tinker) 
    newarray=list(set(newarray))
    poltype.restrainatomgroup1=newarray
    poltype.restrainatomgroup2=array
    outp += ('group 2 %s\n'%(' '.join(['%5d'%(_) for _ in [ligidx[imin] + 1]])))
    outp += ("#r_0=%.3f"%(dmin*10))
    #print('group 2 %s'%(' '.join(['%5d'%(_) for _ in [ligidx[imin] + 1]])))
    #print("#r_0=%.3f"%(dmin*10))
    return dmin*10, outp



def ConvertRangeSyntaxToNumberList(indices):
    newindices=[]
    for i in range(len(indices)):
        index=indices[i]
        if index>=0:
            newindices.append(index)
        else:
            nextindex=indices[i+1]
            rangedlist=list(range(int(np.abs(index)),nextindex+1)) 
            for idx in rangedlist:
                newindices.append(idx)
    return newindices



def error_exit(msg):
    print("ERROR:", msg)
    sys.exit(1)

def write_tinker_idx(idxs):
    idx_out = []
    rs = []
    for i0 in sorted(idxs):
        if len(rs) == 0 or rs[-1][1]+1 < i0:
            rs.append([i0, i0])
        else: 
            rs[-1][1] = i0
    for r in rs:
        if r[0] == r[1]:
            idx_out.append(r[0])
        elif r[0] == r[1] - 1:
            idx_out.append(r[0])
            idx_out.append(r[1])
        else:
            idx_out.append(-r[0])
            idx_out.append(r[1])
    return idx_out
        
def ComputeIdealGroupRestraints(poltype,fxyz):
    ligidx = [_-1 for _ in poltype.ligandindices[0]]
    if len(ligidx) == 0:
        error_exit("ligand keyword not found")
    res = None
    r_thr = 0.5
    for r0 in (0.5, 0.6, 0.7):
        for a in (1.0, 0.5, 0.2):
            rmin, outp = find_grp_idx(poltype,fxyz, ligidx, 'CA C', r0, alpha=a)
            if res is None or rmin < res[0]:
                res = rmin, outp
            if rmin < r_thr:
                break
        if rmin < r_thr:
            #keymods.AddKeyWord(poltype,poltype.configkeyfilename[0],outp)       
            break
    r0_rest = 2.0
    if res[0] >= r0_rest:
        raise ValueError("#Warning: r0 = %.3f > %.3f Angstrom"%(res[0], r0_rest)+" This means the centers of mass of the two groups are too far away,"+" which will lead to poor convergence. Please adjust group definitions.")
    if poltype.flatbotrest==False:
        string="restrain-groups 1 2 %f %.3f %.3f"%(poltype.distancerestraintconstant,min(r0_rest, res[0]+0.5),min(r0_rest, res[0]+0.5))+'\n'

    else:
        value=min(r0_rest, res[0]+0.5)
        string="restrain-groups 1 2 %f 0.0 %.3f"%(poltype.distancerestraintconstant,min(r0_rest, res[0]+0.5))+'\n'
    #keymods.AddKeyWord(poltype,poltype.configkeyfilename[0],string)       

    


