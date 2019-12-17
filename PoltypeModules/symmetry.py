import sys
import os
import openbabel


def CalculateSymmetry(poltype,pmol, frag_atoms, symmetry_classes):
    """
    Intent: Uses and builds on openbabel's 'GetGIVector' method which,
    "Calculates a set of graph invariant indexes using the graph theoretical distance,
    number of connected heavy atoms, aromatic boolean,
    ring boolean, atomic number, and summation of bond orders connected to the atom", to
    find initial symmetry classes.
    Input: 
        pmol: OBMol object 
        frag_atoms: OBBitVec object containing information about the largest fragment in 'pmol' 
        symmetry_classes: the symmetry_classes array which will be filled in
    Output: 
		nclasses: # of classes
        symmetry_classes: array is filled in
    Referenced By: gen_canonicallabels
    Description:
    1. vectorUnsignedInt object is created, 'vgi'
    2. It is filled in with the call to GetGIVector
    3. The 'symmetry_classes' array is initially filled in to match with 'vgi'
    4. These initial invariant classes do not suit our needs perfectly,
       so the ExtendInvariants method is called to find the more 
       refined classes that we need
    """

    vgi = openbabel.vectorUnsignedInt()
    natoms = pmol.NumAtoms()
    nfragatoms = frag_atoms.CountBits()
    pmol.GetGIVector(vgi)
    iteratom = openbabel.OBMolAtomIter(pmol)
    for atom in iteratom:
        idx = atom.GetIdx()
        if(frag_atoms.BitIsOn(idx)):
            symmetry_classes.append([atom, vgi[idx-1]])
    nclasses = ExtendInvariants(poltype,pmol, symmetry_classes,frag_atoms,nfragatoms,natoms)
    return nclasses

def ExtendInvariants(poltype,pmol, symmetry_classes,frag_atoms,nfragatoms,natoms):
    """
    Intent: Refine the invariants found by openbabel's GetGIVector
    Input: 
        pmol: OBMol object
        symmetry_classes: symmetry classes array
        frag_atoms: vector containing information about which atoms belong to the largest fragment
        nfragatoms: # of atoms belonging to the largest fragment
        natoms: # of atoms belonging to the molecule
    Output: 
        nclasses1: # of symmetry classes found 
        symmetry_classes: this array is updated
    Referenced By: CalculateSymmetry
    Description:
    
    Description:
    1. Find the # of current classes found by openbabel, nclasses1, 
       and renumber (relabel) the classes to 1, 2, 3, ...
    2. Begin loop
       a. CreateNewClassVector is called which fills in the 'tmp_classes' array with
          a new set of classes by considering bonding information
       b. The number of classes in tmp_classes is found, 'nclasses2'
       c. If there was no change, nclasses1 == nclasses2, break
       d. If the number of classes changed, set nclasses1 to nclasses2, then continue loop
       e. The loop is continued because now that the symmetry classes have changed,
          CreateNewClassVector may find more classes
    3. Return # of classes found
    """
    nclasses1 = CountAndRenumberClasses(poltype,symmetry_classes)
    tmp_classes = []
    if(nclasses1 < nfragatoms):
        #stops when number of classes don't change
        for i in range(100):
            CreateNewClassVector(poltype,pmol,symmetry_classes, tmp_classes, frag_atoms, natoms)
            nclasses2 = CountAndRenumberClasses(poltype,tmp_classes)
            del symmetry_classes[:]
            symmetry_classes.extend(tmp_classes)
            if(nclasses1 == nclasses2):
                break
            nclasses1 = nclasses2
    return nclasses1

def CountAndRenumberClasses(poltype,symmetry_classes):
    """
    Intent: Counts the number of symmetry classes and renumbers them to 1, 2, 3, ...
    Input: 
        symmetry_classes: Array of symmetry classes
    Output: 
        count: # of symmetry classes
        symmetry_classes array is updated
    Referenced By: ExtendInvariants
    Description: -
    """
    count = 1
    symmetry_classes = sorted(symmetry_classes, key=lambda sym: sym[1])
    if(len(symmetry_classes) > 0):
        idatom = symmetry_classes[0][1]
        symmetry_classes[0][1] = 1
        for i in range(1,len(symmetry_classes)):
            if(symmetry_classes[i][1] != idatom):
                idatom = symmetry_classes[i][1]
                count = count + 1
                symmetry_classes[i][1] = count
            else:
                symmetry_classes[i][1] = count
    return count

def CreateNewClassVector(poltype,pmol,symmetry_classes, tmp_classes, frag_atoms, natoms):
    """
    Intent: Find new symmetry classes if possible
    If two atoms were originally of the same sym class but are bound to atoms of differing
    sym classes, then these two atoms will now belong to two different sym classes
    Input:
        pmol: OBMol object
        symmetry_classes: previous set of symmetry classes
        tmp_classes: tmp array of new symmetry classes
        frag_atoms: atoms in largest fragment
        natoms: number of atoms
    Ouptut:
        tmp_classes is edited
    Referenced By: ExtendInvariants
    Description:
    1. dict idx2index is created which maps atom idx's to an index ranging
       from 0 to # of symmetry classes - 1
    2. For each atom a:
           a. For each bond b that a belongs to:
               i. Find the idx of the other atom, nbratom, in the bond b
               ii. Find and append the symmetry class that nbratom belongs to to vtmp
           b. Using vtmp, create a label for atom a
              i. This label contains information about the symmetry classes of the atoms
              that a is bound to
              ii. This label will be different for two atoms that were originally the same 
              symmetry class but are bound to atoms of differing symmetry classes
    """
    idx2index = dict()
    index = 0
    del tmp_classes[:]
    for s in symmetry_classes:
        idx2index.update({s[0].GetIdx() : index})
        index = index + 1
    for s in symmetry_classes:
        iterbond = openbabel.OBMolBondIter(pmol)
        atom = s[0]
        idatom = s[1]
        nbridx =  0
        vtmp = []
        for b in iterbond:
            #if atom belongs to bond b
            if atom.GetIdx() == b.GetEndAtomIdx() or atom.GetIdx() == b.GetBeginAtomIdx():
                if(atom.GetIdx() == b.GetEndAtomIdx()):
                    nbridx = b.GetBeginAtomIdx()
                elif(atom.GetIdx() == b.GetBeginAtomIdx()):
                    nbridx = b.GetEndAtomIdx()
                if(frag_atoms.BitIsOn(nbridx)):
                    vtmp.append(symmetry_classes[idx2index[nbridx]][1])
        vtmp.sort()
        m = 100
        for v in vtmp:
            idatom = idatom + v * m
            m = 100 * m
        tmp_classes.append([atom,idatom])


def gen_canonicallabels(poltype,mol):
    poltype.symmetryclass = [ 0 ] * mol.NumAtoms()
    """
    Intent: Find the symmetry class that each atom belongs to
    Input: 
        mol: OBMol object 
    Output: 
        The global variable 'symmetryclass' is altered
    Referenced By: main
    Description:
    1. An empty bit vector is created, 'frag_atoms'
    2. OBMol.FindLargestFragment is called to fill in the 'frag_atoms' bit vector (the
    vector is filled with a 1 or 0 depending on whether the atom is part of the largest
    fragment or not)
    3. 'CalculateSymmetry' method is called to find initial symmetry classes
    4. Terminal atoms of the same element are collapsed to one symmetry class
    5. Possibly renumber the symmetry classes
    """
    # Returns symmetry classes for each atom ID
    frag_atoms = openbabel.OBBitVec()
    symmclasslist = []
    mol.FindLargestFragment(frag_atoms)
    CalculateSymmetry(poltype,mol, frag_atoms, symmclasslist)
    for ii in range(len(poltype.symmetryclass)):
        poltype.symmetryclass[ii] = symmclasslist[ii][1]

    # Collapse terminal atoms of same element to one type
    for a in openbabel.OBMolAtomIter(mol):
        for b in openbabel.OBAtomAtomIter(a):
            if b.GetValence() == 1:
                for c in openbabel.OBAtomAtomIter(a):
                    if ((b is not c) and
                        (c.GetValence() == 1) and
                        (b.GetAtomicNum() == c.GetAtomicNum()) and
                        (poltype.symmetryclass[b.GetIdx()-1] !=
                            poltype.symmetryclass[c.GetIdx()-1])):
                        poltype.symmetryclass[c.GetIdx()-1] = \
                            poltype.symmetryclass[b.GetIdx()-1]

    # Renumber symmetry classes
    allcls=list(set(poltype.symmetryclass))
    allcls.sort()
    for ii in range(len(poltype.symmetryclass)):
        poltype.symmetryclass[ii] = allcls.index(poltype.symmetryclass[ii]) + 1

def get_symm_class(poltype,x):
    """
    Intent: For a given atom, output it's symmetry class if it exists
    Input: 
        x: atom in question
    Output:
        symmetry class if it exists, -1 if not
    """
    if (str(type(poltype.mol.GetAtom(x))).find("OBAtom") >= 0):
        return poltype.symmetryclass[x-1]
    else:
        return -1

def getCan(poltype,x):
    """
    Intent: For a given atom, output it's canonical label if it exists
    Input: 
        x: atom in question
    Output:
        canonical label if it exists, -1 if not
    """
    if (str(type(mol.GetAtom(x))).find("OBAtom") >= 0):
        return poltype.canonicallabel[x-1]
    else:
        return -1

def get_class_number(poltype,idx):
    """
    Intent: Given an atom idx, return the atom's class number
    """
    maxidx =  max(poltype.symmetryclass)
    return poltype.prmstartidx + (maxidx - poltype.symmetryclass[idx - 1])


