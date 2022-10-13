import torch
import torchani
from openbabel import openbabel
from ase.optimize import BFGS
from ase import Atoms, Atom
from ase.io.trajectory import Trajectory
from ase.constraints import FixInternals
import getopt
import sys
from ase.io import read, write
import shutil

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = torchani.models.ANI2x(periodic_table_index=True).to(device)

def GenerateMOLObject(xyzfile):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    obConversion.SetInFormat('xyz')
    obConversion.ReadFile(mol,xyzfile)
    return mol


def BabelCoordsAtomicNums(mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomvecls=[]
    atomicnums=[]
    iteratombab = openbabel.OBMolAtomIter(mol)
    for atm in iteratombab:
        atmindex=atm.GetIdx()
        coords=[atm.GetX(),atm.GetY(),atm.GetZ()]
        atomvecls.append(coords)
        atomicnum=atm.GetAtomicNum()
        atomicnums.append(atomicnum)

    return atomvecls,atomicnums



def GenerateSpeciesAndCoordinatesForANI(xyzfile):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    mol=GenerateMOLObject(xyzfile)
    atomicvecls,atomicnums=BabelCoordsAtomicNums(mol)
    coordinates = torch.tensor([atomicvecls],requires_grad=True, device=device)
    species = torch.tensor([atomicnums], device=device)
    return coordinates,species

def GrabANIEnergy(coordinates,species,outputlogname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    Hartree2kcal_mol=627.5095
    energy = model((species, coordinates)).energies
    energy=energy.item()*Hartree2kcal_mol
    GenerateSPLogFile(outputlogname,energy)


def GenerateSPLogFile(logname,energy):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(logname,'w')
    temp.write("Normal termination"+'\n') 
    temp.write('Energy: '+str(energy)+' kcal/mol'+'\n')
    temp.close()


def OptimizeANI(xyzfile,outputfile,outputoptxyz,fmax,listofindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    molecule = read(xyzfile)
    molecule.set_calculator(torchani.models.ANI2x().ase())

    #molecule = Atoms(atomicsymbs,positions=newatomicvecls,calculator=torchani.models.ANI2x().ase())
    ls=[]
    for indices in listofindices: 
        dihedral = [molecule.get_dihedral(*indices), indices] 
        ls.append(dihedral)
    c = FixInternals(dihedrals_deg=ls)
    molecule.set_constraint(c)
    dyn = BFGS(molecule,trajectory=outputfile)
    try:
        dyn.run(fmax=fmax)
        traj = Trajectory(outputfile)
        atoms=traj[-1]
        pos=atoms.get_positions()
        GenerateOptimizedStructureXYZ(atomicsymbs,pos,outputoptxyz)
        GenerateOptimizedLogFile(outputoptxyz.replace('.xyz','.log')) # just need to make fake log so poltype can read termination signal
    except:
        GenerateOptimizedLogFile(outputoptxyz.replace('.xyz','.log'),error=True)


def GenerateOptimizedStructureXYZ(atomicsymbs,pos,outputoptxyz):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(outputoptxyz,'w')
    temp.write(str(len(atomicsymbs))+'\n')
    temp.write('\n')
    for i in range(len(atomicsymbs)):
        atomicsymb=atomicsymbs[i]
        currentpos=pos[i]
        line=str(atomicsymb)+' '+str(currentpos[0])+' '+str(currentpos[1])+' '+str(currentpos[2])+' '+'\n'
        temp.write(line)

    temp.close()

def GenerateOptimizedLogFile(logname,error=None):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(logname,'w')
    if error==None:
        temp.write("Normal termination"+'\n')
    else:
        temp.write("error"+'\n')

    temp.close()



def SanitizeXYZ(inputxyz):
    if '.xyz_2' in inputxyz:
        newxyz=inputxyz.replace('.xyz_2','.xyz')
        shutil.copy(inputxyz,newxyz)
        return newxyz
    else:
        return inputxyz


def ReadRestraintFile(inputres):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    listofindices=[]
    temp=open(inputres,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)>1:
            indices=[int(i) for i in linesplit]
            listofindices.append(indices)


    return listofindices

opt=False
opt_list, args = getopt.getopt(sys.argv[1:], 'p:pdb:',['inputxyz=','inputres=','fmax=','outputname=','opt=','outputlogname='])
for o, a in opt_list:
    if o in ('--inputxyz'):
        inputxyz=a
    if o in ('--inputres'):
        inputres=a
    if o in ('--fmax'):
        fmax=float(a)
    if o in ('--outputname'):
        outputname=a
    if o in ('--opt'):
        opt=bool(a)
    if o in ('--outputlogname'):
        outputlogname=a


inputxyz=SanitizeXYZ(inputxyz)
if opt==True:
    listofindices=ReadRestraintFile(inputres)
    outputfile=inputxyz.replace('.xyz','.traj')
    OptimizeANI(inputxyz,outputfile,outputname,fmax,listofindices)
else:
    coordinates,species=GenerateSpeciesAndCoordinatesForANI(inputxyz)
    GrabANIEnergy(coordinates,species,outputlogname)


