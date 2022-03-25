import os
import submitjobs as submit
import openbabel
import re

def GenerateProteinTinkerXYZFile(poltype):
    if poltype.uncomplexedproteinpdbname==None:
        proteinindextocoordinates=GenerateUncomplexedProteinPDBFromComplexedPDB(poltype)
    poltype.uncomplexedxyzname=poltype.uncomplexedproteinpdbname.replace('.pdb','.xyz')
    poltype.complexedxyzname=poltype.uncomplexedxyzname.replace('.xyz','_comp.xyz')
    poltype.receptorligandxyzfilename=poltype.complexedxyzname
    poltype.ReadReceptorCharge()

    if not os.path.isfile(poltype.complexedxyzname):
        cmdstr=poltype.pdbxyzpath+' '+poltype.uncomplexedproteinpdbname+' '+poltype.prmfilepath
        submit.call_subsystem(poltype,cmdstr,wait=True)    
    atoms,coord,order,types,connections=readTXYZ(poltype,poltype.uncomplexedxyzname)
    uncomplexedatomnum=len(proteinindextocoordinates.keys())
    newuncomplexedatomnum=len(atoms)
    shift= newuncomplexedatomnum- uncomplexedatomnum
    indextocoordinates=GrabLigandCoordinates(poltype,uncomplexedatomnum,shift)
    poltype.ligandindices[0]=list(indextocoordinates.keys())
    GenerateComplexedTinkerXYZFile(poltype,poltype.uncomplexedxyzname,indextocoordinates,newuncomplexedatomnum) 
 
def readTXYZ(poltype,TXYZ):
    temp=open(TXYZ,'r')
    lines = temp.readlines()[1:] #TINKER coordinate starts from second line
    atoms=[];coord=[]
    order=[];types=[];connections=[]
    for line in lines:
        data=line.split()
        order.append(data[0])
        types.append(data[5])
        connections.append(data[6:])
        atoms.append(data[1])
        coord.append([float(data[2]), float(data[3]), float(data[4])])
    return atoms,coord,order, types, connections


def GrabLigandCoordinates(poltype,uncomplexedatomnum,shift): # assumes appended to end of PDB file
    indextocoordinates={}
    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,poltype.complexedproteinpdbname)
    atomiter=openbabel.OBMolAtomIter(pdbmol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        if atomidx>uncomplexedatomnum:
            coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
            indextocoordinates[atomidx+shift]=coords
    if len(indextocoordinates.keys())==0:
        raise ValueError('Complexed PDB missing atoms or ligand in Complexed PDB is not appended to the end of file')
    return indextocoordinates    

def GenerateComplexedTinkerXYZFile(poltype,uncomplexedxyzname,indextocoordinates,uncomplexedatomnum):
    temp=open(uncomplexedxyzname,'r')
    results=temp.readlines()
    temp.close() 
    atoms,coord,order,types,connections=readTXYZ(poltype,poltype.ligandxyzfilename)
    temp=open(poltype.complexedxyzname,'w')
    newatomnum=uncomplexedatomnum+len(indextocoordinates)
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if lineidx==0:
            newline=str(newatomnum)+'\n' 
        else:
            newline=line
        temp.write(newline)
    for idx in range(len(atoms)):
        element=atoms[idx]
        oldindex=order[idx]
        newindex=int(oldindex)+uncomplexedatomnum
        coords=indextocoordinates[newindex]
        x=coords[0]
        y=coords[1]
        z=coords[2]
        typenum=types[idx]
        conns=connections[idx]
        conns=[int(i) for i in conns]
        conns=[i+uncomplexedatomnum for i in conns]
        newline='    '+str(newindex)+'  '+element+'     '+str(x)+'   '+str(y)+'   '+str(z)+'    '+str(typenum)+'     '
        for con in conns:
            newline+=str(con)+'     '
        newline+='\n'
        temp.write(newline)
    temp.close()

def GenerateUncomplexedProteinPDBFromComplexedPDB(poltype):
    uncomplexedatomindices=[]
    indextocoordinates={}
    temp=open(poltype.complexedproteinpdbname,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'ATOM' in line:
            linesplit=line.split()
            index=int(linesplit[1])
            uncomplexedatomindices.append(index) 

    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,poltype.complexedproteinpdbname)
    totalatoms=pdbmol.NumAtoms()
    atomiter=openbabel.OBMolAtomIter(pdbmol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        if atomidx in uncomplexedatomindices:
            coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
            indextocoordinates[atomidx]=coords

    
    indexestodelete=[]
    for i in range(1,totalatoms+1):
        if i not in uncomplexedatomindices:
            indexestodelete.append(i)
    indexestodelete.sort(reverse=True)
    for idx in indexestodelete:
        atom=pdbmol.GetAtom(idx)
        pdbmol.DeleteAtom(atom)
    
    obConversion.SetOutFormat('pdb')
    poltype.uncomplexedproteinpdbname='uncomplexed.pdb'
    obConversion.WriteFile(pdbmol,poltype.uncomplexedproteinpdbname)
    return indextocoordinates



