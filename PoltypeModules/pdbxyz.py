import os
import submitjobs as submit
from openbabel import openbabel
import re
import shutil
import warnings
import time


def GenerateProteinTinkerXYZFile(poltype):
    if poltype.uncomplexedproteinpdbname==None:
        proteinindextocoordinates=GenerateUncomplexedProteinPDBFromComplexedPDB(poltype)
    poltype.uncomplexedxyzname=poltype.uncomplexedproteinpdbname.replace('.pdb','.xyz')
    poltype.complexedxyzname=poltype.uncomplexedxyzname.replace('.xyz','_comp.xyz')
    poltype.receptorligandxyzfilename=poltype.complexedxyzname
    poltype.ReadReceptorCharge()
    chainnum=DetectNumberOfChains(poltype,poltype.uncomplexedproteinpdbname)
    resarray=FindCurrentResidueArray(poltype,poltype.uncomplexedproteinpdbname)
    missingresidues=FindMissingResidues(poltype,resarray)

    if not os.path.isfile(poltype.complexedxyzname):
        if chainnum==1:

            cmdstr=poltype.pdbxyzpath+' '+poltype.uncomplexedproteinpdbname+' '+poltype.prmfilepath
        else:
            cmdstr=poltype.pdbxyzpath+' '+poltype.uncomplexedproteinpdbname+' '+'ALL'+' '+poltype.prmfilepath

        submit.call_subsystem(poltype,cmdstr,wait=True)    
    atoms,coord,order,types,connections=readTXYZ(poltype,poltype.uncomplexedxyzname)
    uncomplexedatomnum=len(proteinindextocoordinates.keys())
    newuncomplexedatomnum=len(atoms)
    shift= newuncomplexedatomnum- uncomplexedatomnum
    if shift!=0:
        string='WARNING! Missing atoms from original PDB have been added by PDBXYZ. Number of atoms added = '+str(shift)
        warnings.warn(string)
        poltype.WriteToLog(string)
    if len(missingresidues)!=0:
        string="WARNING! Residues are missing "+str(missingresidues)
        warnings.warn(string)
        poltype.WriteToLog(string)

    indextocoordinates=GrabLigandCoordinates(poltype,uncomplexedatomnum,shift)
    poltype.ligandindices[0]=list(indextocoordinates.keys())
    GenerateComplexedTinkerXYZFile(poltype,poltype.uncomplexedxyzname,indextocoordinates,newuncomplexedatomnum)
    ligandindices=poltype.ligandindices[0]
    GeneratePDBFileFromXYZ(poltype,poltype.complexedxyzname,ligandindices)
 
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


def GeneratePDBFileFromXYZ(poltype,xyzfile,ligandindices):
    cmdstr=poltype.xyzpdbpath+' '+xyzfile+' '+poltype.prmfilepath 
    submit.call_subsystem(poltype,cmdstr,wait=True)
    newpdb=xyzfile.replace('.xyz','.pdb')
    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,newpdb)
    obConversion.SetOutFormat('xyz')
    tempname='temp.xyz'
    obConversion.WriteFile(pdbmol,tempname)
    obConversion.SetInFormat('xyz')
    obConversion.SetOutFormat('pdb')
    finalname='topologyguess.pdb'
    newpdbmol=openbabel.OBMol()
    obConversion.ReadFile(newpdbmol,tempname)
    obConversion.WriteFile(newpdbmol,finalname)
    tempname=finalname.replace('.pdb','_TEMP.pdb')
    temp=open(finalname,'r')
    results=temp.readlines()
    temp.close()
    temp=open(tempname,'w')
    for line in results:
        linesplit=re.split(r'(\s+)', line)
        if 'HETATM' in line:
            othersplit=line.split()
            index=int(othersplit[1])
            if index in ligandindices:
                if 'UNL' in linesplit:
                    idx=linesplit.index('UNL')
                    linesplit[idx]='LIG'
                    line=''.join(linesplit)
        temp.write(line)
    temp.close()
    os.rename(tempname,finalname)
    if not os.path.isdir(poltype.visfolder):
        os.mkdir(poltype.visfolder)
    newpath=os.path.join(poltype.visfolder,finalname)
    shutil.copy(finalname,newpath)

def DetectNumberOfChains(poltype,pdbfile):
    temp=open(pdbfile,'r')
    results=temp.readlines()
    temp.close()
    chainls=[]
    for line in results:
        if 'ATOM' in line:
            linesplit=line.split()
            chain=linesplit[4]
            if chain not in chainls:
                chainls.append(chain)
    chainnum=len(chainls)
    return chainnum


def WriteSEQFile(poltype,filename,code):
    from modeller import Environ,log,Model,Alignment
    from modeller.automodel import LoopModel,refine    # Loa
    e = Environ()
    m = Model(e, file=filename)
    aln = Alignment(e)
    aln.append_model(m, align_codes=code)
    seqfilename=filename.replace('.pdb','.seq')
    aln.write(file=filename.replace('.pdb','.seq'))
    return seqfilename


def FindCurrentResidueArray(poltype,filename):
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    resarray=[]
    for line in results:
        linesplit=line.split()
        if 'ATOM' in line and linesplit[0]=='ATOM':
            resnum=int(linesplit[5])
            if resnum not in resarray:
                resarray.append(resnum)
            
    return resarray


def FindMissingResidues(poltype,resarray):
    missingresidues=[]
    firstres=1 # dont start at first residue found, start at 1
    lastres=resarray[-1]
    allres=list(range(firstres,lastres+1))
    for res in allres:
        if res not in resarray:
            missingresidues.append(res)



    return missingresidues



def GrabLetterSequence(poltype,seqfilename):
    lettersequence=[]
    chainsequence=[]
    temp=open(seqfilename,'r')
    results=temp.readlines()
    temp.close()
    count=0
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=0:
            if len(linesplit)==1:
                letters=linesplit[0]
                array=[]
                allalpha=True
                for e in letters:
                    if e.isalpha():
                        array.append(e)
                    else:
                        if e!='*':
                            allalpha=False
                            break
                if allalpha==True:
                    lettersequence.extend(array)
                    chainseq=[count]*len(array)
                    chainsequence.extend(chainseq)
                    count+=1



    return lettersequence,chainsequence


def GrabMissingLetterSequence(poltype,filename):
    chainlettertonum={'A':0,'B':1,'C':2,'D':3}
    resnumtomissingthreeletter={}
    resnumtomissingchainseq={}
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    foundmissresline=False
    for line in results:
        linesplit=line.split()
        if 'M RES C SSSEQI' in line:
            foundmissresline=True
        if foundmissresline==True and len(linesplit)==5:
            threeletercode=linesplit[2]
            chain=linesplit[-2]
            chainnum=chainlettertonum[chain]
            resnum=int(linesplit[-1])
            resnumtomissingthreeletter[resnum]=threeletercode
            resnumtomissingchainseq[resnum]=chainnum
        elif foundmissresline==True and len(linesplit)!=5:
            if 'M RES C SSSEQI' not in line:
                foundmissresline==False
                break


    return resnumtomissingthreeletter,resnumtomissingchainseq


def ConvertThreeLetterToSingleLetter(poltype,resnumtomissingthreeletter,threelettercodetosinglelettercode):
    resnumtomissingsingleletter={}
    for resnum,threeletter in resnumtomissingthreeletter.items():
        singleletter=threelettercodetosinglelettercode[threeletter]
        resnumtomissingsingleletter[resnum]=singleletter
    return resnumtomissingsingleletter



def CombineData(poltype,resnumtomissingsingleletter,resnumtocurrentsingleletter,resnumtomissingchainseq,resnumtocurrentchainseq):
    resnumtochainseq={}
    resnumtomissing={}
    resnumtosingleletter={}
    for resnum,missingsingleletter in resnumtomissingsingleletter.items():
        resnumtomissing[resnum]=True
        resnumtosingleletter[resnum]=missingsingleletter
        resnumtochainseq[resnum]=resnumtomissingchainseq[resnum]

    for resnum,currentsingleletter in resnumtocurrentsingleletter.items():
        resnumtomissing[resnum]=False
        resnumtosingleletter[resnum]=currentsingleletter
        chainseq=resnumtocurrentchainseq[resnum]
        resnumtochainseq[resnum]=chainseq
            


    return SortByKey(resnumtochainseq),SortByKey(resnumtomissing),SortByKey(resnumtosingleletter)


def SortByKey(dic):
    newdic={}
    for key in sorted(dic):
        newdic[key]=dic[key]
    return newdic



def GenerateGapAndFilledArrays(poltype,resnumtochainseq,resnumtomissing,resnumtosingleletter):
    gapresiduearrays=[]
    filledresiduearrays=[]
    chainseqtoresnums={}
    for resnum,chainseq in resnumtochainseq.items():
        if chainseq not in chainseqtoresnums.keys():
            chainseqtoresnums[chainseq]=[]
        if resnum not in chainseqtoresnums[chainseq]:
            chainseqtoresnums[chainseq].append(resnum)

    for chainseq,resnums in chainseqtoresnums.items():
        sortedresnums=sorted(resnums)
        gapresarray=[]
        filledresarray=[]
        for resnum in sortedresnums:
            missing=resnumtomissing[resnum]
            singleletter=resnumtosingleletter[resnum]
            if missing==True:
                gaplet='-'
                filledlet=singleletter
            else:
                gaplet=singleletter
                filledlet=singleletter
            gapresarray.append(gaplet)
            filledresarray.append(filledlet)
        gapresiduearrays.append(gapresarray)
        filledresiduearrays.append(filledresarray)
            

    return gapresiduearrays,filledresiduearrays


def GenerateAlignmentFile(poltype,resnumtochainseq,resnumtomissing,resnumtosingleletter,seqfilename,code):
    gapresiduearrays,filledresiduearrays=GenerateGapAndFilledArrays(poltype,resnumtochainseq,resnumtomissing,resnumtosingleletter)
    alignmentfilename='alignment.ali'
    temp=open(seqfilename,'r')
    results=temp.readlines()
    temp.close()
    temp=open(alignmentfilename,'w')
    newcode=code+'_filled'
    WriteAlignmentFileComponent(poltype,temp,results,gapresiduearrays,filledresiduearrays,code,newcode,filled=False)
    temp.flush()
    os.fsync(temp.fileno())
    WriteAlignmentFileComponent(poltype,temp,results,gapresiduearrays,filledresiduearrays,code,newcode,filled=True)
    temp.close()
    return alignmentfilename,newcode

def WriteAlignmentFileComponent(poltype,temp,results,gapresiduearrays,filledresiduearrays,code,newcode,filled=False):
    count=0
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=0:
            if len(linesplit)==1:
                letters=linesplit[0]
                array=[]
                allalpha=True
                for e in letters:
                    if e.isalpha():
                        array.append(e)
                    else:
                        if e!='*':
                            allalpha=False
                            break
                if allalpha==True:
                    gapresarray=gapresiduearrays[count]
                    filledresarray=filledresiduearrays[count]
                    gapstring=''.join(gapresarray)
                    filledstring=''.join(filledresarray)
                    if count==len(gapresiduearrays)-1:
                        gapstring+='*'
                        filledstring+='*'
                    gapstring+='\n'
                    filledstring+='\n'
                    if filled==False:
                        temp.write(gapstring) 
                    else:
                        temp.write(filledstring)
                    count+=1
                else:
                    if filled==True:
                        if code in line:
                            line=line.replace(code,newcode)
                    temp.write(line)
            else:
                temp.write(line)


def GenerateLoops(poltype,alignmentfilename,newcode,code):
    from modeller import Environ,log,Model,Alignment
    from modeller.automodel import LoopModel,refine    # Loa
    log.verbose()
    env = Environ()
    env.io.atom_files_directory = ['.', '../atom_files']
    a = LoopModel(env, alnfile = alignmentfilename,knowns = code, sequence = newcode)
    a.starting_model= 1
    a.ending_model  = 1
    
    a.loop.starting_model = 1
    a.loop.ending_model   = 2
    a.loop.md_level       = refine.fast
    
    a.make()


def GrabLastGeneratedPDB(poltype):
    files=os.listdir()
    timetofilename={}
    for f in files:
        if '.pdb' in f:
            Ftime=os.path.getmtime(f)
            reltime=time.time()-Ftime
            timetofilename[reltime]=f
    mintime=min(timetofilename.keys())
    minfilename=timetofilename[mintime]
    return minfilename


def FillInMissingResidues(poltype,code):
    from modeller import Environ,log,Model,Alignment
    from modeller.automodel import LoopModel,refine    # Load the AutoModel class
    from Bio.PDB.PDBList import PDBList

    pdbl = PDBList()
    filename = pdbl.retrieve_pdb_file(code, file_format='pdb')
    if '.ent' in filename:
        newfilename=filename.replace('.ent','.pdb')
        os.rename(filename,filename.replace('.ent','.pdb'))
        filename=newfilename
    threelettercodetosinglelettercode= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    seqfilename=WriteSEQFile(poltype,filename,code)
    resarray=FindCurrentResidueArray(poltype,filename)
    missingresidues=FindMissingResidues(poltype,resarray)
    lettersequence,chainsequence=GrabLetterSequence(poltype,seqfilename)
    resnumtomissingthreeletter,resnumtomissingchainseq=GrabMissingLetterSequence(poltype,filename)
    resnumtocurrentsingleletter=dict(zip(resarray,lettersequence))
    resnumtocurrentchainseq=dict(zip(resarray,chainsequence))
    resnumtomissingsingleletter=ConvertThreeLetterToSingleLetter(poltype,resnumtomissingthreeletter,threelettercodetosinglelettercode)
    resnumtochainseq,resnumtomissing,resnumtosingleletter=CombineData(poltype,resnumtomissingsingleletter,resnumtocurrentsingleletter,resnumtomissingchainseq,resnumtocurrentchainseq)
    alignmentfilename,newcode=GenerateAlignmentFile(poltype,resnumtochainseq,resnumtomissing,resnumtosingleletter,seqfilename,code)
    GenerateLoops(poltype,alignmentfilename,newcode,code)
    finalpdb=GrabLastGeneratedPDB(poltype)



def CallPDB2PQR(poltype,pdbfilename):
    outputfile=pdbfilename.replace('.pdb','.pqr')
    cmdstr='pdb2pqr30'+' '+pdbfilename+' '+outputfile
    os.system(cmdstr)
    finaloutputfile=outputfile.replace('.pqr','_final.pdb')
    ConvertPQRToPDB(poltype,outputfile,finaloutputfile)


def ConvertPQRToPDB(poltype,outputfile,finaloutputfile):
    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    inFormat = obConversion.FormatFromExt(outputfile)
    obConversion.SetInFormat(inFormat)
    obConversion.ReadFile(mol, outputfile)
    obConversion.SetOutFormat('pdb')
    obConversion.WriteFile(mol,finaloutputfile)

