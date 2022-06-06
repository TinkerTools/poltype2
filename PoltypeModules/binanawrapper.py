import sys
import binana
import py3Dmol
from openbabel import openbabel
import numpy as np


def ExtractLigand(ligandreceptorfilename):
    ligandpdbfilename='ligand.pdb'
    receptorpdbfilename='receptor.pdb'
    receptormol=ExtractMOLObject(ligandreceptorfilename,receptorpdbfilename,receptor=True)
    ligandmol=ExtractMOLObject(ligandreceptorfilename,ligandpdbfilename,receptor=False)


    return ligandpdbfilename,receptorpdbfilename


def ExtractMOLObject(ligandreceptorfilename,newpdbfilename,receptor):
    pdbmol=GenerateMOLObject(ligandreceptorfilename)
    iteratom = openbabel.OBMolAtomIter(pdbmol)
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat('pdb')
    atmindicestodelete=[]

    for atm in iteratom:
        atmindex=atm.GetIdx()
        res=atm.GetResidue()
        reskey=res.GetName()

        if reskey=='LIG' and receptor==True:
            #pdbmol.DeleteAtom(atm)
            atmindicestodelete.append(atmindex)
        elif reskey!='LIG' and receptor==False:
            atmindicestodelete.append(atmindex)
            #pdbmol.DeleteAtom(atm)
    atmindicestodelete.sort(reverse=True)
    for atmindex in atmindicestodelete:
        atm=pdbmol.GetAtom(atmindex)
        pdbmol.DeleteAtom(atm)
    obConversion.WriteFile(pdbmol,newpdbfilename)
    return pdbmol



def GrabAtomValues(listofdics,key):
    values=[]
    for subitm in listofdics:
        for subkey,subvalue in subitm.items():
                if subkey==key and subvalue not in values:
                    values.append(subvalue)

    return values

def GrabFromToAtomIndex(fromatomindices,fromatomnames,toatomindices,listofelements):
    for j in range(len(fromatomindices)):
        atomindex=fromatomindices[j]
        atomname=fromatomnames[j]
        foundele=False
        for ele in listofelements:
            if ele in atomname:
                foundele=True
        if foundele==False:
            fromatomindex=atomindex

        toatomindex=toatomindices[0]
    return fromatomindex,toatomindex


def GrabDirectedBondFromToAtomIndices(allligandatomindices,allligandatomnames,allreceptoratomindices,allreceptoratomnames,listofelements):
    listoffromtoatomindices=[]
    for i in range(len(allligandatomindices)):
        ligandatomindices=allligandatomindices[i]
        ligandatomnames=allligandatomnames[i]
        receptoratomindices=allreceptoratomindices[i]
        receptoratomnames=allreceptoratomnames[i]
        if len(ligandatomindices)>1: # then contains H and start from here
            fromatomindex,toatomindex=GrabFromToAtomIndex(ligandatomindices,ligandatomnames,receptoratomindices,listofelements)
        elif len(receptoratomindices)>1:
            fromatomindex,toatomindex=GrabFromToAtomIndex(receptoratomindices,receptoratomnames,ligandatomindices,listofelements)
        listoffromtoatomindices.append([fromatomindex,toatomindex])


    return listoffromtoatomindices

def GrabLigandReceptorAtomIndices(data,datakey,trueligandatomindices,residcollector):
    allligandatomindices=[]
    allreceptoratomindices=[]
    allligandatomnames=[]
    allreceptoratomnames=[]
    for itm in data[datakey]:
        for key,value in itm.items():
            if key=='ligandAtoms':
                ligandatomindices=GrabAtomValues(value,'atomIndex')
                allligandatomindices.append(ligandatomindices)
                ligandatomnames=GrabAtomValues(value,'atomName')
                allligandatomnames.append(ligandatomnames)
            elif key=='receptorAtoms':
                resid=GrabAtomValues(value,'resID')
                receptoratomindices=GrabAtomValues(value,'atomIndex')
                allreceptoratomindices.append(receptoratomindices)
                receptoratomnames=GrabAtomValues(value,'atomName')
                allreceptoratomnames.append(receptoratomnames)
                for residnum in resid:
                    if residnum not in residcollector:
                        residcollector.append(residnum)


    allligandatomindices=ShiftIndices(allligandatomindices,trueligandatomindices)
    return allligandatomindices,allreceptoratomindices,allligandatomnames,allreceptoratomnames,residcollector

def WritePDBFile(pdb_txt,pdbfilename):
    temp=open(pdbfilename,'w')
    for l in pdb_txt.split("\n") :
        temp.write(l+'\n')
    temp.close()

def GenerateMOLObject(pdbfilename):
    obConversion = openbabel.OBConversion()
    pdbmol = openbabel.OBMol()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,pdbfilename)
    return pdbmol

def GrabAtomVectorsFromIndices(listoffromtoatomindices,pdbmol):
    listoffromtoatomvecs=[]
    for atomindices in listoffromtoatomindices:
        atomveclist=[]
        for atomindex in atomindices:
            atom=pdbmol.GetAtom(atomindex)
            vec=[atom.GetX(),atom.GetY(),atom.GetZ()]
            atomveclist.append(vec)
        listoffromtoatomvecs.append(atomveclist)

    return listoffromtoatomvecs

def GrabLigandAtomIndices(pdbmol):
    ligandatomindices=[]
    iteratom = openbabel.OBMolAtomIter(pdbmol)
    for atm in iteratom:
        atmindex=atm.GetIdx()
        res=atm.GetResidue()
        reskey=res.GetName()
        if reskey=='LIG':
            ligandatomindices.append(atmindex)

    return ligandatomindices



def GrabResidueNames(pdbmol,residcollector):
    residtoresname={}
    iteratom = openbabel.OBMolAtomIter(pdbmol)
    passedlig=False
    for atm in iteratom:
        atmindex=atm.GetIdx()
        res=atm.GetResidue()
        reskey=res.GetName()
        resid=res.GetNum()
        if reskey=='LIG':
            passedlig=True
        if resid in residcollector and passedlig==False:
            residtoresname[resid]=reskey


    return residtoresname


def ShiftIndices(allligandatomindices,ligandatomindices):
    newallligandatomindices=[]
    for atomls in allligandatomindices:
        newatomls=[]
        for atomindex in atomls:
            zeroindex=atomindex-1
            newatomindex=ligandatomindices[zeroindex]
            newatomls.append(newatomindex)
        newallligandatomindices.append(newatomls)


    return newallligandatomindices


def ComputePositionCenter(allatomvecs):
    allcentroids=[]
    for atomvecls in allatomvecs:
        centroid=np.array([0.0,0.0,0.0])
        for vec in atomvecls:
            centroid+=np.array(vec)
        centroid=centroid/len(atomvecls)
        allcentroids.append(centroid)

    return allcentroids



def CombineLigandReceptorVecs(ligandvecs,receptorvecs):
    listoffromtoatomvecs=[]
    for i in range(len(ligandvecs)):
        ligandvec=ligandvecs[i]
        receptorvec=receptorvecs[i]
        listoffromtoatomvecs.append([ligandvec,receptorvec])

    return listoffromtoatomvecs


def GrabToFromVectorsFromCentroid(data,datakey,ligandatomindices,residcollector,pdbmol):
    allligandatomindices,allreceptoratomindices,allligandatomnames,allreceptoratomnames,residcollector=GrabLigandReceptorAtomIndices(data,datakey,ligandatomindices,residcollector)
    allligandatomvecs=GrabAtomVectorsFromIndices(allligandatomindices,pdbmol)
    allreceptoratomvecs=GrabAtomVectorsFromIndices(allreceptoratomindices,pdbmol)
    ligandcentroid=ComputePositionCenter(allligandatomvecs)
    receptorcentroid=ComputePositionCenter(allreceptoratomvecs)
    listoffromtoatomvecs=CombineLigandReceptorVecs(ligandcentroid,receptorcentroid)
    return residcollector,listoffromtoatomvecs




def GrabToFromVectorsFromDirectedBonds(data,datakey,ligandatomindices,residcollector,pdbmol,listofelements):
    allligandatomindices,allreceptoratomindices,allligandatomnames,allreceptoratomnames,residcollector=GrabLigandReceptorAtomIndices(data,datakey,ligandatomindices,residcollector)
    listoffromtoatomindices=GrabDirectedBondFromToAtomIndices(allligandatomindices,allligandatomnames,allreceptoratomindices,allreceptoratomnames,listofelements)
    listoffromtoatomvecs=GrabAtomVectorsFromIndices(listoffromtoatomindices,pdbmol)
    return residcollector,listoffromtoatomvecs



def CombineArrays(listoflists):
    hydrophobicatomindices=[]
    for largels in listoflists:
        for ls in largels:
            for itm in ls:
                if itm not in hydrophobicatomindices:
                    hydrophobicatomindices.append(itm)


    return hydrophobicatomindices



def CreateViewDictionaries(listoffromtoatomvecshbond,listoffromtoatomvecspipi,listoffromtoatomvecststack,listoffromtoatomvecscatpi,listoffromtoatomvecshalbond,listoffromtoatomvecssaltbridge):
    colordic={"CCN":'red','SAL':'red','PIS':'blue','PIT':"aqua","PIC":'navy',"MTL":'purple',"HYD":'grey',"HBN":"black","HALBN":'green'}
    dasheddic={"CCN":False,'SAL':True,'PIS':True,'PIT':True,"PIC":True,"MTL":False,"HYD":False,"HBN":False,'HALBN':False}
    labeldic={"CCN":'Closest contacts','SAL':'Salt bridge','PIS':'pi-pi stacking','PIT':"T-stacking","PIC":'cation-pi',"MTL":'metal coordination',"HYD":'hydrophobic contacts',"HBN":"hydrogen bond",'HALBN':'halogen bond'}
    tofromdic={"CCN":[],'SAL':listoffromtoatomvecssaltbridge,'PIS':listoffromtoatomvecspipi,'PIT':listoffromtoatomvecststack,"PIC":listoffromtoatomvecscatpi,"MTL":[],"HYD":[],"HBN":listoffromtoatomvecshbond,'HALBN':listoffromtoatomvecshalbond}
    return colordic,labeldic,tofromdic,dasheddic


def CreateViewer(listoffromtoatomvecshbond,listoffromtoatomvecspipi,listoffromtoatomvecststack,listoffromtoatomvecscatpi,listoffromtoatomvecshalbond,listoffromtoatomvecssaltbridge,residcollector,hydrophobicatomvecs,residtoresname):
    ligandinputfilename='ligand.pdb'
    receptorinputfilename='receptor.pdb'
    ligand, receptor = binana.load_ligand_receptor.from_files(ligandinputfilename,receptorinputfilename)
    all_inf = binana.interactions.get_all_interactions(ligand, receptor)
    pdbtxt = binana.output.pdb_file.write_all(ligand, receptor,all_inf,None,as_str=True)
    view = py3Dmol.view(width=1000,height=500)
    view.removeAllModels()
    view.addModel(pdbtxt)
    colordic,labeldic,tofromdic,dasheddic=CreateViewDictionaries(listoffromtoatomvecshbond,listoffromtoatomvecspipi,listoffromtoatomvecststack,listoffromtoatomvecscatpi,listoffromtoatomvecshalbond,listoffromtoatomvecssaltbridge)
    view.setStyle({}, {"cartoon": {'color': 'spectrum'}})
    view.setStyle({"resn": "LIG"}, {'stick':{'radius':0.3}})
    for resid in residcollector:
        view.addStyle({'resi':int(resid)},{'stick':{'colorscheme':'greenCarbon'}})
        resname=residtoresname[resid]
        view.addLabel(resname+str(resid),{'fontOpacity':.8,'backgroundColor':'white','fontColor':'black','showBackground':'false'},{'resi':int(resid)})
    for vec in hydrophobicatomvecs:
        view.addSphere({'center':{'x':vec[0],'y':vec[1],'z':vec[2]},'radius':1.0,'color':'grey','opacity': 0.6})
    for key,fromto in tofromdic.items():
        col=colordic[key]
        dashed=dasheddic[key]
        value=labeldic[key]
        if key=='HYD':
            print(value+' '+'= '+col+' spheres ')
        if len(fromto)!=0:
            print(value+' '+'= '+col)
            for fromtols in fromto:
                fromvec=fromtols[0]
                tovec=fromtols[1]
                if dashed==True:
                    view.addCylinder({"start": dict(x=fromvec[0], y=fromvec[1], z=fromvec[2]),"end":   dict(x=tovec[0], y=tovec[1], z=tovec[2]),"color": col,"radius": .15,"dashed": dashed,"fromCap": 1,"toCap": 1})
                else:
                    view.addArrow({"start": dict(x=fromvec[0], y=fromvec[1], z=fromvec[2]),"end":   dict(x=tovec[0], y=tovec[1], z=tovec[2]),"color": col,"radius": .15,"dashed": dashed,"fromCap": 1,"toCap": 1})
    return view

def PrepareInputs(ligandreceptorfilename):
    ligandinputfilename,receptorinputfilename=ExtractLigand(ligandreceptorfilename)
    residcollector=[]
    ligand, receptor = binana.load_ligand_receptor.from_files(ligandinputfilename,receptorinputfilename)
    all_inf = binana.interactions.get_all_interactions(ligand, receptor)
    pdbtxt = binana.output.pdb_file.write_all(ligand, receptor,all_inf,None,as_str=True)
    pdbfilename='output.pdb'
    WritePDBFile(pdbtxt,pdbfilename)
    pdbmol=GenerateMOLObject(pdbfilename)
    ligandatomindices=GrabLigandAtomIndices(pdbmol)
    hbond_inf = binana.interactions.get_hydrogen_bonds(ligand, receptor)
    halbond_inf = binana.interactions.get_halogen_bonds(ligand, receptor)
    pipi_inf = binana.interactions.get_pi_pi(ligand, receptor)
    cationpi_inf = binana.interactions.get_cation_pi(ligand, receptor)
    salt_inf=binana.interactions.get_salt_bridges(ligand, receptor)
    hyd_inf=binana.interactions.get_hydrophobics(ligand, receptor)
    data = binana.output.dictionary.collect(hydrogen_bonds=hbond_inf,pi_pi=pipi_inf,cat_pi=cationpi_inf,halogen_bonds=halbond_inf,salt_bridges=salt_inf,hydrophobics=hyd_inf)
    residcollector,listoffromtoatomvecshbond=GrabToFromVectorsFromDirectedBonds(data,"hydrogenBonds",ligandatomindices,residcollector,pdbmol,['H'])
    residcollector,listoffromtoatomvecspipi=GrabToFromVectorsFromCentroid(data,'piPiStackingInteractions',ligandatomindices,residcollector,pdbmol)
    residcollector,listoffromtoatomvecststack=GrabToFromVectorsFromCentroid(data,'tStackingInteractions',ligandatomindices,residcollector,pdbmol)
    residcollector,listoffromtoatomvecscatpi=GrabToFromVectorsFromCentroid(data,'cationPiInteractions',ligandatomindices,residcollector,pdbmol)
    residcollector,listoffromtoatomvecshalbond=GrabToFromVectorsFromDirectedBonds(data,"halogenBonds",ligandatomindices,residcollector,pdbmol,['F','Cl','I','Br'])
    residcollector,listoffromtoatomvecssaltbridge=GrabToFromVectorsFromCentroid(data,'saltBridges',ligandatomindices,residcollector,pdbmol)
    allligandatomindiceshydpho,allreceptoratomindiceshydpho,allligandatomnameshydpho,allreceptoratomnameshydpho,residcollector=GrabLigandReceptorAtomIndices(data,'hydrophobicContacts',ligandatomindices,residcollector)
    hydrophobicatomindices=CombineArrays([allligandatomindiceshydpho,allreceptoratomindiceshydpho])
    hydrophobicatomvecs=GrabAtomVectorsFromIndices([hydrophobicatomindices],pdbmol)[0]
    residtoresname=GrabResidueNames(pdbmol,residcollector)
    return listoffromtoatomvecshbond,listoffromtoatomvecspipi,listoffromtoatomvecststack,listoffromtoatomvecscatpi,listoffromtoatomvecshalbond,listoffromtoatomvecssaltbridge,residcollector,hydrophobicatomvecs,residtoresname
