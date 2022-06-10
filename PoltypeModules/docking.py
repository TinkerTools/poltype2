from openbabel import openbabel
import re
import os
import numpy as np

'''
conda install -c hcc autodock --yes
conda install -c anaconda pip --yes
pip install vina
conda install -c conda-forge openbabel --yes
'''



def ExtractLigand(ligandreceptorfilename):
    ligandpdbfilename='ligand.pdb'
    receptorpdbfilename='receptor.pdb'
    receptormol=ExtractMOLObject(ligandreceptorfilename,receptorpdbfilename,receptor=True)
    ligandmol=ExtractMOLObject(ligandreceptorfilename,ligandpdbfilename,receptor=False)

    ConvertUNLToLIG(ligandpdbfilename)
    return ligandpdbfilename,receptorpdbfilename

def GenerateMOLObject(pdbfilename):
    obConversion = openbabel.OBConversion()
    pdbmol = openbabel.OBMol()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,pdbfilename)
    return pdbmol

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
            atmindicestodelete.append(atmindex)
        elif reskey!='LIG' and receptor==False:
            atmindicestodelete.append(atmindex)
    atmindicestodelete.sort(reverse=True)
    for atmindex in atmindicestodelete:
        atm=pdbmol.GetAtom(atmindex)
        pdbmol.DeleteAtom(atm)
    obConversion.WriteFile(pdbmol,newpdbfilename)
    return pdbmol


def ConvertUNLToLIG(filename):
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    tempname=filename.replace('.pdb','_TEMP.pdb')
    temp=open(tempname,'w')
    for line in results:
        linesplit=re.split(r'(\s+)', line)
        if 'UNL' in line:
            lineindex=linesplit.index('UNL')
            linesplit[lineindex]='LIG'
            line=''.join(linesplit)
        if 'UNK' in line:
            lineindex=linesplit.index('UNK')
            linesplit[lineindex]='LIG'
            line=''.join(linesplit)

        temp.write(line)
    temp.close()
    os.rename(tempname,filename)


def GenerateListString(ls):
    ls=[str(i) for i in ls]
    string=','.join(ls)
    return string

def GenerateCommandString(python2path,preparescript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename):
    cmdstr=''
    dockinggridstring=GenerateListString(dockgridcenter)
    dockinggridsizestring=GenerateListString(dockgridsize)
    cmdstr+=python2path+' '+preparescript+' '+'--dockgridcenter='+dockinggridstring+' '+'--dockgridsize='+dockinggridsizestring+' '+'--spacing='+str(gridspacing)+' '+'--ligandname='+ligandpdbfilename+' '+'--receptorname='+receptorpdbfilename
    return cmdstr


def GenerateGridFile(receptorgridname):
    cmdstr='autogrid4 '+'-p '+receptorgridname
    os.system(cmdstr)


def RunAutoDock4(ad4parameterfilename):
    cmdstr='autodock4 '+'-p '+ad4parameterfilename
    os.system(cmdstr)

def GeneratePDBQTFiles(python2path,prepdockscript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename):

    cmdstr=GenerateCommandString(python2path,prepdockscript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename)

    os.system(cmdstr)



def AutoDock4(python2path,prepdockscript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename):
    receptorname=receptorpdbfilename.replace('.pdb','.pdbqt')
    ligandname=ligandpdbfilename.replace('.pdb','.pdbqt')
    receptorgridname=receptorpdbfilename.replace('.pdb','.gpf')
    ad4parameterfilename=ligandpdbfilename.replace('.pdb','')+'_'+receptorpdbfilename.replace('.pdb','')+'.dpf'
    GenerateGridFile(receptorgridname)
    RunAutoDock4(ad4parameterfilename)
    ad4logfile=ad4parameterfilename.replace('.dpf','.dlg')
    modeltoscore,modeltostructure=GenerateIndividualLigandFiles(ad4logfile,'ad4',ligandname)
    return modeltoscore,modeltostructure


def GenerateIndividualLigandFiles(outputposes,dockingprogram,ligandname):
    temp=open(outputposes,'r')
    results=temp.readlines()
    temp.close()
    modelcount=1
    modeltoscore={}
    modeltostructure={}
    name=ligandname.replace('.pdbqt','')+'_'+dockingprogram+'_'+str(modelcount)+'.pdbqt'
    modeltostructure[modelcount]=name
    temp=open(name,'w')
    for line in results:
        if 'ENDMDL' in line:
            if dockingprogram=='ad4' and 'DOCK' not in line:
                continue
            if 'DOCK' in line:
                truelinesplit=re.split(r'(\s+)', line)
                truelinesplit=truelinesplit[2:]
                line=''.join(truelinesplit)

            temp.write(line)
            temp.close()
            modelcount+=1
            name=ligandname.replace('.pdbqt','')+'_'+dockingprogram+'_'+str(modelcount)+'.pdbqt'
            temp=open(name,'w')
            modeltostructure[modelcount]=name
        elif "RESULT" in line:
            temp.write(line)
            linesplit=line.split()
            score=float(linesplit[3])
            modeltoscore[modelcount]=score
        elif 'Estimated Free Energy' in line:
            truelinesplit=re.split(r'(\s+)', line)
            truelinesplit=truelinesplit[2:]
            line=''.join(truelinesplit)

            temp.write(line)
            linesplit=line.split()
            score=float(linesplit[-3])
            modeltoscore[modelcount]=score

        else:
            if "REMARK" in line or "ATOM" in line or 'BRANCH' in line or 'ROOT' in line or 'TORSDOF' in line or 'TER' in line or "HETATM" in line:
                if dockingprogram=='ad4':
                    if 'DOCK' in line:
                        truelinesplit=re.split(r'(\s+)', line)
                        truelinesplit=truelinesplit[2:]
                        line=''.join(truelinesplit)
                        temp.write(line)
                else:
                    temp.write(line)
    return modeltoscore,modeltostructure

def AutoDockVinaVinardo(dockingprogram,receptorname,ligandname,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes):
    outputposes=ligandname.replace('.pdbqt','_'+dockingprogram+'_out'+'.pdbqt')
    from vina import Vina
    v = Vina(sf_name=dockingprogram,verbosity=1)
    v.set_receptor(receptorname)
    v.set_ligand_from_file(ligandname)
    v.compute_vina_maps(center=dockgridcenter, box_size=dockgridsize)
    v.dock(exhaustiveness=vinaexhaustiveness, n_poses=nposes)
    v.write_poses(outputposes, n_poses=nposes, overwrite=True)
    modeltoscore,modeltostructure=GenerateIndividualLigandFiles(outputposes,dockingprogram,ligandname)
    return modeltoscore,modeltostructure


def GenerateGOLDConfigurationFile(ligandfilename,receptorfilename,dockgridcenter,dockgridsize,nposes):
    radius=dockgridsize[0]
    conffile='gold.conf'
    temp=open(conffile,'w')
    temp.write('  GOLD CONFIGURATION FILE'+'\n')
    temp.write(' AUTOMATIC SETTINGS'+'\n')
    temp.write('autoscale = 1'+'\n')
    temp.write(' POPULATION'+'\n')
    temp.write('popsiz = auto'+'\n')
    temp.write('select_pressure = auto'+'\n')
    temp.write('n_islands = auto'+'\n')
    temp.write('maxops = auto'+'\n')
    temp.write('niche_siz = auto'+'\n')
    temp.write(' GENETIC OPERATORS'+'\n')
    temp.write('pt_crosswt = auto'+'\n')
    temp.write('allele_mutatewt = auto'+'\n')
    temp.write('migratewt = auto'+'\n')
    temp.write(' FLOOD FILL'+'\n')
    temp.write('radius = '+str(radius)+'\n')
    temp.write('origin = '+str(dockgridcenter[0])+' '+str(dockgridcenter[1])+' '+str(dockgridcenter[2])+'\n')
    temp.write('do_cavity = 1'+'\n')
    temp.write('floodfill_atom_no = 0'+'\n')
    temp.write(' DATA FILES'+'\n')
    temp.write('ligand_data_file '+ligandfilename+' '+str(nposes)+'\n')
    temp.write('param_file = DEFAULT'+'\n')
    temp.write('set_ligand_atom_types = 1'+'\n')
    temp.write('set_protein_atom_types = 0'+'\n')
    temp.write('tordist_file = DEFAULT'+'\n')
    temp.write('make_subdirs = 0'+'\n')
    temp.write('save_lone_pairs = 1'+'\n')
    temp.write('fit_points_file = fit_pts.mol2'+'\n')
    temp.write('read_fitpts = 0'+'\n')
    temp.write(' FLAGS'+'\n')
    temp.write('internal_ligand_h_bonds = 0'+'\n')
    temp.write('flip_free_corners = 1'+'\n')
    temp.write('match_ring_templates = 0'+'\n')
    temp.write('flip_amide_bonds = 0'+'\n')
    temp.write('flip_planar_n = 1 flip_ring_NRR flip_ring_NHR'+'\n')
    temp.write('flip_pyramidal_n = 0'+'\n')
    temp.write('rotate_carboxylic_oh = flip'+'\n')
    temp.write('use_tordist = 1'+'\n')
    temp.write('postprocess_bonds = 1'+'\n')
    temp.write('rotatable_bond_override_file = DEFAULT'+'\n')
    temp.write('solvate_all = 1'+'\n')
    temp.write(' TERMINATION'+'\n')
    temp.write('early_termination = 1'+'\n')
    temp.write('n_top_solutions = 3'+'\n')
    temp.write('rms_tolerance = 1.5'+'\n')
    temp.write(' CONSTRAINTS'+'\n')
    temp.write('force_constraints = 0'+'\n')
    temp.write(' COVALENT BONDING'+'\n')
    temp.write('covalent = 0'+'\n')
    temp.write(' SAVE OPTIONS'+'\n')
    temp.write('save_score_in_file = 1'+'\n')
    temp.write('save_protein_torsions = 1'+'\n')
    temp.write(' FITNESS FUNCTION SETTINGS'+'\n')
    temp.write('initial_virtual_pt_match_max = 3'+'\n')
    temp.write('relative_ligand_energy = 1'+'\n')
    temp.write('gold_fitfunc_path = plp'+'\n')
    temp.write('score_param_file = DEFAULT'+'\n')
    temp.write(' PROTEIN DATA'+'\n')
    temp.write('protein_datafile = '+receptorfilename+'\n')

    return conffile

def GenerateGOLDCommand(goldbin,filename):
    cmdstr=goldbin+' '+filename
    return cmdstr


def GOLDDocking(goldbin,ligandname,receptorname,dockgridcenter,dockgridsize,nposes):
    filename=GenerateGOLDConfigurationFile(ligandname,receptorname,dockgridcenter,dockgridsize,nposes)
    cmdstr=GenerateGOLDCommand(goldbin,filename)
    os.system(cmdstr)
    modeltoscoregold,modeltostructuregold=ParseGOLDOutput(ligandname)
    return modeltoscoregold,modeltostructuregold

def ParseGOLDOutput(ligandname):
    modeltoscoregold={}
    modeltostructuregold={}
    ligandprefix=ligandname.replace('.pdb','')
    outputname=ligandprefix+'_m1'+'.rnk'
    temp=open(outputname,'r')
    results=temp.readlines()
    temp.close()
    count=1
    for line in results:
        linesplit=line.split()
        if len(linesplit)==10:
            score=float(linesplit[1])
            model=int(linesplit[0])
            structurename='gold_soln_'+ligandprefix+'_m1_'+str(model)+'.pdb'
            modeltoscoregold[count]=score
            modeltostructuregold[count]=structurename
            count+=1

    return modeltoscoregold,modeltostructuregold


def GrabAtomPositions(pdbmol):
    atomvecls=[]
    iteratombab = openbabel.OBMolAtomIter(pdbmol)
    for atm in iteratombab:
        atmindex=atm.GetIdx()
        coords=[atm.GetX(),atm.GetY(),atm.GetZ()]
        atomvecls.append(coords)


    return atomvecls


def GrabLigandCentroid(ligandpdbfilename):
    pdbmol=GenerateMOLObject(ligandpdbfilename)
    atomvecls=GrabAtomPositions(pdbmol)
    centroid=np.array([0.0,0.0,0.0])
    for vec in atomvecls:
        centroid+=np.array(vec)
    centroid=centroid/len(atomvecls)

    
    return centroid



def CombineDockingInfo(modeltodockoutput,modeltoscore,modeltostructure,dockingprogram,modeltocomplexedstructure):
    modeltodockoutput[dockingprogram]={}
    modeltodockoutput[dockingprogram]['structure']=modeltostructure
    modeltodockoutput[dockingprogram]['complexedstructure']=modeltocomplexedstructure
    modeltodockoutput[dockingprogram]['score']=modeltoscore

    return modeltodockoutput


def FileConverter(inputfilename,informat,outformat):
    outputfilename=inputfilename.replace('.'+informat,'.'+outformat)
    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    obConversion.SetInFormat(informat)
    obConversion.ReadFile(mol, inputfilename)
    obConversion.SetOutFormat(outformat)
    obConversion.WriteFile(mol,outputfilename)
    return outputfilename


def ConvertPDBQTFilesToPDB(modeltostructure):
    newmodeltostructure={}
    for model,structure in modeltostructure.items():
        newstructure=FileConverter(structure,'pdbqt','pdb')
        newmodeltostructure[model]=newstructure

    return newmodeltostructure


def CombineMolObjects(pdbmol,ligmol,complexedstructure):
    iteratombab = openbabel.OBMolAtomIter(ligmol)
    for atm in iteratombab:
        pdbmol.AddAtom(atm)
    bonditer=openbabel.OBMolBondIter(ligmol)
    for bond in bonditer:
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        bondorder=bond.GetBondOrder()
        diditwork=pdbmol.AddBond(obgnidx,oendidx,bondorder)
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat('pdb')
    obConversion.WriteFile(pdbmol,complexedstructure)
    ConvertUNLToLIG(complexedstructure)



def GenerateComplex(structure,receptorpdbfilename):
    complexedstructure=structure.replace('.pdb','')+'_'+receptorpdbfilename.replace('.pdb','')+'.pdb'
    pdbmol=GenerateMOLObject(receptorpdbfilename)
    ligmol=GenerateMOLObject(structure)
    CombineMolObjects(pdbmol,ligmol,complexedstructure)

    return complexedstructure



def AddDockedLigandBackInComplex(modeltostructure,receptorpdbfilename):
    modeltocomplexedstructure={}
    for model,structure in modeltostructure.items():
        complexedstructure=GenerateComplex(structure,receptorpdbfilename)
        modeltocomplexedstructure[model]=complexedstructure


    return modeltocomplexedstructure

def Docking(ligandpdbfilename,receptorpdbfilename,ligandname,receptorname,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes,goldbin,usevina,usead4,usevinardo,usegold,python2path,prepdockscript,gridspacing):
    modeltodockoutput={}
   
    if usegold==True:
        modeltoscoregold,modeltostructuregold=GOLDDocking(goldbin,ligandpdbfilename,receptorpdbfilename,dockgridcenter,dockgridsize,nposes)
        modeltocomplexedstructuregold=AddDockedLigandBackInComplex(modeltostructuregold,receptorpdbfilename)
        modeltodockoutput=CombineDockingInfo(modeltodockoutput,modeltoscoregold,modeltostructuregold,'GOLD',modeltocomplexedstructuregold)
    GeneratePDBQTFiles(python2path,prepdockscript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename)
    if usead4==True:
        modeltoscoread4,modeltostructuread4=AutoDock4(python2path,prepdockscript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename)
        modeltostructuread4=ConvertPDBQTFilesToPDB(modeltostructuread4)
        modeltocomplexedstructuread4=AddDockedLigandBackInComplex(modeltostructuread4,receptorpdbfilename)
        modeltodockoutput=CombineDockingInfo(modeltodockoutput,modeltoscoread4,modeltostructuread4,'ad4',modeltocomplexedstructuread4)
    if usevina==True:
        modeltoscorevina,modeltostructurevina=AutoDockVinaVinardo('vina',receptorname,ligandname,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes)
        modeltostructurevina=ConvertPDBQTFilesToPDB(modeltostructurevina)
        modeltocomplexedstructurevina=AddDockedLigandBackInComplex(modeltostructurevina,receptorpdbfilename)
        modeltodockoutput=CombineDockingInfo(modeltodockoutput,modeltoscorevina,modeltostructurevina,'vina',modeltocomplexedstructurevina)
    if usevinardo==True:
        modeltoscorevinardo,modeltostructurevinardo=AutoDockVinaVinardo('vinardo',receptorname,ligandname,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes)
        modeltostructurevinardo=ConvertPDBQTFilesToPDB(modeltostructurevinardo)
        modeltocomplexedstructurevinardo=AddDockedLigandBackInComplex(modeltostructurevinardo,receptorpdbfilename)
        modeltodockoutput=CombineDockingInfo(modeltodockoutput,modeltoscorevinardo,modeltostructurevinardo,'vinardo',modeltocomplexedstructurevinardo)
    GenerateDockingReport(modeltodockoutput)

def GenerateDockingReport(modeltodockoutput):
    name='DockingReport.txt'
    temp=open(name,'w')
    for model,dic in modeltodockoutput.items():
        temp.write('Docking model : '+model+'\n')
        for keyword,innerdic in dic.items():
            temp.write('Dictionary '+keyword+'\n')
            for model,value in innerdic.items():
                temp.write('Ranking '+str(model)+' '+str(value)+'\n')
           

    temp.close()

def DockingWrapper(ligandreceptorfilename,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes,goldbin,usevina,usead4,usevinardo,usegold,dockingenvpath,prepdockscript,gridspacing):
    python2path=os.path.join(dockingenvpath,'bin')
    python2path=os.path.join(python2path,'python')
    ligandpdbfilename,receptorpdbfilename=ExtractLigand(ligandreceptorfilename)
    if dockgridcenter==None:
        dockgridcenter=GrabLigandCentroid(ligandpdbfilename)
    receptorname=receptorpdbfilename.replace('.pdb','.pdbqt')
    ligandname=ligandpdbfilename.replace('.pdb','.pdbqt')
    Docking(ligandpdbfilename,receptorpdbfilename,ligandname,receptorname,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes,goldbin,usevina,usead4,usevinardo,usegold,python2path,prepdockscript,gridspacing)


