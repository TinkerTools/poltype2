from openbabel import openbabel
import re
import os
import numpy as np
from pathlib import Path
import sys








def GenerateListString(ls):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    ls=[str(i) for i in ls]
    string=','.join(ls)
    return string

def GenerateCommandString(python2path,preparescript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    cmdstr=''
    dockinggridstring=GenerateListString(dockgridcenter)
    dockinggridsizestring=GenerateListString(dockgridsize)
    cmdstr+=python2path+' '+preparescript+' '+'--dockgridcenter='+dockinggridstring+' '+'--dockgridsize='+dockinggridsizestring+' '+'--spacing='+str(gridspacing)+' '+'--ligandname='+ligandpdbfilename+' '+'--receptorname='+receptorpdbfilename
    return cmdstr


def GenerateGridFile(receptorgridname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    cmdstr='autogrid4 '+'-p '+receptorgridname
    os.system(cmdstr)


def RunAutoDock4(ad4parameterfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    cmdstr='autodock4 '+'-p '+ad4parameterfilename
    os.system(cmdstr)

def GeneratePDBQTFiles(python2path,prepdockscript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """


    cmdstr=GenerateCommandString(python2path,prepdockscript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename)

    os.system(cmdstr)



def AutoDock4(python2path,prepdockscript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

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
    if poltype.covalentdocking==False:
        temp.write('origin = '+str(dockgridcenter[0])+' '+str(dockgridcenter[1])+' '+str(dockgridcenter[2])+'\n')
    else:
        temp.write('floodfill_atom_no = '+str(poltype.prolinkatom)+'\n')
        temp.write('floodfill_center = atom'+'\n')
    temp.write('do_cavity = 1'+'\n')

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
    if poltype.covalentdock==False:
        temp.write(' COVALENT BONDING'+'\n')
        temp.write('covalent = 0'+'\n')
    else:
        temp.write('covalent = 1'+'\n')
        temp.write('covalent_protein_atom_no = '+str(poltype.prolinkatom)+'\n')
        temp.write('covalent_ligand_atom_no = '+str(poltype.liglinkatom)+'\n')


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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    cmdstr=goldbin+' '+filename
    return cmdstr


def GOLDDocking(goldbin,ligandname,receptorname,dockgridcenter,dockgridsize,nposes):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    filename=GenerateGOLDConfigurationFile(ligandname,receptorname,dockgridcenter,dockgridsize,nposes)
    cmdstr=GenerateGOLDCommand(goldbin,filename)
    os.system(cmdstr)
    modeltoscoregold,modeltostructuregold=ParseGOLDOutput(ligandname)
    return modeltoscoregold,modeltostructuregold

def ParseGOLDOutput(ligandname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

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






def CombineDockingInfo(modeltodockoutput,modeltoscore,modeltostructure,dockingprogram,modeltocomplexedstructure,ligandpdbfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    if dockingprogram not in modeltodockoutput.keys():
        modeltodockoutput[dockingprogram]={}
    if ligandpdbfilename not in modeltodockoutput[dockingprogram].keys():
        modeltodockoutput[dockingprogram][ligandpdbfilename]={}
    if 'structure' not in modeltodockoutput[dockingprogram][ligandpdbfilename].keys():
        modeltodockoutput[dockingprogram][ligandpdbfilename]['structure']={}
    if 'complexedstructure' not in modeltodockoutput[dockingprogram][ligandpdbfilename].keys():
        modeltodockoutput[dockingprogram][ligandpdbfilename]['complexedstructure']={}
    if 'score' not in modeltodockoutput[dockingprogram][ligandpdbfilename].keys():
        modeltodockoutput[dockingprogram][ligandpdbfilename]['score']={}
    modeltodockoutput[dockingprogram][ligandpdbfilename]['structure'].update(modeltostructure)
    modeltodockoutput[dockingprogram][ligandpdbfilename]['complexedstructure'].update(modeltocomplexedstructure)
    modeltodockoutput[dockingprogram][ligandpdbfilename]['score'].update(modeltoscore)

    return modeltodockoutput


def FileConverter(inputfilename,informat,outformat):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    outputfilename=inputfilename.replace('.'+informat,'.'+outformat)
    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    obConversion.SetInFormat(informat)
    obConversion.ReadFile(mol, inputfilename)
    obConversion.SetOutFormat(outformat)
    obConversion.WriteFile(mol,outputfilename)
    return outputfilename


def ConvertPDBQTFilesToPDB(modeltostructure):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    newmodeltostructure={}
    for model,structure in modeltostructure.items():
        newstructure=FileConverter(structure,'pdbqt','pdb')
        newmodeltostructure[model]=newstructure

    return newmodeltostructure


def CombineMolObjects(pdbmol,ligmol,complexedstructure):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

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
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    complexedstructure=structure.replace('.pdb','')+'_'+receptorpdbfilename.replace('.pdb','')+'.pdb'
    pdbmol=poltype.GenerateMOLObject(receptorpdbfilename)
    ligmol=poltype.GenerateMOLObject(structure)
    CombineMolObjects(pdbmol,ligmol,complexedstructure)

    return complexedstructure



def AddDockedLigandBackInComplex(modeltostructure,receptorpdbfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    modeltocomplexedstructure={}
    for model,structure in modeltostructure.items():
        complexedstructure=GenerateComplex(structure,receptorpdbfilename)
        modeltocomplexedstructure[model]=complexedstructure


    return modeltocomplexedstructure

def Docking(ligandpdbfilename,receptorpdbfilename,ligandname,receptorname,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes,goldbin,usevina,usead4,usevinardo,usegold,python2path,prepdockscript,gridspacing,modeltodockoutput):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    if usegold==True:
        modeltoscoregold,modeltostructuregold=GOLDDocking(goldbin,ligandpdbfilename,receptorpdbfilename,dockgridcenter,dockgridsize,nposes)
        modeltocomplexedstructuregold=AddDockedLigandBackInComplex(modeltostructuregold,receptorpdbfilename)
        modeltodockoutput=CombineDockingInfo(modeltodockoutput,modeltoscoregold,modeltostructuregold,'GOLD',modeltocomplexedstructuregold,ligandpdbfilename)
    GeneratePDBQTFiles(python2path,prepdockscript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename)
    if usead4==True:
        modeltoscoread4,modeltostructuread4=AutoDock4(python2path,prepdockscript,dockgridcenter,dockgridsize,gridspacing,ligandpdbfilename,receptorpdbfilename)
        modeltostructuread4=ConvertPDBQTFilesToPDB(modeltostructuread4)
        modeltocomplexedstructuread4=AddDockedLigandBackInComplex(modeltostructuread4,receptorpdbfilename)
        modeltodockoutput=CombineDockingInfo(modeltodockoutput,modeltoscoread4,modeltostructuread4,'ad4',modeltocomplexedstructuread4,ligandpdbfilename)
    if usevina==True:
        modeltoscorevina,modeltostructurevina=AutoDockVinaVinardo('vina',receptorname,ligandname,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes)
        modeltostructurevina=ConvertPDBQTFilesToPDB(modeltostructurevina)
        modeltocomplexedstructurevina=AddDockedLigandBackInComplex(modeltostructurevina,receptorpdbfilename)
        modeltodockoutput=CombineDockingInfo(modeltodockoutput,modeltoscorevina,modeltostructurevina,'vina',modeltocomplexedstructurevina,ligandpdbfilename)
    if usevinardo==True:
        modeltoscorevinardo,modeltostructurevinardo=AutoDockVinaVinardo('vinardo',receptorname,ligandname,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes)
        modeltostructurevinardo=ConvertPDBQTFilesToPDB(modeltostructurevinardo)
        modeltocomplexedstructurevinardo=AddDockedLigandBackInComplex(modeltostructurevinardo,receptorpdbfilename)
        modeltodockoutput=CombineDockingInfo(modeltodockoutput,modeltoscorevinardo,modeltostructurevinardo,'vinardo',modeltocomplexedstructurevinardo,ligandpdbfilename)
    return modeltodockoutput

def GenerateDockingReport(modeltodockoutput):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    name='DockingReport.txt'
    temp=open(name,'w')
    for model,dic in modeltodockoutput.items():
        temp.write('Docking model : '+model+'\n')
        for lig,nextdic in dic.items():
            temp.write('Ligand : '+lig+'\n')
            for keyword,innerdic in nextdic.items():
                temp.write('Dictionary '+keyword+'\n')
                for model,value in innerdic.items():
                    temp.write('Pose Ranking '+str(model)+' '+str(value)+'\n')
           

    temp.close()

def DetermineScoreTopPoses(modeltodockoutput):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    modeltoligandtotopposescore={}
    for model,dic in modeltodockoutput.items():
        for lig,nextdic in dic.items():
            for keyword,innerdic in nextdic.items():
                if keyword=='score':
                    topscore=innerdic[1]
                    if model not in modeltoligandtotopposescore.keys():
                        modeltoligandtotopposescore[model]={}
                    if lig not in modeltoligandtotopposescore[model].keys():
                        modeltoligandtotopposescore[model][lig]=topscore

    return modeltoligandtotopposescore


def DetermineRanking(modeltoligandtotopposescore):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    modeltoligandtorank={}
    for model,ligandtotopposescore in modeltoligandtotopposescore.items():
        newligandtotopposescore={} # sometimes scoring function is negative sometimes positive so keep all positive
        for ligand,score in ligandtotopposescore.items():
            if score<0:
                score=score*-1
            newligandtotopposescore[ligand]=score
        sortedligtoscore={k: v for k, v in sorted(newligandtotopposescore.items(), key=lambda item: item[1],reverse=True)}
        count=1
        for sortedlig,score in sortedligtoscore.items():
            if model not in modeltoligandtorank.keys():
                modeltoligandtorank[model]={}

            modeltoligandtorank[model][sortedlig]=count
            count+=1
    return modeltoligandtorank


def ECRScore(rank,ecrexpect,dockingprogram,moleculename,nametomodeltoecrscore):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    ecrscore=(1/ecrexpect)*np.exp(-rank/ecrexpect)
    if moleculename not in nametomodeltoecrscore.keys():
        nametomodeltoecrscore[moleculename]={}
    nametomodeltoecrscore[moleculename][dockingprogram]=ecrscore
    return nametomodeltoecrscore


def ExponentialConsensusRank(modeltodockoutput,ecrexpect):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    

    nametoconsensusrank={}
    modeltoligandtotopposescore=DetermineScoreTopPoses(modeltodockoutput)
    modeltoligandtorank=DetermineRanking(modeltoligandtotopposescore)
    nametomodeltoecrscore={}
    for model,ligandtorank in modeltoligandtorank.items():
        for ligand,rank in ligandtorank.items():
            nametomodeltoecrscore=ECRScore(rank,ecrexpect,model,ligand,nametomodeltoecrscore) 
    for name,modeltoscore in nametomodeltoecrscore.items():
        Sum=0
        for model,score in modeltoscore.items():
            Sum+=score
        nametoconsensusrank[name]=Sum
    sortednametoconsensusrank={k: v for k, v in sorted(nametoconsensusrank.items(), key=lambda item: item[1],reverse=True)}
    newnametoconsensusrank={}
    count=1
    for name,score in sortednametoconsensusrank.items():
        newnametoconsensusrank[name]=count
        count+=1

    return newnametoconsensusrank


def GenerateConsensusDockingReport(nametoconsensusrank):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    temp=open('ConsensusDockingReport.txt','w')
    for name,consensusrank in nametoconsensusrank.items():
        temp.write('Ligand Name : '+name+' , Consensus Score : '+str(consensusrank)+'\n')
    temp.close()



def DockingWrapper(poltype,ligandreceptorfilename,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes,goldbin,usevina,usead4,usevinardo,usegold,dockingenvname,prepdockscript,gridspacing,listofligands,ecrexpect):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    pythonpath=poltype.which('python')
    head,tail=os.path.split(pythonpath)
    pythonpath=Path(head) 
    envdir=pythonpath.parent.absolute()
    envpath=Path(envdir)
    allenvs=envpath.parent.absolute()
    dockingenvpath=os.path.join(allenvs,dockingenvname)
    python2path=os.path.join(dockingenvpath,'bin')
    python2path=os.path.join(python2path,'python')
    ligandpdbfilename,receptorpdbfilename=ExtractLigand(ligandreceptorfilename)
    modeltodockoutput={}
    if dockgridcenter==None:
        dockgridcenter=poltype.GrabLigandCentroid(ligandpdbfilename)
    for ligand in listofligands:
        receptorname=receptorpdbfilename.replace('.pdb','.pdbqt')
        ligandname=ligandpdbfilename.replace('.pdb','.pdbqt')
        ligandsplit=ligand.split('.')
        ligandprefix=ligandsplit[:-1]
        ligandsuffix=ligandsplit[-1]
        if ligandsuffix!='pdb':
            ligandpdbfilename=FileConverter(ligand,ligandsuffix,'pdb')
            ligandname=ligandpdbfilename.replace('.pdb','.pdbqt')
        modeltodockoutput=Docking(ligandpdbfilename,receptorpdbfilename,ligandname,receptorname,dockgridcenter,dockgridsize,vinaexhaustiveness,nposes,goldbin,usevina,usead4,usevinardo,usegold,python2path,prepdockscript,gridspacing,modeltodockoutput)

    GenerateDockingReport(modeltodockoutput)
    nametoconsensusrank=ExponentialConsensusRank(modeltodockoutput,ecrexpect)
    GenerateConsensusDockingReport(nametoconsensusrank)
