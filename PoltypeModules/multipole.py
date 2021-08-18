import electrostaticpotential as esp
import time
import os
import sys
import openbabel
import shutil
import re
from collections import deque


def SanitizeMultipoleFrames(poltype,keyfilename): # pearl script for averging only understands 0 and not empty spaces for parsing
    temp=open(keyfilename,'r')
    results=temp.readlines()
    temp.close()
    tempname=keyfilename.replace('.key','_temp.key')
    temp=open(tempname,'w')
    extraspace='     '
    for line in results:
        if 'multipole' in line:
            linesplit=line.split()
            realsplit=re.split(r'(\s+)', line)
            if len(linesplit)==3:
                realsplit=realsplit[:3]+[extraspace]+['0']+[extraspace]+['0']+realsplit[3:]
            elif len(linesplit)==4:
                realsplit=realsplit[:5]+[extraspace]+['0']+realsplit[5:]
            line=''.join(realsplit)
        temp.write(line)
    temp.close()
    os.remove(keyfilename)
    os.rename(tempname,keyfilename)


def gen_peditinfile(poltype,mol):
    # write out the local frames
    iteratom = openbabel.OBMolAtomIter(mol)
    f = open (poltype.peditinfile, 'w')
    f.write("\n")
    f.write('A'+'\n')

    #Find aromatic carbon, halogens, and bonded hydrogens to correct polarizability
    iteratom = openbabel.OBMolAtomIter(mol)
    writesection=False
    lines=[]
    for a in iteratom:
        if (a.GetAtomicNum() == 6 and a.IsAromatic()):
            lines.append(str(a.GetIdx()) + " " + str(1.750) + "\n")
            writesection=True
        elif (a.GetAtomicNum() == 1):
            iteratomatom = openbabel.OBAtomAtomIter(a)
            for b in iteratomatom:
                if (b.GetAtomicNum() == 6 and b.IsAromatic()):
                    lines.append(str(a.GetIdx()) + " " + str(0.696) + "\n")
                    writesection=True
    if writesection:
        for line in lines:
            f.write(line)

    


    f.write("\n")
    f.flush()
    os.fsync(f.fileno())
    f.write("2\n")
    f.write("N\n")
    f.write("Y\n")


    f.flush()
    os.fsync(f.fileno())

    f.close()
    


def rm_esp_terms_keyfile(poltype,keyfilename):
    """
    Intent: Remove unnecessary terms from the key file
    """
    tmpfname = keyfilename + "_tmp"
    temp=open(tmpfname,'w')
    anothertemp=open(keyfilename,'r')
    results=anothertemp.readlines()
    anothertemp.close()
    passedatomblock=False
    foundatomblock=False
    foundpolarize=False
    for line in results:
        if ('potential-offset' in line or 'none' in line or 'fix-monopole' in line) and 'atom' not in line:
            pass
        else:
            linesplit=line.split()
            if len(linesplit)>0:
                if linesplit[0]=='atom':
                    if foundatomblock==False:
                        foundatomblock=True
                    temp.write(line)

                else:
                    if foundatomblock==True:
                        passedatomblock=True
                    if linesplit[0]=='polarize':
                        foundpolarize=True
                    if passedatomblock==True and foundpolarize==False: # then old multipoles
                        pass
                    else:
                        temp.write(line)
            else:
                temp.write(line)

    temp.close()
    shutil.move(tmpfname, keyfilename)
   
    
def prepend_keyfile(poltype,keyfilename,optmol,dipole=False):
    """
    Intent: Adds a header to the key file given by 'keyfilename'
    """
    while not os.path.isfile(keyfilename):
        time.sleep(5)
        poltype.WriteToLog('Waiting for '+keyfilename)
    tmpfname = keyfilename + "_tmp"
    tmpfh = open(tmpfname, "w")
    keyfh = open(keyfilename, "r")
    tmpfh.write("parameters " + poltype.paramhead + "\n")
    tmpfh.write('OPENMP-THREADS '+str(poltype.numproc)+'\n')
    tmpfh.write("bondterm none\n")
    tmpfh.write("angleterm none\n")
    tmpfh.write("torsionterm none\n")
    tmpfh.write("vdwterm none\n")
    tmpfh.write("fix-monopole\n")
    tmpfh.write("digits 8\n")
    tmpfh.write("potential-offset 1.0\n\n")
    logname=poltype.logespfname
    if poltype.fitqmdipole==True and os.path.exists(logname):
        qmdipole=esp.GrabQMDipoles(poltype,optmol,logname)
        tmpfh.write('TARGET-DIPOLE'+' '+str(qmdipole[0])+' '+str(qmdipole[1])+' '+str(qmdipole[2])+'\n')

    for line in keyfh:
        tmpfh.write(line)
    shutil.move(tmpfname, keyfilename)


def CheckFileForString(poltype,filetocheck):
    temp=open(filetocheck,'r')
    results=temp.readlines()
    temp.close()
    usemp2=False
    for line in results:
        if 'MP2' in line and 'Density' in line:
            usemp2=True
    if usemp2==False:
        string='CC'
    else:
        string='MP2'
    return string


def gen_gdmain(poltype,gdmainfname,molecprefix,fname,dmamethod):
    """
    Intent: Generate GDMA input file for the molecule
    Input:
        gdmainfname: file name for the gdma input file
        molecprefix: name of the molecule
        fname: file name for the *-dma.fchk file
    Output: *.gdmain file is created
    Referenced By: run_gdma
    Description:
        1. create pointer file dma.fchk to *-dma.fchk
        2. create *.gdmain file and write in all the necessary information for the gdma run
    """
    fnamesym = "dma.fchk"
    try:
        os.symlink(fname,fnamesym)
    except Exception as e:
        print(e)



    #punfname = os.path.splitext(fname)[0] + ".punch"
    punfname = "dma.punch"

    try:
        tmpfh = open(gdmainfname, "w")
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))
        sys.exit(e.errno)

    tmpfh.write("Title " + molecprefix + " gdmain\n")
    tmpfh.write("\n")
    if poltype.dmamethod=='MP2':
        densitystring='MP2'
    else:
        densitystring='SCF'
    if poltype.dmamethod=='MP2':
        densitystring=CheckFileForString(poltype,fnamesym)
    tmpfh.write("File " + fnamesym  + " density %s\n"%(densitystring))
    tmpfh.write("Angstrom\n")
    tmpfh.write("AU\n")
    tmpfh.write("Multipoles\n")
    tmpfh.write("Switch 0\n") # v1.3, comment out for v2.2
    tmpfh.write("Limit 2\n")
    tmpfh.write("Punch " + punfname + "\n")
    tmpfh.write("Radius H 0.65\n")
    tmpfh.write("Radius S 0.80\n")
    tmpfh.write("Radius P 0.75\n")
    tmpfh.write("\n")
    tmpfh.write("Start\n")
    tmpfh.write("\n")
    tmpfh.write("Finish\n")
    tmpfh.close()

def run_gdma(poltype):
    """
    Intent: Runs GDMA to find multipole information
    The GDMA program carries out distributed multipole analysis of the wavefunctions
    calculated by Gaussian (using the *-dma.fchk file)
    Input:
    Output: *.gdmaout is created; this is used as input by poledit, a tinker program
    Referenced By: main
    Description:
    1. Generates the gdma input file by calling 'gen_gdmain'
    2. Runs the following command: gdma < *.gdmain > *.gdmaout
    """
    poltype.WriteToLog("NEED DMA: Executing GDMA")

    if not os.path.isfile(poltype.fckdmafname):
        poltype.fckdmafname = os.path.splitext(poltype.fckdmafname)[0]

    try:
        assert os.path.isfile(poltype.fckdmafname), "Error: " + poltype.fckdmafname + " does not exist."+' '+os.getcwd()
    except:
        poltype.DeleteFilesWithString(['dma'])
        poltype.GenerateParameters()

    poltype.gdmainfname = poltype.assign_filenames ( "gdmainfname" , ".gdmain")
    gen_gdmain(poltype,poltype.gdmainfname,poltype.molecprefix,poltype.fckdmafname,poltype.dmamethod)

    cmdstr = poltype.gdmaexe + " < " + poltype.gdmainfname + " > " + poltype.gdmafname
    poltype.call_subsystem(cmdstr,True)

    assert os.path.getsize(poltype.gdmafname) > 0, "Error: " + os.getcwd() +' '+os.path.basename(poltype.gdmaexe) + " cannot create .gdmaout file."
   
def AverageMultipoles(poltype,optmol):
    # gen input file
    gen_avgmpole_groups_file(poltype)
    # call avgmpoles.pl
    avgmpolecmdstr = poltype.avgmpolesexe + " " + poltype.keyfname + " " + poltype.xyzfname + " " + poltype.grpfname + " " + poltype.key2fname + " " + poltype.xyzoutfile + " " + str(poltype.prmstartidx)
    poltype.call_subsystem(avgmpolecmdstr,True)
    prepend_keyfile(poltype,poltype.key2fname,optmol,True)

def gen_avgmpole_groups_file(poltype):
    """
    Intent: Print out *-groups.txt which is a map from symm class to idx
    Also, symm class labels are altered from 1, 2, 3, ... to 401, 402, ...
    (or some other labeling system dependent on 'prmstartidx')
    Input:
    Output:
        *-groups.txt: map from symm group to atom idx
    Referenced By: main
    Description: -
    """
    symclasstoidxlist={}
    for idx,symclass in poltype.idxtosymclass.items():
        if symclass not in symclasstoidxlist.keys():
            symclasstoidxlist[symclass]=[]
        symclasstoidxlist[symclass].append(idx)
    f = open(poltype.grpfname,"w")
    for symclass,idxlist in symclasstoidxlist.items():
        string=str(symclass)+' '
        for idx in idxlist:
            string+=str(idx)+' '
        string+='\n'
        f.write(string)
    f.close()
