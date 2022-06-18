import os
import sys
from pymol import cmd
import json
import getopt


def CallFPocket(bindir,pdbfile):
    binpath=os.path.join(bindir,'fpocket')
    descname='pocket_descriptors.csv'
    outputfolder=pdbfile.replace('.pdb','_out')
    cmdstr=binpath+' '+'-f'+' '+pdbfile+' '+'-d'+' '+'>'+' '+descname
    os.system(cmdstr)
    return descname,outputfolder


def GrabPocketGrid(outputfolder):
    extending=5
    pocketnumtocenter={}
    pocketnumtosize={}
    for file in os.listdir(outputfolder):
        if 'pqr' in file:
            pocket_num=int(file.split('_')[0].replace('pocket',''))
            cmd.load(filename=os.path.join(outputfolder,file),format='pqr',object=pocket_num)
            ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(pocket_num) 
            minX = minX - float(extending)
            minY = minY - float(extending)
            minZ = minZ - float(extending)
            maxX = maxX + float(extending)
            maxY = maxY + float(extending)
            maxZ = maxZ + float(extending)
            SizeX = maxX - minX
            SizeY = maxY - minY
            SizeZ = maxZ - minZ
            CenterX =  (maxX + minX)/2
            CenterY =  (maxY + minY)/2
            CenterZ =  (maxZ + minZ)/2
            center=[CenterX,CenterY,CenterZ]
            size=[SizeX,SizeY,SizeZ]
            pocketnumtocenter[pocket_num]=center
            pocketnumtosize[pocket_num]=size

    return pocketnumtocenter,pocketnumtosize



def WriteDictionaryToFile(dic,filename):
    with open(filename, 'w') as convert_file:
        convert_file.write(json.dumps(dic))




opt_list, args = getopt.getopt(sys.argv[1:], 'p:pdb:',['bindir=','pdbfile='])
for o, a in opt_list:
    if o in ('--bindir'):
        bindir=a
    if o in ('--pdbfile'):
        pdbfile=a


descname,outputfolder=CallFPocket(bindir,pdbfile)
pocketnumtocenter,pocketnumtosize=GrabPocketGrid(outputfolder)
WriteDictionaryToFile(pocketnumtocenter,'pocketnumtocenter.json')
WriteDictionaryToFile(pocketnumtosize,'pocketnumtosize.json')


