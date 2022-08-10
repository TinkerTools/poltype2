import getopt
import sys
import numpy as np
import itertools
from PIL import Image

def Chunks(lst, n):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def ChunksList(gen):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newlst=[]
    for item in gen:
        newlst.append(item)
    return newlst


def FindDimensionsOfMolecule(xyzfile):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    veclist=[]
    temp=open(xyzfile,'r')
    results=temp.readlines()
    temp.close()
    count=0
    for line in results:
        linesplit=line.split()
        if len(linesplit)==1:
            atomnum=int(linesplit[0])
            if 'H20' in xyzfile:
                maxatomnum=atomnum-3
            else:
                maxatomnum=int(atomnum/2)
        if len(linesplit)>1: # not line containing number of atoms
            if count>=maxatomnum:
                break
            vec=np.array([float(linesplit[1]),float(linesplit[2]),float(linesplit[3])])
            veclist.append(vec)
            count+=1

    pairs=list(itertools.combinations(veclist, 2))
    disttodiffvec={}
    
    for pairidx in range(len(pairs)):
        pair=pairs[pairidx]
        progress=(pairidx*100)/len(pairs)
        diff=np.array(pair[0])-np.array(pair[1])
        dist=np.linalg.norm(diff)
        
        disttodiffvec[dist]=diff
    distlist=list(disttodiffvec.keys())
    if len(distlist)!=0:
        mindist=np.amax(np.array(distlist))
        diffvec=disttodiffvec[mindist]
    else:
        mindist=0
    return mindist


def SmallestDivisor(n):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    a=[]
    for i in range(2,n+1):
        if(n%i==0):
            a.append(i)
    a.sort()
    return a[0]




def PlotDimers3D(filenamearray,allindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    from pymol import cmd,preset,util
    from pymol.vfont import plain
    from pymol.cgo import CYLINDER,cyl_text

    molsPerImage=len(filenamearray)
    if (molsPerImage % 2) == 0 or (molsPerImage ** 0.5) % 1==0:
        n=molsPerImage
    else:
        n=molsPerImage+1 
    molsperrow=SmallestDivisor(n)
    if n==2:
        molsperrow=1
    imagesize=1180
    dpi=300
    size=FindDimensionsOfMolecule(filenamearray[0])
    filenamechunks=ChunksList(Chunks(filenamearray,molsPerImage))
    indiceschunks=ChunksList(Chunks(allindices,molsPerImage))
    prevmatslen=len(filenamechunks[0])
    for i in range(len(filenamechunks)):
        filenamesublist=filenamechunks[i]
        indicessublist=indiceschunks[i]
        imagenames=[]
        for j in range(len(filenamesublist)):
            filename=filenamesublist[j]
            indices=indicessublist[j]
            ls=range(len(filenamesublist))
            chunks=ChunksList(Chunks(ls,molsperrow))
            indextorow={}
            for rowidx in range(len(chunks)):
                row=chunks[rowidx]
                for j in row:
                    indextorow[j]=rowidx
            
            fileprefix=filename.split('.')[0]
            imagename=fileprefix+'_3D.'+'png'
            imagenames.append(imagename)
            cmd.delete('all')
            cmd.load(filename)
            preset.ball_and_stick(selection='all', mode=1)
            cmd.bg_color("white")
            cmd.color("grey50","all")
            util.cnc("all")
            cmd.set('label_size',26) 
            cmd.set('depth_cue',0)
            cmd.set('ray_trace_fog',0) 
            firstidx=indices[0]
            secondidx=indices[1]
            firstlab=str(firstidx)
            secondlab=str(secondidx)
            lab=firstlab+'-'+secondlab
            cmd.distance(lab,'index '+firstlab,'index '+secondlab)
            cmd.zoom(lab,size,0,0)
            cmd.png(imagename, imagesize,imagesize,dpi,1)
            cmd.save(fileprefix+'_3D.'+'pse')
        if i>0:
            factor=1
        else:
            factor=0
        firstj=i*prevmatslen+factor
        secondj=firstj+len(filenamesublist)-factor
        prevmatslen=len(filenamesublist)

        basename=fileprefix+'_'+str(firstj)+'-'+str(secondj)
        indextoimage={}
        for index in range(len(filenamesublist)):
            imagename=imagenames[index]
            image=Image.open(imagename)
            indextoimage[index]=image
        cols=len(set(list(indextorow.values())))
        dest = Image.new('RGB', (imagesize*molsperrow,imagesize*cols))
        for j in range(len(filenamesublist)):
            row=indextorow[j]
            x=(j-molsperrow*(row))*imagesize
            y=(row)*imagesize
            dest.paste(indextoimage[j],(x,y))
        dest.show()
        dest.save(basename+'.png')



def PlotESPSurfaces(name,cubefiles):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    from pymol import cmd,preset,util
    from pymol.vfont import plain
    from pymol.cgo import CYLINDER,cyl_text

    imagesize=1180
    dpi=300
    imagename=name+'_ESP.'+'png'
    cmd.delete('all')
    for filename in cubefiles:
        cmd.load(filename)
    preset.ball_and_stick(selection='all', mode=1)
    cmd.bg_color("white")
    cmd.color("grey50","all")
    util.cnc("all")
    cmd.set('label_size',26) 
    cmd.set('depth_cue',0)
    cmd.set('ray_trace_fog',0) 
    cmd.isosurface('Dt2','Dt', 0.001)
    cmd.ramp_new('espcol', 'ESP', [-.01,-.005,0,.005,.01], ['red','orange', 'yellow','green', 'blue'])
    cmd.set('surface_color', 'espcol', 'Dt2')
    cmd.zoom()
    cmd.disable('espcol')
    cmd.png(imagename, imagesize,imagesize,dpi,1)
    cmd.save(name+'_ESP.'+'pse')



filenames=None
allindices=None
name=None
cubefiles=None
opt_list, args = getopt.getopt(sys.argv[1:], 'dgc:dgz:s:ln:rn:',['filenames=','allindices=','cubefiles=','name='])
for o, a in opt_list:
    if o in ('--cubefiles'):
        cubefiles=a.split(',')
    if o in ('--name'):
        name=a
    if o in ('--filenames'):
        filenames=a.split(',')
    if o in ('--allindices'):
        allindices=a.split(',')
        templist=[]
        for ele in allindices:
            nums=ele.lstrip().rstrip().split('_')
            temp=[]
            for e in nums:
                temp.append(int(e))
            templist.append(temp)
        allindices=templist
if filenames!=None and allindices!=None:
    PlotDimers3D(filenames,allindices)

if name!=None and cubefiles!=None:
    PlotESPSurfaces(name,cubefiles)
