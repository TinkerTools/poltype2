import unittest
import os,sys,inspect
from poltypepackage import poltype
import shutil
from parmed.tinker import parameterfile
from scipy.optimize import fmin
import argparse
import getopt

class TestPoltype(unittest.TestCase):

    def __init__(self, method, examplepath):
         super(TestPoltype, self).__init__(method)
         self.examplepath = examplepath

    def MakeTestCasePath(self,testfoldername):

        self.poltypemodulepath = os.path.abspath(os.path.join(__file__, os.pardir))
        self.poltypepath=os.path.abspath(os.path.join(self.poltypemodulepath, os.pardir))+r'/'
        self.testcaseparentpath=self.poltypepath+'TestCases'
        if os.path.isdir(self.testcaseparentpath):
            shutil.rmtree(self.testcaseparentpath)
            os.mkdir(self.testcaseparentpath)
        os.chdir(self.testcaseparentpath)
        testcasepath=self.testcaseparentpath+'/'+testfoldername
        if not os.path.isdir(testcasepath):
            os.makedirs(testcasepath)
        os.chdir(testcasepath)
        return testcasepath

    def GenericCopy(self,testcasepath,examplefolder,examplestructure):
        exampleparentfolder='Examples/'
        examplekey5fname=examplestructure.replace('.sdf','.key_5')
        examplestructurepath=self.poltypepath+exampleparentfolder+examplefolder+examplestructure
        examplekeyfilepath=self.ptypepath+exampleparentfolder+examplefolder+examplekey5fname
        shutil.copy(examplestructurepath,testcasepath+examplestructure)
        return examplekeyfilepath

    def GenericTest(self,exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath):
        testparams=ptype.main()
        exampleparams = parameterfile.AmoebaParameterSet(examplekeyfilepath)
        for key in testparams.bonds.keys():
            testitem=testparams.bonds[key]
            exampleitem=exampleparams.bonds[key]
            testk=testitem.k
            examplek=exampleitem.k
            self.assertEqual(testk,examplek)
            testreq=testitem.req
            examplereq=exampleitem.req
            self.assertEqual(testreq,examplereq)
        for key in testparams.angles.keys():     
            testitem=testparams.angles[key]
            exampleitem=exampleparams.angles[key]
            testk=testitem.k
            examplek=exampleitem.k
            self.assertEqual(testk,examplek)
            testtheteq=testitem.theteq
            exampletheteq=exampleitem.theteq
            self.assertEqual(testtheteq,exampletheteq)
        for key in testparams.stretch_bends.keys():
            testitem=testparams.stretch_bends[key]
            exampleitem=exampleparams.stretch_bends[key]
            testk1=testitem.k1
            examplek1=exampleitem.k1
            self.assertEqual(testk1,examplek1)
            testk2=testitem.k2
            examplek2=exampleitem.k2
            self.assertEqual(testk2,examplek2)

        for key in testparams.dihedrals.keys():
            testitem=testparams.dihedrals[key]
            exampleitem=exampleparams.dihedrals[key]
            testk=testitem.k
            examplek=exampleitem.k
            self.assertEqual(testk,examplek)
            testperiodicity=testitem.periodicity
            exampleperiodicity=exampleitem.periodicity
            self.assertEqual(testperiodicity,exampleperiodicity)
            testphase=testitem.phase
            examplephase=exampleitem.phase
            self.assertEqual(testphase,examplephase)

        for key in testparams.multipoles.keys():
            testitem=testparams.multipoles[key]
            exampleitem=exampleparams.multipoles[key]
            testmpoles=testitem.potential_terms
            examplempoles=exampleitem.potential_terms
            for termidx in range(len(testmpoles)):
                testterm=testmpoles[termidx]
                exampleterm=examplempoles[termidx]
                self.assertEqual(testterm,exampleterm)           

        self.assertEqual(testparams.multipoles,exampleparams.multipoles)
        try: # parmED needs to be pushed to new branch, they had a bug, use exception until main branch has merged changes
            for atomkey in testparams.atoms.keys():
                testatom=testparams.atoms[atomkey]
                exampleatom=exampleparams.atoms[atomkey]
                testvdwr=testatom.size
                examplevdwr=exampleatom.size
                self.assertEqual(testvdwr,examplevdwr)
                testvdwe=testatom.epsilon
                exampleatome=exampleatom.epsilon
                self.assertEqual(testvdwe,exampleatome)
                testpol=testatom.polarizability
                examplepol=exampleatom.polarizability
                self.assertEqual(testpol,examplepol)
                testpolconntypes=testatom.connected_types
                examplepolconntypes=exampleatom.connected_types
                self.assertEqual(testpolconntypes,examplepolconntypes)
        except:
            pass

    def GenericFolderCopy(self,testcasepath):
        head,tail=os.path.split(testcasepath) # remove last folder since it needs to be copied from examples
        os.chdir('..')
        shutil.rmtree(tail)
        exampleparentfolder='Examples/'
        if self.examplepath==None:
            raise ValueError('examplepath not defined')
        head,examplefoldername=os.path.split(self.examplepath)
        examplefoldername='Fragmenter/'+examplefoldername+r'/'
        newexamplepath=self.poltypepath+exampleparentfolder+examplefoldername
        if not os.path.exists(newexamplepath):
            shutil.copytree(examplepath,newexamplepath)
        shutil.copytree(examplepath,testcasepath)
        testcaseqmtorsionpath=testcasepath+r'/'+'qm-torsion'
        if os.path.exists(testcaseqmtorsionpath):
            shutil.rmtree(testcaseqmtorsionpath)
   
    def GrabTorsions(self,folderpath):
        listoftorsions=[]
        listofqmminusmmtxtfiles=[]
        curdir=os.getcwd()
        qmtorsionfolderpath=folderpath+'qm-torsion'
        os.chdir(qmtorsionfolderpath)
        excludedfoldnames=[]
        files=os.listdir()
        for f in files:
            filenamesplit= f.split('.')
            if filenamesplit[1]=='txt':
                if 'fit' in f:
                    newsplit=filenamesplit[0].split('fit-')
                    array=newsplit[1].split('-')
                    torsion=[int(i) for i in array]
                    listoftorsions.append(torsion)
                    listofqmminusmmtxtfiles.append(qmtorsionfolderpath+f)
        os.chdir(curdir)
        return listoftorsions,listofqmminusmmtxtfiles

    def GrabArrayData(self,filepath):
        anglearray=[]
        qmminusmmarray=[]
        with open(filepath) as infile:
            for line in infile:
                linesplit=line.split()
                angle=float(linesplit[0])
                qmminusmm=float(linesplit[1])
                anglearray.append(angle)
                qmminusmmarray.append(qmminusmm)
        return anglearray,qmminusmmarray


    def FindFragmentJobPath(self,parentrotbnd,testcasepath):
        curdir=os.getcwd()
        os.chdir(testcasepath+'Fragmenter/')
        fragfolders=os.listdir()
        fragpath=None
        for f in fragfolders:
            os.chdir(f)
            fragpath=self.RecursiveFolderSearch('qm-torsion',parentrotbnd)
            if fragpath!=None:
                os.chdir(curdir)
                return fragpath
            os.chdir('..')
        os.chdir(curdir)
        return fragpath

    def RecursiveFolderSearch(self,foldername,parentrotbnd):
        path=None
        for root, dirs, files in os.walk(".", topdown=False):
            for name in dirs:
               path=os.path.join(root, name)
               if self.DoesFolderExistInDirectory(path,foldername)==True:
                   if self.ReadTorsionsFile(parentrotbnd)==True:
                       return path
        return path

    def DoesFolderExistInDirectory(self,path,foldername):
        curdir=os.getcwd()
        os.chdir(path)
        foundfoldername=False
        for f in os.listdir():
            if f==foldername:
                foundfoldername=True
        os.chdir(curdir)
        return foundfoldername

    def ReadTorsionsFile(self,parentrotbnd):
        foundparent=False
        if os.path.isfile('torsions.txt'):
            temp=open('torsions.txt','r')
            results=temp.readlines()
            temp.close()
            firstline=results[0]
            if parentrotbnd in firstline:
                foundparent=True
        return foundparent

    def ReadDictionaryFromFile(self,filepath):
        curdir=os.getcwd()
        os.chdir(filepath)
        dictionary=json.load(open("parentindextofragindex.txt"))
        os.chdir(curdir)
        return dictionary


    def test_Fragmenter(self):
        print('Testing Fragmenter')
        os.chdir(self.examplepath)
        ptype=poltype.PolarizableTyper(readinionly=True)
        head,tail=os.path.split(self.examplepath)
        testfoldername='TestFragmenter/'+tail
        testcasepath=self.MakeTestCasePath(testfoldername)
        self.GenericFolderCopy(testcasepath)
        os.chdir(testcasepath)
        poltype.PolarizableTyper(dontfrag=False,structure=ptype.molstructfname,poltypeini=False,suppressdipoleerr=True,optmethod=ptype.optmethod,toroptmethod=ptype.toroptmethod,espmethod=ptype.espmethod,torspmethod=ptype.torspmethod,dmamethod=ptype.dmamethod,torspbasisset=ptype.torspbasisset,espbasisset=ptype.espbasisset,dmabasisset=ptype.dmabasisset,toroptbasisset=ptype.toroptbasisset,optbasisset=ptype.optbasisset,bashrcpath=ptype.bashrcpath,externalapi=ptype.externalapi,use_gaus=ptype.use_gaus,use_gausoptonly=ptype.use_gausoptonly,isfragjob=False,poltypepath=ptype.poltypepath,numproc=ptype.numproc,maxmem=ptype.maxmem,maxdisk=ptype.maxdisk,printoutput=True)
        listofparenttorsions,listofparentqmminusmmtxtfiles=self.GrabTorsions(exampleparentfolder+examplefolder) 
        for i in range(len(listofparenttorsions)):
            torsion=listofparenttorsions=[i]
            parentqmminusmmtxtfile=listofparentqmminusmmtxtfiles[i]
            parentrotbnd=str(torsion[1])+' '+str(torsion[2])
            filepath=self.FindFragmentJobPath(parentrotbnd,testcasepath)
            listoffragtorsions,listoffragqmminusmmtxtfiles=self.GrabTorsions(filepath)
            parentindextofragindex=self.ReadDictionaryFromFile(filepath)
            fragrotbnd=str(parentindextofragindex[torsion[1]])+'-'+str(parentindextofragindex[torsion[2]])
            revfragrotbnd=str(parentindextofragindex[torsion[2]])+'-'+str(parentindextofragindex[torsion[1]])
            for j in range(len(listoffragtorsions)):
                fragtorsion=listoffragtorsions[j]
                fragqmminusmmtxtfile=listoffragqmminusmmtxtfiles[j]
                rotbnd=str(fragtorsion[1])+' '+str(fragtorsion[2])
                if fragrotbnd==rotbnd or revfragrotbnd==rotbnd:
                    parentanglearray,parentqmminusmmarray=self.GrabArrayData(parentqmminusmmtxtfile)
                    fraganglearray,fragqmminusmmarray=self.GrabArrayData(fragqmminusmmtxtfile)
                    if len(parentanglearray)!=len(fraganglearray):
                        print('parentanglearray',parentanglearray,'fraganglearray',fraganglearray)
                        print('length of parent and fragment energy arrays are not same length so cannot compare')
                        continue
                    def RMSDW(c):
                        return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.divide(parentqmminusmmarray,weight),c))))
                    resultW=fmin(RMSDW,.5)
                    minRMSDW=RMSDW(resultW[0])

                    def RMSD(c):
                        return numpy.sqrt(numpy.mean(numpy.square(numpy.add(parentqmminusmmarray,c))))
                    result=fmin(RMSD,.5)
                    minRMSD=RMSD(result[0])
                    if minRMSD>1 and minRMSDW>1:
                        print('parent qm-mm array has bad fit, will not compare to fragment')
                        continue
                    else:

                        def RMSD(c):
                            return numpy.sqrt(numpy.mean(numpy.square(numpy.add(numpy.subtract(parentqmminusmmarray,fragqmminusmmarray,c)))))
                        result=fmin(RMSD,.5)
                        minRMSD=RMSD(result[0])
                        self.assertLessEqual(minRMSD, 1) 




  
    def test_MethylamineCommonInputs(self):
        print('Testing common inputs')
        testcasefolder='TestMethylamineCommonInputs/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='SymmetryMethylamine/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)
    

    def test_MethylaminePsi4(self):
        print('Testing Gaussian')
        testcasefolder='TestMethylamineGaus/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='MethylamineGaus/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,use_gaus=False,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)


    def test_MethylaminePsi4SPOnly(self):
        print('Testing Gaussian OPT only')
        testcasefolder='TestMethylamineGausOptOnly/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='MethylamineGausOptOnly/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,use_gausoptonly=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)

    
    def test_MethylamineDMAESPOptions(self):
        print('Testing DMA and ESP options')
        testcasefolder='TestMethylamineDMAESPOptions/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='Methylamine/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,qmonly=True,espbasisset='aug-cc-pVTZ',dmabasisset='cc-pVTZ',structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
   

 
    def test_MethylamineGeometryOptimizationOptions(self):
        print('Testing Geometry Optimization options ')
        testcasefolder='TestMethylamineGeometryOptimizationOptions/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='MethylamineGeometryOptimizationOptions/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',optbasisset='6-311G*',gausoptcoords='cartesian',freq=True,optpcm=True,optmaxcycle=500,optmethod='HF',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)
    
    def test_MethanolTorsion(self):
        print('Testing torsion options')
        testcasefolder='TestMethanolTorsion/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='MethanolTorsion/'
        examplestructure='methanol.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,rotalltors=True,optpcm=True,torsppcm=True,toroptpcm=True,torsionrestraint=.2,foldnum=4,tordatapointsnum=13,toroptmethod='HF',torspmethod='MP2',toroptbasisset='6-311G*',torspbasisset='6-311++G(2d,2p)',structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)

    def test_MethanolTorsionGaus(self):
        print('Testing torsion options')
        testcasefolder='TestMethanolTorsionGaus/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='MethanolTorsionGaus/'
        examplestructure='methanol.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,use_gaus=True,rotalltors=True,optpcm=True,torsppcm=True,toroptpcm=True,torsionrestraint=.2,foldnum=4,tordatapointsnum=13,toroptmethod='HF',torspmethod='MP2',toroptbasisset='6-311G*',torspbasisset='6-311++G(2d,2p)',structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)
     
 
    def test_MethanolOnlyRotBnd(self):
        print('Testing onlyrotbnd')
        testcasefolder='TestMethanolOnlyRotBnd/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='MethanolOnlyRotBnd/'
        examplestructure='methanol.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,onlyrotbnd='1,2',structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)


    def test_MethanolFitRotBndsList(self):
        print('Testing fitrotbndslist')
        testcasefolder='TestMethanolFitRotBndsList/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='MethanolFitRotBndsList/'
        examplestructure='methanol.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,fitrotbndslist='1 2',structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)


    def test_MethylamineDontDoTorFit(self):
        print('Testing dontdotorfit')
        testcasefolder='TestMethylamineDontDoTorFit/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='Methylamine/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,rotalltors=True,dontdotorfit=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)

    
    def test_MethylamineDontDoTor(self):
        print('Testing dontdotor')
        testcasefolder='TestMethylamineDontDoTor/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='Methylamine/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,dontdotor=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
    

    def test_ModifiedAminoAcid(self):
        print('Testing modified amino acids')
        testcasefolder='TestModifiedAminoAcid/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='ModifiedAminoAcid/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,use_gaus=True,amoebabioprmpath='/home/bdw2292/Tinker-release/params/amoebabio18.prm',unmodifiedproteinpdbname='SNase_WT_H.pdb',mutatedresiduenumber='102',mutatedsidechain='CNC.sdf',modifiedresiduepdbcode='CNC',numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True) 
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)

    def test_Methane(self):
        print('Testing methane symmetry')
        testcasefolder='TestSymmetryMethane/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='SymmetryMethane/'
        examplestructure='methane.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)

    def test_Water(self):
        print('Testing water symmetry')
        testcasefolder='TestSymmetryWater/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='SymmetryWater/'
        examplestructure='water.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)

    def test_Acetamide(self):
        print('Testing acetamide symmetry')
        testcasefolder='TestSymmetryAcetamide/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='SymmetryAcetamide/'
        examplestructure='acetamide.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)

    def test_Ethene(self):
        print('Testing ethene symmetry')
        testcasepath='TestSymmetryEthene/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='SymmetryEthene/'
        examplestructure='ethene.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)

    def test_Ammonia(self):
        print('Testing ammonia symmetry')
        testcasefolder='TestSymmetryAmmonia/'
        testcasepath=self.MakeTestCasePath(testcasefolder)
        examplefolder='SymmetryAmmonia/'
        examplestructure='ammonia.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)


    def test_Aniline(self):
        print('Testing aniline symmetry')
        testcasefolder='TestSymmetryAniline/'
        testcasepath=self.MakeTestCasePath(testcasepath)
        examplefolder='SymmetryAniline/'
        examplestructure='aniline.sdf'
        examplekeyfilepath=self.GenericCopy(testcasepath,examplefolder,examplestructure)
        ptype=poltype.PolarizableTyper(dontfrag=True,structure=examplestructure,numproc=4,maxmem='20GB',maxdisk='100GB',ptypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,ptype,examplekeyfilepath)
    

if __name__ == '__main__':
    command_line_param = sys.argv[1:]
    del sys.argv[1:]
    examplepath=None
    opts, xargs = getopt.getopt(command_line_param,'e:m',["method=","examplepath="])
    for o, a in opts:
        if o in ("--examplepath"):
            examplepath=a
        elif o in ("--method"):
            method=a
 
    suite = unittest.TestSuite()
    suite.addTest(TestPoltype(method,examplepath))
    runner = unittest.TextTestRunner()
    runner.run(suite)
