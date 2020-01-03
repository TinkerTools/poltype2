import unittest
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)+'/'
poltypepath=parentdir
sys.path.insert(0,poltypepath)
from poltype import PolarizableTyper
import shutil
from parmed.tinker import parameterfile


class TestPoltype(unittest.TestCase):
    global testcaseparentpath
    testcaseparentpath=poltypepath+'TestCases'
    if os.path.isdir(testcaseparentpath):
        shutil.rmtree(testcaseparentpath)
    os.mkdir(testcaseparentpath)
    os.chdir(testcaseparentpath)

    def MakeTestCasePath(self,testcasepath):
        if not os.path.isdir(testcasepath):
            os.mkdir(testcasepath)
        os.chdir(testcasepath)

    def GenericCopy(self,exampleparentfolder,testcasepath,examplefolder,examplestructure):
        examplekey5fname=examplestructure.replace('.sdf','.key_5')
        examplestructurepath=poltypepath+exampleparentfolder+examplefolder+examplestructure
        examplekeyfilepath=poltypepath+exampleparentfolder+examplefolder+examplekey5fname
        shutil.copy(examplestructurepath,testcasepath+examplestructure)
        return examplekeyfilepath

    def GenericTest(self,exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath):
        testparams=poltype.main()
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


           
    def test_MethylamineCommonInputs(self):
        print('Testing common inputs')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethylamineCommonInputs/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='SymmetryMethylamine/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)
    

    def test_MethylaminePsi4(self):
        print('Testing Gaussian')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethylamineGaus/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='MethylamineGaus/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(use_gaus=False,structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)


    def test_MethylaminePsi4SPOnly(self):
        print('Testing Gaussian OPT only')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethylamineGausOptOnly/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='MethylamineGausOptOnly/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(use_gausoptonly=True,structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)

    
    def test_MethylamineDMAESPOptions(self):
        print('Testing DMA and ESP options')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethylamineDMAESPOptions/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='Methylamine/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(qmonly=True,espbasisset='aug-cc-pVTZ',dmabasisset='cc-pVTZ',structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
    
    def test_MethylamineGeometryOptimizationOptions(self):
        print('Testing Geometry Optimization options ')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethylamineGeometryOptimizationOptions/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='MethylamineGeometryOptimizationOptions/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',optbasisset='6-311G*',gausoptcoords='cartesian',freq=True,optpcm=True,optmaxcycle=500,optmethod='HF',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)
    
    def test_MethanolTorsion(self):
        print('Testing torsion options')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethanolTorsion/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='MethanolTorsion/'
        examplestructure='methanol.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(rotalltors=True,optpcm=True,torsppcm=True,toroptpcm=True,torsionrestraint=.2,foldnum=4,tordatapointsnum=13,toroptmethod='HF',torspmethod='MP2',toroptbasisset='6-311G*',torspbasisset='6-311++G(2d,2p)',structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)

    def test_MethanolTorsionGaus(self):
        print('Testing torsion options')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethanolTorsionGaus/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='MethanolTorsionGaus/'
        examplestructure='methanol.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(use_gaus=True,rotalltors=True,optpcm=True,torsppcm=True,toroptpcm=True,torsionrestraint=.2,foldnum=4,tordatapointsnum=13,toroptmethod='HF',torspmethod='MP2',toroptbasisset='6-311G*',torspbasisset='6-311++G(2d,2p)',structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)
     
 
    def test_MethanolOnlyRotBnd(self):
        print('Testing onlyrotbndlist')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethanolOnlyRotBnd/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='MethanolOnlyRotBnd/'
        examplestructure='methanol.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(onlyrotbndlist='1,2',structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)


    def test_MethanolFitRotBndsList(self):
        print('Testing fitrotbndslist')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethanolFitRotBndsList/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='MethanolFitRotBndsList/'
        examplestructure='methanol.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(fitrotbndslist='1 2',structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)


    def test_MethylamineDontDoTorFit(self):
        print('Testing dontdotorfit')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethylamineDontDoTorFit/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='Methylamine/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(rotalltors=True,dontdotorfit=True,structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)

    
    def test_MethylamineDontDoTor(self):
        print('Testing dontdotor')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestMethylamineDontDoTor/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='Methylamine/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(dontdotor=True,structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
    

    def test_ModifiedAminoAcid(self):
        print('Testing modified amino acids')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestModifiedAminoAcid/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='ModifiedAminoAcid/'
        examplestructure='methylamine.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(use_gaus=True,amoebabioprmpath='/home/bdw2292/Tinker-release/params/amoebabio18.prm',unmodifiedproteinpdbname='SNase_WT_H.pdb',mutatedresiduenumber='102',mutatedsidechain='CNC.sdf',modifiedresiduepdbcode='CNC',numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True) 
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)

    def test_Methane(self):
        print('Testing methane symmetry')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestSymmetryMethane/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='SymmetryMethane/'
        examplestructure='methane.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)

    def test_Water(self):
        print('Testing water symmetry')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestSymmetryWater/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='SymmetryWater/'
        examplestructure='water.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)

    def test_Acetamide(self):
        print('Testing acetamide symmetry')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestSymmetryAcetamide/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='SymmetryAcetamide/'
        examplestructure='acetamide.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)

    def test_Ethene(self):
        print('Testing ethene symmetry')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestSymmetryEthene/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='SymmetryEthene/'
        examplestructure='ethene.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)

    def test_Ammonia(self):
        print('Testing ammonia symmetry')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestSymmetryAmmonia/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='SymmetryAmmonia/'
        examplestructure='ammonia.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)


    def test_Aniline(self):
        print('Testing aniline symmetry')
        exampleparentfolder='Examples/'
        testcasepath=testcaseparentpath+'/'+'TestSymmetryAniline/'
        self.MakeTestCasePath(testcasepath)
        examplefolder='SymmetryAniline/'
        examplestructure='aniline.sdf'
        examplekeyfilepath=self.GenericCopy(exampleparentfolder,testcasepath,examplefolder,examplestructure)
        poltype=PolarizableTyper(structure=examplestructure,numproc=4,maxmem='5GB',maxdisk='20GB',poltypeini=False,printoutput=True)
        self.GenericTest(exampleparentfolder,testcasepath,examplefolder,examplestructure,poltype,examplekeyfilepath)




if __name__ == '__main__':
    unittest.main() 
