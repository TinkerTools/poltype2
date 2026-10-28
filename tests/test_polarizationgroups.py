import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'PoltypeModules'))
from polarizationgroups import (aromatic_group_bonds, poledit_cut_bonds,
                                prepare_poledit_groups)


class AromaticGroupTests(unittest.TestCase):
    def test_terminal_substituents_on_carbon_and_heteroaromatic_rings(self):
        for smiles in ('Nc1ccccc1', '[NH3+]c1ccccc1', 'Oc1ccccc1',
                       'Sc1ccccc1', 'Nc1ncccc1', 'Nn1cccc1',
                       'Nc1cccs1', 'Nc1ccc2ccccc2c1', 'Nc1ccccccccc1'):
            for explicit_h in (False, True):
                with self.subTest(smiles=smiles, explicit_h=explicit_h):
                    mol = Chem.MolFromSmiles(smiles)
                    if explicit_h:
                        mol = Chem.AddHs(mol)
                    protected = aromatic_group_bonds(mol)
                    self.assertIn((1, 2), protected)
                    self.assertTrue(all(
                        tuple(sorted((b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1)))
                        in protected for b in mol.GetBonds() if b.GetIsAromatic()))

    def test_methyl_and_other_substituents_are_not_protected(self):
        for smiles in ('Cc1ccccc1', 'CNc1ccccc1', 'COc1ccccc1',
                       'CSc1ccccc1', 'NCCc1ccccc1', 'NC1CCCCC1'):
            with self.subTest(smiles=smiles):
                self.assertEqual(aromatic_group_bonds(Chem.MolFromSmiles(smiles)), set())

    def test_methyl_bond_still_cuttable_on_aminotoluene(self):
        protected = aromatic_group_bonds(Chem.MolFromSmiles('Nc1ccc(C)cc1'))
        self.assertIn((1, 2), protected)
        self.assertNotIn((5, 6), protected)

    def test_nucleobases(self):
        for smiles in ('Nc1ncnc2[nH]cnc12', 'Nc1cc[nH]c(=O)n1',
                       'Nc1nc2[nH]cnc2c(=O)[nH]1'):
            with self.subTest(smiles=smiles):
                self.assertIn((1, 2), aromatic_group_bonds(Chem.MolFromSmiles(smiles)))

    def test_read_default_cuts_with_different_damping_formats(self):
        with tempfile.TemporaryDirectory() as directory:
            xyz = Path(directory) / 'probe.xyz'
            key = Path(directory) / 'probe.key'
            xyz.write_text('3\n1 N 0 0 0 10 2\n2 C 1 0 0 20 1 3\n3 C 2 0 0 30 2\n')
            for damping in ('', ' 0.3900', ' 0.3900 0.7000'):
                with self.subTest(damping=damping):
                    key.write_text(f'polarize 10 1.0000{damping}\n'
                                   f'polarize 20 1.0000{damping} 30\n'
                                   f'polarize 30 1.0000{damping} 20\n')
                    self.assertEqual(poledit_cut_bonds(xyz, key), {(1, 2)})

    def test_no_probe_when_no_protected_substituent(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'toluene.mol'
            Chem.MolToMolFile(Chem.MolFromSmiles('Cc1ccccc1'), str(path))
            poltype = SimpleNamespace(molstructfname=str(path), peditinfile='original.in')
            self.assertEqual(prepare_poledit_groups(poltype), 'original.in')


@unittest.skipUnless(os.environ.get('POLEDIT_TEST_EXE'), 'Set POLEDIT_TEST_EXE for Tinker integration')
class PoleditIntegrationTests(unittest.TestCase):
    def test_overrides_are_used_during_poledit(self):
        for smiles in ('Nc1ccccc1', '[NH3+]c1ccccc1', 'Oc1ccccc1',
                       'Sc1ccccc1', 'Nn1cccc1', 'Nc1ccc(C)cc1',
                       'Nc1ncnc2[nH]cnc12', 'Nc1cc[nH]c(=O)n1',
                       'Nc1nc2[nH]cnc2c(=O)[nH]1'):
            with self.subTest(smiles=smiles), tempfile.TemporaryDirectory() as directory:
                workdir = Path(directory)
                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=1234), 0)
                AllChem.UFFOptimizeMolecule(mol)
                molfile = workdir / 'molecule.mol'
                Chem.MolToMolFile(mol, str(molfile))
                gdma = workdir / 'molecule.gdmaout'
                # Zero moments suffice to exercise grouping and polarization.
                with gdma.open('w') as handle:
                    for atom in mol.GetAtoms():
                        point = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                        handle.write(f'{atom.GetSymbol():<11}x ={point.x:10.6f}'
                                     f'  y ={point.y:10.6f}  z ={point.z:10.6f}\n\n')
                    handle.write('Total multipoles\n')
                header = Path(__file__).resolve().parents[1] / 'ParameterFiles/amoebabio18_header.prm'
                pedit = workdir / 'molecule-peditin.txt'
                pedit.write_text(f'{header}\n\nA\n' + ''.join(
                    f'{a.GetIdx()+1} 0.500\n' for a in mol.GetAtoms()) + '\n2\nN\nY\n')

                def run_commands(commands, wait):
                    for command in commands:
                        subprocess.run(command, shell=True, check=True, timeout=30,
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

                poltype = SimpleNamespace(
                    molstructfname=str(molfile), peditinfile=str(pedit),
                    paramhead=str(header), gdmafname=str(gdma),
                    peditexe=os.environ['POLEDIT_TEST_EXE'],
                    WriteToLog=lambda message: None, call_subsystem=run_commands)
                finalinput = prepare_poledit_groups(poltype)
                with open(finalinput) as handle:
                    result = subprocess.run(
                        [poltype.peditexe, '1', gdma.name], cwd=directory,
                        stdin=handle, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        check=True, timeout=30)
                self.assertIn(b'Polarization Groups for Multipole Sites', result.stdout)
                finalcuts = poledit_cut_bonds(workdir / 'molecule.xyz', workdir / 'molecule.key')
                protected = aromatic_group_bonds(mol)
                self.assertFalse(finalcuts & protected)
                self.assertNotIn((1, 2), finalcuts)
                if smiles == 'Nc1ccc(C)cc1':
                    self.assertIn((5, 6), finalcuts)
                # The manual list must preserve every other automatic cut.
                automatic = workdir / 'automatic'
                automatic.mkdir()
                shutil.copyfile(gdma, automatic / gdma.name)
                with pedit.open() as handle:
                    subprocess.run([poltype.peditexe, '1', gdma.name], cwd=automatic,
                                   stdin=handle, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   check=True, timeout=30)
                baseline = poledit_cut_bonds(automatic / 'molecule.xyz',
                                            automatic / 'molecule.key')
                self.assertEqual(finalcuts, baseline - protected)


if __name__ == '__main__':
    unittest.main()
