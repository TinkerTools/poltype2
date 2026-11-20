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


# Expected protected bonds outside aromatic rings (one-based atom indices).
AMINE_CASES = {
    'CNc1ccccc1': {(2, 3)},                         # secondary amine
    'CN(C)c1ccccc1': {(2, 4)},                      # tertiary amine
    'C[NH2+]c1ccccc1': {(2, 3)},                    # protonated secondary amine
    'C[NH+](C)c1ccccc1': {(2, 4)},                  # protonated tertiary amine
    'CNc1ncccc1': {(2, 3)},                        # heteroaromatic ring
    'CN(C)c1cccs1': {(2, 4)},                       # heteroaromatic ring
    'CNn1cccc1': {(2, 3)},                         # aromatic nitrogen attachment
    'CN(C)c1ccc2ccccc2c1': {(2, 4)},                # fused aromatic rings
    'N(c1ccccc1)c1ccccc1': {(1, 2), (1, 8)},        # two directly attached rings
    'CN(c1ccccc1)c1ccccc1': {(2, 3), (2, 9)},       # tertiary diarylamine
    'c1ccccc1N1CCCC1': {(6, 7)},                    # cyclic amine
}


CARBONYL_CASES = {
    'O=Cc1ccccc1': {(1, 2), (2, 3)},                  # benzaldehyde
    'CC(=O)c1ccccc1': {(2, 3), (2, 4)},               # acetophenone
    'O=C(O)c1ccccc1': {(1, 2), (2, 3), (2, 4)},      # benzoic acid
    'O=C([O-])c1ccccc1': {(1, 2), (2, 3), (2, 4)},   # benzoate
    'NC(=O)c1ccccc1': {(1, 2), (2, 3), (2, 4)},      # benzamide
    'COC(=O)c1ccccc1': {(2, 3), (3, 4), (3, 5)},     # methyl benzoate
    'CC(=O)Nc1ccccc1': {(2, 3), (2, 4), (4, 5)},     # acetanilide
    'CC(=O)Oc1ccccc1': {(2, 3), (2, 4), (4, 5)},     # phenyl acetate
    'CC(=O)N(C)c1ccccc1': {(2, 3), (2, 4), (4, 6)},  # N-methylacetanilide
    'CN(C)C(=O)c1ccccc1': {(2, 4), (4, 5), (4, 6)},  # N,N-dimethylbenzamide
    'O=CNc1ccccc1': {(1, 2), (2, 3), (3, 4)},        # formanilide
    'O=COc1ccccc1': {(1, 2), (2, 3), (3, 4)},        # phenyl formate
    'CC(=O)n1cccc1': {(2, 3), (2, 4)},               # N-acetylpyrrole
    'O=Cc1cccs1': {(1, 2), (2, 3)},                  # heteroaromatic aldehyde
    'CC(=O)c1ncccc1': {(2, 3), (2, 4)},               # heteroaromatic ketone
    'O=C(O)c1cccs1': {(1, 2), (2, 3), (2, 4)},       # heteroaromatic acid
    'NC(=O)c1ncccc1': {(1, 2), (2, 3), (2, 4)},      # heteroaromatic amide
    'COC(=O)c1cccs1': {(2, 3), (3, 4), (3, 5)},      # heteroaromatic ester
    'O=Cc1ccc2ccccc2c1': {(1, 2), (2, 3)},           # fused aromatic rings
    'O=Cc1ccccccccc1': {(1, 2), (2, 3)},             # larger aromatic ring
    'O=C(c1ccccc1)c1ccccc1': {(1, 2), (2, 3), (2, 9)},  # benzophenone
    'O=C(Nc1ccccc1)c1ccccc1': {(1, 2), (2, 3), (3, 4), (2, 10)},
    'O=C(Oc1ccccc1)c1ccccc1': {(1, 2), (2, 3), (3, 4), (2, 10)},
}


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
        for smiles in ('Cc1ccccc1', 'COc1ccccc1', 'CSc1ccccc1',
                       'NCCc1ccccc1', 'NC1CCCCC1', 'CNCCc1ccccc1',
                       'CN(C)Cc1ccccc1', 'CNC1CCCCC1', 'CN(C)C1CCCCC1',
                       'N#Cc1ccccc1', 'O=[N+]([O-])c1ccccc1', 'CN=Nc1ccccc1'):
            with self.subTest(smiles=smiles):
                self.assertEqual(aromatic_group_bonds(Chem.MolFromSmiles(smiles)), set())

    def test_substituted_amines_stay_with_rings_but_not_alkyl_substituents(self):
        for smiles, attachments in AMINE_CASES.items():
            for explicit_h in (False, True):
                with self.subTest(smiles=smiles, explicit_h=explicit_h):
                    mol = Chem.MolFromSmiles(smiles)
                    if explicit_h:
                        mol = Chem.AddHs(mol)
                    ring_bonds = {
                        tuple(sorted((b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1)))
                        for b in mol.GetBonds() if b.GetIsAromatic()}
                    self.assertEqual(aromatic_group_bonds(mol), attachments | ring_bonds)

    def test_methyl_bond_still_cuttable_on_aminotoluene(self):
        protected = aromatic_group_bonds(Chem.MolFromSmiles('Nc1ccc(C)cc1'))
        self.assertIn((1, 2), protected)
        self.assertNotIn((5, 6), protected)

    def test_carbonyl_groups_and_directly_attached_rings_stay_together(self):
        for smiles, core_bonds in CARBONYL_CASES.items():
            for explicit_h in (False, True):
                with self.subTest(smiles=smiles, explicit_h=explicit_h):
                    mol = Chem.MolFromSmiles(smiles)
                    if explicit_h:
                        mol = Chem.AddHs(mol)
                    ring_bonds = {
                        tuple(sorted((b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1)))
                        for b in mol.GetBonds() if b.GetIsAromatic()}
                    self.assertEqual(aromatic_group_bonds(mol), core_bonds | ring_bonds)

    def test_carbonyl_groups_with_spacers_or_without_aromatic_rings(self):
        for smiles in ('CC(=O)Cc1ccccc1', 'O=CCc1ccccc1', 'O=C(O)Cc1ccccc1',
                       'NC(=O)Cc1ccccc1', 'COC(=O)Cc1ccccc1',
                       'CC(=O)NCc1ccccc1', 'CC(=O)OCc1ccccc1',
                       'CC(=O)C', 'CCC=O', 'CC(=O)O', 'CC(=O)N', 'CC(=O)OC',
                       'CC(=O)C1CCCCC1', 'O=CC1CCCCC1', 'O=C(O)C1CCCCC1',
                       'NC(=O)C1CCCCC1', 'COC(=O)C1CCCCC1',
                       'O=C(Cl)c1ccccc1', 'O=C(SC)c1ccccc1'):
            for explicit_h in (False, True):
                with self.subTest(smiles=smiles, explicit_h=explicit_h):
                    mol = Chem.MolFromSmiles(smiles)
                    if explicit_h:
                        mol = Chem.AddHs(mol)
                    self.assertEqual(aromatic_group_bonds(mol), set())

    def test_unrelated_ring_and_methyl_bonds_remain_cuttable(self):
        mol = Chem.MolFromSmiles('O=Cc1ccc(-c2ccccc2)cc1C')
        self.assertEqual(aromatic_group_bonds(mol), {
            (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 13), (13, 14), (3, 14)})

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
        for smiles in ('Cc1ccccc1', 'CC(=O)Cc1ccccc1', 'CC(=O)N', 'CC(=O)OC'):
            with self.subTest(smiles=smiles), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'molecule.mol'
                Chem.MolToMolFile(Chem.MolFromSmiles(smiles), str(path))
                poltype = SimpleNamespace(molstructfname=str(path), peditinfile='original.in')
                self.assertEqual(prepare_poledit_groups(poltype), 'original.in')


@unittest.skipUnless(os.environ.get('POLEDIT_TEST_EXE'), 'Set POLEDIT_TEST_EXE for Tinker integration')
class PoleditIntegrationTests(unittest.TestCase):
    def test_overrides_are_used_during_poledit(self):
        for smiles in ('Nc1ccccc1', '[NH3+]c1ccccc1', 'Oc1ccccc1',
                       'Sc1ccccc1', 'Nn1cccc1', 'Nc1ccc(C)cc1',
                       'Nc1ncnc2[nH]cnc12', 'Nc1cc[nH]c(=O)n1',
                       'Nc1nc2[nH]cnc2c(=O)[nH]1', *AMINE_CASES, *CARBONYL_CASES):
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
                if smiles not in CARBONYL_CASES and smiles not in AMINE_CASES:
                    self.assertNotIn((1, 2), finalcuts)
                if smiles in ('CNc1ccccc1', 'CN(C)c1ccccc1'):
                    self.assertIn((1, 2), finalcuts)
                    if smiles == 'CN(C)c1ccccc1':
                        self.assertIn((2, 3), finalcuts)
                if smiles == 'Nc1ccc(C)cc1':
                    self.assertIn((5, 6), finalcuts)
                if smiles == 'O=C(O)c1ccccc1':
                    # The entire COOH, including O-H, stays with the ring.
                    self.assertEqual(finalcuts, set())
                if smiles == 'COC(=O)c1ccccc1':
                    # Only O-CH3 is cut; the methyl hydrogens stay with CH3.
                    self.assertEqual(finalcuts, {(1, 2)})
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
