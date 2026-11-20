"""Preserve selected aromatic substituents in POLEDIT groups."""

import os
import shlex
import shutil
import tempfile

from rdkit import Chem


# Match only the functional-group core. Recursive constraints identify the
# group without including alkyl substituents in the protected bonds.
_CARBONYL_GROUPS = tuple(Chem.MolFromSmarts(smarts) for smarts in (
    '[CX3H1;$([C]-[#6])]=[OX1]',                       # aldehyde
    '[CX3;$([C](-[#6])-[#6])]=[OX1]',                 # ketone
    '[CX3;H1,$([C]-[#6])](=[OX1])-[OX2H1,OX1-]',      # acid/carboxylate
    '[CX3;H1,$([C]-[#6])](=[OX1])-[#7X3]',             # amide
    '[CX3;H1,$([C]-[#6])](=[OX1])-[OX2H0;$([O]-[#6])]',  # ester
))


def aromatic_group_bonds(mol):
    """Return protected bonds as sorted pairs of one-based atom indices.

    Keep amine N/NH/NH2 (including protonated forms), terminal OH/SH,
    and ketone, aldehyde, carboxylic acid,
    amide, and ester cores with directly attached aromatic rings. Attachment
    through an amide N or ester O also qualifies. Aromatic atoms may have
    any element and belong to rings of any size. Include the aromatic ring
    bonds so the whole ring stays together even if Tinker's aromaticity
    heuristic does not recognize it. Alkyl substituent bonds retain the
    installed POLEDIT's default treatment.
    """
    protected = set()
    ring_atoms = set()
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            continue
        hydrogens = atom.GetTotalNumHs(includeNeighbors=True)
        heavy = [a for a in atom.GetNeighbors() if a.GetAtomicNum() != 1]
        amine = (atom.GetAtomicNum() == 7 and
                 atom.GetTotalDegree() in (3, 4) and
                 all(b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds()))
        terminal_oh_sh = (atom.GetAtomicNum() in (8, 16) and
                          hydrogens == 1 and len(heavy) == 1)
        if not (amine or terminal_oh_sh):
            continue
        for aromatic in heavy:
            if aromatic.GetIsAromatic() and aromatic.IsInRing():
                protected.add(tuple(sorted((atom.GetIdx()+1, aromatic.GetIdx()+1))))
                ring_atoms.add(aromatic.GetIdx())

    for pattern in _CARBONYL_GROUPS:
        for match in mol.GetSubstructMatches(pattern):
            attachments = set()
            for index in match:
                for neighbor in mol.GetAtomWithIdx(index).GetNeighbors():
                    if neighbor.GetIsAromatic() and neighbor.IsInRing():
                        attachments.add(tuple(sorted((index+1, neighbor.GetIdx()+1))))
                        ring_atoms.add(neighbor.GetIdx())
            if not attachments:
                continue
            protected.update(attachments)
            for bond in pattern.GetBonds():
                protected.add(tuple(sorted((match[bond.GetBeginAtomIdx()]+1,
                                            match[bond.GetEndAtomIdx()]+1))))

    pending = list(ring_atoms)
    while pending:
        atom = mol.GetAtomWithIdx(pending.pop())
        for bond in atom.GetBonds():
            if not (bond.GetIsAromatic() and bond.IsInRing()):
                continue
            other = bond.GetOtherAtom(atom)
            protected.add(tuple(sorted((atom.GetIdx()+1, other.GetIdx()+1))))
            if other.GetIdx() not in ring_atoms:
                ring_atoms.add(other.GetIdx())
                pending.append(other.GetIdx())
    return protected


def poledit_cut_bonds(xyzfile, keyfile):
    """Recover the installed POLEDIT's cuts from its uncondensed output."""
    with open(xyzfile) as handle:
        count = int(handle.readline().split()[0])
        atoms = {}
        for line in handle:
            fields = line.split()
            if fields:
                atoms[int(fields[0])] = (int(fields[5]), list(map(int, fields[6:])))
    if len(atoms) != count:
        raise ValueError('Incomplete POLEDIT coordinates for polarization groups')

    groups = {}
    with open(keyfile) as handle:
        for line in handle:
            fields = line.split('#', 1)[0].split()
            if fields and fields[0].lower() == 'polarize':
                # POLEDIT prints real-valued polarizability/damping fields
                # with decimal points; trailing integers are group types.
                groups[int(fields[1])] = {
                    int(value) for value in fields[3:] if value.isdigit()}
    if any(typenum not in groups for typenum, _ in atoms.values()):
        raise ValueError('Missing POLEDIT polarize parameters for group detection')

    cuts = set()
    for index, (typenum, neighbors) in atoms.items():
        for neighbor in neighbors:
            if index < neighbor and atoms[neighbor][0] not in groups[typenum]:
                cuts.add((index, neighbor))
    return cuts


def prepare_poledit_groups(poltype):
    """Build a manual cut list before POLEDIT removes induced dipoles.

    A temporary automatic run supplies the installed Tinker's default cuts;
    only protected aromatic bonds are removed from that list. The final run
    still starts from the original GDMA moments. Molecules without matching
    substituents use the original input and incur no extra POLEDIT run.
    """
    mol = Chem.MolFromMolFile(poltype.molstructfname, removeHs=False)
    if mol is None:
        raise ValueError('Cannot read molecule for polarization group rules')
    protected = aromatic_group_bonds(mol)
    if not protected:
        return poltype.peditinfile

    with open(poltype.peditinfile) as handle:
        lines = handle.read().splitlines()
    if [line.strip() for line in lines[-3:]] != ['2', 'N', 'Y']:
        raise ValueError('Aromatic polarization group rules require the '
                         'Poltype POLEDIT input ending in 2, N, Y')
    # Resolve the parameter file before entering the temporary directory.
    lines[0] = os.path.abspath(lines[0].strip() or poltype.paramhead)
    poltype.WriteToLog('Preserving aromatic amine N/NH/NH2 (including protonated forms), OH/SH and '
                      'ketone/aldehyde/acid/amide/ester polarization groups')
    with tempfile.TemporaryDirectory(prefix='poltype-polarization-') as workdir:
        shutil.copyfile(poltype.gdmafname, os.path.join(workdir, 'probe.gdmaout'))
        probeinput = os.path.join(workdir, 'probe.in')
        with open(probeinput, 'w') as handle:
            handle.write('\n'.join(lines) + '\n')
        executable = poltype.peditexe
        if os.path.dirname(executable):
            executable = os.path.abspath(executable)
        command = (f'cd {shlex.quote(workdir)} && '
                   f'{shlex.quote(executable)} 1 probe.gdmaout < probe.in')
        poltype.call_subsystem([command], True)
        cuts = poledit_cut_bonds(os.path.join(workdir, 'probe.xyz'),
                                os.path.join(workdir, 'probe.key'))

    cuts.difference_update(protected)
    output = os.path.splitext(poltype.peditinfile)[0] + '-groups.txt'
    with open(output, 'w') as handle:
        handle.write('\n'.join(lines[:-3]) + '\n3\n')
        for first, second in sorted(cuts):
            handle.write(f'{first} {second}\n')
        handle.write('\nN\nY\n')
    return output
