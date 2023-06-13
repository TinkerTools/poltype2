#!/usr/bin/env python

import sys
import argparse
import logging
from collections import defaultdict
import re
import numpy as np
import mdtraj
from rdkit.Chem import rdmolfiles, rdmolops
from rdkit import Chem
try:
    from openbabel import openbabel
except (ImportError, ModuleNotFoundError):
    import openbabel
logging.basicConfig(level=20)

def txyz_to_xyzblock(ftxyz):
    """ Convert txyz to xyz string
    """
    with open(ftxyz, 'r') as fh:
        lines = fh.readlines()
        try:
            natom = int(lines[0].split()[0])
        except (ValueError, IndexError):
            logging.warning(f"Unable to read number of atoms from Tinker xyz file {ftxyz}")
            return
        if len(lines) < natom+1:
            logging.warning(f"Wrong format for Tinker xyz file {ftxyz}")
            return
        outlines = ['%d\n'%natom, '\n']
        for iatom in range(natom):
            w = lines[iatom+1].split()
            try: 
                outlines.append( '%s %10.4f %10.4f %10.4f\n'%(w[1], float(w[2]), float(w[3]), float(w[4])))
            except (ValueError, IndexError):
                logging.warning(f"Unable to read atom coordinate in Tinker xyz file '{ftxyz}': {lines[iatom+1]}")
                return
        return ''.join(outlines)

def txyz_to_xyz(ftxyz, fxyz):
    traj = mdtraj.load_arc(ftxyz)
    traj.save_xyz(fxyz)

def txyz_to_pdb(ftxyz, fpdb):
    traj = mdtraj.load_arc(ftxyz)
    traj.save_pdb(fpdb)

def read_txyz_to_rdmol(ftxyz):
    xyz_str = txyz_to_xyzblock(ftxyz)
    if xyz_str is None:
        raise ValueError(f'Could not read Tinker xyz file {ftxyz}')
    obconvert = openbabel.OBConversion()
    obconvert.SetInFormat('xyz')
    obconvert.SetOutFormat('mol')
    obmol = openbabel.OBMol()
    obconvert.ReadString(obmol, xyz_str)
    mol_str = obconvert.WriteString(obmol)

    m1 = rdmolfiles.MolFromMolBlock(mol_str, removeHs=False)
    logging.info("Read molecule %s"%(rdmolfiles.MolToSmiles(m1)))
    return m1

    #tmpfile = "tmp_convert.pdb"
    #txyz_to_pdb(ftxyz, tmpfile)
    #m1 = rdmolfiles.MolFromPDBFile(tmpfile, removeHs=False)
    #return m1

def read_txyz_atomtypes(ftxyz):
    atypes = []
    with open(ftxyz) as fh:
        lines = fh.readlines()
        if len(lines) == 0:
            return atypes
        w = lines[0].split()
        assert len(w) > 0
        assert w[0].isdigit(), f"First line ('{w[0]}') in Tinker xyz should specify natoms"
        natoms = int(w[0])

        assert len(lines) >= natoms + 1

        for iatom in range(natoms):
            w = lines[iatom+1].split()
            assert len(w) >= 6
            assert w[0] == '%d'%(iatom+1)
            assert w[5].isdigit()

            atype = int(w[5])
            atypes.append(atype)
    return atypes

def read_key_atomclass(fkey):
    atype_map = {}
    with open(fkey) as fh:
        for line in fh:
            if not line.startswith('atom'):
                continue
            w = line.split()
            if len(w) < 7:
                continue
            try:
                atype = int(w[1])
                aclass = int(w[2])
                atype_map[atype] = aclass
            except ValueError:
                logging.warning(f"Wrong format. Skip line: '{line.strip()}'")
    return atype_map

class ParmMod(object):
    def __init__(self):
        self._terms = defaultdict(list)
        # force field terms that use atom type; by default atom class is used
        self.ATOM_TYPE_TERMS = ('polarize', 'quadrupole-scale', 'multipole')
        self.IMPLEMENTED_TERMS = ('quadrupole-scale', 'vdw', 'vdw-scale', 'polarize')
        self.PRM_FORMAT = {}
        self.PRM_FORMAT['vdw'] = ('%11d', '%20.4f', '%10.4f', '%10.4f')
        self.PRM_FORMAT['vdwpair'] = ('%7d', '%4d', '%15.4f', '%10.4f', '10.4f')

    @staticmethod
    def convert_data(string):
        for datatype in (int, float, str):
            try:
                return datatype(string)
            except ValueError:
                pass

    def read_inpfile(self, inpfile):
        """Read prms from an input file
        """
        with open(inpfile, 'r') as fh:
            for line in fh:
                # skip comment line 
                if line.startswith('#'):
                    continue
                # remove in-line comment
                words = line.split(' #')[0].split()
                if len(words) < 3:
                    continue
                term = words[0]
                smarts = words[1]
                prm = list(map(self.convert_data, words[2:]))
                self._terms[term].append((smarts, prm))

    def read_list(self, smarts_list, x_list, term=''):
        """Read prms from list
        """
        for smarts, x in zip(smarts_list, x_list):
            if hasattr(x, '__iter__'):
                prm = list(x)
            else:
                prm = [x]
            self._terms[term].append((smarts, prm))

    def reset_terms(self):
        self._terms = defaultdict(list)

    def match_atomtypes(self, ftxyz, fsdf=None):
        """Find atom types that match prm definition
        """
        atypes = read_txyz_atomtypes(ftxyz)
        prm_types = defaultdict(dict)
        if fsdf is None:
            m1 = read_txyz_to_rdmol(ftxyz)
        else:
            m1 = Chem.MolFromMolFile(fsdf, removeHs=False)

        for term, prmlist in self._terms.items():
            if term not in self.IMPLEMENTED_TERMS:
                logging.warning("prm term '%s' not implemented"%term)
                continue
            for smarts, prm in prmlist:
                rd_smarts = Chem.MolFromSmarts(smarts)
                matches = m1.GetSubstructMatches(rd_smarts)
                if not matches:
                    continue
                for match in matches:
                    prm_def = tuple(atypes[_] for _ in match)
                    prm_types[term].update({prm_def:prm})
        return prm_types
    def print_prm_types(self, prm_types):
        outlines = []
        for term, prmdict in prm_types.items():
            for atomtype_def, prm in prmdict.items():
                w = (term,) + tuple(atomtype_def) + tuple(prm)
                outlines.append(' '.join(map(str, w)) + '\n')
        return ''.join(outlines)

    def convert_prmtype_to_prmclass(self, prm_types, atomclass_map):
        """Convert atom types in prm def to atom classes
        """
        prm_classes = defaultdict(dict)
        for term, prmdict in prm_types.items():
            if term in self.ATOM_TYPE_TERMS:
                prm_classes[term] = prmdict
                continue
            for atomtype_def, prm in prmdict.items():
                atomclass_def = tuple(atomclass_map.get(_, _) for _ in atomtype_def)
                prm_classes[term].update({atomclass_def:prm})
        return prm_classes

    def format_prm_line(self, term, prm_def, prm, keyword=None):
        """ Get tinker-formatted prm line

        Args:
            term: the term to look up in the format dictionary
            prm_def: tuple of integer atom types/classes
            prm: list of floats
            keyword: the tinker keyword to be written to the key file. Default is the same as 'term'
        """
        if keyword is None:
            keyword = term
        data = tuple(prm_def) + tuple(prm)
        if term in self.PRM_FORMAT and len(data) < len(self.PRM_FORMAT[term]):
            outlist = [keyword] + [fmt%val for fmt,val in zip(self.PRM_FORMAT[term], data)]
        else:
            outlist = [keyword] + [str(_) for _ in data]
        return ' '.join(outlist) + '\n'

    def modify_key_file(self, ftxyz, inpkey, outkey, fsdf=None):
        '''Modify key file according to prm patch file

        Args:
            ftxyz: tinker xyz file
            inpkey: input key file
            outkey: output key file
            fsdf: sdf file for the ligand. If not provided, the structure
                will be read from xyz file, which will be less accurate
        '''
        aclass_map = read_key_atomclass(inpkey)
        # atom type def
        prm_types1 = self.match_atomtypes(ftxyz, fsdf)
        # atom type/class def
        prm_types = self.convert_prmtype_to_prmclass(prm_types1, aclass_map)
        mod_msg = self.print_prm_types(prm_types)
        if mod_msg:
            logging.info("Parameter modifications\n%s"%(mod_msg))
        outlines = []
        with open(inpkey, 'r') as fh:
            lines = fh.readlines()
            outlines = list(lines)

            iline_mpole = -5
            curr_type = None
            curr_prm = None
            curr_term = ''
            qpole_prm = []
            qpole_fmt = []
            for iline, line1 in enumerate(lines):
                line = line1.strip()
                w = re.findall('(\s+\S+)', line1)
                words = line.split()
                newline = line1
                if line.startswith('multipole'):
                    curr_type = (int(words[1]), )
                    term = 'quadrupole-scale'
                    if curr_type in prm_types[term]:
                        curr_term = term
                        curr_prm = prm_types[term][curr_type]
                        iline_mpole = iline
                        qpole_prm = []
                        qpole_fmt = []
                elif line.startswith('vdw '):
                    curr_type = (int(words[1]), )
                    if curr_type in prm_types['vdw']:
                        curr_prm = prm_types['vdw'][curr_type]
                        newline = self.format_prm_line('vdw', curr_type, curr_prm)
                    elif curr_type in prm_types['vdw-scale']:
                        curr_prm = prm_types['vdw-scale'][curr_type]
                        if len(curr_prm) == 1:
                            curr_prm = (curr_prm[0], 1)
                        old_prm = tuple(float(_) for _ in words[2:])
                        # scale the first parameter, duplicate the rest
                        new_prm = (old_prm[0]*curr_prm[0],old_prm[1]*curr_prm[1]) + tuple(old_prm[2:])
                        newline = self.format_prm_line('vdw', curr_type, new_prm)
                elif line.startswith('polarize'):
                    curr_type = (int(words[1]), )
                    if curr_type in prm_types['polarize']:
                        curr_prm = prm_types['polarize'][curr_type]
                        #old_prm = tuple(float(_) for _ in words[2:])
                        #new_prm = (curr_prm[0],) + tuple(old_prm[1:])
                        polarize_fmt = line[:(len('polarize'+w[0]))] + '%' + str(len(w[1])) + '.3f' + line[(len('polarize'+w[0]+w[1])):] + '\n'
                        newline = polarize_fmt%(curr_prm[0])
                elif iline <= iline_mpole + 4 and iline > iline_mpole + 1:
                    # number of of qpole parameters
                    n_pole = iline - iline_mpole - 1
                    try:
                        qpole_prm.append([float(_) for _ in w[:n_pole]])
                        qpole_fmt.append(''.join(['%'+str(len(_))+'.5f' for _ in w[:n_pole]]+w[n_pole:])+'\n')
                        #_scale = curr_prm[0]
                        #w_new = [('%' + str(len(_)) + '.5f')%(float(_)*_scale) for _ in w[:n_pole]] + w[n_pole:]
                        #newline = ''.join(w_new) + '\n'
                    except ValueError or IndexError as e:
                        logging.warning(f"Failed to read quadrupole at line {iline} of {inpkey}: '{line}'. Expected {n_pole} floats")
                    if iline == iline_mpole + 4:
                        _scale = curr_prm[0]
                        qpole_prm2 = [[np.round(_val*_scale, 5) for _val in _vals] for _vals in qpole_prm]
                        qpole_prm2[0][0] = -(qpole_prm2[1][1] + qpole_prm2[2][2])
                        newline = ''.join(_fmt%tuple(_val) for (_fmt, _val) in zip(qpole_fmt, qpole_prm2))
                    else:
                        newline = ''
                outlines[iline] = newline
        if len(outlines) == 0:
            return
        with open(outkey, 'w') as fh:
            fh.write(''.join(outlines))


def modify_key(xyzfile, keyfile, keyout, sdffile=None, inpfile=None, sm_list=None, x_list=None, term=''):
    pm = ParmMod()
    if inpfile is not None:
        if hasattr(inpfile, '__iter__') and not isinstance(inpfile, str):
            inplist = inpfile
        else:
            inplist = [inpfile]
        for f1 in inplist:
            pm.read_inpfile(f1)
    elif sm_list is not None and x_list is not None and term != '':
        pm.read_list(sm_list, x_list, term=term)
    else:
        logging.warning("No input provided")
        return
    pm.modify_key_file(xyzfile, keyfile, keyout, sdffile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='parmmod.py')
    parser.add_argument('prm', nargs='*', help='smarts scale_factor [smarts scale_factor ...]')
    parser.add_argument('-x', nargs=1, metavar='txyz', required=True, help='txyz')
    parser.add_argument('-s', nargs=1, metavar='sdf', required=False, help='sdf')
    parser.add_argument('-i', nargs='+', metavar='input', required=False, help='input file for list of smarts')
    parser.add_argument('-k', nargs=1, metavar='key', required=False, help='input key file smarts')
    parser.add_argument('-o', nargs=1, metavar='key2', required=False, help='output key file')
    parser.add_argument('-t', nargs=1, metavar='term', required=False, default=['vdw'], help='FF term to modify')

    v = (parser.parse_args())

    fxyz = v.x[0]
    fsdf = None if v.s is None else v.s[0]
    finpkey = None if v.k is None else v.k[0]
    foutkey = None if v.o is None else v.o[0]
    term = v.t[0]
    pm = ParmMod()
    if v.i is None:
        sm_list = v.prm[0::2]
        scale_list = [(float(_), ) for _ in v.prm[1::2]]
        pm.read_list(sm_list, scale_list, term=term)
    else:
        for inpfile in v.i:
            pm.read_inpfile(inpfile)
    if finpkey and foutkey:
        pm.modify_key_file(fxyz, finpkey, foutkey, fsdf)

