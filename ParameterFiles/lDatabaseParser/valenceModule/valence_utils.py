
# valence utils to be used
try:
  from openbabel import pybel 
except:
  import pybel
import networkx as nx
from rdkit import Chem

def findSp2AtomTypes(txyz, sdffile):
  sp2AtomTypes = []
  atom2type = {}
  g = nx.Graph()
  nodes = []
  edges = []
  lines = open(txyz).readlines()
  for line in lines[1:]:
    d = line.split()
    if int(d[0]) not in atom2type:
      atom2type[int(d[0])] = int(d[5])
    if d[0] not in nodes: 
      nodes.append(int(d[0]))
    for c in d[6:]:
      s = sorted([int(d[0]), int(c)])
      if s not in edges:
        edges.append(s)
  g.add_nodes_from(nodes)
  g.add_edges_from(edges)
 
  # rdkit is more robust
  if sdffile.endswith('mol2'):
    rdkitmol = Chem.MolFromMol2File(sdffile,removeHs=False)
  else:
    rdkitmol = Chem.MolFromMolFile(sdffile,removeHs=False)
  num_atoms = rdkitmol.GetNumAtoms() 
  for i in range(num_atoms):
    atom = rdkitmol.GetAtomWithIdx(i)
    if (atom.GetHybridization() == Chem.HybridizationType.SP2) and (g.degree[i+1] == 3):
      if str(atom2type[i+1]) not in sp2AtomTypes:
        sp2AtomTypes.append(str(atom2type[i+1]))
  
  # Detect Nitrogen of Amide
  amide_nitrogens = [] 
  pattern = Chem.MolFromSmarts('[#6](=O)[#7X3]')
  matches = rdkitmol.GetSubstructMatches(pattern, uniquify=False)
  for match in matches:
    amide_nitrogens.append(match[2])
  
  amide_nitrogen_types = []
  for n in amide_nitrogens:
    atype = atom2type[n+1]
    if atype not in amide_nitrogen_types:
      amide_nitrogen_types.append(atype)

  # Detect Nitrogen of Aniline-like molecule
  aniline_like_nitrogens = [] 
  pattern = Chem.MolFromSmarts('[a][#7X3;h3]')
  matches = rdkitmol.GetSubstructMatches(pattern, uniquify=False)
  for match in matches:
    aniline_like_nitrogens.append(match[1])
  
  aniline_like_nitrogen_types = []
  for n in aniline_like_nitrogens:
    atype = atom2type[n+1]
    if atype not in aniline_like_nitrogen_types:
      aniline_like_nitrogen_types.append(atype)
  
  # N of Aniline-like molecule DONOT need opbend/anglep
  for n in aniline_like_nitrogen_types:
    if str(n) in sp2AtomTypes:
      sp2AtomTypes.remove(str(n))
      print(f'Removing type {str(n)} from SP2 atoms since it is an aniline_like_nitrogen')
  
  # N of Amide DO need opbend/anglep
  for n in amide_nitrogen_types:
    if str(n) not in sp2AtomTypes:
      sp2AtomTypes.append(str(n))
      print(f'Adding type {str(n)} to SP2 atoms since it is an amide_nitrogen')
  return sp2AtomTypes

