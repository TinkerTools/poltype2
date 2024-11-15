
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
  rdkitmol = Chem.MolFromMolFile(sdffile,removeHs=False)
  num_atoms = rdkitmol.GetNumAtoms() 
  for i in range(num_atoms):
    atom = rdkitmol.GetAtomWithIdx(i)
    if (atom.GetHybridization() == Chem.HybridizationType.SP2) and (g.degree[i+1] == 3):
      if str(atom2type[i+1]) not in sp2AtomTypes:
        sp2AtomTypes.append(str(atom2type[i+1]))
  return sp2AtomTypes

