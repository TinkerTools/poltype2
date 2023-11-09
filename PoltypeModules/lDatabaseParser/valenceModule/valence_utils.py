
# valence utils to be used
try:
  from openbabel import pybel 
except:
  import pybel
import networkx as nx

def findSp2AtomTypes(txyz):
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
  
  for mol in pybel.readfile("txyz", txyz):
    natoms = len(mol.atoms)
    for i in range(natoms):
      hyb = mol.atoms[i].hyb
      if (hyb == 2) and (g.degree[i+1] == 3):
        if atom2type[i+1] not in sp2AtomTypes:
          sp2AtomTypes.append(str(atom2type[i+1]))
  return sp2AtomTypes

