# scale factor for quadrupole
quadrupole-scale [$([OD2]([#6])[#6])]            0.70 # ether
quadrupole-scale [$([OH1D2][CD4])]               0.70 # non-aromatic hydroxyl
quadrupole-scale [$([NH2D3][CD4])]               0.80 # primary amine
quadrupole-scale [$([NH1D3]([CD4])[CD4])]        0.80 # secondary amine
quadrupole-scale [$([NH0D3]([CD4])([CD4])[CD4])] 0.80 # tertiary amine
quadrupole-scale [c]                             0.90 # aromatic carbon
quadrupole-scale [n]                             0.90 # aromatic nitrogen
# polarizability parameters
polarize         [n]                             1.073 # aromatic nitrogen. Use amoeba09 parameter
# updated vdw parameters
vdw [$([O]=[C][OD1])]             3.6000 0.0950 # carboxylate. r0 increased by 1.5% from amoebabio18 (3.550 0.095)
vdw [$([O]=[C]([#6])[#6,#7,OD2])] 3.4000 0.1120 # carbonyl. r0 increased by 3% from amoeba09 (3.300 0.112)
vdw [SD2]                         4.2000 0.3550 # thiol, sulfide, disulfide. r0 increased by 5% from amoeba09 (4.005 0.355)
vdw [n]                           3.7100 0.1260 # aromatic nitrogen. epsilon increased by 20% from amoeba09 (3.710 0.105)
