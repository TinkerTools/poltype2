


class PySCF_init_setup():

    def __init__(self,poltype_obj,cure_dir,COM_file,mol,\
                 modred,b_a_restraint,skip_error,charge,\
                 opt_conv,torsion_restraints):

        self.pol_obj = poltype_obj
        self.mol_obj = mol
        self.init_data = {'topdir': cure_dir, 'COM_file': COM_file,
                          'modred': modred,'b_a_rest': b_a_restraint,
                          'skip_error': skip_error,'opt_conv': opt_conv,
                          'tor_rest': torsion_restraints}
        self.mol_data_init = {'charge': charge}


    

    def read_Gauss_inp_coord(self):

        self.mol_data_init['Atoms_list'] = []
        self.mol_data_init['Atoms_coord'] = []
        with open(f"{self.init_data['topdir']}/{self.init_data['COM_file']}", 'r') as f:
            for lines in f:
                L = lines.strip('\n').split()
                if len(L)==4 and '#' not in L:
                    self.mol_data_init['Atoms_list'].append(L[0])
                    self.mol_data_init['Atoms_coord'].append([float(i) for i in L[1:]])


        return

    def write_xyz(self,name):

        self.PySCF_inp_xyz = f"{self.init_data['COM_file'].strip('.com')}_{name}.xyz"

        with open(f"{self.init_data['topdir']}/{self.PySCF_inp_xyz}", 'w') as f:
            f.write(f"{len(self.mol_data_init['Atoms_coord'])}\n\n")
            for A,C in zip(self.mol_data_init['Atoms_list'],self.mol_data_init['Atoms_coord']):
                f.write(f'{A:3s} {C[0]:9.6f} {C[1]:9.6f} {C[2]:9.6f}\n')

        return

    def write_geometric_tor_const(self,name):

        self.PySCF_inp_const = f"PySCF_{name}_constraint.txt"

        with open(f"{self.init_data['topdir']}/{self.PySCF_inp_const}", 'w') as f:
            f.write('$set\n')
            for tor in self.init_data['tor_rest']:
                angle = self.mol_obj.GetTorsion(tor[0], tor[1], tor[2], tor[3]) % 360
                f.write(f'dihedral {tor[0]} {tor[1]} {tor[2]} {tor[3]} {angle:3.4f}\n')

        return


    def write_PySCF_input(self):

        self.PySCF_inp_file = f"{self.init_data['COM_file'].strip('.com')}.py" 
        self.PySCF_out_file = f"{self.init_data['COM_file'].strip('.com')}_pyscf.out" 

        self.write_PySCF_header()

        self.write_PySCF_params()

        self.write_PySCF_main()

    def write_PySCF_header(self):

        with open(f"{self.init_data['topdir']}/{self.PySCF_inp_file}", 'w') as f:
            f.write('import pyscf\n')
            f.write('from pyscf import gto, scf\n')
            f.write('from pyscf.dft import rks,uks\n')
            f.write('from pyscf import dft\n')
            f.write('from pyscf.geomopt.geometric_solver import optimize\n\n')

        return

    def write_PySCF_params(self):

        with open(f"{self.init_data['topdir']}/{self.PySCF_inp_file}", 'a+') as f:
            f.write(f'params = {{"constraints": "{self.PySCF_inp_const}",\n')
            if self.pol_obj.optconvergence == 'LOOSE':
                f.write('            "convergence_energy": 1e-4,\n') 
                f.write('            "convergence_grms": 1.7e-3,\n') 
                f.write('            "convergence_gmax": 2.5e-3,\n') 
                f.write('            "convergence_drms": 6.7e-3,\n') 
                f.write('            "convergence_dmax": 1.0e-2,\n') 

            f.write('}\n\n') 
        return 


    def write_PySCF_main(self):

        atm_s = f'atom="{self.PySCF_inp_xyz}"'
        basis_s = f'basis = "{self.pol_obj.optbasisset}"'
        spin_s = f'spin = {int(self.mol_obj.GetTotalSpinMultiplicity())-1}'
        if self.mol_data_init['charge'] == None:
            chg_s = f'charge = {self.mol_obj.GetTotalCharge()}'
        else:
            chg_s = f'charge = {self.mol_data_init["charge"]}'
        unit_s = 'unit = "A"'

        with open(f"{self.init_data['topdir']}/{self.PySCF_inp_file}", 'a+') as f:
            f.write(f'mol = pyscf.M({atm_s}, {basis_s}, {spin_s}, {chg_s}, {unit_s})\n')
            f.write('mol.build()\n')
            f.write(f'mf = pyscf.solvent.PCM(mol.UKS(xc = "{self.pol_obj.pyscf_opt_meth}").density_fit())\n')
            f.write(f'mf.with_solvent.method = "{self.pol_obj.pyscf_sol_imp}"\n')
            f.write(f'mf.with_solvent.eps = {self.pol_obj.pyscf_sol_eps}\n\n')
            f.write(f'mol_eq = optimize(mf, maxsteps={self.pol_obj.optmaxcycle}, **params)\n')
            f.write('print(f"Final optimized structure:")\n')
            f.write('print(mol_eq.tostring())\n\n')


    

