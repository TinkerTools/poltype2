import numpy as np
import os, sys


def fchk_info(fname, num_atom):
  path = os.getcwd()
  log = "%s.log" %(fname)
  fchk = "%s.fchk" %(fname)
  bond_list = []
  angle_list = []
  if(not(os.path.exists(os.path.join(path, log))) or not(os.path.exists(os.path.join(path, fchk)))):
    raise FileError("log or fchk file does not exist")
    sys.exit(0)

  f_log = open(self.log)
  lines = f_log.readlines()
  f_log.close()

  i = 0
  j = 0
  for line in lines:
    if(j != 1 and len(line) > 80 and line[:81] == ' ! Name  Definition              Value          Derivative Info.                !'):
      i = 1
      j = 1
      continue

    if(i == 1):
      terms = line.split()
      if(len(terms) > 6 and terms[1][0] == 'R' and terms[-1] == '!'):
        bond_list.append([int(terms[2].split(',')[0][2:]), int(terms[2].split(',')[1][:-1])])
      elif(len(terms) > 6 and terms[1][0] == 'A' and terms[-1] == '!'):
        if(len(terms[2].split(',')) == 3):
          angle_list.append([int(terms[2].split(',')[0][2:]), int(terms[2].split(',')[1]), int(terms[2].split(',')[2][:-1])])
        else:
          angle_list.append([int(terms[2].split(',')[0][2:]), int(terms[2].split(',')[1]), int(terms[2].split(',')[2])])

  crd = os.popen("grep 'Current cartesian coordinates' %s" %(fchk)).read()
  crd_num = int(crd.split()[5])
  list_crd = os.popen("grep -A%d 'Current cartesian coordinates' %s" %(np.ceil(crd_num/5.), fchk)).read()
  coords = [float(i) for i in list_crd.split()[6:]]
  coords = np.array(coords) * 0.529  # Bohr to Ang
  coords = coords.reshape(num_atom, 3)

  hes = os.popen("grep 'Cartesian Force Constants' %s" %(fchk)).read()
  hes_num = int(hes.split()[5])
  list_hes = os.popen("grep -A%d 'Cartesian Force Constants' %s" %(np.ceil(hes_num/5.), self.fchk)).read()
  unprocessed_hessian = [float(i) for i in list_hes.split()[6:]]
  unprocessed_hessian = np.array(unprocessed_hessian) * 627.509391 / (np.square(0.529)) # Hartree/bohr to kcal/mol/ang
  hessian_mat = np.zeros((3*num_atom, 3*num_atom))
  for i in range(3*num_atom):
    for j in range(i+1):
      hessian_mat[i][j] = unprocessed_hessian[int(i*(i+1)/2+j)]
      hessian_mat[j][i] = unprocessed_hessian[int(i*(i+1)/2+j)]
  return bond_list, angle_list, coords, hessian_mat

def eigen(num_atom, hessian_mat):
  N = num_atom
  eigenvectors = np.zeros((N,N,3,3))
  eigenvalues = np.zeros((N,N,3))
  for i in range(N):
    for j in range(N):
      v1, v2 = np.linalg.eig(hessian_mat[(i*3):((i+1)*3), (j*3):((j+1)*3)])
      eigenvalues[i, j, :] = v1
      eigenvectors[i, j, :, :] = v2
  return eigenvalues, eigenvectors

def bond_projection(bond_list, coords, eigenvalues, eigenvectors, scale):
  # vectors for bond
  bond_num = len(bond_list)
  k_b = np.zeros(bond_num) # array for all force constants
  k_b_dict = {}
  for i in range(bond_num):
    bond_AB = bond_list[i]
    # coords index from 0
    vec_AB = coords[bond_AB[1]-1] - coords[bond_AB[0]-1]
    eigenvalues_AB = eigenvalues[bond_AB[0]-1, bond_AB[1]-1, :]
    eigenvectors_AB = eigenvectors[bond_AB[0]-1, bond_AB[1]-1, 0:3, 0:3]
    eigenvalues_BA = eigenvalues[bond_AB[1]-1, bond_AB[0]-1, :]
    eigenvectors_BA = eigenvectors[bond_AB[1]-1, bond_AB[0]-1, 0:3, 0:3]
    unit_vec_AB = vec_AB / b_a.bond_value[i]
    unit_vec_BA = - vec_AB / b_a.bond_value[i]

    for j in range(3):
      k_b[i] = k_b[i] - 0.5 * eigenvalues_AB[j] * np.abs(unit_vec_AB.dot(eigenvectors_AB[:,j]))
      k_b[i] = k_b[i] - 0.5 * eigenvalues_BA[j] * np.abs(unit_vec_BA.dot(eigenvectors_BA[:,j]))
    k_b[i] *= 0.5 * (scale**2)
    k_b_dict[bond_list[i]] = k_b[i]

  return k_b_dict

def angle_projection(angle_list, coords, eigenvalues, eigenvectors, scale):
  # vector for angle
  angle_num = len(angle_list)
  k_a = np.zeros((angle_num,))
  k_a_dict = {}
  p_dict = {}
  p_list_A = []
  p_list_C = []
  bond_AB = []
  bond_CB = []

  for i in range(angle_num):
    angle_ABC = angle_list[i]
    tmp_1 = "%d-%d" %(angle_ABC[1], angle_ABC[0]) # atom B - atom A
    tmp_2 = "%d-%d" %(angle_ABC[1], angle_ABC[2]) # atom B - atom C
    p_dict[tmp_1] = []
    p_dict[tmp_2] = []

  for i in range(angle_num):
    angle_ABC = angle_list[i]
    vec_AB = coords[angle_ABC[0]-1] - coords[angle_ABC[1]-1]
    vec_CB = coords[angle_ABC[2]-1] - coords[angle_ABC[1]-1]
    unit_AB = vec_AB / np.linalg.norm(vec_AB)
    unit_CB = vec_CB / np.linalg.norm(vec_CB)

    unit_N = np.cross(unit_CB, unit_AB)
    unit_N = unit_N / np.linalg.norm(unit_N)

    unit_PA = np.cross(unit_AB, unit_N)
    unit_PC = np.cross(unit_CB, unit_N)
    unit_PA = unit_PA / np.linalg.norm(unit_PA)
    unit_PC = unit_PC / np.linalg.norm(unit_PC)

    p_list_A.append(unit_PA)
    p_list_C.append(unit_PC)

    tmp_1 = "%d-%d" %(angle_ABC[1], angle_ABC[0]) # atom B - atom A
    tmp_2 = "%d-%d" %(angle_ABC[1], angle_ABC[2]) # atom B - atom C
    p_dict[tmp_1].append(unit_PA)
    p_dict[tmp_2].append(unit_PC)
    bond_AB.append(np.linalg.norm(vec_AB))
    bond_CB.append(np.linalg.norm(vec_CB))

  # modified force constant for angles
  for i in range(angle_num):
    angle_ABC = angle_list[i]
    sum_1 = np.zeros((3,))
    sum_2 = np.zeros((3,))
    eigenvalues_AB = eigenvalues[angle_ABC[0]-1, angle_ABC[1]-1, :]
    eigenvalues_CB = eigenvalues[angle_ABC[2]-1, angle_ABC[1]-1, :]
    eigenvectors_AB = eigenvectors[angle_ABC[0]-1, angle_ABC[1]-1, 0:3, 0:3]
    eigenvectors_CB = eigenvectors[angle_ABC[2]-1, angle_ABC[1]-1, 0:3, 0:3]
    unit_PA = p_list_A[i]
    unit_PC = p_list_C[i]
    for j in range(3):
      eig_AB = eigenvectors_AB[:,j]
      sum_1[j] = eigenvalues_AB[j] * abs(unit_PA.dot(eig_AB.T))
      eig_CB = eigenvectors_CB[:,j]
      sum_2[j] = eigenvalues_CB[j] * abs(unit_PC.dot(eig_CB.T))
    sum_1 = np.sum(sum_1)
    sum_2 = np.sum(sum_2)

    scale_A = 0.
    scale_C = 0.
    # scaling factor
    tmp_1 = "%d-%d" %(angle_ABC[1], angle_ABC[0]) # atom B - atom A
    if(len(p_dict[tmp_1]) == 1):
      scale_A = 1.
    else:
      for j in p_dict[tmp_1]:
        scale_A += np.dot(unit_PA, j.T) ** 2
      scale_A =  1 + ((scale_A - 1) / (len(p_dict[tmp_1]) - 1))
    tmp_2 = "%d-%d" %(angle_ABC[1], angle_ABC[2]) # atom B - atom C
    if(len(p_dict[tmp_2]) == 1):
      scale_C = 1.
    else:
      for j in p_dict[tmp_2]:
        scale_C += np.dot(unit_PC, j.T) ** 2
      scale_C =  1 + ((scale_C - 1) / (len(p_dict[tmp_2]) - 1))

    sum_1 = sum_1/scale_A
    sum_2 = sum_2/scale_C
    k_theta = 1. / ((bond_AB[i] ** 2) * sum_1) + 1. / ((bond_CB[i] ** 2) * sum_2)
    k_theta = 1./k_theta
    k_theta *= -0.5 * (scale**2)
    k_a[i] = k_theta
    k_a_dict[angle_ABC] = k_a[i]

  return k_a_dict

def modified_Seminario(mode, comb, ttypes, dict_):
  if(mode == 'b'):
    terms = comb.split('_')
    atom_1 = int(terms[0]) - int(ttypes[0]) + 1
    atom_2 = int(terms[1]) - int(ttypes[0]) + 1
    if([atom_1, atom_2] in dict_.keys()):
      return dict_[[atom_1, atom_2]]
    elif([atom_2, atom_1] in dict_.keys()):
      return dict_[[atom_2, atom_1]]
  elif(mode == 'a'):
    terms = comb.split('_')
    atom_1 = int(terms[0]) - int(ttypes[0]) + 1
    atom_2 = int(terms[1]) - int(ttypes[0]) + 1
    atom_3 = int(terms[2]) - int(ttypes[0]) + 1
    if([atom_1, atom_2, atom_3] in dict_.keys()):
      return dict_[[atom_1, atom_2, atom_3]]
    elif([atom_3, atom_2, atom_1] in dict_.keys()):
      return dict_[[atom_3, atom_2, atom_1]]
