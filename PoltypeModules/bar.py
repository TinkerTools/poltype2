import submitjobs
import tables
import terminate as term
import os
import sys
import time
import numpy as np
import csv
import restraints
import pandas as pd
import pymbar
import warnings

def convert_1d(arr):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: convert nx2 array to nx1 array of the difference
  """

  if len(arr.shape) == 1:
    arr1 = arr
  elif len(arr.shape) ==2 and arr.shape[1] == 2:
    arr1 = arr[:, 0] - arr[:, 1]
  else:
    arr1 = arr[:, 0]
  return arr1
  
def calc_eff_size(arr, equil=False):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  arr1 = convert_1d(arr)
  if len(arr1) == 0:
    return 0, 1
  NBLOCK = 71
  nskip = max(1, len(arr1)//NBLOCK)
  t0 = 0
  indices = np.arange(len(arr1), dtype=int)
  if equil:
    [t0, g, Neff_max] = pymbar.timeseries.detectEquilibration(arr1, nskip=nskip)
    #indices = pymbar.timeseries.subsampleCorrelatedData(arr1[t0:], g=g)
  else:
    try:
        indices = pymbar.timeseries.subsampleCorrelatedData(arr1)
        length=len(indices)
    except:
        length=len(arr1)
    g = len(arr1)/max(1, length)
  return t0, g

def resample(arr, n, uniform=False):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  if uniform:
    idx = np.array((np.arange(n, dtype=int)*(len(arr)-1))//(n-1), dtype=int)
    return arr[idx]
  else:
    return arr[np.random.randint(0, len(arr), n, dtype=int)]
    
def subsample1(arr1, equil=False, corr=True, skip=0, last=None, n_aug=1, uniform=True):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  '''
  equil: discard equilibration samples
  corr: compute correlation
  skip: first sample
  last: last sample
  n_aug: multiply sample size by this factor
  '''
  kwargs = {'equil':equil}
  if skip is None:
    skip = 0
  if last is None:
    last = len(arr1)
  elif last < 0:
    last = last % len(arr1)
  if (skip-1)*skip<=0 and (last-1)*last<=0:
    skip = int(skip*len(arr1))
    last = int(last*len(arr1))
  last = min(last, len(arr1))

  arr1b = arr1[skip:last]
  t1, g1 = calc_eff_size(arr1b, **kwargs)
  if not corr:
    g1 = 1
  n1 = len(arr1b)*n_aug//max(1, g1)
  arr1c = resample(arr1b, n1, uniform=True)
  return arr1c, t1 + skip, last, g1

def subsample_arrays(arr1s, **kwargs):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  arr2s = []
  # t0, t1, g, n
  tgn = np.zeros((len(arr1s), 4))
  for i, arr1 in enumerate(arr1s):
    res = subsample1(arr1, **kwargs)
    arr2s.append(res[0])
    tgn[i, 0] = res[1]
    tgn[i, 1] = res[2]
    tgn[i, 2] = res[3]
    tgn[i, 3] = len(res[0])

  arr2 = np.concatenate(arr2s, axis=0)
  t1 = np.sum(tgn[:, 0])
  t2 = np.sum(tgn[:, 1])
  g1 = np.sum(tgn[:, 3] * tgn[:, 2])/max(1, np.sum(tgn[:, 3]))
  return arr2, t1, t2, g1

def subsample(arr, equil=False, corr=True):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  arr1 = convert_1d(arr)
  NBLOCK = 71
  nskip = max(1, len(arr1)//NBLOCK)
  t0 = 0
  indices = np.arange(len(arr1), dtype=int)
  if equil and corr:
    [t0, g, Neff_max] = pymbar.timeseries.detectEquilibration(arr1, nskip=nskip)
    indices = pymbar.timeseries.subsampleCorrelatedData(arr1[t0:], g=g)
  elif corr:
    indices = pymbar.timeseries.subsampleCorrelatedData(arr1)
  return [_+t0 for _ in indices]

def get_index_summary(indices):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  return '(%d,%d,%d)'%(indices[0], indices[-1], indices[1]-indices[0])

def get_index_summary2(idx1, idx2):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  msg = 'A:%s B:%s'%(get_index_summary(idx1), get_index_summary(idx2))
  return msg

def read_tinker_bars(inplist, temp=298, press=1.0, ibegin=0, iend=-1):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  u1s = []
  u2s = []
  temps = []
  for inp in inplist:
    res = read_tinker_bar(inp)
    if res is not None:
      res0 = res[0][ibegin:iend]
      res1 = res[1][ibegin:iend]
      #if len(res[0][ibegin:iend]) + len(res[1][ibegin:iend]) > 0:
      if len(res0) + len(res1) > 0:
        u1s.append(res0)
        u2s.append(res1)
        temps.append(res[2])
  if len(u1s) == 0:
    return None
  return u1s, u2s, temps

def read_tinker_bar(inp, temp=298, press=1.0):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  #beta = 1.0/(8.314*temp/4184)
  # bar * ang^3 in kcal/mol
  PV_UNIT_CONV = 1e5*1e-30/4184*6.02e23
  with open(inp, 'r') as fh:
    lines = fh.readlines()
    if len(lines) == 0:
      return
    w = lines[0].split()
    if len(w) == 0 or (not w[0].isdigit()):
      return
    n1 = int(w[0])
    temp1 = float(w[1]) if len(w) >= 2 else temp
    w = lines[n1+1].split()
    if len(w) == 0 or (not w[0].isdigit()):
      return
    n2 = int(w[0])
    temp2 = float(w[1]) if len(w) >= 2 else temp

    beta1 = 1.0/(8.314*temp1/4184)
    beta2 = 1.0/(8.314*temp2/4184)

    if len(lines) != n1+n2+2:
      print("nlines (%d) != %d + %d"%(len(lines), n1, n2))
      return

    arr1 = np.fromstring(''.join(lines[1:n1+1]), sep=' ').reshape(n1, -1)
    arr2 = np.fromstring(''.join(lines[n1+2:]), sep=' ').reshape(n2, -1)

    u1 = calc_tinker_reduced(arr1, beta=beta1, press=press)
    u2 = calc_tinker_reduced(arr2, beta=beta2, press=press)
    #return (np.concatenate((u1[idx1], u2[idx2]), axis=0).transpose(), [len(idx1), len(idx2)], msg)
    return u1, u2, temp1

def calc_tinker_pv(arr1, press=1.0):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  PV_UNIT_CONV = 1e5*1e-30/4184*6.02e23
  if arr1.shape[1] >= 4:
    pv1 = press*arr1[:, 3:4]*PV_UNIT_CONV
  else:
    pv1 = 0
  return pv1

def calc_tinker_reduced(arr1, beta, press=1.0):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  pv1 = calc_tinker_pv(arr1, press=press)
  u1 = beta*(arr1[:, 1:3] + pv1)
  return u1

def concat_arr(arr1, arr2, idx1, idx2):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  msg = get_index_summary2(idx1, idx2)
  return (np.concatenate((arr1[idx1], arr2[idx2]), axis=0).transpose(), [len(idx1), len(idx2)], msg)

def tinker_to_mbar(arr1, arr2, equil=False, corr=True):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  assert len(arr1.shape) == 2 and arr1.shape[1] == 2
  idx1 = subsample(arr1[:, 1] - arr1[:, 0], equil=equil, corr=corr)
  idx2 = subsample(arr2[:, 1] - arr2[:, 0], equil=equil, corr=corr)
  return concat_arr(arr1, arr2, idx1, idx2)

def arr_to_mbar(arr1, arr2):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  return (np.concatenate((arr1, arr2), axis=0).transpose(), [len(arr1), len(arr2)])

def exp_ave(ener, return_sd=True, ener_unit=1):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  if len(ener) == 0:
    return np.nan, 0
  emin = np.min(ener)
  eave = -np.log(np.mean(np.exp(-(ener-emin)))) + emin
  if return_sd:
    _sd = np.std(ener)
    return eave*ener_unit, _sd*ener_unit
  return eave*ener_unit

def calc_mbar(u_kn, N_k, teff=1, temp=298):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  #data = concat_arr(arr1, arr2, idx1, idx2)
  #cols = 'dF(kcal/mol) se_dF(kcal/mol) dH se_dH TdS se_TdS dF_fwd dF_bwd dE_fwd dE_bwd sd_dE_fwd sd_dE_bwd p_A(traj_B) p_B(traj_A) eig_overlap'.split()
  cols = 'dF(kcal/mol) se_dF(kcal/mol) dH TdS se_dH dF_fwd dF_bwd dE_fwd dE_bwd sd_dE_fwd sd_dE_bwd overlap'.split()
  if u_kn is None:
    return cols
  kt = temp*8.314/4184
  mbar = pymbar.MBAR(u_kn, N_k)
  #results = mbar.getFreeEnergyDifferences()
  results = mbar.computeEntropyAndEnthalpy()
  overlap = mbar.computeOverlap()
  #msg1 = get_index_summary(idx1)
  #msg2 = get_index_summary(idx2)
  es_fwd = (u_kn[1, 0:N_k[0]] - u_kn[0, 0:N_k[0]])
  es_bwd = (u_kn[0, N_k[0]:sum(N_k[0:2])] - u_kn[1, N_k[0]:sum(N_k[0:2])])
  de_fwd = np.mean(es_fwd)*kt
  de_bwd = np.mean(es_bwd)*kt

  fwd, sd_fwd = exp_ave(es_fwd, ener_unit=kt)
  bwd, sd_bwd = exp_ave(es_bwd, ener_unit=kt)
  ghs = [] # free energy, enthalpy, entropy
  for i in range(3):
    ghs.append(kt*results[i*2][0, 1])
    ghs.append(kt*results[i*2+1][0, 1]*np.sqrt(teff))
  #return ghs[0], ghs[1], ghs[2], ghs[3], ghs[4], ghs[5], fwd, -bwd, de_fwd, -de_bwd, sd_fwd, sd_bwd, overlap['matrix'][0, 1], overlap['matrix'][1, 0], overlap['scalar']
  return ghs[0], ghs[1], ghs[2], ghs[4], ghs[3], fwd, -bwd, de_fwd, -de_bwd, sd_fwd, sd_bwd, overlap['scalar']


def get_bar_res(data, ishift=0, **kwargs):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  '''
  data: [array1, array2, temp]
  return 1 row of data frame
  '''
  col0 = 'start_A end_A start_B end_B g_A g_B'.split()
  if data is None:
    return col0 + calc_mbar(None, None)
  opt1 = kwargs
  arr1, t1, t1b, g1 = subsample_arrays(data[0], **opt1)
  arr2, t2, t2b, g2 = subsample_arrays(data[1], **opt1)
  if len(arr1) + len(arr2) == 0:
    print("WARNING: empty array")
    return np.nan
  u_kn, N_k = arr_to_mbar(arr1, arr2)
  if len(u_kn) == 0:
    print("WARNING: empty array")
    return np.nan
  res = calc_mbar(u_kn, N_k, temp=data[2][0])
  return [_+ishift for _ in [t1, t1b, t2, t2b]] + [g1, g2] + list(res)


def compute_total_sd(df1):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  ''' compute sum on df columns
  sqrt(sum of squares) for 'sd_XXX', 'se_XXX'
  mininum for 'overlap*'
  '''
  arr1 = df1.sum(axis=0)
  for i, col in enumerate(df1.columns):
    if col.startswith('se_') or col.startswith('sd_'):
      arr1.iloc[i] = np.sqrt(np.sum(df1[col]**2.0))
    elif col.startswith('overlap') or col.endswith('overlap'):
      arr1.iloc[i] = df1[col].min()
    elif col.startswith('g_'):
      arr1.iloc[i] = df1[col].mean()
  return arr1

def compute_block_ave(df1):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  ''' compute average on df columns
  if df contains se_XXX and XXX, se_XXX will be calculated by the stderr of XXX
  '''
  arr1 = df1.mean(axis=0)
  nrow = len(df1.index)
  for i, col in enumerate(df1.columns):
    if col.startswith('se_') and col[3:] in df1.columns:
      arr1.iloc[i] = (df1[col[3:]].std())/np.sqrt(nrow)
  return arr1

def sum_df_list(dflist):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  ''' compute sum for list of dfs
  sqrt(sum of squares) for 'sd_XXX', 'se_XXX'
  mininum for 'overlap*' '*overlap'
  '''
  if len(dflist) == 0:
    return None
  elif len(dflist) == 1:
    return dflist[0]

  idx0 = dflist[0].index
  for df1 in dflist[1:]:
    #idx0 = idx0 & df1.index
    idx0 = idx0.intersection(df1.index)
  df2list = []
  df3 = pd.DataFrame(columns=dflist[0].columns)
  for row in idx0:
    df2 = pd.concat([_.loc[[row], :] for _ in dflist], axis=0, ignore_index=True)
    df3.loc[row, :] = (compute_total_sd(df2))
  return df3

def check_data_size(data, ibegin, iend):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  data1 = [[], [], []]
  for i in range(len(data[0])):
    if len(data[0][i][ibegin:iend]) + len(data[1][i][ibegin:iend]) > 0:
      for j in range(len(data)):
        data1[j].append(data[j][i])
  return data1

def get_summary(df1, verbose=True):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  name_val = []
  name_val.append(('dF_equ', df1.loc['equ', 'dF(kcal/mol)']))
  name_val.append(('dF_all', df1.loc['all', 'dF(kcal/mol)']))
  if 'blockave_all' in df1.index:
    name_val.append(('se_dF_block', df1.loc['blockave_all', 'se_dF(kcal/mol)']))
  else:
    name_val.append(('se_dF_block', np.nan))
  nA = (df1.loc['equ', 'end_A']-df1.loc['equ', 'start_A'])//df1.loc['equ', 'g_A']
  nB = (df1.loc['equ', 'end_B']-df1.loc['equ', 'start_B'])//df1.loc['equ', 'g_B']

  sd_de_fwd = df1.loc['all', 'sd_dE_fwd']
  sd_de_bwd = df1.loc['all', 'sd_dE_bwd']

  name_val.append(('se_dF_bar', df1.loc['equ', 'se_dF(kcal/mol)']))
  name_val.append(('dF_fwd', df1.loc['all', 'dF_fwd']))
  name_val.append(('dF_bwd', df1.loc['all', 'dF_bwd']))
  name_val.append(('ddF_fwd_bwd', abs(df1.loc['all', 'dF_fwd']-df1.loc['all', 'dF_bwd'])))
  name_val.append(('dH', df1.loc['equ', 'dH']))
  name_val.append(('se_dH', df1.loc['equ', 'se_dH']))
  name_val.append(('sd_dE', np.sqrt(np.sum(np.array([sd_de_fwd, sd_de_bwd])**2.0))))
  name_val.append(('overlap_all', df1.loc['all', 'overlap']))
  name_val.append(('overlap_equ', df1.loc['equ', 'overlap']))
  name_val.append(('nA', nA))
  name_val.append(('nB', nB))
  names, vals = zip(*name_val)
  df2 = pd.DataFrame([vals], columns=names, index=['SUMMARY'])
  CHECK_THR = (('se_dF_block', 0.5, 1), ('se_dF_bar', 0.5, 1), ('ddF_fwd_bwd', 1.0, 1), ('sd_dE', 2.0, 1), ('overlap_all', 0.3, -1))
  name_desc = {}
  name_desc.update({'dF_equ':'Free energy difference (equilibrated) (kcal/mol)'})
  name_desc.update({'dF_all':'Free energy difference (all data) (kcal/mol) '})
  name_desc.update({'se_dF_block':'Standard error (block average) (kcal/mol) '})
  name_desc.update({'ddF_fwd_bwd':'Difference between forward and backward FEP (kcal/mol) '})
  name_desc.update({'sd_dE':'Standard deviation of energy difference (kcal/mol)'})
  name_desc.update({'overlap_all':'Connectivity between two states'})
  if verbose:
    for (_row, _thr, _sign) in CHECK_THR:
      if _row not in df2.columns:
        continue
      val = df2.loc['SUMMARY', _row] 
      if (val - _thr)*_sign > 0:
        _name = name_desc[_row] if _row in name_desc else _row
        if _sign > 0:
          _direct = 'above'
        else:
          _direct = 'below'
        #print("WARNING: %s %s %s (%.3f vs. %.3f)"%(_name,  _direct, 'threshold', val, _thr))
  return df2

def calc_dg(flist, ibegin=0, iend=-1, nblocks=5, summary=True, verbose=True):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  NBLOCK = nblocks
  NMIN = 100
  names = 'all uncorr equ'.split()
  opts = [{'equil':False, 'corr':False}, 
          {'equil':False, 'corr':True},
          {'equil':True,  'corr':True},
          ]

  df_out = pd.DataFrame(columns=get_bar_res(None))
  namesb = names[:1]
  optsb = opts[:1]
  df_blocks = [pd.DataFrame(columns=get_bar_res(None)) for _ in namesb]
  n_name_blocks = 1
  data0 = read_tinker_bars(flist, ibegin=ibegin, iend=iend)
  if data0 is None:
    print("ERROR: empty data")
    return
  teff = 1
  opt1 = {'skip':0, 'last':None, 'uniform':True, "ishift":ibegin}
  for opt, name in zip(opts, names):
    opt1.update(opt)
    res = get_bar_res(data0, **opt1)
    if res is not None:
      df_out.loc[name, :] = res

  if len(data0[0]) == 1:
    # calculate results of blocks
    #n1 = len(data0[0][0][ibegin:iend])
    #n2 = len(data0[1][0][ibegin:iend])
    n1 = len(data0[0][0])
    n2 = len(data0[1][0])
    if min(n1, n2) >= NBLOCK*NMIN:
      #bn1 = n1//NBLOCK
      #bn2 = n2//NBLOCK
      iblk = -1
      for opt, name1 in zip(optsb, namesb):
        iblk += 1
        for i in range(NBLOCK):
          name = "block%d%s"%(i+1, name1)
          opt1.update({'skip':i/NBLOCK, 'last':(i+1)/NBLOCK})
          opt1.update(opt)
          res = get_bar_res(data0, **opt1)
          if res is not None:
            #df_out.loc[name, :] = res
            df_blocks[iblk].loc[name, :] = res
        break
  elif len(data0[0]) > 1:
    #opt1 = {'skip':ibegin, 'last':iend, 'uniform':True}
    opt1 = {'skip':0, 'last':None, 'uniform':True, "ishift":ibegin}
    iblk = -1
    for opt, name1 in zip(optsb, namesb):
      iblk += 1
    # calculate results of files
      for i in range(len(data0[0])):
        name = "file%d%s"%(i+1, name1)
        opt1.update(opt)
        data1 = [data0[0][i:i+1], data0[1][i:i+1], data0[2][i:i+1]]
        res = get_bar_res(data1, **opt1)
        if res is not None:
          df_blocks[iblk].loc[name, :] = res

  for name, df2 in zip(namesb, df_blocks):
    if len(df2.index) > 1:
      df_out.loc['blockave_%s'%name, :] = compute_block_ave(df2)
  df_out = pd.concat([df_out]+df_blocks)
  #if verbose:
  #  print(df_out.to_string(max_rows=len(df_out.index), max_cols=len(df_out.columns)))
  df_sum = get_summary(df_out,verbose)
  #if summary:
    #print()
    #print(df_sum.to_string(max_rows=len(df_sum.index), max_cols=len(df_sum.columns)))
  return df_out, df_sum, (namesb, df_blocks)

def add_line(inp):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  return '\n'+'*'*8+' '+inp+' '+'*'*8+'\n'

def list_to_string(alist, maxl=12):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  if len(alist) == 0:
    return '[]'
  def trim(s):
    if len(s) > maxl:
      return '..'+s[-maxl:]
    else:
      return s

  return '[%s ... %s]'%(trim(str(alist[0])), trim(str(alist[-1])))

def calc_dg_summary(flist, seperate=False, outfile=None, **kwargs):
  """
  Intent:
  Input:
  Output:
  Referenced By: 
  Description: 
  """
  '''
  wrapper for calc_dg

  if seperate, analyze each file independently
  otherwise combine all files
  '''
  if seperate and len(flist) > 1:
    dflist = []
    df1list = []
    df_block_list = []
    names_block = []
    if outfile is None:
      fout = sys.__stdout__
    else:
      fout = open(outfile, 'w')
    for f1 in flist:
      #print(add_line('Analysis of %s'%f1))
      df1, df1sum, data_block = calc_dg([f1], **kwargs)
      dflist.append(df1sum)
      df1list.append(df1)
      df_block_list.append(data_block[1])
      names_block = data_block[0]
    df_blocks = [sum_df_list(_) for _ in zip(*df_block_list)]
    #print(add_line('Sum of %d files %s'%(len(flist), list_to_string(flist))))
    df3 = sum_df_list(df1list)
    for name, df_block in zip(names_block, df_blocks):
        df3.loc['blockave_%s'%name, :] = compute_block_ave(df_block)
    #print(df3.to_string())
    df2 = pd.concat(dflist, axis=0)
    #df2 = pd.DataFrame(df2.to_numpy(), index=flist, columns=df2.columns)
    df2.set_axis(flist, axis=0, inplace=True)
    #df2.reindex(flist)
    #df2.loc['Total', :] = compute_total_sd(df2)
    df2.loc['Total', :] = get_summary(df3, verbose=False).iloc[0, :].to_numpy()
    for col in 'nA nB'.split():
      if col in (df2.columns & df3.columns):
        df2.loc['Total', col] = df2.iloc[:-1, col].sum()

    #print(add_line('SUMMARY'))
    #print(df2.to_string(), file=fout)
    if outfile is not None:
      fout.close()
  else:
    df1, df1sum, _ = calc_dg(flist, **kwargs)
  matrix_sum = df2.to_numpy() 
  dF_equ_list=matrix_sum[0:len(flist),0]   
  dF_all_list=matrix_sum[0:len(flist),1] 
  se_dF_block_list=matrix_sum[0:len(flist),2] 
  se_dF_bar_list=matrix_sum[0:len(flist),3]   
  dF_fwd_list=matrix_sum[0:len(flist),4]   
  dF_bwd_list=matrix_sum[0:len(flist),5] 
  ddF_fwd_bwd_list=matrix_sum[0:len(flist),6]       
  dH_list=matrix_sum[0:len(flist),7]    
  se_dH_list=matrix_sum[0:len(flist),8]    
  sd_dE_list=matrix_sum[0:len(flist),9] 
  overlap_all_list=matrix_sum[0:len(flist),10] 
  overlap_equ_list=matrix_sum[0:len(flist),11]       
  nA_list=matrix_sum[0:len(flist),12]      
  nB_list=matrix_sum[0:len(flist),13]
  
  return dF_equ_list,dF_all_list,se_dF_block_list,se_dF_bar_list,dF_fwd_list,dF_bwd_list,ddF_fwd_bwd_list,dH_list,se_dH_list,sd_dE_list,overlap_all_list,overlap_equ_list,nA_list,nB_list


def GrabBarFileName(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    files=os.listdir()
    for f in files:
        if '.' in f:
            fsplit=f.split('.')
            ext=fsplit[1]
            if 'bar' in f:
                barfile=f 
    barfilepath=os.path.join(os.getcwd(),barfile)
    return barfilepath

def ComputeThermoProperties(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    curdir=os.getcwd()
    
    poltype.freeenergyconv=[]
    poltype.freeenergyerrorconv=[]
    for i in range(len(poltype.baroutputfilepath)):
        totalframes=int(poltype.proddynframenum[i])
        divisor=int(totalframes/poltype.barinterval)
        framechunks=list(np.arange(divisor,totalframes,divisor))
        framechunks.append(None)
        baroutputfilepathlist=poltype.baroutputfilepath[i]
        freeenergylistoflist=[]
        freeenergyerrorlistoflist=[]
        enthalpylistoflist=[]
        enthalpyerrorlistoflist=[]
        entropylistoflist=[]
        entropyerrorlistoflist=[]
        freeenergyfwdlistoflist=[]
        freeenergybwdlistoflist=[]
        totalfreeenergy=[]
        totalenthalpy=[]
        totalentropy=[]
        totalerrorfreeenergy=[]
        totalerrorenthalpy=[]
        totalerrorentropy=[]
        totalfreeenergyfwd=[]
        totalfreeenergybwd=[]
        overlaplistoflist=[]
        totalfreeenergyconv=[]
        totalerrorfreeenergyconv=[]
        poltype.freeenergyconv.append([])
        poltype.freeenergyerrorconv.append([])
        for j in range(len(baroutputfilepathlist)):
            totalfreeenergyconv.append([])
            totalerrorfreeenergyconv.append([])
            baroutputfilepath=baroutputfilepathlist[j]
            barfiles=[]
            for path in baroutputfilepath:
                head,tail=os.path.split(path)
                os.chdir(head)
                barfilepath=GrabBarFileName(poltype) 
                barfiles.append(barfilepath) 
            for framenumidx in range(len(framechunks)):
                framenum=framechunks[framenumidx]
                verbose=False
                if framenumidx==len(framechunks)-1:
                    verbose=True
                dF_equ_list,dF_all_list,se_dF_block_list,se_dF_bar_list,dF_fwd_list,dF_bwd_list,ddF_fwd_bwd_list,dH_list,se_dH_list,sd_dE_list,overlap_all_list,overlap_equ_list,nA_list,nB_list=calc_dg_summary(barfiles, seperate=True, outfile='table.txt', ibegin=0, iend=framenum, nblocks=5, summary=True, verbose=verbose)
                if framenumidx==len(framechunks)-1:
                    overlaplistoflist.append(overlap_all_list)
                    freeenergylistoflist.append(dF_all_list)
                    freeenergyerrorlistoflist.append(se_dF_bar_list)
                    enthalpylistoflist.append(dH_list)
                    enthalpyerrorlistoflist.append(se_dH_list)
                    TdS_list=np.array(dH_list)-np.array(dF_all_list)
                    entropylistoflist.append(TdS_list)
                    freeenergyfwdlistoflist.append(dF_fwd_list)
                    freeenergybwdlistoflist.append(dF_bwd_list)
                    dF_all=np.sum(dF_all_list)
                    totalfreeenergy.append(dF_all)
                    dH=np.sum(dH_list)
                    totalenthalpy.append(dH)
                    valuelist=np.add(np.square(se_dH_list),np.square(se_dF_bar_list))
                    valuelist=valuelist.astype(float)
                    se_TdS_list=np.sqrt(valuelist)
                    entropyerrorlistoflist.append(se_TdS_list)
                    TdS=(dF_all)-(dH)
                    totalentropy.append(TdS)
                    dF_fwd=np.sum(dF_fwd_list)
                    dF_bwd=np.sum(dF_bwd_list)
                    totalfreeenergyfwd.append(dF_fwd)
                    totalfreeenergybwd.append(dF_bwd)
                    se_dF_bar=np.sqrt(np.sum(np.square(se_dF_bar_list)))
                    se_dH=np.sqrt(np.sum(np.square(se_dH_list)))
                    totalerrorfreeenergy.append(se_dF_bar)
                    totalerrorenthalpy.append(se_dH)
                    se_TdS=np.sqrt(se_dH**2+se_dF_bar**2)
                    totalerrorentropy.append(se_TdS)
                dF_all=np.sum(dF_all_list)
                se_dF_bar=np.sqrt(np.sum(np.square(se_dF_bar_list)))
                totalfreeenergyconv[j].append(dF_all)
                totalerrorfreeenergyconv[j].append(se_dF_bar)
        poltype.overlaplist[i]=overlaplistoflist
        poltype.freeenergylist[i]=freeenergylistoflist
        poltype.freeenergy[i]=sum(totalfreeenergy)
        poltype.freeenergyerrorlist[i]=freeenergyerrorlistoflist
        poltype.enthalpylist[i]=enthalpylistoflist
        poltype.enthalpy[i]=sum(totalenthalpy)
        poltype.enthalpyerrorlisttotal[i]=enthalpyerrorlistoflist
        poltype.entropylist[i]=entropylistoflist
        poltype.entropy[i]=sum(totalentropy)
        poltype.entropyerrorlisttotal[i]=entropyerrorlistoflist
        poltype.freeenergylistfwd[i]=freeenergyfwdlistoflist
        poltype.freeenergyfwd[i]=sum(totalfreeenergyfwd)
        poltype.freeenergylistbwd[i]=freeenergybwdlistoflist
        poltype.freeenergybwd[i]=sum(totalfreeenergybwd)
        poltype.enthalpyerror[i]=np.sqrt(np.sum(np.square(totalerrorenthalpy)))
        poltype.entropyerror[i]=np.sqrt(np.sum(np.square(totalerrorentropy)))
        poltype.freeenergyerror[i]=np.sqrt(np.sum(np.square(totalerrorfreeenergy)))
        totalfreeenergyconvsum=[]
        totalerrorfreeenergyconvsum=[]
        totalfreeenergyconv=np.transpose(np.array(totalfreeenergyconv))
        totalerrorfreeenergyconv=np.transpose(np.array(totalerrorfreeenergyconv))

        for idx in range(len(totalfreeenergyconv)):
            freeg=sum(totalfreeenergyconv[idx])
            freegerrls=totalerrorfreeenergyconv[idx]
            Sum=0
            for val in freegerrls:
                Sum+=val**2
            freegerr=np.sqrt(Sum)
            totalfreeenergyconvsum.append(freeg)
            totalerrorfreeenergyconvsum.append(freegerr)
        poltype.freeenergyconv[i]=totalfreeenergyconvsum
        poltype.freeenergyerrorconv[i]=totalerrorfreeenergyconvsum


        if poltype.solvation==True and poltype.complexation==False:
            poltype.tabledict[i][u'ΔGˢᵒˡᵛ']=poltype.freeenergy[i]
            poltype.tabledict[i][u'ΔGˢᵒˡᵛᵉʳʳ']=poltype.freeenergyerror[i]
            poltype.tabledict[i][u'ΔHˢᵒˡᵛ']=poltype.enthalpy[i]
            poltype.tabledict[i][u'ΔHˢᵒˡᵛᵉʳʳ']=poltype.enthalpyerror[i]
            poltype.tabledict[i][u'ΔSˢᵒˡᵛ']=poltype.entropy[i]
            poltype.tabledict[i][u'ΔSˢᵒˡᵛᵉʳʳ']=poltype.entropyerror[i]
        elif poltype.complexation==True and poltype.solvation==False:
            poltype.tabledict[i][u'ΔGᶜᵒᵐᵖᵘⁿᶜᵒʳʳ']=poltype.freeenergy[i]
            poltype.tabledict[i][u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ']=poltype.freeenergyerror[i]
            poltype.tabledict[i][u'ΔHᶜᵒᵐᵖ']=poltype.enthalpy[i]
            poltype.tabledict[i][u'ΔHᶜᵒᵐᵖᵉʳʳ']=poltype.enthalpyerror[i]
            poltype.tabledict[i][u'ΔSᶜᵒᵐᵖ']=poltype.entropy[i]
            poltype.tabledict[i][u'ΔSᶜᵒᵐᵖᵉʳʳ']=poltype.entropyerror[i]
        elif poltype.complexation==True and poltype.solvation==True:
            if i==0:
                poltype.tabledict[i][u'ΔGᶜᵒᵐᵖᵘⁿᶜᵒʳʳ']=poltype.freeenergy[i]
                poltype.tabledict[i][u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ']=poltype.freeenergyerror[i]
                poltype.tabledict[i][u'ΔHᶜᵒᵐᵖ']=poltype.enthalpy[i]
                poltype.tabledict[i][u'ΔHᶜᵒᵐᵖᵉʳʳ']=poltype.enthalpyerror[i]
                poltype.tabledict[i][u'ΔSᶜᵒᵐᵖ']=poltype.entropy[i]
                poltype.tabledict[i][u'ΔSᶜᵒᵐᵖᵉʳʳ']=poltype.entropyerror[i]
            elif i==1:
                poltype.tabledict[i][u'ΔGˢᵒˡᵛ']=poltype.freeenergy[i]
                poltype.tabledict[i][u'ΔGˢᵒˡᵛᵉʳʳ']=poltype.freeenergyerror[i]
                poltype.tabledict[i][u'ΔHˢᵒˡᵛ']=poltype.enthalpy[i]
                poltype.tabledict[i][u'ΔHˢᵒˡᵛᵉʳʳ']=poltype.enthalpyerror[i]
                poltype.tabledict[i][u'ΔSˢᵒˡᵛ']=poltype.entropy[i]
                poltype.tabledict[i][u'ΔSˢᵒˡᵛᵉʳʳ']=poltype.entropyerror[i]

    tables.WriteTableUpdateToLog(poltype,verbose=False)
    comptable=[]
    solvtable=[]

    
    vdwlambdapairlist,vdwlambdapairlistoflist=FlattenListOfListGenerateBARPairs(poltype,poltype.vdwlambdascheme)
    elelambdapairlist,elelambdapairlistoflist=FlattenListOfListGenerateBARPairs(poltype,poltype.estatlambdascheme)
    comptable.append(['Vdw-Lambda']+vdwlambdapairlist)
    comptable.append(['Ele-Lambda']+elelambdapairlist)
    comptable.append(['Rest-Lambda']+FlattenListOfListGenerateBARPairs(poltype,poltype.restlambdascheme)[0])
    solvtable.append(['Vdw-Lambda']+vdwlambdapairlist)
    solvtable.append(['Ele-Lambda']+elelambdapairlist)
    solvtable.append(['Rest-Lambda']+FlattenListOfListGenerateBARPairs(poltype,poltype.restlambdascheme)[0])

    if poltype.solvation==True and poltype.complexation==False:
        barpairs=[]
        for ls in poltype.lambdafolderlist[0]:
            pairs=GenerateBARPairs(poltype,ls)
            barpairs.extend(pairs)
        solvtable.append(['Folder']+barpairs)
        solvfreeenergylist=FlattenListOfList(poltype,poltype.freeenergylist[0])
        totalvdw,totalele,vdwfreeenergynoioncorrection,elefreeenergynoioncorrection,totalfreeenergynoioncorrection,vdwfreeenergynogas,elefreeenergynogas,totalfreeenergynogas,vdwfreeenergygas,elefreeenergygas,totalfreeenergygas=GenerateEleVdwSums(poltype,solvfreeenergylist,vdwlambdapairlist,elelambdapairlist,baroutputfilepathlist)
        fwdenergylist=FlattenListOfList(poltype,poltype.freeenergylistfwd[0])
        bwdenergylist=FlattenListOfList(poltype,poltype.freeenergylistbwd[0])
        fwdbwddiff=list(np.abs(np.array(fwdenergylist)-np.array(bwdenergylist)))
        fwdbwddiff=[str(i) for i in fwdbwddiff]
        CheckIfForwardBackwardPertubationsAreSimilar(poltype,fwdenergylist,bwdenergylist)
        poltype.tabledict[0][u'ΔGˢᵒˡᵛᵉˡᵉ']=totalele
        poltype.tabledict[0][u'ΔGˢᵒˡᵛᵛᵈʷ']=totalvdw
        poltype.tabledict[0][u'ΔG']=totalfreeenergynoioncorrection
        poltype.tabledict[0][u'ΔGᵉˡᵉ']=elefreeenergynoioncorrection
        poltype.tabledict[0][u'ΔGᵛᵈʷ']=vdwfreeenergynoioncorrection
        poltype.tabledict[0][u'ΔGᵍᵃˢ']=totalfreeenergygas
        poltype.tabledict[0][u'ΔGˢᵒˡ']=totalfreeenergynogas
        poltype.tabledict[0][u'ΔGᵉˡᵉᵍᵃˢ']=elefreeenergygas
        poltype.tabledict[0][u'ΔGᵉˡᵉˢᵒˡ']=elefreeenergynogas
        poltype.tabledict[0][u'ΔGᵛᵈʷᵍᵃˢ']=vdwfreeenergygas
        poltype.tabledict[0][u'ΔGᵛᵈʷˢᵒˡ']=vdwfreeenergynogas
        solvtable.append([u'ΔGˢᵒˡᵛ']+solvfreeenergylist)
        solvtable.append([u'ΔGˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlist[0]))
        solvtable.append([u'ΔGˢᵒˡᵛᶠʷᵈ']+fwdenergylist)
        solvtable.append([u'ΔGˢᵒˡᵛᵇʷᵈ']+bwdenergylist)
        solvtable.append([u'ΔGˢᵒˡᵛᶠʷᵈᵇʷᵈ']+fwdbwddiff)
        solvtable.append([u'ΔHˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.enthalpylist[0]))
        solvtable.append([u'ΔHˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.enthalpyerrorlisttotal[0]))
        solvtable.append([u'TΔSˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.entropylist[0]))
        solvtable.append([u'TΔSˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.entropyerrorlisttotal[0]))
        solvtable.append(['SolvOverlap']+FlattenListOfList(poltype,poltype.overlaplist[0]))
        CheckOverLapValues(poltype,poltype.overlaplist[0],barpairs)
        newsolvtable=list(map(list, zip(*solvtable)))
        tempname='SolvBARResults.csv'
        with open(poltype.outputpath+tempname, mode='w') as energy_file:
            energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for row in newsolvtable:
                energy_writer.writerow(row)

    elif poltype.complexation==True and poltype.solvation==True:
        barpairs=[]
        for ls in poltype.lambdafolderlist[1]:
            pairs=GenerateBARPairs(poltype,ls)
            barpairs.extend(pairs)
        solvtable.append(['Folder']+barpairs)
        compfreeenergylist=FlattenListOfList(poltype,poltype.freeenergylist[0])
        totalvdw,totalele,vdwfreeenergynoioncorrection,elefreeenergynoioncorrection,totalfreeenergynoioncorrection,vdwfreeenergynogas,elefreeenergynogas,totalfreeenergynogas,vdwfreeenergygas,elefreeenergygas,totalfreeenergygas=GenerateEleVdwSums(poltype,compfreeenergylist,vdwlambdapairlist,elelambdapairlist,baroutputfilepathlist)
        poltype.tabledict[0][u'ΔGᶜᵒᵐᵖᵉˡᵉ']=totalele
        poltype.tabledict[0][u'ΔGᶜᵒᵐᵖᵛᵈʷ']=totalvdw
        compfwdenergylist=FlattenListOfList(poltype,poltype.freeenergylistfwd[0])
        compbwdenergylist=FlattenListOfList(poltype,poltype.freeenergylistbwd[0])
        CheckIfForwardBackwardPertubationsAreSimilar(poltype,compfwdenergylist,compbwdenergylist)
        solvfwdenergylist=FlattenListOfList(poltype,poltype.freeenergylistfwd[1])
        solvbwdenergylist=FlattenListOfList(poltype,poltype.freeenergylistbwd[1])
        solvfwdbwddiff=list(np.abs(np.array(solvfwdenergylist)-np.array(solvbwdenergylist)))
        solvfwdbwddiff=[str(i) for i in solvfwdbwddiff]
        compfwdbwddiff=list(np.abs(np.array(compfwdenergylist)-np.array(compbwdenergylist)))
        compfwdbwddiff=[str(i) for i in compfwdbwddiff]
        CheckIfForwardBackwardPertubationsAreSimilar(poltype,solvfwdenergylist,solvbwdenergylist)
        solvfreeenergylist=FlattenListOfList(poltype,poltype.freeenergylist[1])
        totalvdw,totalele,vdwfreeenergynoioncorrection,elefreeenergynoioncorrection,totalfreeenergynoioncorrection,vdwfreeenergynogas,elefreeenergynogas,totalfreeenergynogas,vdwfreeenergygas,elefreeenergygas,totalfreeenergygas=GenerateEleVdwSums(poltype,solvfreeenergylist,vdwlambdapairlist,elelambdapairlist,baroutputfilepathlist)
        poltype.tabledict[1][u'ΔGˢᵒˡᵛᵉˡᵉ']=totalele
        poltype.tabledict[1][u'ΔGˢᵒˡᵛᵛᵈʷ']=totalvdw

        if poltype.addgas==True:
            solvlists=poltype.freeenergylist[1]
            liq=solvlists[0]
            gas=solvlists[1]
            empty=[None for i in gas]
            inter=liq+gas[::-1]
            inter=[round(i,2) for i in inter]
            inter=[str(i) for i in inter]+empty

        solvtable.append([u'ΔGˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.freeenergylist[1]))
        if poltype.addgas==True:
            solvtable.append([u'ΔGˢᵒˡᵛⁱⁿᵗᵉʳ']+inter)
        solvtable.append([u'ΔGˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlist[1]))
        solvtable.append([u'ΔGˢᵒˡᵛᶠʷᵈ']+solvfwdenergylist)
        solvtable.append([u'ΔGˢᵒˡᵛᵇʷᵈ']+solvbwdenergylist)
        solvtable.append([u'ΔGˢᵒˡᵛᶠʷᵈᵇʷᵈ']+solvfwdbwddiff)

        solvtable.append([u'ΔHˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.enthalpylist[1]))
        solvtable.append([u'ΔHˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.enthalpyerrorlisttotal[1]))
        solvtable.append([u'TΔSˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.entropylist[1]))
        solvtable.append([u'TΔSˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.entropyerrorlisttotal[1]))
        solvtable.append(['SolvOverlap']+FlattenListOfList(poltype,poltype.overlaplist[1]))
        CheckOverLapValues(poltype,poltype.overlaplist[1],barpairs)
        barpairs=[]
        for ls in poltype.lambdafolderlist[0]:
            pairs=GenerateBARPairs(poltype,ls)
            barpairs.extend(pairs)
        if poltype.addgas==True:
            complists=poltype.freeenergylist[0]
            liq=solvlists[0]
            gas=solvlists[1]
            empty=[None for i in gas]
            inter=liq+gas[::-1]
            inter=[round(i,2) for i in inter]
            inter=[str(i) for i in inter]+empty
        comptable.append(['Folder']+barpairs)
        comptable.append([u'ΔGᶜᵒᵐᵖ']+compfreeenergylist)
        if poltype.addgas==True:
            comptable.append([u'ΔGᶜᵒᵐᵖⁱⁿᵗᵉʳ']+inter)
        comptable.append([u'ΔGᶜᵒᵐᵖᵉʳʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlist[0]))
        comptable.append([u'ΔGᶜᵒᵐᵖᶠʷᵈ']+compfwdenergylist)
        comptable.append([u'ΔGᶜᵒᵐᵖᵇʷᵈ']+compbwdenergylist)
        comptable.append([u'ΔGᶜᵒᵐᵖᶠʷᵈᵇʷᵈ']+compfwdbwddiff)
        comptable.append([u'ΔHᶜᵒᵐᵖ']+FlattenListOfList(poltype,poltype.enthalpylist[0]))
        comptable.append([u'ΔHᶜᵒᵐᵖᵉʳʳ']+FlattenListOfList(poltype,poltype.enthalpyerrorlisttotal[0]))
        comptable.append([u'TΔSᶜᵒᵐᵖ']+FlattenListOfList(poltype,poltype.entropylist[0]))
        comptable.append([u'TΔSᶜᵒᵐᵖᵉʳʳ']+FlattenListOfList(poltype,poltype.entropyerrorlisttotal[0]))
        comptable.append(['CompOverlap']+FlattenListOfList(poltype,poltype.overlaplist[0]))
        CheckOverLapValues(poltype,poltype.overlaplist[0],barpairs)
        newcomptable=list(map(list, zip(*comptable)))
        newsolvtable=list(map(list, zip(*solvtable)))

        tempname='CompBARResults.csv'
        with open(poltype.outputpath+tempname, mode='w') as energy_file:
            energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for row in newcomptable:
                energy_writer.writerow(row)

        tempname='SolvBARResults.csv'
        with open(poltype.outputpath+tempname, mode='w') as energy_file:
            energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for row in newsolvtable:
                energy_writer.writerow(row)

    os.chdir(curdir)



def CheckOverLapValues(poltype,overlaplist,barpairs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tol=.3
    for i in range(len(overlaplist)):
        valuelist=overlaplist[i]
        barpair=barpairs[i]
        for value in valuelist:
            if float(value)<tol:
                string='WARNING: Overlap between neighboring states is too small '+str(value)+' neighboring states are '+barpair+' tolerance is '+str(tol)
                warnings.warn(string)
                poltype.WriteToLog(string)


def CheckIfForwardBackwardPertubationsAreSimilar(poltype,fwdenergylist,bwdenergylist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tol=1
    for i in range(len(fwdenergylist)):
        fwdenergy=fwdenergylist[i]
        bwdenergy=bwdenergylist[i]
        diff=np.abs(fwdenergy-bwdenergy)
        if diff>tol:

            string='Forward and Backward pertubations for have detected pertubation greather than '+str(tol)+' kcal/mol '+'difference of '+str(diff)+' kcal/mol was detected'
            warnings.warn(string)
            poltype.WriteToLog(string)

def DecomposeLambdaArray(poltype,lambdals):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indiceschange=[]
    indicesconstant=[]
    most_common = max(lambdals, key = lambdals.count)
    for lamidx in range(len(lambdals)):
        lam=lambdals[lamidx]
        if lam==most_common:
            indicesconstant.append(lamidx)
        else:
            indiceschange.append(lamidx)

    return indiceschange,indicesconstant


def GenerateEleVdwSums(poltype,freeenergylist,vdwlambdapairlist,elelambdapairlist,baroutputfilepathlist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indiceslist=[]
    count=0
    eleboundaryindices=[]
    vdwboundaryindices=[]
    lenprevarray=0
    for i in range(len(poltype.vdwlambdascheme)):
        vdwlambdals=poltype.vdwlambdascheme[i]
        elelambdals=poltype.estatlambdascheme[i]
        elelambdapairlist=[]
        vdwlambdapairlist=[]
        for j in range(len(vdwlambdals)-1):
            vdwlambda=vdwlambdals[j]
            nextvdwlambda=vdwlambdals[j+1]
            elelambda=elelambdals[j]
            nextelelambda=elelambdals[j+1]
            elelambdapairlist.append(str(elelambda)+'-'+str(nextelelambda))
            vdwlambdapairlist.append(str(vdwlambda)+'-'+str(nextvdwlambda))
        eleindiceschange,eleindicesconstant=DecomposeLambdaArray(poltype,elelambdapairlist)
        vdwindiceschange,vdwindicesconstant=DecomposeLambdaArray(poltype,vdwlambdapairlist)
        eleindiceschange=[k+lenprevarray for k in eleindiceschange]
        vdwindiceschange=[k+lenprevarray for k in vdwindiceschange]
        eleboundaryindices.append([eleindiceschange[0],eleindiceschange[-1]])
        vdwboundaryindices.append([vdwindiceschange[0],vdwindiceschange[-1]])
        lenprevarray=len(vdwlambdapairlist)+lenprevarray
          
    vdwenergies=[]
    eleenergies=[]
    for i in range(len(vdwboundaryindices)):
        vdwindices=vdwboundaryindices[i]
        eleindices=eleboundaryindices[i]
        vdw1,vdw2=vdwindices[:]
        ele1,ele2=eleindices[:]
        vdwenergyls=freeenergylist[vdw1:vdw2+1]
        eleenergyls=freeenergylist[ele1:ele2+1]
        vdwenergies.append(vdwenergyls)
        eleenergies.append(eleenergyls)
    vdwenergiessum=[sum(ls) for ls in vdwenergies]
    eleenergiessum=[sum(ls) for ls in eleenergies]
    if len(vdwenergiessum)==3 and poltype.annihilatevdw==True:
        vdwfreeenergynoioncorrection=sum(vdwenergiessum[:1])
        elefreeenergynoioncorrection=sum(eleenergiessum[:1])
        totalfreeenergynoioncorrection=elefreeenergynoioncorrection+vdwfreeenergynoioncorrection


    elif len(vdwenergiessum)==2 and poltype.annihilatevdw==False:
        vdwfreeenergynoioncorrection=vdwenergiessum[0]
        elefreeenergynoioncorrection=eleenergiessum[0]
        totalfreeenergynoioncorrection=elefreeenergynoioncorrection+vdwfreeenergynoioncorrection
    else:
        vdwfreeenergynoioncorrection=None
        elefreeenergynoioncorrection=None
        totalfreeenergynoioncorrection=None

    vdwfreeenergynogas=vdwenergiessum[0]
    elefreeenergynogas=eleenergiessum[0]
    totalfreeenergynogas=elefreeenergynogas+vdwfreeenergynogas
    if len(vdwenergiessum)>=2:
        vdwfreeenergygas=vdwenergiessum[1]
        elefreeenergygas=eleenergiessum[1]
        totalfreeenergygas=elefreeenergygas+vdwfreeenergygas
    else:
        vdwfreeenergygas=None
        elefreeenergygas=None
        totalfreeenergygas=None


    totalvdw=sum(vdwenergiessum)
    totalele=sum(eleenergiessum)
    return totalvdw,totalele,vdwfreeenergynoioncorrection,elefreeenergynoioncorrection,totalfreeenergynoioncorrection,vdwfreeenergynogas,elefreeenergynogas,totalfreeenergynogas,vdwfreeenergygas,elefreeenergygas,totalfreeenergygas


def ExecuteBARSecondOption(poltype,joblistname):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    jobtolog={}
    jobtojobpath={}
    jobtoinputfilepaths={}
    jobtooutputfiles={}
    jobtoabsolutebinpath={}

    for j in range(len(poltype.simfoldname)):
        simfoldname=poltype.simfoldname[j]
        os.chdir(poltype.outputpath+simfoldname)
        proddynarcboxfilename=poltype.proddynarcboxfilename[j]
        proddynboxfilename=poltype.proddynboxfilename[j]
        for k in range(len(poltype.thermooutputfilepath[j])):
            thermooutputfilepath=poltype.thermooutputfilepath[j][k]
            barfilepathlist=poltype.barfilepath[j][k]
            for i in range(len(thermooutputfilepath)): 
                outputfilepath=thermooutputfilepath[i]
                barfilepath=barfilepathlist[i]
                path=os.path.split(outputfilepath)[0]
                head,tail=os.path.split(outputfilepath)
                if 'Gas' in tail: 
                    barpath=poltype.barpath
                else:
                    barpath=poltype.truebarpath

                maxframenum=int(os.path.getsize(path+'/'+proddynarcboxfilename)/os.path.getsize(path+'/'+proddynboxfilename))
                firstframenum=1
                cmdstr=BARSecondOptionCommand(poltype,barfilepath,firstframenum,maxframenum,outputfilepath,barpath)
                terminate,deletefile,error=term.CheckFileTermination(poltype,outputfilepath)
                if terminate==False:
                    jobtolog[cmdstr]=poltype.outputpath+joblistname
                    jobtojobpath[cmdstr]=path
                    head,tail=os.path.split(outputfilepath)
                    outputfilenames=[tail]
                    absbinpath=poltype.which(barpath)
                    inputfilepaths=[poltype.outputpath+poltype.prmfilepath,barfilepath]
                    jobtoinputfilepaths[cmdstr]=inputfilepaths
                    jobtooutputfiles[cmdstr]=outputfilenames
                    jobtoabsolutebinpath[cmdstr]=absbinpath
        os.chdir('..')
    if len(jobtolog.keys())!=0:
        poltype.WriteToLog('Computing free energy from .bar files ',prin=True)
        submitjobs.SubmitJobs(poltype,jobtolog,jobtojobpath,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,poltype.outputpath+joblistname)

def BARSecondOptionCommand(poltype,barfilepath,firstframenum,maxframenum,outputname,barpath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    head,tail=os.path.split(outputname)
    cmdstr=barpath+' '+'2'+' '+barfilepath+' '+str(firstframenum)+' '+str(maxframenum)+' '+'1'+' '+str(firstframenum)+' '+str(maxframenum)+' '+'1'+' > '+tail
    return cmdstr

def ExecuteBAR(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    jobtolog={}
    jobtojobpath={}
    jobtoinputfilepaths={}
    jobtooutputfiles={}
    jobtoabsolutebinpath={}
    jobtooutputpath={}

    for j in range(len(poltype.simfoldname)):
        simfoldname=poltype.simfoldname[j]
        os.chdir(poltype.outputpath+simfoldname)
        for k in range(len(poltype.baroutputfilepath[j])):
            baroutputfilepath=poltype.baroutputfilepath[j][k]
            secondarcpaths=poltype.secondarcpaths[j][k]
            firstarcpaths=poltype.firstarcpaths[j][k]
            for i in range(len(baroutputfilepath)): 
                outputfilepath=baroutputfilepath[i]
                path=os.path.split(outputfilepath)[0]
                head,tail=os.path.split(outputfilepath)
                if 'Gas' in tail:
                    barpath=poltype.barpath
                else:
                    if poltype.externalapi!=None:
                        barpath=poltype.barommpath
                    else:
                        barpath=poltype.truebarpath
                secondarcpath=secondarcpaths[i]
                firstarcpath=firstarcpaths[i]
                cmdstr=BARCommand(poltype,secondarcpath,firstarcpath,outputfilepath,barpath)
                terminate,deletefile,error=term.CheckFileTermination(poltype,outputfilepath)
                if poltype.redobar==True:
                    terminate=False 
                if terminate==False:
                    jobtolog[cmdstr]=poltype.outputpath+poltype.barjobsfilename
                    jobtojobpath[cmdstr]=path
                    secondarcname=os.path.split(secondarcpath)[1]
                    firstarcname=os.path.split(firstarcpath)[1]
                    secondkeypath=secondarcpath.replace('.arc','.key')
                    firstkeypath=firstarcpath.replace('.arc','.key')
                    barname=secondarcname.replace('.arc','.bar')
                    inputfilepaths=[secondarcpath,firstarcpath,poltype.outputpath+poltype.prmfilepath,secondkeypath,firstkeypath]
                    head,tail=os.path.split(outputfilepath)
                    outputfilenames=[tail,barname]
                    absbinpath=poltype.which(barpath)
                    jobtoinputfilepaths[cmdstr]=inputfilepaths
                    jobtooutputfiles[cmdstr]=outputfilenames
                    jobtoabsolutebinpath[cmdstr]=absbinpath
                    jobtooutputpath[cmdstr]=head
                    curdir=os.getcwd()
                    os.chdir(path)
                    files=os.listdir()
                    for f in files:
                        if '.bar' in f or '.err' in f or '.end' in f:
                            thepath=os.path.join(os.getcwd(),f)
                            poltype.WriteToLog('Removing file '+thepath,prin=True)
                            os.remove(f) # sometimes BAR fails and want to delete bad .bar files
                        if "BAR" in f and '.out' in f and poltype.redobar==True:
                            os.remove(f)
                    os.chdir(curdir)
        os.chdir('..')
    if len(jobtolog.keys())!=0:
        poltype.WriteToLog('Running BAR',prin=True)
        submitjobs.SubmitJobs(poltype,jobtolog,jobtojobpath,jobtoinputfilepaths,jobtooutputfiles,jobtoabsolutebinpath,poltype.outputpath+poltype.barjobsfilename,jobtooutputpath)

def BARCommand(poltype,secondarcpath,firstarcpath,outputfilepath,barpath):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    head,tail=os.path.split(outputfilepath)
    cmdstr=barpath+' '+'1'+' '+secondarcpath+' '+str(poltype.equilibriatescheme[-1])+' '+firstarcpath+' '+str(poltype.equilibriatescheme[-1])+' N'+' > '+tail
    return cmdstr

def SumTheFreeEnergyStepsFromBAR(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for j in range(len(poltype.simfoldname)):
        simfoldname=poltype.simfoldname[j]
        os.chdir(poltype.outputpath+simfoldname)
        freeenergylistoflist=[]
        freeenergyerrorlistoflist=[]
        freeenergyviabootstraplistoflist=[]
        freeenergyviabootstraperrorlistoflist=[]
        enthalpylistoflist=[]
        enthalpyerrorlistoflist=[]
        entropylistoflist=[]
        entropyerrorlistoflist=[]
        freeenergyfwdlistoflist=[]
        freeenergyfwderrorlistoflist=[]
        freeenergybwdlistoflist=[]
        freeenergybwderrorlistoflist=[]
        freeenergyviabariterlistoflist=[]
        freeenergyviabaritererrorlistoflist=[]

        for k in range(len(poltype.thermooutputfilepath[j])):
            thermooutputfilepath=poltype.thermooutputfilepath[j][k]
            totalfreeenergylist=[]
            totalfreeenergyviabootstraplist=[]
            totalenthalpylist=[]
            totalentropylist=[]
            totalfreeenergyfwdlist=[]
            totalfreeenergybwdlist=[]


            freeenergylist=[]
            freeenergyerrorlist=[]
            freeenergyviabootstraplist=[]
            freeenergyviabootstraperrorlist=[]
            enthalpylist=[]
            enthalpyerrorlist=[]
            entropylist=[]
            entropyerrorlist=[]
            freeenergyfwdlist=[]
            freeenergyfwderrorlist=[]
            freeenergybwdlist=[]
            freeenergybwderrorlist=[]
            freeenergyviabariterlist=[]
            freeenergyviabaritererrorlist=[]


            for i in range(len(thermooutputfilepath)):
                outputfilepath=thermooutputfilepath[i]
                temp=open(outputfilepath,'r')
                foundfreeenergyline=False
                freeenergy=0
                freeenergyviabootstrap=0
                enthalpy=0
                entropy=0
                freeenergyfwd=0
                freeenergybwd=0
                foundbariter=False
                for line in temp.readlines():
                    if 'Free Energy via BAR Iteration' in line and foundbariter==False:
                        foundbariter=True
                        foundfreeenergyline=True
                        linesplit=line.split()
                        freeenergyvalue=float(linesplit[5])
                        index=line.find('+/-')
                        newstring=line[index+3:]
                        newstringlinesplit=newstring.split()
                        freenergyerror=float(newstringlinesplit[0])
                        freeenergy+=freeenergyvalue
                        freeenergylist.append(freeenergyvalue)
                        freeenergyerrorlist.append(freenergyerror)
                        freeenergyviabariterlist.append(freeenergyvalue)
                        freeenergyviabaritererrorlist.append(freenergyerror)
                        

                    elif 'Free Energy via BAR Bootstrap' in line: 

                        foundfreeenergyline=True
                        linesplit=line.split()
                        freeenergyvalue=float(linesplit[5])
                        index=line.find('+/-')
                        newstring=line[index+3:]
                        newstringlinesplit=newstring.split()
                        freenergyerror=float(newstringlinesplit[0])
                        freeenergyviabootstrap+=freeenergyvalue
                        freeenergyviabootstraplist.append(freeenergyvalue)
                        freeenergyviabootstraperrorlist.append(freenergyerror)

                    elif 'Enthalpy via BAR Estimate' in line:
                        linesplit=line.split()
                        enthalpyvalue=float(linesplit[4])
                        index=line.find('+/-')
                        newstring=line[index+3:]
                        newstringlinesplit=newstring.split()
                        enthalpyerror=float(newstringlinesplit[0])
                        enthalpy+=enthalpyvalue
                        enthalpylist.append(enthalpyvalue)
                        enthalpyerrorlist.append(enthalpyerror)
                        

                    elif 'Entropy via BAR Estimate' in line:
                        linesplit=line.split()
                        entropyvalue=float(linesplit[4])
                        entropyerror=np.sqrt((enthalpyerror)**2+(freenergyerror)**2)
                        entropy+=entropyvalue
                        entropylist.append(entropyvalue)
                        entropyerrorlist.append(entropyerror)
                        

                    elif 'Free Energy via Forward FEP' in line:
                        linesplit=line.split()
                        freeenergyvalue=float(linesplit[5])
                        index=line.find('+/-')
                        newstring=line[index+3:]
                        newstringlinesplit=newstring.split()
                        freenergyerror=float(newstringlinesplit[0])
                        freeenergyfwd+=freeenergyvalue
                        
                        freeenergyfwdlist.append(freeenergyvalue)
                        freeenergyfwderrorlist.append(freenergyerror)


                    elif 'Free Energy via Backward FEP' in line:
                        linesplit=line.split()
                        freeenergyvalue=float(linesplit[5])
                        index=line.find('+/-')
                        newstring=line[index+3:]
                        newstringlinesplit=newstring.split()
                        freenergyerror=float(newstringlinesplit[0])
                        freeenergybwd+=freeenergyvalue
                        
                        freeenergybwdlist.append(freeenergyvalue)
                        freeenergybwderrorlist.append(freenergyerror)

                totalfreeenergylist.append(freeenergy)
                totalfreeenergyviabootstraplist.append(freeenergyviabootstrap)
                totalenthalpylist.append(enthalpy)
                totalentropylist.append(entropy)
                totalfreeenergyfwdlist.append(freeenergyfwd)
                totalfreeenergybwdlist.append(freeenergybwd)


                temp.close()
                

            freeenergylistoflist.append(freeenergylist)
            freeenergyerrorlistoflist.append(freeenergyerrorlist)
            freeenergyviabootstraplistoflist.append(freeenergyviabootstraplist)
            freeenergyviabootstraperrorlistoflist.append(freeenergyviabootstraperrorlist)
            enthalpylistoflist.append(enthalpylist)
            enthalpyerrorlistoflist.append(enthalpyerrorlist)
            entropylistoflist.append(entropylist)
            entropyerrorlistoflist.append(entropyerrorlist)
            freeenergyfwdlistoflist.append(freeenergyfwdlist)
            freeenergyfwderrorlistoflist.append(freeenergyfwderrorlist)
            freeenergybwdlistoflist.append(freeenergybwdlist)
            freeenergybwderrorlistoflist.append(freeenergybwderrorlist)
            freeenergyviabariterlistoflist.append(freeenergyviabariterlist)
            freeenergyviabaritererrorlistoflist.append(freeenergyviabaritererrorlist)


        
        poltype.freeenergylist[j]=freeenergylistoflist
        poltype.freeenergy[j]=sum(FlattenListOfList(poltype,freeenergylistoflist))
        poltype.freeenergyerrorlist[j]=freeenergyerrorlistoflist
        poltype.freeenergylistviabariter[j]=freeenergyviabariterlistoflist
        poltype.freeenergyerrorlistviabariter[j]=freeenergyviabaritererrorlistoflist
        poltype.freeenergylistviabootstrap[j]=freeenergyviabootstraplistoflist
        poltype.freeenergyviabootstrap[j]=sum(FlattenListOfList(poltype,freeenergyviabootstraplistoflist))
        poltype.freeenergyerrorlistviabootstrap[j]=freeenergyviabootstraperrorlistoflist
        poltype.enthalpylist[j]=enthalpylistoflist
        poltype.enthalpy[j]=sum(FlattenListOfList(poltype,enthalpylistoflist))
        poltype.enthalpyerrorlisttotal[j]=enthalpyerrorlistoflist
        poltype.entropylist[j]=entropylistoflist
        poltype.entropy[j]=sum(FlattenListOfList(poltype,entropylistoflist))
        poltype.entropyerrorlisttotal[j]=entropyerrorlistoflist
        poltype.freeenergylistfwd[j]=freeenergyfwdlistoflist
        poltype.freeenergyfwd[j]=sum(FlattenListOfList(poltype,freeenergyfwdlistoflist))
        poltype.freeenergyerrorlistfwd[j]=freeenergyfwderrorlistoflist
        poltype.freeenergylistbwd[j]=freeenergybwdlistoflist
        poltype.freeenergybwd[j]=sum(FlattenListOfList(poltype,freeenergybwdlistoflist))
        poltype.freeenergyerrorlistbwd[j]=freeenergybwderrorlistoflist
        poltype.enthalpyerror[j]=[]
        poltype.entropyerror[j]=[]
        poltype.freeenergyerror[j]=[]
        for m in range(len(poltype.enthalpyerrorlisttotal[j])):
            enthalpyerrorlist=poltype.enthalpyerrorlisttotal[j][m]
            entropyerrorlist=poltype.entropyerrorlisttotal[j][m]
            freeenergyerrorlist=poltype.freeenergyerrorlist[j][m]
            enthalpyerror=np.sqrt(np.sum(np.square(enthalpyerrorlist)))
            entropyerror=np.sqrt(np.sum(np.square(entropyerrorlist)))
            freeenergyerror=np.sqrt(np.sum(np.square(freeenergyerrorlist)))
            poltype.enthalpyerror[j].append(enthalpyerror)
            poltype.entropyerror[j].append(entropyerror)
            poltype.freeenergyerror[j].append(freeenergyerror)
        poltype.enthalpyerror[j]=sum(map(lambda i : i * i, poltype.enthalpyerror[j]))

        poltype.entropyerror[j]=sum(map(lambda i : i * i, poltype.entropyerror[j]))

        poltype.freeenergyerror[j]=sum(map(lambda i : i * i, poltype.freeenergyerror[j]))



        if poltype.solvation==True and poltype.complexation==False:
            poltype.tabledict[j][u'ΔGˢᵒˡᵛ']=poltype.freeenergy[j]
            poltype.tabledict[j][u'ΔGˢᵒˡᵛᵉʳʳ']=poltype.freeenergyerror[j]
            poltype.tabledict[j][u'ΔHˢᵒˡᵛ']=poltype.enthalpy[j]
            poltype.tabledict[j][u'ΔHˢᵒˡᵛᵉʳʳ']=poltype.enthalpyerror[j]
            poltype.tabledict[j][u'ΔSˢᵒˡᵛ']=poltype.entropy[j]
            poltype.tabledict[j][u'ΔSˢᵒˡᵛᵉʳʳ']=poltype.entropyerror[j]
        elif poltype.complexation==True and poltype.solvation==False:
            poltype.tabledict[j][u'ΔGᶜᵒᵐᵖᵘⁿᶜᵒʳʳ']=poltype.freeenergy[j]
            poltype.tabledict[j][u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ']=poltype.freeenergyerror[j]
            poltype.tabledict[j][u'ΔHᶜᵒᵐᵖ']=poltype.enthalpy[j]
            poltype.tabledict[j][u'ΔHᶜᵒᵐᵖᵉʳʳ']=poltype.enthalpyerror[j]
            poltype.tabledict[j][u'ΔSᶜᵒᵐᵖ']=poltype.entropy[j]
            poltype.tabledict[j][u'ΔSᶜᵒᵐᵖᵉʳʳ']=poltype.entropyerror[j]
        elif poltype.complexation==True and poltype.solvation==True:
            if j==0:
                poltype.tabledict[j][u'ΔGᶜᵒᵐᵖᵘⁿᶜᵒʳʳ']=poltype.freeenergy[j]
                poltype.tabledict[j][u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ']=poltype.freeenergyerror[j]
                poltype.tabledict[j][u'ΔHᶜᵒᵐᵖ']=poltype.enthalpy[j]
                poltype.tabledict[j][u'ΔHᶜᵒᵐᵖᵉʳʳ']=poltype.enthalpyerror[j]
                poltype.tabledict[j][u'ΔSᶜᵒᵐᵖ']=poltype.entropy[j]
                poltype.tabledict[j][u'ΔSᶜᵒᵐᵖᵉʳʳ']=poltype.entropyerror[j]
            elif j==1:
                poltype.tabledict[j][u'ΔGˢᵒˡᵛ']=poltype.freeenergy[j]
                poltype.tabledict[j][u'ΔGˢᵒˡᵛᵉʳʳ']=poltype.freeenergyerror[j]
                poltype.tabledict[j][u'ΔHˢᵒˡᵛ']=poltype.enthalpy[j]
                poltype.tabledict[j][u'ΔHˢᵒˡᵛᵉʳʳ']=poltype.enthalpyerror[j]
                poltype.tabledict[j][u'ΔSˢᵒˡᵛ']=poltype.entropy[j]
                poltype.tabledict[j][u'ΔSˢᵒˡᵛᵉʳʳ']=poltype.entropyerror[j]

    tables.WriteTableUpdateToLog(poltype,verbose=False)
    tempname='BARResults.csv'
    with open(poltype.outputpath+tempname, mode='w') as energy_file:
        energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        table.append(['Vdw-Lambda']+FlattenListOfListGenerateBARPairs(poltype,poltype.vdwlambdascheme))
        table.append(['Ele-Lambda']+FlattenListOfListGenerateBARPairs(poltype,poltype.estatlambdascheme))
        table.append(['Rest-Lambda']+FlattenListOfListGenerateBARPairs(poltype,poltype.restlambdascheme))
        if poltype.solvation==True and poltype.complexation==False:
            table.append([u'ΔGˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.freeenergylist[0]))
            table.append([u'ΔGˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlist[0]))
            table.append([u'ΔGˢᵒˡᵛᵇᵃʳᶦᵗᵉʳ']+FlattenListOfList(poltype,poltype.freeenergylistviabariter[0]))
            table.append([u'ΔGˢᵒˡᵛᵉʳʳᵇᵃʳᶦᵗᵉʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistviabariter[0]))
            table.append([u'ΔGˢᵒˡᵛᵇᵒᵒᵗˢᵗʳᵃᵖ']+FlattenListOfList(poltype,poltype.freeenergylistviabootstrap[0]))
            table.append([u'ΔGˢᵒˡᵛᵉʳʳᵇᵒᵒᵗˢᵗʳᵃᵖ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistviabootstrap[0]))
            table.append([u'ΔGˢᵒˡᵛᶠʷᵈ']+FlattenListOfList(poltype,poltype.freeenergylistfwd[0]))
            table.append([u'ΔGˢᵒˡᵛᵉʳʳᶠʷᵈ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistfwd[0]))
            table.append([u'ΔGˢᵒˡᵛᵇʷᵈ']+FlattenListOfList(poltype,poltype.freeenergylistbwd[0]))
            table.append([u'ΔGˢᵒˡᵛᵉʳʳᵇʷᵈ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistbwd[0]))
            table.append([u'ΔHˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.enthalpylist[0]))
            table.append([u'ΔHˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.enthalpyerrorlisttotal[0]))
            table.append([u'ΔSˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.entropylist[0]))
            table.append([u'ΔSˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.entropyerrorlisttotal[0]))
        elif poltype.complexation==True and poltype.solvation==False:
            table.append([u'ΔGᶜᵒᵐᵖ']+FlattenListOfList(poltype,poltype.freeenergylist[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵉʳʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlist[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵇᵃʳᶦᵗᵉʳ']+FlattenListOfList(poltype,poltype.freeenergylistviabariter[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵉʳʳᵇᵃʳᶦᵗᵉʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistviabariter[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵇᵒᵒᵗˢᵗʳᵃᵖ']+FlattenListOfList(poltype,poltype.freeenergylistviabootstrap[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵉʳʳᵇᵒᵒᵗˢᵗʳᵃᵖ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistviabootstrap[0]))
            table.append([u'ΔGᶜᵒᵐᵖᶠʷᵈ']+FlattenListOfList(poltype,poltype.freeenergylistfwd[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵉʳʳᶠʷᵈ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistfwd[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵇʷᵈ']+FlattenListOfList(poltype,poltype.freeenergylistbwd[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵉʳʳᵇʷᵈ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistbwd[0]))
            table.append([u'ΔHᶜᵒᵐᵖ']+FlattenListOfList(poltype,poltype.enthalpylist[0]))
            table.append([u'ΔHᶜᵒᵐᵖᵉʳʳ']+FlattenListOfList(poltype,poltype.enthalpyerrorlisttotal[0]))
            table.append([u'ΔSᶜᵒᵐᵖ']+FlattenListOfList(poltype,poltype.entropylist[0]))
            table.append([u'ΔSᶜᵒᵐᵖᵉʳʳ']+FlattenListOfList(poltype,poltype.entropyerrorlisttotal[0]))
        elif poltype.complexation==True and poltype.solvation==True:
            table.append([u'ΔGˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.freeenergylist[1]))
            table.append([u'ΔGˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlist[1]))
            table.append([u'ΔGˢᵒˡᵛᵇᵃʳᶦᵗᵉʳ']+FlattenListOfList(poltype,poltype.freeenergylistviabariter[1]))
            table.append([u'ΔGˢᵒˡᵛᵉʳʳᵇᵃʳᶦᵗᵉʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistviabariter[1]))
            table.append([u'ΔGˢᵒˡᵛᵇᵒᵒᵗˢᵗʳᵃᵖ']+FlattenListOfList(poltype,poltype.freeenergylistviabootstrap[1]))
            table.append([u'ΔGˢᵒˡᵛᵉʳʳᵇᵒᵒᵗˢᵗʳᵃᵖ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistviabootstrap[1]))
            table.append([u'ΔGˢᵒˡᵛᶠʷᵈ']+FlattenListOfList(poltype,poltype.freeenergylistfwd[1]))
            table.append([u'ΔGˢᵒˡᵛᵉʳʳᶠʷᵈ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistfwd[1]))
            table.append([u'ΔGˢᵒˡᵛᵇʷᵈ']+FlattenListOfList(poltype,poltype.freeenergylistbwd[1]))
            table.append([u'ΔGˢᵒˡᵛᵉʳʳᵇʷᵈ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistbwd[1]))
            table.append([u'ΔHˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.enthalpylist[1]))
            table.append([u'ΔHˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.enthalpyerrorlisttotal[1]))
            table.append([u'ΔSˢᵒˡᵛ']+FlattenListOfList(poltype,poltype.entropylist[1]))
            table.append([u'ΔSˢᵒˡᵛᵉʳʳ']+FlattenListOfList(poltype,poltype.entropyerrorlisttotal[1]))
            table.append([u'ΔGᶜᵒᵐᵖ']+FlattenListOfList(poltype,poltype.freeenergylist[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵉʳʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlist[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵇᵃʳᶦᵗᵉʳ']+FlattenListOfList(poltype,poltype.freeenergylistviabariter[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵉʳʳᵇᵃʳᶦᵗᵉʳ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistviabariter[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵇᵒᵒᵗˢᵗʳᵃᵖ']+FlattenListOfList(poltype,poltype.freeenergylistviabootstrap[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵉʳʳᵇᵒᵒᵗˢᵗʳᵃᵖ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistviabootstrap[0]))
            table.append([u'ΔGᶜᵒᵐᵖᶠʷᵈ']+FlattenListOfList(poltype,poltype.freeenergylistfwd[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵉʳʳᶠʷᵈ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistfwd[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵇʷᵈ']+FlattenListOfList(poltype,poltype.freeenergylistbwd[0]))
            table.append([u'ΔGᶜᵒᵐᵖᵉʳʳᵇʷᵈ']+FlattenListOfList(poltype,poltype.freeenergyerrorlistbwd[0]))
            table.append([u'ΔHᶜᵒᵐᵖ']+FlattenListOfList(poltype,poltype.enthalpylist[0]))
            table.append([u'ΔHᶜᵒᵐᵖᵉʳʳ']+FlattenListOfList(poltype,poltype.enthalpyerrorlisttotal[0]))
            table.append([u'ΔSᶜᵒᵐᵖ']+FlattenListOfList(poltype,poltype.entropylist[0]))
            table.append([u'ΔSᶜᵒᵐᵖᵉʳʳ']+FlattenListOfList(poltype,poltype.entropyerrorlisttotal[0]))

        newtable=list(map(list, zip(*table)))
        for row in newtable:
            energy_writer.writerow(row)
            
def FlattenListOfList(poltype,listoflist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    flat=[item for sublist in listoflist for item in sublist]
    newflat=[]
    for item in flat:
        try:
            float(item)
            item=round(item,3)
        except:
            pass
        newflat.append(item)
    return newflat

def FlattenListOfListGenerateBARPairs(poltype,listoflist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """   """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newflat=[]
    newlistoflist=[]
    for j in range(len(listoflist)):
        ls=listoflist[j]
        temp=[]
        for i in range(len(ls)-1):
            first=str(ls[i]).replace('Sim','').replace('Vdw','V').replace('Ele','E').replace('Res','R')
            second=str(ls[i+1]).replace('Sim','').replace('Vdw','V').replace('Ele','E').replace('Res','R')
            new=str(second)+'-'+str(first)
            newflat.append(new)
            temp.append(new)
        newlistoflist.append(temp)
    return newflat,newlistoflist


def GenerateBARPairs(poltype,ls):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newflat=[]
    for i in range(len(ls)-1):
        first=ls[i]
        second=ls[i+1]
        new=str(second)+'-'+str(first)
        newflat.append(new)
    return newflat

def DeleteBARFiles(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for j in range(len(poltype.simfoldname)):
       simfoldname=poltype.simfoldname[j]
       os.chdir(poltype.outputpath+simfoldname)
       proddynoutfilepathlistoflist=poltype.proddynoutfilepath[j]
       for k in range(len(proddynoutfilepathlistoflist)):
           proddynoutfilepath=proddynoutfilepathlistoflist[k]
           for i in range(len(proddynoutfilepath)):
               outputfilepath=proddynoutfilepath[i]
               path,tail=os.path.split(outputfilepath)
               os.chdir(path)
               if not os.path.isdir('Backup'):
                   os.mkdir('Backup')
               files=os.listdir()
               for f in files:
                   if '.bar' in f or 'BAR1' in f or 'BAR2' in f:
                       os.replace(f,os.path.join('Backup',f)) 
       os.chdir('..')

   

def BARProtocol(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    if poltype.barfilesfinished==False:
        for i in range(len(poltype.baroutputfilepath)):
            baroutputfilepathlist=poltype.baroutputfilepath[i]
            for j in range(len(baroutputfilepathlist)):
                baroutputfilepath=baroutputfilepathlist[j]
                if term.CheckFilesTermination(poltype,baroutputfilepath)[0]==False or poltype.redobar==True:         
                    ExecuteBAR(poltype)
        messages=[]
        for i in range(len(poltype.baroutputfilepath)):
            baroutputfilepathlist=poltype.baroutputfilepath[i]
            for j in range(len(baroutputfilepathlist)):
                baroutputfilepath=baroutputfilepathlist[j]
                checkfin=term.CheckFilesTermination(poltype,baroutputfilepath)
                finished=checkfin[0]
                percentfinished=checkfin[1]
                while finished==False:
                    msg='BAR is not complete, '+str(percentfinished)+'% of jobs finished out of a total = '+str(len(baroutputfilepath))+ ' jobs '
                    if msg not in messages:
                        poltype.WriteToLog(msg,prin=True)
                        messages.append(msg)
                    time.sleep(poltype.waitingtime)
                    checkfin=term.CheckFilesTermination(poltype,baroutputfilepath)
                    finished=checkfin[0]
                    percentfinished=checkfin[1]

    if poltype.usetinkerforthermoprops==True:
        for i in range(len(poltype.thermooutputfilepath)):
            thermooutputfilepathlist=poltype.thermooutputfilepath[i]
            for j in range(len(thermooutputfilepathlist)):
                thermooutputfilepath=thermooutputfilepathlist[j]
                if term.CheckFilesTermination(poltype,thermooutputfilepath)[0]==False:
                    ExecuteBARSecondOption(poltype,poltype.freeenergyjobsfilename)
        messages=[]
        for i in range(len(poltype.thermooutputfilepath)):
            thermooutputfilepathlist=poltype.thermooutputfilepath[i]
            for j in range(len(thermooutputfilepathlist)):
                thermooutputfilepath=thermooutputfilepathlist[j]
                finished=False
                percentfinished=0
                while finished==False:
                    msg='BAR option 2 is not complete ,'+str(percentfinished)+'% of jobs finished'
                    if msg not in messages:
                        poltype.WriteToLog(msg,prin=True)
                        messages.append(msg)
                    time.sleep(poltype.waitingtime)
                    checkfin=term.CheckFilesTermination(poltype,thermooutputfilepath)
                    finished=checkfin[0]
                    percentfinished=checkfin[1]

    poltype.WriteToLog('Generating .bar files is complete',prin=True)
    if poltype.usetinkerforthermoprops==True:
        SumTheFreeEnergyStepsFromBAR(poltype)
    else:
        poltype.WriteToLog('Computing BAR from .bar files',prin=True)
        ComputeThermoProperties(poltype)
    if poltype.complexation==True: 
        poltype.freeenergy[0]=poltype.freeenergy[0]+poltype.rescorrection
        for idx in range(len(poltype.freeenergyconv[0])):
            poltype.freeenergyconv[0][idx]+=poltype.rescorrection
        poltype.tabledict[0][u'ΔGᵃⁿᵃᶜᵒᵐᵖᶜᵒʳʳ']=poltype.rescorrection
        poltype.tabledict[0][u'ΔGᶜᵒᵐᵖᶜᵒʳʳ']=poltype.freeenergy[0]
        if poltype.solvation==True: 
            poltype.tabledict[1][u'ΔGᵃⁿᵃᶜᵒᵐᵖᶜᵒʳʳ']=poltype.rescorrection
            poltype.tabledict[1][u'ΔGᶜᵒᵐᵖᶜᵒʳʳ']=poltype.freeenergy[0]
    #DeleteBARFiles(poltype) # if want to add more time to MD, then delete BAR files, so when run again, BAR is recomputed from new files
    poltype.WriteToLog('Job is complete your molecule has been annihilated!',prin=True)


