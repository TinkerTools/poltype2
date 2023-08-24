import time
import os,sys
import numpy as np
import subprocess
from multiprocessing import Pool


def single_point(params, fname, savetxt=True):
  node = []
  keys = [line.split()[0] for line in open("p0.txt").readlines()]
  subprocess.run("rm -rf ./outs/*.out", shell=True)
  oFile = open("%s_out.key" %(fname), 'w')
  lines = open("%s.key").readlines()
  for line in lines:
    if "PRM" in line:
      for i in range(len(keys)):	
        search = "PRM"+ str(i) + '_'
        if search in line and (line[0:1] != "#"):
          params[i] = "%20.10f"%(float(params[i]))
          line=line.replace(search, str(params[i]))
    oFile.write(line)
  oFile.close()

  freq_QM = []
  freq_MM = []
  #Check whether all analyze jobs finished!
  files = [line.split(".xyz")[0] for line in open(filelist).readlines()]
  nFile = len(files)
  polarizeexe = "/home/xy3866/Tinker/bin/vibrate" 
  cmdstr_1 ="%s %s.txyz %s_new.key CR > %s.out &\n" %(polarizeexe, fname, fname, fname)
  subprocess.run(cmdstr_1, shell=True)
  os.system("grep 'result: ' > result_1.p")
  
  freq_QM = np.loadtxt('frequency.log',usecols=(-1,))
  freq_MM = np.loadtxt('result_1.p',usecols=(-1,))
  if savetxt:
    np.savetxt("param.dat", np.transpose(params), fmt="%20.10f")
  QM = np.array(freq_QM)
  MM = np.array(freq_MM)
  return QM,MM

def costFUNC(params, x0):
  QM, MM = single_point(params, 'filelist_0.txt', "nodelist") 
  QM_1 = []
  MM_1 = []
  for i, j in enumerate(QM):
    if(j>1000):
      QM_1.append(j)
      MM_1.append(MM[i])
      
  
  cost = MSE(QM_1,MM_1)
  print("MSE: %15.8f "%(np.sqrt(cost)))
  print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

  return cost

def fitting(fname):
  x0 = np.loadtxt("p0.txt",usecols=(-1,))
  ret = least_squares(costFUNC, x0, jac='3-point', max_nfev=100)
  np.savetxt("p1.txt", ret.x,fmt='%13.5f')
  subprocess.run("paste p0.txt p1.txt >temp && mv temp p0.txt", shell=True)
