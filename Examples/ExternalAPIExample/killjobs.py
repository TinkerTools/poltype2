import os
import sys
import subprocess
cpunodepath='cpunodes.txt'

def ReadCPUNodes(cpunodepath):
    cpunodes=[]
    temp=open(cpunodepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if '#' not in line:
            linesplit=line.split()
            cpunodes.append(linesplit[0])
    return cpunodes


def KillSignal(cpunodes,killstring):
    for node in cpunodes:
        cmdstr='ssh %s "%s"'%(node,killstring)
        p = subprocess.Popen(cmdstr, stdout=subprocess.PIPE,shell=True)


killstring='pkill -U bdw2292'
cpunodes=ReadCPUNodes(cpunodepath)
KillSignal(cpunodes,killstring)
