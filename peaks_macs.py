import os, sys 
import commands
import read
import write
import numpy as np

tmp=sys.argv
fl=tmp[1]
#vs=str(tmp[3])
mark=map(str,tmp[2].split('-'))
mark=mark[1]
vs='input'
data=read.read_dat(fl)
out='table'+'-'+vs+'.txt'
f=open(out,'a')
tmp=fl[fl.index('/')+1:]
solution=tmp.split('_')[1]
rep=tmp.split('_')[2]

tmp = str(data[-1][-1])
tmp=tmp.split('_')
rank=tmp[2]
coverage=tmp[-1]
vec=tmp[0]
try:
	print >> f, rep,mark,vec,solution,vs,coverage
except:
	print vec
f.close()

