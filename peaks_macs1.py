import os, sys 
import commands
import read
import write
import numpy as np

tmp=sys.argv
fl=tmp[1]
table=read.read_dat("peaks.txt",'\t')
data=read.read_dat(fl,'\t')
mark=map(str,tmp[2].split('_'))
#mark=mark[mark.index('_')+1:]
mark=mark[-1]
vecs=['pcDNA','K4E','Y69H','D83V']
#print data[1]D83V
for row in table:
	for tmp in row:
		if "bwa-0." in tmp:
			if tmp in data[1][5]:
				#print data[1][5]
				
				cellline = data[-1][-1]
				cellline = cellline [:cellline.index('_')]
				print cellline

vs='WT'

out='table'+'-'+vs+'_all_vectors.txt'
f=open(out,'a')
tmp=fl[fl.index('/')+1:]
rep=int(tmp[-len('2_peaks.xls'):-len('2_peaks.xls')+1])
vec=tmp.split('_')[1]
tmp = (data[-1][-1])

#print tmp


tmp=tmp.split('_')
coverage=int(tmp[-1])
mark=str(tmp[2])
solution=str(tmp[1])
print >> f, mark,rep,solution,vs,cellline,coverage
#print  mark,rep,solution,vs,cellline,coverage


f.close()

