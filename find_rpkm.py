#!/usr/local/bin/env python
#coding: utf8 
import sys, os, numpy as np
sys.path.append('/Users/ssaberim/epigenomics/code')
import read,write,bedtools,analyse



ID=str(sys.argv[1])
drrna=str(sys.argv[2])
sample=str(sys.argv[3])
outfile=str(sys.argv[4])
field=int(sys.argv[5])

libsfile=str(sys.argv[6])

mark='RNA'
libsdata=read.read_dat(libsfile,'\t')

indmark=libsdata[0].index(mark)

list1=[]
for i in range(1,len(libsdata)):

	if libsdata[i][field] in sample and len(libsdata[i][field])>0:
		list1.append(libsdata[i][indmark])

hugo="/Users/ssaberim/epigenomics/resources/list.genes2.txt"
hugo=read.read_dat(hugo)
#print hugo[0]
ID=read.read_dat(ID,'\t')
for i in range(len(ID)):
	for j  in range(len(hugo)):
		#print ID[i][0]
		if ID[i][0] in hugo[j][0]+hugo[j][1]:
			ID[i][0]=hugo[j][1]
			break

#print ID
f=open(outfile,'a')
#print outfile
for lib in list1:
	st='ls '+drrna+lib+'*A.rpkm.pc'
	direct=os.popen(st)
	fl=direct.readline()[:-1]
	#print fl
	fl=read.read_dat(fl,'\t')
	
	ind=0
	nfound=[]
	for i in ID:
		found=False
		#print i[0]
		for j in fl:
			if i[0] in j[0]:
				found= True
				a=float(j[2])
				print >> f, lib, sample, i[0],  a
				#print i[0]
				break

f.close()
