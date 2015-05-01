#!/usr/local/bin/env python
#coding: utf8 
import os,sys,read, write,bedtools, analyse


mark=sys.argv[1]
sample1=sys.argv[2]
sample2=sys.argv[3]
pival=float(sys.argv[4])

dirout=sys.argv[5]
field=int(sys.argv[6])
libsfile=str(sys.argv[7])
TSSdir=str(sys.argv[8])
flagstat='/Users/ssaberim/epigenomics/rhabdoid/bamstat.Aligned.17'

beds={}
libsdata=read.read_dat(libsfile,'\t')
indmark=libsdata[0].index(mark)
list1=[]
list2=[]
for i in xrange(1,len(libsdata)):
	if sample1 in libsdata[i][field] and len(libsdata[i][indmark])>1 and '#' not in libsdata[i][indmark]:
		list1.append(libsdata[i][indmark])

	elif sample2 in libsdata[i][field] and len(libsdata[i][indmark])>1 and '#' not in libsdata[i][indmark]:
		list2.append(libsdata[i][indmark])


beds={}
#print list1
#print indmark
files1,beds[sample1]=read.readall_bed(TSSdir,mark,list1,flagstat)
files2,beds[sample2]=read.readall_bed(TSSdir,mark,list2,flagstat)
direct=os.popen('ls '+TSSdir+'/*.coverage')
genefile=direct.readlines()[0].strip()
genes=read.read_gene_pos(genefile)

compared=bedtools.Compare_all_bed_pi(beds[sample1],beds[sample2],pival)
commongenes=bedtools.intersect_bed_gene(compared,genes)

print 'Num# of common genes for', mark+',', sample1, 'vs.',sample2, 'are:', len(commongenes.keys())
outfile=dirout+mark+'genes'+sample1+'-'+sample2+'.txt'
f=open(outfile,'w')


for i in commongenes.keys():
    print >> f, i
    
f.close()







