#!/usr/local/bin/env python
#coding: utf8 
import sys,read, write,bedtools, analyse


mark=sys.argv[1]
sample1=sys.argv[2]
sample2=sys.argv[3]
num=int(sys.argv[4])
T=float(sys.argv[5])
dirout=sys.argv[6]
field=int(sys.argv[7])
libsfile=str(sys.argv[8])
TSSdir=str(sys.argv[9])


beds={}
libsdata=read.read_dat(libsfile,'\t')
indmark=libsdata[0].index(mark)
list1=[]
list2=[]
for i in xrange(1,len(libsdata)):
	if sample1 in libsdata[i][field] and len(libsdata[i][indmark])>1:
		list1.append(libsdata[i][indmark])

	elif sample2 in libsdata[i][field] and len(libsdata[i][indmark])>1:
		list2.append(libsdata[i][indmark])


beds={}
files1,beds[sample1]=read.readall_bed(TSSdir,mark,list1)
files2,beds[sample2]=read.readall_bed(TSSdir,mark,list2)
genes=read.read_gene_pos('rhabdoid/coverage/TSS_2000_all/hg19v69_genes.TSS_2000.pc.A03480.H3K27me3.GE02.coverage')

compared=bedtools.Compare_all_bed(beds[sample1],beds[sample2])
ranked=[]
ind=-1
for j in compared:
   for i in j:
     tmp=bedtools.filter_bed_num(i,num)
     tmp2=bedtools.intersect_bed_gene(tmp,genes)
     ranked.append(tmp2)
     ind+=1


allgenes=bedtools.intersect_allgenes(ranked)

thresh=int(T*float(ind))
commongenes=bedtools.filter_common_genes(allgenes,thresh)
print 'Num# of common genes for', mark+',', sample1, 'vs.',sample2, 'are:', len(commongenes.keys())
outfile=dirout+mark+'genes'+sample1+'-'+sample2+'.txt'
f=open(outfile,'w')


for i in commongenes.keys():
    print >> f, i
    
f.close()







