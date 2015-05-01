#!python2.7.5

import time, os, sys
dr=os.getcwd()
sys.path.insert(0, '/home/ssaberi/code')
import read,write,bedtools

dir1="/projects/epigenomics/nci_rt/ChIPseq_RT/macs2/TUM_NS_H3K27me3-05/"
dir2="/projects/epigenomics/nci_rt/ChIPseq_RT/macs2/TUM_NS_H3K27me3-05-mine/"
direct1=os.popen('ls '+dir1 + '*narr*')
direct2=os.popen('ls '+dir2 + '*')

f=open('inetrsect.out','w')
ind=0
files1=[]
for tmp in direct1.readlines():
 files1.append(tmp)
files2=[]
for tmp in direct2.readlines():
 files2.append(tmp)



for fl1 in files1:
	tmp=fl1.split('/')[-1]
	libs=tmp.split('_')[:2]
	#print libs
	ind+=1	
	for fl2 in files2:
		fl2=str(fl2)
		if libs[1] in fl2 and  libs[0] in fl2:
			bed1=read.read_bed(fl1[:-1])
			print fl2, fl1
			bed2=read.read_bed(fl2[:-1])
			inter=bedtools.intersect_beds(bed2,bed1,1000)
			a=bedtools.print_bed(inter)
			print >> f, libs[0], libs[1], bedtools.print_bed(bed1), bedtools.print_bed(bed2), a, a/float(bedtools.print_bed(bed1))

print ind
f.close()
