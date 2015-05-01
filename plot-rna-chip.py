#!/usr/local/bin/env python

#python code/plot-rna-chip.py rhabdoid/coverage/TSS_2000_all/
#rhabdoid/genes/fc/ sig_fc.FDR5.rpkm H3K4me3 TUM Norm 5 rhabdoid/alllibsvs.TUM.txt

import sys, os.path
sys.path.append('/Users/ssaberim/epigenomics/code')
import read,write,bedtools, analyse
import matplotlib.pyplot as plt

dirPD=str(sys.argv[1])
dirin=str(sys.argv[2])
fl=str(sys.argv[3])
fl=dirin+fl
mark=str(sys.argv[4])

sample1=str(sys.argv[5])
sample2=str(sys.argv[6])
field=int(sys.argv[7])
libsfile=str(sys.argv[8])
hugo=str(sys.argv[9])



libsdata=read.read_dat(libsfile,'\t')

genelist=[]
genelist=read.read_dat(hugo)

if os.path.isfile(fl) and os.path.getsize(fl) > 0:
	print fl,field,libsfile
	nm,libsRNA,gene,genex,libs=analyse.heat_rna(fl,genelist,field,libsfile)
	rpkmmat=nm*1.
	
	plt.xlabel(sample1+' vs. '+sample2+' rpkm')
	plt.savefig(dirin+sample1+'-'+sample2+'-rpkm.pdf', bbox_inches='tight')
	print dirin+sample1+'-'+sample2+'-rpkm.pdf'
	plt.close()
#	try:
	
	#print genex
	#genex=genex[::-1]
	a=analyse.heatmap_vec(libsRNA,gene,nm)
	plt.savefig(dirin+sample1+'-'+sample2+'-rpkm-dend.pdf', bbox_inches='tight')
	write.write_mat(libsRNA,genex,nm,dirin+sample1+'-'+sample2+'-rpkm-dend.txt')


	indmark=libsdata[0].index(mark)
	tmp=[]
	tmp2=[]
	libsdata=read.read_dat(libsfile,'\t')
	for i in libs:
		for j in libsdata:
			if i.split('-')[0] in j and len(j[indmark])>0 :
				tmp.append(j[indmark]+'-'+j[field])
				tmp2.append(j[0]+'-'+j[field])
				break
	if len(gene)>0:
		#try:
		nm=analyse.heat_density(dirPD,mark,tmp,gene,genelist)
		#except:
		#  print 'pd aha'

		#print mark, 'correlation=%3.4f std=%3.4f' % analyse.rpkm_pd_corr(rpkmmat,nm,libsRNA,tmp2)
		plt.xticks(range(len(tmp2)),tmp2)
		plt.xlabel(mark+' Promoter Density')
		plt.savefig(dirin+mark+sample1+'-'+sample2+'-PD.pdf', bbox_inches='tight')
		try:
		  a=analyse.heatmap_vec(tmp2,genex,nm)
		  plt.savefig(dirin+mark+sample1+'-'+sample2+'-dend.pdf', bbox_inches='tight')
		  write.write_mat(tmp2,genex,nm,dirin+mark+sample1+'-'+sample2+'-dend.txt')
		except:
		  pass
	else:
		print dirin+mark+sample1+'-'+sample2+'-PD.pdf'
	plt.close()
	print 'number of rpkm plotted genes are:', len(genex)

