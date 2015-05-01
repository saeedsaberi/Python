#!python2.7.1

import matplotlib.pyplot as plt
import sys, os, bx
sys.path.append('/home/ssaberi/code')
import read,write,bedtools, analyse
import numpy as np
fl=sys.argv[1]
dist=int(sys.argv[2])
from bx.bbi.bigwig_file import BigWigFile

genes=read.dat("/home/ssaberi/resources/list.genes.txt",'\t')
table=read.dat("/projects/epigenomics/MarcoJuliaPon/peaks.txt",'\t')
mygenes=read.dat("/projects/epigenomics/MarcoJuliaPon/mygenes.txt",'\t')
ens=[]
for i in mygenes:
	for gn in genes:
		if i in gn[0]:
			ens.append(gn[1])
			break

genespos=read.read_gene_pos('/home/ssaberi/resources/hg19v69_genes.TSS_2000.pc.A03480.H3K27me3.GE02.coverage')
genesbed=bedtools.makebed_genpos(ens,genespos,100000)
              


f = open(fl)
bw = BigWigFile(file=f)
mat=[]
for bed_i in genesbed:
   vals = bw.get( bed_i[0], bed_i[1], bed_i[2])
   mat.append(np.array(vals))
mat=np.array(mat)
plt.matshow(mat,aspect='auto',cmap='YlOrBr')
fl=fl[-fl[::-1].index('/'):-fl[::-1].index('.')]
plt.save(fl+".pdf")
