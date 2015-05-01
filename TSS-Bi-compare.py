import sys,read, write,bedtools, analyse

T=0.8
mark=sys.argv[1]
sample1=sys.argv[2]
sample2=sys.argv[3]
dirout=sys.argv[6]
num=int(sys.argv[4])
T=float(sys.argv[5])
if mark=='H3K27me3':
	tum=["A33463","A33882","A33889","A33896","A34008","A34031","A34070","A34190","A34457","A34464"]
	NS=['A03480','A03284']
	cell=["A32397","A32404","A32411","A32418"]
	ES=["A32425","A32432","A33456"]
if mark=='H3K4me3':
	tum=["A33880","A34029","A33887","A33894","A34006","A34068","A33461","A34188","A34455","A34462"]
	NS=['A03478','A03282']
	ES=["A32423","A32430","A33454"]

	cell=["A32395","A32402","A32409","A32416"]

beds={}
filesTUM,bedsTUM=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,tum)
beds['TUM']=bedsTUM
filesES,bedsES=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,ES)
beds['ES']=bedsES
filesCELL,bedsCELL=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,cell)
beds['CELL']=bedsCELL
filesNS,bedsNS=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,NS)
beds['NS']=bedsNS
beds['TUMCELL']=beds['TUM']+beds['CELL']
genes=read.read_gene_pos('rhabdoid/coverage/TSS_2000_all/hg19v69_genes.TSS_2000.pc.A03480.H3K27me3.GE02.coverage')

#sample e.i. 'NS' or 'TUM'

compared=bedtools.Compare_all_bed(beds[sample1])
ranked=[]
ind=-1
for j in compared:
   for i in j:
     tmp=bedtools.filter_bed_num(i,num)
     bedtools.print_bed(tmp)
     tmp2=bedtools.intersect_bed_gene(tmp,genes)
     ranked.append(tmp2)
     ind+=1
     #print ind
   
allgenes=bedtools.intersect_allgenes(ranked)
thresh=int(T*float(ind))
commongenes=bedtools.filter_common_genes(allgenes,thresh)
ES2=["A32423","A32430","A33454"]



#Bivalent ones
#H3K27me3 libs:
tum1=["A33463","A33882","A33889","A33896","A34008","A34031","A34070","A34190","A34457","A34464"]
NS1=['A03480','A03284']
cell1=["A32397","A32404","A32411","A32418"]
ES1=["A32425","A32432","A33456"]
#H3K4me3 libs:
tum2=["A33880","A34029","A33887","A33894","A34006","A34068","A33461","A34188","A34455","A34462"]
NS2=['A03478','A03282']
ES2=["A32423","A32430","A33454"]
cell2=["A32395","A32402","A32409","A32416"]
mark1=str(mark)
mark='H3K27me3'
beds1={}
filesTUM,bedsTUM1=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,tum1)
beds1['TUM']=bedsTUM1
filesES,bedsES1=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,ES1)
beds1['ES']=bedsES1
filesCELL,bedsCELL1=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,cell1)
beds1['CELL']=bedsCELL1
filesNS,bedsNS1=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,NS1)
beds1['NS']=bedsNS1
beds1['TUMCELL']=beds1['TUM']+beds1['CELL']

mark='H3K4me3'
beds2={}
filesTUM,bedsTUM2=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,tum2)
beds2['TUM']=bedsTUM2
filesES,bedsES2=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,ES2)
beds2['ES']=bedsES2
filesCELL,bedsCELL2=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,cell2)
beds2['CELL']=bedsCELL2
filesNS,bedsNS2=read.readall_bed('rhabdoid/coverage/TSS_2000_all/',mark,NS2)
beds2['NS']=bedsNS2
beds2['TUMCELL']=beds2['TUM']+beds2['CELL']
print 'Starting bivalent'
prod2=bedtools.bedsall_prod(beds1[sample2],beds2[sample2])
prod1=bedtools.bedsall_prod(beds1[sample1],beds2[sample1])


compared_bi=bedtools.Compare_all_bed(prod1,prod2) # all biv in sample 2 not existance in sample1
ranked_bi=[]
ind=-1
for j in compared_bi:
   for i in j:
     tmp=bedtools.filter_bed_num(i,num)
     tmp2=bedtools.intersect_bed_gene(tmp,genes)
     ranked_bi.append(tmp2)
     ind+=1

   
allgenes_bi=bedtools.intersect_allgenes(ranked_bi)
thresh=int(T*float(ind))

commongenes_bi=bedtools.filter_common_genes(allgenes_bi,thresh)
print 'Num# of common genes for,', mark1, 'vs. BiV are:', len(commongenes_bi.keys())
print 'Num# of common genes for,', sample1, 'are:', len(commongenes.keys())


commonall=bedtools.intersect_gene_names(commongenes,commongenes_bi)

print 'Num# of common genes for histon modification', mark1,'ending with bivalent state for,', sample1, 'vs. ', sample2, 'are:', len(commonall.keys())
f=open(dirout+mark1+'genes'+sample1+'-'+sample2+'.txt','w')
for i in commonall.keys():
    print >> f, i
f.close()







