#!python2.7.1

import sys,read,write,bedtools
import time
#import psyco
#psyco.full()


start = time.time()
wigfile1=str(sys.argv[1])
wigfile2=str(sys.argv[2])
res=int(sys.argv[3])
peaknum=int(sys.argv[4])
outfile=str(sys.argv[5])
region=1000
gap=100
# data is imported


print 'runing with: wigfile1=', wigfile1, 'runing with: wigfile2=', wigfile2,
print ', res=',res,', outfile=',outfile,', peaknum=',peaknum
wig1=read.read_wig(wigfile1)
wig2=read.read_wig(wigfile2)
print 'reading is over, calculating domains'
end = time.time()
print 'reading time is:',end - start
wig1=bedtools.filter_bed_chrom(wig1)
wig2=bedtools.filter_bed_chrom(wig2)

binned1=bedtools.bin_wig(wig1,res)
binned2=bedtools.bin_wig(wig2,res)

compared=bedtools.compare_bed(binned1,binned2)
#a,b=histogram(compared['chr1'][:,2])
#plt.plot(b,a[:len(a)])

filtered=bedtools.filter_bed_num(compared,peaknum)
write.write_bed(filtered,outfile)

print ' calculating FDR and domains'
#domain,FDR=bedtools.FDR_bed(filtered,10,gap,region)

#write.write_bed(domain,outfile)
end = time.time()
print 'running time is:',end - start
