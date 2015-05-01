import os, sys 
import commands
import read
import numpy as np
import time
import subprocess
import write
table=read.read_dat("peaks.txt",'\t')
data=table
cells=["D83V","K4E", "pcDNA","WT","Y69H"]
marks=["H3K27ac","H3K4me3","H3K9me3","V5","input"]
ran=0 
for irow in xrange(2,len(data)):
	row=data[irow]
	vector=row[-3]
	solution=row[-1][:-2]
	mymark=""
	#print row[-2]
	if "input" in row:
				mymark="input"
	else:
		mymark=row[-2]
	if '1' in row[-4]:
		rep='1'
	else:
		rep='2'
	bam=row[5]
	if 'WT' not in vector and 'input' not in mymark and mymark in marks:
		for ind in xrange(-1,4):
			try:
				ii=data[irow+ind][-3]
				#print ii
				if 'WT' in ii:
					wtfile=data[irow+ind][5]
					#print data[irow+ind][5]
					p=ii
			except:
				print 'end '

		bam=os.popen("ls -f " + bam +"*.bam")
		bam=bam.readline()[:-1]
		wt=os.popen("ls -f " + wtfile +"*.bam")
		wt=wt.readline()[:-1]
		if not os.path.isdir('WT-'+mymark):
			os.makedirs('WT-'+mymark)
		os.chdir('WT-'+mymark)
		#print os.getcwd()
		if  "H3K9me3" not in mymark:
			string=" nohup macs2 callpeak -t "+wt+" -c " + bam
			string+=" -f BAM -g hs -n "
			tmp= str(vector)+'_'+(solution)+'_'+(mymark)+'_'+(rep)
			string+= tmp + " -B -q 0.01  --SPMR --nomodel " 
			#string+= "> "  +tmp+'.log ' #" 2>&1&"
		elif "H3K9" in mymark:
			string=" nohup macs2 callpeak -t "+wt +" -c " + bam
			string+=" -f BAM -g hs -n "
			tmp= str(vector)+'_'+(solution)+'_'+(mymark)+'_'+(rep)
			string+= tmp+ "  -B -q 0.01 --broad --nomodel  "
			#string+= "> " +tmp+'.log ' #+ " 2>&1&"
		print ran,p,data[irow][-3]
		fl=open(tmp+'.log','a')
		#print os.getcwd()
		#print os.path.isfile(tmp+'.log')
		#print string
		process = subprocess.Popen(args=map(str,string.split()), stdout=fl,stderr=fl)

		fl.close()
		os.chdir('../')
		ran+=1
		if np.mod(ran,20)==0:
			time.sleep(6600)

		#print string, ran





