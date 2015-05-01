import os, sys 
import commands
import read
import write
import numpy as np
import time
import subprocess

table=read.read_dat("peaks.txt")
data=table
cells=["D83V","K4E", "pcDNA","WT","Y69H"]
marks=["H3K27ac","H3K4me3","H3K9me3","V5","input"]
ran=0 
for irow in xrange(1,len(data)):
	row=data[irow]
	vector=row[-3]
	solution=row[-1][:-2]
	mymark=""
	if "input" in row[-2]:
				mymark="input"
	else:
		mymark=row[-2]
	if '1' in row[-4]:
		rep='1'
	else:
		rep='2'
	bam=row[5]
	if 'WT' not in vector and 'input' not in mymark and 'H3K4' not in mymark:
		for ind in xrange(-1,4):
			try:
				ii=data[irow+ind][-3]
				if 'WT' in ii:
					wt=data[irow+ind][5]
			except:
				print 'end '
			if 'WT' in ii:
				wt=data[irow+ind][5]
		bam=os.popen("ls -f " + bam +"*.bam")
		bam=bam.readline()[:-1]
		wt=os.popen("ls -f " + wt +"*.bam")
		wt=wt.readline()[:-1]
		if not os.path.isdir(mymark):
			os.makedirs(mymark)
		os.chdir(mymark)
		#print os.getcwd()
		if "V5" in mymark :
			string=" nohup macs2 callpeak -t "+bam +" -c " + wt
			string+=" -f BAM -g hs -n "
			tmp= str(vector)+'_'+(solution)+'_'+(mymark)+'_'+(rep)
			string+= tmp + "  -q 0.05  --SPMR --nomodel " 
			#string+= "> "  +tmp+'.log ' #" 2>&1&"
		#elif "H3K9" in mymark:
		#	string=" nohup macs2 callpeak -t "+bam +" -c " + wt
		#	string+=" -f BAM -g hs -n "
		#	tmp= str(vector)+'_'+(solution)+'_'+(mymark)+'_'+(rep)
		#	string+= tmp+ "  -B -q 0.01 --broad --nomodel  "
		#	#string+= "> " +tmp+'.log ' #+ " 2>&1&"
		print ran
		fl=open(tmp+'.log','a')
		print os.getcwd()
		print os.path.isfile(tmp+'.log')
		print string
		
		
		process = subprocess.Popen(args=map(str,string.split()), stdout=fl,stderr=fl)
		fl.close()
		os.chdir('../')
		ran+=1
		if np.mod(ran,20)==0:
			time.sleep(3600)

		#print string, ran





