#python profile-heatmap.py <profile1> <profile2> <rpkmfile> <outputname>

import sys, os
rnafile=str(sys.argv[-2])
print str(os.path.basename(rnafile))
import numpy as np
import matplotlib.pyplot as plt
def write_data(fl,data):
	fwrite=open(fl,'w')

	for i in data:
		tmpstr=''
		tmpstr+= str(i)+'\t'
	#	print tmpstr
		print >> fwrite, tmpstr[:-1]

	fwrite.close()

def data2arr(data):
   for row in data:
      for i in row:
            i=float(i)
      row=np.array(row)
   data=np.array(data)

   return data

def read_dat(file,*spl):
    if len(spl)!=1:
       spl=" "
    f=open(file,'r')
    a=[]
    if len(spl)>0:
        spl=str(spl[0])
    ind=0
    for  ln  in f:
       s = ln.split(spl)
       for i in xrange(len(s)):
            s[i]=s[i].strip()
            try :
                 s[i]=float(s[i])
            except:
               pass
       a.append(s)
    f.close()
    return  a

CGI=[]
cgifile=read_dat('/Users/ssaberim/profiles-novo/hg19v65_genes.CGI','\t')
for i in cgifile:
    CGI.append(i[0])
profile=[]
for j in sys.argv[1:-2]:
  profile.append(str(j))
rnafile=str(sys.argv[-2])
outfile=str(sys.argv[-1])
rna0=read_dat(rnafile,'\t')
data=[]
for i in profile:
   data.append(read_dat(i)[:-3])

gnrna=[]
rna=[]
for i in range(len(rna0[0])): #:###
  #if rna0[0][i] in CGI:
    gnrna.append(str(rna0[0][i]))
    rna.append(rna0[1][i])
rna=np.array(rna)
gnrna=np.array(gnrna)
ind=rna.argsort()[::-1]
rna=rna[ind]
#print rna[:10]
gnrna=gnrna[ind]
#for i in gnrna[:100]:
#        print i
ind=[]
mat=[]
mat0=[]
gn=[]


num=len(profile)
for j in xrange(num):
    mat.append([])
    gn.append([])
    for i in data[j]:
         mat[j].append(i[1:])
         gntmp=str(i[0])
         gn[j].append(gntmp)
    mat[j]=np.array(mat[j])
    #print 'aha', len(mat[j]),len(gn[j])

a=[]
for k in xrange(num):
    inds=[]
    for i in gnrna:
        try:
            inds.append(gn[k].index(i))
        except:
            pass
    inds=np.array(inds)
    gn[k]=np.array(gn[k])
    #print np.max(mat[k]),len(inds),k
    mat0.append(mat[k][inds,:]) #np.random.randint(low=0,high=len(mat[k]),size=10000)])
    mat0[k]=np.array(mat0[k])
    #plt.figure()
    #plt.plot(mat0[k][100])
    #plt.plot(np.mean(mat[k][:100],axis=0))
    #plt.plot(np.mean(mat0[k],axis=0))
    #plt.plot(np.mean(mat0[k]*len(mat0[k]),axis=0))
    #print np.mean(mat[k]), np.mean(mat0[k]),  len(inds),inds
    #plt.plot(np.mean(mat0[k],axis=0))
    #plt.show()

for k in range(num):
    a=np.amax(mat0[k],axis=1)
    for i in xrange(len(mat0)):
        if a[i]!=0:
          mat0[k][i]/=a[i]
sub=num+1
tmp=0
sz=8

plt.figure()
for i in range(num):
   tmp+=1
   plt.subplot(100+sub*10+tmp)
   plt.imshow(mat0[i],aspect='auto',cmap='RdYlBu',interpolation='nearest',vmin=0,vmax=20.)
   vec=np.mean(mat0[i],axis=0)
   plt.xticks(rotation=90)
   plt.xticks((range(0,201,20)),range(-2000,2001,400),fontsize=sz)
   plt.yticks([])
   plt. xlabel('bp',fontsize=sz)
   t= str(sys.argv[tmp]).split('/')[-1].split('_')[0]
   if str('nuc') in str(sys.argv[tmp]).split('/')[-1]:
       t+='_nucfree'
   plt.title( t,fontsize=10)
   plt.ylabel("All genes ranked",fontsize=sz)
   if 'CpG' in str(os.path.basename(str(sys.argv[-1]))) or 'CPG' in str(os.path.basename(str(sys.argv[-1]))):
       t+='_CGI'
   else:
       t+='_All'
   #plt.figure()
   #plt.plot(vec)
   #plt.show()
   write_data(t+'.txt',vec)
tmp+=1
#plt.colorbar()
plt.subplot(100+sub*10+tmp)
print rna[0], len(rna)
plt.plot(rna[::-1],range(len(rna)),'-')
plt.ylim(0,len(rna))
plt.yticks([])
plt.xscale('log')
plt.xlabel("RNA expression (RPKM)",fontsize=sz)
plt.xticks(fontsize=sz-2)
plt.grid()
plt.savefig(outfile+'.pdf', bbox_inches='tight')
plt.close()

