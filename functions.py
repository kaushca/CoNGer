#! /usr/bin/python2.7

import os
import sys
import subprocess
import shlex
from multiprocessing import Process, Manager, Pool
import random

scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))

def libsize(bam,bed,outDir,sample_names):
	tempDir=os.path.join(outDir,"temp_lib")
	os.mkdir(tempDir)
	inF_count=len(bam)
	sizeF=open(os.path.join(outDir,"libsize.txt"),"w")
	sizeF.write("\t".join(sample_names)+"\n")
	
	if (inF_count%4 >0):
		it=(inF_count/4)+1
	else:
		it=(inF_count/4)
	
	def temp_libsize(bam,bed,outF,name):
		#subprocess.call("intersectBed -abam %s -b %s > %s" %(bam,bed,outF+".bam"),shell=True) # ontarget reads are selected
		subprocess.call("samtools view -F 4 %s |wc -l > %s" %(bam,outF+".size"),shell=True) # skipping unmapped reads from the bam when viewing (-F 4)
		

	i=0
	k=0
	while(k<it):
		if(inF_count<((k+1)*4)):
			j = 0
			while(j<inF_count%4):
				temp_libsize(bam[(k*4)+j],bed,os.path.join(tempDir,sample_names[(k*4)+j]),sample_names[(k*4)+j])
				j=j+1
		else:
			pr1 = Process(target = temp_libsize, args=(bam[(k*4)+i],bed,os.path.join(tempDir,sample_names[(k*4)+i]),sample_names[(k*4)+i]))
			i=i+1
			pr2 = Process(target = temp_libsize, args=(bam[(k*4)+i],bed,os.path.join(tempDir,sample_names[(k*4)+i]),sample_names[(k*4)+i]))
			i=i+1
			pr3 = Process(target = temp_libsize, args=(bam[(k*4)+i],bed,os.path.join(tempDir,sample_names[(k*4)+i]),sample_names[(k*4)+i]))
			i=i+1
			pr4 = Process(target = temp_libsize, args=(bam[(k*4)+i],bed,os.path.join(tempDir,sample_names[(k*4)+i]),sample_names[(k*4)+i]))
			i=0
		
			pr1.start()
			pr2.start()
			pr3.start()
			pr4.start()
		
			pr1.join()
			pr2.join()
			pr3.join()
			pr4.join()
		
		k = k+1
	
	size=map(str,range(0,inF_count))
	for i in range(0,inF_count):
		f=open(tempDir+"/"+sample_names[i]+".size")
		size[i]=f.readline().rstrip()
		f.close()
	sizeF.write("\t".join(size)+"\n")
	
	sizeF.close()
	
	subprocess.call("rm -rf %s" %(tempDir),shell=True)
	
	return size

def getMeanCoverage(bam,bed,outDir,sample,gm,libsize):
	print("Calculate mean coverage...")
	
	fileName = sample+".coverage"
	
	tempF=os.path.join(outDir,fileName)

	subprocess.call("coverageBed -abam %s -b %s -d > %s" %(bam, bed,tempF),shell=True)	
	inF = open(tempF)
	outF = open(outDir+"/"+sample+".mean","w")

	tot = 0
	count = 0
	eps = 10

	for line in inF:
		l = line.rstrip().split("\t")
	
		tot = tot + float(l[5])
		if (float(l[5]) ==0):
			count = count + 1
	
		if (int(l[2]) - int(l[1]) == int(l[4])):
			#if count > 10:
			if count > -1:
				chr = str(l[0])
				start  = str(l[1])
				end = str(l[2])
				if count<int(l[4]):
					m=round(float(tot)/(int(l[4])-count),0)
				else:
					m=0
				c = m*gm/libsize
				outF.write(chr+"\t"+start+"\t"+end+"\t"+l[3]+"\t"+str(c)+"\n")
			tot = 0
			count =0
	
			
	inF.close()
	outF.close()
	
	subprocess.call("rm %s" %(tempF),shell=True)

def calculateCoverage(inF,bed,outDir,all_samples,tot_gm,size):
	inF_count=len(inF)
	
	if (inF_count%4 >0):
		it=(inF_count/4)+1
	else:
		it=(inF_count/4)
		
	i=0
	k=0
	outD=os.path.join(outDir,"temp")
	os.mkdir(outD)
	
	while(k<it):
		if(inF_count<((k+1)*4)):
			j = 0
			while(j<inF_count%4):
				getMeanCoverage(inF[(k*4)+j],bed,outD,all_samples[(k*4)+j],tot_gm,size[(k*4)+j])
				j=j+1
		else:
			pr1 = Process(target = getMeanCoverage, args=(inF[(k*4)+i],bed,outD,all_samples[(k*4)+i],tot_gm,size[(k*4)+i]))
			i=i+1
			pr2 = Process(target = getMeanCoverage, args=(inF[(k*4)+i],bed,outD,all_samples[(k*4)+i],tot_gm,size[(k*4)+i]))
			i=i+1
			pr3 = Process(target = getMeanCoverage, args=(inF[(k*4)+i],bed,outD,all_samples[(k*4)+i],tot_gm,size[(k*4)+i]))
			i=i+1
			pr4 = Process(target = getMeanCoverage, args=(inF[(k*4)+i],bed,outD,all_samples[(k*4)+i],tot_gm,size[(k*4)+i]))
			i=0
		
			pr1.start()
			pr2.start()
			pr3.start()
			pr4.start()
		
			pr1.join()
			pr2.join()
			pr3.join()
			pr4.join()
		
		k = k+1

def create_ratios(inF,control,outF):
	inF=open(inF)
	allControl=open(os.path.join(control,"gc_normalize.txt"))
	control=open(os.path.join(control,"control.txt"))
	rControl=open(os.path.join(outF,"control_ratios.txt"),"w")
	rTumour=open(os.path.join(outF,"tumour_ratios.txt"),"w")
	
	hT=inF.readline() #header
	hC=allControl.readline()
	control.readline()
	
	rControl.write(hC)
	rTumour.write(hT)
	
	for line in control:
		l=line.rstrip().split("\t")
		if float(l[4])==0:
			#removing zero locations of control from tumour
			inF.readline()
			allControl.readline()
		else:
			allT=map(float,inF.readline().rstrip().split("\t")[4:len(hT.split("\t"))]) #Coverage values of Tumours
			allC=map(float,allControl.readline().rstrip().split("\t")[4:len(hC.split("\t"))]) #coverage values of controls
			
			def divide(x,y):
				return (x/y)
			
			rT = map(divide,allT,[float(l[4])]*len(allT))
			rC = map(divide,allC,[float(l[4])]*len(allC))
			
			s = "\t".join(l[0:4])+"\t"+"\t".join(map(str,rT))+"\n"
			rTumour.write(s)
			
			s = "\t".join(l[0:4])+"\t"+"\t".join(map(str,rC))+"\n"
			rControl.write(s)
	
	inF.close()
	allControl.close()
	control.close()
	rControl.close()
	rTumour.close()

def divideExon(start,end,minLen=100):
	tempStart=start
	tempEnd=tempStart+minLen
	
	origins=[]
	terminals=[]
	
	while(tempEnd<=end):
		origins.append(tempStart)
		terminals.append(tempEnd)
		
		tempStart=tempEnd
		tempEnd=tempStart+minLen

	if(tempStart<end):
		origins.append(tempStart)
		terminals.append(end)
	
	return(origins,terminals)

def createBed(inFile,outFile):
	inFile= open(inFile)
	outFile = open(outFile,"w")
	
	for line in inFile:
		l = line.rstrip().split("\t")
		chrom = l[0]
		start = l[1]
		end = l[2]
		ID = l[3]
	
		origins,terminals = divideExon(int(start),int(end),100)
		
		rep = len(origins)
		
		for i in range(0,rep):
			outFile.write(chrom+"\t"+str(origins[i])+"\t"+str(terminals[i])+"\t"+ID+"\n")
	
	inFile.close()
	outFile.close()


def calculateGC(fasta,bed,outF):
	
	print("Calculating GC content of the targets...")
	
	subprocess.call("fastaFromBed -fi %s -bed %s -fo %s -tab" %(fasta,bed,os.path.join(outF,"tempGC.txt")),shell=True)
	
	bed =open(bed)
	fastabed=open(os.path.join(outF,"tempGC.txt"))
	gc=open(os.path.join(outF,"gc.txt"),"w")
	
	bedline=bed.readline().rstrip()
	seq=fastabed.readline().rstrip()
	if len(seq)>0:
		seq=seq.split("\t")[1]
	
	while 1:
		g = seq.count("g") + seq.count("G")
		a = seq.count("a") + seq.count("A")
		t = seq.count("t") + seq.count("T")
		c = seq.count("c") + seq.count("C")
		up = g+c
		down = max(up+a+t,1)
		ratio = up*100.0/down
		gc.write(bedline+"\t"+str(up)+"\t"+str(down)+"\t"+str(ratio)+"\n")
		bedline = bed.readline().rstrip()
		seq=fastabed.readline().rstrip()
		if len(seq)>0:
			seq=seq.split("\t")[1]
		else:
			break
				
	subprocess.call("rm %s" %(os.path.join(outF,"tempGC.txt")),shell=True)
	
	bed.close()
	fastabed.close()
	gc.close()
	
def getFileNames(tumor,control=""):
	tumor=open(tumor)

	tumorFiles=[]
	tumorNames=[]
	for line in tumor:
		l = line.rstrip()
		tumorFiles.append(l)
		fname=os.path.split(l)[1].split("_merged.tmp.realign.recal.bam")[0]
		tumorNames.append(fname)

	tumor.close()
	if control!="":
		control=open(control)
		controlFiles = []
		controlNames = []
		for line in control:
			l = line.rstrip()
			controlFiles.append(l)
			fname = os.path.split(l)[1].split("_merged.tmp.realign.recal.bam")[0]
			controlNames.append(fname)
	
			control.close()
			return(tumorFiles,tumorNames,controlFiles,controlNames)
		
	return(tumorFiles,tumorNames)

def gc_normalize(gc,libF,outF,scriptPath):	
	rScriptName = os.path.join(scriptPath, "gc_loess.R")
	args = shlex.split("Rscript %s %s %s %s" %(rScriptName,gc,libF,outF))
	rscr = subprocess.call(args)
#	subprocess.call("Rscript %s %s %s %s" %(rScriptName,gc,libF,outF),shell=True)

def addExtension (x,xx):
	s = x+".mean"
	y = os.path.join(xx,s)
	
	return y

def createDirPath(x,y): return os.path.join(x,y)

def getUniCnt(details):
	nCol=len(details[0])
	codes=[]
	for row in details:
		codes.append("".join(row[2:nCol]))
	
	#for i in range(2,nCol):
	#	uniqueVal.append(set[row[i]])
	
	uniqueVal=list(set(codes))
	codeCnt=[]
	
	for i in uniqueVal:
		codeCnt.append(codes.count(i))

	return (uniqueVal,codeCnt,codes)

def fragCounts(x,y,val):
	values=[j for j,i in enumerate(x) if i==val and y[j]]
	
	return(values,len(values))
	
def extractInfo(info,fragHigh=340,fragLow=325,minSamples=5,nPool=15):
	info = open(info)
	infoMat = []
	colNames=info.readline().rstrip().split("\t")
	sampleIndex=colNames.index("sample")
	fragIndex=colNames.index("fragSize")
	batchIndex=colNames.index("flowCell")
	prepInd=colNames.index("prepPlate")
	laneInd=colNames.index("seqPool")
	
	for line in info:
		l = line.rstrip().split("\t")
		infoMat.append(l)
	
	largeFrag=[int(row[fragIndex])>=fragHigh for row in infoMat] #select >340 samples as good samples
	lowFrag=[int(row[fragIndex])>=fragLow and int(row[fragIndex])<fragHigh for row in infoMat]  # select > 325 samples as an alternative
	badFrag=[int(row[fragIndex])<fragLow for row in infoMat]
	
	lanes=[row[laneInd] for row in infoMat]
	uniLanes=list(set(lanes))
	
	lanePlates=["___".join([row[prepInd],row[laneInd]]) for row in infoMat]
	uniLanePlate=list(set(lanePlates))
	
	lP=[ x.split("___") for i,x in enumerate(uniLanePlate)]
	s=[[row[1] for row in lP].count(i) for i in uniLanes]
	
	samples=[]
	for i in uniLanePlate:
		#longFrag=[j for j,x in enumerage(lanePlates) if x==i and largeFrag[j]]
		#longCnt=len(longFrag)
		k=i.split("___")
		for j,x in enumerate(uniLanes):
			if x==k[1]:
				l=j
		
		if (s[l]==1):
			longFrag,longCnt=fragCounts(lanePlates,largeFrag,i)
			shortFrag,shortCnt=fragCounts(lanePlates,lowFrag,i)
			restFrag,restCnt=fragCounts(lanePlates,badFrag,i)			
		else:
			plates=[row[0] for row in lP if row[1]==k[1]]
			for p in plates:
				pLabel="___".join([p,k[1]])
				if (pLabel==i and lanePlates.count(i)>=minSamples):
					longFrag,longCnt=fragCounts(lanePlates,largeFrag,pLabel)
					shortFrag,shortCnt=fragCounts(lanePlates,lowFrag,pLabel)
					restFrag,restCnt=fragCounts(lanePlates,badFrag,pLabel)
				elif pLabel==i:
					longFrag,longCnt=fragCounts(lanes,largeFrag,k[1])
					shortFrag,shortCnt=fragCounts(lanes,lowFrag,k[1])
					restFrag,restCnt=fragCounts(lanes,badFrag,k[1])	
				else: continue	
		
		if longCnt>=nPool:
			r=random.sample(longFrag,nPool)
		elif shortCnt+longCnt>=nPool:
			rr=[longFrag,random.sample(shortFrag,nPool-longCnt)]
			r = [e for row in rr for e in row]
		elif lanePlates.count(i)>=nPool:
			rr=[longFrag,shortFrag,random.sample(restFrag,nPool-longCnt-shortCnt)]
			r = [e for row in rr for e in row]
		else:
			if lanePlates.count(i)<2:
				print "Not enough samples to pool in plate "+k[0]+" and lane "+k[1]+" exiting from software"
				sys.exit(1)
			elif lanePlates.count(i)<minSamples:
				print "Less than "+str(minSamples)+" samples in plate "+k[0]+" and lane "+k[1]
			rr = [longFrag,shortFrag,restFrag]
			r = [e for row in rr for e in row]
			
		samples.append(r)
		
	return (uniLanePlate,samples,infoMat)


	
