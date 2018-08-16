#! /usr/bin/python

import os
import argparse
import sys
import subprocess
import shlex
from functions import *
from multiprocessing import Process, Manager, Pool

#absolute script path
scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))

class Options:

	def __init__(self):
		self.parser = argparse.ArgumentParser("Normalizing controls")
		self.parser.add_argument('-c','--control',help = 'File containing collection of paths to  control BAM files, each in new line',required=True)
		self.parser.add_argument('-b','--bed',help = 'Bed definition of covered regions',required=True)
		self.parser.add_argument('--out',help = 'Output directory',required=True)
		self.parser.add_argument('--fasta',help = 'FASTA file of the reference genome',required=True)
		self.parser.add_argument('--info',help="Tab delimited file containing sample information (sample name, library average fragment size and flow cell are mandatory)")
		self.parser.add_argument('--gc',help = 'GC content %% in targets [optional]',default =False)
		self.parser.add_argument('--lowFragments',help='Lower threshold for fragment count (used in creating the pooled normal and frequency cutoffs for CNV calls) [optional]',default=325)
		self.parser.add_argument('--highFragments',help='Lower threshold for fragment count (used in creating the pooled normal) [optional]',default=340)
		self.parser.add_argument('--minSamples',help='Minimum number of samples to be pooled per group',default=5)
		self.parser.add_argument('--nPool',help='Number of samples to be pooled per group',default=15)
		
	
		args = self.parser.parse_args()
		
		if args.control:
			self.control=args.control
		if args.info:
			self.info_status=True
			self.info=args.info
		else:
			self.info_status=False
		if args.bed:
			self.bed=args.bed
		if args.out:
			self.outF=args.out
		if args.fasta:
			self.fasta=args.fasta
		if args.gc:
			self.gcT=True
			self.gc=args.gc
		else:
			self.gcT=False
		if args.lowFragments:
			self.fragLow=args.lowFragments
		if args.highFragments:
			self.fragHigh=args.highFragments
		if args.nPool:
			self.nPool=args.nPool
		if args.minSamples:
			self.minSamples=args.minSamples

def createPool(inF,bed,outDir,all_samples,tot_gm,gc):
	size = map(float,libsize(inF,bed,outDir,all_samples))	
		
	print("Calculating mean coverage...")
	
	calculateCoverage(inF,bed,outDir,all_samples,tot_gm,size)
	
	##moving control mean coverage files
	controlF = os.path.join(outDir,"controlRaw.txt")
	controlTemp = os.path.join(outDir,"tempCon")
	subprocess.call("mkdir %s" %(controlTemp),shell=True)
	x=map(addExtension,all_samples,[os.path.join(outDir,"temp")]*len(all_samples))
	subprocess.call("mv %s %s" %(" ".join(x),controlTemp),shell=True)
	s = "\t".join(['Chr','Start','End','ID']+all_samples)
	args = shlex.split("bash %s %s %s %r" %(os.path.join(scriptPath,"concat_files.sh"),controlTemp,controlF,s))
	subprocess.call(args)
	
	pr2 = Process(target = gc_normalize, args=(gc,controlF,os.path.join(outDir,"controlGC.txt"),scriptPath))####GC Normalization of control
	pr2.start()
	pr2.join()
	
	subprocess.call("rm -r %s %s" %(controlTemp,os.path.join(outDir,"temp")),shell=True)

def main():
	options = Options()
	
	control = options.control
	outDir = options.outF
	bed = options.bed
	fasta = options.fasta
	
	if (options.info_status==True):
		info = options.info
		uniLanePlates,sIndex,infoMat=extractInfo(info)
	
	os.mkdir(outDir)
	
	controlFiles,controlNames = getFileNames(control)
	createBed(bed,os.path.join(outDir,"targets.bed"))
	bed = os.path.join(outDir,"targets.bed")
	
	tot_gm=1000000000 #normalizing based on million reads instead of geometric mean (this makes sure control and tumour samples have same libsize)
	#and calculating RPKM i.e. (# reads * 10^9)/(target length*total mapped reads) = (average reads per targe base*10^9)/total mapped reads
	
	if (options.gcT==True):
		gc = options.gc
	else:
		calculateGC(fasta,bed,outDir)
		gc = os.path.join(outDir,"gc.txt")
	
	if (options.info_status==True):
		for j,x in enumerate(uniLanePlates):
			indices = sIndex[j]
			samplesX = [infoMat[i][0] for i in indices]
			idx = []
			for y in samplesX:
				for i, xx in enumerate(controlNames):
					if xx==y:
					 	idx.append(i)
			inF = [controlFiles[i] for i in idx]
			all_samples=[controlNames[i] for i in idx]
			x=x.replace("/","")
			x=x.replace(" ","-")
			if len(inF)>0:
				subDir=os.path.join(outDir,x)
				os.mkdir(subDir)
				f=open(os.path.join(subDir,"samples.txt"),"w")
				f.write("\n".join(all_samples))
				f.close()
				print("Creating Pool for "+x)
				createPool(inF,bed,subDir,all_samples,tot_gm,gc)
			
	else:
		all_samples = controlNames
		inF = controlFiles
		print("No batch based pooling")
		createPool(inF,bed,outDir,all_samples,tot_gm,gc)
	
if __name__ == "__main__":
	main()	

#END
