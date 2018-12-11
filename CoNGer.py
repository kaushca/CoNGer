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
		self.parser = argparse.ArgumentParser("Germline Copy Number detection tool for Haloplex data")
		self.parser.add_argument('-t','--tumor',help = 'File containing collection of paths to  tumour BAM files, each in new line',required=True)
		self.parser.add_argument('--info',help="Tab delimited file containing sample information (sample name, library average fragment size, flow cell, library prep plate and sequencing pool)")
		self.parser.add_argument('-c','--control',help = 'Path to normalized samples to create the required control -either a directory containing all subdirectories for each batch or a file with normalized samples (with no batch based info)',required=True)
		self.parser.add_argument('-b','--bed',help = 'Bed definition of covered regions',required=True)
		self.parser.add_argument('--out',help = 'Output directory',required=True)
		self.parser.add_argument('--fasta',help = 'FASTA file of the reference genome',required=True)
		self.parser.add_argument('--gc',help = 'GC content %% in targets [optional]',default =False)
		self.parser.add_argument('--samples',help='Provide tumour sample names separated by spaces [optional]',nargs='+')
		self.parser.add_argument('--lowFragments',help='Lower threshold for fragment count (used in creating the pooled normal and frequency cutoffs for CNV calls) [optional]',default=325)

		try:
			args = self.parser.parse_args()
		except:
			self.parser.print_help()
			sys.exit(0)
		
		if args.tumor:
			self.tumor=args.tumor
		if args.control:
			if os.path.isfile(args.control):
				#cmd=shlex.split("awk '{print NF}' %s | tail -n 1" %(args.control))
				#nC,err=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
				#self.normalizeC=False
				#if(int(nC)>3):
				print("No batch based control...")
				self.normalizeC=True
			else:
				if(args.info):
					self.normalizeC=False
					self.info=args.info
				else:
					print("Information on batches of the tumour need to be supplied")
					sys.exit(1)
			self.control=args.control
		if args.bed:
			self.bed=args.bed
		if args.out:
			self.outF=args.out
		if args.fasta:
			self.fasta=args.fasta
		if args.gc:
			self.gcT = True
			self.gc=args.gc
		else:
			self.gcT= False
		if args.samples:
			if len(self.tumor)!=len(args.samples):
				self.parser.error('Number of samples provided should equal to the number of sample names')
			self.samples=args.samples
		else:
			self.samples = []
		if args.lowFragments:
			self.fragLow=args.lowFragments


def main():
	options = Options()
	
	tumor = options.tumor
	control = options.control
	normalizeC=options.normalizeC
	outDir = options.outF
	bed = options.bed
	sample_names = options.samples
	fasta = options.fasta
	fragThresh=options.fragLow
	
	
	try:
		os.mkdir(outDir)
	except:
		print("cannot create folder '%s'" %outDir)
		print("if folder already exist, please specify other folder")
		sys.exit(1)
		
	tumorFiles,tumorNames=getFileNames(tumor)	

	if len(sample_names)==0:
		sample_names=tumorNames
	
	print("Analyzing following tumor samples...")
	print("\n".join(sample_names))
	
	all_samples = tumorNames
	inF = tumorFiles
	

	size = map(float,libsize(inF,bed,outDir,all_samples))
	print("sizes calculated")
	tot_gm=1000000000 #normalizing based on million reads instead of geometric mean (this makes sure control and tumour samples have same libsize)
	#and calculating RPKM i.e. (# reads * 10^9)/(target length*total mapped reads) = (average reads per targe base*10^9)/total mapped reads
	
	createBed(bed,os.path.join(outDir,"targets.bed"))
	
	bed = os.path.join(outDir,"targets.bed")
	
	if (options.gcT==True):
		gc = options.gc
	else:
		calculateGC(fasta,bed,outDir)
		gc = os.path.join(outDir,"gc.txt")
		
	print("Calculating mean coverage...")
	
	calculateCoverage(inF,bed,outDir,all_samples,tot_gm,size)
	
	
	####moving tumor mean coverage files
	tumorF = os.path.join(outDir,"tumorRaw.txt")
	tumorTemp = os.path.join(outDir,"tempTumor")
	subprocess.call("mkdir %s" %(tumorTemp),shell=True)
	x=map(addExtension,sample_names,[os.path.join(outDir,"temp")]*len(sample_names))
	subprocess.call("mv %s %s" %(" ".join(x),tumorTemp),shell=True)
	s = "\t".join(['Chr','Start','End','ID']+sample_names)
	args = shlex.split("bash %s %s %s %r" %(os.path.join(scriptPath,"concat_files.sh"),tumorTemp,tumorF,s))
	subprocess.call(args)
	
	
	
	pr1 = Process(target = gc_normalize, args=(gc,tumorF,os.path.join(outDir,"tumorGC.txt"),scriptPath))###GC Normalization of tumor
	pr1.start()
	pr1.join()
	subprocess.call("rm -r %s %s" %(tumorTemp,os.path.join(outDir,"temp")),shell=True)
	
	if normalizeC==False:
		info = options.info	
		infoFile=open(info)
		line=infoFile.readline().rstrip().split("\t")
		prepInd=line.index("prepPlate")
		laneInd=line.index("seqPool")
		fragIndex=line.index("fragSize")
		tmpFile=open(os.path.join(outDir,"temp.txt"),"w")
		tmpFile.write(control)
		tmpFile.close()
		control=os.path.join(outDir,"temp.txt")
		
	for sample in all_samples:
		if normalizeC==False:
			print("Selecting control based on batch")
			cmd=shlex.split("grep %s %s" %(sample,info))
			i=(subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate())
			i=str(i[0]).rstrip().split("\t")
			x="___".join([i[prepInd],i[laneInd]])
			x=x.replace("/","")
			x=x.replace(" ","-")
			control=open(control).readline().rstrip()
			controlPath=os.path.join(control,x,"controlGC.txt")
			print(controlPath)
			subprocess.call("cp %s %s" %(controlPath,os.path.join(outDir,x+".txt")),shell=True)
			controlGC=controlPath
		else:
			controlGC=control
			subprocess.call("cp %s %s" %(control, os.path.join(outDir,"controlGC.txt")),shell=True)
		
	
	x=map(createDirPath,[outDir]*len(sample_names),sample_names)
	subprocess.call("mkdir %s" %(" ".join(x)),shell=True)
	
	rScriptName = os.path.join(scriptPath,"analyseTumor.R")
	args = shlex.split("Rscript %s %s %s %s" %(rScriptName,os.path.join(outDir,"tumorGC.txt"),controlGC,outDir))
	rscr = subprocess.call(args)

if __name__ == "__main__":
	main()	

#END
