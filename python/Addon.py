import os
import sys
sys.argv.append("-b")
print '\nLoading ROOT ... \n\n'
#import ROOT
#from ROOT import TFile, TTree, TLorentzVector, kTRUE, TMath, TNtuple, gRandom, TCanvas, TH2F
from ROOT import *
import math
print 'ROOT loaded.'

import numpy
import array
import random

SampleList7 = []
SampleList8 = []

SampleList7.append(['ZH105',0.5449,0.100974])
SampleList7.append(['ZH115',0.4107,0.100974])
SampleList7.append(['ZH125',0.3158,0.100974])
SampleList7.append(['ZH135',0.2453,0.100974])
SampleList7.append(['ZH145',0.1930,0.100974])

SampleList8.append(['ZH105',0.6750,0.100974])
SampleList8.append(['ZH115',0.5117,0.100974])
SampleList8.append(['ZH125',0.3943,0.100974])
SampleList8.append(['ZH135',0.3074,0.100974])
SampleList8.append(['ZH145',0.2424,0.100974])

def scale2data(afile,chan,energy8): #finstate
	EE = (chan > 1.5)
	MM = (chan < 1.5)
	wscale = 1.0
	zscale = 1.0
	if ('SingleT' in afile) or ('TTJets' in afile) or (('W' and 'Jets') in afile) or ('WW' in afile):
		if EE:
			wscale = 2.95
			if energy8:
				wscale = 0.72
		if MM:
			wscale = 2.01
			if energy8:
				wscale = 0.70
	else:

		if ('WZ' in afile) or ('ZZ' in afile):
			if EE:
				zscale = 1.042
				if energy8:
					zscale = 1.025
			if MM:
				zscale = 1.073
				if energy8:
					zscale = 1.034

	return [wscale,zscale]

def EWcorrection(afile,Zpt,Hpt):
	Tpt = 0.0
	Hpt = 0.0
	Zpt = 0.0
	corr = 1.0

	if (Zpt > Hpt):#trailing boson
		Tpt = Hpt
	else:
		Tpt = Zpt

	if (Hpt >= 50.0) and ("ZH" in afile):
		corr = 0.94 - 1.0*(0.2-0.068)/400.0*(Hpt-50.0)
	if (Tpt >= 50.0):
		if ("ZZ" in afile): #trailing boson pt
			corr = 1.0 + (0.55 - 0.071*Tpt)/100.0
		if ("WZ" in afile): #trailing boson pt
			corr = 1.0 + (1.9 - 0.037*Tpt)/100.0

	return corr


def XSarc(afile,XS,BR,samplelist):
	out = [XS,BR]
	for ss in samplelist:
		if ss[0] in afile:
			out = [ss[1],ss[2]]
	return out


def AddMore(inputdir,inputfile,inputtree,outputdir,TeV8,SampleList):
	FInA = TFile.Open(inputdir + inputfile,"")
	TInA= FInA.Get(inputtree)
	tdir = FInA.Get(inputtree)
	x = tdir.GetListOfBranches() #retrieve all old variables
	tagalongvariables = []
	for y in x:
		tagalongvariables.append(y.GetName())
	print tagalongvariables

	addons = ['Wscale','Zscale','EWCorr','XS2','BR2']

	for var in (tagalongvariables + addons):
  		exec(var+' = array.array(\'f\',[0])')
	
	FOut = TFile.Open(outputdir + inputfile,"RECREATE")
	TOut = TTree("tmvatree","tmvatree")
	for var in (tagalongvariables + addons): #add tag-alongs to output tree
		exec("TOut.Branch('"+var+"', "+var+",'"+var+"/F')")

	N = TInA.GetEntries()
	print N, "INITIAL N"

	for n in range(N):
		if (n%10000 == 1) or (n+1==N):
			print inputfile+":   "+str(n+1) +' of '+str(N) +' events evaluated.'
		TInA.GetEntry(n)

		for var in tagalongvariables:
			exec(var+'[0] = TInA.'+var)

		Wscale[0] = scale2data(inputfile,finstate[0],TeV8)[0]
		Zscale[0] = scale2data(inputfile,finstate[0],TeV8)[1]
		if ("ZH" in inputfile) or ("ZZ" in inputfile) or ("WZ" in inputfile):
			EWCorr[0] = EWcorrection(inputfile,zpt,met)
		else:
			EWCorr[0] = 1.0

		XSa = XSarc(inputfile,XS[0],BR[0],SampleList) #use old ZH cross sections, but keep new ones, so add XS2 and BR2
		XS2[0] = XSa[0]
		BR2[0] = XSa[1]

		TOut.Fill()
	NN = TOut.GetEntries()
	print NN, "FINAL N"
	FOut.Write("", TObject.kOverwrite)
	FOut.Close()


inputdir = "/afs/cern.ch/work/c/chasco/WW_7/"
En8 = False

if En8:
	SampleList = SampleList8
else:
	SampleList = SampleList7

outputdir = inputdir + "Addon2/"
os.system("mkdir "+outputdir)

inputdirlist = os.listdir(inputdir)

inputdirlistroot = []
for x in inputdirlist:
	if (".root" in x) and ('BKGD.root' not in x) and ('BKGDandZZ.root' not in x) and ('ZHcombo.root' not in x):
		inputdirlistroot.append(x)

print inputdirlistroot

for x in inputdirlistroot:
	AddMore(inputdir,x,"tmvatree",outputdir,En8,SampleList)