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

syst = ["","_jerup","_jerdown","_jesup","_jesdown","_umetup","_umetdown","_lesup","_lesdown","_puup","_pudown","_btagup","_btagdown"]

#inputdir = "/afs/cern.ch/work/c/chasco/WW_8/Addon/OUT_v8_WR/"
inputdir = "/afs/cern.ch/work/c/chasco/OCT19_p66_8/OUT_v2/"
TeV8 = True
os.system("mkdir v2")
#outputdir = inputdir+"v8/"
inputdirlist = os.listdir(inputdir)
inputdirlistroot = []
for ii in inputdirlist:
	if (".root" in ii) and ('BKGD.root' not in ii) and ('BKGDandZZ.root' not in ii) and ('ZHcombo.root' not in ii):
		inputdirlistroot.append(ii)
print inputdirlistroot

WEIGHTING = "Eweight*XS*BR*LUM*(1/NGE)*(B2/B3)"#*Wscale*Zscale"
WEIGHTINGT = "WT4" #application sample only
CUTTING = "(Zmetphi > 2.6)*(REDmet > 110)*((met/zpt)>0.8)*((met/zpt)<1.2)*(mass > 76)*(mass < 106)*(pBveto>0.0)" #application sample only
CUTTINGT = "(training4<0.1)" #application sample only
EE = "(finstate > 1.5)"
MM = "(finstate < 1.5)"

VARS = ['LikelihoodZH125vsBKGDandZZ','BDTZH125vsBKGDandZZ','MLPZH125vsBKGDandZZ','SVMZH125vsBKGDandZZ','CFMlpANNZH125vsBKGDandZZ','TransMass3','l1pt']
# bin = 10
# _min = 0.0
# _max = 1.0

treeNameLEP = [["muons","MM",MM],["electrons","EE",EE]]

for LEP in treeNameLEP:

	treeName = LEP[0]

	for v in VARS:
		bin = 20
		_min = 0.0
		_max = 1000.0
		if ("Likelihood" in v):
			print "Likelihood"
			bin = 20
			_min = 0.0
			_max = 1.0
		if ("BDT" in v):
			print "BDT"
			bin = 20
			_min = -1.0
			_max = 1.0
		if ("SVM" in v):
			print "SVM"
			bin = 20
			_min = 0.268
			_max = 0.277
		if ("MLP" in v):
			print "MLP"
			bin = 20
			_min = -0.7
			_max = 1.3
		if ("CFMlpANN" in v):
			print "CFMlpANN"
			bin = 20
			_min = 0.0
			_max = 0.9




		FileFF=TFile.Open("v2/"+v+"_"+LEP[1]+str(TeV8)+".root","RECREATE")
		FallaDirectory = FileFF.mkdir(treeName,"histogramcrackers")
		FallaDirectory.cd()

		HData = TH1F("data_obs","",bin,_min,_max)
		HData.Sumw2()

		for f in inputdirlistroot:

			n = f.replace(".root","")
			fin = TFile.Open(inputdir+f,"READ")
			tin = fin.Get("tmvatree")
			#print tin.GetEntries()
			FallaDirectory.cd()

			if ("Data" in f):

				exec(n+'=TH1F("'+n+'","",bin,_min,_max)')
				exec(n+'.Sumw2()')
				tin.Draw(v+">>"+n,WEIGHTING+"*"+CUTTING+"*(sys>0.0)"+"*"+LEP[2])
				exec('HData.Add('+n+')')

			else:
				for s in syst:
					#print "*"*20
					#print(n+s+'=TH1F("'+n+s+'","",bin,_min,_max)') 
					exec(n+s+'=TH1F("'+n+s.replace("up","Up").replace("down","Down")+'","",bin,_min,_max)')  #make histograms for each file
					exec(n+s+'.Sumw2()')
					#if (s == ""):
					if "Likelihood" in v:
					#print "Likelihood! "*20
						tin.Draw(v+">>"+n+s.replace("up","Up").replace("down","Down"),WEIGHTING+"*"+WEIGHTINGT+"*"+CUTTING+"*(sys"+s+">0.5)*"+LEP[2]+"*"+CUTTINGT)
					else:
						tin.Draw(v+">>"+n+s.replace("up","Up").replace("down","Down"),WEIGHTING+"*"+CUTTING+"*(sys"+s+">0.5)*"+LEP[2])
					#else:
						#tin.Draw(v+">>"+n+s,WEIGHTING+"*"+CUTTING+"*(sys"+s+">0.5)*"+LEP[2])
					#print WEIGHTING+"*"+CUTTING+"*sys"+s+"*"+LEP[2]
					exec(n+s+'.Write()')

			fin.Close()

		HData.Write()
		FileFF.Close()

