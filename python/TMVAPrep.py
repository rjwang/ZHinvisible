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
ROOT.gStyle.SetOptStat(0)

hadder = str(sys.argv[1])

whataboutsystematics = str(sys.argv[2]) #options: nonvar *nothing*=both

def scale2data(file,chan,energy8): #finstate
	EE = (chan > 1.5)
	MM = (chan < 1.5)
	wscale = 1.0
	zscale = 1.0
	if ('SingleT' in file) or ('TTJets' in file) or (('W' and 'Jets') in file) or ('WW' in file):
		if EE:
			wscale = 2.95
			if energy8:
				wscale = 0.72
		if MM:
			wscale = 2.01
			if energy8:
				wscale = 0.70
	else:

		if ('WZ' in file) or ('ZZ' in file):
			if EE:
				zscale = 1.042
				if energy8:
					zscale = 1.025
			if MM:
				zscale = 1.073
				if energy8:
					zscale = 1.034

	return [wscale,zscale]

# def MT(pt1,pt2,phi):
# 	MTout= math.sqrt(2*pt1*pt2*(1-math.cos(abs(phi))))
# 	return MTout

# def MT2(pt1,eta1,phi1,pt2,eta2,phi2,met,metphi):
# 	mu1 = TLorentzVector()
# 	mu2 = TLorentzVector()
# 	misg = TLorentzVector()
# 	v1 = TLorentzVector()
# 	v2 = TLorentzVector()
# 	Pi = math.pi
# 	mu1.SetPtEtaPhiM(pt1, eta1, phi1, 0)
# 	mu2.SetPtEtaPhiM(pt2, eta2, phi2, 0)
# 	misg.SetPtEtaPhiM(met,0,metphi,0)
# 	qq = 0.0
# 	n=8 #number of iterations in phi
# 	m=9 #number of iterations in met
# 	k=3.0 #how large can energy get w.r.t. met
# 	s11=0.0
# 	s22=0.0
# 	s12=0.0
# 	s21=0.0
# 	smin=0.0
# 	smax=0.0
# 	store = []
# 	for q in range(m+1):
# 		qq = q*k*met/(1.0*m)
# 		for a in range(n+1):
# 			v1.SetPtEtaPhiM(qq,0,a*Pi/(1.0*n),0)
# 			v2 = misg - v1
# 			s11 = MT(v1.Pt(),mu1.Pt(),v1.DeltaPhi(mu1))
# 			s22 = MT(v2.Pt(),mu2.Pt(),v2.DeltaPhi(mu2))
# 			s12 = MT(v1.Pt(),mu2.Pt(),v1.DeltaPhi(mu2))
# 			s21 = MT(v2.Pt(),mu1.Pt(),v2.DeltaPhi(mu1))
# 			smin = min(s11,s22,s12,s21)
# 			if (smin > smax):
# 				smax = smin
# 				store = [smax,s11,s22,s12,s21,q,a]
# 	return store

# def MT2(pt1,eta1,phi1,pt2,eta2,phi2,met,metphi):
# 	mu1 = TLorentzVector()
# 	mu2 = TLorentzVector()
# 	misg = TLorentzVector()
# 	v1 = TLorentzVector()
# 	v2 = TLorentzVector()
# 	Pi = math.pi
# 	mu1.SetPtEtaPhiM(pt1, eta1, phi1, 0)
# 	mu2.SetPtEtaPhiM(pt2, eta2, phi2, 0)
# 	misg.SetPtEtaPhiM(met,0,metphi,0)
# 	qq = 0.0
# 	n=8 #number of iterations in phi
# 	m=9 #number of iterations in met
# 	k=3.0 #how large can energy get w.r.t. met
# 	s11=0.0
# 	s22=0.0
# 	s12=0.0
# 	s21=0.0
# 	smin=9999999.0
# 	smax=0.0
# 	store = []
# 	for q in range(m+1):
# 		qq = q*k*met/(1.0*m)
# 		for a in range(n+1):
# 			v1.SetPtEtaPhiM(qq,0,a*Pi/(1.0*n),0)
# 			v2 = misg - v1
# 			s11 = MT(v1.Pt(),mu1.Pt(),v1.DeltaPhi(mu1))
# 			s22 = MT(v2.Pt(),mu2.Pt(),v2.DeltaPhi(mu2))
# 			s12 = MT(v1.Pt(),mu2.Pt(),v1.DeltaPhi(mu2))
# 			s21 = MT(v2.Pt(),mu1.Pt(),v2.DeltaPhi(mu1))
# 			smax = max(s11,s22,s12,s21)
# 			if (smin > smax):
# 				smin = smax
# 				store = [smin,s11,s22,s12,s21,q,a]
# 	return store

def costheta_CS(pt1,eta1,phi1,pt2,eta2,phi2,l1id):
	mu1 = TLorentzVector()
	mu2 = TLorentzVector()
	Q = TLorentzVector()
	
	if (l1id < 0):
		mu1.SetPtEtaPhiM(pt1, eta1, phi1, 0)
		mu2.SetPtEtaPhiM(pt2, eta2, phi2, 0)
	else:
		mu1.SetPtEtaPhiM(pt2, eta2, phi2, 0)
		mu2.SetPtEtaPhiM(pt1, eta1, phi1, 0)
		
	Q = mu1 + mu2

	mu1plus  = ((2.0)**(-0.5))*(mu1.E() + mu1.Pz())
	mu1minus  = ((2.0)**(-0.5))*(mu1.E() - mu1.Pz())

	mu2plus  = ((2.0)**(-0.5))*(mu2.E() + mu2.Pz())
	mu2minus = ((2.0)**(-0.5))*(mu2.E() - mu2.Pz())

	costheta = ((2.0/Q.Mag())/((Q.Mag()**2.0 + Q.Pt()**2.0)**(0.5)))*(mu1plus*mu2minus - mu1minus*mu2plus)

	return costheta
	
def sintheta_CM(pt1,eta1,phi1,ptz,etaz,phiz,mass):
	mu1 = TLorentzVector()
	zb = TLorentzVector()
	mu1.SetPtEtaPhiM(pt1,eta1,phi1,0)
	zb.SetPtEtaPhiM(ptz,etaz,phiz,mass)
	Sintheta = 2.0*(pt1/mass)*math.sin(zb.Angle(mu1.Vect()))
	#if (	zb.Angle(mu1.Vect()) < 0.0 ):
		#print zb.Angle(mu1.Vect()), "******%%%%%%%%%%*"
	return Sintheta
	
def Boosted_Angle(pt1,eta1,phi1,pt2,eta2,phi2,ptz,etaz,phiz,mass):
	mu1 = TLorentzVector()
	mu2 = TLorentzVector()
	smu1 = TLorentzVector()
	smu2 = TLorentzVector()
	zb = TLorentzVector()
	mu1.SetPtEtaPhiM(pt1,eta1,phi1,0)
	mu2.SetPtEtaPhiM(pt2,eta2,phi2,0)
	smu1.SetPtEtaPhiM(pt1,eta1,phi1,0)
	smu2.SetPtEtaPhiM(pt2,eta2,phi2,0)
	angle = mu1.Angle(mu2.Vect())
	zb.SetPtEtaPhiM(ptz,etaz,phiz,mass)
	angle_Z1 = zb.Angle(mu1.Vect())
	angle_Z2 = zb.Angle(mu2.Vect())
	mu1.Boost(-zb.Px()/zb.E(),-zb.Py()/zb.E(),-zb.Pz()/zb.E())
	mu2.Boost(-zb.Px()/zb.E(),-zb.Py()/zb.E(),-zb.Pz()/zb.E())
	angleBoost = mu1.Angle(mu2.Vect())
	angleBoost_Z1 = zb.Angle(mu1.Vect())
	angleBoost_Z2 = zb.Angle(mu2.Vect())
	Boostdiff11 = mu1.Angle(smu1.Vect())
	Boostdiff22 = mu2.Angle(smu2.Vect())
	Boostdiff12 = mu1.Angle(smu2.Vect())
	Boostdiff21 = mu2.Angle(smu1.Vect())
	#print "******&&&&******", angle, angleBoost
	return [angleBoost,angle,angleBoost_Z1,angle_Z1,angle_Z2,Boostdiff11,Boostdiff22,Boostdiff12,Boostdiff21]

def Convert(f,filesANDsamples,filesANDsamples2,indir,lumi,TeV8):

	fin = TFile.Open(f,"")
	print f, "F"*20
	H=fin.Get('all_cutflow')
	NUMGENEVENT = H.GetBinContent(1)
	BIN2 = H.GetBinContent(2)
	BIN3 = H.GetBinContent(3)

	RND_CUT = 0.66

	L1 = TLorentzVector() ####### causes function to be out of scope?!
	L2 = TLorentzVector()
	ZB = TLorentzVector()
	MET = TLorentzVector()

	aa = fin.GetListOfKeys()
	TREES = [] #array of all the input systematic trees
	for a in aa:
		#print a.GetName()
		if "finalTree" in a.GetName():
			TREES.append(a.GetName())

	if ("Data" in f) or ("nonvar" in whataboutsystematics): # no systematics for Data!!!!
		TREES = [TREES[0]]
	print TREES

	tdir = fin.Get(TREES[0]) #get list of branches to put into new tree
	#x = tdir.GetListOfKeys()
	#print tdir.GetName()
	x = tdir.GetListOfBranches()
	tagalongvariables = []
	for y in x:
		#print y.GetName()
		tagalongvariables.append(y.GetName())
	#print tagalongvariables

	addinvariables = ['XS','BR','LUM','NGE','B2','B3','RND']#,'WT']
	addinvariables += ['WT0','WT1','WT2','WT3','WT4']
	addinvariables += ['training0','training1','training2','training3','training4']#,'Wscale','Zscale']
	addinvariables += ['XS2','BR2']
	addinvariables += ['phil1met','phil2met','Thrust','DeltaPz','DeltaPhi_ZH','TransMass3','TransMass4','CScostheta','CMsintheta']
	addinvariables += ['Theta_lab','ZL1_Boost','ZL1_lab','ZL2_lab','ZRapidity','Lep2Dover3D','ZMEToverLep3D','ZMEToverLep2D','l1l2metPt','l1l2minusmetPt']
	addinvariables += ['Boost11','Boost22','Boost12','Boost21']
	#addinvariables += ['mt_max','mt_11','mt_22','mt_12','mt_21','qmet','qang']


	for var in (tagalongvariables + addinvariables + TREES): #declare tag-alongs and new variables for new tree
		var = var.replace('finalTree','sys')
		exec(var+' = array.array(\'f\',[0])')

	TN = 0
	for SYS in TREES: #loop over systematic variations
		exec(SYS+'tin=fin.Get("'+SYS+'")')
		exec('N='+SYS+'tin.GetEntries()')
		print N, SYS
		TN = TN + N

	print TN, "TOTAL IN"

	FOut = TFile.Open(f.replace(indir,outdir),"RECREATE")
  	TOut = TTree("tmvatree","tmvatree")

	for var in (tagalongvariables+addinvariables+TREES): #makes branches in new tree for each tag-along
		var = var.replace('finalTree','sys')
  		exec("TOut.Branch('"+var+"',"+var+",'"+var+"/F')")


  	for SYS in TREES: #loop over each systematic variation, including "no variation"
		exec('N='+SYS+'tin.GetEntries()')
		print N, SYS

		for n in range(N): #loop over every event
			if (n%10000 == 1) or (n+1==N):
				print f+"_"+SYS+":   "+str(n+1) +' of '+str(N) +' events evaluated.'
			exec(SYS+'tin.GetEntry(n)')

			for var in tagalongvariables:
				exec(var+'[0] = '+SYS+'tin.'+var)

			for var in TREES:
				if (var == SYS):
					var = var.replace('finalTree','sys')
					exec(var+'[0] = 1.0')
				else:
					var = var.replace('finalTree','sys')
					exec(var+'[0] = 0.0')

			if "Data" in f:
				XS[0]=1.0
				BR[0]=1.0
				XS2[0]=1.0
				BR2[0]=1.0
				LUM[0]=1.0
				NGE[0]= 1.0
				B2[0] = 1.0
				B3[0] = 1.0
				RND[0] = 1.0
				#WT[0] = 1.0
				for rr in range(5):
					exec('training'+str(rr)+'[0]=0.0')
					exec('WT'+str(rr)+'[0]=1.0')
				# Wscale[0] = 1.0
				# Zscale[0] = 1.0
			else:
				XS[0]=filesANDsamples[f.replace(indir,'')][1]
				BR[0]=filesANDsamples[f.replace(indir,'')][2]
				XS2[0]=filesANDsamples2[f.replace(indir,'')][1]
				BR2[0]=filesANDsamples2[f.replace(indir,'')][2]
				LUM[0]=lumi
				NGE[0]= NUMGENEVENT
				B2[0] = BIN2
				B3[0] = BIN3
				# Wscale[0] = scale2data(f,finstate[0],TeV8)[0]
				# Zscale[0] = scale2data(f,finstate[0],TeV8)[1]
				if SYS.replace("finalTree","sys") == "sys":
					RND[0] = 1.0*random.uniform(0,1) #only do training/testing on non-varied
					# if (RND[0] < RND_CUT):
					# 	WT[0] = 1.0/(RND_CUT) #testing-training
					# 	#training0[0]=1.0
					# else:
					# 	WT[0] = 1.0/(1.0-RND_CUT) #application
						#training0[0]=0.0
					for rr in range(5):
						if (1.0*random.uniform(0,1) < RND_CUT):
							exec('training'+str(rr)+'[0]=1.0')
							exec('WT'+str(rr)+'[0]=1.0/(RND_CUT)')
						else:
							exec('training'+str(rr)+'[0]=0.0')
							exec('WT'+str(rr)+'[0]=1.0/(1.0-RND_CUT)')
				else:
					RND[0] = 1.0
					#WT[0]=1.0

					for rr in range(5):
						exec('training'+str(rr)+'[0]=0.0')
						exec('WT'+str(rr)+'[0]=1.0')


			L1.SetPtEtaPhiM(l1pt[0],l1eta[0],l1phi[0],0)
			L2.SetPtEtaPhiM(l2pt[0],l2eta[0],l2phi[0],0)
			ZB.SetPtEtaPhiM(zpt[0],zeta[0],zphi[0],mass[0])
			MET.SetPtEtaPhiM(met[0],0,metphi[0],0)

			phil2met[0] = abs(L2.DeltaPhi(MET))
			phil1met[0] = abs(L1.DeltaPhi(MET))
			Thrust[0] = (L1-L2).Pt()
			DeltaPz[0] = abs((L1-L2).Pz())
			DeltaPhi_ZH[0] = abs(ZB.DeltaPhi(MET))
			#TransMass[0] = (L1 + L2 + MET).Mt()
			#TransMass_Eff[0] = L1.Et() + L2.Et() + MET.Et()
			#TransMass2[0]= math.sqrt(2*(L1.Pt()*L2.Pt()*(1-math.cos(abs(L1.DeltaPhi(L2)))) + L1.Pt()*MET.Pt()*(1-math.cos(abs(L1.DeltaPhi(MET)))) + L2.Pt()*MET.Pt()*(1-math.cos(abs(L2.DeltaPhi(MET)) ))))
			TransMass3[0]= math.sqrt(2*(L1+L2).Pt()*MET.Pt()*(1-math.cos(abs(MET.DeltaPhi(L1+L2)))))
			TransMass4[0]= math.sqrt( ( math.sqrt( ((L1+L2).Pt())**2 + ((L1+L2).M())**2 ) + math.sqrt( (MET.Pt())**2 + ((L1+L2).M())**2 ) )**2 - ( (L1 + L2 + MET).Px() )**2 - ( (L1 + L2 + MET).Py() )**2 )

			CScostheta[0] = costheta_CS(l1pt[0],l1eta[0],l1phi[0],l2pt[0],l2eta[0],l2phi[0],l1id[0])
			CMsintheta[0] = sintheta_CM(l1pt[0],l1eta[0],l1phi[0],zpt[0],zeta[0],zphi[0],mass[0])
		
			Boost = Boosted_Angle(l1pt[0],l1eta[0],l1phi[0],l2pt[0],l2eta[0],l2phi[0],zpt[0],zeta[0],zphi[0],mass[0])
			#Theta_Boost[0] = Boost[0] #should always be ~pi
			Theta_lab[0] = Boost[1]
			ZL1_Boost[0] = Boost[2]
			ZL1_lab[0] = Boost[3]
			ZL2_lab[0] = Boost[4]

			Boost11[0] = Boost[-4]
			Boost22[0] = Boost[-3]
			Boost12[0] = Boost[-2]
			Boost21[0] = Boost[-1]

			ZRapidity[0] = ZB.Rapidity() #0.5*math.log((ZB.E()+ZB.Pz())/(ZB.E()-ZB.Pz()))

			Lep2Dover3D[0] = abs((L1.DeltaPhi(L2))/(Theta_lab[0] + 0.0001))
			ZMEToverLep3D[0] = abs((ZB.DeltaPhi(MET))/(Theta_lab[0] + 0.0001))
			ZMEToverLep2D[0] = abs((ZB.DeltaPhi(MET))/((L1.DeltaPhi(L2)) + 0.0001))

			l1l2metPt[0] = (L1+L2+MET).Pt()
			l1l2minusmetPt[0] = (L1+L2).Pt() - (MET).Pt()

			# mt2v = MT2(l1pt[0],l1eta[0],l1phi[0],l2pt[0],l2eta[0],l2phi[0],met[0],metphi[0])
			# mt_max[0] = mt2v[0]
			# mt_11[0] = mt2v[1]
			# mt_22[0] = mt2v[2]
			# mt_12[0] = mt2v[3]
			# mt_21[0] = mt2v[4]
			# qmet[0] = mt2v[5]
			# qang[0] = mt2v[6]

			TOut.Fill()

			# if (n==1) or (n==1000):
		 #    	print TInA.x, "X"
		 #    	print TOut.varx, "VARX"

		# MVAOutput[0] = reader.EvaluateMVA('BDT')
		#TOut.Fill()
		TNO = TOut.GetEntries()
	print TNO, "total out"
	if TNO != TN:
		sys.exit("INPUT YIELD DIFFERENT FROM OUTPUT YIELD: TREE PROBLEM!!!")

	FOut.Write("", TObject.kOverwrite)
	FOut.Close()

# WHICH ENERGY?
TeV8 = True
outputdir="/afs/cern.ch/work/c/chasco/OCT19_p66_8/"

SampleList = []
SampleList1 = []
SampleList2 = []
if (TeV8): #organize this better
	lumi = 19700
	 #sample, cross section, branching ratio
	SampleList.append(['DYJetsToLL_10to50',860.5,1])
	SampleList.append(['DYJetsToLL_50toInf',3532.8,1])
	SampleList.append(['SingleTbar_s',1.76,1])
	SampleList.append(['SingleTbar_to.',30.7,1])
	SampleList.append(['SingleTbar_tW.',11.1,1])
	SampleList.append(['SingleT_s',3.79,1])
	SampleList.append(['SingleT_to.',56.4,1])
	SampleList.append(['SingleT_tW.',11.1,1])
	SampleList.append(['TTJets',225.197,0.10608049])
	SampleList.append(['W1Jets',5400.0,1])
	SampleList.append(['W2Jets',1750.0,1])
	SampleList.append(['W3Jets',519.0,1])
	SampleList.append(['W4Jets',214.0,1])
	SampleList.append(['WW',57.1097,0.104976])
	SampleList.append(['WZ',32.3,0.032715576])
	SampleList.append(['ZZ',0.355036,1])
	SampleList2.append(['ZH105',0.6750,0.100974])
	SampleList2.append(['ZH115',0.5117,0.100974])
	SampleList2.append(['ZH125',0.3943,0.100974])
	SampleList2.append(['ZH135',0.3074,0.100974])
	SampleList2.append(['ZH145',0.2424,0.100974])
	SampleList1.append(['ZH105',0.7022,0.100974])
	SampleList1.append(['ZH115',0.5358,0.100974])
	SampleList1.append(['ZH125',0.4153,0.100974])
	SampleList1.append(['ZH135',0.3259,0.100974])
	SampleList1.append(['ZH145',0.2583,0.100974])
	SampleList.append(['Data',1,1])

	Processes = []
	XSections = []
	BranchRat = []
	for s in (SampleList+SampleList1):
		Processes.append(s[0])
		XSections.append(s[1])
		BranchRat.append(s[2])

	Lumi =[]
	for x in Processes:
		if 'Data' in x:
			Lumi.append(1)
		else:
			Lumi.append(lumi)

else:
	lumi = 5035
	SampleList.append(['DYJetsToLL',3048,1])
	SampleList.append(['SingleTbar_s',1.44,1])
	SampleList.append(['SingleTbar_to.',22.65,1])
	SampleList.append(['SingleTbar_tW.',7.87,1])
	SampleList.append(['SingleT_s',3.19,1])
	SampleList.append(['SingleT_to.',41.92,1])
	SampleList.append(['SingleT_tW.',7.87,1])
	SampleList.append(['TTJets',165,1])
	SampleList.append(['WJetsToLNu',31314,1])
	SampleList.append(['WW',5.5,1])
	SampleList.append(['WZ',0.856,1])
	SampleList.append(['ZZ',0.260094,1]) #"br":[0.038647521,0.994055467]
	SampleList2.append(['ZH105',0.5449,0.100974])
	SampleList2.append(['ZH115',0.4107,0.100974])
	SampleList2.append(['ZH125',0.3158,0.100974])
	SampleList2.append(['ZH135',0.2453,0.100974])
	SampleList2.append(['ZH145',0.1930,0.100974])
	SampleList1.append(['ZH105',0.5724,0.100974])
	SampleList1.append(['ZH115',0.4345,0.100974])
	SampleList1.append(['ZH125',0.3351,0.100974])
	SampleList1.append(['ZH135',0.2616,0.100974])
	SampleList1.append(['ZH145',0.2068,0.100974])
	#SampleList.append(['ZH150',0.01171,0.100974])
	SampleList.append(['Data',1,1])

	Processes = []
	XSections = []
	BranchRat = []
	for s in (SampleList+SampleList1):
		Processes.append(s[0])
		XSections.append(s[1])
		BranchRat.append(s[2])

	Lumi =[]
	for x in Processes:
		if 'Data' in x:
			Lumi.append(1)
		else:
			Lumi.append(lumi)

#RND_CUT = 1.0/3.0
SampleList_1 = SampleList + SampleList1
SampleList_2 = SampleList + SampleList2
SampleList_1.sort()
SampleList_2.sort()
filesANDsamples = dict()
filesANDsamples2 = dict()
for samp in SampleList_1:
	if 'Data' not in samp[0]:
		filesANDsamples[samp[0].replace('.','')+'.root']=samp
for samp in SampleList_2:
	if 'Data' not in samp[0]:
		filesANDsamples2[samp[0].replace('.','')+'.root']=samp

#print filesANDsamples, lumi, filesANDsamples['ZZ.root'],  filesANDsamples['ZZ.root'][1]
#sys.exit("ok")


if len(Processes) != len(XSections):
	sys.exit('Processes not equal to XSections!!') #if arrays unequal, error
if len(Processes) != len(BranchRat):
	sys.exit('Processes not equal to BranchingRat!!')
if len(Processes) != len(Lumi):
	sys.exit('Processes not equal to Lumi!!')

initdir = '/tmp/chasco/INIT/' #make this directory by hand and put analyzer output in it
indir = initdir + 'HADD/' 
outdir = indir + 'TMVA/'

if "hadd" in hadder:
	os.system('rm -r /tmp/chasco/INIT/HADD')

	rfiles = os.listdir(initdir)
	for rr in rfiles:
		if "SingleT_t." in rr:
			os.system("mv "+initdir+rr+" "+initdir+rr.replace("_t.","_to.")) #erase name ambiguity for later stages
		if "SingleTbar_t." in rr:
			os.system("mv "+initdir+rr+" "+initdir+rr.replace("_t.","_to."))

	os.system('mkdir '+indir)#creates new hadd output directories
	os.system('mkdir '+outdir)

	# for suff in SystematicSuffixList:
	# 	os.system('mkdir '+outdir+'Trees'+suff) #creates directories for each systematic variation, including "no variation"

	files_string = str(os.listdir(initdir))
	print files_string

	for p in range(len(Processes)): #hadd together MC to normalize to 1/pb correctly in tree
		if ("Data" not in Processes[p]) and (Processes[p] in files_string):
			os.system('hadd ' + indir + Processes[p].replace('.','') + '.root ' + initdir + '*' + Processes[p].replace('.','') + '*.root')
	if ("Data" in str(Processes)) and ("Data" in files_string):
		os.system('cp '+initdir+ '*Data*.root ' + indir) #don't hadd data
else:
	os.system('mkdir '+outdir)
		
files = os.listdir(indir)

infiles=[]
#outfiles=[]

for f in files:
	
	if '.root' not in f:
		continue
	
	infiles.append(indir+f)
	#outfiles.append(outdir+f)

print "files"
print infiles
#print outfiles

# SystematicSuffixList = ["","_jerup","_jerdown","_jesup","_jesdown","_umetup","_umetdown","_lesup","_lesdown","_puup","_pudown","_btagup","_btagdown"]#,"_sherpaup","_sherpadown"]
# DataSuff = [""] #don't do systematics for data!!!

infiles.sort()

MCinfiles = []
for i in infiles:
	if "Data" not in i:
		MCinfiles.append(i)

print len(MCinfiles), len(SampleList+SampleList1), "compare "*10

#infiles = [infiles[-2]] #one file run, COMMENT THIS!!!!

for f in infiles: #loop over files
	Convert(f,filesANDsamples,filesANDsamples2,indir,lumi,TeV8)

os.system("mkdir "+outputdir)
os.system("cp "+outdir+"*.root "+outputdir)


