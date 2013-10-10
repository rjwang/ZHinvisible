#include <iostream>
#include <boost/shared_ptr.hpp>

#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuSummaryHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuPhysicsEvent.h"
#include "CMGTools/HtoZZ2l2nu/interface/METUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/ZHUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/GammaEventHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/plotter.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/interface/SmartSelectionMonitor.h"
#include "CMGTools/HtoZZ2l2nu/interface/TMVAUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/MacroUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/EventCategory.h"
#include "CMGTools/HtoZZ2l2nu/interface/EfficiencyMap.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "CMGTools/HtoZZ2l2nu/src/MuScleFitCorrector_v4_1/MuScleFitCorrector.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TMath.h"

using namespace std;

int main(int argc, char* argv[])
{
    //##############################################
    //########    GLOBAL INITIALIZATION     ########
    //##############################################

    // check arguments
    if(argc<2) {
        std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
        exit(0);
    }

    // load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    // configure the process
    const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

    bool use2011Id = runProcess.getParameter<bool>("is2011");
    bool useCHS(true);
    bool useJERsmearing(true);
    bool useJetsOnlyInTracker(false);
    bool usePUsubJetId(true);
    // invert ajets <-> jets ...
    // 2011:  ajets -> non CHS,  jets -> CHS
    // 2012:  ajets -> CHS,      jets -> non CHS
    if(use2011Id) useCHS = !useCHS;  // May17. RJ TEST
    cout << "Note: will apply " << (use2011Id ? 2011 : 2012) << " version of the id's" << endl;

    bool isMC       = runProcess.getParameter<bool>("isMC");
    int mctruthmode = runProcess.getParameter<int>("mctruthmode");

    TString url=runProcess.getParameter<std::string>("input");
    TString outFileUrl(gSystem->BaseName(url));
    outFileUrl.ReplaceAll(".root","");
    if(mctruthmode!=0) {
        outFileUrl += "_filt";
        outFileUrl += mctruthmode;
    }
    TString outdir=runProcess.getParameter<std::string>("outdir");
    TString outUrl( outdir );
    gSystem->Exec("mkdir -p " + outUrl);
    int fType(0);
    if(url.Contains("DoubleEle")) fType=EE;
    if(url.Contains("DoubleMu"))  fType=MUMU;
    if(url.Contains("MuEG"))      fType=EMU;
    if(url.Contains("SingleMu"))  fType=MUMU;
    bool isSingleMuPD(!isMC && url.Contains("SingleMu"));

    TString outTxtUrl= outUrl + "/" + outFileUrl + ".txt";
    FILE* outTxtFile = NULL;
    //if(!isMC)
    outTxtFile = fopen(outTxtUrl.Data(), "w");
    printf("TextFile URL = %s\n",outTxtUrl.Data());

    //tree info
    int evStart     = runProcess.getParameter<int>("evStart");
    int evEnd       = runProcess.getParameter<int>("evEnd");
    TString dirname = runProcess.getParameter<std::string>("dirName");

    //jet energy scale uncertainties
    TString uncFile = runProcess.getParameter<std::string>("jesUncFileName");
    gSystem->ExpandPathName(uncFile);
    JetCorrectionUncertainty jecUnc(uncFile.Data());

    //systematics
    bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
	runSystematics = true;    
	TString varNames[]= {"",
                         "_jerup","_jerdown",
                         "_jesup","_jesdown",
                         "_umetup","_umetdown",
                         "_lesup","_lesdown",
                         "_puup","_pudown",
                         "_btagup","_btagdown",
                         "_sherpaup","_sherpadown"
                        };
    size_t nvarsToInclude(1);
    if(runSystematics) {
        cout << "Systematics will be computed for this analysis" << endl;
        nvarsToInclude=sizeof(varNames)/sizeof(TString);
    }

    // Muon scale/resolution corrections
    TString fitParametersFile = "/afs/cern.ch/work/c/chasco/RJS/CMSSW_5_3_3_patch2/src/CMGTools/HtoZZ2l2nu/src/MuScleFitCorrector_v4_1/";
    if(use2011Id) {
        if(isMC) fitParametersFile += "MuScleFit_2011_MC_44X.txt";
        else     fitParametersFile += "MuScleFit_2011_DATA_44X.txt";
    } else {
        if(isMC) fitParametersFile += "MuScleFit_2012_MC_53X.txt"; // CHANGE!!!
        else {
            if( url.Contains("2012D") ) fitParametersFile += "MuScleFit_2012D_DATA_53X.txt";
            else                        fitParametersFile += "MuScleFit_2012ABC_DATA_53X.txt";
        }
    }
    MuScleFitCorrector *corrector_ = new MuScleFitCorrector(fitParametersFile);


  TTree *stmkr=new TTree("finalTree","finalTree");
  Double_t finstate = -99999.0;
  stmkr->Branch("finstate",&finstate,"finstate/D");
  Double_t Eweight = 1.0;
  stmkr->Branch("Eweight",&Eweight,"Eweight/D");
  Double_t mass = -99999.0;
  stmkr->Branch("mass",&mass,"mass/D");
  Double_t zpt = -99999.0;
  stmkr->Branch("zpt",&zpt,"zpt/D");
  Double_t zeta = -99999.0;
  stmkr->Branch("zeta",&zeta,"zeta/D");
  Double_t zphi = -99999.0;
  stmkr->Branch("zphi",&zphi,"zphi/D");
  Double_t l1pt = -99999.0;
  stmkr->Branch("l1pt",&l1pt,"l1pt/D");
  Double_t l2pt = -99999.0;
  stmkr->Branch("l2pt",&l2pt,"l2pt/D");
  Double_t l1eta = -99999.0;
  stmkr->Branch("l1eta",&l1eta,"l1eta/D");
  Double_t l2eta = -99999.0;
  stmkr->Branch("l2eta",&l2eta,"l2eta/D");
  Double_t l1phi = -99999.0;
  stmkr->Branch("l1phi",&l1phi,"l1phi/D");
  Double_t l2phi = -99999.0;
  stmkr->Branch("l2phi",&l2phi,"l2phi/D");
  Double_t met = -99999.0;
  stmkr->Branch("met",&met,"met/D");
  Double_t metphi = -99999.0;
  stmkr->Branch("metphi",&metphi,"metphi/D");
  Double_t l1Err = -99999.0;
  stmkr->Branch("l1Err",&l1Err,"l1Err/D");
  Double_t l2Err = -99999.0;
  stmkr->Branch("l2Err",&l2Err,"l2Err/D");
  Int_t l1id = -1;
  stmkr->Branch("l1id",&l1id,"l1id/I");
  Int_t l2id = -1;
  stmkr->Branch("l2id",&l2id,"l2id/I");
  Double_t REDmet = -99999.0;
  stmkr->Branch("REDmet",&REDmet,"REDmet/D");
  Int_t pZmass = 0;
  Int_t pZpt = 0;
  Int_t pBveto = 0;
  Int_t pLepVeto = 0;
  Int_t pDphijmet = 0;
  Int_t pBalance = 0;
  Int_t predMet = 0;
  Int_t pJetVeto = 0;
  Double_t REDmetmetphi = -999.0;
  stmkr->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  Int_t pMT80 = 0;
  Int_t pMT105 = 0;
  Int_t pMT120 = 0;
  Int_t pMT135 = 0;
  Int_t pET65 = 0;
  Int_t pET80 = 0;
  Double_t baldiff = -99999.0; //difference in METs
  stmkr->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr->Branch("pET65",&pET65,"pET65/I");
  stmkr->Branch("pET80",&pET80,"pET80/I");
  stmkr->Branch("baldiff",&baldiff,"baldiff/D");
  Double_t Zmetphi = -99999.0;
  stmkr->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  Double_t llphi = -99999.0;
  stmkr->Branch("llphi",&llphi,"llphi/D");
  Int_t nj15 = 0;
  stmkr->Branch("nj15",&nj15,"nj15/I");
    Int_t nj20 = 0;
  stmkr->Branch("nj20",&nj20,"nj20/I");
    Int_t nj25 = 0;
  stmkr->Branch("nj25",&nj25,"nj25/I");
    Int_t nj30 = 0;
  stmkr->Branch("nj30",&nj30,"nj30/I");
  Double_t mtzh3 = -99999.0;
  stmkr->Branch("mtzh3",&mtzh3,"mtzh3/D");
    Double_t mtzh = -99999.0;
  stmkr->Branch("mtzh",&mtzh,"mtzh/D");
    Double_t ColinSoper = -99999.0;
  stmkr->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  Double_t AxialMet = -99999.0;
  stmkr->Branch("AxialMet",&AxialMet,"AxialMet/D");
  Double_t ZREDmetphi = -99999.0;
  stmkr->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

    //Int_t nj10 = 0;
  //stmkr->Branch("nj10",&nj10,"nj10/I");


  TTree *stmkr_jesup=new TTree("finalTree_jesup","finalTree_jesup");
  stmkr_jesup->Branch("finstate",&finstate,"finstate/D");
  stmkr_jesup->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_jesup->Branch("mass",&mass,"mass/D");
  stmkr_jesup->Branch("zpt",&zpt,"zpt/D");
  stmkr_jesup->Branch("zeta",&zeta,"zeta/D");
  stmkr_jesup->Branch("zphi",&zphi,"zphi/D");
  stmkr_jesup->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_jesup->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_jesup->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_jesup->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_jesup->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_jesup->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_jesup->Branch("met",&met,"met/D");
  stmkr_jesup->Branch("metphi",&metphi,"metphi/D");
  stmkr_jesup->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_jesup->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_jesup->Branch("l1id",&l1id,"l1id/I");
  stmkr_jesup->Branch("l2id",&l2id,"l2id/I");
  stmkr_jesup->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_jesup->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_jesup->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_jesup->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_jesup->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_jesup->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_jesup->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_jesup->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_jesup->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_jesup->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_jesup->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_jesup->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_jesup->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_jesup->Branch("pET65",&pET65,"pET65/I");
  stmkr_jesup->Branch("pET80",&pET80,"pET80/I");
  stmkr_jesup->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_jesup->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_jesup->Branch("llphi",&llphi,"llphi/D");
  stmkr_jesup->Branch("nj15",&nj15,"nj15/I");
  stmkr_jesup->Branch("nj20",&nj20,"nj20/I");
  stmkr_jesup->Branch("nj25",&nj25,"nj25/I");
  stmkr_jesup->Branch("nj30",&nj30,"nj30/I");
  stmkr_jesup->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_jesup->Branch("mtzh",&mtzh,"mtzh/D");
    stmkr_jesup->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_jesup->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_jesup->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

  TTree *stmkr_jesdown=new TTree("finalTree_jesdown","finalTree_jesdown");
  stmkr_jesdown->Branch("finstate",&finstate,"finstate/D");
  stmkr_jesdown->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_jesdown->Branch("mass",&mass,"mass/D");
  stmkr_jesdown->Branch("zpt",&zpt,"zpt/D");
  stmkr_jesdown->Branch("zeta",&zeta,"zeta/D");
  stmkr_jesdown->Branch("zphi",&zphi,"zphi/D");
  stmkr_jesdown->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_jesdown->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_jesdown->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_jesdown->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_jesdown->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_jesdown->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_jesdown->Branch("met",&met,"met/D");
  stmkr_jesdown->Branch("metphi",&metphi,"metphi/D");
  stmkr_jesdown->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_jesdown->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_jesdown->Branch("l1id",&l1id,"l1id/I");
  stmkr_jesdown->Branch("l2id",&l2id,"l2id/I");
  stmkr_jesdown->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_jesdown->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_jesdown->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_jesdown->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_jesdown->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_jesdown->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_jesdown->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_jesdown->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_jesdown->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_jesdown->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_jesdown->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_jesdown->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_jesdown->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_jesdown->Branch("pET65",&pET65,"pET65/I");
  stmkr_jesdown->Branch("pET80",&pET80,"pET80/I");
  stmkr_jesdown->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_jesdown->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_jesdown->Branch("llphi",&llphi,"llphi/D");
  stmkr_jesdown->Branch("nj15",&nj15,"nj15/I");
  stmkr_jesdown->Branch("nj20",&nj20,"nj20/I");
  stmkr_jesdown->Branch("nj25",&nj25,"nj25/I");
  stmkr_jesdown->Branch("nj30",&nj30,"nj30/I");
  stmkr_jesdown->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_jesdown->Branch("mtzh",&mtzh,"mtzh/D");
    stmkr_jesdown->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_jesdown->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_jesdown->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

  TTree *stmkr_jerup=new TTree("finalTree_jerup","finalTree_jerup");
  stmkr_jerup->Branch("finstate",&finstate,"finstate/D");
  stmkr_jerup->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_jerup->Branch("mass",&mass,"mass/D");
  stmkr_jerup->Branch("zpt",&zpt,"zpt/D");
  stmkr_jerup->Branch("zeta",&zeta,"zeta/D");
  stmkr_jerup->Branch("zphi",&zphi,"zphi/D");
  stmkr_jerup->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_jerup->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_jerup->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_jerup->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_jerup->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_jerup->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_jerup->Branch("met",&met,"met/D");
  stmkr_jerup->Branch("metphi",&metphi,"metphi/D");
  stmkr_jerup->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_jerup->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_jerup->Branch("l1id",&l1id,"l1id/I");
  stmkr_jerup->Branch("l2id",&l2id,"l2id/I");
  stmkr_jerup->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_jerup->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_jerup->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_jerup->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_jerup->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_jerup->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_jerup->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_jerup->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_jerup->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_jerup->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_jerup->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_jerup->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_jerup->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_jerup->Branch("pET65",&pET65,"pET65/I");
  stmkr_jerup->Branch("pET80",&pET80,"pET80/I");
  stmkr_jerup->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_jerup->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_jerup->Branch("llphi",&llphi,"llphi/D");
  stmkr_jerup->Branch("nj15",&nj15,"nj15/I");
  stmkr_jerup->Branch("nj20",&nj20,"nj20/I");
  stmkr_jerup->Branch("nj25",&nj25,"nj25/I");
  stmkr_jerup->Branch("nj30",&nj30,"nj30/I");
  stmkr_jerup->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_jerup->Branch("mtzh",&mtzh,"mtzh/D");
  stmkr_jerup->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_jerup->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_jerup->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

  TTree *stmkr_jerdown=new TTree("finalTree_jerdown","finalTree_jerdown");
  stmkr_jerdown->Branch("finstate",&finstate,"finstate/D");
  stmkr_jerdown->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_jerdown->Branch("mass",&mass,"mass/D");
  stmkr_jerdown->Branch("zpt",&zpt,"zpt/D");
  stmkr_jerdown->Branch("zeta",&zeta,"zeta/D");
  stmkr_jerdown->Branch("zphi",&zphi,"zphi/D");
  stmkr_jerdown->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_jerdown->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_jerdown->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_jerdown->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_jerdown->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_jerdown->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_jerdown->Branch("met",&met,"met/D");
  stmkr_jerdown->Branch("metphi",&metphi,"metphi/D");
  stmkr_jerdown->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_jerdown->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_jerdown->Branch("l1id",&l1id,"l1id/I");
  stmkr_jerdown->Branch("l2id",&l2id,"l2id/I");
  stmkr_jerdown->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_jerdown->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_jerdown->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_jerdown->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_jerdown->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_jerdown->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_jerdown->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_jerdown->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_jerdown->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_jerdown->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_jerdown->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_jerdown->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_jerdown->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_jerdown->Branch("pET65",&pET65,"pET65/I");
  stmkr_jerdown->Branch("pET80",&pET80,"pET80/I");
  stmkr_jerdown->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_jerdown->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_jerdown->Branch("llphi",&llphi,"llphi/D");
  stmkr_jerdown->Branch("nj15",&nj15,"nj15/I");
  stmkr_jerdown->Branch("nj20",&nj20,"nj20/I");
  stmkr_jerdown->Branch("nj25",&nj25,"nj25/I");
  stmkr_jerdown->Branch("nj30",&nj30,"nj30/I");
  stmkr_jerdown->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_jerdown->Branch("mtzh",&mtzh,"mtzh/D");
  stmkr_jerdown->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_jerdown->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_jerdown->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

  TTree *stmkr_umetup=new TTree("finalTree_umetup","finalTree_umetup");
  stmkr_umetup->Branch("finstate",&finstate,"finstate/D");
  stmkr_umetup->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_umetup->Branch("mass",&mass,"mass/D");
  stmkr_umetup->Branch("zpt",&zpt,"zpt/D");
  stmkr_umetup->Branch("zeta",&zeta,"zeta/D");
  stmkr_umetup->Branch("zphi",&zphi,"zphi/D");
  stmkr_umetup->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_umetup->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_umetup->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_umetup->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_umetup->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_umetup->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_umetup->Branch("met",&met,"met/D");
  stmkr_umetup->Branch("metphi",&metphi,"metphi/D");
  stmkr_umetup->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_umetup->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_umetup->Branch("l1id",&l1id,"l1id/I");
  stmkr_umetup->Branch("l2id",&l2id,"l2id/I");
  stmkr_umetup->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_umetup->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_umetup->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_umetup->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_umetup->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_umetup->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_umetup->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_umetup->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_umetup->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_umetup->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_umetup->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_umetup->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_umetup->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_umetup->Branch("pET65",&pET65,"pET65/I");
  stmkr_umetup->Branch("pET80",&pET80,"pET80/I");
  stmkr_umetup->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_umetup->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_umetup->Branch("llphi",&llphi,"llphi/D");
  stmkr_umetup->Branch("nj15",&nj15,"nj15/I");
  stmkr_umetup->Branch("nj20",&nj20,"nj20/I");
  stmkr_umetup->Branch("nj25",&nj25,"nj25/I");
  stmkr_umetup->Branch("nj30",&nj30,"nj30/I");
  stmkr_umetup->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_umetup->Branch("mtzh",&mtzh,"mtzh/D");
  stmkr_umetup->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_umetup->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_umetup->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

  TTree *stmkr_umetdown=new TTree("finalTree_umetdown","finalTree_umetdown");
  stmkr_umetdown->Branch("finstate",&finstate,"finstate/D");
  stmkr_umetdown->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_umetdown->Branch("mass",&mass,"mass/D");
  stmkr_umetdown->Branch("zpt",&zpt,"zpt/D");
  stmkr_umetdown->Branch("zeta",&zeta,"zeta/D");
  stmkr_umetdown->Branch("zphi",&zphi,"zphi/D");
  stmkr_umetdown->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_umetdown->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_umetdown->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_umetdown->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_umetdown->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_umetdown->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_umetdown->Branch("met",&met,"met/D");
  stmkr_umetdown->Branch("metphi",&metphi,"metphi/D");
  stmkr_umetdown->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_umetdown->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_umetdown->Branch("l1id",&l1id,"l1id/I");
  stmkr_umetdown->Branch("l2id",&l2id,"l2id/I");
  stmkr_umetdown->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_umetdown->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_umetdown->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_umetdown->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_umetdown->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_umetdown->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_umetdown->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_umetdown->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_umetdown->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_umetdown->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_umetdown->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_umetdown->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_umetdown->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_umetdown->Branch("pET65",&pET65,"pET65/I");
  stmkr_umetdown->Branch("pET80",&pET80,"pET80/I");
  stmkr_umetdown->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_umetdown->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_umetdown->Branch("llphi",&llphi,"llphi/D");
  stmkr_umetdown->Branch("nj15",&nj15,"nj15/I");
  stmkr_umetdown->Branch("nj20",&nj20,"nj20/I");
  stmkr_umetdown->Branch("nj25",&nj25,"nj25/I");
  stmkr_umetdown->Branch("nj30",&nj30,"nj30/I");
  stmkr_umetdown->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_umetdown->Branch("mtzh",&mtzh,"mtzh/D");
  stmkr_umetdown->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_umetdown->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_umetdown->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

  TTree *stmkr_lesup=new TTree("finalTree_lesup","finalTree_lesup");
  stmkr_lesup->Branch("finstate",&finstate,"finstate/D");
  stmkr_lesup->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_lesup->Branch("mass",&mass,"mass/D");
  stmkr_lesup->Branch("zpt",&zpt,"zpt/D");
  stmkr_lesup->Branch("zeta",&zeta,"zeta/D");
  stmkr_lesup->Branch("zphi",&zphi,"zphi/D");
  stmkr_lesup->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_lesup->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_lesup->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_lesup->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_lesup->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_lesup->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_lesup->Branch("met",&met,"met/D");
  stmkr_lesup->Branch("metphi",&metphi,"metphi/D");
  stmkr_lesup->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_lesup->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_lesup->Branch("l1id",&l1id,"l1id/I");
  stmkr_lesup->Branch("l2id",&l2id,"l2id/I");
  stmkr_lesup->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_lesup->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_lesup->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_lesup->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_lesup->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_lesup->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_lesup->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_lesup->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_lesup->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_lesup->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_lesup->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_lesup->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_lesup->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_lesup->Branch("pET65",&pET65,"pET65/I");
  stmkr_lesup->Branch("pET80",&pET80,"pET80/I");
  stmkr_lesup->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_lesup->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_lesup->Branch("llphi",&llphi,"llphi/D");
  stmkr_lesup->Branch("nj15",&nj15,"nj15/I");
  stmkr_lesup->Branch("nj20",&nj20,"nj20/I");
  stmkr_lesup->Branch("nj25",&nj25,"nj25/I");
  stmkr_lesup->Branch("nj30",&nj30,"nj30/I");
  stmkr_lesup->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_lesup->Branch("mtzh",&mtzh,"mtzh/D");
  stmkr_lesup->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_lesup->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_lesup->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

  TTree *stmkr_lesdown=new TTree("finalTree_lesdown","finalTree_lesdown");
  stmkr_lesdown->Branch("finstate",&finstate,"finstate/D");
  stmkr_lesdown->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_lesdown->Branch("mass",&mass,"mass/D");
  stmkr_lesdown->Branch("zpt",&zpt,"zpt/D");
  stmkr_lesdown->Branch("zeta",&zeta,"zeta/D");
  stmkr_lesdown->Branch("zphi",&zphi,"zphi/D");
  stmkr_lesdown->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_lesdown->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_lesdown->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_lesdown->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_lesdown->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_lesdown->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_lesdown->Branch("met",&met,"met/D");
  stmkr_lesdown->Branch("metphi",&metphi,"metphi/D");
  stmkr_lesdown->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_lesdown->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_lesdown->Branch("l1id",&l1id,"l1id/I");
  stmkr_lesdown->Branch("l2id",&l2id,"l2id/I");
  stmkr_lesdown->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_lesdown->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_lesdown->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_lesdown->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_lesdown->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_lesdown->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_lesdown->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_lesdown->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_lesdown->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_lesdown->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_lesdown->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_lesdown->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_lesdown->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_lesdown->Branch("pET65",&pET65,"pET65/I");
  stmkr_lesdown->Branch("pET80",&pET80,"pET80/I");
  stmkr_lesdown->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_lesdown->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_lesdown->Branch("llphi",&llphi,"llphi/D");
  stmkr_lesdown->Branch("nj15",&nj15,"nj15/I");
  stmkr_lesdown->Branch("nj20",&nj20,"nj20/I");
  stmkr_lesdown->Branch("nj25",&nj25,"nj25/I");
  stmkr_lesdown->Branch("nj30",&nj30,"nj30/I");
  stmkr_lesdown->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_lesdown->Branch("mtzh",&mtzh,"mtzh/D");
  stmkr_lesdown->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_lesdown->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_lesdown->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

  TTree *stmkr_puup=new TTree("finalTree_puup","finalTree_puup");
  stmkr_puup->Branch("finstate",&finstate,"finstate/D");
  stmkr_puup->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_puup->Branch("mass",&mass,"mass/D");
  stmkr_puup->Branch("zpt",&zpt,"zpt/D");
  stmkr_puup->Branch("zeta",&zeta,"zeta/D");
  stmkr_puup->Branch("zphi",&zphi,"zphi/D");
  stmkr_puup->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_puup->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_puup->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_puup->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_puup->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_puup->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_puup->Branch("met",&met,"met/D");
  stmkr_puup->Branch("metphi",&metphi,"metphi/D");
  stmkr_puup->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_puup->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_puup->Branch("l1id",&l1id,"l1id/I");
  stmkr_puup->Branch("l2id",&l2id,"l2id/I");
  stmkr_puup->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_puup->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_puup->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_puup->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_puup->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_puup->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_puup->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_puup->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_puup->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_puup->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_puup->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_puup->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_puup->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_puup->Branch("pET65",&pET65,"pET65/I");
  stmkr_puup->Branch("pET80",&pET80,"pET80/I");
  stmkr_puup->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_puup->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_puup->Branch("llphi",&llphi,"llphi/D");
  stmkr_puup->Branch("nj15",&nj15,"nj15/I");
  stmkr_puup->Branch("nj20",&nj20,"nj20/I");
  stmkr_puup->Branch("nj25",&nj25,"nj25/I");
  stmkr_puup->Branch("nj30",&nj30,"nj30/I");
  stmkr_puup->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_puup->Branch("mtzh",&mtzh,"mtzh/D");
  stmkr_puup->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_puup->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_puup->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

  TTree *stmkr_pudown=new TTree("finalTree_pudown","finalTree_pudown");
  stmkr_pudown->Branch("finstate",&finstate,"finstate/D");
  stmkr_pudown->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_pudown->Branch("mass",&mass,"mass/D");
  stmkr_pudown->Branch("zpt",&zpt,"zpt/D");
  stmkr_pudown->Branch("zeta",&zeta,"zeta/D");
  stmkr_pudown->Branch("zphi",&zphi,"zphi/D");
  stmkr_pudown->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_pudown->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_pudown->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_pudown->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_pudown->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_pudown->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_pudown->Branch("met",&met,"met/D");
  stmkr_pudown->Branch("metphi",&metphi,"metphi/D");
  stmkr_pudown->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_pudown->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_pudown->Branch("l1id",&l1id,"l1id/I");
  stmkr_pudown->Branch("l2id",&l2id,"l2id/I");
  stmkr_pudown->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_pudown->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_pudown->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_pudown->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_pudown->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_pudown->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_pudown->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_pudown->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_pudown->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_pudown->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_pudown->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_pudown->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_pudown->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_pudown->Branch("pET65",&pET65,"pET65/I");
  stmkr_pudown->Branch("pET80",&pET80,"pET80/I");
  stmkr_pudown->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_pudown->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_pudown->Branch("llphi",&llphi,"llphi/D");
  stmkr_pudown->Branch("nj15",&nj15,"nj15/I");
  stmkr_pudown->Branch("nj20",&nj20,"nj20/I");
  stmkr_pudown->Branch("nj25",&nj25,"nj25/I");
  stmkr_pudown->Branch("nj30",&nj30,"nj30/I");
  stmkr_pudown->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_pudown->Branch("mtzh",&mtzh,"mtzh/D");
  stmkr_pudown->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_pudown->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_pudown->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");
/*
  TTree *stmkr_sherpaup=new TTree("finalTree_sherpaup","finalTree_sherpaup");
  stmkr_sherpaup->Branch("finstate",&finstate,"finstate/D");
  stmkr_sherpaup->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_sherpaup->Branch("mass",&mass,"mass/D");
  stmkr_sherpaup->Branch("zpt",&zpt,"zpt/D");
  stmkr_sherpaup->Branch("zeta",&zeta,"zeta/D");
  stmkr_sherpaup->Branch("zphi",&zphi,"zphi/D");
  stmkr_sherpaup->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_sherpaup->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_sherpaup->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_sherpaup->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_sherpaup->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_sherpaup->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_sherpaup->Branch("met",&met,"met/D");
  stmkr_sherpaup->Branch("metphi",&metphi,"metphi/D");
  stmkr_sherpaup->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_sherpaup->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_sherpaup->Branch("l1id",&l1id,"l1id/I");
  stmkr_sherpaup->Branch("l2id",&l2id,"l2id/I");
  stmkr_sherpaup->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_sherpaup->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_sherpaup->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_sherpaup->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_sherpaup->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_sherpaup->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_sherpaup->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_sherpaup->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_sherpaup->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_sherpaup->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_sherpaup->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_sherpaup->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_sherpaup->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_sherpaup->Branch("pET65",&pET65,"pET65/I");
  stmkr_sherpaup->Branch("pET80",&pET80,"pET80/I");
  stmkr_sherpaup->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_sherpaup->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_sherpaup->Branch("llphi",&llphi,"llphi/D");
  stmkr_sherpaup->Branch("nj15",&nj15,"nj15/I");

  TTree *stmkr_sherpadown=new TTree("finalTree_sherpadown","finalTree_sherpadown");
  stmkr_sherpadown->Branch("finstate",&finstate,"finstate/D");
  stmkr_sherpadown->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_sherpadown->Branch("mass",&mass,"mass/D");
  stmkr_sherpadown->Branch("zpt",&zpt,"zpt/D");
  stmkr_sherpadown->Branch("zeta",&zeta,"zeta/D");
  stmkr_sherpadown->Branch("zphi",&zphi,"zphi/D");
  stmkr_sherpadown->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_sherpadown->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_sherpadown->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_sherpadown->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_sherpadown->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_sherpadown->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_sherpadown->Branch("met",&met,"met/D");
  stmkr_sherpadown->Branch("metphi",&metphi,"metphi/D");
  stmkr_sherpadown->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_sherpadown->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_sherpadown->Branch("l1id",&l1id,"l1id/I");
  stmkr_sherpadown->Branch("l2id",&l2id,"l2id/I");
  stmkr_sherpadown->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_sherpadown->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_sherpadown->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_sherpadown->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_sherpadown->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_sherpadown->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_sherpadown->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_sherpadown->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_sherpadown->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_sherpadown->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_sherpadown->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_sherpadown->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_sherpadown->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_sherpadown->Branch("pET65",&pET65,"pET65/I");
  stmkr_sherpadown->Branch("pET80",&pET80,"pET80/I");
  stmkr_sherpadown->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_sherpadown->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_sherpadown->Branch("llphi",&llphi,"llphi/D");
  stmkr_sherpadown->Branch("nj15",&nj15,"nj15/I");*/

  TTree *stmkr_btagup=new TTree("finalTree_btagup","finalTree_btagup");
  stmkr_btagup->Branch("finstate",&finstate,"finstate/D");
  stmkr_btagup->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_btagup->Branch("mass",&mass,"mass/D");
  stmkr_btagup->Branch("zpt",&zpt,"zpt/D");
  stmkr_btagup->Branch("zeta",&zeta,"zeta/D");
  stmkr_btagup->Branch("zphi",&zphi,"zphi/D");
  stmkr_btagup->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_btagup->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_btagup->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_btagup->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_btagup->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_btagup->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_btagup->Branch("met",&met,"met/D");
  stmkr_btagup->Branch("metphi",&metphi,"metphi/D");
  stmkr_btagup->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_btagup->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_btagup->Branch("l1id",&l1id,"l1id/I");
  stmkr_btagup->Branch("l2id",&l2id,"l2id/I");
  stmkr_btagup->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_btagup->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_btagup->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_btagup->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_btagup->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_btagup->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_btagup->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_btagup->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_btagup->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_btagup->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_btagup->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_btagup->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_btagup->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_btagup->Branch("pET65",&pET65,"pET65/I");
  stmkr_btagup->Branch("pET80",&pET80,"pET80/I");
  stmkr_btagup->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_btagup->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_btagup->Branch("llphi",&llphi,"llphi/D");
  stmkr_btagup->Branch("nj15",&nj15,"nj15/I");
  stmkr_btagup->Branch("nj20",&nj20,"nj20/I");
  stmkr_btagup->Branch("nj25",&nj25,"nj25/I");
  stmkr_btagup->Branch("nj30",&nj30,"nj30/I");
  stmkr_btagup->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_btagup->Branch("mtzh",&mtzh,"mtzh/D");
  stmkr_btagup->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_btagup->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_btagup->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");

  TTree *stmkr_btagdown=new TTree("finalTree_btagdown","finalTree_btagdown");
  stmkr_btagdown->Branch("finstate",&finstate,"finstate/D");
  stmkr_btagdown->Branch("Eweight",&Eweight,"Eweight/D");
  stmkr_btagdown->Branch("mass",&mass,"mass/D");
  stmkr_btagdown->Branch("zpt",&zpt,"zpt/D");
  stmkr_btagdown->Branch("zeta",&zeta,"zeta/D");
  stmkr_btagdown->Branch("zphi",&zphi,"zphi/D");
  stmkr_btagdown->Branch("l1pt",&l1pt,"l1pt/D");
  stmkr_btagdown->Branch("l2pt",&l2pt,"l2pt/D");
  stmkr_btagdown->Branch("l1eta",&l1eta,"l1eta/D");
  stmkr_btagdown->Branch("l2eta",&l2eta,"l2eta/D");
  stmkr_btagdown->Branch("l1phi",&l1phi,"l1phi/D");
  stmkr_btagdown->Branch("l2phi",&l2phi,"l2phi/D");
  stmkr_btagdown->Branch("met",&met,"met/D");
  stmkr_btagdown->Branch("metphi",&metphi,"metphi/D");
  stmkr_btagdown->Branch("l1Err",&l1Err,"l1Err/D");
  stmkr_btagdown->Branch("l2Err",&l2Err,"l2Err/D");
  stmkr_btagdown->Branch("l1id",&l1id,"l1id/I");
  stmkr_btagdown->Branch("l2id",&l2id,"l2id/I");
  stmkr_btagdown->Branch("REDmet",&REDmet,"REDmet/D");
  stmkr_btagdown->Branch("pZmass",&pZmass,"pZmass/I");
  stmkr_btagdown->Branch("pZpt",&pZpt,"pZpt/I");
  stmkr_btagdown->Branch("pBveto",&pBveto,"pBveto/I");
  stmkr_btagdown->Branch("pJetVeto",&pJetVeto,"pJetVeto/I");
  stmkr_btagdown->Branch("pLepVeto",&pLepVeto,"pLepVeto/I");
  stmkr_btagdown->Branch("pDphijmet",&pDphijmet,"pDphijmet/I");
  stmkr_btagdown->Branch("pBalance",&pBalance,"pBalance/I");
  stmkr_btagdown->Branch("REDmetmetphi",&REDmetmetphi,"REDmetmetphi/D");
  stmkr_btagdown->Branch("pMT80",&pMT80,"pMT80/I");
  stmkr_btagdown->Branch("pMT105",&pMT105,"pMT105/I");
  stmkr_btagdown->Branch("pMT120",&pMT120,"pMT120/I");
  stmkr_btagdown->Branch("pMT135",&pMT135,"pMT135/I");
  stmkr_btagdown->Branch("pET65",&pET65,"pET65/I");
  stmkr_btagdown->Branch("pET80",&pET80,"pET80/I");
  stmkr_btagdown->Branch("baldiff",&baldiff,"baldiff/D");
  stmkr_btagdown->Branch("Zmetphi",&Zmetphi,"Zmetphi/D");
  stmkr_btagdown->Branch("llphi",&llphi,"llphi/D");
  stmkr_btagdown->Branch("nj15",&nj15,"nj15/I");
  stmkr_btagdown->Branch("nj20",&nj20,"nj20/I");
  stmkr_btagdown->Branch("nj25",&nj25,"nj25/I");
  stmkr_btagdown->Branch("nj30",&nj30,"nj30/I");
  stmkr_btagdown->Branch("mtzh3",&mtzh3,"mtzh3/D");
  stmkr_btagdown->Branch("mtzh",&mtzh,"mtzh/D");
  stmkr_btagdown->Branch("ColinSoper",&ColinSoper,"ColinSoper/D");
  stmkr_btagdown->Branch("AxialMet",&AxialMet,"AxialMet/D");
  stmkr_btagdown->Branch("ZREDmetphi",&ZREDmetphi,"ZREDmetphi/D");



    //##############################################
    //########    INITIATING HISTOGRAMS     ########
    //##############################################
    SmartSelectionMonitor mon;

    // pileup control
    mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtx_dy",";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F("npfjets_dy",  ";Jet multiplicity (p_{T}>30 GeV);Events",5,0,5) );
    mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) );
    mon.addHistogram( new TH1F( "rho25",";#rho(#eta<2.5);Events",50,0,25) );
    mon.addHistogram(  new TProfile("metvsrun"    ,      "Run number",     600, 190000,196000) ) ;
    mon.addHistogram(  new TProfile("metvsavginstlumi",  "Avg. inst lumi", 60,  400,1000));
    mon.addHistogram(  new TProfile("nvtxvsrun",         "Run number",     600, 190000,196000) ) ;
    mon.addHistogram(  new TProfile("nvtxvsavginstlumi", "Avg. inst lumi", 60,  400,1000));
    mon.addHistogram( new TH1F( "RunDep_Yields", ";Run;Events",4000,170000,210000) );
    mon.addHistogram( new TProfile( "RunDep_Met", ";Run;<Met>",4000,170000,210000) );

    TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;

    TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 12,0,12) );
    h->GetXaxis()->SetBinLabel(1,"Trigger");
    h->GetXaxis()->SetBinLabel(2,"#geq 2 leptons");
    h->GetXaxis()->SetBinLabel(3,"#geq 2 iso leptons");
    h->GetXaxis()->SetBinLabel(4,"|M_{ll}-M_{Z}|<15");
    h->GetXaxis()->SetBinLabel(5,"p_{T}^{ll}>50");
    h->GetXaxis()->SetBinLabel(6,"3^{rd}-lepton veto");
    h->GetXaxis()->SetBinLabel(7,"Jet veto");
    h->GetXaxis()->SetBinLabel(8,"b-veto");
    h->GetXaxis()->SetBinLabel(9,"#Delta #phi(Z,E_{T}^{miss})>2.6");
    h->GetXaxis()->SetBinLabel(10,"Reduced E_{T}^{miss}>110");
    h->GetXaxis()->SetBinLabel(11,"0.8<E_{T}^{miss}/p_{T}^{ll}<1.2");
    h->GetXaxis()->SetBinLabel(12,"M_{T}>220");

    h=(TH1F*) mon.addHistogram( new TH1F ("eventflow_gamma", ";;Events", 5,0,5) );
    //h->GetXaxis()->SetBinLabel(1,"Trigger");
    //h->GetXaxis()->SetBinLabel(2,"#geq 2 leptons");
    //h->GetXaxis()->SetBinLabel(3,"#geq 2 iso leptons");
    //h->GetXaxis()->SetBinLabel(4,"|M-M_{Z}|<15");
    //h->GetXaxis()->SetBinLabel(5,"p_{T}^{ll}>55");
    //h->GetXaxis()->SetBinLabel(6,"3^{rd}-lepton veto");
    //h->GetXaxis()->SetBinLabel(7,"Jet veto");
    //h->GetXaxis()->SetBinLabel(1,"b-veto");
    //h->GetXaxis()->SetBinLabel(2,"#Delta #phi(jet,E_{T}^{miss})>0.5");
    //h->GetXaxis()->SetBinLabel(3,"E_{T}^{miss}>70");
    //h->GetXaxis()->SetBinLabel(4,"0.4<E_{T}^{miss}/p_{T}^{ll}<1.8");
    h->GetXaxis()->SetBinLabel(1,"PreSelection");
    h->GetXaxis()->SetBinLabel(2,"#Delta #phi(Z,E_{T}^{miss})>2.6");
    h->GetXaxis()->SetBinLabel(3,"Reduced E_{T}^{miss}>110"); //RJ
    h->GetXaxis()->SetBinLabel(4,"0.8<E_{T}^{miss}/p_{T}^{ll}<1.2");
    h->GetXaxis()->SetBinLabel(5,"M_{T}>220");


    mon.addHistogram(new TH1F("ereliso",           ";RelIso;Leptons",50,0,2) );
    mon.addHistogram(new TH1F("mureliso",           ";RelIso;Leptons",50,0,2) );
    mon.addHistogram( new TH1F( "leadpt", ";p_{T}^{l} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "leadeta", ";#eta^{l};Events", 50,-2.6,2.6) );
    mon.addHistogram( new TH1F( "trailerpt", ";p_{T}^{l} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "trailereta", ";#eta^{l};Events", 50,-2.6,2.6) );

    Double_t qtaxis[100];
    for(size_t i=0; i<40; i++)  qtaxis[i]=2.5*i;       //0-97.5
    for(size_t i=0; i<20; i++)  qtaxis[40+i]=100+5*i;  //100-195
    for(size_t i=0; i<15; i++)  qtaxis[60+i]=200+10*i; //200-340
    for(size_t i=0; i<25; i++)  qtaxis[75+i]=350+25*i; //350-976

    Double_t QTaxis[54];
    for(size_t i=0; i<20; i++)  QTaxis[i]=4*i;       //0-100
    for(size_t i=0; i<10; i++)  QTaxis[20+i]=80+12*i;  //100-200
    for(size_t i=0; i<10; i++)  QTaxis[30+i]=200+10*i; //200-300
    for(size_t i=0; i<14; i++)  QTaxis[40+i]=300+50*i; //300-900

    //mon.addHistogram( new TH1D( "qt_rebin",        ";p_{T}^{#gamma} [GeV];Events",53,QTaxis));
    mon.addHistogram( new TH1D( "qt_rebin",        ";p_{T}^{#gamma} [GeV];Events",99,qtaxis));
    mon.addHistogram( new TH1F( "qt",        ";p_{T}^{#gamma} [GeV];Events",1000,0,1000));
    mon.addHistogram( new TH1F( "qmass", ";M^{ll};Events", 15,76,106) );
    mon.addHistogram( new TH1F( "qCosl1Z_CS", ";cos#theta*;Events", 80,-1.,1) );
    mon.addHistogram( new TH1F( "qdphi2lep", ";#Delta#phi_{ll} [rad];Events", 50,0,TMath::Pi()) );

    //mon.addHistogram( new TH2F("nvtxVSqtrebin",";Vertices;p_{T}^{#gamma} [GeV];Events",50,0,50,99,QTaxis));
    //mon.addHistogram( new TH2F("nvtxVSqt",";Vertices;p_{T}^{#gamma} [GeV];Events",50,0,50,1000,0,1000));

    mon.addHistogram( new TH1F( "zpt", ";p_{T}^{ll} [GeV];Events", 50,0,500) );
    //double zpt_bins[] = {0., 20., 40., 60., 80., 100., 150., 200., 300., 400., 600., 800.};
    double zpt_bins[] = {0., 20., 40., 60., 80., 100., 200., 400., 800.};
    const int n_zpt_bins = sizeof(zpt_bins)/sizeof(double) - 1;
    mon.addHistogram( new TH1F( "zpt_rebin", ";p_{T}^{ll} [GeV];Events", n_zpt_bins, zpt_bins) );
    mon.addHistogram( new TH1F( "zpt_dy", ";p_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "zpt_noJet", ";p_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "zptGen_noJet", ";p_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "zptNM1", ";p_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "zeta", ";#eta^{ll};Events", 50,-10,10) );
    mon.addHistogram( new TH1F( "zmass", ";M^{ll} [GeV];Events", 100,40,250) );
    mon.addHistogram( new TH1F( "zmassNM1", ";M^{ll} [GeV];Events", 100,40,250) );

    // 3rd lepton control
    mon.addHistogram( new TH1F( "thirdleptonpt", ";p_{T}^{l} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "thirdleptoneta", ";#eta^{l};Events", 50,-2.6,2.6) );
    h = (TH1F*) mon.addHistogram( new TH1F( "nleptons", ";Lepton multiplicity;Events", 3,2,4) );
    for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin+1);
        h->GetXaxis()->SetBinLabel(ibin,label);
    }

    mon.addHistogram( new TH1F( "nleptons_final", ";Lepton multiplicity;Events", 3,2,4) );
    mon.addHistogram( new TH1F( "nleptonsNM1", ";Lepton multiplicity;Events", 3,2,4) );
    mon.addHistogram( new TH1F( "mt3" , ";M_{T}^{3rd lepton} [GeV];Events",20,0,200) );

    // background control
    mon.addHistogram( new TH1F( "Ctrl_WZ_RedMet",";RedMet [GeV];Events",100,0,300) );
    mon.addHistogram( new TH1F( "Ctrl_WZ_zmass","Mass [GeV];Ctrl_WZ;Events",100,40,170) );
    mon.addHistogram( new TH1F( "Ctrl_WZ_Mt","Trasv. Mass [GeV];Ctrl_WZ;Events",100,0.,130) );

    //RJ for preApproval request
    mon.addHistogram( new TH1F( "WZ_contrl_Met_final",";E_{T}^{miss} [GeV];Events",100,0,500) );
    mon.addHistogram( new TH1F( "WZ_contrl_Mt_final",";M_{T}(Z,H) [GeV];Events",100,0.,1000) );
    mon.addHistogram( new TH1F( "WZ_contrl_zpt", ";p_{T}^{ll} [GeV];Events", 100,0,500) );
    //mon.addHistogram( new TH1F( "WZ_contrl_redmet", ";Reduced E_{T}^{miss} [GeV];Events", 100,0,500) );
    mon.addHistogram( new TH1F( "WZ_contrl_met", ";E_{T}^{miss} [GeV];Events", 100,0,500) );
    mon.addHistogram( new TH1F( "WZ_contrl_met_extra_e", ";E_{T}^{miss} [GeV];Events", 100,0,500) );
    mon.addHistogram( new TH1F( "WZ_contrl_met_extra_mu", ";E_{T}^{miss} [GeV];Events", 100,0,500) );
    //mon.addHistogram( new TH1F( "WZ_contrl_balance", ";E_{T}^{miss}/q_{T};Events", 25,0,5) );
    //mon.addHistogram( new TH1F( "WZ_contrl_dphiZMET",";#Delta#phi(Z,E_{T}^{miss});Events", 50,0,TMath::Pi()) );
    //mon.addHistogram( new TH1F( "WZ_contrl_dphiextrLepMET",";#Delta#phi(3rd Lepton,E_{T}^{miss});Events", 50,0,TMath::Pi()) );
    //mon.addHistogram( new TH1F( "WZ_contrl_dphiZextrLep",";#Delta#phi(Z,3rd Lepton);Events", 50,0,TMath::Pi()) );
    //mon.addHistogram( new TH1F( "WZ_contrl_dphiZLepMET",";#Delta#phi(Z,3rd Lepton+E_{T}^{miss});Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "WZ_contrl_dilepmass",";M_{ll} [GeV];Events",100,50,250) );
    mon.addHistogram( new TH1F( "WZ_contrl_trilepmass",";M_{lll} [GeV];Events",100,50,250) );
    mon.addHistogram( new TH1F( "WZ_contrl_lepmass_pair1",";M_{ll} [GeV];Events",100,0,200) );
    mon.addHistogram( new TH1F( "WZ_contrl_lepmass_pair2",";M_{ll} [GeV];Events",100,0,200) );
    mon.addHistogram( new TH1F( "WZ_contrl_lepmass_pair3",";M_{ll} [GeV];Events",100,0,200) );
    //mon.addHistogram( new TH1D( "WZ_contrl_qt_rebin",";p_{T}^{#gamma} [GeV];Events",99,qtaxis));
    mon.addHistogram( new TH1F( "WZ_contrl_newMt","Trasv. Mass [GeV];M_{T}(E_{T}^{miss},3rd Lepton) [GeV];Events",100,0.,130) );
    mon.addHistogram( new TH1F( "WZ_contrl_pt_lep1", ";Z daughter Leading Lepton p_{T} [GeV];Events", 70,0,140) );
    mon.addHistogram( new TH1F( "WZ_contrl_pt_lep2", ";Z daughter 2nd Leading Lepton p_{T} [GeV];Events", 70,0,140) );
    mon.addHistogram( new TH1F( "WZ_contrl_pt_lep3", ";W daughter Lepton p_{T} [GeV];Events", 70,0,140) );
    mon.addHistogram( new TH1F( "WZ_contrl_pt_3rdlepton_e", ";W daughter Lepton e p_{T} [GeV];Events", 70,0,140) );
    mon.addHistogram( new TH1F( "WZ_contrl_pt_3rdlepton_mu", ";W daughter Lepton #mu p_{T} [GeV];Events", 70,0,140) );



    mon.addHistogram( new TH1F( "Ctrl_WZ_Mt1","Trasv. Mass [GeV];M_{T}(E_{T}^{miss},3rd Lepton) [GeV];Events",100,0.,130) );
    mon.addHistogram( new TH1F( "Ctrl_WZ_Mt2","Trasv. Mass [GeV];Ctrl_WZ;Events",100,0.,130) );
    mon.addHistogram( new TH1F( "Ctrl_WZ_Mt3","Trasv. Mass [GeV];Ctrl_WZ;Events",100,0.,130) );
    mon.addHistogram( new TH1F( "Ctrl_WZ_Mt4","Trasv. Mass [GeV];Ctrl_WZ;Events",100,0.,130) );
    mon.addHistogram( new TH1F( "Ctrl_T_zmass",";Mass [GeV];Events",50,40,300) );
    mon.addHistogram( new TH1F( "Ctrl_Tside_RedMet",";RedMet [GeV];Events",50,0,300) );
    mon.addHistogram( new TH1F( "Ctrl_WW_RedMet",";RedMet [GeV];Events",50,0,300) );
    mon.addHistogram( new TH1F( "Ctrl_WW_RedMet1",";RedMet [GeV];Events",50,0,300) );
    mon.addHistogram( new TH1F( "Ctrl_WW_RedMet2",";RedMet [GeV];Events",50,0,300) );
    mon.addHistogram( new TH1F( "Ctrl_WW_RedMet3",";RedMet [GeV];Events",50,0,300) );
    mon.addHistogram( new TH1F( "Ctrl_WW_RedMet4",";RedMet [GeV];Events",50,0,300) );
    mon.addHistogram( new TH1F( "Ctrl_WW_RedMet5",";RedMet [GeV];Events",50,0,300) );
    mon.addHistogram( new TH1F( "Ctrl_WW_zmass",";M^{ll} [GeV];Events",45,30,250) );


    //jet control
    mon.addHistogram( new TH1F("pfjetpt"       , ";p_{T} [GeV];Events",50,0,250) );
    mon.addHistogram( new TH1F("pfjeteta"       , ";|#eta|;Events",25,0,5) );
    mon.addHistogram( new TH2F("npfjetsvspu",          ";Pileup interactions;Jet multiplicity (p_{T}>30 GeV);Events",50,0,50,5,0,5) );
    h=(TH1F *)mon.addHistogram( new TH1F("npfjets",  ";Jet multiplicity (p_{T}>30 GeV);Events",5,0,5) );
    for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h->GetXaxis()->SetBinLabel(ibin,label);
    }
    mon.addHistogram( new TH1F("npfjets_pres",  ";Jet multiplicity after preselection (p_{T}>30 GeV);Events",5,0,5) );
    mon.addHistogram( new TH1F("npfjetsNM1",    ";Jet multiplicity (p_{T}>30 GeV);Events",5,0,5) );
    mon.addHistogram( new TH1F("pfjetbeta"    , ";#beta;Events",50,0,1) );
    mon.addHistogram( new TH1F("pfjetmva"     , ";MVA;Events",50,-1,1) );
    mon.addHistogram( new TH1F("bpfjetstags",     ";b tags;Events",50,0,1) );
    mon.addHistogram( new TH1F("otherpfjetstags", ";udscg tags;Events",50,0,1) );
    mon.addHistogram( new TH1F("pfjetstags",         ";Discriminator;Events",50,0,1) );
    mon.addHistogram( new TH1F("npfjetsbtags",    ";b-tag multiplicity;Events",5,0,5) );
    mon.addHistogram( new TH1F("npfjetsbtagsNM1",    ";b-tag multiplicity;Events",5,0,5) );


    //MET control
    mon.addHistogram( new TH1F( "mindphilmet", ";min(#Delta#phi(lepton,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1F( "mindphilmet_aftPhijPt", ";min(#Delta#phi(lepton,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1F( "mindphilmet_aftPhij", ";min(#Delta#phi(lepton,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1F( "mindphilmet_aftMet", ";min(#Delta#phi(lepton,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1F( "maxdphilmet", ";max(#Delta#phi(lepton,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1F( "mindphijmet_0", ";min #Delta#phi(jet,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1F( "mindphijmet_25", ";min #Delta#phi(jet,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1F( "mindphijmet_50", ";min #Delta#phi(jet,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1F( "mindphijmet", ";min #Delta#phi(jet,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1F( "mindphijmet_aftPhil", ";min #Delta#phi(jet,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1F( "mindphijmetNM1", ";min #Delta#phi(jet,E_{T}^{miss});Events",40,0,4) );
    mon.addHistogram( new TH1D( "balancedif",    ";|E_{T}^{miss}-q_{T}|/q_{T};Events", 5,0,1.0) );
    mon.addHistogram( new TH1D( "balance",    ";E_{T}^{miss}/q_{T};Events", 25,0,5) );
    mon.addHistogram( new TH2F( "balancevsmt",          ";E_{T}^{miss}/q_{T};M_{T};Events", 4,0.8,1.2,8,20,1000) );
    mon.addHistogram( new TH1D( "balanceNM1", ";E_{T}^{miss}/q_{T};Events", 20,0,2) );
    mon.addHistogram( new TH2D( "met_mindphilmet"  , ";E_{T}^{miss};min(#Delta#phi(lepton,E_{T}^{miss});Events", 50,0,250,40,0,4) );
    mon.addHistogram( new TH2D( "metoverlpt_mindphilmet"  , ";E_{T}^{miss}/p_{T}^{lepton};min(#Delta#phi(lepton,E_{T}^{miss});Events", 50,0,2,40,0,4) );
    mon.addHistogram( new TH1F( "met_metSB",        ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_redMetSB",     ";Reduced E_{T}^{miss} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "met_met",          ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_metL",          ";Axial-E_{T}^{miss} [GeV];Events", 50,-50,150) );
    mon.addHistogram( new TH1F( "met_metCompL",          ";Long. PF E_{T}^{miss} [GeV];Events", 75,-250,500) );
    mon.addHistogram( new TH1F( "met_metCompT",          ";Trans. PF E_{T}^{miss} [GeV];Events", 75,-250,500) );
    mon.addHistogram( new TH1F( "met_metNM1",          ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_redMetNM1",       ";Reduced E_{T}^{miss} [GeV];Events", 50,0,500) );

    //RJ
    mon.addHistogram( new TH1F( "met_met_ROC",          ";E_{T}^{miss} [GeV];Events", 100,0,500) ); //RJ
    mon.addHistogram( new TH1F( "met_redMet_ROC",       ";Reduced E_{T}^{miss} [GeV];Events", 100,0,500) );
    mon.addHistogram( new TH1F( "met_met_ROC2",          ";E_{T}^{miss} [GeV];Events", 100,0,500) );
    mon.addHistogram( new TH1F( "met_redMet_ROC2",       ";Reduced E_{T}^{miss} [GeV];Events", 100,0,500) );

    mon.addHistogram( new TH2F( "metvsmt",          ";E_{T}^{miss};M_{T};Events", 100,0,500,200,0,1000) );
    mon.addHistogram( new TH2F( "typeImetvstypeImt",          ";type I E_{T}^{miss};type I M_{T};Events", 100,0,500,200,0,1000) );
    mon.addHistogram( new TH1F( "met_metbeforedphilmet",          ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_met_blind",    ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_mvamet",       ";MVA E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_mvamet_blind",       ";MVA E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_typeImet",     ";Type I E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_typeImet_blind",     ";Type I E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_fulltypeImet", ";Type I E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_fulltypeImetbeforedphilmet", ";Type I E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_fulltypeImet_blind", ";Type I E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_redMet",       ";Reduced E_{T}^{miss} [GeV];Events", 50,0,500) );
    //double redMet_bins[] = {0., 20., 40., 60., 80., 100., 200., 400., 800.};
    double redMet_bins[] = {0., 20., 40., 60., 80., 100., 150., 400., 800.};
    const int n_redMet_bins = sizeof(redMet_bins)/sizeof(double) - 1;
    mon.addHistogram( new TH1F( "met_redMet_rebin",       ";Reduced E_{T}^{miss} [GeV];Events", n_redMet_bins,redMet_bins) );
    mon.addHistogram( new TH1F( "met_redMet_phi_dy"  , ";#phi(Reduced E_{T}^{miss});Events", 60,-TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "met_met_phi_dy"     , ";#phi(E_{T}^{miss});Events", 60,-TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "met_redMet_blind",       ";Reduced E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_redMetL",      ";Long. Reduced E_{T}^{miss} [GeV];Events", 75,-250,500) );
    mon.addHistogram( new TH1F( "met_redMetT",      ";Trans. Reduced E_{T}^{miss} [GeV];Events", 75,-250,500) );
    mon.addHistogram( new TH1F( "met_clusteredMetL",	";Long. clustered E_{T}^{miss} [GeV];Events", 75,-250,500) );
    mon.addHistogram( new TH1F( "met_clusteredMetT",	";Trans. clustered E_{T}^{miss} [GeV];Events", 75,-250,500) );
    mon.addHistogram( new TH1F( "RecoilL", ";U_{L} [GeV];Events", 75,-250,500) );
    mon.addHistogram( new TH1F( "RecoilT", ";U_{T} [GeV];Events", 75,-250,500) );

    mon.addHistogram( new TH1F( "met_met3leptons",  ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    //mon.addHistogram( new TH2F( "met_met_vspu",       ";Pileup events; E_{T}^{miss};Events", 50,0,50,50,0,250) );
    //mon.addHistogram( new TH2F( "met_mvamet_vspu",    ";Pileup events; MVA E_{T}^{miss};Events", 50,0,50,50,0,250) );
    //mon.addHistogram( new TH2F( "met_redMet_vspu",    ";Pileup events; Reduced E_{T}^{miss};Events", 50,0,50,50,0,250) );
    //mon.addHistogram( new TH2F( "met_typeImet_vspu",  ";Pileup events; Type I E_{T}^{miss};Events", 50,0,50,50,0,250) );
    //mon.addHistogram( new TH2F( "met_redMetCHS_vspu", ";Pileup events; Type I E_{T}^{miss};Events", 50,0,50,50,0,250) );
    mon.addHistogram( new TH1F( "mt"  , ";M_{T}(Z,H) [GeV];Events", 85,150,1000) );
    mon.addHistogram( new TH1F( "mtNM1"  , ";M_{T}(Z,H) [GeV];Events", 80,200.,1000) );
    mon.addHistogram( new TH1F( "mt_blind"  , ";M_{T}(Z,H) [GeV];Events", 100,0,1000) );
    mon.addHistogram( new TH1F( "typeImt"  , ";M_{T}(Z,H) [GeV];Events", 100,0,1000) );

    // Final distributions
    mon.addHistogram( new TH1F( "mt_final"  , ";M_{T}(Z,H) [GeV];Events", 8,200,1000) );
    mon.addHistogram( new TH1F( "new_mt_final"  , ";M_{T}(Z,H) [GeV];Events", 10,0,1000) );
    mon.addHistogram( new TH1F( "zpt_final", ";p_{T}^{ll} [GeV];Events", 25,0,500) );
    mon.addHistogram( new TH1F( "zpt_rebin_final", ";p_{T}^{ll} [GeV];Events", n_zpt_bins, zpt_bins) );
    mon.addHistogram( new TH1F( "met_met_final"  , ";E_{T}^{miss} [GeV];Events", 25,0,250) );
    mon.addHistogram( new TH1F( "met_redMet_final"  , ";Reduced E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_redMetL_final"  , ";Longitudinal Reduced E_{T}^{miss} [GeV];Events", 25,-100,400) );
    mon.addHistogram( new TH1F( "met_redMet_rebin_final", ";Reduced E_{T}^{miss} [GeV];Events", n_redMet_bins,redMet_bins) );

    //RJ
    mon.addHistogram( new TH1F( "dphi_Zll_met_met",";#Delta#phi(Z,E_{T}^{miss});Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "dphi_Zll_met_met_final",";#Delta#phi(Z,E_{T}^{miss});Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "dphi_Zll_met_metNM1",";#Delta#phi(Z,E_{T}^{miss});Events", 50,2,TMath::Pi()) );
    mon.addHistogram( new TH1F( "dphi_Zll_met_redMet_final",";#Delta#phi(Z,Reduced E_{T}^{miss});Events", 50,0,TMath::Pi()) );

    //Collins Soper Frame
    mon.addHistogram( new TH1F( "Cos_lep1Z_CS",";cos(#theta*(l,Z));Events", 80,-1.,1.) );
    mon.addHistogram( new TH1F( "Cos_lep1Z_CS_final",";cos(#theta*(l,Z));Events", 80,-1.,1.) );

    //RJ. Looser preSelection (PRS)
    mon.addHistogram( new TH1F( "dphi_Zll_met_metPRS",";#Delta#phi(Z,E_{T}^{miss});Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "mtPRS",";M_{T} [GeV];Events", 100,0,1000) );
    mon.addHistogram( new TH1F( "zmassPRS",";M^{ll} [GeV];Events", 60,60,120) );
    mon.addHistogram( new TH1F( "zptPRS",";p_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "nleptonsPRS",";Lepton multiplicity;Events", 3,2,4) );
    mon.addHistogram( new TH1F( "npfjetsPRS",";Jet multiplicity (p_{T}>30 GeV);Events",5,0,5) );
    mon.addHistogram( new TH1F( "npfjetsbtagsPRS",";b-tag multiplicity;Events",5,0,5) );
    mon.addHistogram( new TH1D( "balancePRS",";E_{T}^{miss}/q_{T};Events", 25,0,5) );
    mon.addHistogram( new TH1F( "met_metPRS",";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "met_redMetPRS",";Reduced E_{T}^{miss} [GeV];Events", 50,0,500) );


    //NLO EWK Correction study
    mon.addHistogram( new TH1F( "v_pt",";p_{T} [GeV];Events", 500,0,500) );
    mon.addHistogram( new TH1F( "v_pt_ewkcorr",";p_{T} EWK Corr [GeV];Events", 500,0,500) );

    //##############################################
    //######## STUFF FOR CUTS OPTIMIZATION  ########
    //##############################################

/* 
      //run full optimization
       std::vector<double> optim_Cuts1_met;
       std::vector<double> optim_Cuts1_balance;
       std::vector<double> optim_Cuts1_dphi;
       std::vector<double> optim_Cuts1_zmass;
       for(double met=90;met<140;met+=10) {
         for(double balance=0.05;balance<=0.4;balance+=0.05) {
           for(double dphi=2.5;dphi<=3.0;dphi+=0.1) {
             for(double zm=10;zm<=20;zm+=2.5)
    	 	{
    	   		optim_Cuts1_met	.push_back(met);
    	  		optim_Cuts1_balance	.push_back(balance);
    	   		optim_Cuts1_dphi	.push_back(dphi);
    	   		optim_Cuts1_zmass	.push_back(zm);
    	  	}
           }
         }
       }
*/  


    //optimization
    std::vector<double> optim_Cuts1_met;
    std::vector<double> optim_Cuts1_balance;
    std::vector<double> optim_Cuts1_dphi;
    std::vector<double> optim_Cuts1_zmass;

    //fast run fixed limit
    for(double met=90; met<=110; met+=10) {
        for(double balance=0.2; balance<=0.2; balance+=0.05) {
            for(double dphi=2.6; dphi<=2.6; dphi+=0.1) {
                for(double zm=15; zm<=17.5; zm+=2.5) {
		//for(double zm=0.; zm<=3.14; zm+=0.1) {
                    optim_Cuts1_met	.push_back(met);
                    optim_Cuts1_balance	.push_back(balance);
                    optim_Cuts1_dphi	.push_back(dphi);
                    optim_Cuts1_zmass	.push_back(zm);
                }
            }
        }
    }

//Fast check hists


    //make it as a TProfile so hadd does not change the value
    TProfile* Hoptim_cuts1_met    =  (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_met",    ";cut index;met",     optim_Cuts1_met.size(),0,optim_Cuts1_met.size()) ) ;
    TProfile* Hoptim_cuts1_balance=  (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_balance",";cut index;dphi",    optim_Cuts1_met.size(),0,optim_Cuts1_met.size()) ) ;
    TProfile* Hoptim_cuts1_dphi    =  (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_dphi",    ";cut index;dphi",  optim_Cuts1_met.size(),0,optim_Cuts1_met.size()) ) ;
    TProfile* Hoptim_cuts1_zmass  =  (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_zm",     ";cut index;zmass",   optim_Cuts1_met.size(),0,optim_Cuts1_met.size()) ) ;
    for(unsigned int index=0; index<optim_Cuts1_met.size(); index++) {
        Hoptim_cuts1_met		->Fill(index, optim_Cuts1_met[index]);
        Hoptim_cuts1_balance		->Fill(index, optim_Cuts1_balance[index]);
        Hoptim_cuts1_dphi		->Fill(index, optim_Cuts1_dphi[index]);
        Hoptim_cuts1_zmass		->Fill(index, optim_Cuts1_zmass[index]);
    }

    TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;



    double BalanceXaxis[4], DphiLLXaxis[4], CosThetaXaxis[4];
    BalanceXaxis[0] = 0.8;
    BalanceXaxis[1] = 0.9;
    BalanceXaxis[2] = 1.1;
    BalanceXaxis[3] = 1.2;
    DphiLLXaxis[0] = 0.;
    DphiLLXaxis[1] = 0.6;
    DphiLLXaxis[2] = 1.4;
    DphiLLXaxis[3] = TMath::Pi();
    CosThetaXaxis[0] = -1.;
    CosThetaXaxis[1] = -0.5;
    CosThetaXaxis[2] = 0.;
    CosThetaXaxis[3] = 1.;

    for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
        Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);

        mon.addHistogram( new TH1F( "zptvar"+varNames[ivar], ";p_{T}^{ll} [GeV];Events", 100,0,1000) );

        Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);

        // for ZZ analysis (use Zpt or MET shape)
        //mon.addHistogram( new TH2F (TString("redMet_shapes")+varNames[ivar],";cut index;red-MET [GeV];#events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), 160,0,800) );
        mon.addHistogram( new TH2F (TString("redMet_shapes")+varNames[ivar],";cut index;red-MET [GeV];#Events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),80,0,800) );
        mon.addHistogram( new TH2F (TString("redMet_rebin_shapes")+varNames[ivar],";cut index;red-MET [GeV];#Events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), n_redMet_bins,redMet_bins) );
        //RJ
        mon.addHistogram( new TH2F (TString("mt_shapes")+varNames[ivar],";cut index;M_{T}(Z,H) [GeV];#Events (/100GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),12,0,1200) );
        mon.addHistogram( new TH2F (TString("new_mt_shapes")+varNames[ivar],";cut index;M_{T}(Z,H) [GeV];#Events (/100GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),12,0,1200) );

        //2D shapes
        TH2F *hh=(TH2F *) mon.addHistogram( new TH2F (TString("balance_mt_shapes")+varNames[ivar],";cut index;M_{T}(Z,H) [GeV];#Events (/100GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),36,0,3600) );
        for(int j=1; j<=36; j++) hh->GetYaxis()->SetBinLabel(j,"");
        hh->GetYaxis()->SetBinLabel(1 ,"0.8<E_{T}^{miss}/q_{T}<0.9");
        hh->GetYaxis()->SetBinLabel(13,"0.9<E_{T}^{miss}/q_{T}<1.1");
        hh->GetYaxis()->SetBinLabel(25,"1.1<E_{T}^{miss}/q_{T}<1.2");
        hh->GetYaxis()->SetBinLabel(5 ,"0<M_{T}(Z,H)<1200");
        hh->GetYaxis()->SetBinLabel(17,"0<M_{T}(Z,H)<1200");
        hh->GetYaxis()->SetBinLabel(29,"0<M_{T}(Z,H)<1200");

        hh=(TH2F *) mon.addHistogram( new TH2F (TString("dphiLL_mt_shapes")+varNames[ivar],";cut index;M_{T}(Z,H) [GeV];#Events (/100GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),36,0,3600) );
        for(int j=1; j<=36; j++) hh->GetYaxis()->SetBinLabel(j,"");
        hh->GetYaxis()->SetBinLabel(1 ,"0<#Delta#phi_{ll}<0.6");
        hh->GetYaxis()->SetBinLabel(13,"0.6<#Delta#phi_{ll}<1.4");
        hh->GetYaxis()->SetBinLabel(25,"1.4<#Delta#phi_{ll}<#pi");
        hh->GetYaxis()->SetBinLabel(5 ,"0<M_{T}(Z,H)<1200");
        hh->GetYaxis()->SetBinLabel(17,"0<M_{T}(Z,H)<1200");
        hh->GetYaxis()->SetBinLabel(29,"0<M_{T}(Z,H)<1200");

        hh=(TH2F *) mon.addHistogram( new TH2F (TString("coslZ_mt_shapes")+varNames[ivar],";cut index;M_{T}(Z,H) [GeV];#Events (/100GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),36,0,3600) );
        for(int j=1; j<=36; j++) hh->GetYaxis()->SetBinLabel(j,"");
        hh->GetYaxis()->SetBinLabel(1 ,"-1.0<cos#theta*<-0.5");
        hh->GetYaxis()->SetBinLabel(13,"-0.5<cos#theta*<0.0");
        hh->GetYaxis()->SetBinLabel(25,"0.0<cos#theta*<1.0");
        hh->GetYaxis()->SetBinLabel(5 ,"0<M_{T}(Z,H)<1200");
        hh->GetYaxis()->SetBinLabel(17,"0<M_{T}(Z,H)<1200");
        hh->GetYaxis()->SetBinLabel(29,"0<M_{T}(Z,H)<1200");

        /* control plots */
        mon.addHistogram( new TH2F (TString("balance1Dshape")+varNames[ivar],";cut index;E_{T}^{miss}/q_{T};Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),3,BalanceXaxis) );
        mon.addHistogram( new TH2F (TString("dphiLL1Dshape")+varNames[ivar],";cut index;#Delta#phi_{ll} [rad];Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),3,DphiLLXaxis) );
        mon.addHistogram( new TH2F (TString("coslZ1Dshape")+varNames[ivar],";cut index;cos#theta* [rad];Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),3,CosThetaXaxis) );

	//add 2D shape plots
    	//mon.addHistogram( new TH2F ("dphiLLvsMT",";M_{T} [GeV];#Delta#phi_{ll} [rad];Events",12,0,1200,3,DphiLLXaxis) );



        mon.addHistogram( new TH2F (TString("dphi_shapes")+varNames[ivar],";cut index; #Delta#phi(Z,E_{T}^{miss}) [rad];#events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),5,2.5,TMath::Pi()) );
        mon.addHistogram( new TH2F (TString("met_shapes")+varNames[ivar],";cut index;MET [GeV];#events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), 160,0,800) );
        mon.addHistogram( new TH2F (TString("met_rebin_shapes")+varNames[ivar],";cut index;MET [GeV];#events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), n_redMet_bins,redMet_bins) );

        mon.addHistogram( new TH2F (TString("zpt_shapes")+varNames[ivar],";cut index;Z p_{T} [GeV];#events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), 160,0,800) );
        mon.addHistogram( new TH2F (TString("zpt_rebin_shapes")+varNames[ivar],";cut index;Z p_{T} [GeV];#events (/5GeV)",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), n_zpt_bins, zpt_bins) );


        //Non-Resonant Bkgs

        TH2F *h2=(TH2F *) mon.addHistogram( new TH2F ("redMet_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        h2=(TH2F *) mon.addHistogram( new TH2F ("redMet_rebin_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        //RJ
        h2=(TH2F *) mon.addHistogram( new TH2F ("mt_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        h2=(TH2F *) mon.addHistogram( new TH2F ("new_mt_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        //2D shapes
        h2=(TH2F *) mon.addHistogram( new TH2F ("balance_mt_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        h2=(TH2F *) mon.addHistogram( new TH2F ("dphiLL_mt_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        h2=(TH2F *) mon.addHistogram( new TH2F ("coslZ_mt_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        /*control plots */
        h2=(TH2F *) mon.addHistogram( new TH2F ("balance1Dshape_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        h2=(TH2F *) mon.addHistogram( new TH2F ("dphiLL1Dshape_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        h2=(TH2F *) mon.addHistogram( new TH2F ("coslZ1Dshape_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");




        h2=(TH2F *) mon.addHistogram( new TH2F ("dphi_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        h2=(TH2F *) mon.addHistogram( new TH2F ("met_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        h2=(TH2F *) mon.addHistogram( new TH2F ("met_rebin_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        h2=(TH2F *) mon.addHistogram( new TH2F ("zpt_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");

        h2=(TH2F *) mon.addHistogram( new TH2F ("zpt_rebin_shapes_NRBctrl"+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
        h2->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
        h2->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
        h2->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");
    }


    //RJ
    //Phi modulation correction
    //Hadronic recoil correction
    mon.addHistogram( new TH2F ("clusteredMetL_vsPhi",";Z-#phi; Long. clustered E_{T}^{miss} [GeV];", 60,-TMath::Pi(),TMath::Pi(), 75,-250,500) );
    mon.addHistogram( new TH2F ("clusteredMetT_vsPhi",";Z-#phi; Trans. clustered E_{T}^{miss} [GeV];", 60,-TMath::Pi(),TMath::Pi(), 75,-250,500) );
    mon.addHistogram( new TH2F ("RecoilMetL_vsPhi",";photon-#phi; U_{L} [GeV];", 60,-TMath::Pi(),TMath::Pi(), 75,-250,500) );
    mon.addHistogram( new TH2F ("RecoilMetT_vsPhi",";photon-#phi; -U_{T} [GeV];", 60,-TMath::Pi(),TMath::Pi(), 75,-250,500) );

    //mon.addHistogram( new TH1F( "RecoilL", ";U_{L} [GeV];Events", 75,-250,500) );
    //mon.addHistogram( new TH1F( "RecoilT", ";-U_{T} [GeV];Events", 75,-250,500) );

    mon.addHistogram( new TH2F ("RecoilT_vsqT",";q_{T} [GeV]; -U_{T} [GeV];Events",99,qtaxis,75,-250,500) );
    mon.addHistogram( new TH2F ("RecoilL_vsqT",";q_{T} [GeV]; U_{L} [GeV];Events",99,qtaxis,75,-250,500) );





    //##############################################
    //######## GET READY FOR THE EVENT LOOP ########
    //##############################################

    //open the file and get events tree
    ZZ2l2nuSummaryHandler evSummaryHandler;
    TFile *file = TFile::Open(url);
    printf("Looping on %s\n",url.Data());
    if(file==0) return -1;
    if(file->IsZombie()) return -1;
    if( !evSummaryHandler.attachToTree( (TTree *)file->Get(dirname) ) ) {
        file->Close();
        return -1;
    }


    //check run range to compute scale factor (if not all entries are used)
    const Int_t totalEntries= evSummaryHandler.getEntries();
    float rescaleFactor( evEnd>0 ?  float(totalEntries)/float(evEnd-evStart) : -1 );
    if(evEnd<0 || evEnd>evSummaryHandler.getEntries() ) evEnd=totalEntries;
    if(evStart > evEnd ) {
        file->Close();
        return -1;
    }

    //MC normalization (to 1/pb)
    float cnorm=1.0;
    if(isMC) {
        TH1F* cutflowH = (TH1F *) file->Get("evAnalyzer/h2zz/cutflow");
        if(cutflowH) cnorm=cutflowH->GetBinContent(1);
        if(rescaleFactor>0) cnorm /= rescaleFactor;
        printf("cnorm = %f\n",cnorm);
    }
    Hcutflow->SetBinContent(1,cnorm);


    //pileup weighting: based on vtx for now...
    std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
    std::vector<float> dataPileupDistribution;
    for(unsigned int i=0; i<dataPileupDistributionDouble.size(); i++) {
        dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);
    }
    std::vector<float> mcPileupDistribution;
    bool useObservedPU(true);
    //bool useObservedPU(use2011Id);
    if(!use2011Id && url.Contains("toZZto2L")) useObservedPU=true;
    if(isMC) {
        TString puDist("evAnalyzer/h2zz/pileuptrue");
        if(useObservedPU) puDist="evAnalyzer/h2zz/pileup";
        TH1F* histo = (TH1F *) file->Get(puDist);
        if(!histo)std::cout<<"pileup histogram is null!!!\n";
        for(int i=1; i<=histo->GetNbinsX(); i++) {
            mcPileupDistribution.push_back(histo->GetBinContent(i));
        }
        delete histo;
        if(dataPileupDistribution.size()==0) dataPileupDistribution=mcPileupDistribution;
    }
    while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);

    gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
    edm::LumiReWeighting *LumiWeights=0;
    PuShifter_t PuShifters;
    if(isMC) {
        LumiWeights= new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
        PuShifters=getPUshifters(dataPileupDistribution,0.05);
    }

    //event Categorizer
    //EventCategory eventCategoryInst(4); //jet(0,>=1)+vbf binning
    EventCategory eventCategoryInst(1); //jet(0,1,>=2) binning //RJ



    //##############################################
    //########           EVENT LOOP         ########
    //##############################################
    //loop on all the events
    printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
    printf("Scanning the ntuple :");
    int treeStep = (evEnd-evStart)/50;
    if(treeStep==0)treeStep=1;
    DuplicatesChecker duplicatesChecker;
    int nDuplicates(0);
    for( int iev=evStart; iev<evEnd; iev++) {
        if((iev-evStart)%treeStep==0) {
            printf(".");
            fflush(stdout);
        }

        //##############################################   EVENT LOOP STARTS   ##############################################
        //load the event content from tree
        evSummaryHandler.getEntry(iev);
        ZZ2l2nuSummary_t &ev=evSummaryHandler.getEvent();
        if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) {
            nDuplicates++;
            continue;
        }
        PhysicsEvent_t phys=getPhysicsEventFrom(ev);

        //event category
        bool isSameFlavor(ev.cat==MUMU || ev.cat==EE);
        TString tag_cat;
        switch(ev.cat) {
        case MUMU :
            tag_cat = "mumu";
            break;
        case EE   :
            tag_cat = "ee";
            break;
        case EMU  :
            tag_cat = "emu";
            break;
        default   :
            continue;
        }
        //      if(isMC && mctruthmode==1 && !isDYToLL(ev.mccat) && !isZZ2l2nu(ev.mccat) ) continue;
        if(isMC && mctruthmode==1 && !isDYToLL(ev.mccat) ) continue;
        if(isMC && mctruthmode==2 && !isDYToTauTau(ev.mccat) ) continue;

        //require compatibilitiy of the event with the PD
        bool hasTrigger(false);
        bool hasEEtrigger = ev.triggerType & 0x1;
        bool hasMMtrigger = (ev.triggerType >> 1 ) & 0x1;
        bool hasEMtrigger = (ev.triggerType >> 2 ) & 0x1;
        bool hasMtrigger  = (ev.triggerType >> 3 ) & 0x1;
        if(!isMC) {
            if(ev.cat!=fType) continue;

            if(ev.cat==EE   && !hasEEtrigger) continue;
            if(ev.cat==MUMU && !(hasMMtrigger||hasMtrigger) ) continue;
            if(ev.cat==EMU  && !hasEMtrigger) continue;

            //this is a safety veto for the single mu PD
            if(isSingleMuPD) {
                if(!hasMtrigger) continue;
                if(hasMtrigger && hasMMtrigger) continue;
            }

            hasTrigger=true;
        } else {
            if(ev.cat==EE   && hasEEtrigger) hasTrigger=true;
            if(ev.cat==MUMU && (hasMMtrigger || hasMtrigger) ) hasTrigger=true;
            if(ev.cat==EMU  && hasEMtrigger) hasTrigger=true;
            if(use2011Id) hasTrigger = true;
        }

        //prepare the tag's vectors for histo filling
        std::vector<TString> tags(1,"all");

        //pileup weight
        float weight = 1.0;
        double TotalWeight_plus = 1.0;
        double TotalWeight_minus = 1.0;
        if(isMC) {
            weight            = LumiWeights->weight(useObservedPU ? ev.ngenITpu : ev.ngenTruepu);
            TotalWeight_plus  = PuShifters[PUUP]->Eval(useObservedPU ? ev.ngenITpu : ev.ngenTruepu);
            TotalWeight_minus = PuShifters[PUDOWN]->Eval(useObservedPU ? ev.ngenITpu : ev.ngenTruepu);
        }
        Hcutflow->Fill(1,1);
        Hcutflow->Fill(2,weight);
        Hcutflow->Fill(3,weight*TotalWeight_minus);
        Hcutflow->Fill(4,weight*TotalWeight_plus);


	////////////////////////////////////////////////// what?
	if(isMC) {
		TString input_File = url.Data();
		bool isZH(false);
		bool isWZ(false);
		bool isZZ(false);
		if(input_File.Contains("TeV_ZH")) { /*fprintf(outTxtFile,"ZH \n");*/ isZH=true;}
		if(input_File.Contains("TeV_WZ")) { /*fprintf(outTxtFile,"WZ \n");*/ isWZ=true;}
		if(input_File.Contains("TeV_ZZ")) { /*fprintf(outTxtFile,"ZZ \n");*/ isZZ=true;}

        	//NLO EWK Signal correction
		//double h_genPt = phys.genhiggs[0].pt();
/*
		int ngenleps = phys.genleptons.size();
		fprintf(outTxtFile,"Nleps:%d [",ngenleps);
		for(int j=0; j<ngenleps; j++) fprintf(outTxtFile,"%d ",phys.genleptons[j].id);
		fprintf(outTxtFile,"]\n");

		int Nneutrinos = phys.genneutrinos.size();
		fprintf(outTxtFile,"Nv:%d [",Nneutrinos);
		for(int j=0; j<Nneutrinos; j++) fprintf(outTxtFile,"%d ",phys.genneutrinos[j].id);
		fprintf(outTxtFile,"]\n");		
*/
		//for(int j=0; j<Nneutrinos; j++) fprintf(outTxtFile,"px: %f, py: %f, pz: %f \n",phys.genneutrinos[j].px(),
		//					phys.genneutrinos[j].py(), phys.genneutrinos[j].pz());



		if(isZH){
			//double pt_zll = (phys.genleptons[0]+phys.genleptons[1]).pt();
			double pt_zll = phys.genhiggs[0].pt();
			float ewk_nloweightsZH = ZHUtils::weightNLOEWKsignal(pt_zll);

			mon.fillHisto("v_pt",tags,pt_zll,weight);
			weight *= ewk_nloweightsZH;
			//fprintf(outTxtFile,"ZH:%f \n",ewk_nloweightsZH);
			mon.fillHisto("v_pt_ewkcorr",tags,pt_zll,weight);
		}

		if(isZZ){
			double pt_zll = (phys.genleptons[0]+phys.genleptons[1]).pt();
			//double pt_zvv = phys.genmet[0].pt();
			double pt_zvv = (phys.genneutrinos[0]+phys.genneutrinos[1]).pt();
			double trailing_zzpt = pt_zvv;
			if(pt_zll<pt_zvv) trailing_zzpt = pt_zll;
			double ewk_nloweightsZZ = ZHUtils::weightNLOEWKzz(trailing_zzpt);

			mon.fillHisto("v_pt",tags,trailing_zzpt,weight);
			weight *= ewk_nloweightsZZ;
			//fprintf(outTxtFile,"ZZ:%f \n",ewk_nloweightsZZ);
			mon.fillHisto("v_pt_ewkcorr",tags,trailing_zzpt,weight);
		}

		if(isWZ){
			int ngenleps = phys.genleptons.size();
			int Nneutrinos = phys.genneutrinos.size();
			if(ngenleps==3 && Nneutrinos==1){
				int id_lep0 = fabs(phys.genleptons[0].id);
				int id_lep1 = fabs(phys.genleptons[1].id);
				int id_lep2 = fabs(phys.genleptons[2].id);

				int id_Lep0 = phys.genleptons[0].id;
				int id_Lep1 = phys.genleptons[1].id;
				int id_Lep2 = phys.genleptons[2].id;
				bool hassameids = (id_lep0==id_lep1)&&(id_lep1==id_lep2);	
				bool isWminus(false);

				double pt_zll, pt_wlv;
				if(!hassameids){
					if(id_lep0==id_lep1){
						pt_zll = (phys.genleptons[0]+phys.genleptons[1]).pt();
						pt_wlv = (phys.genleptons[2]+phys.genneutrinos[0]).pt();
						//fprintf(outTxtFile,"Z candidates: %d - %d\n",id_Lep0,id_Lep1);
						if(id_Lep2<0) isWminus=true;
					}
					else if(id_lep0==id_lep2){
						pt_zll = (phys.genleptons[0]+phys.genleptons[2]).pt();
						pt_wlv = (phys.genleptons[1]+phys.genneutrinos[0]).pt();
						//fprintf(outTxtFile,"Z candidates: %d - %d\n",id_Lep0,id_Lep2);
						if(id_Lep1<0) isWminus=true;
					}
					else if(id_lep1==id_lep2){
						pt_zll = (phys.genleptons[1]+phys.genleptons[2]).pt();
						pt_wlv = (phys.genleptons[0]+phys.genneutrinos[0]).pt();
						//fprintf(outTxtFile,"Z candidates: %d - %d\n",id_Lep1,id_Lep2);
						if(id_Lep0<0) isWminus=true;
					}
				}
				else{ // if 3 leptons has same ids

					//find the one different lepton
					int sum_lepid = id_Lep0+id_Lep1+id_Lep2;
					if(id_Lep0 == -sum_lepid){
						double mass1 = fabs((phys.genleptons[0]+phys.genleptons[1]).mass()-91.);
						double mass2 = fabs((phys.genleptons[0]+phys.genleptons[2]).mass()-91.);
						if(mass1 < mass2) {
							pt_zll = (phys.genleptons[0]+phys.genleptons[1]).pt();
							pt_wlv = (phys.genleptons[2]+phys.genneutrinos[0]).pt();
							//fprintf(outTxtFile,"Z candidates: %d - %d\n",id_Lep0,id_Lep1);
							if(id_Lep2<0) isWminus=true;
						}else{
							pt_zll = (phys.genleptons[0]+phys.genleptons[2]).pt();
							pt_wlv = (phys.genleptons[1]+phys.genneutrinos[0]).pt();
							//fprintf(outTxtFile,"Z candidates: %d - %d\n",id_Lep0,id_Lep2);
							if(id_Lep1<0) isWminus=true;
						} 
					}
					else if(id_Lep1 == -sum_lepid){
						double mass1 = fabs((phys.genleptons[1]+phys.genleptons[0]).mass()-91.);
						double mass2 = fabs((phys.genleptons[1]+phys.genleptons[2]).mass()-91.);
						if(mass1 < mass2) {
							pt_zll = (phys.genleptons[1]+phys.genleptons[0]).pt();
							pt_wlv = (phys.genleptons[2]+phys.genneutrinos[0]).pt();
							//fprintf(outTxtFile,"Z candidates: %d - %d\n",id_Lep1,id_Lep0);
							if(id_Lep2<0) isWminus=true;
						}else{
							pt_zll = (phys.genleptons[1]+phys.genleptons[2]).pt();
							pt_wlv = (phys.genleptons[0]+phys.genneutrinos[0]).pt();
							//fprintf(outTxtFile,"Z candidates: %d - %d\n",id_Lep1,id_Lep2);
							if(id_Lep0<0) isWminus=true;
						}
					}
					else if(id_Lep2 == -sum_lepid){
						double mass1 = fabs((phys.genleptons[2]+phys.genleptons[0]).mass()-91.);
						double mass2 = fabs((phys.genleptons[2]+phys.genleptons[1]).mass()-91.);
						if(mass1 < mass2) {
							pt_zll = (phys.genleptons[2]+phys.genleptons[0]).pt();
							pt_wlv = (phys.genleptons[1]+phys.genneutrinos[0]).pt();
							//fprintf(outTxtFile,"Z candidates: %d - %d\n",id_Lep2,id_Lep0);
							if(id_Lep1<0) isWminus=true;
						}else{
							pt_zll = (phys.genleptons[2]+phys.genleptons[1]).pt();
							pt_wlv = (phys.genleptons[0]+phys.genneutrinos[0]).pt();
							//fprintf(outTxtFile,"Z candidates: %d - %d\n",id_Lep2,id_Lep1);
							if(id_Lep0<0) isWminus=true;
						}
					}
				}//if 3 leptons has same ids

				double trailing_pt = pt_wlv;
                        	if(pt_zll<pt_wlv) trailing_pt = pt_zll;
				//fprintf(outTxtFile,"isWminus:%d \n",isWminus);
				double ewk_nloweightsWZ = ZHUtils::weightNLOEWKwz(trailing_pt);
  	                      	//double ewk_nloweightsWZ = ZHUtils::weightNLOEWKwplusz(trailing_pt);
				//if(isWminus) ewk_nloweightsWZ = ZHUtils::weightNLOEWKwminuz(trailing_pt);

                        	mon.fillHisto("v_pt",tags,trailing_pt,weight);
                        	weight *= ewk_nloweightsWZ;
                        	//fprintf(outTxtFile,"WZ:%f \n",ewk_nloweightsWZ);
                        	mon.fillHisto("v_pt_ewkcorr",tags,trailing_pt,weight);
				

			} //ngenleps==3 && Nneutrinos==1
		} //isWZ

	}
	//////////////////////////////////////////////////


        //MET variables
        LorentzVector rawMetP4=phys.met[2];
        if(use2011Id) rawMetP4=phys.met[0]; //V3 Test
        LorentzVector fullTypeIMetP4=phys.met[0];
        LorentzVector mvaMetP4=phys.met[7];

        //apply JER base corrections to jets (and compute associated variations on the MET variable)
        // std PF
        std::vector<PhysicsObjectJetCollection> variedAJets;
        LorentzVectorCollection zvvs;
        PhysicsObjectJetCollection &recoJets = ( useCHS ? phys.ajets : phys.jets);
        METUtils::computeVariation(recoJets, phys.leptons, rawMetP4, variedAJets, zvvs, &jecUnc);
        if(!useJERsmearing) zvvs[0] = rawMetP4; // this stinks a bit...

        //
        // LEPTON ANALYSIS
        //
        LorentzVector lep1=phys.leptons[0];
        LorentzVector lep2=phys.leptons[1];
        //Muon scale and resolution corrections
        if(fabs(phys.leptons[0].id)==13) {
            TLorentzVector tmpLep1(lep1.Px(), lep1.Py(), lep1.Pz(), lep1.E());
            int tmpCh1 = phys.leptons[0].id < 0.0 ? -1 : 1; // "id" has same sign as charge here (opposite to "genid")
            corrector_->applyPtCorrection(tmpLep1, tmpCh1);
            if( isMC && (!use2011Id) ) corrector_->applyPtSmearing(tmpLep1, tmpCh1);
            lep1 = LorentzVector(tmpLep1.Px(), tmpLep1.Py(), tmpLep1.Pz(), tmpLep1.E());
        }
        if(fabs(phys.leptons[1].id)==13) {
            TLorentzVector tmpLep2(lep2.Px(), lep2.Py(), lep2.Pz(), lep2.E());
            int tmpCh2 = phys.leptons[0].id < 0.0 ? -1 : 1; // "id" has same sign as charge here (opposite to "genid")
            corrector_->applyPtCorrection(tmpLep2, tmpCh2);
            if( isMC && (!use2011Id) ) corrector_->applyPtSmearing(tmpLep2, tmpCh2);
            lep2 = LorentzVector(tmpLep2.Px(), tmpLep2.Py(), tmpLep2.Pz(), tmpLep2.E());
        }
        LorentzVector zll(lep1+lep2);
        double dphi2l = fabs(deltaPhi(lep1.phi(),lep2.phi()));
	bool passdphi2l(true);
	//bool passdphi2l(cos(dphi2l)>0);

        //two opposite sign leptons
        int  id1=phys.leptons[0].id;
        int  id2=phys.leptons[1].id;
        bool passOppositeSign((id1*id2)<0);

        bool passId(true);
        passId &= passOppositeSign;
        bool passIdAndIso(true);
        passIdAndIso &= passOppositeSign;
        bool passZmass(fabs(zll.mass()-91)<15); //RJ changed

        bool isZsideBand( (zll.mass()>40 && zll.mass()<70) || (zll.mass()>110 && zll.mass()<200));
        bool isZsideBandPlus( (zll.mass()>110 && zll.mass()<200));
        bool passZpt(zll.pt()>50); //for Gamma+Jets, can also be used for control plots

        //check alternative selections for the dilepton
        double llScaleFactor(1.0),llTriggerEfficiency(1.0);
        LorentzVector genZP4(0,0,0,0); // for checks on Sherpa ZZ
        int genmatchid[2] = {-1, -1};
        double genmatchdr[2] = {0.1, 0.1};
        for(int ilep=0; ilep<2; ilep++) {
            TString lepStr( fabs(phys.leptons[ilep].id)==13 ? "mu" : "e");

            //generator level matching
            //  int matchid(0);
            //LorentzVector genP4(0,0,0,0);
            for(size_t igl=0; igl<phys.genleptons.size(); igl++) {
                //if(deltaR(phys.genleptons[igl],phys.leptons[ilep])>0.1) continue;
                if(ilep==1 && int(igl)==genmatchid[0]) continue;
                if(deltaR(phys.genleptons[igl],phys.leptons[ilep])<genmatchdr[ilep]) {
                    genmatchdr[ilep] = deltaR(phys.genleptons[igl],phys.leptons[ilep]);
                    genmatchid[ilep] = igl;
                }
            }
            if(genmatchid[0]>-1 && genmatchid[1]>-1) {
                //genP4=phys.genleptons[igl];
                //matchid=phys.genleptons[igl].id;
                genZP4 = phys.genleptons[0] + phys.genleptons[1];
            }

            //id and isolation
            int lpid=phys.leptons[ilep].pid;
            float relIso2011    = phys.leptons[ilep].relIsoRho(ev.rho);
            float relIso = (lepStr=="mu") ?
                           phys.leptons[ilep].pfRelIsoDbeta() :
                           phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho,ev.en_sceta[lpid]); //RENJIE
            std::vector<int> passIds;
            std::map<int,bool> passIsos;
            bool hasGoodId(false), isIso(false);
            if(fabs(phys.leptons[ilep].id)==13) {
                if( hasObjectId(ev.mn_idbits[lpid], MID_LOOSE) )    {
                    passIds.push_back(0);
                    passIsos[0]=(relIso<0.2);
                }
                if( hasObjectId(ev.mn_idbits[lpid], MID_TIGHT) )    {
                    passIds.push_back(1);
                    passIsos[1]=(relIso<0.2);
                    if(!use2011Id) {
                        hasGoodId=true;
                        isIso=passIsos[0];
                    }
                }
                if( hasObjectId(ev.mn_idbits[lpid], MID_VBTF2011) ) {
                    passIds.push_back(2);
                    passIsos[2]=(relIso2011<0.15);
                    if(use2011Id) {
                        hasGoodId=true;
                        isIso=passIsos[2];
                    }
                }
                if( hasObjectId(ev.mn_idbits[lpid], MID_SOFT) )     {
                    passIds.push_back(3);
                    passIsos[3]=true;
                }
                if(use2011Id) {
                    try {
                        llScaleFactor *= muonScaleFactor(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()));
                        llTriggerEfficiency *= muonTriggerEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()));
                    } catch(std::string &e) {
                    }
                } else {
                    llScaleFactor *= 1;
                    llTriggerEfficiency *= 1.0; //muonTriggerEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),2012);
                }
            } else {
                int wps[]= {EgammaCutBasedEleId::LOOSE,EgammaCutBasedEleId::MEDIUM, EID_VBTF2011, EgammaCutBasedEleId::VETO};
                for(int iwp=0; iwp<4; iwp++) {
                    if(iwp==2 && hasObjectId(ev.en_idbits[lpid], EID_VBTF2011)) {
                        passIds.push_back(2);
                        passIsos[2]=(relIso2011<0.10);
                        if(use2011Id) {
                            hasGoodId=true;
                            isIso=passIsos[2];
                            try {
                                llScaleFactor *= electronScaleFactor(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()));
                                llTriggerEfficiency *= electronTriggerEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()));
                            } catch(std::string &e) {
                            }
                        }
                    } else {
                        bool passWp = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::WorkingPoint(wps[iwp]),
                                      (fabs(phys.leptons[ilep].eta())<1.4442),
                                      phys.leptons[ilep].pt(), phys.leptons[ilep].eta(),
                                      ev.en_detain[lpid],  ev.en_dphiin[lpid], ev.en_sihih[lpid], ev.en_hoe[lpid],
                                      ev.en_ooemoop[lpid], phys.leptons[ilep].d0, phys.leptons[ilep].dZ,
                                      0., 0., 0.,
                                      !hasObjectId(ev.en_idbits[lpid], EID_CONVERSIONVETO),0,ev.rho);
                        if(passWp) {
                            passIds.push_back(iwp);
                            passIsos[iwp]=(relIso<0.15);
                            if(wps[iwp]==EgammaCutBasedEleId::MEDIUM && !use2011Id) {
                                hasGoodId=true;
                                isIso=passIsos[iwp];
                            }
                        }
                        if(!use2011Id) {
                            llScaleFactor *= 1;
                            llTriggerEfficiency *= 1.0; //electronTriggerEfficiency(phys.leptons[ilep].pt(),fabs(phys.leptons[ilep].eta()),2012);
                        }
                    }
                }
            }
            if(!hasGoodId)  {
                passId &=false;    //RJ
                passIdAndIso &=false;
            }
            if(!isIso) passIdAndIso &=false; //RJ

            //        fill control histograms (constrained to the Z mass)
            if(passZmass && isSameFlavor) {
                if(hasGoodId) {
                    mon.fillHisto(lepStr+"reliso",     tags, use2011Id? relIso2011 : relIso,   weight);
                }
            }
        }

        tags.push_back(tag_cat);
        if(tag_cat=="mumu" || tag_cat=="ee") tags.push_back("ll");


        //
        // 3rd LEPTON ANALYSIS
        // 1.for lepton veto
        // 2.for WZ Control Sample
        bool pass3dLeptonVeto(true);
        bool has3dLepton(false);
        int nextraleptons(0), nextraleptons_WZ(0);
        std::vector<int> nextraleptonsid;
        std::vector<LorentzVector> extraLeptonsP4, extraLeptonsP4_WZ;
        for(size_t ilep=2; ilep<phys.leptons.size(); ilep++) {
            //lepton type
            bool isGood(false), isGood_WZ(false);
            int lpid=phys.leptons[ilep].pid;
            if(fabs(phys.leptons[ilep].id)==13) {
                if(!use2011Id) {
                    isGood = (hasObjectId(ev.mn_idbits[lpid], MID_LOOSE) && phys.leptons[ilep].pfRelIsoDbeta()<0.2);
                    isGood |= (hasObjectId(ev.mn_idbits[lpid], MID_SOFT) && phys.leptons[ilep].pt()>3);
                    //RENJIE
                    isGood_WZ = (hasObjectId(ev.mn_idbits[lpid], MID_LOOSE) && phys.leptons[ilep].pfRelIsoDbeta()<0.2 && phys.leptons[ilep].pt()>20);
                } else {
                    isGood = (hasObjectId(ev.mn_idbits[lpid], MID_VBTF2011) && phys.leptons[ilep].relIsoRho(ev.rho)<0.15 && phys.leptons[ilep].pt()>10);
                    isGood |= (hasObjectId(ev.mn_idbits[lpid], MID_SOFT2011) && phys.leptons[ilep].pt()>3);
                    //RENJIE
                    isGood_WZ = (hasObjectId(ev.mn_idbits[lpid], MID_VBTF2011) && phys.leptons[ilep].relIsoRho(ev.rho)<0.15 && phys.leptons[ilep].pt()>20);
                }
            } else {
                if(!use2011Id) { //RENJIE
                    isGood = ( hasObjectId(ev.en_idbits[lpid],EID_VETO) && phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho,ev.en_sceta[lpid])<0.15 && phys.leptons[ilep].pt()>10);
                    //RENJIE
                    isGood_WZ=(hasObjectId(ev.en_idbits[lpid],EID_VETO) && phys.leptons[ilep].ePFRelIsoCorrected2012(ev.rho,ev.en_sceta[lpid])<0.15 && phys.leptons[ilep].pt()>20);
                } else {
                    isGood = ( hasObjectId(ev.en_idbits[lpid],EID_VBTF2011) && phys.leptons[ilep].relIsoRho(ev.rho)<0.1 && phys.leptons[ilep].pt()>10);
                    //RENJIE
                    isGood_WZ=(hasObjectId(ev.en_idbits[lpid],EID_VBTF2011) && phys.leptons[ilep].relIsoRho(ev.rho)<0.1 && phys.leptons[ilep].pt()>20);
                }
            }
            nextraleptons += isGood;
            nextraleptons_WZ += isGood_WZ;
            if(isGood_WZ) nextraleptonsid.push_back(phys.leptons[ilep].id);

            //if(!isGood) continue;
            //Muon scale and resolution corrections
            // probably unnecessary
            LorentzVector tmpLep = phys.leptons[ilep];
            if(fabs(phys.leptons[ilep].id)==13) {
                TLorentzVector tmpTLep(tmpLep.Px(), tmpLep.Py(), tmpLep.Pz(), tmpLep.E());
                int tmpCh = phys.leptons[ilep].id < 0.0 ? -1 : 1; // "id" has same sign as charge here (opposite to "genid")
                corrector_->applyPtCorrection(tmpTLep, tmpCh);
                if( isMC && (!use2011Id) ) corrector_->applyPtSmearing(tmpTLep, tmpCh);
                tmpLep = LorentzVector(tmpTLep.Px(), tmpTLep.Py(), tmpTLep.Pz(), tmpTLep.E());
            }
            if(isGood) extraLeptonsP4.push_back(tmpLep);
            if(isGood_WZ) extraLeptonsP4_WZ.push_back(tmpLep);
            //extraLeptonsP4.push_back( phys.leptons[ilep] );
        }
        pass3dLeptonVeto=(nextraleptons==0);
        has3dLepton=(nextraleptons_WZ>0);

        //
        //STD PF JET ANALYSIS
        //
        bool passJetveto(true);
        bool passBveto(true);
        bool passRedMet(true);
        bool passDphijmet(true);
        bool passBalanceCut(true);
        PhysicsObjectJetCollection &aJets = ( useJERsmearing ? variedAJets[0] : recoJets );
        PhysicsObjectJetCollection aGoodIdJets;
        LorentzVector aClusteredMetP4(zll);
        aClusteredMetP4 *= -1;
        LorentzVector Recoil(zvvs[0]);
        Recoil *= -1;
        Recoil -= zll;
        int nABtags(0),nAJetsGood30(0), nCSVMtags(0), nCSVTtags(0);
        float mindphijmet(999999.),mindphijmet15(999999.);
        for(size_t ijet=0; ijet<aJets.size(); ijet++) {
            if(aJets[ijet].pt()<15) continue;

            float idphijmet( fabs(deltaPhi(aJets[ijet].phi(),zvvs[0].phi()) ) );
            if(aJets[ijet].pt()>15) if(idphijmet<mindphijmet15)  mindphijmet15=idphijmet;
            if(aJets[ijet].pt()>30) if(idphijmet<mindphijmet)  mindphijmet=idphijmet;

            bool isGoodJet = hasObjectId(aJets[ijet].pid,JETID_LOOSE);
            if(usePUsubJetId) isGoodJet = hasObjectId(aJets[ijet].pid,JETID_CUTBASED_LOOSE);
            mon.fillHisto("pfjetbeta",     tags,aJets[ijet].beta,     weight);
            mon.fillHisto("pfjetmva",      tags,aJets[ijet].pumva,    weight);
            if(!isGoodJet) continue;


            aClusteredMetP4 -= aJets[ijet];

            if(useJetsOnlyInTracker && fabs(aJets[ijet].eta())>2.5) continue;

            aGoodIdJets.push_back(aJets[ijet]);
            if(aJets[ijet].pt()>30)nAJetsGood30++;

            if(aJets[ijet].pt()>20 && fabs(aJets[ijet].eta())<2.5)  nABtags += (aJets[ijet].btag2>0.244);
            if(aJets[ijet].pt()>20 && fabs(aJets[ijet].eta())<2.5)  nCSVMtags += (aJets[ijet].btag2>0.679);
            if(aJets[ijet].pt()>20 && fabs(aJets[ijet].eta())<2.5)  nCSVTtags += (aJets[ijet].btag2>0.898);
        }
        //passJetveto=(nAJetsGood30==0);
        passBveto=(nABtags==0);
        //bool passMBveto=(nCSVMtags==0);
        bool passTBveto=(nCSVTtags==0);
        //passDphijmet=(mindphijmet>0.5);
        //if(nAJetsGood30==0) passDphijmet=(mindphijmet15>0.5);
        passBalanceCut=(zvvs[0].pt()/zll.pt()>0.8 && zvvs[0].pt()/zll.pt()<1.2);
        bool passBalanceCutWWCtrl=(zvvs[0].pt()/zll.pt()>0.4 && zvvs[0].pt()/zll.pt()<1.8);


        //ad-hoc cut for obvious correlations between MET and a lepton
        double dphil1met=fabs(deltaPhi(lep1.phi(),zvvs[0].phi()));
        double dphil2met=fabs(deltaPhi(lep2.phi(),zvvs[0].phi()));
        bool passLMetVeto(true);
        if(!use2011Id && zvvs[0].pt()>60 && min(dphil1met,dphil2met)<0.2) passLMetVeto=false;

        //other mets
        METUtils::stRedMET aRedMetOut;
        LorentzVector null(0,0,0,0);
        LorentzVector aRedMet=METUtils::redMET(METUtils::INDEPENDENTLYMINIMIZED, zll, 0, null, 0, aClusteredMetP4, zvvs[0], true, &aRedMetOut);
        double aRedMetL=aRedMetOut.redMET_l;
        double aRedMetT=aRedMetOut.redMET_t;
        TVector2 aClusteredMet2(aClusteredMetP4.px(),aClusteredMetP4.py());
        double clusteredMetL=aRedMetOut.a_l*aClusteredMet2;
        double clusteredMetT=aRedMetOut.a_t*aClusteredMet2;
        passRedMet=(aRedMet.pt()>110); //ReducedMET: aRedMet.pt(); PFMET: zvvs[0].pt()
        bool passRedMetWWCtrl=(aRedMet.pt()>65);

        TVector2 RecoilMet2(Recoil.px(),Recoil.py());
        double RecoilMetL=aRedMetOut.a_l*RecoilMet2;
        double RecoilMetT=aRedMetOut.a_t*RecoilMet2;


        //RJ compute dphi(zpt,redMet)
        double dphiZllmet=fabs(deltaPhi(zll.phi(),zvvs[0].phi()));
        double dphiZllredMet=fabs(deltaPhi(zll.phi(),aRedMet.phi()));
        //RJ
        bool passdphiZllmetCut(dphiZllmet>2.6);

        //transverse masses
        double aMT=METUtils::transverseMass(zll,zvvs[0],true);
        double aMTmassless=METUtils::transverseMass(zll,zvvs[0],false);
        double balanceDif=fabs(zvvs[0].pt()-zll.pt())/zll.pt();
        TVector2 dil2(zll.px(),zll.py());
        TVector2 met2(zvvs[0].px(),zvvs[0].py());
        double axialMet=dil2*met2;
        axialMet /= -zll.pt();
        double pfMetCompL = aRedMetOut.a_l*met2;
        double pfMetCompT = aRedMetOut.a_t*met2;

        //Collins-Soper Frame
        double colin_soper = ZHUtils::Collins_Soper(lep1,lep2);


        //RJ
        bool passMTcut(aMT>220 && aMT<1200);
	//passMTcut &= passdphi2l;

        int selWord1(0);
        int selWord2(0);
        selWord1 |= hasTrigger;
        selWord1 |= (passId<<1);
        selWord1 |= (passIdAndIso<<2);
        selWord1 |= (passZmass<<3);
        selWord1 |= (passZpt<<4);
        selWord1 |= (pass3dLeptonVeto<<5);
        selWord1 |= (passBveto<<6);

        selWord2 |= (passLMetVeto<<1);
        selWord2 |= (passdphiZllmetCut<<2);
        selWord2 |= (passRedMet<<3);
        selWord2 |= (passBalanceCut<<4);
        selWord2 |= (passMTcut<<5);

        TString evCat(tag_cat);
        {
            int eventSubCat  = eventCategoryInst.Get(phys,&aGoodIdJets);
            TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);
            evCat+=tag_subcat;
        }
        //fprintf(outTxtFile,"%d %d %d %d %d %f %f %s\n",ev.run,ev.lumi,ev.event,selWord,nAJetsGood30,zvvs[0].pt(),aMT,evCat.Data());
        //fprintf(outTxtFile,"%d %d %d %d %d %s\n",ev.run,ev.lumi,ev.event,selWord1,selWord2,evCat.Data());
        fprintf(outTxtFile,"%d %d %d %d%d%d%d%d%d%d%d%d%d%d %s %d\n",ev.run,ev.lumi,ev.event,
                passId,passIdAndIso,passZmass,passZpt,pass3dLeptonVeto,passBveto,passLMetVeto,passdphiZllmetCut,
                passRedMet,passBalanceCut,passMTcut,
                evCat.Data(),hasTrigger);

        //
        //RUN PRESELECTION AND CONTROL PLOTS
        //
        if(isMC && use2011Id) weight *= llScaleFactor*llTriggerEfficiency;
        if(hasTrigger)                 {
            mon.fillHisto("eventflow",tags,0,weight);
        }
        if(hasTrigger && passId)       {
            mon.fillHisto("eventflow",tags,1,weight);
        }
        if(hasTrigger && passIdAndIso) {
            mon.fillHisto("eventflow",tags,2,weight);
        } else continue;

        // Control Regions
        bool WZ_ctrl_pt ( passZpt && nextraleptons==1 && passZmass && passBveto );
        bool WZ_ctrl_m  ( passZpt && nextraleptons==1 && passBveto );
        bool T_ctrl     ( passZpt && nABtags>0 && isZsideBand );
        bool T_ctrl_side( passZpt && nABtags>0 );
        bool WW_ctrl    ( passZpt && passBveto && (nAJetsGood30==0) && passBalanceCut && passDphijmet && pass3dLeptonVeto && isZsideBand );
        bool WW_ctrl_1  ( passZpt && isZsideBand );
        bool WW_ctrl_2  ( passZpt && isZsideBand && passBveto );
        bool WW_ctrl_3  ( passZpt && isZsideBand && passBveto && passBalanceCut && (nAJetsGood30==0) );
        bool WW_ctrl_4  ( passZpt && isZsideBand && passBveto && passBalanceCut && (nAJetsGood30==0) && passDphijmet);
        bool WW_ctrl_5  ( passZpt && isZsideBand && passBveto && passBalanceCut && (nAJetsGood30==0) && passDphijmet && pass3dLeptonVeto);
        //bool WW_ctrl_m  ( passZpt && passRedMet  && passBveto && passBalanceCut && (nAJetsGood30==0) && passDphijmet && pass3dLeptonVeto);
        bool WWCtrl_RJ  ( passZpt && passRedMetWWCtrl  && passTBveto && passBalanceCutWWCtrl && (nAJetsGood30==0) && pass3dLeptonVeto);
        if(WZ_ctrl_pt)  mon.fillHisto("Ctrl_WZ_RedMet",tags,aRedMet.pt(),weight);
        if(WZ_ctrl_m) {
            mon.fillHisto("Ctrl_WZ_zmass",tags,zll.mass(),weight);
            if(passZmass) mon.fillHisto("Ctrl_WZ_Mt", tags, sqrt((2*zvvs[0].pt()*extraLeptonsP4[0].pt())*(1-cos(deltaPhi(zvvs[0].phi(),extraLeptonsP4[0].phi())))), weight);
            if(passZmass) mon.fillHisto("Ctrl_WZ_Mt1", tags, sqrt((2*zvvs[0].pt()*extraLeptonsP4[0].pt())*(1-cos(deltaPhi(zvvs[0].phi(),extraLeptonsP4[0].phi())))), weight);
        }
        if(passZmass && passRedMet && nextraleptons==1) mon.fillHisto("Ctrl_WZ_Mt2", tags, sqrt((2*zvvs[0].pt()*extraLeptonsP4[0].pt())*(1-cos(deltaPhi(zvvs[0].phi(),extraLeptonsP4[0].phi())))), weight);
        if(passZmass  && aRedMet.pt()>50 && nextraleptons==1) mon.fillHisto("Ctrl_WZ_Mt3", tags, sqrt((2*zvvs[0].pt()*extraLeptonsP4[0].pt())*(1-cos(deltaPhi(zvvs[0].phi(),extraLeptonsP4[0].phi())))), weight);
        if(passZmass  && aRedMet.pt()>70 && nextraleptons==1) mon.fillHisto("Ctrl_WZ_Mt4", tags, sqrt((2*zvvs[0].pt()*extraLeptonsP4[0].pt())*(1-cos(deltaPhi(zvvs[0].phi(),extraLeptonsP4[0].phi())))), weight);



        //RENJIE's WZ Control Plots
        bool passSameFlavor(fabs(id1)==fabs(id2));
        double new3rdMT=0.0;
        if(has3dLepton && passSameFlavor) {
            int id3=nextraleptonsid[0];
            double relZrange = fabs((lep1+lep2).mass()-91);

            LorentzVector relLep1 = lep1;
            LorentzVector relLep2 = lep2;
            LorentzVector rel3rdLep = extraLeptonsP4_WZ[0];
            /*
            	//checking if 3 leptons have same flavor
            	if(id3==id1){
            		LorentzVector fakeZ = lep2+extraLeptonsP4_WZ[0];
            		double fakeZrange = fabs(fakeZ.mass()-91);
            		if(fakeZrange<relZrange){
            			relLep1=extraLeptonsP4_WZ[0];
            			rel3rdLep=lep1;
            		}
            	}else
            	if(id3==id2){
            		LorentzVector fakeZ = lep1+extraLeptonsP4_WZ[0];
            		double fakeZrange = fabs(fakeZ.mass()-91);
            		if(fakeZrange<relZrange){
            			relLep2=extraLeptonsP4_WZ[0];
            			rel3rdLep=lep2;
            		}
            	}
            */
            LorentzVector trueZll=relLep1+relLep2;
            bool passZmass_WZ(trueZll.mass()>60 && trueZll.mass()<120);
            //bool passZmass_WZ(fabs(trueZll.mass()-91)<15);

            //|M(lll)-91|>|M(ll)-91|
            LorentzVector Trileptons = trueZll+rel3rdLep;
            double trilepmass = fabs(Trileptons.mass()-91);
            double dilepmass = fabs(trueZll.mass()-91);
            bool passtriLepmass(trilepmass>dilepmass);

            //manually add 3rd lepton to pfMET
            LorentzVector newMet = zvvs[0]+rel3rdLep;
            float newMT = METUtils::transverseMass(trueZll,newMet,false);

            //for old WZ control plots
            new3rdMT = METUtils::transverseMass(rel3rdLep,zvvs[0],false);

            if(!passZmass_WZ) continue;
            if(!passBveto) continue;
            if(!passtriLepmass) continue;

            //if(relLep1.pt()<20) continue;
            //if(relLep2.pt()<20) continue;
            //if(rel3rdLep.pt()<20) continue;

            mon.fillHisto("WZ_contrl_pt_lep1",tags,relLep1.pt(),weight);
            mon.fillHisto("WZ_contrl_pt_lep2",tags,relLep2.pt(),weight);
            mon.fillHisto("WZ_contrl_pt_lep3",tags,rel3rdLep.pt(),weight);
            if(fabs(id3)==11)
                mon.fillHisto("WZ_contrl_pt_3rdlepton_e",tags,rel3rdLep.pt(),weight);
            if(fabs(id3)==13)
                mon.fillHisto("WZ_contrl_pt_3rdlepton_mu",tags,rel3rdLep.pt(),weight);

            mon.fillHisto("WZ_contrl_lepmass_pair1",tags,(relLep1+relLep2).mass(),weight);
            mon.fillHisto("WZ_contrl_lepmass_pair2",tags,(relLep1+rel3rdLep).mass(),weight);
            mon.fillHisto("WZ_contrl_lepmass_pair3",tags,(rel3rdLep+relLep2).mass(),weight);


            mon.fillHisto("WZ_contrl_dilepmass",tags,trueZll.mass(),weight);
            mon.fillHisto("WZ_contrl_trilepmass",tags,Trileptons.mass(),weight);


            mon.fillHisto("WZ_contrl_zpt",tags,trueZll.pt(),weight);
            //mon.fillHisto("WZ_contrl_qt_rebin",tags,trueZll.pt(),weight);
            //mon.fillHisto("WZ_contrl_redmet",tags,aRedMet.pt(),weight);
            mon.fillHisto("WZ_contrl_met",tags,zvvs[0].pt(),weight);
            mon.fillHisto("WZ_contrl_newMt",tags,new3rdMT,weight);

            if(fabs(id3) == 11) mon.fillHisto("WZ_contrl_met_extra_e",tags,zvvs[0].pt(),weight);
            if(fabs(id3) == 13) mon.fillHisto("WZ_contrl_met_extra_mu",tags,zvvs[0].pt(),weight);

            //double dphiZmet=fabs(deltaPhi(zll.phi(),zvvs[0].phi()));
            //double dphiZextrlep=fabs(deltaPhi(zll.phi(),extraLeptonsP4[0].phi()));
            //double dphiextrlepMET=fabs(deltaPhi(zvvs[0].phi(),extraLeptonsP4[0].phi()));
            //double dphiZnewMET=fabs(deltaPhi(zll.phi(),newMet.phi()));

            if((zvvs[0].pt()>40)) {

                //if(!(dphiZextrlep<2.6)) continue;
                //mon.fillHisto("WZ_contrl_dphiZMET",tags,dphiZmet,weight);
                //mon.fillHisto("WZ_contrl_dphiextrLepMET",tags,dphiextrlepMET,weight);
                //mon.fillHisto("WZ_contrl_dphiZextrLep",tags,dphiZextrlep,weight);
                //mon.fillHisto("WZ_contrl_dphiZLepMET",tags,dphiZnewMET,weight);
                //mon.fillHisto("WZ_contrl_balance",tags,newMet.pt()/zll.pt(),weight);

                mon.fillHisto("WZ_contrl_Met_final",tags,newMet.pt(),weight);
                mon.fillHisto("WZ_contrl_Mt_final", tags,newMT, weight);
            }

        }



        if(T_ctrl)      mon.fillHisto("Ctrl_T_zmass",tags,zll.mass(),weight);
        if(T_ctrl_side) mon.fillHisto("Ctrl_Tside_RedMet",tags,aRedMet.pt(),weight);
        if(WW_ctrl)     mon.fillHisto("Ctrl_WW_RedMet",tags,aRedMet.pt(),weight);
        if(WW_ctrl_1)   mon.fillHisto("Ctrl_WW_RedMet1",tags,aRedMet.pt(),weight);
        if(WW_ctrl_2)   mon.fillHisto("Ctrl_WW_RedMet2",tags,aRedMet.pt(),weight);
        if(WW_ctrl_3)   mon.fillHisto("Ctrl_WW_RedMet3",tags,aRedMet.pt(),weight);
        if(WW_ctrl_4)   mon.fillHisto("Ctrl_WW_RedMet4",tags,aRedMet.pt(),weight);
        if(WW_ctrl_5)   mon.fillHisto("Ctrl_WW_RedMet5",tags,aRedMet.pt(),weight);
        if(WWCtrl_RJ)   mon.fillHisto("Ctrl_WW_zmass",tags,zll.mass(),weight);

        //Sherpa control
        // N.B. 7 TeV:  Sherpa  Zpt>12,  MG Zpt>50  --> cut Zpt>50 to fill the plot
        //      8 TeV:  Sherpa  Zpt>12,  MG Zpt>12  --> no need to cut
        if( (!use2011Id) || zll.mass()>50. ) {
            mon.fillHisto("npfjets_pres", tags, nAJetsGood30,weight);
        }
        mon.fillHisto("zmass",       tags, zll.mass(), weight);

        //##############################################
        //########  Main Event Selection        ########
        //##############################################
        //event category
        int eventSubCat  = eventCategoryInst.Get(phys,&aGoodIdJets);
        TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);

        if(tag_subcat != "geq2jets") {
            tags.push_back(tag_cat+tag_subcat);
            if(tag_cat=="mumu" || tag_cat=="ee") tags.push_back("ll"+tag_subcat);
        }


        if(passZmass) {
            mon.fillHisto("eventflow",   tags, 3,            weight);

            mon.fillHisto("nvtx"     ,   tags, ev.nvtx,      weight);
            mon.fillHisto("nvtxraw"  ,   tags, ev.nvtx,      1);
            mon.fillHisto("rho"      ,   tags, ev.rho,       weight);
            mon.fillHisto("rho25"    ,   tags, ev.rho25Neut, weight);

            mon.fillHisto("zpt"      ,   tags, zll.pt(),     weight);
            mon.fillHisto("zpt_rebin",   tags, zll.pt(),     weight);
            if(nAJetsGood30==0) mon.fillHisto("zpt_noJet", tags, zll.pt(), weight);
            if(nAJetsGood30==0 && genmatchid[0]>-1 && genmatchid[1]>-1 && genmatchid[0]!=genmatchid[1]) mon.fillHisto("zptGen_noJet", tags, genZP4.pt(), weight);
            mon.fillHisto("zeta"     ,   tags, zll.eta(),    weight);

            if(passZpt) {
                mon.fillHisto("eventflow",tags,4,weight);

                mon.fillHisto("Cos_lep1Z_CS",tags,colin_soper,weight);

                //analyze lepton kinematics
                LorentzVector leadingLep(phys.leptons[0].pt()>phys.leptons[1].pt() ? phys.leptons[0]: phys.leptons[1]);
                LorentzVector trailerLep(phys.leptons[0].pt()>phys.leptons[1].pt() ? phys.leptons[1]: phys.leptons[0]);
                mon.fillHisto("leadeta"     ,   tags, leadingLep.eta()   ,weight);
                mon.fillHisto("leadpt"      ,   tags, leadingLep.pt()    ,weight);
                mon.fillHisto("trailereta"  ,   tags, trailerLep.eta()   ,weight);
                mon.fillHisto("trailerpt"   ,   tags, trailerLep.pt()    ,weight);
                mon.fillHisto("nleptons",tags,2+nextraleptons,weight);
                for(size_t iel=0; iel<extraLeptonsP4.size(); iel++) {
                    mon.fillHisto("thirdleptonpt" ,   tags,extraLeptonsP4[iel].pt()     ,weight);
                    mon.fillHisto("thirdleptoneta",   tags,extraLeptonsP4[iel].eta()   ,weight);
                }

                if(pass3dLeptonVeto) {
                    mon.fillHisto("eventflow",tags,5,weight);

                    //pre-tagged jet control
                    for(size_t ij=0; ij<aGoodIdJets.size(); ij++) {
                        mon.fillHisto("pfjetpt",  tags, aGoodIdJets[ij].pt(),weight);
                        mon.fillHisto("pfjeteta",  tags, fabs(aGoodIdJets[ij].eta()),weight);
                    }

                    //final jet control
                    mon.fillHisto("npfjets",              tags, nAJetsGood30,weight);
                    mon.fillHisto("npfjetsvspu",          tags, ev.ngenITpu, nAJetsGood30,weight);

                    // jet veto
                    if(passJetveto) { //set passJetveto always true

                        mon.fillHisto("eventflow",tags,6,weight);

                        //pre-tagged jet control
                        for(size_t ij=0; ij<aGoodIdJets.size(); ij++) {
                            if(aGoodIdJets[ij].pt()<30 || fabs(aGoodIdJets[ij].eta())>2.5) continue;
                            if(isMC) {
                                if(fabs(aGoodIdJets[ij].flavid)==5) mon.fillHisto("bpfjetstags",     tags, aGoodIdJets[ij].btag2, weight);
                                else                                mon.fillHisto("otherpfjetstags", tags, aGoodIdJets[ij].btag2, weight);
                            }
                            mon.fillHisto("pfjetstags",  tags, aGoodIdJets[ij].btag2,weight);
                        }
                        mon.fillHisto("npfjetsbtags",  tags, nABtags ,weight);

                        //b-veto
                        if(passBveto) {
                            mon.fillHisto("eventflow",tags,/*6*/ 7,weight);
                            mon.fillHisto("eventflow_gamma",tags,0,weight);

                            //Data-Driven DY template //RJ
                            mon.fillHisto("qt", tags, zll.pt(), weight, true);
                            mon.fillHisto("qmass",tags, zll.mass(),weight);
                            mon.fillHisto("qt_rebin", tags, zll.pt(), weight, true);
                            mon.fillHisto("nvtx_dy"     ,   tags, ev.nvtx,      weight);
                            mon.fillHisto("met_redMet_phi_dy",tags,aRedMet.phi(),weight);
                            mon.fillHisto("met_met_phi_dy",   tags,zvvs[0].phi(),weight);
                            mon.fillHisto("qdphi2lep",  tags, dphi2l, weight);
                            mon.fillHisto("qCosl1Z_CS",  tags, colin_soper, weight);

                            //mon.fillHisto("nvtxVSqtrebin",tags,zll.pt(),ev.nvtx,weight,true);
                            //mon.fillHisto("nvtxVSqt",tags,zll.pt(),ev.nvtx,weight,true);

                            mon.fillHisto("zpt_dy"      ,   tags, zll.pt(),     weight);
                            mon.fillHisto("npfjets_dy",              tags, nAJetsGood30,weight);

                            if(zvvs[0].pt()>0)  mon.fillHisto("mindphijmet_0",  tags, nAJetsGood30==0 ? mindphijmet15:mindphijmet,weight);
                            if(zvvs[0].pt()>25) mon.fillHisto("mindphijmet_25", tags, nAJetsGood30==0 ? mindphijmet15:mindphijmet,weight);
                            if(zvvs[0].pt()>50) mon.fillHisto("mindphijmet_50", tags, nAJetsGood30==0 ? mindphijmet15:mindphijmet,weight);
                            if(zvvs[0].pt()>50) mon.fillHisto("mindphijmet",    tags, nAJetsGood30==0 ? mindphijmet15:mindphijmet,weight);
                            if(zvvs[0].pt()>50 && passLMetVeto) mon.fillHisto("mindphijmet_aftPhil",    tags, nAJetsGood30==0 ? mindphijmet15:mindphijmet,weight);
                            mon.fillHisto("mindphilmet",tags, min(dphil1met,dphil2met) ,weight);
                            if(passDphijmet) mon.fillHisto("mindphilmet_aftPhij",tags, min(dphil1met,dphil2met) ,weight);
                            if(passDphijmet && zvvs[0].pt()>60 ) mon.fillHisto("mindphilmet_aftPhijPt",tags, min(dphil1met,dphil2met) ,weight);
                            mon.fillHisto("maxdphilmet",tags, max(dphil1met,dphil2met) ,weight);
                            mon.fillHisto("met_metbeforedphilmet",         tags,  zvvs[0].pt(),  weight);
                            mon.fillHisto("met_mindphilmet",tags,zvvs[0].pt(),min(dphil1met,dphil2met),weight);


                            if(passDphijmet && passLMetVeto) {
                                mon.fillHisto("eventflow",tags, 8, weight);
                                mon.fillHisto("eventflow_gamma",tags,1,weight);

                                //RJ, Data-Driven DY sample
                                mon.fillHisto("met_met",         tags,  zvvs[0].pt(),  weight);
                                mon.fillHisto("met_metL",    tags,  axialMet,  weight);
                                mon.fillHisto("met_mvamet",      tags,  mvaMetP4.pt(), weight);
                                mon.fillHisto("met_typeImet",    tags,  fullTypeIMetP4.pt(),  weight);
                                mon.fillHisto("met_redMet",tags,aRedMet.pt(),weight);
                                mon.fillHisto("met_redMet_rebin",tags,aRedMet.pt(),weight);
                                //mon.fillHisto("met_redMet_phi_dy",tags,aRedMet.phi(),weight);
                                //mon.fillHisto("met_met_phi_dy",	  tags,zvvs[0].phi(),weight);
                                mon.fillHisto("met_redMetL",tags,aRedMetT,weight);
                                mon.fillHisto("met_redMetT",tags,aRedMetL,weight);
                                mon.fillHisto("met_clusteredMetL",tags,clusteredMetL,weight);
                                mon.fillHisto("met_clusteredMetT",tags,clusteredMetT,weight);
                                mon.fillHisto("met_metCompL",tags,pfMetCompL,weight);
                                mon.fillHisto("met_metCompT",tags,pfMetCompT,weight);
                                mon.fillHisto("balance",tags, zvvs[0].pt()/zll.pt(),weight);
                                mon.fillHisto("dphi_Zll_met_met",tags, dphiZllmet, weight);
                                mon.fillHisto("mt",tags, aMT, weight);

                                mon.fillHisto("met_met_ROC",tags,  zvvs[0].pt(),  weight);
                                mon.fillHisto("met_redMet_ROC",tags,aRedMet.pt(), weight);

                                //clusteredMetL_vsPhi
                                mon.fillHisto("clusteredMetL_vsPhi",tags, zll.phi(),clusteredMetL,weight);
                                mon.fillHisto("clusteredMetT_vsPhi",tags, zll.phi(),clusteredMetT,weight);
                                mon.fillHisto("RecoilL",tags,RecoilMetL,weight);
                                mon.fillHisto("RecoilT",tags,RecoilMetT,weight);
                                mon.fillHisto("RecoilMetL_vsPhi",tags, zll.phi(),RecoilMetL,weight);
                                mon.fillHisto("RecoilMetT_vsPhi",tags, zll.phi(),RecoilMetT,weight);
                                mon.fillHisto("RecoilT_vsqT",tags,zll.pt(),RecoilMetT,weight);
                                mon.fillHisto("RecoilL_vsqT",tags,zll.pt(),RecoilMetL,weight);



                                if(passdphiZllmetCut) {
                                    mon.fillHisto("deltaleptonpt",tags, leadingLep.pt()-trailerLep.pt()    ,weight);
                                    mon.fillHisto("deltazpt",tags, zll.pt()-zvvs[0].pt(),weight);
                                    mon.fillHisto("balancedif",tags, balanceDif,weight);

                                    if(passRedMet) {
                                        mon.fillHisto("mindphilmet_aftMet",tags, min(dphil1met,dphil2met) ,weight);
                                        mon.fillHisto("eventflow",tags,/*8*/ 9,weight);
                                        mon.fillHisto("eventflow_gamma",tags,2,weight);

                                        if(passBalanceCut) {   //RJ
                                            mon.fillHisto("eventflow",tags,10,weight);
                                            mon.fillHisto("eventflow_gamma",tags,3,weight);

                                            if(passMTcut) { //RJ
                                                mon.fillHisto("eventflow",tags,11,weight);
                                                mon.fillHisto("eventflow_gamma",tags,4,weight);

                                                // Final distributions
                                                mon.fillHisto("mt_final",          tags, aMT,          weight);
                                                mon.fillHisto("new_mt_final",	   tags, aMTmassless,  weight);
                                                mon.fillHisto("zpt_final",         tags, zll.pt(),     weight);
                                                mon.fillHisto("zpt_rebin_final",   tags, zll.pt(),     weight);
                                                mon.fillHisto("met_met_final",     tags, zvvs[0].pt(), weight);
                                                mon.fillHisto("met_redMet_final",  tags, aRedMet.pt(), weight);
                                                mon.fillHisto("met_redMet_rebin_final", tags, aRedMet.pt(), weight);
                                                mon.fillHisto("met_redMetL_final", tags, aRedMetL,     weight);

                                                //RJ
                                                mon.fillHisto("dphi_Zll_met_met_final",tags, dphiZllmet, weight);
                                                mon.fillHisto("dphi_Zll_met_redMet_final",tags, dphiZllredMet, weight);
                                                mon.fillHisto("balancevsmt",tags,zvvs[0].pt()/zll.pt(),aMT,weight);

                                                mon.fillHisto("Cos_lep1Z_CS_final",tags,colin_soper,weight);

						//mon.fillHisto("dphiLLvsMT",tags,aMTmassless,dphi2l,weight);

                                                /*
                                                TString evCat(tag_cat);
                                                	{
                                                		int eventSubCat  = eventCategoryInst.Get(phys,&aGoodIdJets);
                                                		TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);
                                                		evCat+=tag_subcat;
                                                	}
                                                fprintf(outTxtFile,"%s %d %d %d\n",evCat.Data(),ev.run,ev.lumi,ev.event);
                                                */

                                                //Print out event informaiton
                                                //cout << "--------------------------"<< endl;
                                                //cout << "ev.cat: " << ev.cat << endl;
                                                //cout << "ev.run: " << ev.run << endl;
                                                //cout << "ev.lumi: " << ev.lumi << endl;
                                                //cout << "ev.event: " << ev.event << endl;
                                                //cout << "--------------------------"<< endl;


                                            } //end MT cut
                                        }//end passBalanceCut
                                    } // end passRedMet
                                }//end passdphiZllmetCut
                            }//end passLMetVeto,
                        }//end passBveto
                    }//end passJetveto
                }//3lept
            }//end passZpt
        }//end passZmass

        if(passMTcut 			&& passZmass && passZpt && pass3dLeptonVeto && passJetveto && passBveto && passDphijmet && passBalanceCut && passRedMet) {
            mon.fillHisto("dphi_Zll_met_metNM1", tags, dphiZllmet,weight);
        }
        if(             passdphiZllmetCut && passZmass && passZpt && pass3dLeptonVeto && passJetveto && passBveto && passDphijmet && passBalanceCut && passRedMet) {
            mon.fillHisto("mtNM1", tags, aMT ,weight);
        }
        if(passMTcut && passdphiZllmetCut &&             passZpt && pass3dLeptonVeto && passJetveto && passBveto && passDphijmet && passBalanceCut && passRedMet) {
            mon.fillHisto("zmassNM1",tags,zll.mass(),weight);
        }
        if(passMTcut && passdphiZllmetCut && passZmass            && pass3dLeptonVeto && passJetveto && passBveto && passDphijmet && passBalanceCut && passRedMet) {
            mon.fillHisto("zptNM1",tags,zll.pt(),weight);
        }
        if(passMTcut && passdphiZllmetCut && passZmass && passZpt                     && passJetveto && passBveto && passDphijmet && passBalanceCut && passRedMet) {
            mon.fillHisto("nleptonsNM1",tags,2+nextraleptons,weight);
        }
        if(passMTcut && passdphiZllmetCut && passZmass && passZpt && pass3dLeptonVeto                && passBveto && passDphijmet && passBalanceCut && passRedMet) {
            mon.fillHisto("npfjetsNM1", tags,nAJetsGood30,weight);
        }
        if(passMTcut && passdphiZllmetCut && passZmass && passZpt && pass3dLeptonVeto && passJetveto              && passDphijmet && passBalanceCut && passRedMet) {
            mon.fillHisto("npfjetsbtagsNM1", tags, nABtags ,weight);
        }
        if(passMTcut && passdphiZllmetCut && passZmass && passZpt && pass3dLeptonVeto && passJetveto && passBveto                 && passBalanceCut && passRedMet) {
            mon.fillHisto("mindphijmetNM1", tags, nAJetsGood30==0 ? mindphijmet15:mindphijmet ,weight);
        }
        if(passMTcut && passdphiZllmetCut && passZmass && passZpt && pass3dLeptonVeto && passJetveto && passBveto && passDphijmet                   && passRedMet) {
            mon.fillHisto("balanceNM1",tags,zvvs[0].pt()/zll.pt(),weight);
        }
        if(passMTcut && passdphiZllmetCut && passZmass && passZpt && pass3dLeptonVeto && passJetveto && passBveto && passDphijmet && passBalanceCut              ) {
            mon.fillHisto("met_metNM1", tags, zvvs[0].pt() ,weight);
            mon.fillHisto("met_redMetNM1", tags, aRedMet.pt(),weight);
            //RJ
            mon.fillHisto("met_met_ROC2", 	 tags, zvvs[0].pt(),weight);
            mon.fillHisto("met_redMet_ROC2", tags, aRedMet.pt(),weight);
        }

        //RJ. Looser preSelection.
        bool passZmassPRS = (fabs(zll.mass()-91)<25);
        bool passZptPRS = (zll.pt()>30);
        bool passRedMetPRS = (aRedMet.pt()>50);

        if(passZmassPRS && passZptPRS && passRedMetPRS) {
            mon.fillHisto("dphi_Zll_met_metPRS", tags, dphiZllmet,weight);
            mon.fillHisto("mtPRS", tags, aMT ,weight);
            mon.fillHisto("zmassPRS",tags,zll.mass(),weight);
            mon.fillHisto("zptPRS",tags,zll.pt(),weight);
            mon.fillHisto("nleptonsPRS",tags,2+nextraleptons,weight);
            mon.fillHisto("npfjetsPRS", tags,nAJetsGood30,weight);
            mon.fillHisto("npfjetsbtagsPRS", tags, nABtags ,weight);
            mon.fillHisto("balancePRS",tags,zvvs[0].pt()/zll.pt(),weight);
            mon.fillHisto("met_metPRS", tags, zvvs[0].pt() ,weight);
            mon.fillHisto("met_redMetPRS", tags, aRedMet.pt(),weight);
        }


        //RENJIE, Top Control, MET control in the sideband
        bool passSB( ((zll.mass()>40 && zll.mass()<70) || (zll.mass()>110 && zll.mass()<200)) && zll.pt()>55 );
        if(passSB && pass3dLeptonVeto && passDphijmet && !passBveto && passLMetVeto) {
            mon.fillHisto("met_metSB",tags,zvvs[0].pt(),weight);
            mon.fillHisto("met_redMetSB",tags,aRedMet.pt(),weight);
        }

        //
        // HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
        //
        //Fill histogram for posterior optimization, or for control regions
        for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
            float iweight = weight;                                               //nominal
            if(ivar==9)                         iweight *=TotalWeight_plus;        //pu up
            if(ivar==10)                        iweight *=TotalWeight_minus;       //pu down
            float Sherpa_weight(1.);
            if(use2011Id && (ivar==13 || ivar==14)) {
                if( zll.pt() < 84) Sherpa_weight = 0.976138;
                else if( zll.pt()< 98) Sherpa_weight = 1.01469;
                else if( zll.pt()< 112) Sherpa_weight = 1.00612;
                else if( zll.pt()< 126) Sherpa_weight = 0.981601;
                else if( zll.pt()< 140) Sherpa_weight = 0.969641;
                else if( zll.pt()< 154) Sherpa_weight = 0.980155;
                else if( zll.pt()< 168) Sherpa_weight = 1.00758;
                else if( zll.pt()< 182) Sherpa_weight = 1.04195;
                else if( zll.pt()< 196) Sherpa_weight = 1.07586;
                else if( zll.pt()< 210) Sherpa_weight = 1.10551;
                else if( zll.pt()< 224) Sherpa_weight = 1.12925;
                else if( zll.pt()< 238) Sherpa_weight = 1.14636;
                else if( zll.pt()< 252) Sherpa_weight = 1.15645;
                else if( zll.pt()< 266) Sherpa_weight = 1.15935;
                else if( zll.pt()< 280) Sherpa_weight = 1.15498;
                else if( zll.pt()< 294) Sherpa_weight = 1.14343;
                else if( zll.pt()< 308) Sherpa_weight = 1.12491;
                else if( zll.pt()< 322) Sherpa_weight = 1.09977;
                else if( zll.pt()< 336) Sherpa_weight = 1.06847;
                else                    Sherpa_weight = 1.03156;
                if(ivar==20) Sherpa_weight = (1./Sherpa_weight);
            }
            if(!use2011Id && (ivar==13 || ivar==14)) {
                if( zll.pt() < 84) Sherpa_weight = 1.0724;
                else if( zll.pt()< 112) weight = 0.981944;
                else if( zll.pt()< 126) weight = 0.945408;
                else if( zll.pt()< 140) weight = 0.947075;
                else if( zll.pt()< 154) weight = 0.974582;
                else if( zll.pt()< 168) weight = 1.00435;
                else if( zll.pt()< 182) weight = 1.01752;
                else if( zll.pt()< 196) weight = 1.00596;
                else if( zll.pt()< 210) weight = 0.970417;
                else if( zll.pt()< 224) weight = 0.915903;
                else if( zll.pt()< 238) weight = 0.848171;
                else if( zll.pt()< 252) weight = 0.772162;
                else if( zll.pt()< 266) weight = 0.691807;
                else if( zll.pt()< 280) weight = 0.610271;
                else if( zll.pt()< 294) weight = 0.530154;
                else if( zll.pt()< 308) weight = 0.453577;
                else if( zll.pt()< 322) weight = 0.382189;
                else if( zll.pt()< 336) weight = 0.317164;
                else if( zll.pt()< 350) weight = 0.25922;
            }

            //recompute MET/MT if JES/JER was varied
            LorentzVector zvv = zvvs[ivar>8 ? 0 : ivar];
            PhysicsObjectJetCollection &varJets = ( ivar<=4 ? variedAJets[ivar] : ( useJERsmearing ? variedAJets[0] : recoJets ) );
            PhysicsObjectJetCollection tightVarJets;
            LorentzVector clusteredMetP4(zll);
            clusteredMetP4 *= -1;
            bool passLocalJetveto(true);
            bool passLocalBveto(true); //bool passLocalBveto(passBveto);
            bool passLocalDphijmet(true);
            bool passLocalBalanceCut(true);
            float localmindphijmet(999999.),localmindphijmet15(999999.);
            int localNAJetsGood30(0);
            int localNAJetsGood15(0);
            int localNAJetsGood20(0);
            int localNAJetsGood25(0);
            for(size_t ijet=0; ijet<varJets.size(); ijet++) {
                //if(varJets[ijet].pt()<15) continue;
                if(varJets[ijet].pt()<15){
                          continue;
                        }
                        else{
                          localNAJetsGood15++;
                        }

                // dphi
                float idphijmet( fabs(deltaPhi(varJets[ijet].phi(),zvv.phi()) ) );
                if(varJets[ijet].pt()>15) if(idphijmet<localmindphijmet15) localmindphijmet15 = idphijmet;
                if(varJets[ijet].pt()>30) if(idphijmet<localmindphijmet)   localmindphijmet   = idphijmet;

                bool isGoodJet = hasObjectId(aJets[ijet].pid,JETID_LOOSE);
                if(usePUsubJetId) isGoodJet = hasObjectId(aJets[ijet].pid,JETID_CUTBASED_LOOSE);
                if(!isGoodJet) continue;

                clusteredMetP4 -= varJets[ijet];

                if(useJetsOnlyInTracker && fabs(varJets[ijet].eta())>2.5) continue;

                tightVarJets.push_back( varJets[ijet] );
                if(varJets[ijet].pt()>30)localNAJetsGood30++;
                if(varJets[ijet].pt()>25)localNAJetsGood25++;
                if(varJets[ijet].pt()>20)localNAJetsGood20++;

                if(varJets[ijet].pt()>20 && fabs(varJets[ijet].eta())<2.5) {
                    if(ivar==11)      passLocalBveto &= (varJets[ijet].btag2<0.250);
                    else if(ivar==12) passLocalBveto &= (varJets[ijet].btag2<0.240);
                    else              passLocalBveto &= (varJets[ijet].btag2<0.244);
                }
            }
            passLocalJetveto=(localNAJetsGood30==0);
            //passLocalDphijmet=(localmindphijmet>0.5); //RJ
            //if(localNAJetsGood30==0) passLocalDphijmet=(localmindphijmet15>0.5);
            passLocalBalanceCut=(zvv.pt()/zll.pt()>0.8 && zvv.pt()/zll.pt()<1.2);

            double dphil1met=fabs(deltaPhi(lep1.phi(),zvv.phi()));
            double dphil2met=fabs(deltaPhi(lep2.phi(),zvv.phi()));
              double dphi12 = fabs(deltaPhi(lep1.phi(),lep2.phi()));
              double dphiZmet = fabs(deltaPhi((lep1+lep2).phi(),zvv.phi()));
              double dphiZH = fabs(deltaPhi(zll.phi(),zvv.phi()));
              //double MTZH0 = sqrt(2*zll.pt()*zvv.pt()*(1-cos(dphiZH)));
              //double MTZH = (zll + zvv).Mt();
            bool passLocalLMetVeto(true);
            if(!use2011Id && zvv.pt()>60 && min(dphil1met,dphil2met)<0.2) passLMetVeto=false;

            float mt = METUtils::transverseMass(zll,zvv,true);
            LorentzVector nullP4(0,0,0,0);
            LorentzVector redMet = METUtils::redMET(METUtils::INDEPENDENTLYMINIMIZED, zll, 0, nullP4, 0, clusteredMetP4, zvv, true);

            float new_mt = METUtils::transverseMass(zll,zvv,false); //RJ

            double dphiredMetMet=fabs(deltaPhi(redMet.phi(),zvv.phi()));
  ///////////////////////////////////// TREE VARIABLES

      //bool passFullSelection (isSameFlavor && passZmass && passZpt && passLocalBveto && passLocalJetveto && pass3dLeptonVeto && passLocalDphijmet && passLocalBalanceCut && (zvv.pt() > 50.0));// && passRedMet);//(redMet.pt() > 70));
      //bool passFullSelection (isSameFlavor && passZmass15 && passZpt && passLocalBveto && passLocalJetveto && pass3dLeptonVeto && passLocalBalanceCut && (zvv.pt() > 50.0));// && passRedMet);//(redMet.pt() > 70));
      bool passFullSelection (isSameFlavor && passZmass && passLocalJetveto && pass3dLeptonVeto && (redMet.pt() > 90) && passLocalBalanceCut);// && (zll.pt() > 50.0));// && passRedMet);//(redMet.pt() > 70));

      bool passMT135 (new_mt > 135.0);
      bool passMT120 (new_mt > 120.0);
      bool passMT105 (new_mt > 105.0);
      bool passMT80 (new_mt > 80.0);
      bool passET65 (zvv.pt() > 65.0);
      bool passET80 (zvv.pt() > 80.0);
      bool passET95 (zvv.pt() > 95.0);

      pZmass = 0;
      pZpt = 0;
      pBveto = 0;
      pLepVeto = 0;
      pDphijmet = 0;
      pBalance = 0;
      predMet = 0;
      pJetVeto = 0;

      pMT80 = 0;
      pMT105 = 0;
      pMT120 = 0;
      pMT135 = 0;
      pET65 = 0;
      pET80 = 0;
      //metdiff = fabs(redMet.pt()-redMet_cor.pt());

      // if (passZmass) pZmass = 1;
      // if (passZpt) pZpt = 1;
      // if (passLocalBveto)  pBveto = 1;
      // if (pass3dLeptonVeto)  pLepVeto = 1;
      // if (passDphijmet) pDphijmet = 1;
      // if (passLocalBalanceCut) pBalance = 1;
      // if (redMet.pt() > 70) predMet = 1;

      if (passZmass)    pZmass = 1;
      if (passZpt)  pZpt = 1;
      if (passLocalBveto)   pBveto = 1;
      if (passLocalJetveto) pJetVeto = 1;
      if (pass3dLeptonVeto) pLepVeto = 1;
      if (passLocalDphijmet) pDphijmet = 1;
      if (passLocalBalanceCut) pBalance = 1;
      if (passRedMet) predMet = 1;

      if(passMT135) pMT135 = 1;
      if(passMT120) pMT120 = 1;
      if(passMT105) pMT105 = 1;
      if(passMT80) pMT80 = 1;
      //if(passET95) metdiff = 1;
      if(passET80) pET80 = 1;
      if(passET65) pET65 = 1;

      //bool passFullSelection (isSameFlavor && passZmass && passZpt);

      Eweight = 0;
      mass = 0;
      zpt = 0;
      zeta = 0;
      zphi = 0;
      l1pt = 0;
      l2pt = 0;
      l1eta = 0;
      l2eta = 0;
      l1phi = 0;
      l2phi = 0;
      met = 0;
      metphi = 0;
      l1Err = 0;
      l2Err = 0;
      l1id = 0;
      l2id = 0;
      REDmet = 0;
      Zmetphi = 0;
      if(passFullSelection)
      {
        finstate = Double_t(ev.cat);
        Eweight = iweight;
        mass = zll.mass();
        zpt = zll.pt();
        zeta = zll.eta();
        zphi = zll.phi();
        l1pt = lep1.pt();
        l2pt = lep2.pt();
        l1eta = lep1.eta();
        l2eta = lep2.eta();
        l1phi = lep1.phi();
        l2phi = lep2.phi();
        met = zvv.pt();  //zvv NOT zvvs?
        metphi = zvv.phi();
        l1Err = phys.leptons[0].ptErr;
        l2Err = phys.leptons[1].ptErr;
        l1id = phys.leptons[0].id;
        l2id = phys.leptons[1].id;
        REDmet = redMet.pt();
        REDmetmetphi = dphiredMetMet;
        Zmetphi = dphiZmet;
        llphi = dphi12;
        nj15 = localNAJetsGood15;
        nj30 = localNAJetsGood30;
        nj25 = localNAJetsGood25;
        nj20 = localNAJetsGood20;
        mtzh = new_mt;
        mtzh3 = new3rdMT;
        baldiff = balanceDif;
        ColinSoper = colin_soper;
        AxialMet = axialMet;
        ZREDmetphi = dphiZllredMet;
        //nj10 = localNAJetsGood10;


        if (ivar == 0) stmkr->Fill();
        if (runSystematics) {
        if (ivar == 1) stmkr_jerup->Fill();
        if (ivar == 2) stmkr_jerdown->Fill();
        if (ivar == 3) stmkr_jesup->Fill();
        if (ivar == 4) stmkr_jesdown->Fill();
        if (ivar == 5) stmkr_umetup->Fill();
        if (ivar == 6) stmkr_umetdown->Fill();
        if (ivar == 7) stmkr_lesup->Fill();
        if (ivar == 8) stmkr_lesdown->Fill();
        if (ivar == 9) stmkr_puup->Fill();
        if (ivar == 10) stmkr_pudown->Fill();
        if (ivar == 11) stmkr_btagup->Fill();
        if (ivar == 12) stmkr_btagdown->Fill();


/*                                 "_jerup","_jerdown",
                         "_jesup","_jesdown",
                         "_umetup","_umetdown",
                         "_lesup","_lesdown",
                         "_puup","_pudown",
                         "_btagup","_btagdown",
                         "_sherpaup","_sherpadown"
*/
        //ivar+1, varNames[ivar]

        }
    }


            // with standard Z mass, Z pt, RedMet (will be variated later)
            bool passPreselection(passZmass && passZpt && pass3dLeptonVeto && passLocalJetveto && passLocalBveto && passLocalDphijmet && passLocalLMetVeto && passLocalBalanceCut && passRedMet);
            //bool passPreselectionMbvetoMzmass(passZpt && pass3dLeptonVeto && passLocalJetveto && passLocalDphijmet && passLocalLMetVeto && passLocalBalanceCut && passRedMet);

            //RJ compute dphi(zpt,redMet) for statistics
            double LocaldphiZllmet=fabs(deltaPhi(zll.phi(),zvv.phi()));
            //bool passLocaldphiZllmetCut(LocaldphiZllmet>2.6);

            //RJ
            bool passLocalMTcut(mt>220 && mt<1200);
	    //passLocalMTcut &= passdphi2l;

            //re-assign the event category if jets were varied
            int eventSubCat  = eventCategoryInst.Get(phys,&tightVarJets);
            TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);
            tags.clear();

            if(tag_subcat != "geq2jets") {
                //tags.push_back(tag_cat);
                tags.push_back(tag_cat+tag_subcat);
                if(tag_cat=="mumu" || tag_cat=="ee") tags.push_back(string("ll")+tag_subcat);
            }

            //fill shapes
            for(unsigned int index=0; index<optim_Cuts1_met.size(); index++) {

                float minMet=optim_Cuts1_met[index];
                float minBalance=optim_Cuts1_balance[index];
                float minDphi=optim_Cuts1_dphi[index];
                float deltaZ=optim_Cuts1_zmass[index];

                bool passLocalRedMet(redMet.pt()>minMet); //Reduced MET: redMet.pt() ; PFMET: zvv.pt()
                bool passLocalRedMet_RJ(redMet.pt()>65);  //Reduced MET: redMet.pt() ; PFMET: zvv.pt()  for passPreselectionMjvetoMbvetoMzmass
                bool passLocalZmass(fabs(zll.mass()-91)<deltaZ);

                bool passLocalZpt(zll.pt()>50.); //fix Zpt cut
                passLocalBalanceCut=(zvv.pt()/zll.pt()>(1.-minBalance) && zvv.pt()/zll.pt()<(1.+minBalance));
                bool passLocalBalanceCut_RJ=(zvv.pt()/zll.pt()>0.4 && zvv.pt()/zll.pt()<1.8);
                bool passLocaldphiZllmetCut(LocaldphiZllmet>minDphi);

                passPreselection = (passLocalMTcut && passLocaldphiZllmetCut && passLocalZmass && passLocalZpt && pass3dLeptonVeto && passLocalJetveto && passLocalBveto && passLocalDphijmet && passLocalBalanceCut && passLocalRedMet);
                //passPreselectionMbvetoMzmass = (passLocalMTcut && passLocaldphiZllmetCut && passLocalZpt && pass3dLeptonVeto && passLocalJetveto && passLocalDphijmet && passLocalBalanceCut && passLocalRedMet); //not use
                bool passPreselectionMjvetoMbvetoMzmass = ( passLocalRedMet_RJ && passLocalZpt && pass3dLeptonVeto && passLocalDphijmet && passLocalBalanceCut_RJ );
                //bool passPreselectionMjvetoMbvetoMzmass = (passLocalRedMet && passLocalZpt && pass3dLeptonVeto && passLocalDphijmet && passLocalBalanceCut && passLocaldphiZllmetCut);

                if( passPreselection ) {
                    mon.fillHisto(TString("redMet_shapes")+varNames[ivar],tags,index, redMet.pt(),iweight);
                    mon.fillHisto(TString("redMet_rebin_shapes")+varNames[ivar],tags,index, redMet.pt(),iweight);
                    //RJ
                    mon.fillHisto(TString("mt_shapes")+varNames[ivar],tags,index,mt,iweight);
                    mon.fillHisto(TString("new_mt_shapes")+varNames[ivar],tags,index,new_mt,iweight);

                    /* balance_mt_shapes */
                    float i_balance = zvv.pt()/zll.pt();
                    if(i_balance>=BalanceXaxis[0] && i_balance< BalanceXaxis[1]) mon.fillHisto(TString("balance_mt_shapes")+varNames[ivar],tags,index,new_mt,iweight);
                    if(i_balance>=BalanceXaxis[1] && i_balance< BalanceXaxis[2]) mon.fillHisto(TString("balance_mt_shapes")+varNames[ivar],tags,index,new_mt+1200.,iweight);
                    if(i_balance>=BalanceXaxis[2] && i_balance<=BalanceXaxis[3]) mon.fillHisto(TString("balance_mt_shapes")+varNames[ivar],tags,index,new_mt+2400.,iweight);

                    /* dphiLL_mt_shapes */
                    if(dphi2l>=DphiLLXaxis[0] && dphi2l< DphiLLXaxis[1]) mon.fillHisto(TString("dphiLL_mt_shapes")+varNames[ivar],tags,index,new_mt,iweight);
                    if(dphi2l>=DphiLLXaxis[1] && dphi2l< DphiLLXaxis[2]) mon.fillHisto(TString("dphiLL_mt_shapes")+varNames[ivar],tags,index,new_mt+1200.,iweight);
                    if(dphi2l>=DphiLLXaxis[2] && dphi2l<=DphiLLXaxis[3]) mon.fillHisto(TString("dphiLL_mt_shapes")+varNames[ivar],tags,index,new_mt+2400.,iweight);

                    /* coslZ_mt_shapes */
                    if(colin_soper>=CosThetaXaxis[0] && colin_soper< CosThetaXaxis[1]) mon.fillHisto(TString("coslZ_mt_shapes")+varNames[ivar],tags,index,new_mt,iweight);
                    if(colin_soper>=CosThetaXaxis[1] && colin_soper< CosThetaXaxis[2]) mon.fillHisto(TString("coslZ_mt_shapes")+varNames[ivar],tags,index,new_mt+1200.,iweight);
                    if(colin_soper>=CosThetaXaxis[2] && colin_soper<=CosThetaXaxis[3]) mon.fillHisto(TString("coslZ_mt_shapes")+varNames[ivar],tags,index,new_mt+2400.,iweight);

                    /* control plots */
                    mon.fillHisto(TString("balance1Dshape")+varNames[ivar],tags,index,i_balance,iweight);
                    mon.fillHisto(TString("dphiLL1Dshape")+varNames[ivar],tags,index,dphi2l,iweight);
                    mon.fillHisto(TString("coslZ1Dshape")+varNames[ivar],tags,index,colin_soper,iweight);


                    mon.fillHisto(TString("dphi_shapes")+varNames[ivar],tags,index,LocaldphiZllmet,iweight);
                    mon.fillHisto(TString("met_shapes")+varNames[ivar],tags,index,zvv.pt(),iweight);
                    mon.fillHisto(TString("met_rebin_shapes")+varNames[ivar],tags,index,zvv.pt(),iweight);
                    mon.fillHisto(TString("zpt_shapes")+varNames[ivar],tags,index, zll.pt(),iweight);
                    mon.fillHisto(TString("zpt_rebin_shapes")+varNames[ivar],tags,index, zll.pt(),iweight);
                }
                if( passPreselectionMjvetoMbvetoMzmass && passLocalZmass && passLocalJetveto && passLocalBveto ) {
                    mon.fillHisto("redMet_shapes_NRBctrl"+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto("redMet_rebin_shapes_NRBctrl"+varNames[ivar],tags,index,0,iweight);
                    //RJ
                    mon.fillHisto(TString("mt_shapes_NRBctrl")+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto(TString("new_mt_shapes_NRBctrl")+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto(TString("balance_mt_shapes_NRBctrl")+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto(TString("dphiLL_mt_shapes_NRBctrl")+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto(TString("coslZ_mt_shapes_NRBctrl")+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto(TString("balance1Dshape_NRBctrl")+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto(TString("dphiLL1Dshape_NRBctrl")+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto(TString("coslZ1Dshape_NRBctrl")+varNames[ivar],tags,index,0,iweight);

                    mon.fillHisto(TString("dphi_shapes_NRBctrl")+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto(TString("met_shapes_NRBctrl")+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto(TString("met_rebin_shapes_NRBctrl")+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto("zpt_shapes_NRBctrl"+varNames[ivar],tags,index,0,iweight);
                    mon.fillHisto("zpt_rebin_shapes_NRBctrl"+varNames[ivar],tags,index,0,iweight);
                }
                if( passPreselectionMjvetoMbvetoMzmass && isZsideBand && passLocalJetveto && passLocalBveto ) {
                    mon.fillHisto("redMet_shapes_NRBctrl"+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto("redMet_rebin_shapes_NRBctrl"+varNames[ivar],tags,index,1,iweight);
                    //RJ
                    mon.fillHisto(TString("mt_shapes_NRBctrl")+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto(TString("new_mt_shapes_NRBctrl")+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto(TString("balance_mt_shapes_NRBctrl")+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto(TString("dphiLL_mt_shapes_NRBctrl")+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto(TString("coslZ_mt_shapes_NRBctrl")+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto(TString("balance1Dshape_NRBctrl")+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto(TString("dphiLL1Dshape_NRBctrl")+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto(TString("coslZ1Dshape_NRBctrl")+varNames[ivar],tags,index,1,iweight);

                    mon.fillHisto(TString("dphi_shapes_NRBctrl")+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto(TString("met_shapes_NRBctrl")+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto(TString("met_rebin_shapes_NRBctrl")+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto("zpt_shapes_NRBctrl"+varNames[ivar],tags,index,1,iweight);
                    mon.fillHisto("zpt_rebin_shapes_NRBctrl"+varNames[ivar],tags,index,1,iweight);
                }
                if( passPreselectionMjvetoMbvetoMzmass && isZsideBandPlus && passLocalJetveto && passLocalBveto ) {
                    mon.fillHisto("redMet_shapes_NRBctrl"+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto("redMet_rebin_shapes_NRBctrl"+varNames[ivar],tags,index,2,iweight);
                    //RJ
                    mon.fillHisto(TString("mt_shapes_NRBctrl")+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto(TString("new_mt_shapes_NRBctrl")+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto(TString("balance_mt_shapes_NRBctrl")+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto(TString("dphiLL_mt_shapes_NRBctrl")+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto(TString("coslZ_mt_shapes_NRBctrl")+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto(TString("balance1Dshape_NRBctrl")+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto(TString("dphiLL1Dshape_NRBctrl")+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto(TString("coslZ1Dshape_NRBctrl")+varNames[ivar],tags,index,2,iweight);

                    mon.fillHisto(TString("dphi_shapes_NRBctrl")+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto(TString("met_shapes_NRBctrl")+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto(TString("met_rebin_shapes_NRBctrl")+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto("zpt_shapes_NRBctrl"+varNames[ivar],tags,index,2,iweight);
                    mon.fillHisto("zpt_rebin_shapes_NRBctrl"+varNames[ivar],tags,index,2,iweight);
                }


                //RJ, changing eventCategory for NRB control sample, i.e., no jetveto
                //prepare the tag's vectors for histo filling
                std::vector<TString> NRBtags(1,"all");
                NRBtags.clear();
                //NRBtags.push_back(tag_cat);
                NRBtags.push_back(tag_cat+"eq0jets");
                NRBtags.push_back(tag_cat+"eq1jets");
                if(tag_cat=="mumu" || tag_cat=="ee") {
                    NRBtags.push_back(string("ll")+"eq0jets");
                    NRBtags.push_back(string("ll")+"eq1jets");
                }

                if( passPreselectionMjvetoMbvetoMzmass && passLocalZmass && !passLocalBveto ) {
                    mon.fillHisto("redMet_shapes_NRBctrl"+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto("redMet_rebin_shapes_NRBctrl"+varNames[ivar],NRBtags,index,3,iweight);
                    //RJ
                    mon.fillHisto(TString("mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto(TString("new_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto(TString("balance_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto(TString("dphiLL_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto(TString("coslZ_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto(TString("balance1Dshape_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto(TString("dphiLL1Dshape_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto(TString("coslZ1Dshape_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);

                    mon.fillHisto(TString("dphi_shapes_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto(TString("met_shapes_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto(TString("met_rebin_shapes_NRBctrl")+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto("zpt_shapes_NRBctrl"+varNames[ivar],NRBtags,index,3,iweight);
                    mon.fillHisto("zpt_rebin_shapes_NRBctrl"+varNames[ivar],NRBtags,index,3,iweight);
                }
                if( passPreselectionMjvetoMbvetoMzmass && isZsideBand && !passLocalBveto ) {
                    mon.fillHisto("redMet_shapes_NRBctrl"+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto("redMet_rebin_shapes_NRBctrl"+varNames[ivar],NRBtags,index,4,iweight);
                    //RJ
                    mon.fillHisto(TString("mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto(TString("new_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto(TString("balance_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto(TString("dphiLL_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto(TString("coslZ_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto(TString("balance1Dshape_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto(TString("dphiLL1Dshape_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto(TString("coslZ1Dshape_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);

                    mon.fillHisto(TString("dphi_shapes_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto(TString("met_shapes_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto(TString("met_rebin_shapes_NRBctrl")+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto("zpt_shapes_NRBctrl"+varNames[ivar],NRBtags,index,4,iweight);
                    mon.fillHisto("zpt_rebin_shapes_NRBctrl"+varNames[ivar],NRBtags,index,4,iweight);
                }
                if( passPreselectionMjvetoMbvetoMzmass && isZsideBandPlus && !passLocalBveto ) {
                    mon.fillHisto("redMet_shapes_NRBctrl"+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto("redMet_rebin_shapes_NRBctrl"+varNames[ivar],NRBtags,index,5,iweight);
                    //RJ
                    mon.fillHisto(TString("mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto(TString("new_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto(TString("balance_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto(TString("dphiLL_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto(TString("coslZ_mt_shapes_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto(TString("balance1Dshape_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto(TString("dphiLL1Dshape_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto(TString("coslZ1Dshape_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);

                    mon.fillHisto(TString("dphi_shapes_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto(TString("met_shapes_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto(TString("met_rebin_shapes_NRBctrl")+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto("zpt_shapes_NRBctrl"+varNames[ivar],NRBtags,index,5,iweight);
                    mon.fillHisto("zpt_rebin_shapes_NRBctrl"+varNames[ivar],NRBtags,index,5,iweight);
                }

                //changing eventCategory for NRB control sample


            }
        }
    }

    printf("\n");
    file->Close();

    //##############################################
    //########     SAVING HISTO TO FILE     ########
    //##############################################
    //save control plots to file
    outUrl += "/";
    outUrl += outFileUrl + ".root";
    printf("Results save in %s\n", outUrl.Data());

    //save all to the file
    TFile *ofile=TFile::Open(outUrl, "recreate");
    mon.Write();
    stmkr->Write();
    if (runSystematics){
    stmkr_jerup->Write();
    stmkr_jerdown->Write();
    stmkr_jesup->Write();
    stmkr_jesdown->Write();
    stmkr_umetup->Write();
    stmkr_umetdown->Write();
    stmkr_lesup->Write();
    stmkr_lesdown->Write();
    stmkr_puup->Write();
    stmkr_pudown->Write();
    stmkr_btagup->Write();
    stmkr_btagdown->Write();
    }
    ofile->Close();

    if(outTxtFile)fclose(outTxtFile);
}





