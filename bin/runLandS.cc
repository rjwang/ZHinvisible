#include <iostream>
#include <boost/shared_ptr.hpp>
#include "Math/GenVector/Boost.h"

#include <sstream>
#include "CMGTools/HtoZZ2l2nu/src/tdrstyle.C"
#include "CMGTools/HtoZZ2l2nu/src/JSONWrapper.cc"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/MacroUtils.h"
#include "CMGTools/HtoZZ2l2nu/interface/plotter.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>

using namespace std;
double NonResonnantSyst = 0.25; //0.1;//0.25;
double GammaJetSyst = 1.0; //0.5;//0.5, 1.0;


// Map for non-shape systematics (lnN)
std::map<TString, double> normSysts;
void initNormalizationSysts();

//wrapper for a projected shape for a given set of cuts
class Shape_t {
public:
    TH1* data, *totalBckg;
    std::vector<TH1 *> bckg, signal, bckgInSignal;
    //the key corresponds to the proc name
    //the key is the name of the variation: e.g. jesup, jesdown, etc.
    std::map<TString,std::vector<std::pair<TString, TH1*> > > bckgVars, signalVars, bckgInSignalVars;

    std::map<TString, double> xsections;
    std::map<TString, double> BRs;

    Shape_t() {}
    ~Shape_t() {}


    void clear() {
        std::cout<<"shape is destructed...";
        if(data)delete data;
        if(totalBckg)totalBckg->Delete();
        for(unsigned int i=0; i<bckg.  size(); i++) {
            delete bckg  [i];
        }
        bckg  .clear();
        for(unsigned int i=0; i<signal.size(); i++) {
            delete signal[i];
        }
        signal.clear();
        for(std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator it=bckgVars  .begin(); it!=bckgVars  .end(); it++) {
            for(unsigned int i=0; i<(*it).second.size(); i++) {
                delete (*it).second[i].second;
            }
        }
        bckgVars  .clear();
        for(std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator it=signalVars.begin(); it!=signalVars.end(); it++) {
            for(unsigned int i=0; i<(*it).second.size(); i++) {
                delete (*it).second[i].second;
            }
        }
        signalVars.clear();
        std::cout<<"done\n";

    }


};

typedef std::pair<TString,TString> RateKey_t;
struct DataCardInputs {
    TString shapesFile;
    std::vector<TString> ch;
    std::vector<TString> procs;
    std::map<RateKey_t, Double_t> obs;
    std::map<RateKey_t, Double_t> rates;
    std::map<TString, std::map<RateKey_t,Double_t> > systs;
    int nsignalproc;
};


void printHelp();
Shape_t getShapeFromFile(TFile* inF, TString ch, TString shapeName, int cutBin,JSONWrapper::Object &Root,double minCut=0, double maxCut=9999, bool onlyData=false);
void showShape(std::vector<TString>& selCh ,map<TString, Shape_t>& allShapes, TString mainHisto, TString SaveName);
void Draw1DHistogram(TH1* mc, THStack *stack, TH1 *mcPlusRelUnc, TGraphErrors *errors, std::vector<TH1 *>& spimpose, TH1 *data, TLegend *legA, bool noLog, TString finstate, TString AnalysisBins);
void Draw2DHistogram(std::map<TString, TH1*>& mapbkgs, std::vector<TH1 *>& spimpose, TH1 *data, TString finstate, TString AnalysisBins);

void getYieldsFromShape(std::vector<TString> ch, const map<TString, Shape_t> &allShapes, TString shName);
void getEffFromShape(std::vector<TString> ch, const map<TString, Shape_t> &allShapes, TString shName);




void convertHistosForLimits_core(DataCardInputs& dci, TString& proc, TString& bin, TString& ch, std::vector<TString>& systs, std::vector<TH1*>& hshapes);
DataCardInputs convertHistosForLimits(Int_t mass,TString histo="finalmt",TString url="plotter.root",TString Json="");
DataCardInputs convertHistosForLimits(TString atgcpar,Int_t mass,TString histo="finalmt",TString url="plotter.root",TString Json="");
std::vector<TString> buildDataCard(Int_t mass, TString histo="finalmt", TString url="plotter.root",TString Json="");
std::vector<TString> buildDataCard(TString atgcpar,Int_t mass, TString histo="finalmt", TString url="plotter.root",TString Json="");
void doBackgroundSubtraction(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t> &allShapes, TString mainHisto, TString sideBandHisto, TString url, JSONWrapper::Object &Root);
void doDYReplacement(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString metHistoForRescale);
void doWZSubtraction(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString sideBandHisto);
void BlindData(std::vector<TString>& selCh, map<TString, Shape_t>& allShapes, TString mainHisto, bool addSignal);

void RescaleForInterference(std::vector<TString>& selCh,map<TString, Shape_t>& allShapes, TString mainHisto);
void SignalInterpolation(std::vector<TString>& selCh,map<TString, Shape_t>& allShapesL, map<TString, Shape_t>& allShapes, map<TString, Shape_t>& allShapesR, TString mainHisto);


void setTGraph(TString proc, TString suffix);
void initializeTGraph();

TGraph *ggH7TG_xsec=NULL, *ggH7TG_errp=NULL, *ggH7TG_errm=NULL, *ggH7TG_scap=NULL, *ggH7TG_scam=NULL, *ggH7TG_pdfp=NULL, *ggH7TG_pdfm=NULL;
TGraph *qqH7TG_xsec=NULL, *qqH7TG_errp=NULL, *qqH7TG_errm=NULL, *qqH7TG_scap=NULL, *qqH7TG_scam=NULL, *qqH7TG_pdfp=NULL, *qqH7TG_pdfm=NULL;
TGraph *ggH8TG_xsec=NULL, *ggH8TG_errp=NULL, *ggH8TG_errm=NULL, *ggH8TG_scap=NULL, *ggH8TG_scam=NULL, *ggH8TG_pdfp=NULL, *ggH8TG_pdfm=NULL;
TGraph *qqH8TG_xsec=NULL, *qqH8TG_errp=NULL, *qqH8TG_errm=NULL, *qqH8TG_scap=NULL, *qqH8TG_scam=NULL, *qqH8TG_pdfp=NULL, *qqH8TG_pdfm=NULL;
TGraph *    TG_xsec=NULL, *    TG_errp=NULL, *    TG_errm=NULL, *    TG_scap=NULL, *    TG_scam=NULL, *    TG_pdfp=NULL, *    TG_pdfm=NULL;
// Maps with fit values for ATGC
std::map<TString, std::vector<double> > extrvalpoints;
std::map<TString, std::vector<double> > extrerrpoints;

void initializeAtgcMaps();

TGraph* TG_QCDScaleK0ggH0=NULL, *TG_QCDScaleK0ggH1=NULL, *TG_QCDScaleK1ggH1=NULL, *TG_QCDScaleK1ggH2=NULL, *TG_QCDScaleK2ggH2=NULL;
TGraph* TG_UEPSf0=NULL, *TG_UEPSf1=NULL, *TG_UEPSf2=NULL;

bool subNRB2011 = false;
bool subNRB2012 = false;
bool MCclosureTest = false;

bool mergeWWandZZ = false;
bool skipWW = true;
std::vector<TString> Channels;
std::vector<TString> AnalysisBins;
bool fast = false;
bool skipGGH = false;
bool skipQQH = false;
bool subDY = false;
bool subWZ = false;
double DDRescale = 1.0;
double MCRescale = 1.0;
bool blindData = false;
bool blindWithSignal = false;
TString DYFile ="";
TString inFileUrl(""),jsonFile(""), histo("");
TString postfix="";
TString systpostfix="";
double shapeMin = 0;
double shapeMax = 9999;
double shapeMinVBF = 0;
double shapeMaxVBF = 9999;
bool doInterf = false;

int indexvbf = -1;
int indexcut   = -1, indexcutL=-1, indexcutR=-1;
int mass=-1, massL=-1, massR=-1;
TString atgcpar="", atgcpar2="";
bool runSystematics = false;
bool shape = false;
float sysSherpa=1.;

void initNormalizationSysts()
{
    normSysts["lumi_7TeV"] = 0.022;
    normSysts["lumi_8TeV"] = 0.026;
    normSysts["accept_7TeV"] = 0.;//0.02;//0.003; //RJ
    normSysts["accept_8TeV"] = 0.;//0.02;//0.018; //RJ
    normSysts["sherpa_kin_syst"] = sysSherpa-1.0;
    normSysts["CMS_eff_e"] = 0.03;
    normSysts["CMS_eff_m"] = 0.04;
    normSysts["CMS_scale_e"] = 0.01; // do we need it? There's shape uncertainty "les"...
    normSysts["CMS_scale_m"] = 0.01; // do we need it? There's shape uncertainty "les"...
    normSysts["QCDscale_VV_zz_7Tev"] = 0.07021;
    normSysts["QCDscale_VV_wz_7Tev"] = 0.059;
    normSysts["QCDscale_VV_zz_8Tev"] = 0.0944;
    normSysts["QCDscale_VV_wz_8Tev"] = 0.054;
    normSysts["QCDscale_VV1in_zz_7Tev"] = 0.91657-1.0; // s'ha da fa' cosi'...
    normSysts["QCDscale_VV1in_zz_8Tev"] = 0.92617-1.0; // s'ha da fa' cosi'...
    normSysts["pdf_VV_zz_7TeV"] = 0.0115;
    normSysts["pdf_VV_wz_7TeV"] = 0.0116;
    normSysts["pdf_VV_zz_8TeV"] = 0.0112;
    normSysts["pdf_VV_wz_8TeV"] = 0.0120;
    normSysts["sys_zlldata_7TeV"] = GammaJetSyst;
    normSysts["sys_topwwwjetsdata_7TeV"] = NonResonnantSyst;
    normSysts["sys_zlldata_8TeV"] = GammaJetSyst;
    normSysts["sys_topwwwjetsdata_8TeV"] = NonResonnantSyst;
    //
    normSysts["CMS_zh2l2v_mumueq0jets_leptonVeto"] = 0.01;
    normSysts["CMS_zh2l2v_eeeq0jets_leptonVeto"] = 0.01;
    normSysts["CMS_zh2l2v_mumueq1jets_leptonVeto"] = 0.013;
    normSysts["CMS_zh2l2v_eeeq1jets_leptonVeto"] = 0.013;
}

void printHelp()
{
    printf("Options\n");
    printf("--in        --> input file with from plotter\n");
    printf("--json      --> json file with the sample descriptor\n");
    printf("--histo     --> name of histogram to be used\n");
    printf("--shapeMin  --> left cut to apply on the shape histogram\n");
    printf("--shapeMax  --> right cut to apply on the shape histogram\n");
    printf("--shapeMinVBF  --> left cut to apply on the shape histogram for Vbf bin\n");
    printf("--shapeMaxVBF  --> right cut to apply on the shape histogram for Vbf bin\n");
    printf("--indexvbf  --> index of selection to be used for the vbf bin (if unspecified same as --index)\n");
    printf("--index     --> index of selection to be used (Xbin in histogram to be used)\n");
    printf("--indexL    --> index of selection to be used (Xbin in histogram to be used) used for interpolation\n");
    printf("--indexR    --> index of selection to be used (Xbin in histogram to be used) used for interpolation\n");
    printf("--m         --> higgs mass to be considered\n");
    printf("--mL        --> higgs mass on the left  of the mass to be considered (used for interpollation\n");
    printf("--mR        --> higgs mass on the right of the mass to be considered (used for interpollation\n");
    printf("--atgc      --> aTGC parameter (ex. string format: \"f4z=-0.01\")\n");
    printf("--syst      --> use this flag if you want to run systematics, default is no systematics\n");
    printf("--shape     --> use this flag if you want to run shapeBased analysis, default is cut&count\n");
    printf("--subNRB    --> use this flag if you want to subtract non-resonant-backgounds similarly to what was done in 2011 (will also remove H->WW)\n");
    printf("--subNRB12  --> use this flag if you want to subtract non-resonant-backgounds using a new technique that keep H->WW\n");
    printf("--subDY     --> histogram that contains the Z+Jets background estimated from Gamma+Jets)\n");
    printf("--subWZ     --> use this flag if you want to subtract WZ background by the 3rd lepton SB)\n");
    printf("--DDRescale --> factor to be used in order to multiply/rescale datadriven estimations\n");
    printf("--closure   --> use this flag if you want to perform a MC closure test (use only MC simulation)\n");
    printf("--bins      --> list of bins to be used (they must be comma separated without space)\n");
    printf("--HWW       --> use this flag to consider HWW signal)\n");
    printf("--skipGGH   --> use this flag to skip GGH signal)\n");
    printf("--skipQQH   --> use this flag to skip GGH signal)\n");
    printf("--blind     --> use this flag to replace observed data by total predicted background)\n");
    printf("--blindWithSignal --> use this flag to replace observed data by total predicted background+signal)\n");
    printf("--fast      --> use this flag to only do assymptotic prediction (very fast but inaccurate))\n");
    printf("--postfix    --> use this to specify a postfix that will be added to the process names)\n");
    printf("--systpostfix    --> use this to specify a syst postfix that will be added to the process names)\n");
    printf("--MCRescale    --> use this to specify a syst postfix that will be added to the process names)\n");
    printf("--interf     --> use this to rescale xsection according to WW interferences)\n");
    printf("--aTGC_Syst    --> use this to specify that you want to add and extra systematic in the shape)\n");
}

//
int main(int argc, char* argv[])
{
    setTDRStyle();
    gStyle->SetPadTopMargin   (0.06);
    gStyle->SetPadBottomMargin(0.12);
    //gStyle->SetPadRightMargin (0.16);
    gStyle->SetPadRightMargin (0.06);
    gStyle->SetPadLeftMargin  (0.14);
    gStyle->SetTitleSize(0.04, "XYZ");
    gStyle->SetTitleXOffset(1.1);
    gStyle->SetTitleYOffset(1.45);
    gStyle->SetPalette(1);
    gStyle->SetNdivisions(505);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    //init the TGraphs
    initializeTGraph();

    //init the ATGC maps
    initializeAtgcMaps();

    //get input arguments
    for(int i=1; i<argc; i++) {
        string arg(argv[i]);
        if(arg.find("--help")          !=string::npos) {
            printHelp();
            return -1;
        } else if(arg.find("--subNRB12") !=string::npos) {
            subNRB2012=true;
            skipWW=false;
            printf("subNRB2012 = True\n");
        } else if(arg.find("--subNRB")   !=string::npos) {
            subNRB2011=true;
            skipWW=true;
            printf("subNRB2011 = True\n");
        } else if(arg.find("--subDY")    !=string::npos) {
            subDY=true;
            DYFile=argv[i+1];
            i++;
            printf("Z+Jets will be replaced by %s\n",DYFile.Data());
        } else if(arg.find("--subWZ")    !=string::npos) {
            subWZ=true;
            printf("WZ will be estimated from 3rd lepton SB\n");
        } else if(arg.find("--DDRescale")!=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&DDRescale);
            i++;
        } else if(arg.find("--MCRescale")!=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&MCRescale);
            i++;
        } else if(arg.find("--HWW")      !=string::npos) {
            skipWW=false;
            printf("HWW = True\n");
        } else if(arg.find("--skipGGH")  !=string::npos) {
            skipGGH=true;
            printf("skipGGH = True\n");
        } else if(arg.find("--skipQQH")  !=string::npos) {
            skipQQH=true;
            printf("skipQQH = True\n");
        } else if(arg.find("--blindWithSignal")    !=string::npos) {
            blindData=true;
            blindWithSignal=true;
            printf("blindData = True; blindWithSignal = True\n");
        } else if(arg.find("--blind")    !=string::npos) {
            blindData=true;
            printf("blindData = True\n");
        } else if(arg.find("--closure")  !=string::npos) {
            MCclosureTest=true;
            printf("MCclosureTest = True\n");
        } else if(arg.find("--shapeMinVBF") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&shapeMinVBF);
            i++;
            printf("Min cut on shape for VBF = %f\n", shapeMinVBF);
        } else if(arg.find("--shapeMaxVBF") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&shapeMaxVBF);
            i++;
            printf("Max cut on shape for VBF = %f\n", shapeMaxVBF);
        } else if(arg.find("--shapeMin") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&shapeMin);
            i++;
            printf("Min cut on shape = %f\n", shapeMin);
        } else if(arg.find("--shapeMax") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&shapeMax);
            i++;
            printf("Max cut on shape = %f\n", shapeMax);
        } else if(arg.find("--interf")    !=string::npos) {
            doInterf=true;
            printf("doInterf = True\n");
        } else if(arg.find("--indexvbf") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexvbf);
            i++;
            printf("indexVBF = %i\n", indexvbf);
        } else if(arg.find("--indexL")    !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexcutL);
            i++;
            printf("indexL = %i\n", indexcutL);
        } else if(arg.find("--indexR")    !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexcutR);
            i++;
            printf("indexR = %i\n", indexcutR);
        } else if(arg.find("--index")    !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexcut);
            i++;
            printf("index = %i\n", indexcut);
        } else if(arg.find("--in")       !=string::npos && i+1<argc)  {
            inFileUrl = argv[i+1];
            i++;
            printf("in = %s\n", inFileUrl.Data());
        } else if(arg.find("--json")     !=string::npos && i+1<argc)  {
            jsonFile  = argv[i+1];
            i++;
            printf("json = %s\n", jsonFile.Data());
        } else if(arg.find("--histo")    !=string::npos && i+1<argc)  {
            histo     = argv[i+1];
            i++;
            printf("histo = %s\n", histo.Data());
        } else if(arg.find("--mL")       !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&massL );
            i++;
            printf("massL = %i\n", massL);
        } else if(arg.find("--mR")       !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&massR );
            i++;
            printf("massR = %i\n", massR);
        } else if(arg.find("--m")        !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&mass );
            i++;
            printf("mass = %i\n", mass);
        } else if(arg.find("--atgc")     !=string::npos && i+1<argc)  {
            atgcpar = argv[i+1];
            i++;
            printf("aTGC parameter: %s\n", atgcpar.Data());
        } else if(arg.find("--bins")     !=string::npos && i+1<argc)  {
            char* pch = strtok(argv[i+1],",");
            printf("bins are : ");
            while (pch!=NULL) {
                printf(" %s ",pch);
                AnalysisBins.push_back(pch);
                pch = strtok(NULL,",");
            }
            printf("\n");
            i++;
        } else if(arg.find("--channels") !=string::npos && i+1<argc)  {
            char* pch = strtok(argv[i+1],",");
            printf("channels are : ");
            while (pch!=NULL) {
                printf(" %s ",pch);
                Channels.push_back(pch);
                pch = strtok(NULL,",");
            }
            printf("\n");
            i++;
        } else if(arg.find("--fast")     !=string::npos) {
            fast=true;
            printf("fast = True\n");
        } else if(arg.find("--postfix")   !=string::npos && i+1<argc)  {
            postfix = argv[i+1];
            systpostfix = argv[i+1];
            i++;
            printf("postfix '%s' will be used\n", postfix.Data());
        } else if(arg.find("--systpostfix")   !=string::npos && i+1<argc)  {
            systpostfix = argv[i+1];
            i++;
            printf("systpostfix '%s' will be used\n", systpostfix.Data());
        } else if(arg.find("--syst")     !=string::npos) {
            runSystematics=true;
            printf("syst = True\n");
        } else if(arg.find("--shape")    !=string::npos) {
            shape=true;
            printf("shapeBased = True\n");
        } else if(arg.find("--aTGC_Syst")      !=string::npos) {
            sscanf(argv[i+1],"%f",&sysSherpa);
            i++;
            printf("Additional systematic on the shape setted: %f\n", sysSherpa);
        }
    }
    if(jsonFile.IsNull() || inFileUrl.IsNull() || histo.IsNull() || indexcut == -1 || mass==-1) {
        printHelp();
        return -1;
    }
    if(AnalysisBins.size()==0)AnalysisBins.push_back("");
    if(Channels.size()==0) {
        Channels.push_back("ee");
        Channels.push_back("mumu");
        Channels.push_back("ll"); //RENJIE, add all
    }

    atgcpar2 = atgcpar;
    atgcpar2.ReplaceAll("-", "M");
    atgcpar2.ReplaceAll("+", "P");
    atgcpar2.ReplaceAll("=0.", "=P0.");
    atgcpar2.ReplaceAll(".", "p");
    atgcpar2.ReplaceAll("=", "_");

    //if(systpostfix.Contains('8')) {NonResonnantSyst = 0.25; GammaJetSyst = 0.40;}
    if(atgcpar.Length()>0) GammaJetSyst = 0.0;

    initNormalizationSysts();

    //build the datacard for this mass point
    //std::vector<TString> dcUrls = buildDataCard(mass,histo,inFileUrl, jsonFile);
    std::vector<TString> dcUrls = buildDataCard(atgcpar,mass,histo,inFileUrl, jsonFile);
}




//
Shape_t getShapeFromFile(TFile* inF, TString ch, TString shapeName, int cutBin, JSONWrapper::Object &Root, double minCut, double maxCut, bool onlyData)
{
    gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

    Shape_t shape;
    shape.totalBckg=NULL;
    shape.data=NULL;

    std::vector<TString> BackgroundsInSignal;

    //iterate over the processes required
    std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
    for(unsigned int i=0; i<Process.size(); i++) {
        TString procCtr("");
        procCtr+=i;
        TString proc=(Process[i])["tag"].toString();
        TDirectory *pdir = (TDirectory *)inF->Get(proc);
        if(pdir==0) {
            /*printf("Skip Proc=%s because its directory is missing in root file\n", proc.Data());*/ continue;
        }

        bool isData(Process[i]["isdata"].toBool());
        if(onlyData && !isData)continue; //just here to speedup the NRB prediction

        bool isSignal(Process[i].isTag("issignal") && Process[i]["issignal"].toBool());
        if(Process[i]["spimpose"].toBool() && (proc.Contains("ggH") || proc.Contains("qqH")))isSignal=true;
        bool isInSignal(Process[i].isTag("isinsignal") && Process[i]["isinsignal"].toBool());
        int color(1);
        if(Process[i].isTag("color" ) ) color  = (int)Process[i]["color" ].toInt();
        int lcolor(color);
        if(Process[i].isTag("lcolor") ) lcolor = (int)Process[i]["lcolor"].toInt();
        int mcolor(color);
        if(Process[i].isTag("mcolor") ) mcolor = (int)Process[i]["mcolor"].toInt();
        int fcolor(color);
        if(Process[i].isTag("fcolor") ) fcolor = (int)Process[i]["fcolor"].toInt();
        int lwidth(1);
        if(Process[i].isTag("lwidth") ) lwidth = (int)Process[i]["lwidth"].toInt();
        int lstyle(1);
        if(Process[i].isTag("lstyle") ) lstyle = (int)Process[i]["lstyle"].toInt();
        int fill(1001);
        if(Process[i].isTag("fill"  ) ) fill   = (int)Process[i]["fill"  ].toInt();
        int marker(20);
        if(Process[i].isTag("marker") ) marker = (int)Process[i]["marker"].toInt();

        TH1* syst = (TH1*)pdir->Get("optim_systs");
        if(syst==NULL) {
            syst=new TH1F("optim_systs","optim_systs",1,0,1);
            syst->GetXaxis()->SetBinLabel(1,"");
        }
        for(int ivar = 1; ivar<=syst->GetNbinsX(); ivar++) {
            TH1D* hshape   = NULL;

            TString varName = syst->GetXaxis()->GetBinLabel(ivar);
            TString histoName = ch+"_"+shapeName+varName ;
            TH2* hshape2D = (TH2*)pdir->Get(histoName );
            if(!hshape2D) {
//            if(varName==""){
                //replace by empty histogram (take inclusive histo to make sure it has same binning)
                hshape2D = (TH2*)pdir->Get(shapeName+varName);
                if(hshape2D)hshape2D->Reset();
//            }else{
//               continue;
//            }
            }

            if(hshape2D) {
                histoName.ReplaceAll(ch,ch+"_proj"+procCtr);
                hshape   = hshape2D->ProjectionY(histoName,cutBin,cutBin);
                if(hshape->Integral()<=0 && varName=="" && !isData) {
                    hshape->Reset();
                    hshape->SetBinContent(1, 1E-10);
                }

                if(isnan((float)hshape->Integral())) {
                    hshape->Reset();
                }
                hshape->SetDirectory(0);
                hshape->SetTitle(proc);
                fixExtremities(hshape,true,true);
                hshape->SetFillColor(color);
                hshape->SetLineColor(lcolor);
                hshape->SetMarkerColor(mcolor);
                hshape->SetFillStyle(fill);
                hshape->SetLineWidth(lwidth);
                hshape->SetMarkerStyle(marker);
                hshape->SetLineStyle(lstyle);
            } else {
                printf("Histo %s does not exist for syst:%s\n", histoName.Data(), varName.Data());
                continue;
            }


            //if current shape is the one to cut on, then apply the cuts
            if(shapeName == histo) {
                for(int x=0; x<=hshape->GetXaxis()->GetNbins()+1; x++) {
                    if(hshape->GetXaxis()->GetBinCenter(x)<=minCut || hshape->GetXaxis()->GetBinCenter(x)>=maxCut) {
                        hshape->SetBinContent(x,0);
                        hshape->SetBinError(x,0);
                    }
                }
                //hshape->Rebin(2);
                hshape->GetYaxis()->SetTitle("Entries (/25GeV)");
            }

            hshape->Scale(MCRescale);



            //save in structure
            if(isData) {
                if(varName=="")  shape.data=hshape;
                else continue;
            } else if(isSignal) {
                if(skipGGH && proc.Contains("ggH"))continue;
                if(skipQQH && proc.Contains("qqH"))continue;

                if(skipWW && string(proc.Data()).find("WW")!=string::npos )continue;
                if(!skipWW && mergeWWandZZ) {
                    proc.ReplaceAll("WW","VV");
                    proc.ReplaceAll("ZZ","VV");
                }

                if(varName=="") {
                    if(Process[i]["data"].daughters()[0].isTag("xsec"))shape.xsections[proc] = Process[i]["data"].daughters()[0]["xsec"].toDouble();
                    if(Process[i]["data"].daughters()[0].isTag("br")) {
                        std::vector<JSONWrapper::Object> BRs = Process[i]["data"].daughters()[0]["br"].daughters();
                        double totalBR=1.0;
                        for(size_t ipbr=0; ipbr<BRs.size(); ipbr++) {
                            totalBR*=BRs[ipbr].toDouble();
                        }
                        shape.BRs[proc] = totalBR;
                    }

                    int procIndex = -1;
                    for(unsigned int i=0; i<shape.signal.size(); i++) {
                        if(string(proc.Data())==shape.signal[i]->GetTitle() ) {
                            procIndex=i;
                            break;
                        }
                    }
                    if(procIndex>=0) shape.signal[procIndex]->Add(hshape);
                    else             {
                        hshape->SetTitle(proc);
                        shape.signal.push_back(hshape);
                    }

                    //printf("Adding signal %s\n",proc.Data());
                } else {
                    std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator it = shape.signalVars.find(proc);

                    bool newVar = true;
                    if(it!=shape.signalVars.end()) {
                        for(unsigned int i=0; i<it->second.size(); i++) {
                            if( string(it->second[i].first.Data()) == varName ) {
                                it->second[i].second->Add(hshape);
                                newVar=false;
                                break;
                            }
                        }
                    }

                    if(newVar) {
                        shape.signalVars[proc].push_back( std::pair<TString,TH1*>(varName,hshape) );
                    }
                }
            } else {
                if(isInSignal) {
                    if(varName=="")  BackgroundsInSignal.push_back(proc);
                    if(varName=="")  shape.bckgInSignal.push_back(hshape);
                    else             shape.bckgInSignalVars[proc].push_back( std::pair<TString,TH1*>(varName,hshape) );
                } else {
                    if(varName=="")  shape.bckg.push_back(hshape);
                    else             shape.bckgVars[proc].push_back( std::pair<TString,TH1*>(varName,hshape) );
                }

                //printf("histoName = B %i -- %i  -- %s - %s --> %s\n", i, int(varName==""), proc.Data(), histoName.Data(), hshape->GetTitle());
            }
        }
    }

    //compute the total
    for(size_t i=0; i<shape.bckg.size(); i++) {
        if(i==0) {
            shape.totalBckg = (TH1 *)shape.bckg[i]->Clone(ch+"_"+shapeName+"_total");
            shape.totalBckg->SetDirectory(0);
        } else     {
            shape.totalBckg->Add(shape.bckg[i]);
        }
    }

    if(MCclosureTest) {
        if(shape.totalBckg) {
            if(!shape.data) {
                shape.data=(TH1F*)shape.totalBckg->Clone("data");
                shape.data->SetDirectory(0);
                shape.data->SetTitle("data");
            } else {
                shape.data->Reset();
                shape.data->Add(shape.totalBckg, 1);
            }
        }
    }


    //subtract background that are included to signal sample
    //if(shape.signal.size()>1 && BackgroundsInSignal.size()>0)printf("YOU ARE TRYING TO SUBTRACT BACKGROUNDS FROM SIGNAL SAMPLE WHILE YOU HAVE MORE THAN ONE SIGNAL GIVEN!  ASSUME THAT FIRST SIGNAL SAMPLE IS THE ONE CONTAINING THE BACKGROUND\n");
    if(shape.signal.size()>1 && BackgroundsInSignal.size()>0) printf("YOU ARE TRYING TO SUBTRACT BACKGROUNDS FROM SIGNAL SAMPLE WHILE YOU HAVE MORE THAN ONE SIGNAL GIVEN! ASSUME THAT ALL SIGNAL SAMPLES CONTAIN THE BACKGROUND, AND SUBTRACT THEM ALL!\n");
    //for(size_t i=0;i<BackgroundsInSignal.size()  && shape.signal.size()>0  ;i++){
    //deal with central value
    //for(size_t j=0;j<shape.bckg.size();j++){
    for(size_t j=0; j<shape.bckgInSignal.size(); j++) {
        //printf("Compare %s and %s\n",shape.bckg[j]->GetTitle(), BackgroundsInSignal[i].Data());
        //if(BackgroundsInSignal[i]!=shape.bckg[j]->GetTitle())continue;
        //printf("Subtract %s to %s\n",shape.signal[0]->GetTitle(), BackgroundsInSignal[i].Data());
        //shape.signal[0]->Add(shape.bckg[j],-1);
        std::map<TString, std::vector<double> > atgcShiftMap;
        for(size_t k=0; k<shape.signal.size(); ++k) {
            printf("Subtract %s to %s\n", shape.signal[k]->GetTitle(), shape.bckgInSignal[j]->GetTitle());
            shape.signal[k]->Add(shape.bckgInSignal[j],-1);
            // Get bin values for interpolated points
            TString atgcCorrIdx = ch+"_"+shape.signal[k]->GetTitle();
            if(systpostfix.Contains("7"))      atgcCorrIdx.Prepend("7TeV_");
            else if(systpostfix.Contains("8")) atgcCorrIdx.Prepend("8TeV_");
            std::vector<double> atgcCorrVect; //, atgcErrVect;
            if(extrvalpoints.count(atgcCorrIdx.Data())>0) atgcCorrVect = extrvalpoints[atgcCorrIdx.Data()];
            else std::cout << "No ATGC correction found with label " << atgcCorrIdx.Data() << "!" << std::endl;
            //if(extrerrpoints.count(atgcCorrIdx.Data())>0) atgcErrVect = extrerrpoints[atgcCorrIdx.Data()];
            if(atgcCorrVect.size()!=((unsigned int)shape.signal[k]->GetNbinsX())) {
                std::cout << "ATGC correction has different size than signal shape!" << std::endl;
                std::cout << "           atgcCorrVect.size() = " << atgcCorrVect.size() << "  (" << atgcCorrIdx.Data() << ")" << std::endl;
                std::cout << "  shape.signal[k]->GetNbinsX() = " << shape.signal[k]->GetNbinsX() << "  (" << shape.signal[k]->GetTitle() << ")" << std::endl;
                atgcCorrVect.clear(); //atgcErrVect.clear();
            }
            std::vector<double> atgcShiftVec(atgcCorrVect.size(), 0.);
            for(int x=0; x<=shape.signal[k]->GetNbinsX()+1; x++) {
                if(shape.signal[k]->GetBinContent(x)<0)shape.signal[k]->SetBinContent(x,0.0);
                if(atgcCorrVect.size()==0) continue;
                if(x==0 || x==shape.signal[k]->GetNbinsX()+1) continue;
                if(atgcCorrVect[x-1]<0.) continue;
                atgcShiftVec[x-1] = atgcCorrVect[x-1] - shape.signal[k]->GetBinContent(x); // before changing, save the difference (the syst plots will be shifted by the same amount)
                shape.signal[k]->SetBinContent(x, atgcCorrVect[x-1]); // replace the content of the bin with the expected value for the extrapolated point
            }
            if(atgcShiftVec.size()>0) atgcShiftMap[atgcCorrIdx] = atgcShiftVec;
        }
        //}

        //deal with systematics
        //const std::vector<std::pair<TString, TH1*> >& bckgSysts = shape.bckgVars  [BackgroundsInSignal[i] ];
        const std::vector<std::pair<TString, TH1*> >& bckgSysts = shape.bckgInSignalVars[shape.bckgInSignal[j]->GetTitle()];
        for(size_t k=0; k<shape.signal.size(); ++k) {
            const std::vector<std::pair<TString, TH1*> >& signSysts = shape.signalVars[shape.signal[k]->GetTitle() ];
            if(bckgSysts.size()!=signSysts.size())printf("Problem the two vectors have different size!\n");
            TString atgcCorrIdx = ch+"_"+shape.signal[k]->GetTitle();
            if(systpostfix.Contains("7"))      atgcCorrIdx.Prepend("7TeV_");
            else if(systpostfix.Contains("8")) atgcCorrIdx.Prepend("8TeV_");
            for(size_t jsys=0; jsys<bckgSysts.size(); jsys++) {
                signSysts[jsys].second->Add(bckgSysts[jsys].second, -1);
                for(int x=0; x<=signSysts[jsys].second->GetNbinsX()+1; x++) {
                    if(signSysts[jsys].second->GetBinContent(x)<0)signSysts[jsys].second->SetBinContent(x,0.0);
                    if(atgcShiftMap.count(atgcCorrIdx)==0) continue;
                    if(x==0 || x==shape.signal[k]->GetNbinsX()+1) continue;
                    signSysts[jsys].second->SetBinContent( x, signSysts[jsys].second->GetBinContent(x) + atgcShiftMap[atgcCorrIdx][x-1] ); // apply to the syst plots the same shift used for the main distribution
                }
            }
        }
    }



    //all done
    return shape;
}

//
void showShape(std::vector<TString>& selCh ,map<TString, Shape_t>& allShapes, TString mainHisto, TString SaveName)
{
    std::map<TString, TH1*> CCHistos;


    TCanvas* c1 = new TCanvas("c1","c1",300*AnalysisBins.size(),300*selCh.size());
    c1->SetTopMargin(0.00);
    c1->SetRightMargin(0.00);
    c1->SetBottomMargin(0.00);
    c1->SetLeftMargin(0.00);
    TPad* t1 = new TPad("t1","t1", 0.03, 0.03, 1.0, 1.00);
    t1->Draw();
    t1->cd();
    t1->Divide(AnalysisBins.size(), selCh.size(), 0, 0);

    t1->cd(selCh.size()*AnalysisBins.size());
    TLegend* legA  = new TLegend(0.6,0.5,0.99,0.85, "NDC");
    legA->SetNColumns(1); //RJ
    t1->cd(0);

    bool haserrorsleg(false);

    for(size_t s=0; s<selCh.size(); s++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            TVirtualPad* pad = t1->cd(1+s*AnalysisBins.size()+b);
            pad->SetTopMargin(0.06);
            pad->SetRightMargin(0.03);
            pad->SetBottomMargin(0.07);
            pad->SetLeftMargin(0.06);

            TH1* allbkg=NULL;
            std::map<TString, TH1*> mapbkg;
            std::map<TString, TH1*> mapsig;
            TH1* alldata=NULL;

            Shape_t& shape = allShapes.find(selCh[s]+AnalysisBins[b]+mainHisto)->second;
            cout << "selCh[s]+AnalysisBins[b]+mainHisto: " << selCh[s]+AnalysisBins[b]+mainHisto << endl;

            if(!allbkg) {
                allbkg=(TH1*)shape.totalBckg->Clone("mc");
            } else {
                allbkg->Add(shape.totalBckg);
            }
            if(!alldata) {
                alldata=(TH1*)shape.data->Clone("data");
            } else {
                alldata->Add(shape.data);
            }

            for(size_t i=0; i<shape.bckg.size(); i++) {
                //cout << "----------------------------------------------" << endl;
                //cout << "shape.bckg[i]->GetTitle(): " << shape.bckg[i]->GetTitle() << "  Integral: " << shape.bckg[i]->Integral() << endl;
                //cout << "----------------------------------------------" << endl;
                if(shape.bckg[i]->Integral()<=1E-6) continue;
                cout << "\n ishape.bckg[i]->GetTitle(): " << shape.bckg[i]->GetTitle() << "\n";
                if(mapbkg.find(shape.bckg[i]->GetTitle())!=mapbkg.end()) {
                    mapbkg[shape.bckg[i]->GetTitle()]->Add(shape.bckg[i]);
                } else {
                    mapbkg[shape.bckg[i]->GetTitle()]=(TH1*)shape.bckg[i]->Clone(shape.bckg[i]->GetTitle());
                }

                //cut & count plot
                double VAL, ERR;
                VAL = shape.bckg[i]->IntegralAndError(1,shape.bckg[i]->GetXaxis()->GetNbins(),ERR);
                if(CCHistos.find(shape.bckg[i]->GetTitle())==CCHistos.end()) {
                    CCHistos[shape.bckg[i]->GetTitle()] = new TH1D(TString(shape.bckg[i]->GetTitle())+"CC", "", 25,0,25);
                    CCHistos[shape.bckg[i]->GetTitle()]->SetLineColor(shape.bckg[i]->GetLineColor());
                    CCHistos[shape.bckg[i]->GetTitle()]->SetLineStyle(shape.bckg[i]->GetLineStyle());
                    CCHistos[shape.bckg[i]->GetTitle()]->SetLineWidth(shape.bckg[i]->GetLineWidth());
                    CCHistos[shape.bckg[i]->GetTitle()]->SetFillColor(shape.bckg[i]->GetFillColor());
                    CCHistos[shape.bckg[i]->GetTitle()]->SetFillStyle(shape.bckg[i]->GetFillStyle());
                }
                CCHistos[shape.bckg[i]->GetTitle()]->SetBinContent(1+s*AnalysisBins.size()+b,VAL);
                CCHistos[shape.bckg[i]->GetTitle()]->SetBinError(1+s*AnalysisBins.size()+b,VAL);
            }

            TString massStr("");
            if(mass>0) massStr += mass;
            if(atgcpar.Length()>0) massStr = atgcpar;
            TString massStr2 = atgcpar2;
            for(size_t ip=0; ip<shape.signal.size(); ip++) {
                TString proc(shape.signal[ip]->GetTitle());
                if(mass>0 && !proc.Contains(massStr))continue;
                if(mass>0 && proc.Contains("ggH") && proc.Contains("ZZ"))proc = "ggHZZ2l2v";
                else if(mass>0 && proc.Contains("qqH") && proc.Contains("ZZ"))proc = "qqHZZ2l2v";
                else if(mass>0 && proc.Contains("ggH") && proc.Contains("WW"))proc = "ggHWW2l2v";
                else if(mass>0 && proc.Contains("qqH") && proc.Contains("WW"))proc = "qqHWW2l2v";
                else if(mass>0 && proc.Contains("ZH")                        )proc = "ZH"+massStr+"2lMET";

                if(proc=="qqHZZ") {
                    shape.signal[ip]->SetLineStyle(2);
                }

                //if(atgcpar.Length()>0 && !proc.Contains(massStr)) continue;
                if(atgcpar.Length()>0 && !proc.EndsWith(massStr)) continue;
                if(atgcpar.Length()>0) proc = "ZZ2l2nu_"+massStr2;

                if(mapsig.find(shape.signal[ip]->GetTitle())!=mapsig.end()) {
                    mapsig[shape.signal[ip]->GetTitle()]->Add(shape.signal[ip]);
                } else {
                    mapsig[shape.signal[ip]->GetTitle()]=(TH1*)shape.signal[ip]->Clone(shape.signal[ip]->GetTitle());
                }

                //cut & count plot
                double VAL, ERR;
                VAL = shape.signal[ip]->IntegralAndError(1,shape.signal[ip]->GetXaxis()->GetNbins(),ERR);
                if(CCHistos.find(shape.signal[ip]->GetTitle())==CCHistos.end()) {
                    CCHistos[shape.signal[ip]->GetTitle()] = new TH1D(TString(shape.signal[ip]->GetTitle())+"CC", "", 25,0,25);
                    CCHistos[shape.signal[ip]->GetTitle()]->SetLineColor(shape.signal[ip]->GetLineColor());
                    CCHistos[shape.signal[ip]->GetTitle()]->SetLineStyle(shape.signal[ip]->GetLineStyle());
                    CCHistos[shape.signal[ip]->GetTitle()]->SetLineWidth(shape.signal[ip]->GetLineWidth());
                    CCHistos[shape.signal[ip]->GetTitle()]->SetFillColor(shape.signal[ip]->GetFillColor());
                    CCHistos[shape.signal[ip]->GetTitle()]->SetFillStyle(shape.signal[ip]->GetFillStyle());
                }
                CCHistos[shape.signal[ip]->GetTitle()]->SetBinContent(1+s*AnalysisBins.size()+b,VAL);
                CCHistos[shape.signal[ip]->GetTitle()]->SetBinError(1+s*AnalysisBins.size()+b,ERR);
            }

            THStack *stack=0;
            TH1* mc=allbkg;
            TH1* mcPlusRelUnc=0;
            TGraphErrors* errors=0;

            if(allbkg) {
                mcPlusRelUnc = (TH1 *) allbkg->Clone("totalmcwithunc");
                mcPlusRelUnc->SetDirectory(0);
                mcPlusRelUnc->Reset();
                stack = new THStack("stack","stack");
                for(std::map<TString, TH1*>::iterator it=mapbkg.begin(); it!=mapbkg.end(); ++it) {
                    it->second->SetLineColor( it->second->GetFillColor());
                    stack->Add(it->second,"HIST");
                    if(s==0 && b==0)legA->AddEntry(it->second,it->second->GetTitle(),"F");

                    double baseRelUnc = it->second->GetBinError(0)/it->second->Integral();
                    //if(TString(it->second->GetTitle()).Contains("rightarrow ll (data)")){printf("replace uncertainty %g",baseRelUnc); baseRelUnc=1.0;}
                    for(int ibin=1; ibin<=mcPlusRelUnc->GetXaxis()->GetNbins(); ibin++) {
                        double val = it->second->GetBinContent(ibin);
                        double err = it->second->GetBinError(ibin);
                        double value = mcPlusRelUnc->GetBinContent(ibin) + val;
                        double error = sqrt(pow(mcPlusRelUnc->GetBinError(ibin),2) + pow(err,2) + pow(val*baseRelUnc,2));

                        //Syst Err from Shape
                        std::vector<float> ErrShape_tot;
                        ErrShape_tot.clear();
                        // Loop on proc
                        std::cout << " ************** " << it->first.Data() << " ************** " << std::endl;
                        for(std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator jt=shape.bckgVars.begin(); jt!=shape.bckgVars.end(); ++jt) {
                            if(it->first.Contains("data")) continue;
                            if(it->first.CompareTo(jt->first)!=0) continue;
                            std::cout << " ************** " << jt->first.Data() << " ************** " << std::endl;
                            float Err(0);
                            //cout<<"proc: "<<jt->first<<endl; //Z#rightarrow ll
                            std::vector<float> ErrShape_up;
                            ErrShape_up.clear();
                            std::vector<float> ErrShape_down;
                            ErrShape_down.clear();
                            // Loop on syst
                            for(int isel=0; isel<14 ; isel++) {
                                if((jt->second)[isel].first.Contains("up")) {
                                    TH1F *h1 = (TH1F*) (jt->second)[isel].second;
                                    float err_tmp = h1->GetBinContent(ibin);
                                    err_tmp = (err_tmp-value)/value;
                                    //cout<<"Syst Name: "<<(jt->second)[isel].first<<" : "<<h1->GetBinContent(ibin)<<" "<<value<<" -> "<<err_tmp<<endl;
                                    if(err_tmp!=-1 && value!=0) ErrShape_up.push_back(err_tmp);
                                    else                        ErrShape_up.push_back(0.);
                                } else if((jt->second)[isel].first.Contains("down")) {
                                    TH1F *h1 = (TH1F*) (jt->second)[isel].second;
                                    float err_tmp = h1->GetBinContent(ibin);
                                    err_tmp = (err_tmp-value)/value;
                                    //cout<<"Syst Name: "<<(jt->second)[isel].first<<" : "<<h1->GetBinContent(ibin)<<" "<<value<<" -> "<<err_tmp<<endl;
                                    if(err_tmp!=-1 && value!=0) ErrShape_down.push_back(err_tmp);
                                    else                        ErrShape_down.push_back(0.);
                                } else {
                                    cout<<"WARNING: Syst not computed properly"<<endl;
                                }
                            }//syst
                            std::vector<float> ErrShape;
                            ErrShape.clear();
                            for(unsigned int isel=0; isel<ErrShape_down.size() ; isel++) {
                                if(ErrShape_down.size() != ErrShape_up.size()) cout<<"WARNING: Syst not computed properly"<<endl;
                                ErrShape.push_back( std::max(fabs(ErrShape_down[isel]), fabs(ErrShape_up[isel])) );
                            }
                            // Non-shape uncertainties ("lnN")
                            //for(unsigned int inormsys=0; inormsys<normSysts.size(); ++inormsys) {
                            for(std::map<TString, double>::const_iterator kt=normSysts.begin(); kt!=normSysts.end(); ++kt) {
                                if(kt->first.Contains("7TeV") && !systpostfix.Contains("7TeV")) continue;
                                if(kt->first.Contains("8TeV") && !systpostfix.Contains("8TeV")) continue;
                                if(kt->first.EndsWith("_e") && !selCh[s].Contains("ee")) continue;
                                if(kt->first.EndsWith("_m") && !selCh[s].Contains("mumu")) continue;
                                if(kt->first.Contains("_zz_") && !it->first.Contains("ZZ")) continue;
                                if(kt->first.Contains("_wz_") && !it->first.Contains("WZ")) continue;
                                if(kt->first.Contains("_zlldata_") && !it->first.Contains("Z#rightarrow ll (data)")) continue;
                                if(kt->first.Contains("_topwwwjetsdata_") && !it->first.Contains("Top/WW/W+Jets  (data)")) continue;

                                ErrShape.push_back(kt->second);
                            }
                            float TotErrShape(0);
                            for(unsigned int iEr=0; iEr<ErrShape.size(); iEr++) TotErrShape+=pow(ErrShape[iEr],2);
                            TotErrShape=sqrt(TotErrShape);
                            //cout<<"Sing process Unc:"<<TotErrShape<<endl;
                            ErrShape_tot.push_back(TotErrShape);
                        }//process
                        float TotErrShape_tot(0.);
                        for(unsigned int iEr=0; iEr<ErrShape_tot.size(); iEr++) TotErrShape_tot+=pow(ErrShape_tot[iEr],2);
                        TotErrShape_tot=sqrt(TotErrShape_tot);
                        cout<<"BIN unc:"<<TotErrShape_tot<<" , the other is: "<<error<<" summed are: "<<sqrt(TotErrShape_tot*TotErrShape_tot+error*error)<<endl;

                        //Syst Err not form datacard (harcoded)

                        mcPlusRelUnc->SetBinContent(ibin,value);
                        mcPlusRelUnc->SetBinError(ibin,error);
                    }
                }


                double mcErorMax=0;
                errors = new TGraphErrors(mcPlusRelUnc->GetXaxis()->GetNbins());
                int icutg=0;
                for(int ibin=1; ibin<=mcPlusRelUnc->GetXaxis()->GetNbins(); ibin++) {
                    if(mcPlusRelUnc->GetBinContent(ibin)>0)
                        errors->SetPoint(icutg,mcPlusRelUnc->GetXaxis()->GetBinCenter(ibin), mcPlusRelUnc->GetBinContent(ibin));
                    errors->SetPointError(icutg,mcPlusRelUnc->GetXaxis()->GetBinWidth(ibin)/2.0, mcPlusRelUnc->GetBinError(ibin));
                    icutg++;
                    mcErorMax = std::max(mcErorMax, mcPlusRelUnc->GetBinContent(ibin) + mcPlusRelUnc->GetBinError(ibin));
                }

                TH1* axis = (TH1*)allbkg->Clone("axis");
                axis->Reset();
                axis->GetXaxis()->SetTitle(mc->GetXaxis()->GetTitle());
                axis->GetYaxis()->SetTitle(b==0?mc->GetYaxis()->GetTitle():"");
                //axis->SetMinimum(mc->GetMinimum());
                //axis->SetMaximum(1.1*std::max(mcErorMax, alldata->GetMaximum()));
                axis->SetMaximum(1.9*std::max(mcErorMax, alldata->GetMaximum()));
                //axis->GetXaxis()->SetRangeUser(150,700);//
                axis->GetXaxis()->SetRangeUser(0,400);//
                axis->Draw();
                stack->Draw("same");


                errors->Set(icutg);
                //errors->SetFillStyle(3427);
                //errors->SetFillStyle(3002);
                errors->SetFillStyle(3254); //RJ
                errors->SetFillColor(kGray+3);
                errors->SetLineStyle(1);
                errors->SetLineColor(0);
                errors->Draw("2 same");
                if(!haserrorsleg) {
                    legA->AddEntry(errors,"Syst.+Stat. Unc.","F"); //RJ
                    haserrorsleg=true;
                }
            }


            std::vector<TH1 *> sigvec;
            for(std::map<TString, TH1*>::iterator it=mapsig.begin(); it!=mapsig.end(); it++) {
                it->second->SetFillStyle(0);
                it->second->Draw("hist same");
                if(s==0 && b==0)legA->AddEntry(it->second,it->second->GetTitle(),"L");
                sigvec.push_back(it->second);
            }



            if(alldata) {
                alldata->Draw("E1 same");
                if(s==0 && b==0)legA->AddEntry(alldata,alldata->GetTitle(),"PL");
            }


            TPaveText* Label = new TPaveText(0.1,0.81,0.94,0.89, "NDC");
            Label->SetFillColor(0);
            Label->SetFillStyle(0);
            Label->SetLineColor(0);
            Label->SetBorderSize(0);
            Label->SetTextAlign(31);
            TString LabelText = selCh[s]+"  -  "+AnalysisBins[b];
            LabelText.ReplaceAll("mumu","#mu#mu");
            LabelText.ReplaceAll("geq2jets","#geq2jets");
            LabelText.ReplaceAll("eq0jets","0jet");
            LabelText.ReplaceAll("eq1jets","1jet");
            Label->AddText(LabelText);
            Label->Draw();



            if(s==selCh.size()-1 && b==AnalysisBins.size()-1) {
                legA->SetFillColor(0);
                legA->SetFillStyle(0);
                legA->SetLineColor(0);
                legA->SetBorderSize();
                legA->SetHeader("");
                legA->Draw("same");
                legA->SetTextFont(42);
            }

            Draw1DHistogram(mc, stack, mcPlusRelUnc, errors, sigvec, alldata, legA, true, selCh[s], AnalysisBins[b]);
            if(histo.Contains("dphiLL_mt_shapes"))
                Draw2DHistogram(mapbkg, sigvec, alldata, selCh[s], AnalysisBins[b]);


            /*
              TH1 *ratio=0;
              if(allbkg && alldata){
                  c1->cd();
                  TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.2);     t2->Draw();
                  t2->cd();
                  t2->SetGridy(true);
                  t2->SetTopMargin(0);   t2->SetBottomMargin(0.5);
                  float yscale = (1.0-0.2)/(0.18-0);
                  TH1 *ratio = (TH1*)alldata->Clone("RatioHistogram");
                  ratio->SetDirectory(0);
                  ratio->Divide(allbkg);
                  ratio->GetYaxis()->SetTitle("Obs/Ref");
                  ratio->GetXaxis()->SetTitle("");
                  ratio->SetMinimum(0);
                  ratio->SetMaximum(2.2);
                  ratio->GetXaxis()->SetTitleOffset(1.3);
                  ratio->GetXaxis()->SetLabelSize(0.033*yscale);
                  ratio->GetXaxis()->SetTitleSize(0.036*yscale);
                  ratio->GetXaxis()->SetTickLength(0.03*yscale);
                  ratio->GetYaxis()->SetTitleOffset(0.3);
                  ratio->GetYaxis()->SetNdivisions(5);
                  ratio->GetYaxis()->SetLabelSize(0.033*yscale);
                  ratio->GetYaxis()->SetTitleSize(0.036*yscale);
                  ratio->GetXaxis()->SetRangeUser(175,450);
                  ratio->Draw("E1");
              }
            */
            pad->SetLogy();
            pad->Update();
        }
    }

    t1->cd();
    t1->Update();
    c1->cd();
    TPaveText* T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetBorderSize(0);
    T->SetTextAlign(22);
    if(systpostfix.Contains('8')) {
        T->AddText("CMS preliminary, ZH #rightarrow ll+MET, #sqrt{s}=8.0 TeV, #scale[0.5]{#int} L=19.6  fb^{-1}");
    } else {
        T->AddText("CMS preliminary, ZH #rightarrow ll+MET, #sqrt{s}=7.0 TeV, #scale[0.5]{#int} L=5.1  fb^{-1}");
    }
    T->Draw();
    c1->Update();
    c1->SaveAs(SaveName+"_Shape.eps");
    c1->SaveAs(SaveName+"_Shape.png");
    c1->SaveAs(SaveName+"_Shape.pdf");
    c1->SaveAs(SaveName+"_Shape.C");
    delete c1;
//  if(alldata)delete alldata;
//  if(allbkg)delete allbkg;
//  if(stack) delete stack;
//  if(ratio) delete ratio;
}


void Draw1DHistogram(TH1* mc, THStack *stack, TH1 *mcPlusRelUnc, TGraphErrors *errors, std::vector<TH1 *>& spimpose, TH1 *data, TLegend *legA, bool noLog, TString finstate, TString AnalysisBins)
{

    bool showsigvsBkg(true);
    bool isDataBlind = true;

    TCanvas* c1d = new TCanvas("c1d","c1d",750,800);//,800,800);

    TPad* t1d = new TPad();
    if(!isDataBlind) t1d = new TPad("t1d","t1d", 0.0, 0.3, 1.0, 1.0);
    else t1d = new TPad("t1d","t1d", 0.0, 0.0, 1.0, 1.0);

    t1d->Draw();
    t1d->cd();
    if(!isDataBlind) t1d->SetBottomMargin(0); //RJ
    TString name_denRelUncH;
    if(!noLog) t1d->SetLogy(true);
    if(histo.Contains("redMet_rebin_shapes") ||
            histo.Contains("zpt_rebin_shapes") ) t1d->SetLogx(true);
    float maximumFound(noLog);

    if(mc &&  maximumFound<mc->GetMaximum()) maximumFound=mc->GetMaximum() * (noLog ? (systpostfix.Contains('8') ? 1.2 : 1.4) : 1.1);
    if(data && maximumFound<data->GetMaximum()) maximumFound=data->GetMaximum() * (noLog ? (systpostfix.Contains('8') ? 1.2 : 1.4) : 1.1);

    bool canvasIsFilled(false);
    if(stack && stack->GetStack() && stack->GetStack()->GetEntriesFast()>0) {
        stack->Draw("");
        TH1 *hist=(TH1*)stack->GetStack()->At(0);
        //Set YTitle, RJ
        float binsize = hist->GetBinWidth(1);
        std::ostringstream strs;
        strs << binsize;
        TString binSize = strs.str();
        if(stack->GetXaxis()) {
            stack->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
            if(histo.Contains("new_mt_shapes")) stack->GetXaxis()->SetTitle("M_{T}(ll,#slash{E}_{T}) [GeV]");
            if(isDataBlind) stack->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
            //stack->GetYaxis()->SetTitle("Entries");
            TString xtitle = hist->GetXaxis()->GetTitle();
            if(xtitle.Contains("GeV") )
                stack->GetYaxis()->SetTitle("Events / "+binSize+" GeV");
            else	stack->GetYaxis()->SetTitle("Events / "+binSize);

            stack->GetXaxis()->SetTitleSize(0.05);
            stack->GetYaxis()->SetTitleSize(0.05);
            stack->GetXaxis()->SetMoreLogLabels(true);
            name_denRelUncH = hist->GetXaxis()->GetTitle(); //RJ
            if(histo.Contains("redMet_rebin_shapes")) stack->GetXaxis()->SetLimits(40., stack->GetXaxis()->GetXmax());
            if(histo.Contains("zpt_rebin_shapes")) stack->GetXaxis()->SetLimits(40., stack->GetXaxis()->GetXmax());
            if(histo.Contains("redMet_shapes")) {
                stack->GetXaxis()->SetLimits(40., 400.);
            }
            if(histo.Contains("zpt_shapes")) {
                stack->GetXaxis()->SetLimits(40., 400.);
            }

            stack->SetMinimum(hist->GetMinimum());
            stack->SetMaximum(maximumFound);
            if(noLog) {
                stack->SetMaximum(1.2*maximumFound);
                //stack->GetXaxis()->SetRangeUser(hist->GetMinimum(),maximumFound);
            }
        }
        canvasIsFilled = true;
    }

    stack->GetYaxis()->SetTitleOffset(1.3);
    stack->GetYaxis()->SetLabelSize(0.04);
    stack->GetYaxis()->SetTitleSize(0.04);

    stack->GetXaxis()->SetTitleOffset(1.3);
    stack->GetXaxis()->SetLabelSize(0.04);
    stack->GetXaxis()->SetTitleSize(0.04);



    if(errors) {
        errors->Draw("2 same");
    }

    if(data) {
        //if(!isDataBlind)
        //data->Draw(canvasIsFilled ? "E1 same" : "E1"); //RJ, blind data point for Higgs-Exotic Meeting
        canvasIsFilled = true;
    }

    for(size_t ip=0; ip<spimpose.size(); ip++) {
        //double mrksz = spimpose[ip]->GetMarkerSize(); spimpose[ip]->SetMarkerSize(mrksz*0.1); spimpose[ip]->Draw((canvasIsFilled ? "e1same": "e1") ); canvasIsFilled = true; //spimpose[ip]->SetMarkerSize(mrksz);
        TString opt = "hist";
        spimpose[ip]->Draw( (opt + (canvasIsFilled ? "same": "")).Data() );
        canvasIsFilled = true;
    }

    TPaveText* T = new TPaveText(0.1,0.995,0.9,0.95, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);
    char Buffer[1024];
    bool isSim = data ? false : true;
    double iEcm  = systpostfix.Contains('8') ? 8.0 : 7.0;
    double iLumi = systpostfix.Contains('8') ? 19577 : 5051;


    T = new TPaveText(0.18,0.9,0.33,0.8, "NDC");
    sprintf(Buffer, "#splitline{#bf{CMS}}{Preliminary}");
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);

    T = new TPaveText(0.35,0.9,0.5,0.8, "NDC");
    if(finstate.Contains("mumu"))   sprintf(Buffer, "#splitline{#it{ZH #rightarrow l^{+}l^{-}+#slash{E}_{T}}}{#it{#mu#mu channel}}");
    if(finstate.Contains("ee"))     sprintf(Buffer, "#splitline{#it{ZH #rightarrow l^{+}l^{-}+#slash{E}_{T}}}{#it{ee channel}}");
    if(finstate.Contains("ll"))     sprintf(Buffer, "#splitline{#it{ZH #rightarrow l^{+}l^{-}+#slash{E}_{T}}}{#it{ll channel}}");
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);

    T = new TPaveText(0.35,0.8,0.5,0.75, "NDC");
    sprintf(Buffer, "#sqrt{s} = %.1f TeV", iEcm);
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);

    T = new TPaveText(0.35,0.75,0.52,0.7, "NDC");
    sprintf(Buffer, "#scale[1.1]{#font[12]{#int}}Ldt = %.1f fb^{-1}",iLumi/1000);
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);


    legA->SetX1(0.55);
    legA->SetY1(0.65);
    legA->SetX2(0.95);
    legA->SetY2(0.92);
    legA->SetFillColor(0);
    legA->SetFillStyle(0);
    legA->SetLineColor(0);
    legA->SetHeader("");
    legA->Draw("same");
    legA->SetTextFont(42);

    std::vector<TH1 *> compDists;
    if(data)                   compDists.push_back(data);
    else if(spimpose.size()>0) compDists = spimpose;

    if(mc && compDists.size() && !isDataBlind) {
        c1d->cd();
        TPad* t2d = new TPad("t2d","t2d", 0.0, 0.0, 1.0, 0.3);
        t2d->Draw();
        t2d->cd();
        t2d->SetGridy(true);
        t2d->SetTopMargin(0);
        t2d->SetBottomMargin(0.5);
        if(histo.Contains("redMet_rebin_shapes") ||
                histo.Contains("zpt_rebin_shapes") ) t2d->SetLogx(true);

        // MC stats + syst
        TH1D *denRelUncH=0;
        if(mcPlusRelUnc) denRelUncH=(TH1D *) mcPlusRelUnc->Clone("mcrelunc");
        else             denRelUncH=(TH1D *) mc->Clone("mcrelunc");
        for(int xbin=1; xbin<=denRelUncH->GetXaxis()->GetNbins(); xbin++) {
            if(denRelUncH->GetBinContent(xbin)==0) continue;
            Double_t err=denRelUncH->GetBinError(xbin)/denRelUncH->GetBinContent(xbin);
            denRelUncH->SetBinContent(xbin,1);
            denRelUncH->SetBinError(xbin,err);
        }
        TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH);
        denRelUnc->SetLineColor(1);
        denRelUnc->SetFillStyle(3001);
        denRelUnc->SetFillColor(kGray);
        denRelUnc->SetMarkerColor(1);
        denRelUnc->SetMarkerStyle(1);
        denRelUncH->Reset("ICE");
        denRelUncH->Draw();
        denRelUnc->Draw("3");
        float yscale = (1.0-0.2)/(0.18-0);
        denRelUncH->GetYaxis()->SetTitle("Data/#Sigma MC");
        denRelUncH->SetMinimum(0.1);//0.4);
        denRelUncH->SetMaximum(1.9);//1.6);
        denRelUncH->GetXaxis()->SetTitle("");
        //denRelUncH->SetMinimum(0);
        //denRelUncH->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.10);
        denRelUncH->GetXaxis()->SetTitleOffset(1.3);
        denRelUncH->GetXaxis()->SetLabelSize(0.033*yscale);
        denRelUncH->GetXaxis()->SetTitleSize(0.036*yscale);
        denRelUncH->GetXaxis()->SetTickLength(0.03*yscale);
        denRelUncH->GetXaxis()->SetTitle(name_denRelUncH); //RJ
        denRelUncH->GetYaxis()->SetTitleOffset(0.3);
        denRelUncH->GetYaxis()->SetNdivisions(5);
        denRelUncH->GetYaxis()->SetLabelSize(0.033*yscale);
        denRelUncH->GetYaxis()->SetTitleSize(0.036*yscale);
        denRelUncH->GetXaxis()->SetMoreLogLabels(true);
        //denRelUncH->GetXaxis()->SetNdivisions(9);
        if(histo.Contains("redMet_rebin_shapes")) denRelUncH->GetXaxis()->SetRangeUser(40., denRelUncH->GetXaxis()->GetXmax());
        if(histo.Contains("zpt_rebin_shapes")) denRelUncH->GetXaxis()->SetRangeUser(40, denRelUncH->GetXaxis()->GetXmax());
        if(histo.Contains("redMet_shapes")) denRelUncH->GetXaxis()->SetRangeUser(40., 400.);
        if(histo.Contains("zpt_shapes")) denRelUncH->GetXaxis()->SetRangeUser(40, 400.);

        // Add comparisons
        for(size_t icd=0; icd<compDists.size(); icd++) {
            TString name("CompHistogram");
            name+=icd;
            TH1D *dataToObsH = (TH1D*)compDists[icd]->Clone(name);
            dataToObsH->Divide(mc);
            if(!isDataBlind) dataToObsH->Draw("same"); //RJ, Blind data point
        }
    } else {
        //if not comparison resize the canvas
        //c1d->SetWindowSize(600,400);
        //c1d->SetCanvasSize(600,400);
        t1d->SetPad(0,0,1,1);
    }

    c1d->Modified();
    c1d->Update();
    c1d->cd();
    //cout << "AnalysisBins: " << AnalysisBins << endl;

    string SavePath = (finstate + AnalysisBins + "_" + histo + systpostfix + ".eps").Data();
    while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
    while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
    while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
    while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
    while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
    while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
    while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
    while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
    system(string(("rm -f ") + SavePath).c_str());
    c1d->SaveAs(SavePath.c_str());
    while(SavePath.find(".eps")!=std::string::npos)SavePath.replace(SavePath.find(".eps"),4,".png");
    c1d->SaveAs(SavePath.c_str());
    while(SavePath.find(".png")!=std::string::npos)SavePath.replace(SavePath.find(".png"),4,".pdf");
    c1d->SaveAs(SavePath.c_str());
    delete c1d;
    delete T;


    /******************/
    if(showsigvsBkg) {
        TCanvas* cc = new TCanvas("cc","cc",800,600);//,800,800);
        //cc->cd();
        //TPad* t2 = new TPad("t2d","t2d", 0.0, 0.0, 1.0, 0.3);
        //t2d->Draw();
        //t2d->cd();
        //t2d->SetGridy(true);
        //t2d->SetTopMargin(0);
        //t2d->SetBottomMargin(0.5);

        // Add comparisons
        for(size_t ip=0; ip<spimpose.size(); ip++) {
            TH1D *sigTobkg = (TH1D*)spimpose[ip]->Clone();
            sigTobkg->Divide(mc);
            sigTobkg->Draw();

            float binsize = sigTobkg->GetBinWidth(1);
            std::ostringstream strs;
            strs << binsize;
            TString binSize = strs.str();

            sigTobkg->GetYaxis()->SetTitle("Sig/Bkg / "+binSize);
        }


        //TPaveText* T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
        TPaveText* T = new TPaveText(0.1,0.995,0.9,0.95, "NDC");
        T->SetFillColor(0);
        T->SetFillStyle(0);
        T->SetLineColor(0);
        T->SetTextAlign(22);
        char Buffer[1024];
        bool isSim = data ? false : true;
        double iEcm  = systpostfix.Contains('8') ? 8.0 : 7.0;
        double iLumi = systpostfix.Contains('8') ? 19577 : 5051;
        //isSim = true; //RJ
        if(isSim) sprintf(Buffer, "CMS simulation, ZH #rightarrow ll+MET, #sqrt{s}=%.1f TeV, #scale[0.7]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
        else      sprintf(Buffer, "CMS preliminary, ZH #rightarrow ll+MET, #sqrt{s}=%.1f TeV, #scale[0.7]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
        T->AddText(Buffer);
        T->Draw("same");
        T->SetBorderSize(0);

        //RJ adding ee or mumu channel number!
        TPaveText* T2 = new TPaveText(0.2,0.91,0.3,0.81, "NDC");
        T2->SetFillColor(0);
        T2->SetFillStyle(0);
        T2->SetLineColor(0);
        T2->SetTextAlign(22);
        char Buffer2[1024];
        if(finstate.Contains("ee"))       sprintf(Buffer2,"(ee)");
        if(finstate.Contains("mumu"))     sprintf(Buffer2,"(#mu#mu)");
        if(finstate.Contains("lleq"))	  sprintf(Buffer2,"(ee+#mu#mu)");
        T2->AddText(Buffer2);
        T2->Draw("same");

        string mysavePath_eps = (finstate + AnalysisBins + "_" + histo + systpostfix + "_sigToBkg.eps").Data();
        string mysavePath_png = (finstate + AnalysisBins + "_" + histo + systpostfix + "_sigToBkg.png").Data();
        string mysavePath_pdf = (finstate + AnalysisBins + "_" + histo + systpostfix + "_sigToBkg.pdf").Data();
        cc->SaveAs(mysavePath_eps.c_str());
        cc->SaveAs(mysavePath_png.c_str());
        cc->SaveAs(mysavePath_pdf.c_str());
        delete cc;
    }
    /******************/
}


void Draw2DHistogram(std::map<TString, TH1*>& mapbkgs, std::vector<TH1 *>& spimpose, TH1 *data, TString finstate, TString AnalysisBins)
{

    double DphiLLXaxis[4];//BalanceXaxis[4], DphiLLXaxis[4], CosThetaXaxis[4];
    //BalanceXaxis[0] = 0.8;
    //BalanceXaxis[1] = 0.9;
    //BalanceXaxis[2] = 1.1;
    //BalanceXaxis[3] = 1.2;
    DphiLLXaxis[0] = 0.;
    DphiLLXaxis[1] = 0.6;
    DphiLLXaxis[2] = 1.4;
    DphiLLXaxis[3] = TMath::Pi();
    //CosThetaXaxis[0] = -1.;
    //CosThetaXaxis[1] = -0.5;
    //CosThetaXaxis[2] = 0.;
    //CosThetaXaxis[3] = 1.;

    std::map<TString, TH1*> toDraw;
    TH1* totBkg=NULL;
    for(std::map<TString, TH1*>::iterator it=mapbkgs.begin(); it!=mapbkgs.end(); ++it) {
        toDraw[it->first]=it->second ;
        if(!totBkg) totBkg=(TH1*)it->second->Clone();
        else	totBkg->Add(it->second);
    }
    for(size_t ip=0; ip<spimpose.size(); ip++) {
        toDraw[spimpose[ip]->GetTitle()]=spimpose[ip];
    }
    toDraw["Data"]=data;
    toDraw["Background"]=totBkg;

    for(std::map<TString, TH1*>::iterator it=toDraw.begin(); it!=toDraw.end(); ++it) {
        TCanvas* c_1 = new TCanvas("c_1","c_1",800,800);
        c_1->SetTopMargin(0.08);
        c_1->SetRightMargin(0.15);
        //c_1->SetBottomMargin(0.00);
        //c_1->SetLeftMargin(0.00);
        //TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.0);
        //t1->Draw();
        //t1->cd();

        //t1->SetPadTopMargin   (0.06);
        //t1->SetPadBottomMargin(0.12);
        //t1->SetPadRightMargin (0.1);
        //t1->SetPadLeftMargin  (0.14);

        TH2F *hist = new TH2F("",";M_{T} [GeV];#Delta#phi_{ll} [rad];Events",12,0,1200,3,DphiLLXaxis);
        hist->SetTitleOffset(1.2,"yz");
        TH1* hist_id = it->second;
        TString iTitle = it->first;
        for(int ybin=1; ybin<=3; ybin++) {
            for(int xbin=1; xbin<=12; xbin++) {
                double value = hist_id->GetBinContent(xbin+12*(ybin-1));
                //printf("val: \t%f\n",value);
                hist->SetBinContent(xbin,ybin,value);
            }
        }
        gStyle->SetPalette(1);
        hist->SetMinimum(-1.0e-10);
        hist->Draw("COLZ");


        TPaveText* T = new TPaveText(0.4,0.98,0.86,0.90, "NDC");
        T->SetFillColor(0);
        T->SetFillStyle(0);
        T->SetLineColor(0);
        T->SetTextAlign(22);
        T->SetTextFont(42);
        char Buffer[1024];
        bool isSim = data ? false : true;
        double iEcm  = systpostfix.Contains('8') ? 8.0 : 7.0;
        double iLumi = systpostfix.Contains('8') ? 19577 : 5051;
        if(finstate.Contains("mumu")) sprintf(Buffer, "#it{ZH #rightarrow ll+#slash{E}_{T}}, #it{#mu#mu channel}, #sqrt{s} = %.1f TeV, L = %.1f fb^{-1}", iEcm, iLumi/1000);
        if(finstate.Contains("ee")) sprintf(Buffer, "#it{ZH #rightarrow ll+#slash{E}_{T}}, #it{ee channel}, #sqrt{s} = %.1f TeV, L = %.1f fb^{-1}", iEcm, iLumi/1000);
        T->AddText(Buffer);
        T->Draw("same");
        T->SetBorderSize(0);

        T = new TPaveText(0.58,0.995,0.86,0.95, "NDC");
        isSim = true;
        if(isSim) sprintf(Buffer, "#bf{CMS Preliminary}");
        T->AddText(Buffer);
        T->Draw("same");
        T->SetBorderSize(0);

        T = new TPaveText(0.12,0.98,0.3,0.90, "NDC");
        if(iTitle.Contains("ZH(105)"))		sprintf(Buffer, "SM M_{H}=105GeV");
        else if(iTitle.Contains("ZH(115)"))	sprintf(Buffer, "SM M_{H}=115GeV");
        else if(iTitle.Contains("ZH(125)"))	sprintf(Buffer, "SM M_{H}=125GeV");
        else if(iTitle.Contains("ZH(135)"))     sprintf(Buffer, "SM M_{H}=135GeV");
        else if(iTitle.Contains("ZH(145)"))     sprintf(Buffer, "SM M_{H}=145GeV");
        else if(iTitle.Contains("Data")) 	sprintf(Buffer, "  Data         ");
        else if(iTitle.Contains("WZ"))  	sprintf(Buffer, " WZ#rightarrow 3l#nu     ");
        else if(iTitle.Contains("ZZ")) 		sprintf(Buffer, " ZZ#rightarrow 2l2#nu    ");
        else if(iTitle.Contains("WW")) {
            T = new TPaveText(0.12,0.98,0.37,0.90, "NDC");
            sprintf(Buffer, iTitle);
        } else sprintf(Buffer, iTitle);
        T->AddText(Buffer);
        T->Draw("same");
        T->SetBorderSize(0);

        string SavePath = (finstate + AnalysisBins + "_" + histo + systpostfix + "_" + iTitle + ".eps").Data();
        while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
        while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
        while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
        while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
        while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
        while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
        while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
        while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"_");
        while(SavePath.find(" ")!=std::string::npos)SavePath.replace(SavePath.find(" "),1,"_");
        c_1->SaveAs(SavePath.c_str());
        while(SavePath.find(".eps")!=std::string::npos)SavePath.replace(SavePath.find(".eps"),4,".png");
        c_1->SaveAs(SavePath.c_str());
        while(SavePath.find(".png")!=std::string::npos)SavePath.replace(SavePath.find(".png"),4,".pdf");
        c_1->SaveAs(SavePath.c_str());
        delete c_1;
    }

}





//
void getYieldsFromShape(std::vector<TString> ch, const map<TString, Shape_t> &allShapes, TString shName)
{
    FILE* pFile = fopen("Yields.tex","w");
    fprintf(pFile,"\\begin{table}[htbp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n\\footnotesize\n");

    string Ccol   = "\\begin{tabular}{|c|";
    string Cname  = "channel";
    string Cval   = "";

    TString massStr("");
    if(mass>0)massStr += mass;
    if(atgcpar.Length()>0) massStr = atgcpar;
    TString massStr2 = atgcpar2;


    TH1* h;
    Double_t valerr, val, syst;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t ich=0; ich<ch.size(); ich++) {
            TString icol(ch[ich]+" "+AnalysisBins[b]);
            icol.ReplaceAll("mu","\\mu");
            icol.ReplaceAll("_"," ");
            Cval = "$ "+icol+" $";

            //bckg
            size_t nbckg=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.bckg.size();
            for(size_t ibckg=0; ibckg<nbckg; ibckg++) {
                TH1* h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.bckg[ibckg];
                TString procTitle(h->GetTitle());
                if(procTitle.Contains("QCD"))continue;
                if(procTitle.Contains("W#rightarrow l#nu"))continue;
                if(procTitle.Contains("Z#rightarrow #tau#tau"))continue;
                procTitle.ReplaceAll("#","\\");
                if(b==0&&ich==0)Ccol  += "c|";
                if(b==0&&ich==0)Cname += "&$" + procTitle + "$";

                val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
                syst = h->GetBinError(0)<=0 ? -1 : h->GetBinError(0);
                if(val<1E-6) {
                    val=0.0;
                    valerr=0.0;
                    syst=-1;
                }
                Cval += "&" + toLatexRounded(val,valerr, syst);
            }

            //total bckg
            if(b==0&&ich==0)Ccol  += "c|";
            if(b==0&&ich==0)Cname += "&$Total$";
            h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.totalBckg;
            val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
            syst = h->GetBinError(0)<=0 ? -1 : h->GetBinError(0);
            if(val<1E-6) {
                val=0.0;
                valerr=0.0;
                syst=-1;
            }
            Cval += "&\\boldmath " + toLatexRounded(val,valerr,syst);

            //signal
            size_t nsig=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.signal.size();
            for(size_t isig=0; isig<nsig; isig++) {
                h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.signal[isig];
                TString procTitle(h->GetTitle());
                procTitle.ReplaceAll("#","\\");

                if(mass>0 && !procTitle.Contains(massStr))continue;
                if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("ZZ"))procTitle = "ggH("+massStr+")";
                else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("ZZ"))procTitle = "qqH("+massStr+")";
                else if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("WW"))procTitle = "ggH("+massStr+")WW";
                else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("WW"))procTitle = "qqH("+massStr+")WW";
                else if(mass>0 && procTitle.Contains("ZH")                             )procTitle = "ZH("+massStr+")2lMET";

                //if(atgcpar.Length()>0 && !procTitle.Contains(massStr)) continue;
                if(atgcpar.Length()>0 && !procTitle.EndsWith(massStr)) continue;
                if(atgcpar.Length()>0) procTitle = "ZZ("+massStr+")";

                if(b==0&&ich==0)Ccol  += "c|";
                if(b==0&&ich==0)Cname += "&$" + procTitle+"$";

                val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
                if(val<1E-6) {
                    val=0.0;
                    valerr=0.0;
                }
                Cval += "&" + toLatexRounded(val,valerr);
            }

            //data
            //if(b==0&&ich==0)Ccol  += "c|";
            //if(b==0&&ich==0)Cname += "&$Data$";
            h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.data;
            val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
            if(val<1E-6) {
                val=0.0;
                valerr=0.0;
            }
            char tmpchar[255];
            ////sprintf(tmpchar,"%.0f",val);
            //Cval += "&\\boldmath $" + string(tmpchar)+"$";
//    Cval += "&\\boldmath " + toLatexRounded(val,0.0);

            //endline
            if(b==0&&ich==0)fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
            if(b==0&&ich==0)fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
            fprintf(pFile,"%s\\\\\n", Cval.c_str());
        }
    }
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
    fprintf(pFile,"\n\n\n\n");


//SAME THING BUT SUMMING UP ALL CHANNELS



//    fprintf(pFile,"\\begin{table}[H]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");

    Ccol   = "\\begin{tabular}{|c|";
    Cname  = "channel";
    Cval   = "";

    double VAL = 0;
    double VALERR = 0;
    double SYST = 0;

    for(size_t b=0; b<AnalysisBins.size(); b++) {
        //TString icol(AnalysisBins[b]);
        TString icol("ll");
        icol.ReplaceAll("mu","\\mu");
        icol.ReplaceAll("_"," ");
        Cval = "$ "+icol+" $";

        //bckg
        size_t nbckg=allShapes.find(ch[0]+AnalysisBins[b]+shName)->second.bckg.size();
        for(size_t ibckg=0; ibckg<nbckg; ibckg++) {
            VAL = 0;
            VALERR = 0;
            SYST = 0;
            for(size_t ich=0; ich<ch.size(); ich++) {
                TH1* h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.bckg[ibckg];
                TString procTitle(h->GetTitle());
                if(procTitle.Contains("QCD"))continue;
                if(procTitle.Contains("W#rightarrow l#nu"))continue;
                if(procTitle.Contains("Z#rightarrow #tau#tau"))continue;
                procTitle.ReplaceAll("#","\\");
                if(b==0&&ich==0)Ccol  += "c|";
                if(b==0&&ich==0)Cname += "&$" + procTitle + "$";

                val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
                syst = h->GetBinError(0)<=0 ? -1 : h->GetBinError(0);
                if(val<1E-6) {
                    val=0.0;
                    valerr=0.0;
                    syst=-1;
                }
                VAL += val;


                cout << "procTitle: " << procTitle << endl;
                if(procTitle.Contains("rightarrow ll (data)") ||
                        procTitle.Contains("Top/WW/W+Jets (data)")) { //RENJIE
                    VALERR=VALERR+valerr;
                    SYST+=syst;
                } //RENJIE, stat. & syst. 100% correlated
                else 	VALERR=sqrt(VALERR*VALERR+valerr*valerr);
                SYST=syst;
            }
            if(toLatexRounded(VAL,VALERR, SYST) == "") continue;
            Cval += "&" + toLatexRounded(VAL,VALERR, SYST);
        }

        //total bckg
        if(b==0)Ccol  += "c|";
        if(b==0)Cname += "&$Total$";
        VAL = 0;
        VALERR = 0;
        SYST = 0;
        for(size_t ich=0; ich<ch.size(); ich++) {
            h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.totalBckg;
            val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
            syst = h->GetBinError(0)<=0 ? -1 : h->GetBinError(0);
            if(val<1E-6) {
                val=0.0;
                valerr=0.0;
                syst=-1;
            }
            VAL += val;
            VALERR=sqrt(VALERR*VALERR+valerr*valerr);
            SYST=syst;
        }
        Cval += "&\\boldmath " + toLatexRounded(VAL,VALERR,SYST);

        //signal
        size_t nsig=allShapes.find(ch[0]+AnalysisBins[b]+shName)->second.signal.size();
        for(size_t isig=0; isig<nsig; isig++) {
            h=allShapes.find(ch[0]+AnalysisBins[b]+shName)->second.signal[isig];
            TString procTitle(h->GetTitle());
            procTitle.ReplaceAll("#","\\");

            if(mass>0 && !procTitle.Contains(massStr))continue;
            if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("ZZ"))procTitle = "ggH("+massStr+")";
            else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("ZZ"))procTitle = "qqH("+massStr+")";
            else if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("WW"))procTitle = "ggH("+massStr+")WW";
            else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("WW"))procTitle = "qqH("+massStr+")WW";
            else if(mass>0 && procTitle.Contains("ZH")                             )procTitle = "ZH("+massStr+")2lMET";

            //if(atgcpar.Length()>0 && !procTitle.Contains(massStr)) continue;
            if(atgcpar.Length()>0 && !procTitle.EndsWith(massStr)) continue;
            if(atgcpar.Length()>0) procTitle = "ZZ("+massStr+")";

            if(b==0)Ccol  += "c|";
            if(b==0)Cname += "&$" + procTitle+"$";

            VAL = 0;
            VALERR = 0;
            SYST = 0;
            for(size_t ich=0; ich<ch.size(); ich++) {
                h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.signal[isig];

                val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
                if(val<1E-6) {
                    val=0.0;
                    valerr=0.0;
                }
                VAL += val;
                VALERR=sqrt(VALERR*VALERR+valerr*valerr);
            }
            Cval += "&" + toLatexRounded(VAL, VALERR);
        }

        //data
        if(b==0)Ccol  += "c|";
        if(b==0)Cname += "&$Data$";
        VAL = 0;
        VALERR = 0;
        SYST = 0;
        for(size_t ich=0; ich<ch.size(); ich++) {
            h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.data;
            val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
            if(val<1E-6) {
                val=0.0;
                valerr=0.0;
            }
            VAL += val;
        }
        char tmpchar[255];
        sprintf(tmpchar,"%.0f",VAL);
        Cval += "&\\boldmath $" + string(tmpchar)+"$";
//    Cval += "&\\boldmath " + toLatexRounded(val,0.0);

        //endline
        //if(b==0)fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
        //if(b==0)fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
        //fprintf(pFile,"%s\\\\\n", Cval.c_str());
    }
    //fprintf(pFile,"\\hline\n");
    //fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
    //fprintf(pFile,"\n\n\n\n");






    fclose(pFile);
}


void getEffFromShape(std::vector<TString> ch, map<TString, Shape_t> &allShapes, TString shName)
{
    FILE* pFile = fopen("Efficiency.tex","w");
//  fprintf(pFile,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");

    string Ccol   = "\\begin{tabular}{|c|";
    string Cname  = "channel";
    string Cval   = "";

    TString massStr("");
    if(mass>0)massStr += mass;
    if(atgcpar.Length()>0) massStr = atgcpar;
    TString massStr2 = atgcpar2;


    TH1* h;
    Double_t valerr, val, syst;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t ich=0; ich<ch.size(); ich++) {
            TString icol(ch[ich]+"-"+AnalysisBins[b]);
            icol.ReplaceAll("mu","\\mu");
            icol.ReplaceAll("_"," ");
            Cval = "$ "+icol+" $";

            //signal
            size_t nsig=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.signal.size();
            for(size_t isig=0; isig<nsig; isig++) {
                Shape_t& shape  = allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second;
                h=shape.signal[isig];
                TString procTitle(h->GetTitle());
                procTitle.ReplaceAll("#","\\");

                if(mass>0 && !procTitle.Contains(massStr))continue;
                if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("ZZ"))procTitle = "ggH("+massStr+")";
                else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("ZZ"))procTitle = "qqH("+massStr+")";
                else if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("WW"))procTitle = "ggH("+massStr+")WW";
                else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("WW"))procTitle = "qqH("+massStr+")WW";
                else if(mass>0 && procTitle.Contains("ZH")                             )procTitle = "ZH("+massStr+")2lMET";

                //if(atgcpar.Length()>0 && !procTitle.Contains(massStr)) continue;
                if(atgcpar.Length()>0 && !procTitle.EndsWith(massStr)) continue;
                if(atgcpar.Length()>0) procTitle = "ZZ("+massStr+")";

                if(b==0&&ich==0)Ccol  += "c|";
                if(b==0&&ich==0)Cname += "&$" + procTitle+"$";

//       printf("%s --> Xsec=%E x %E\n",h->GetTitle(), shape.xsections[h->GetTitle()], shape.BRs[h->GetTitle()]);
                double xsecXbr = shape.xsections[h->GetTitle()] * shape.BRs[h->GetTitle()];

                val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
                if(val<1E-6) {
                    val=0.0;
                    valerr=0.0;
                }
                Cval += "&" + toLatexRounded(val/xsecXbr,valerr/xsecXbr);

                fprintf(pFile,"%30s %30s %4.0f %6.2E %6.2E %6.2E %6.2E\n",icol.Data(), procTitle.Data(), (double)mass, shape.xsections[h->GetTitle()], shape.BRs[h->GetTitle()], val/xsecXbr, valerr/xsecXbr);
            }

            //endline
//    if(b==0&&ich==0)fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
//    if(b==0&&ich==0)fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
//    fprintf(pFile,"%s\\\\\n", Cval.c_str());
        }
    }
//  fprintf(pFile,"\\hline\n");
//  fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");
//  fprintf(pFile,"\n\n\n\n");
    fclose(pFile);
}



std::vector<TString>  buildDataCard(Int_t mass, TString histo, TString url, TString Json)
{
    return buildDataCard("", mass, histo, url, Json);
}

std::vector<TString>  buildDataCard(TString atgcpar, Int_t mass, TString histo, TString url, TString Json)
{
    std::vector<TString> dcUrls;

    //get the datacard inputs
    DataCardInputs dci = convertHistosForLimits(atgcpar,mass,histo,url,Json);

    TString eecard = "";
    TString mumucard = "";
    TString combinedcard = "";
    TString combinedcardLL = "";

    //build the datacard separately for each channel
    for(size_t i=1; i<=dci.ch.size(); i++) {
        TString dcName=dci.shapesFile;
        dcName.ReplaceAll(".root","_"+dci.ch[i-1]+".dat");
        FILE* pFile = fopen(dcName.Data(),"w");

        if(!dci.ch[i-1].Contains("lleq")) combinedcard += dci.ch[i-1]+"="+dcName+" ";
        if(dci.ch[i-1].Contains("lleq")) combinedcardLL += dci.ch[i-1]+"="+dcName+" ";
        if(dci.ch[i-1].Contains("ee"))eecard += dci.ch[i-1]+"="+dcName+" ";
        if(dci.ch[i-1].Contains("mumu"))mumucard += dci.ch[i-1]+"="+dcName+" ";

        //header
        fprintf(pFile, "imax 1\n");
        fprintf(pFile, "jmax *\n");
        fprintf(pFile, "kmax *\n");
        fprintf(pFile, "-------------------------------\n");
        if(shape) {
            fprintf(pFile, "shapes * * %s %s/$PROCESS %s/$PROCESS_$SYSTEMATIC\n",dci.shapesFile.Data(), dci.ch[i-1].Data(), dci.ch[i-1].Data());
            fprintf(pFile, "-------------------------------\n");
        }
        //observations
        fprintf(pFile, "bin 1\n");
        fprintf(pFile, "Observation %f\n",dci.obs[RateKey_t("obs",dci.ch[i-1])]);
        fprintf(pFile, "-------------------------------\n");

        fprintf(pFile,"%45s ", "bin");
        for(size_t j=1; j<=dci.procs.size(); j++) {
            if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            fprintf(pFile,"%6i ", 1);
        }
        fprintf(pFile,"\n");

        fprintf(pFile,"%45s ", "process");
        for(size_t j=1; j<=dci.procs.size(); j++) {
            if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            fprintf(pFile,"%6s ", (postfix+dci.procs[j-1]).Data());
        }
        fprintf(pFile,"\n");

        fprintf(pFile,"%45s ", "process");
        int procCtr(1-dci.nsignalproc);
        for(size_t j=1; j<=dci.procs.size(); j++) {
            if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            fprintf(pFile,"%6i ", procCtr );
            procCtr++;
        }
        fprintf(pFile,"\n");

        fprintf(pFile,"%45s ", "rate");
        for(size_t j=1; j<=dci.procs.size(); j++) {
            if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            fprintf(pFile,"%6f ", dci.rates[RateKey_t(dci.procs[j-1],dci.ch[i-1])] );
        }
        fprintf(pFile,"\n");
        fprintf(pFile, "-------------------------------\n");


        //systematics
        char sFile[2048];
        bool isSyst;
        if(runSystematics) {
            if(systpostfix.Contains('8')) {
                fprintf(pFile,"%35s %10s ", "lumi_8TeV", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("data")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["lumi_8TeV"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            } else {
                fprintf(pFile,"%35s %10s ", "lumi_7TeV", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("data")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["lumi_7TeV"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }
            if(systpostfix.Contains('8')) {
                fprintf(pFile,"%35s %10s ", "accept_8TeV", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("data")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["accept_8TeV"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            } else {
                fprintf(pFile,"%35s %10s ", "accept_7TeV", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("data")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["accept_7TeV"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }
            if(sysSherpa != 1.) {
                fprintf(pFile,"%35s %10s ", "sherpa_kin_syst", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("zz2l2nu") && dci.procs[j-1].Contains("_") ) {
                        fprintf(pFile,"%6f ",sysSherpa);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }

//         if(mass>0){
//         fprintf(pFile,"%35s %10s ", "theoryUncXS_HighMH", "lnN");
//         for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
//            if((int)j<=dci.nsignalproc){fprintf(pFile,"%6f ",std::min(1.0+1.5*pow((mass/1000.0),3),2.0));}else{fprintf(pFile,"%6s ","-");}
//         }fprintf(pFile,"\n");
//         }

            //Id+Trigger efficiencies combined
            if(dci.ch[i-1].Contains("ee")) {
                fprintf(pFile,"%35s %10s ", "CMS_eff_e", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("data")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_eff_e"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            } else {
                fprintf(pFile,"%35s %10s ", "CMS_eff_mu", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("data")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_eff_m"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }



            //Oct 30, 2013

            if(dci.ch[i-1].Contains("mumueq0jets")) {
                fprintf(pFile,"%35s %10s ", "CMS_zh2l2v_mumueq0jets_leptonVeto", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("wz3lnu")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_zh2l2v_mumueq0jets_leptonVeto"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }

            if(dci.ch[i-1].Contains("mumueq1jets")) {
                fprintf(pFile,"%35s %10s ", "CMS_zh2l2v_mumueq1jets_leptonVeto", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("wz3lnu")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_zh2l2v_mumueq1jets_leptonVeto"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }



            if(dci.ch[i-1].Contains("eeeq0jets")) {
                fprintf(pFile,"%35s %10s ", "CMS_zh2l2v_eeeq0jets_leptonVeto", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("wz3lnu")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_zh2l2v_eeeq0jets_leptonVeto"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }

            if(dci.ch[i-1].Contains("eeeq1jets")) {
                fprintf(pFile,"%35s %10s ", "CMS_zh2l2v_eeeq1jets_leptonVeto", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("wz3lnu")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_zh2l2v_eeeq1jets_leptonVeto"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }



            // Scale-e and scale-mu replaced by shape-unncertainty "les"
//          if(dci.ch[i-1].Contains("ee")){
//             fprintf(pFile,"%35s %10s ", "CMS_scale_e", "lnN");
//             for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
//                if(!dci.procs[j-1].Contains("data")){fprintf(pFile,"%6f ",1.0+normSysts["CMS_scale_e"]);}else{fprintf(pFile,"%6s ","-");}
//             }fprintf(pFile,"\n");
//          }else{
//             fprintf(pFile,"%35s %10s ", "CMS_scale_m", "lnN");
//             for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
//                if(!dci.procs[j-1].Contains("data")){fprintf(pFile,"%6f ",1.0+normSysts["CMS_scale_m"]);}else{fprintf(pFile,"%6s ","-");}
//             }fprintf(pFile,"\n");
//          }

            if(mass>0) {
//         fprintf(pFile,"%35s %10s ", "Signal_rescaling_8TeV", "lnN");
//         for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
////            if(systpostfix.Contains('8') && (dci.procs[j-1].Contains("ggh") || dci.procs[j-1].Contains("qqh"))){fprintf(pFile,"%6f ",1.25);
//            if(systpostfix.Contains('8') && (dci.procs[j-1].Contains("qqh"))){fprintf(pFile,"%6f ",1.25);
//            }else{fprintf(pFile,"%6s ","-");}
//         }fprintf(pFile,"\n");


//verify cross section!
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("ggh") || dci.procs[j-1].Contains("qqh")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        printf("850 XSECTION=%6f ",TG_xsec->Eval(850,NULL,"S"));
                    }
                    if(dci.procs[j-1].Contains("ggh") || dci.procs[j-1].Contains("qqh")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        printf("900 XSECTION=%6f ",TG_xsec->Eval(900,NULL,"S"));
                    }
                    if(dci.procs[j-1].Contains("ggh") || dci.procs[j-1].Contains("qqh")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        printf("950 XSECTION=%6f ",TG_xsec->Eval(950,NULL,"S"));
                    }

                }


                fprintf(pFile,"%35s %10s ", "pdf_gg", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("ggh")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",1+0.01*sqrt(pow(TG_pdfp->Eval(mass,NULL,"S"),2) + pow(TG_pdfm->Eval(mass,NULL,"S"),2)));
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");

                fprintf(pFile,"%35s %10s ", "pdf_qqbar", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("zh")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",1+0.01*sqrt(pow(TG_pdfp->Eval(mass,NULL,"S"),2) + pow(TG_pdfm->Eval(mass,NULL,"S"),2)));
                    } else if(dci.procs[j-1].BeginsWith("zz")) {
                        if(systpostfix.Contains('8')) {
                            fprintf(pFile,"%6f ",1.0312);
                        } else {
                            fprintf(pFile,"%6f ",1.0360);
                        }
                    } else if(dci.procs[j-1].BeginsWith("wz")) {
                        if(systpostfix.Contains('8')) {
                            fprintf(pFile,"%6f ",1.0455);
                        } else {
                            fprintf(pFile,"%6f ",1.0502);
                        }
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");

                fprintf(pFile,"%35s %10s ", "UEPS", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("ggh") && dci.ch[i-1].Contains("eq0jet")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",TG_UEPSf0->Eval(mass,NULL,"S"));
                    } else if(dci.procs[j-1].Contains("ggh") && dci.ch[i-1].Contains("eq1jet")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",TG_UEPSf1->Eval(mass,NULL,"S"));
                    } else if(dci.procs[j-1].Contains("ggh") && (dci.ch[i-1].Contains("eq2jet") || dci.ch[i-1].Contains("vbf"))) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",TG_UEPSf2->Eval(mass,NULL,"S"));
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");

                fprintf(pFile,"%35s %10s ", "QCDscale_ggH", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("ggh") && dci.ch[i-1].Contains("eq0jet")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",TG_QCDScaleK0ggH0->Eval(mass,NULL,"S"));
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");

                fprintf(pFile,"%35s %10s ", "QCDscale_ggH1in", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("ggh") && dci.ch[i-1].Contains("eq0jet")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",TG_QCDScaleK0ggH1->Eval(mass,NULL,"S"));
                    } else if(dci.procs[j-1].Contains("ggh") && dci.ch[i-1].Contains("eq1jet")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",TG_QCDScaleK1ggH1->Eval(mass,NULL,"S"));
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");

                fprintf(pFile,"%35s %10s ", "QCDscale_ggH2in", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("ggh") && dci.ch[i-1].Contains("eq1jet")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",TG_QCDScaleK1ggH2->Eval(mass,NULL,"S"));
                    } else if(dci.procs[j-1].Contains("ggh") && (dci.ch[i-1].Contains("eq2jet") || dci.ch[i-1].Contains("vbf"))) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",TG_QCDScaleK2ggH2->Eval(mass,NULL,"S"));
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");

                fprintf(pFile,"%35s %10s ", "QCDscale_ZH", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("ggh")) {
                        setTGraph(dci.procs[j-1], systpostfix );
                        fprintf(pFile,"%6f ",1+0.01*sqrt(pow(TG_scap->Eval(mass,NULL,"S"),2) + pow(TG_scam->Eval(mass,NULL,"S"),2)));
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");

            }

            fprintf(pFile,"%35s %10s ", "QCDscale_ggVV", "lnN");
            for(size_t j=1; j<=dci.procs.size(); j++) {
                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                fprintf(pFile,"%6s ","-");
            }
            fprintf(pFile,"\n");
            /*
                     fprintf(pFile,"%35s %10s ", "QCDscale_VV", "lnN");
                     for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            	    if(dci.procs[j-1].BeginsWith("zz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.0+normSysts["QCDscale_VV_zz_8Tev"]);}else{fprintf(pFile,"%6f ",1.0+normSysts["QCDscale_VV_zz_7Tev"]);}  // from k(>=0)
                        //temporary removed to avoid double counts with uncertainty applied on ZZ xsection itself --> should be reintegrated for Higgs computation
                        }else if(dci.procs[j-1].BeginsWith("wz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.0+normSysts["QCDscale_VV_wz_8Tev"]);}else{fprintf(pFile,"%6f ",1.0+normSysts["QCDscale_VV_wz_7Tev"]);}
                        }else{fprintf(pFile,"%6s ","-");}
                     }fprintf(pFile,"\n");
            */
            /*   //RJ
                     fprintf(pFile,"%35s %10s ", "QCDscale_VV1in", "lnN");
                     for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                              if(dci.procs[j-1].BeginsWith("zz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.0+normSysts["QCDscale_VV1in_zz_8Tev"]);}else{fprintf(pFile,"%6f ",1.0+normSysts["QCDscale_VV1in_zz_7Tev"]);}  // from k(>=1)
                        //temporary removed to avoid double counts with uncertainty applied on ZZ xsection itself --> should be reintegrated for Higgs computation
                        }else{fprintf(pFile,"%6s ","-");}
                     }fprintf(pFile,"\n");

                     fprintf(pFile,"%35s %10s ", "pdf_VV", "lnN");
                     for(size_t j=1; j<=dci.procs.size(); j++) { if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            	          if(dci.procs[j-1].BeginsWith("zz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.0+normSysts["pdf_VV_zz_8TeV"]);}else{fprintf(pFile,"%6f ",1.0+normSysts["pdf_VV_zz_7TeV"]);}
            	    }else if(dci.procs[j-1].BeginsWith("wz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.0+normSysts["pdf_VV_wz_8TeV"]);}else{fprintf(pFile,"%6f ",1.0+normSysts["pdf_VV_wz_7TeV"]);}
            	    }else{fprintf(pFile,"%6s ","-");}
                     }fprintf(pFile,"\n");
            */


            ///////////////////////////////////////////////
            // RJ, for count and cut
            ///////////////////////////////////////////////

            fprintf(pFile,"%35s %10s ", "pdf_qqbar", "lnN");
            for(size_t j=1; j<=dci.procs.size(); j++) {
                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                if(dci.procs[j-1].BeginsWith("zh")) {
                    if(systpostfix.Contains('8')) {
                        fprintf(pFile,"%6f ",1.055);
                    } else {
                        fprintf(pFile,"%6f ",1.055);
                    }
                } else if(dci.procs[j-1].BeginsWith("zz")) {
                    if(systpostfix.Contains('8')) {
                        fprintf(pFile,"%6f ",1.057);
                    } else {
                        fprintf(pFile,"%6f ",1.057);
                    }
                } else if(dci.procs[j-1].BeginsWith("wz")) {
                    if(systpostfix.Contains('8')) {
                        fprintf(pFile,"%6f ",1.048);
                    } else {
                        fprintf(pFile,"%6f ",1.048);
                    }
                } else {
                    fprintf(pFile,"%6s ","-");
                }
            }
            fprintf(pFile,"\n");

            fprintf(pFile,"%35s %10s ", "QCDScale_VV", "lnN");
            for(size_t j=1; j<=dci.procs.size(); j++) {
                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                if(dci.procs[j-1].BeginsWith("zz")) {
                    if(systpostfix.Contains('8')) {
                        fprintf(pFile,"%6f ",1.067);
                    } else {
                        fprintf(pFile,"%6f ",1.067);
                    }
                } else if(dci.procs[j-1].BeginsWith("wz")) {
                    if(systpostfix.Contains('8')) {
                        fprintf(pFile,"%6f ",1.077);
                    } else {
                        fprintf(pFile,"%6f ",1.077);
                    }
                } else {
                    fprintf(pFile,"%6s ","-");
                }
            }
            fprintf(pFile,"\n");

            fprintf(pFile,"%35s %10s ", "QCDScale_ZH", "lnN");
            for(size_t j=1; j<=dci.procs.size(); j++) {
                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                if(dci.procs[j-1].BeginsWith("zh")) {
                    if(dci.ch[i-1].Contains("eq0jets")) {
                        if(systpostfix.Contains('8')) {
                            fprintf(pFile,"%6f ",1.072);
                        } else {
                            fprintf(pFile,"%6f ",1.072);
                        }
                    } else if(dci.ch[i-1].Contains("eq1jets")) {
                        if(systpostfix.Contains('8')) {
                            fprintf(pFile,"%6f ",1.105);
                        } else {
                            fprintf(pFile,"%6f ",1.105);
                        }
                    }
                } else {
                    fprintf(pFile,"%6s ","-");
                }
            }
            fprintf(pFile,"\n");

            /*
            	 if(dci.ch[i-1].Contains("ee"))
                     	fprintf(pFile,"%35s %10s ", "scale_e", "lnN");
            	 else fprintf(pFile,"%35s %10s ", "scale_mu", "lnN");
                     for(size_t j=1; j<=dci.procs.size(); j++) { if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            	          if(dci.procs[j-1].BeginsWith("zh")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.01);}else{fprintf(pFile,"%6f ",1.01);}}
            		  else if(dci.procs[j-1].BeginsWith("zz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.01);}else{fprintf(pFile,"%6f ",1.01);}}
            		  else if(dci.procs[j-1].BeginsWith("zll") && !dci.procs[j-1].BeginsWith("zlldata")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.01);}else{fprintf(pFile,"%6f ",1.01);}}
            		  else if(dci.procs[j-1].BeginsWith("wz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.01);}else{fprintf(pFile,"%6f ",1.01);}}
            		  else {fprintf(pFile,"%6s ","-");}
                     }fprintf(pFile,"\n");
            */

//RJ May 17
            /*
                     fprintf(pFile,"%35s %10s ", "jes", "lnN");
                     for(size_t j=1; j<=dci.procs.size(); j++) { if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            	          if(dci.procs[j-1].BeginsWith("zh")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.02);}else{fprintf(pFile,"%6f ",1.02);}}
            		  else if(dci.procs[j-1].BeginsWith("zz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.02);}else{fprintf(pFile,"%6f ",1.02);}}
            		  else if(dci.procs[j-1].BeginsWith("wz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.03);}else{fprintf(pFile,"%6f ",1.03);}}
            		  else {fprintf(pFile,"%6s ","-");}
                     }fprintf(pFile,"\n");
            */

            /*
            	//remove it if use data-driven DY sample
                     fprintf(pFile,"%35s %10s ", "err_zll", "lnN");
                     for(size_t j=1; j<=dci.procs.size(); j++) { if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            	          if(dci.procs[j-1].BeginsWith("zll")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",2.0);}else{fprintf(pFile,"%6f ",2.0);}}
            		  else {fprintf(pFile,"%6s ","-");}
                     }fprintf(pFile,"\n");
            */
            ///////////////////////////////////////////////








            ///////////////////////////////////////////////
            /*
                fprintf(pFile,"%35s %10s ", "XSec_sys_WW", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                   if(dci.procs[j-1].BeginsWith("ww")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.097);}else{fprintf(pFile,"%6f ",1.097);}
                   }else{fprintf(pFile,"%6s ","-");}
                }fprintf(pFile,"\n");

                fprintf(pFile,"%35s %10s ", "XSec_sys_WZ", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                   if(dci.procs[j-1].BeginsWith("wz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.056);}else{fprintf(pFile,"%6f ",1.056);}
                   }else{fprintf(pFile,"%6s ","-");}
                }fprintf(pFile,"\n");
            */

            ///////////////////////////////////////////////


            for(std::map<TString, std::map<RateKey_t,Double_t> >::iterator it=dci.systs.begin(); it!=dci.systs.end(); it++) {
                if(!runSystematics && string(it->first.Data()).find("stat")>0 )continue;

                isSyst=false;
                if(it->first.Contains("_sys_") || it->first.Contains("_interpol_")) {
                    sprintf(sFile,"%35s %10s ", it->first.Data(), "lnN");
                    for(size_t j=1; j<=dci.procs.size(); j++) {
                        if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                        if(it->second.find(RateKey_t(dci.procs[j-1],dci.ch[i-1])) != it->second.end()) {
                            Double_t systUnc = it->second[RateKey_t(dci.procs[j-1],dci.ch[i-1])];
                            if(systUnc<=0) {
                                sprintf(sFile,"%s%6s ",sFile,"-");
                            } else {
                                sprintf(sFile,"%s%6f ",sFile,(1.0+ (systUnc / dci.rates[RateKey_t(dci.procs[j-1],dci.ch[i-1])]) ));
                                isSyst=true;
                            }
                        } else {
                            sprintf(sFile,"%s%6s ",sFile,"-");
                        }
                    }
                    if(isSyst)fprintf(pFile,"%s\n",sFile);

                } else {
                    if(shape) {
                        sprintf(sFile,"%35s %10s ", it->first.Data(), "shapeN2");
                    } else {
                        sprintf(sFile,"%35s %10s ", it->first.Data(), "lnN");
                    }
                    for(size_t j=1; j<=dci.procs.size(); j++) {
                        if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                        if(it->first.Contains("sherpa")) continue; //RJ
                        //if(it->first.Contains("CMS_scale_j")) continue; //RJ
                        //if(it->first.Contains("CMS_res_j")) continue; //RJ
                        if(it->first.Contains("CMS_zh2l2v_les") && dci.ch[i-1].Contains("mumu")) continue; //RJ

                        if(it->second.find(RateKey_t(dci.procs[j-1],dci.ch[i-1])) != it->second.end()) {
                            if(it->first.Contains("sherpa") && (!dci.procs[j-1].Contains("zz2l2nu"))) {
                                sprintf(sFile,"%s%6s ",sFile,"-");
                            } else {
                                sprintf(sFile,"%s%6.3f ",sFile,it->second[RateKey_t(dci.procs[j-1],dci.ch[i-1])]);
                                isSyst=true;
                            }
                        } else {
                            sprintf(sFile,"%s%6s ",sFile,"-");
                        }
                    }
                    if(isSyst)fprintf(pFile,"%s\n",sFile);
                }
            }
        }


        fclose(pFile);
        cout << "Data card for " << dci.shapesFile << " and " << dci.ch[i-1] << " channel @ " << dcName << endl;
        dcUrls.push_back(dcName);
    }

    FILE* pFile = fopen("combineCards.sh","w");
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + combinedcard + " > " + "card_combined.dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + combinedcardLL + " > " + "card_combinedLL.dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + eecard       + " > " + "card_ee.dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + mumucard     + " > " + "card_mumu.dat").Data());
    fclose(pFile);

//  gSystem->Exec(TString("combineCards.py ") + combinedcard + " > " + outDir+"/card_combined.dat");
//  gSystem->Exec(TString("combineCards.py ") + eecard       + " > " + outDir+"/card_ee.dat");
//  gSystem->Exec(TString("combineCards.py ") + mumucard     + " > " + outDir+"/card_mumu.dat");

    //all done
    return dcUrls;
}

//
DataCardInputs convertHistosForLimits(Int_t mass,TString histo,TString url,TString Json)
{
    return convertHistosForLimits("",mass,histo,url,Json);
}

DataCardInputs convertHistosForLimits(TString atgcpar,Int_t mass,TString histo,TString url,TString Json)
{
    DataCardInputs dci;

    //init the json wrapper
    JSONWrapper::Object Root(Json.Data(), true);

    //init globalVariables
    TString massStr("");
    if(mass>0)massStr += mass;
    if(atgcpar.Length()>0) massStr = atgcpar;
    TString massStr2 = atgcpar2;
    if(massStr2.CompareTo("")==0) massStr2 = massStr;
//  std::set<TString> allCh,allProcs;
    std::vector<TString> allCh,allProcs;

    //open input file
    TFile* inF = TFile::Open(url);
    if( !inF || inF->IsZombie() ) {
        cout << "Invalid file name : " << url << endl;
        return dci;
    }

    //get the shapes for each channel
    map<TString, Shape_t> allShapes;
    map<TString, Shape_t> allShapesL;
    map<TString, Shape_t> allShapesR;
    TString ch[]= {"mumu","ee","emu","ll"};
    const size_t nch=sizeof(ch)/sizeof(TString);
    std::vector<TString> sh;
    sh.push_back(histo);
//  if(subNRB2011 || subNRB2012)sh.push_back("nonresbckg_ctrl");
    if(subNRB2011 || subNRB2012)sh.push_back(histo+"_NRBctrl");
    if(subNRB2012)sh.push_back(histo+"BTagSB");
    if(subWZ)sh.push_back(histo+"_3rdLepton");
    const size_t nsh=sh.size();
    for(size_t i=0; i<nch; i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            int indexcut_ = indexcut;
            double cutMin=shapeMin;
            double cutMax=shapeMax;
            if(indexvbf>=0 && AnalysisBins[b].Contains("vbf")) {
                printf("use vbf index(%i) for bin %s\n", indexvbf, AnalysisBins[b].Data());
                indexcut_ = indexvbf;
                cutMin=shapeMinVBF;
                cutMax=shapeMaxVBF;
            }
            for(size_t j=0; j<nsh; j++) {
                printf("i=%i b=%i j=%i\n",(int)i,(int)b,(int)j);
                allShapes[ch[i]+AnalysisBins[b]+sh[j]]=getShapeFromFile(inF, ch[i]+AnalysisBins[b],sh[j],indexcut_,Root,cutMin, cutMax);
                cout << "allShapes: " << ch[i]+AnalysisBins[b]+sh[j] << endl; //RJ
                if(indexcutL>=0 && indexcutR>=0) {
                    if(indexvbf>=0 && AnalysisBins[b].Contains("vbf")) {
                        allShapesL[ch[i]+AnalysisBins[b]+sh[j]]=getShapeFromFile(inF, ch[i]+AnalysisBins[b],sh[j],indexvbf,Root,cutMin, cutMax);
                        allShapesR[ch[i]+AnalysisBins[b]+sh[j]]=getShapeFromFile(inF, ch[i]+AnalysisBins[b],sh[j],indexvbf,Root,cutMin, cutMax);
                    } else {
                        allShapesL[ch[i]+AnalysisBins[b]+sh[j]]=getShapeFromFile(inF, ch[i]+AnalysisBins[b],sh[j],indexcutL,Root,cutMin, cutMax);
                        allShapesR[ch[i]+AnalysisBins[b]+sh[j]]=getShapeFromFile(inF, ch[i]+AnalysisBins[b],sh[j],indexcutR,Root,cutMin, cutMax);
                    }
                }
            }
        }
    }

    //all done with input file
    inF->Close();

    printf("done loading all shapes\n");

    //define vector for search
    std::vector<TString>& selCh = Channels;
//  selCh.push_back("ee"); selCh.push_back("mumu");

    //non-resonant background estimation
    //estimateNonResonantBackground(selCh,"emu",allShapes,"nonresbckg_ctrl");

    //remove the non-resonant background from data
    if(subNRB2011 || subNRB2012)doBackgroundSubtraction(selCh,"emu",allShapes,histo,histo+"_NRBctrl", url, Root);

    //replace WZ by its estimate from 3rd Lepton SB
    if(subWZ)doWZSubtraction(selCh,"emu",allShapes,histo,histo+"_3rdLepton");

    //replace Z+Jet background by Gamma+Jet estimates
    if(subDY)doDYReplacement(selCh,"gamma",allShapes,histo,"met_met");
    //if(subDY)doDYReplacement(selCh,"gamma",allShapes,histo,"met_redMet");

    //replace data by total MC background
    if(blindData)BlindData(selCh,allShapes,histo, blindWithSignal);

    //interpollate signal sample if desired mass point is not available
    SignalInterpolation(selCh,allShapesL, allShapes, allShapesR, histo);


    if(doInterf)RescaleForInterference(selCh,allShapes, histo);



    //print event yields from the mt shapes
    if(!fast)getYieldsFromShape(selCh,allShapes,histo);
    if(!fast)getEffFromShape(selCh,allShapes,histo);

    if(!fast)showShape(selCh,allShapes,histo,"plot");


    //prepare the output
    dci.shapesFile="zh2l2v_"+massStr2+systpostfix+".root";
    TFile *fout=TFile::Open(dci.shapesFile,"recreate");

    //loop on channel/proc/systematics
    for(size_t ich=0; ich<selCh.size(); ich++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            TString chbin = selCh[ich]+AnalysisBins[b];
            fout->mkdir(chbin);
            fout->cd(chbin);
            allCh.push_back(chbin);
            Shape_t shapeSt = allShapes.find(chbin+histo)->second;

            //signals
            dci.nsignalproc = 0;
            size_t nsignal=allShapes.find(chbin+histo)->second.signal.size();
            for(size_t isignal=0; isignal<nsignal; isignal++) {
                TH1* h=shapeSt.signal[isignal];

                TString proc(h->GetTitle());
                if(mass>0 && !proc.Contains(massStr))continue;
                if(mass>0 && proc.Contains("ggH") && proc.Contains("ZZ"))proc = "ggHZZ2l2v";
                else if(mass>0 && proc.Contains("qqH") && proc.Contains("ZZ"))proc = "qqHZZ2l2v";
                else if(mass>0 && proc.Contains("ggH") && proc.Contains("WW"))proc = "ggHWW2l2v";
                else if(mass>0 && proc.Contains("qqH") && proc.Contains("WW"))proc = "qqHWW2l2v";
                else if(mass>0 && proc.Contains("ZH")                        )proc = "ZH"+massStr+"2lMET";

                //if(atgcpar.Length()>0 && !proc.Contains(massStr)) continue;
                if(atgcpar.Length()>0 && !proc.EndsWith(massStr)) continue;
                if(atgcpar.Length()>0) proc = "ZZ2l2nu_"+massStr2;

                std::vector<std::pair<TString, TH1*> > vars = shapeSt.signalVars[h->GetTitle()];
                std::vector<TString> systs;
                std::vector<TH1*>    hshapes;
                systs.push_back("");
                hshapes.push_back(shapeSt.signal[isignal]);
                for(size_t v=0; v<vars.size(); v++) {
                    //printf("SYSTEMATIC FOR SIGNAL %s : %s\n",h->GetTitle(), vars[v].first.Data());
                    systs.push_back(vars[v].first);
                    hshapes.push_back(vars[v].second);
                }

                convertHistosForLimits_core(dci, proc, AnalysisBins[b], chbin, systs, hshapes);
                if(ich==0 && b==0)allProcs.push_back(proc);
                dci.nsignalproc++;
            }

            //backgrounds
            size_t nbckg=allShapes.find(chbin+histo)->second.bckg.size();
            size_t nNonNullBckg=0;
            for(size_t ibckg=0; ibckg<nbckg; ibckg++) {
                TH1* h=shapeSt.bckg[ibckg];
                std::vector<std::pair<TString, TH1*> > vars = shapeSt.bckgVars[h->GetTitle()];

                std::vector<TString> systs;
                std::vector<TH1*>    hshapes;
                systs.push_back("");
                hshapes.push_back(shapeSt.bckg[ibckg]);
                for(size_t v=0; v<vars.size(); v++) {
                    systs.push_back(vars[v].first);
                    hshapes.push_back(vars[v].second);
                }

                TString proc(h->GetTitle());
                convertHistosForLimits_core(dci, proc, AnalysisBins[b], chbin, systs, hshapes);
                if(ich==0 && b==0)allProcs.push_back(proc);

                //remove backgrounds with rate=0 (but keep at least one background)
                std::map<RateKey_t, Double_t>::iterator it = dci.rates.find(RateKey_t(proc,chbin));
                if(it==dci.rates.end()) {
                    printf("proc=%s not found --> THIS SHOULD NEVER HAPPENS.  PLEASE CHECK THE CODE\n",proc.Data());
                } else {
                    if(it->second>0) {
                        nNonNullBckg++;
                    } else if(ibckg<nbckg-1 || nNonNullBckg>0) {
                        dci.rates.erase(dci.rates.find(RateKey_t(proc,chbin)));
                    }
                }

            }
            /*
                 //remove all backgrounds with rate=0
                 size_t nNonNullBckg=0;
                 for(size_t ibckg=0; ibckg<nbckg; ibckg++){
                    TString proc(shapeSt.bckg[ibckg]->GetTitle());
               proc.ReplaceAll("#bar{t}","tbar");
               proc.ReplaceAll("Z-#gamma^{*}+jets#rightarrow ll","dy");
               proc.ReplaceAll("(","");    proc.ReplaceAll(")","");    proc.ReplaceAll("+","");    proc.ReplaceAll(" ","");
               proc.ToLower();

                    std::map<RateKey_t, Double_t>::iterator it = dci.rates.find(RateKey_t(proc,chbin));
                    if(it==dci.rates.end()){printf("proc=%s not found\n",proc.Data());continue;}
                    printf("Rates for %s-%s is %6.2E  - NonNull=%i\n",chbin.Data(), proc.Data(), it->second,(int) nNonNullBckg);
                    if(it->second>0){
                       nNonNullBckg++;
            	   printf("nNonNullBckg++\n");
                    }else{
                       printf("ibckg<nbckg-1=%i || nNonNullBckg=%i --> %i\n",(int)(ibckg<nbckg-1), (int)(nNonNullBckg>0), (int)(ibckg<nbckg-1 || nNonNullBckg>0));
            //           if(ibckg<nbckg-1 || nNonNullBckg>0){printf("Erase Rate %i --> ", (int)dci.rates.size());  dci.rates.erase(dci.rates.find(RateKey_t(proc,chbin))); printf("%i\n", (int)dci.rates.size());}
                       printf("Erase Rate %i --> ", (int)dci.rates.size());  dci.rates.erase(it); printf("%i\n", (int)dci.rates.size());
                    }
                 }
            */

            //data
            TH1* h=shapeSt.data;
            std::vector<TString> systs;
            std::vector<TH1*>    hshapes;
            systs.push_back("");
            hshapes.push_back(h);
            TString proc(h->GetTitle());
            convertHistosForLimits_core(dci, proc, AnalysisBins[b], chbin, systs, hshapes);

            //return to parent dir
            fout->cd("..");
        }
    }
    dci.ch.resize(allCh.size());
    std::copy(allCh.begin(), allCh.end(),dci.ch.begin());
    dci.procs.resize(allProcs.size());
    std::copy(allProcs.begin(), allProcs.end(),dci.procs.begin());


// DEBUGGING
//  printf("-----------------------\n");
//  printf("shapesFile=%s\n",dci.shapesFile.Data());
//  for(unsigned int i=0;i<dci.ch.size();i++){printf("%s - ",dci.ch[i].Data());}printf("\n");
//  for(unsigned int i=0;i<dci.procs.size();i++){printf("%s - ",dci.procs[i].Data());}printf("\n");
//  for(std::map<TString, std::map<RateKey_t,Double_t> >::iterator iter = dci.systs.begin();   iter != dci.systs.end(); ++iter ){
//       printf("%s : ", iter->first.Data());
//       for(std::map<RateKey_t, Double_t>::iterator it = iter->second.begin();   it != iter->second.end(); ++it ){
//                 printf("%s_%s (%f) ", it->first.first.Data(), it->first.second.Data(), it->second);
//       }
//       printf("\n");
//  }
//  printf("-----------------------\n");


    /*
      //################# START BACKGROUND SUBTRACTION CODE

        TH1* proj_em = ((TH2*)fin->Get(TString("data/emu_" ) + "nonresbckg_ctrl" ))->ProjectionY("_py", indexcut,indexcut);
        TH1* proj_mm = ((TH2*)fin->Get(TString("data/mumu_") + "nonresbckg_ctrl" ))->ProjectionY("_py", indexcut,indexcut);
        TH1* proj_ee = ((TH2*)fin->Get(TString("data/ee_"  ) + "nonresbckg_ctrl" ))->ProjectionY("_py", indexcut,indexcut);

        printf("Bin %f %f %f %f %f %f\n", proj_em->GetBinContent(1), proj_em->GetBinContent(2), proj_em->GetBinContent(3), proj_em->GetBinContent(4), proj_em->GetBinContent(5), proj_em->GetBinContent(6) );
        printf("Bin %f %f %f %f %f %f\n", proj_ee->GetBinContent(1), proj_ee->GetBinContent(2), proj_ee->GetBinContent(3), proj_ee->GetBinContent(4), proj_ee->GetBinContent(5), proj_ee->GetBinContent(6) );
        printf("Bin %f %f %f %f %f %f\n", proj_mm->GetBinContent(1), proj_mm->GetBinContent(2), proj_mm->GetBinContent(3), proj_mm->GetBinContent(4), proj_mm->GetBinContent(5), proj_mm->GetBinContent(6) );
        double alpha_e     = proj_ee->GetBinContent(5) / proj_em->GetBinContent(5);
        double alpha_e_err = ( fabs( proj_ee->GetBinContent(5) * proj_em->GetBinError(5) ) + fabs(proj_ee->GetBinError(5) * proj_em->GetBinContent(5) )  ) / pow(proj_em->GetBinContent(5), 2);

        double alpha_m     = proj_mm->GetBinContent(5) / proj_em->GetBinContent(5);
        double alpha_m_err = ( fabs( proj_mm->GetBinContent(5) * proj_em->GetBinError(5) ) + fabs(proj_mm->GetBinError(5) * proj_em->GetBinContent(5) )  ) / pow(proj_em->GetBinContent(5), 2);

        printf("alpha e=%f+-%f mu=%f+-%f\n", alpha_e, alpha_e_err, alpha_m, alpha_m_err);





      //################# END   BACKGROUND SUBTRACTION CODE

    */


    //all done
    fout->Close();

    return dci;
}



void convertHistosForLimits_core(DataCardInputs& dci, TString& proc, TString& bin, TString& ch, std::vector<TString>& systs, std::vector<TH1*>& hshapes)
{
    proc.ReplaceAll("#bar{t}","tbar");
    proc.ReplaceAll("Z-#gamma^{*}+jets#rightarrow ll","dy");
    proc.ReplaceAll("#rightarrow","");
    proc.ReplaceAll("(","");
    proc.ReplaceAll(")","");
    proc.ReplaceAll("+","");
    proc.ReplaceAll(" ","");
    proc.ReplaceAll("/","");
    proc.ReplaceAll("#","");
    proc.ReplaceAll("=","");
    proc.ReplaceAll(".","");
    proc.ReplaceAll("^","");
    proc.ReplaceAll("}","");
    proc.ReplaceAll("{","");
    proc.ToLower();

    for(unsigned int i=0; i<systs.size(); i++) {
        TString syst   = systs[i];
        TH1*    hshape = hshapes[i];
        hshape->SetDirectory(0);

        //Do Renaming and cleaning
        syst.ReplaceAll("down","Down");
        syst.ReplaceAll("up","Up");

        if(syst=="") {
            syst="";
        } else if(syst.BeginsWith("_jes")) {
            syst.ReplaceAll("_jes","_CMS_scale_j");
        } else if(syst.BeginsWith("_jer")) {
            syst.ReplaceAll("_jer","_CMS_res_j"); // continue;//skip res for now
        } else if(syst.BeginsWith("_btag")) {
            syst.ReplaceAll("_btag","_CMS_eff_b");
        } else if(syst.BeginsWith("_pu" )) {
            syst.ReplaceAll("_pu", "_CMS_zh2l2v_pu");
        } else if(syst.BeginsWith("_ren" )) {
            continue;   //already accounted for in QCD scales
        } else if(syst.BeginsWith("_fact" )) {
            continue; //skip this one
        } else if(syst.BeginsWith("_interf" )) {
        } else {
            syst="_CMS_zh2l2v"+syst;
        }

        double systUncertainty = hshape->GetBinError(0);
        hshape->SetBinError(0,0);

        //If cut&count keep only 1 bin in the histo
        if(!shape) {
//          hshape = hshape->Rebin(hshape->GetXaxis()->GetNbins(), TString(hshape->GetName())+"_Rebin");
            hshape = hshape->Rebin(hshape->GetXaxis()->GetNbins());
            //make sure to also count the underflow and overflow
            double bin  = hshape->GetBinContent(0) + hshape->GetBinContent(1) + hshape->GetBinContent(2);
            double bine = sqrt(hshape->GetBinError(0)*hshape->GetBinError(0) + hshape->GetBinError(1)*hshape->GetBinError(1) + hshape->GetBinError(2)*hshape->GetBinError(2));
            hshape->SetBinContent(0,0);
            hshape->SetBinError  (0,0);
            hshape->SetBinContent(1,bin);
            hshape->SetBinError  (1,bine);
            hshape->SetBinContent(2,0);
            hshape->SetBinError  (2,0);
        }

        if(syst=="") {
            //central shape (for data call it data_obs)
            hshape->SetName(proc);
            if(proc=="data") {
                hshape->Write("data_obs");
            } else {
                hshape->Write(proc+postfix);

                if(hshape->Integral()>0) {
                    hshape->SetName(proc+syst);
                    TH1* statup=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"Up");
                    TH1* statdown=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"Down");
                    /*	       //RENJIE
                    	       if(proc=="topwwwjetsdata" || proc=="zlldata") {
                    			statup=(TH1 *)hshape->Clone(proc+"_stat"+proc+"Up");
                    			statdown=(TH1 *)hshape->Clone(proc+"_stat"+proc+"Down");
                    	       }
                    */
                    for(int ibin=1; ibin<=statup->GetXaxis()->GetNbins(); ibin++) {
                        statup  ->SetBinContent(ibin,std::min(2*hshape->GetBinContent(ibin), std::max(0.01*hshape->GetBinContent(ibin), statup  ->GetBinContent(ibin) + statup  ->GetBinError(ibin))));
                        statdown->SetBinContent(ibin,std::min(2*hshape->GetBinContent(ibin), std::max(0.01*hshape->GetBinContent(ibin), statdown->GetBinContent(ibin) - statdown->GetBinError(ibin))));
                    }
                    //RENJIE
//	       if(proc=="topwwwjetsdata" || proc=="zlldata") {
//               	statup  ->Write(proc+postfix+"_CMS_zh2l2v_stat_"+proc+systpostfix+"Up");
//               	statdown->Write(proc+postfix+"_CMS_zh2l2v_stat_"+proc+systpostfix+"Down");
//	       }else{
                    statup  ->Write(proc+postfix+"_CMS_zh2l2v_stat_"+ch+"_"+proc+systpostfix+"Up");
                    statdown->Write(proc+postfix+"_CMS_zh2l2v_stat_"+ch+"_"+proc+systpostfix+"Down");
//	       }


                    if(shape) {
                        dci.systs["CMS_zh2l2v_stat_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=1.0;
                    } else {
                        dci.systs["CMS_zh2l2v_stat_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                    }

                    /*
                                   if(shape){ //RENJIE  Jun15
                    			if(proc=="topwwwjetsdata" || proc=="zlldata")
                    				dci.systs["CMS_zh2l2v_stat_"+proc+systpostfix][RateKey_t(proc,ch)]=1.0;
                    			else
                    				dci.systs["CMS_zh2l2v_stat_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=1.0;
                                   }else{
                    			if(proc=="topwwwjetsdata" || proc=="zlldata")
                    				dci.systs["CMS_zh2l2v_stat_"+proc+systpostfix][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                    			else
                    				dci.systs["CMS_zh2l2v_stat_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                                   }
                    */

                    if(systUncertainty>0) {
                        if(proc.Contains("ggh") || proc.Contains("qqh")) {
                            dci.systs["CMS_zh2l2v_interpol_"+bin+"_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;
                        } else {
                            printf("SYST in %s - %s - %s = %f\n",bin.Data(), ch.Data(), proc.Data(), systUncertainty);
                            //makesure that syst+stat error is never bigger than 100%
                            //double valerr, val  = hshape->IntegralAndError(1,hshape->GetXaxis()->GetNbins(),valerr);
                            //if(sqrt(pow(valerr,2)+pow(systUncertainty,2))>val){systUncertainty = sqrt(std::max(0.0, pow(val,2) - pow(valerr,2)));}

                            //add syst uncertainty as bin dependent or not
//                     dci.systs["CMS_zh2l2v_sys_"+bin+"_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;
                            dci.systs["CMS_zh2l2v_sys_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;
                        }
                    }
                }
            }
        } else if(runSystematics && proc!="data" && (syst.Contains("Up") || syst.Contains("Down"))) {
            //if empty histogram --> no variation is applied
            if(hshape->Integral()<hshapes[0]->Integral()*0.01 || isnan((float)hshape->Integral())) {
                hshape->Reset();
                hshape->Add(hshapes[0],1);
            }

            //write variation to file
            hshape->SetName(proc+syst);
            hshape->Write(proc+postfix+syst);
        } else if(runSystematics) {
            //for one sided systematics the down variation mirrors the difference bin by bin
            hshape->SetName(proc+syst);
            hshape->Write(proc+postfix+syst+"Up");
            TH1 *hmirrorshape=(TH1 *)hshape->Clone(proc+syst+"Down");
            for(int ibin=1; ibin<=hmirrorshape->GetXaxis()->GetNbins(); ibin++) {
//            double bin = hmirrorshape->GetBinContent(ibin);
                double bin = 2*hshapes[0]->GetBinContent(ibin)-hmirrorshape->GetBinContent(ibin);
                if(bin<0)bin=0;
                hmirrorshape->SetBinContent(ibin,bin);
            }
            if(hmirrorshape->Integral()<=0)hmirrorshape->SetBinContent(1, 1E-10);
            hmirrorshape->Write(proc+postfix+syst+"Down");
        }

        if(runSystematics && syst!="") {
            TString systName(syst);
            systName.ReplaceAll("Up","");
            systName.ReplaceAll("Down","");//  systName.ReplaceAll("_","");
            if(systName.First("_")==0)systName.Remove(0,1);


            TH1 *temp=(TH1*) hshape->Clone();
            temp->Add(hshapes[0],-1);
            if(temp->Integral()!=0 || syst.BeginsWith("_interf" )) {
                if(shape) {
                    dci.systs[systName][RateKey_t(proc,ch)]=1.0;
                } else {
                    double Unc = 1 + fabs(temp->Integral()/hshapes[0]->Integral());
                    if(dci.systs.find(systName)==dci.systs.end() || dci.systs[systName].find(RateKey_t(proc,ch))==dci.systs[systName].end() ) {
                        dci.systs[systName][RateKey_t(proc,ch)]=Unc;
                    } else {
                        dci.systs[systName][RateKey_t(proc,ch)]=(dci.systs[systName][RateKey_t(proc,ch)] + Unc)/2.0;
                    }
                }
            }


            delete temp;
//        }else if(proc=="asignal" && syst==""){dci.rates[RateKey_t(proc,ch)]=hshape->Integral();
//        }else if(proc!="data" && syst==""){if(hshape->Integral()>1E-6)dci.rates[RateKey_t(proc,ch)]=hshape->Integral();
        } else if(proc!="data" && syst=="") {
            dci.rates[RateKey_t(proc,ch)]= hshape->Integral()>1E-6 ? hshape->Integral() : 0.0;
        } else if(proc=="data" && syst=="") {
            dci.obs[RateKey_t("obs",ch)]=hshape->Integral();
        }
    }
}

void doBackgroundSubtraction(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString sideBandHisto, TString url, JSONWrapper::Object &Root)
{
    string Lcol   = "\\begin{tabular}{|l";
    string Lchan  = "channel";
    string Lalph1 = "$\\alpha$ measured";
    string Lalph2 = "$\\alpha$ used";
    string Lyield   = "yield data";
    string LyieldMC = "yield mc";

    string Ccol   = "\\begin{tabular}{|l|c|c|c|c|";
    string Cname  = "channel & $\\alpha$ measured & $\\alpha$ used & yields $e\\mu$ & yield data & yield mc";
    string Cval   = "";
    FILE* pFile = NULL;
    if(!fast) {
        pFile = fopen("NonResonnant.tex","w");
        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Non resonant background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
        fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
    }

    for(size_t i=0; i<selCh.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Lcol += " |c";
            Lchan += string(" &")+selCh[i]+string(" - ")+AnalysisBins[b];
            Cval   = selCh[i]+string(" - ")+AnalysisBins[b];


            Shape_t& shapeCtrl_SB = allShapes.find(ctrlCh+AnalysisBins[b]+sideBandHisto)->second;
            TH1* hCtrl_SB=shapeCtrl_SB.data;
            Shape_t& shapeCtrl_SI = allShapes.find(ctrlCh+AnalysisBins[b]+mainHisto)->second;
            TH1* hCtrl_SI=shapeCtrl_SI.data;
            Shape_t& shapeChan_SB = allShapes.find(selCh[i]+AnalysisBins[b]+sideBandHisto)->second;
            TH1* hChan_SB=shapeChan_SB.data;
            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
//        TH1* hChan_SI=shapeChan_SI.data;

            /*
                    //IF HISTO IS EMPTY... LOWER THE CUT AND TAKE THIS WITH 100% UNCERTAINTY
                    if(subNRB2011 && hCtrl_SI->Integral()<=0){
                       printf("LOWERING THE CUTS for %s\n",(selCh[i]+string(" - ")+AnalysisBins[b]).Data());
                       TFile* inF = TFile::Open(url);
                       if( !inF || inF->IsZombie() ){break;}

                       int indexcut_ = indexcut; double cutMin=shapeMin; double cutMax=shapeMax; double cutMin_;
                       if(indexvbf>=0 && AnalysisBins[b].Contains("vbf")){printf("use vbf index(%i) for bin %s\n", indexvbf, AnalysisBins[b].Data()); indexcut_ = indexvbf; cutMin=shapeMinVBF; cutMax=shapeMaxVBF;}
                       while(hCtrl_SI->Integral()<=0 && indexcut_>=0){
                          cutMin_=cutMin;
                          while(hCtrl_SI->Integral()<=0 && cutMin_>=0){
                             hCtrl_SB=getShapeFromFile(inF, ctrlCh+AnalysisBins[b],sideBandHisto,indexcut_,Root,cutMin_, cutMax, true).data;
                             hCtrl_SI=getShapeFromFile(inF, ctrlCh+AnalysisBins[b],mainHisto,indexcut_,Root,cutMin_, cutMax, true).data;
                             hChan_SB=getShapeFromFile(inF, selCh[i]+AnalysisBins[b],sideBandHisto,indexcut_,Root,cutMin_, cutMax, true).data;
                             //printf("indexcut_ = %i cutMin=%f --> Integral = %f\n",indexcut_, cutMin_, hCtrl_SI->Integral());
            //                 cutMin_-=5;
                             cutMin_-=50;
                          }
                          indexcut_--;
                       }
                       inF->Close();

                       printf("Cuts are %i (>%f) --> EMU events=\n", indexcut_, cutMin_, hCtrl_SI->Integral());


                       //set stat error to 100%
                       for(int b=1;b<=hCtrl_SB->GetXaxis()->GetNbins()+1;b++){
                          hCtrl_SB->SetBinContent(b, hCtrl_SB->GetBinContent(b) *0.5 );
                          hCtrl_SI->SetBinContent(b, hCtrl_SI->GetBinContent(b) *0.5 );
                          hChan_SB->SetBinContent(b, hChan_SB->GetBinContent(b) *0.5 );
                          hCtrl_SB->SetBinError  (b, hCtrl_SB->GetBinContent(b) );
                          hCtrl_SI->SetBinError  (b, hCtrl_SI->GetBinContent(b) );
                          hChan_SB->SetBinError  (b, hChan_SB->GetBinContent(b) );
                       }
                    }
            */


            //printf("Channel = %s\n", selCh[i].Data());
            //printf("Bin %f %f %f %f %f %f\n", hCtrl_SB->GetBinContent(1), hCtrl_SB->GetBinContent(2), hCtrl_SB->GetBinContent(3), hCtrl_SB->GetBinContent(4), hCtrl_SB->GetBinContent(5), hCtrl_SB->GetBinContent(6) );
            //printf("Bin %f %f %f %f %f %f\n", hChan_SB->GetBinContent(1), hChan_SB->GetBinContent(2), hChan_SB->GetBinContent(3), hChan_SB->GetBinContent(4), hChan_SB->GetBinContent(5), hChan_SB->GetBinContent(6) );
            double alpha=0 ,alpha_err=0;
            if(hCtrl_SB->GetBinContent(5)>0) {
                alpha     = hChan_SB->GetBinContent(5) / hCtrl_SB->GetBinContent(5);
                alpha_err = ( fabs( hChan_SB->GetBinContent(5) * hCtrl_SB->GetBinError(5) ) + fabs(hChan_SB->GetBinError(5) * hCtrl_SB->GetBinContent(5) )  ) / pow(hCtrl_SB->GetBinContent(5), 2);
                //alpha     = hChan_SB->GetBinContent(2) / hCtrl_SB->GetBinContent(2);
                //alpha_err = ( fabs( hChan_SB->GetBinContent(2) * hCtrl_SB->GetBinError(2) ) + fabs(hChan_SB->GetBinError(2) * hCtrl_SB->GetBinContent(2) )  ) / pow(hCtrl_SB->GetBinContent(2), 2);
            }

            Lalph1 += string(" &") + toLatexRounded(alpha,alpha_err);
            Cval   += string(" &") + toLatexRounded(alpha,alpha_err);

//        if(selCh[i].First("ee"  )!=kNPOS){alpha = 0.339286; alpha_err=0.043549;}
//        if(selCh[i].First("mumu")!=kNPOS){alpha = 0.529018; alpha_err=0.059357;}

            // Temporarily excluded (Daniele)
            //if(selCh[i].First("ee"  )!=kNPOS){alpha = 0.34; alpha_err=0.03;}
            //if(selCh[i].First("mumu")!=kNPOS){alpha = 0.61; alpha_err=0.04;}


            //add 100% syst uncertainty on alpha
            //alpha_err = sqrt(pow(alpha*1.0,2)+pow(alpha_err,2));

            Lalph2 += string(" &") + toLatexRounded(alpha,alpha_err);
            Cval   += string(" &") + toLatexRounded(alpha,alpha_err);


            TH1* NonResonant = NULL;
            if(subNRB2011) {
                NonResonant = (TH1*)hCtrl_SI->Clone("Top/WW/W+Jets (data)");
                NonResonant->SetTitle("Top/WW/W+Jets (data)");
            } else if(subNRB2012) {
                Shape_t& shapeChan_BTag = allShapes.find(selCh[i]+mainHisto+"BTagSB")->second;
                NonResonant = (TH1*)shapeChan_BTag.data->Clone("Top (data)");
                NonResonant->SetTitle("Top (data)");
            } else {
                return;
            }

            double valvalerr, valval;
            valval = NonResonant->IntegralAndError(1,NonResonant->GetXaxis()->GetNbins(),valvalerr);

            Cval   += string(" &") + toLatexRounded(valval,valvalerr);

            printf("VALVAL=%f\n",valval);
            for(int b=1; b<=NonResonant->GetXaxis()->GetNbins()+1; b++) {
                double val = NonResonant->GetBinContent(b);
                double err = NonResonant->GetBinError(b);
                double newval = val*alpha;
                double newerr = sqrt(pow(err*alpha,2) + pow(val*alpha_err,2));
                NonResonant->SetBinContent(b, newval );
                NonResonant->SetBinError  (b, newerr );
            }
            NonResonant->Scale(DDRescale);

            Double_t valerr;
            Double_t val = NonResonant->IntegralAndError(1,NonResonant->GetXaxis()->GetNbins(),valerr);
            Double_t systError = val*NonResonnantSyst;
            NonResonant->SetBinError(0,systError);//save syst error in underflow bin that is always empty
            if(val<1E-6) {
                val=0.0;
                valerr=0.0;
                systError=-1;
            }
            Lyield += string(" &") + toLatexRounded(val,valerr,systError);
            Cval   += string(" &") + toLatexRounded(val,valerr,systError);

            //Clean background collection
            TH1* MCNRB = (TH1*)shapeChan_SI.totalBckg->Clone("MCNRB");
            MCNRB->Reset();
            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                if(( subNRB2011 && (proc.Contains("t#bar{t}") || proc.Contains("Single top") || proc.Contains("WW") || proc.Contains("Z#rightarrow #tau#tau") || proc.Contains("W#rightarrow l#nu")) ) ||
                        ( subNRB2012 && (proc.Contains("t#bar{t}") || proc.Contains("Single top") ) ) ) {
                    MCNRB->Add(shapeChan_SI.bckg[ibckg], 1);
                    NonResonant->SetFillColor(shapeChan_SI.bckg[ibckg]->GetFillColor());
                    NonResonant->SetLineColor(shapeChan_SI.bckg[ibckg]->GetLineColor());
                    shapeChan_SI.bckg.erase(shapeChan_SI.bckg.begin()+ibckg);
                    ibckg--;
                }
            }

            NonResonant->SetFillStyle(1001);
            NonResonant->SetFillColor(393);//592);
            NonResonant->SetLineColor(393);//592);


            //add the background estimate
            shapeChan_SI.bckg.push_back(NonResonant);


            //recompute total background
            shapeChan_SI.totalBckg->Reset();
            for(size_t i=0; i<shapeChan_SI.bckg.size(); i++) {
                shapeChan_SI.totalBckg->Add(shapeChan_SI.bckg[i]);
            }

            val = MCNRB->IntegralAndError(1,MCNRB->GetXaxis()->GetNbins(),valerr);
            if(val<1E-6) {
                val=0.0;
                valerr=0.0;
            }
            LyieldMC += string(" &") + toLatexRounded(val,valerr);
            Cval     += string(" &") + toLatexRounded(val,valerr);

            if(pFile) {
                fprintf(pFile,"%s\\\\\n", Cval.c_str());
            }
        }
    }

    if(pFile) {
        fprintf(pFile,"\\hline\n");
        fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
        fprintf(pFile,"\n\n\n\n");

        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Non resonant background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s|}\\hline\n", Lcol.c_str());
        fprintf(pFile,"%s\\\\\n", Lchan.c_str());
        fprintf(pFile,"%s\\\\\n", Lalph1.c_str());
        fprintf(pFile,"%s\\\\\n", Lalph2.c_str());
        fprintf(pFile,"%s\\\\\n", Lyield.c_str());
        fprintf(pFile,"%s\\\\\n", LyieldMC.c_str());
        fprintf(pFile,"\\hline\n");
        fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
        fclose(pFile);
    }
}



void doDYReplacement(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString metHistoForRescale)
{
    TString DYProcName = "Z#rightarrow ll";
    TString GammaJetProcName = "Instr. background (data)";
    std::map<TString, double> LowMetIntegral;


    string Ccol   = "\\begin{tabular}{|l|c|c|c|";
    string Cname  = "channel & rescale & yield data & yield mc";
    string Cval   = "";
    FILE* pFile = NULL;
    if(!fast) {
        pFile = fopen("GammaJets.tex","w");
        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Instrumental background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
        fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
    }

    FILE* gjFile = NULL;
    gjFile = fopen("DataDY_Yields.log","w");


    //open input file
    TFile* inF = TFile::Open(inFileUrl);
    if( !inF || inF->IsZombie() ) {
        cout << "Invalid file name : " << inFileUrl << endl;
        return;
    }

    TDirectory *pdir = (TDirectory *)inF->Get(DYProcName);
    if(!pdir) {
        printf("Skip Z+Jet estimation because %s directory is missing in root file\n", DYProcName.Data());
        return;
    }

    /*
      //Just to redo what we did for the HIGHMASS pape
      for(size_t i=0;i<selCh.size();i++){
      for(size_t b=0; b<AnalysisBins.size(); b++){
         TH1* met = (TH1*)pdir->Get(selCh[i]+AnalysisBins[b]+"_"+metHistoForRescale);
         LowMetIntegral[selCh[i]+AnalysisBins[b]] = met->Integral(1,met->GetXaxis()->FindBin(50));
      }}
    */
    //all done with input file
    inF->Close();


    //open gamma+jet file
    inF = TFile::Open(DYFile);
    if( !inF || inF->IsZombie() ) {
        cout << "Invalid file name : " << DYFile << endl;
        return;
    }

    pdir = (TDirectory *)inF->Get(GammaJetProcName);
    if(!pdir) {
        printf("Skip Z+Jet estimation because %s directory is missing in Gamma+Jet file\n", GammaJetProcName.Data());
        return;
    }

    for(size_t i=0; i<selCh.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Cval   = selCh[i]+string(" - ")+AnalysisBins[b];
            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;

            //find DY background
            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                if( proc.Contains(DYProcName) ) {

                    std::cout<<"DYTEST >>>>>  " << selCh[i]+string(" - ")+AnalysisBins[b] << "  \n";

                    double rescaleFactor = 1.0;
                    /*
                               //compute rescale factor using low MET events
                               TH1* met = (TH1*)pdir->Get(selCh[i]+AnalysisBins[b]+"_"+metHistoForRescale);
                               double integral = met->Integral(1,met->GetXaxis()->FindBin(50));
                               double rescaleFactor = LowMetIntegral[selCh[i]+AnalysisBins[b]] / integral; //Just to redo what we did for the HIGHMASS paper
                               printf("Rescale in %s = %f/%f = %f\n",  (selCh[i]+AnalysisBins[b]).Data(), LowMetIntegral[selCh[i]+AnalysisBins[b]], integral, rescaleFactor);
                    */
                    char buffer[255];
                    sprintf(buffer,"%6.3f",rescaleFactor);
                    Cval   += string(" &") + buffer;

                    Double_t valerr, val = shapeChan_SI.bckg[ibckg]->IntegralAndError(1,shapeChan_SI.bckg[ibckg]->GetXaxis()->GetNbins(),valerr);
                    if(val<1E-6) {
                        val=0.0;
                        valerr=0.0;
                    }
                    string MCTruth = toLatexRounded(val,valerr);
                    double systError = val;

                    //replace DY histogram by G+Jets data
                    int indexcut_ = indexcut;
                    double cutMin=shapeMin;
                    double cutMax=shapeMax;
                    if(indexvbf>=0 && AnalysisBins[b].Contains("vbf")) {
                        printf("use vbf index(%i) for bin %s\n", indexvbf, AnalysisBins[b].Data());
                        indexcut_ = indexvbf;
                        cutMin=shapeMinVBF;
                        cutMax=shapeMaxVBF;
                    }

                    std::cout<<"DYTEST >>>>>  " << selCh[i]+AnalysisBins[b]+"_"+mainHisto <<"\n";
                    fprintf(gjFile,"[%-10.10s - %s]\n", (selCh[i]+AnalysisBins[b]).Data(), mainHisto.Data());

                    TH2* gjets2Dshape  = (TH2*)pdir->Get(selCh[i]+AnalysisBins[b]+"_"+mainHisto);
                    if(!gjets2Dshape)printf("Can't find histo: %s in g+jets template\n",(selCh[i]+AnalysisBins[b]+"_"+mainHisto).Data());

                    TH1* gjets1Dshape  = gjets2Dshape->ProjectionY("tmpName",indexcut_,indexcut_);

                    //apply the cuts
                    for(int x=0; x<=gjets1Dshape->GetXaxis()->GetNbins()+1; x++) {
                        if(gjets1Dshape->GetXaxis()->GetBinCenter(x)<=cutMin || gjets1Dshape->GetXaxis()->GetBinCenter(x)>=cutMax) {
                            gjets1Dshape->SetBinContent(x,0);
                            gjets1Dshape->SetBinError(x,0);
                        }
                    }
                    //gjets1Dshape->Rebin(2);

                    shapeChan_SI.bckg[ibckg]->SetTitle(DYProcName + " (data)");
                    //for(int i=0; i<shapeChan_SI.bckg[ibckg]->GetNbinsX(); i++) { //THIS IS A BUG!
                    for(int i=1; i<=shapeChan_SI.bckg[ibckg]->GetNbinsX(); i++) {
                        double val = gjets1Dshape->GetBinContent(i);
                        cout << "i: " << i << " \tval: " << val << endl;
                        fprintf(gjFile,"bin: %d, \tval: %f\n", i, val);
                        double err = gjets1Dshape->GetBinError(i);
                        //if(val<0) val=0;
                        //if(AnalysisBins[b]=="eq1jets" && val<0) val=0;
                        double newval = rescaleFactor*val;
                        double newerr = rescaleFactor*err;
                        shapeChan_SI.bckg[ibckg]->SetBinContent(i, newval);
                        shapeChan_SI.bckg[ibckg]->SetBinError  (i, newerr);
                    }
                    delete gjets1Dshape;

                    shapeChan_SI.bckg[ibckg]->Scale(DDRescale);

                    val  = shapeChan_SI.bckg[ibckg]->IntegralAndError(1,shapeChan_SI.bckg[ibckg]->GetXaxis()->GetNbins(),valerr);
                    //systError = std::min(GammaJetSyst * val, fabs(systError - val) ); //syst is difference between G+Jet estimate and MC
                    systError = GammaJetSyst * val; //syst completely coming from G+Jet estimate
                    shapeChan_SI.bckg[ibckg]->SetBinError(0,systError);//save syst error in underflow bin that is always empty
                    if(val<1E-6) {
                        val=0.0;
                        valerr=0.0;
                        systError=-1;
                    }
                    Cval+= string(" &") + toLatexRounded(val,valerr,systError) +" & " + MCTruth;

                    //erase Systematic relatated to DY background
                    std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator dyvars = shapeChan_SI.bckgVars.find(DYProcName);
                    if(dyvars!=shapeChan_SI.bckgVars.end()) {
                        shapeChan_SI.bckgVars.erase(dyvars);
                    }

                    //recompute total background
                    shapeChan_SI.totalBckg->Reset();
                    for(size_t i=0; i<shapeChan_SI.bckg.size(); i++) {
                        shapeChan_SI.totalBckg->Add(shapeChan_SI.bckg[i]);
                    }

                    if(pFile) {
                        fprintf(pFile,"%s\\\\\n", Cval.c_str());
                    }
                }
            }


        }
    }

    if(pFile) {
        fprintf(pFile,"\\hline\n");
        fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
        fprintf(pFile,"\n\n\n\n");
        fclose(pFile);
    }

    //all done with gamma+jet file
    inF->Close();
}



void doWZSubtraction(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString sideBandHisto)
{
    string Ccol   = "\\begin{tabular}{|l|c|c|c|c|";
    string Cname  = "channel & $\\alpha$ measured & $\\alpha$ used & yield data & yield mc";
    string Cval   = "";
    FILE* pFile = NULL;
    if(!fast) {
        pFile = fopen("WZ.tex","w");
        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Non resonant background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
        fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
    }

    for(size_t i=0; i<selCh.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Cval   = selCh[i]+string(" - ")+AnalysisBins[b];

            Shape_t& shapeCtrl_SB = allShapes.find(ctrlCh+AnalysisBins[b]+sideBandHisto)->second;
            Shape_t& shapeCtrl_SI = allShapes.find(ctrlCh+AnalysisBins[b]+mainHisto)->second;
            Shape_t& shapeChan_SB = allShapes.find(selCh[i]+AnalysisBins[b]+sideBandHisto)->second;
            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;

            fprintf(pFile,"#############%s:\n",(string(" &")+selCh[i]+string(" - ")+AnalysisBins[b]).Data());
            fprintf(pFile,"MC: em 3leptons=%6.2E  em 2leptons=%6.2E  ll 3leptons=%6.2E  ll 2leptons=%6.2E\n",shapeCtrl_SB.totalBckg->Integral(), shapeCtrl_SI.totalBckg->Integral(), shapeChan_SB.totalBckg->Integral(), shapeChan_SI.totalBckg->Integral());

            TH1* histo1=NULL, *histo2=NULL, *histo3=NULL, *histo4=NULL;
            for(size_t ibckg=0; ibckg<shapeCtrl_SB.bckg.size(); ibckg++) {
                if(TString(shapeCtrl_SB.bckg[ibckg]->GetTitle()).Contains("WZ"))histo1=shapeCtrl_SB.bckg[ibckg];
            }
            for(size_t ibckg=0; ibckg<shapeCtrl_SI.bckg.size(); ibckg++) {
                if(TString(shapeCtrl_SI.bckg[ibckg]->GetTitle()).Contains("WZ"))histo2=shapeCtrl_SI.bckg[ibckg];
            }
            for(size_t ibckg=0; ibckg<shapeChan_SB.bckg.size(); ibckg++) {
                if(TString(shapeChan_SB.bckg[ibckg]->GetTitle()).Contains("WZ"))histo3=shapeChan_SB.bckg[ibckg];
            }
            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                if(TString(shapeChan_SI.bckg[ibckg]->GetTitle()).Contains("WZ"))histo4=shapeChan_SI.bckg[ibckg];
            }
            fprintf(pFile,"WZ: em 3leptons=%6.2E  em 2leptons=%6.2E  ll 3leptons=%6.2E  ll 2leptons=%6.2E\n",histo1->Integral(), histo2->Integral(), histo3->Integral(), histo4->Integral());

            double Num, Denom, NumError, DenomError;
            Num = histo4->IntegralAndError(1,histo4->GetXaxis()->GetNbins(),NumError);
            Denom = histo3->IntegralAndError(1,histo3->GetXaxis()->GetNbins(),DenomError);
            double ratio = Num/Denom;
            double ratio_err  = sqrt(pow(Num*DenomError,2) + pow(Denom*NumError,2))/ pow(Denom,2);
            double ratio_syst = fabs(histo3->Integral() - shapeChan_SB.totalBckg->Integral())/shapeChan_SB.totalBckg->Integral();
            fprintf(pFile,"Ratio = %s\n",toLatexRounded(ratio,ratio_err,ratio_syst).c_str());

            Double_t valerr;
            Double_t val = histo4->IntegralAndError(1,histo4->GetXaxis()->GetNbins(),valerr);

            Double_t valerr2, valsyst2;
            Double_t val2 = shapeChan_SB.totalBckg->IntegralAndError(1,shapeChan_SB.data->GetXaxis()->GetNbins(),valerr2);
            valerr2= sqrt(pow(valerr2*ratio,2) + pow(val2*ratio_err,2) );
            valsyst2 = val2*ratio_syst;
            val2=val2*ratio;
            fprintf(pFile,"WZ (MC closure test): %s --> %s\n",toLatexRounded(val,valerr).c_str(), toLatexRounded(val2,valerr2,valsyst2).c_str());

            bool noDataObserved=false;
            Double_t valerr3, valsyst3;
            Double_t val3 = shapeChan_SB.data->IntegralAndError(1,shapeChan_SB.data->GetXaxis()->GetNbins(),valerr3);
            if(val3<=0) {
                noDataObserved=true;
                val3=1.0;
                valerr3=1.0;
            }
            valerr3= sqrt(pow(valerr3*ratio,2) + pow(val3*ratio_err,2) );
            valsyst3 = val3*ratio_syst;
            val3=val3*ratio;
            if(!noDataObserved) {
                fprintf(pFile,"WZ (from data)      : %s --> %s\n",toLatexRounded(val,valerr).c_str(), toLatexRounded(val3,valerr3,valsyst3).c_str());
            } else {
                fprintf(pFile,"WZ (from data)      : %s --> smaller than %s (because no data was observed in 3dlepton SideBand--> assume 1+-1 observed data for rescale)\n",toLatexRounded(val,valerr).c_str(), toLatexRounded(val3,valerr3,valsyst3).c_str());
            }


        }
    }

    if(pFile) {
        fclose(pFile);
    }
}

void BlindData(std::vector<TString>& selCh, map<TString, Shape_t>& allShapes, TString mainHisto, bool addSignal)
{
    for(size_t i=0; i<selCh.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
            shapeChan_SI.data->Reset();
            shapeChan_SI.data->Add(shapeChan_SI.totalBckg,1);
            if(addSignal) {
                for(unsigned int s=0; s<shapeChan_SI.signal.size(); s++) {
                    shapeChan_SI.data->Add(shapeChan_SI.signal[s], 1);
                }
            }
        }
    }
}






void SignalInterpolation(std::vector<TString>& selCh,map<TString, Shape_t>& allShapesL, map<TString, Shape_t>& allShapes, map<TString, Shape_t>& allShapesR, TString mainHisto)
{
    if(massL<0 || massR<0 || massL==massR)return;

    std::vector<TString> signalTypes;
    signalTypes.push_back("ggH");
    signalTypes.push_back("qqH");
    for(unsigned int t=0; t<signalTypes.size(); t++) {
        printf("MASS: %i will be interpolated from %i and %i for %s\n",mass, massL,massR,signalTypes[t].Data());

        TString nameL = signalTypes[t]+"(";
        nameL+=massL;
        TString nameR = signalTypes[t]+"(";
        nameR+=massR;

        for(size_t i=0; i<selCh.size(); i++) {
            for(size_t b=0; b<AnalysisBins.size(); b++) {
                Shape_t& shapeL = allShapesL.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
                Shape_t& shape  = allShapes .find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
                Shape_t& shapeR = allShapesR.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;


                int indexL=0, indexR=0;
                for(unsigned int i=0; i<shapeL.signal.size(); i++) {
                    if(TString(shapeL.signal[i]->GetTitle()).BeginsWith(nameL))indexL=i;
                }
                for(unsigned int i=0; i<shapeR.signal.size(); i++) {
                    if(TString(shapeR.signal[i]->GetTitle()).BeginsWith(nameR))indexR=i;
                }

                double Ratio = (mass - massL);
                Ratio/=(massR - massL);

                //centralValue
                TH1* histoL = (TH1*) shapeL.signal[indexL]->Clone("tmpLeft");
                histoL->Scale(1.0/shapeL.xsections[histoL->GetTitle()]);
                TH1* histoR = (TH1*) shapeR.signal[indexR]->Clone("tmpRight");
                histoR->Scale(1.0/shapeR.xsections[histoR->GetTitle()]);
                TString newName = signalTypes[t]+"(";
                newName+= mass;
                TH1* histoNew = (TH1*)histoL->Clone(TString(histoL->GetTitle()).ReplaceAll(nameL,newName));
                histoNew->SetTitle(histoNew->GetName());
                histoNew->Scale(1-Ratio);
                histoNew->Add(histoR,Ratio);

                setTGraph(histoNew->GetTitle(), systpostfix );
                double XSection = TG_xsec->Eval(mass,NULL,"S");
                histoNew->Scale(XSection);
                histoNew->SetBinError(0,0.15*histoNew->Integral());
                shape.signal.insert(shape.signal.begin()+indexL,histoNew);


                //systematics
                std::vector<std::pair<TString, TH1*> > varsLeft  = shapeL.signalVars[histoL->GetTitle()];
                std::vector<std::pair<TString, TH1*> > varsRight = shapeR.signalVars[histoR->GetTitle()];
                std::vector<std::pair<TString, TH1*> > varsNew;
                for(size_t v=0; v<varsLeft.size(); v++) {
                    TH1* histoSystL = (TH1*) varsLeft [v].second->Clone("tmpSystLeft");
                    histoSystL->Scale(1.0/shapeL.xsections[histoL->GetTitle()]);
                    TH1* histoSystR = (TH1*) varsRight[v].second->Clone("tmpSystRight");
                    histoSystR->Scale(1.0/shapeR.xsections[histoR->GetTitle()]);
                    TH1* histoSystNew   = (TH1*)histoSystL->Clone(TString(histoSystL->GetTitle()).ReplaceAll(nameL,newName));
                    histoSystNew->SetTitle(histoSystNew->GetName());
                    histoSystNew->Scale(1-Ratio);
                    histoSystNew->Add(histoSystR,Ratio);
                    histoSystNew->Scale(XSection);
                    //varsNew.push_back(std::make_pair<TString, TH1*>(varsLeft[v].first,histoSystNew));
                    varsNew.push_back(std::pair<TString, TH1*>(varsLeft[v].first,histoSystNew));
                    delete histoSystL;
                    delete histoSystR;
                }
                shape.signalVars[histoNew->GetTitle()] = varsNew;

                delete histoL;
                delete histoR;
            }
        }

    }
}


void RescaleForInterference(std::vector<TString>& selCh,map<TString, Shape_t>& allShapes, TString mainHisto)
{
    TString massStr("");
    if(mass>0)massStr += mass;
    for(size_t i=0; i<selCh.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Shape_t& shape  = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
            //signals
            size_t nsignal=shape.signal.size();
            for(size_t isignal=0; isignal<nsignal; isignal++) {
                TString proc(((TH1D*)shape.signal[isignal])->GetTitle());
                if(mass>0 && !proc.Contains(massStr))continue;
                if(!proc.Contains("ggH"))continue;
                if(mass<400)continue;
                double scaleFactor = 1.45 - 0.00218 * mass + 0.000002625 * mass * mass;
                double scaleFactorDown = 1.0;
                double scaleFactorUp = 1 + (scaleFactor-1)*2;
                printf("Scale Factor : %f [%f,%f] applied on %s\n",scaleFactor, scaleFactorDown, scaleFactorUp, proc.Data());

                ((TH1D*)shape.signal[isignal])->Scale(scaleFactor);
                std::vector<std::pair<TString, TH1*> >& vars = shape.signalVars[proc];
                for(size_t v=0; v<vars.size(); v++) {
                    ((TH1D*)vars[v].second)->Scale(scaleFactor);
                }

                //((TH1D*)shape.signal[isignal])->SetBinContent(0,scaleFactor);

                TH1* down = (TH1D*)shape.signal[isignal]->Clone(proc+"interf_ggHDown");
                down->Scale(scaleFactorDown);
                TH1* up   = (TH1D*)shape.signal[isignal]->Clone(proc+"interf_ggHUp"  );
                up  ->Scale(scaleFactorUp  );
                vars.push_back(std::make_pair("_interf_ggHDown", down) );
                vars.push_back(std::make_pair("_interf_ggHUp"  , up  ) );
            }
        }
    }
}


//TGraph *ggH7TG_xsec=NULL, *ggH7TG_errp=NULL, *ggH7TG_errm=NULL, *ggH7TG_scap=NULL, *ggH7TG_scam=NULL, *ggH7TG_pdfp=NULL, *ggH7TG_pdfm=NULL;
//TGraph *qqH7TG_xsec=NULL, *qqH7TG_errp=NULL, *qqH7TG_errm=NULL, *qqH7TG_scap=NULL, *qqH7TG_scam=NULL, *qqH7TG_pdfp=NULL, *qqH7TG_pdfm=NULL;
//TGraph *ggH8TG_xsec=NULL, *ggH8TG_errp=NULL, *ggH8TG_errm=NULL, *ggH8TG_scap=NULL, *ggH8TG_scam=NULL, *ggH8TG_pdfp=NULL, *ggH8TG_pdfm=NULL;
//TGraph *qqH8TG_xsec=NULL, *qqH8TG_errp=NULL, *qqH8TG_errm=NULL, *qqH8TG_scap=NULL, *qqH8TG_scam=NULL, *qqH8TG_pdfp=NULL, *qqH8TG_pdfm=NULL;
//TGraph *    TG_xsec=NULL, *    TG_errp=NULL, *    TG_errm=NULL, *    TG_scap=NULL, *    TG_scam=NULL, *    TG_pdfp=NULL, *    TG_pdfm=NULL;
void setTGraph(TString proc, TString suffix)
{
    if( suffix.Contains('8') &&  proc.Contains("qq")) {
        TG_xsec=qqH8TG_xsec;
        TG_errp=qqH8TG_errp;
        TG_errm=qqH8TG_errm;
        TG_scap=qqH8TG_scap;
        TG_scam=qqH8TG_scam;
        TG_pdfp=qqH8TG_pdfp;
        TG_pdfm=qqH8TG_pdfm;
    } else if( suffix.Contains('8') && !proc.Contains("qq")) {
        TG_xsec=ggH8TG_xsec;
        TG_errp=ggH8TG_errp;
        TG_errm=ggH8TG_errm;
        TG_scap=ggH8TG_scap;
        TG_scam=ggH8TG_scam;
        TG_pdfp=ggH8TG_pdfp;
        TG_pdfm=ggH8TG_pdfm;
    } else if(!suffix.Contains('8') &&  proc.Contains("qq")) {
        TG_xsec=qqH7TG_xsec;
        TG_errp=qqH7TG_errp;
        TG_errm=qqH7TG_errm;
        TG_scap=qqH7TG_scap;
        TG_scam=qqH7TG_scam;
        TG_pdfp=qqH7TG_pdfp;
        TG_pdfm=qqH7TG_pdfm;
    } else if(!suffix.Contains('8') && !proc.Contains("qq")) {
        TG_xsec=ggH7TG_xsec;
        TG_errp=ggH7TG_errp;
        TG_errm=ggH7TG_errm;
        TG_scap=ggH7TG_scap;
        TG_scam=ggH7TG_scam;
        TG_pdfp=ggH7TG_pdfp;
        TG_pdfm=ggH7TG_pdfm;
    }
}
void initializeTGraph()
{
    double ggH7_mass [] = {90.0,95.0,100.0,105.0,110.0,110.5,111.0,111.5,112.0,112.5,113.0,113.5,114.0,114.5,115.0,115.5,116.0,116.5,117.0,117.5,118.0,118.5,119.0,119.5,120.0,120.5,121.0,121.5,122.0,122.5,123.0,123.5,124.0,124.5,125.0,125.5,126.0,126.5,127.0,127.5,128.0,128.5,129.0,129.5,130.0,130.5,131.0,131.5,132.0,132.5,133.0,133.5,134.0,134.5,135.0,135.5,136.0,136.5,137.0,137.5,138.0,138.5,139.0,139.5,140.0,141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0,151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,160.0,162.0,164.0,166.0,168.0,170.0,172.0,174.0,176.0,178.0,180.0,182.0,184.0,186.0,188.0,190.0,192.0,194.0,196.0,198.0,200.0,202.0,204.0,206.0,208.0,210.0,212.0,214.0,216.0,218.0,220.0,222.0,224.0,226.0,228.0,230.0,232.0,234.0,236.0,238.0,240.0,242.0,244.0,246.0,248.0,250.0,252.0,254.0,256.0,258.0,260.0,262.0,264.0,266.0,268.0,270.0,272.0,274.0,276.0,278.0,280.0,282.0,284.0,286.0,288.0,290.0,295.0,300.0,305.0,310.0,315.0,320.0,325.0,330.0,335.0,340.0,345.0,350.0,360.0,370.0,380.0,390.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0};
    double ggH7_xsec [] = {29.51,26.51,24.00,21.77,19.84,19.66,19.48,19.31,19.13,18.96,18.79,18.63,18.46,18.30,18.14,17.98,17.83,17.67,17.52,17.37,17.22,17.08,16.93,16.79,16.65,16.51,16.37,16.23,16.10,15.97,15.84,15.71,15.58,15.45,15.32,15.20,15.08,14.96,14.85,14.73,14.62,14.50,14.38,14.27,14.16,14.05,13.94,13.83,13.72,13.62,13.51,13.41,13.31,13.21,13.11,13.01,12.91,12.81,12.72,12.62,12.53,12.44,12.35,12.26,12.18,12.00,11.82,11.65,11.49,11.33,11.18,11.02,10.87,10.72,10.58,10.43,10.29,10.16,10.02,9.886,9.754,9.624,9.487,9.349,9.202,8.830,8.519,8.246,8.009,7.786,7.578,7.389,7.212,7.041,6.869,6.696,6.522,6.349,6.179,6.017,5.865,5.725,5.598,5.483,5.377,5.277,5.188,5.106,5.009,4.922,4.833,4.758,4.695,4.608,4.528,4.449,4.381,4.321,4.245,4.177,4.114,4.056,3.990,3.924,3.854,3.789,3.726,3.667,3.611,3.555,3.501,3.449,3.398,3.349,3.301,3.255,3.211,3.167,3.125,3.083,3.044,3.006,2.970,2.934,2.900,2.866,2.833,2.803,2.773,2.744,2.677,2.616,2.563,2.516,2.478,2.443,2.418,2.403,2.398,2.407,2.431,2.428,2.408,2.362,2.283,2.175,2.049,1.776,1.507,1.263,1.050,0.8708,0.7211,0.5976,0.4960,0.4126,0.3444,0.2883,0.2422,0.2042,0.1728,0.1468,0.1252,0.1071,0.09190,0.07930,0.06850,0.05950,0.05180,0.04520,0.03970,0.03480,0.03070,0.02710,0.02410,0.02140,0.01900};
    double ggH7_errp [] = {+16.0,+15.8,+15.5,+15.4,+15.2,+15.2,+15.2,+15.2,+15.2,+15.1,+15.1,+15.1,+15.1,+15.1,+15.1,+15.1,+15.1,+15.0,+15.0,+15.0,+14.9,+14.9,+14.9,+14.8,+14.8,+14.8,+14.8,+14.8,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.6,+14.6,+14.6,+14.6,+14.6,+14.6,+14.6,+14.6,+14.6,+14.6,+14.5,+14.5,+14.5,+14.5,+14.5,+14.5,+14.5,+14.5,+14.5,+14.5,+14.4,+14.4,+14.4,+14.4,+14.4,+14.4,+14.4,+14.3,+14.3,+14.3,+14.3,+14.3,+14.2,+14.2,+14.2,+14.2,+14.1,+14.1,+14.1,+14.0,+14.0,+14.0,+13.9,+13.9,+13.9,+13.9,+13.9,+13.9,+13.9,+13.8,+13.8,+13.7,+13.7,+13.7,+13.7,+13.7,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.5,+13.5,+13.6,+13.7,+13.8,+13.9,+14.1,+14.1,+14.0,+13.9,+13.8,+13.6,+13.6,+13.5,+13.5,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.7,+13.7,+13.7,+13.7,+13.7,+13.7,+13.7,+13.7,+13.7,+13.7,+13.8,+13.8,+13.8,+13.8,+13.8,+13.8,+13.9,+13.9,+14.0,+14.0,+14.1,+14.1,+14.1,+14.1,+14.2,+14.2,+14.3,+14.3,+14.5,+14.7,+15.0,+15.1,+15.2,+15.3,+15.5,+15.6,+15.7,+15.8,+16.0,+16.2,+16.4,+16.5,+16.7,+16.8,+17.0,+17.1,+17.3,+17.4,+17.5,+17.7,+17.9,+18.2,+18.5,+18.8,+19.3,+19.7,+20.1,+20.5,+20.9,+21.2};
    double ggH7_errm [] = {-15.4,-15.3,-15.2,-15.2,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.8,-14.8,-14.8,-14.8,-14.8,-14.8,-14.8,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.8,-14.8,-14.8,-14.8,-14.8,-14.8,-14.8,-14.8,-14.9,-14.9,-14.9,-14.9,-14.9,-14.8,-14.8,-14.8,-14.8,-14.8,-14.7,-14.7,-14.7,-14.7,-14.7,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.5,-14.5,-14.5,-14.5,-14.5,-14.5,-14.5,-14.5,-14.5,-14.5,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.3,-14.3,-14.3,-14.3,-14.3,-14.3,-14.3,-14.3,-14.3,-14.4,-14.4,-14.4,-14.4,-14.3,-14.3,-14.3,-14.3,-14.4,-14.4,-14.4,-14.3,-14.3,-14.3,-14.3,-14.4,-14.3,-14.2,-14.1,-14.0,-13.9,-13.9,-14.0,-14.0,-14.1,-14.1,-14.2,-14.3,-14.4,-14.6,-14.7,-14.9,-15.0,-15.1,-15.2,-15.3,-15.4,-15.6,-15.7,-15.8,-16.0,-16.3,-16.7,-17.0,-17.4,-17.8,-18.2,-18.6,-18.9,-19.2};
    double ggH7_scap [] = {+8.2,+8.0,+7.8,+7.7,+7.5,+7.5,+7.5,+7.5,+7.5,+7.4,+7.4,+7.4,+7.4,+7.4,+7.4,+7.4,+7.4,+7.4,+7.3,+7.3,+7.3,+7.3,+7.3,+7.2,+7.2,+7.2,+7.2,+7.2,+7.2,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.0,+7.0,+7.0,+7.0,+7.0,+7.0,+7.0,+7.0,+7.0,+7.0,+6.9,+6.9,+6.9,+6.9,+6.9,+6.9,+6.9,+6.9,+6.9,+6.9,+6.8,+6.8,+6.8,+6.8,+6.8,+6.8,+6.8,+6.7,+6.7,+6.7,+6.7,+6.7,+6.6,+6.6,+6.6,+6.6,+6.6,+6.5,+6.5,+6.5,+6.5,+6.5,+6.4,+6.4,+6.4,+6.4,+6.4,+6.4,+6.4,+6.3,+6.3,+6.2,+6.2,+6.2,+6.2,+6.2,+6.1,+6.1,+6.1,+6.1,+6.1,+6.1,+6.1,+6.1,+6.0,+6.0,+6.0,+6.0,+6.0,+6.0,+6.1,+6.1,+6.2,+6.4,+6.5,+6.5,+6.4,+6.3,+6.1,+5.9,+5.9,+5.8,+5.8,+5.9,+5.9,+5.9,+5.9,+5.9,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.9,+5.9,+5.9,+5.9,+5.9,+5.9,+5.9,+6.0,+6.0,+6.0,+6.0,+6.0,+6.1,+6.1,+6.2,+6.2,+6.3,+6.3,+6.3,+6.4,+6.4,+6.5,+6.5,+6.5,+6.5,+6.5,+6.6,+6.7,+6.8,+6.8,+6.8,+6.9,+7.0};
    double ggH7_scam [] = {-8.7,-8.6,-8.4,-8.3,-8.1,-8.1,-8.1,-8.1,-8.1,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.7,-7.7,-7.7,-7.7,-7.7,-7.7,-7.7,-7.7,-7.7,-7.7,-7.6,-7.6,-7.6,-7.6,-7.6,-7.6,-7.6,-7.6,-7.6,-7.6,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.4,-7.4,-7.4,-7.4,-7.3,-7.3,-7.3,-7.3,-7.3,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.1,-7.1,-7.0,-7.0,-7.0,-7.0,-7.0,-6.9,-6.9,-6.9,-6.9,-6.9,-6.8,-6.8,-6.8,-6.8,-6.8,-6.8,-6.8,-6.7,-6.7,-6.7,-6.7,-6.6,-6.6,-6.6,-6.6,-6.6,-6.5,-6.5,-6.5,-6.5,-6.5,-6.4,-6.4,-6.4,-6.4,-6.4,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.2,-6.2,-6.2,-6.2,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.0,-6.0,-6.0,-6.0,-6.0,-6.0,-6.0,-5.9,-5.9,-5.9,-5.9,-5.9,-5.8,-5.6,-5.5,-5.4,-5.3,-5.3,-5.3,-5.2,-5.2,-5.2,-5.2,-5.2,-5.2,-5.2,-5.2,-5.2,-5.2,-5.3,-5.3,-5.3,-5.4,-5.4,-5.4,-5.4,-5.4,-5.5,-5.5,-5.6,-5.6,-5.6,-5.7,-5.7,-5.7,-5.7};
    double ggH7_pdfp [] = {+7.8,+7.8,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.5,+7.5,+7.5,+7.5,+7.6,+7.6,+7.6,+7.6,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.8,+7.8,+7.8,+7.8,+7.8,+7.8,+7.8,+7.8,+7.8,+7.8,+7.9,+7.9,+7.9,+7.9,+7.9,+7.9,+7.9,+7.9,+7.9,+7.9,+8.0,+8.0,+8.0,+8.0,+8.0,+8.0,+8.1,+8.1,+8.2,+8.2,+8.3,+8.3,+8.3,+8.3,+8.4,+8.4,+8.4,+8.4,+8.6,+8.8,+9.1,+9.2,+9.3,+9.4,+9.5,+9.6,+9.7,+9.8,+9.9,+10.1,+10.2,+10.4,+10.5,+10.6,+10.7,+10.8,+10.9,+11.0,+11.1,+11.2,+11.4,+11.7,+11.9,+12.3,+12.6,+13.0,+13.3,+13.7,+14.0,+14.2};
    double ggH7_pdfm [] = {-6.7,-6.7,-6.8,-6.9,-6.9,-6.9,-6.9,-6.9,-6.9,-6.9,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.1,-7.1,-7.1,-7.1,-7.1,-7.1,-7.1,-7.1,-7.1,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.4,-7.4,-7.4,-7.4,-7.4,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.6,-7.6,-7.6,-7.6,-7.7,-7.7,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.2,-8.2,-8.2,-8.2,-8.2,-8.3,-8.3,-8.3,-8.3,-8.3,-8.3,-8.3,-8.4,-8.4,-8.4,-8.4,-8.4,-8.4,-8.4,-8.4,-8.5,-8.6,-8.6,-8.6,-8.6,-8.6,-8.7,-8.7,-8.8,-8.9,-9.0,-9.0,-9.1,-9.2,-9.4,-9.5,-9.7,-9.8,-9.8,-9.9,-10.0,-10.1,-10.2,-10.3,-10.4,-10.6,-10.9,-11.1,-11.5,-11.8,-12.2,-12.5,-12.9,-13.2,-13.5};


    /*
       double ggH7_mass [] = {90.0, 95.0, 100.0, 105.0, 110.0, 110.5, 111.0, 111.5, 112.0, 112.5, 113.0, 113.5, 114.0, 114.5, 115.0, 115.5, 116.0, 116.5, 117.0, 117.5, 118.0, 118.5, 119.0, 119.5, 120.0, 120.5, 121.0, 121.5, 122.0, 122.5, 123.0, 123.5, 124.0, 124.5, 125.0, 125.5, 126.0, 126.5, 127.0, 127.5, 128.0, 128.5, 129.0, 129.5, 130.0, 130.5, 131.0, 131.5, 132.0, 132.5, 133.0, 133.5, 134.0, 134.5, 135.0, 135.5, 136.0, 136.5, 137.0, 137.5, 138.0, 138.5, 139.0, 139.5, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 162.0, 164.0, 166.0, 168.0, 170.0, 172.0, 174.0, 176.0, 178.0, 180.0, 182.0, 184.0, 186.0, 188.0, 190.0, 192.0, 194.0, 196.0, 198.0, 200.0, 202.0, 204.0, 206.0, 208.0, 210.0, 212.0, 214.0, 216.0, 218.0, 220.0, 222.0, 224.0, 226.0, 228.0, 230.0, 232.0, 234.0, 236.0, 238.0, 240.0, 242.0, 244.0, 246.0, 248.0, 250.0, 252.0, 254.0, 256.0, 258.0, 260.0, 262.0, 264.0, 266.0, 268.0, 270.0, 272.0, 274.0, 276.0, 278.0, 280.0, 282.0, 284.0, 286.0, 288.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0, 720.0, 740.0, 760.0, 780.0, 800.0, 820.0, 840.0, 860.0, 880.0, 900.0, 920.0, 940.0, 960.0, 980.0, 1000.0};
       double ggH7_xsec [] = {29.47, 26.58, 24.02, 21.78, 19.84, 19.65, 19.48, 19.30, 19.13, 18.95, 18.79, 18.62, 18.45, 18.29, 18.13, 17.97, 17.82, 17.66, 17.51, 17.36, 17.21, 17.06, 16.92, 16.78, 16.63, 16.49, 16.36, 16.22, 16.08, 15.94, 15.82, 15.69, 15.56, 15.43, 15.31, 15.18, 15.06, 14.94, 14.82, 14.70, 14.58, 14.46, 14.35, 14.23, 14.12, 14.01, 13.90, 13.80, 13.69, 13.58, 13.48, 13.37, 13.28, 13.18, 13.08, 12.98, 12.88, 12.78, 12.68, 12.60, 12.50, 12.40, 12.31, 12.22, 12.13, 11.95, 11.78, 11.60, 11.44, 11.27, 11.12, 10.96, 10.80, 10.65, 10.50, 10.36, 10.21, 10.07, 9.934, 9.795, 9.655, 9.514, 9.367, 9.225, 9.080, 8.770, 8.465, 8.190, 7.952, 7.729, 7.515, 7.310, 7.112, 6.923, 6.739, 6.553, 6.379, 6.212, 6.050, 5.896, 5.751, 5.616, 5.487, 5.364, 5.249, 5.136, 5.027, 4.924, 4.822, 4.723, 4.630, 4.539, 4.454, 4.369, 4.288, 4.207, 4.128, 4.053, 3.980, 3.908, 3.839, 3.771, 3.707, 3.643, 3.581, 3.523, 3.468, 3.414, 3.362, 3.312, 3.261, 3.212, 3.164, 3.118, 3.072, 3.028, 2.984, 2.944, 2.903, 2.864, 2.828, 2.793, 2.760, 2.728, 2.696, 2.664, 2.633, 2.603, 2.574, 2.546, 2.480, 2.422, 2.369, 2.322, 2.281, 2.247, 2.221, 2.204, 2.195, 2.198, 2.225, 2.306, 2.361, 2.341, 2.266, 2.158, 2.032, 1.756, 1.482, 1.237, 1.026, 0.8491, 0.7006, 0.5782, 0.4771, 0.3944, 0.3267, 0.2713, 0.2257, 0.1883, 0.1574, 0.1320, 0.1109, 0.09335, 0.07883, 0.06668, 0.05655, 0.04806, 0.04089, 0.03490, 0.02982, 0.02555, 0.02193, 0.01885, 0.01624, 0.01400, 0.01210};
       double ggH7_errp [] = {22.9, 21.9, 21.2, 20.8, 20.4, 20.4, 20.3, 20.3, 20.3, 20.3, 20.2, 20.2, 20.2, 20.1, 20.0, 20.0, 20.0, 19.9, 19.9, 19.9, 19.8, 19.8, 19.7, 19.7, 19.7, 19.7, 19.7, 19.7, 19.7, 19.6, 19.5, 19.5, 19.5, 19.5, 19.5, 19.5, 19.4, 19.4, 19.4, 19.3, 19.3, 19.3, 19.3, 19.2, 19.2, 19.2, 19.2, 19.1, 19.1, 19.0, 19.0, 19.0, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.7, 18.8, 18.8, 18.8, 18.8, 18.8, 18.8, 18.8, 18.7, 18.7, 18.7, 18.7, 18.7, 18.7, 18.7, 18.7, 18.7, 18.6, 18.6, 18.5, 18.5, 18.6, 18.6, 18.6, 18.6, 18.6, 18.5, 18.2, 18.0, 17.9, 17.9, 17.9, 17.9, 18.0, 18.0, 18.1, 17.9, 17.5, 17.3, 17.3, 17.3, 17.3, 17.2, 17.2, 17.2, 17.2, 17.1, 17.1, 16.9, 16.9, 16.9, 16.9, 16.9, 16.9, 16.9, 16.8, 16.8, 16.7, 16.7, 16.6, 16.6, 16.6, 16.6, 16.7, 16.7, 16.7, 16.7, 16.7, 16.6, 16.6, 16.5, 16.5, 16.4, 16.3, 16.3, 16.2, 16.2, 16.3, 16.2, 16.2, 16.2, 16.2, 16.1, 16.1, 16.1, 16.0, 16.0, 16.0, 16.0, 16.1, 16.1, 16.6, 16.8, 16.8, 16.8, 16.8, 16.9, 17.1, 17.2, 17.4, 18.1, 19.5, 19.7, 19.4, 18.2, 17.1, 16.4, 15.7, 14.8, 15.6, 16.4, 17.1, 17.7, 18.1, 18.4, 18.7, 19.1, 19.4, 19.6, 20.0, 20.3, 20.7, 20.9, 21.1, 21.5, 21.7, 22.0, 22.2, 22.5, 23.1, 23.4, 24.1, 24.6, 25.1, 25.6, 26.0, 26.6, 27.1};
       double ggH7_errm [] = {-15.6, -15.9, -15.6, -15.5, -15.3, -15.3, -15.3, -15.3, -15.3, -15.2, -15.3, -15.3, -15.3, -15.3, -15.3, -15.3, -15.3, -15.3, -15.3, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.0, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.2, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -15.1, -14.9, -15.0, -15.0, -15.0, -15.0, -15.0, -15.0, -15.0, -15.0, -15.0, -15.0, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -14.9, -15.0, -15.0, -15.0, -15.0, -15.0, -15.0, -14.9, -14.8, -14.7, -14.9, -14.9, -14.9, -14.8, -14.7, -14.7, -14.7, -14.6, -14.9, -15.0, -15.1, -15.0, -15.1, -15.1, -15.2, -15.1, -15.2, -15.2, -15.2, -15.3, -15.3, -15.3, -15.3, -15.3, -15.3, -15.3, -15.3, -15.3, -15.4, -15.4, -15.5, -15.5, -15.5, -15.4, -15.5, -15.5, -15.4, -15.4, -15.5, -15.6, -15.6, -15.6, -15.7, -15.7, -15.8, -15.8, -15.9, -15.8, -15.8, -15.9, -15.8, -15.8, -15.9, -15.9, -16.1, -16.1, -16.2, -16.2, -16.2, -16.2, -16.1, -16.1, -15.8, -15.5, -15.6, -15.6, -15.8, -15.8, -15.8, -15.8, -15.7, -15.3, -14.8, -14.8, -14.6, -15.1, -15.6, -16.0, -16.3, -16.9, -17.2, -17.2, -17.5, -17.6, -17.6, -17.7, -17.8, -17.8, -17.9, -18.0, -18.0, -18.3, -18.3, -18.4, -18.6, -18.6, -18.8, -18.9, -19.1, -19.4, -19.6, -19.9, -20.2, -20.6, -21.1, -21.4, -21.9, -22.2, -22.6};
       double ggH7_scap [] = {14.8, 13.9, 13.3, 12.8, 12.5, 12.5, 12.4, 12.4, 12.4, 12.3, 12.3, 12.3, 12.2, 12.1, 12.1, 12.1, 12.1, 12.0, 12.0, 12.0, 12.0, 12.0, 11.9, 11.9, 11.9, 11.9, 11.9, 11.8, 11.8, 11.8, 11.7, 11.7, 11.7, 11.7, 11.7, 11.7, 11.6, 11.6, 11.6, 11.5, 11.5, 11.5, 11.5, 11.3, 11.3, 11.3, 11.3, 11.3, 11.3, 11.2, 11.2, 11.2, 11.1, 11.1, 11.1, 11.1, 11.1, 11.1, 11.1, 10.9, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.8, 10.9, 10.8, 10.8, 10.9, 10.9, 10.9, 10.9, 10.9, 10.8, 10.5, 10.3, 10.2, 10.2, 10.2, 10.2, 10.3, 10.3, 10.3, 10.2, 9.8, 9.6, 9.6, 9.6, 9.6, 9.5, 9.5, 9.4, 9.4, 9.3, 9.3, 9.2, 9.2, 9.2, 9.2, 9.2, 9.1, 9.1, 9.0, 9.0, 8.9, 8.9, 8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 8.9, 8.8, 8.8, 8.7, 8.6, 8.6, 8.5, 8.5, 8.3, 8.3, 8.3, 8.3, 8.3, 8.2, 8.2, 8.2, 8.2, 8.1, 8.1, 8.0, 7.9, 7.9, 7.9, 7.9, 8.0, 8.0, 8.4, 8.6, 8.6, 8.5, 8.5, 8.6, 8.7, 8.8, 9.0, 9.6, 11.0, 11.1, 10.8, 9.5, 8.5, 7.6, 6.8, 5.8, 6.4, 7.0, 7.6, 8.0, 8.3, 8.5, 8.7, 8.9, 9.0, 9.1, 9.4, 9.5, 9.8, 9.9, 10.0, 10.3, 10.4, 10.6, 10.7, 10.8, 10.9, 11.0, 11.3, 11.5, 11.6, 11.9, 11.9, 12.1, 12.3};
       double ggH7_scam [] = {-8.7, -9.0, -8.6, -8.5, -8.2, -8.2, -8.2, -8.2, -8.2, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.0, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.7, -7.7, -7.7, -7.7, -7.7, -7.6, -7.6, -7.7, -7.7, -7.5, -7.5, -7.5, -7.6, -7.6, -7.6, -7.6, -7.6, -7.5, -7.5, -7.6, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.2, -7.2, -7.2, -7.2, -7.3, -7.2, -7.2, -7.3, -7.3, -7.3, -7.2, -7.2, -7.2, -7.1, -6.9, -6.8, -6.9, -7.0, -6.9, -6.9, -6.7, -6.7, -6.7, -6.6, -6.9, -7.0, -7.1, -7.0, -7.1, -7.1, -7.2, -7.1, -7.2, -7.2, -7.2, -7.3, -7.3, -7.2, -7.3, -7.3, -7.2, -7.2, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.2, -7.3, -7.3, -7.3, -7.4, -7.4, -7.4, -7.5, -7.5, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.5, -7.6, -7.6, -7.7, -7.8, -7.9, -7.9, -7.9, -7.8, -7.7, -7.7, -7.3, -7.1, -7.1, -7.2, -7.3, -7.2, -7.3, -7.3, -7.1, -6.7, -6.1, -6.2, -5.9, -6.3, -6.9, -7.3, -7.7, -8.2, -8.4, -8.4, -8.5, -8.6, -8.6, -8.5, -8.5, -8.4, -8.3, -8.3, -8.2, -8.3, -8.2, -8.3, -8.3, -8.2, -8.3, -8.3, -8.4, -8.4, -8.3, -8.4, -8.4, -8.4, -8.5, -8.5, -8.6, -8.6, -8.6};
       double ggH7_pdfp [] = {8.1, 8.0, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.7, 7.8, 7.8, 7.8, 7.8, 7.7, 7.7, 7.7, 7.7, 7.7, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.1, 8.1, 8.1, 8.1, 8.1, 8.1, 8.1, 8.1, 8.2, 8.2, 8.2, 8.2, 8.3, 8.3, 8.3, 8.4, 8.4, 8.5, 8.5, 8.6, 8.6, 8.6, 8.6, 8.8, 8.9, 9.1, 9.3, 9.4, 9.5, 9.7, 9.8, 9.9, 10.0, 10.2, 10.4, 10.5, 10.6, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.8, 12.1, 12.4, 12.7, 13.0, 13.4, 13.8, 14.1, 14.4, 14.8};
       double ggH7_pdfm [] = {-6.9, -6.9, -7.0, -7.1, -7.1, -7.1, -7.1, -7.1, -7.1, -7.1, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.2, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.4, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.6, -7.6, -7.6, -7.6, -7.6, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.8, -7.8, -7.8, -7.8, -7.9, -7.9, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.2, -8.2, -8.2, -8.2, -8.2, -8.2, -8.2, -8.2, -8.2, -8.2, -8.3, -8.3, -8.3, -8.3, -8.3, -8.3, -8.3, -8.3, -8.3, -8.3, -8.3, -8.3, -8.3, -8.3, -8.3, -8.4, -8.4, -8.4, -8.4, -8.4, -8.4, -8.5, -8.5, -8.5, -8.5, -8.6, -8.6, -8.6, -8.7, -8.7, -8.7, -8.8, -8.7, -8.7, -8.7, -8.7, -8.8, -8.8, -8.9, -9.0, -9.1, -9.2, -9.3, -9.4, -9.6, -9.7, -9.8, -10.0, -10.1, -10.2, -10.3, -10.4, -10.5, -10.6, -10.7, -10.9, -11.3, -11.6, -11.9, -12.2, -12.6, -12.9, -13.3, -13.6, -14.0};
    */

    double qqH7_mass [] = {90.0, 95.0, 100.0, 105.0, 110.0, 110.5, 111.0, 111.5, 112.0, 112.5, 113.0, 113.5, 114.0, 114.5, 115.0, 115.5, 116.0, 116.5, 117.0, 117.5, 118.0, 118.5, 119.0, 119.5, 120.0, 120.5, 121.0, 121.5, 122.0, 122.5, 123.0, 123.5, 124.0, 124.5, 125.0, 125.5, 126.0, 126.5, 127.0, 127.5, 128.0, 128.5, 129.0, 129.5, 130.0, 130.5, 131.0, 131.5, 132.0, 132.5, 133.0, 133.5, 134.0, 134.5, 135.0, 135.5, 136.0, 136.5, 137.0, 137.5, 138.0, 138.5, 139.0, 139.5, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 162.0, 164.0, 166.0, 168.0, 170.0, 172.0, 174.0, 176.0, 178.0, 180.0, 182.0, 184.0, 186.0, 188.0, 190.0, 192.0, 194.0, 196.0, 198.0, 200.0, 202.0, 204.0, 206.0, 208.0, 210.0, 212.0, 214.0, 216.0, 218.0, 220.0, 222.0, 224.0, 226.0, 228.0, 230.0, 232.0, 234.0, 236.0, 238.0, 240.0, 242.0, 244.0, 246.0, 248.0, 250.0, 252.0, 254.0, 256.0, 258.0, 260.0, 262.0, 264.0, 266.0, 268.0, 270.0, 272.0, 274.0, 276.0, 278.0, 280.0, 282.0, 284.0, 286.0, 288.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0, 720.0, 740.0, 760.0, 780.0, 800.0, 820.0, 840.0, 860.0, 880.0, 900.0, 920.0, 940.0, 960.0, 980.0, 1000.0};
    double qqH7_xsec [] = {1.710, 1.628, 1.546, 1.472, 1.398, 1.391, 1.384, 1.378, 1.371, 1.364, 1.358, 1.351, 1.345, 1.339, 1.332, 1.326, 1.319, 1.313, 1.307, 1.300, 1.294, 1.288, 1.282, 1.276, 1.269, 1.263, 1.257, 1.251, 1.246, 1.240, 1.234, 1.228, 1.222, 1.216, 1.211, 1.205, 1.199, 1.193, 1.188, 1.182, 1.176, 1.171, 1.165, 1.159, 1.154, 1.148, 1.143, 1.137, 1.132, 1.126, 1.121, 1.115, 1.110, 1.105, 1.100, 1.095, 1.090, 1.085, 1.080, 1.076, 1.071, 1.066, 1.062, 1.057, 1.052, 1.043, 1.033, 1.023, 1.013, 1.004, 0.9951, 0.9866, 0.9782, 0.9699, 0.9617, 0.9529, 0.9441, 0.9353, 0.9266, 0.9180, 0.9095, 0.9013, 0.8934, 0.8859, 0.8787, 0.8676, 0.8571, 0.8453, 0.8316, 0.8173, 0.8029, 0.7885, 0.7744, 0.7609, 0.7480, 0.7361, 0.7248, 0.7139, 0.7032, 0.6925, 0.6812, 0.6699, 0.6587, 0.6478, 0.6371, 0.6267, 0.6164, 0.6064, 0.5965, 0.5869, 0.5775, 0.5684, 0.5594, 0.5506, 0.5420, 0.5335, 0.5252, 0.5170, 0.5089, 0.5011, 0.4934, 0.4859, 0.4785, 0.4712, 0.4641, 0.4572, 0.4503, 0.4436, 0.4369, 0.4304, 0.4239, 0.4174, 0.4111, 0.4049, 0.3988, 0.3931, 0.3875, 0.3821, 0.3767, 0.3715, 0.3663, 0.3611, 0.3560, 0.3510, 0.3461, 0.3413, 0.3365, 0.3318, 0.3271, 0.3226, 0.3116, 0.3011, 0.2908, 0.2809, 0.2716, 0.2627, 0.2539, 0.2453, 0.2368, 0.2286, 0.2206, 0.2132, 0.2018, 0.1910, 0.1808, 0.1712, 0.1620, 0.1451, 0.1304, 0.1171, 0.1054, 0.09497, 0.08568, 0.07746, 0.07010, 0.06353, 0.05771, 0.05246, 0.04776, 0.04356, 0.03977, 0.03637, 0.03330, 0.03052, 0.02805, 0.02580, 0.02373, 0.02188, 0.02018, 0.01864, 0.01724, 0.01597, 0.01479, 0.01375, 0.01275, 0.01186, 0.01104};
    double qqH7_errp [] = {2.7, 2.5, 2.6, 2.5, 2.8, 2.8, 2.8, 2.7, 2.7, 2.7, 2.6, 2.6, 2.6, 2.6, 2.5, 2.6, 2.6, 2.6, 2.6, 2.7, 2.7, 2.7, 2.7, 2.7, 2.8, 2.7, 2.8, 2.8, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.9, 2.9, 2.9, 2.9, 2.9, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 2.9, 2.9, 2.9, 2.9, 2.9, 2.8, 2.8, 2.9, 2.9, 3.0, 3.1, 3.1, 3.1, 3.0, 3.0, 3.0, 2.9, 3.0, 3.0, 3.1, 3.1, 3.1, 3.1, 3.0, 3.0, 2.9, 2.9, 3.0, 3.1, 3.1, 3.1, 3.1, 3.2, 3.2, 3.2, 3.1, 3.1, 3.2, 3.4, 3.4, 3.3, 3.3, 3.3, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.6, 3.6, 3.7, 3.7, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.9, 3.9, 4.0, 4.0, 4.1, 4.2, 4.2, 4.3, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.3, 4.3, 4.3, 4.4, 4.4, 4.4, 4.5, 4.5, 4.5, 4.6, 4.7, 4.8, 4.8, 4.9, 5.0, 5.0, 5.0, 5.1, 5.2, 5.2, 5.3, 5.5, 5.7, 5.8, 5.9, 6.2, 6.5, 6.7, 7.0, 7.2, 7.5, 7.8, 8.0, 8.3, 8.6, 8.9, 9.2, 9.4, 9.7, 9.9, 10.2, 10.5, 10.8, 11.0, 11.3, 11.5, 11.8, 12.1, 12.3, 12.6, 12.9, 13.3, 13.6, 13.9, 14.2};
    double qqH7_errm [] = {-2.3, -2.5, -2.4, -2.4, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.5, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.1, -2.1, -2.1, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.1, -2.1, -2.1, -2.1, -2.2, -2.2, -2.2, -2.3, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.1, -2.2, -2.3, -2.4, -2.3, -2.2, -2.2, -2.2, -2.2, -2.3, -2.4, -2.4, -2.4, -2.3, -2.3, -2.3, -2.3, -2.4, -2.4, -2.4, -2.4, -2.4, -2.4, -2.5, -2.5, -2.5, -2.5, -2.4, -2.4, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.6, -2.6, -2.6, -2.6, -2.6, -2.5, -2.5, -2.4, -2.4, -2.4, -2.5, -2.5, -2.5, -2.6, -2.6, -2.6, -2.6, -2.7, -2.7, -2.7, -2.7, -2.7, -2.6, -2.6, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.8, -2.8, -2.9, -2.9, -2.9, -3.0, -3.0, -3.0, -3.0, -3.0, -3.1, -3.1, -3.2, -3.3, -3.4, -3.4, -3.5, -3.5, -3.7, -3.8, -3.8, -3.8, -3.9, -4.0, -4.0, -4.1, -4.2, -4.2, -4.3, -4.3, -4.4, -4.5, -4.5, -4.6, -4.6, -4.7, -4.7, -4.8, -4.8, -4.9};
    double qqH7_scap [] = {0.6, 0.4, 0.4, 0.3, 0.5, 0.5, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.5, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.0, 0.1, 0.3, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.5, 0.6, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 1.0, 1.0, 1.1, 1.1, 1.2, 1.2, 1.3, 1.4, 1.4, 1.4, 1.5, 1.5, 1.5, 1.6, 1.7, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2};
    double qqH7_scam [] = {-0.2, -0.4, -0.3, -0.3, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.4, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.2, -0.2, -0.1, -0.1, -0.1, -0.1, -0.0, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.0, -0.1, -0.1, -0.1, -0.2, -0.2, -0.2, -0.1, -0.1, -0.1, -0.2, -0.1, -0.1, -0.1, -0.2, -0.3, -0.2, -0.2, -0.1, -0.2, -0.2, -0.3, -0.4, -0.4, -0.3, -0.2, -0.3, -0.3, -0.3, -0.3, -0.3, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.5, -0.5, -0.5, -0.5, -0.5, -0.6, -0.6, -0.6, -0.6, -0.6, -0.5, -0.5, -0.4, -0.4, -0.4, -0.5, -0.5, -0.5, -0.6, -0.6, -0.6, -0.6, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.8, -0.8, -0.8, -0.8, -0.7, -0.8, -0.9, -0.9, -0.9, -1.0, -1.0, -1.1, -1.1, -1.1, -1.1, -1.2, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.7, -1.8, -1.9, -2.0, -2.1, -2.1, -2.2, -2.3, -2.4, -2.5, -2.6, -2.7, -2.7, -2.8, -2.8, -2.9, -3.0, -3.1, -3.2, -3.2, -3.3, -3.3, -3.4, -3.5};
    double qqH7_pdfp [] = {2.1, 2.1, 2.2, 2.2, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.9, 2.9, 2.9, 2.9, 3.0, 3.0, 3.0, 3.0, 3.0, 3.1, 3.1, 3.1, 3.1, 3.2, 3.2, 3.2, 3.2, 3.2, 3.3, 3.3, 3.3, 3.3, 3.3, 3.4, 3.4, 3.4, 3.4, 3.5, 3.5, 3.5, 3.5, 3.5, 3.6, 3.6, 3.6, 3.6, 3.7, 3.7, 3.7, 3.7, 3.7, 3.8, 3.8, 3.8, 3.8, 3.8, 3.9, 3.9, 3.9, 3.9, 4.0, 4.0, 4.0, 4.0, 4.0, 4.1, 4.1, 4.1, 4.1, 4.2, 4.2, 4.2, 4.2, 4.2, 4.3, 4.3, 4.4, 4.4, 4.5, 4.5, 4.6, 4.6, 4.7, 4.8, 4.8, 4.9, 4.9, 5.0, 5.1, 5.2, 5.3, 5.5, 5.7, 5.9, 6.1, 6.3, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.9, 8.1, 8.3, 8.5, 8.7, 8.9, 9.2, 9.4, 9.6, 9.8, 10.0, 10.3, 10.5, 10.7, 10.9, 11.1, 11.4, 11.6, 11.8, 12.0};
    double qqH7_pdfm [] = {-2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -1.9, -1.9, -1.9, -1.9, -1.9, -1.9, -1.9, -1.9, -1.9, -1.9, -1.9, -1.9, -1.9, -1.9, -1.8, -1.8, -1.8, -1.8, -1.8, -1.8, -1.7, -1.7, -1.7, -1.7, -1.7, -1.7, -1.6, -1.6, -1.6, -1.6, -1.6, -1.6, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.4, -1.4, -1.4};


    double ggH8_mass [] = {80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 110.5, 111.0, 111.5, 112.0, 112.5, 113.0, 113.5, 114.0, 114.5, 115.0, 115.5, 116.0, 116.5, 117.0, 117.5, 118.0, 118.5, 119.0, 119.5, 120.0, 120.5, 121.0, 121.5, 122.0, 122.5, 123.0, 123.5, 124.0, 124.5, 125.0, 125.5, 126.0, 126.5, 127.0, 127.5, 128.0, 128.5, 129.0, 129.5, 130.0, 130.5, 131.0, 131.5, 132.0, 132.5, 133.0, 133.5, 134.0, 134.5, 135.0, 135.5, 136.0, 136.5, 137.0, 137.5, 138.0, 138.5, 139.0, 139.5, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 162.0, 164.0, 165.0, 166.0, 168.0, 170.0, 172.0, 174.0, 175.0, 176.0, 178.0, 180.0, 182.0, 184.0, 185.0, 186.0, 188.0, 190.0, 192.0, 194.0, 195.0, 196.0, 198.0, 200.0, 202.0, 204.0, 206.0, 208.0, 210.0, 212.0, 214.0, 216.0, 218.0, 220.0, 222.0, 224.0, 226.0, 228.0, 230.0, 232.0, 234.0, 236.0, 238.0, 240.0, 242.0, 244.0, 246.0, 248.0, 250.0, 252.0, 254.0, 256.0, 258.0, 260.0, 262.0, 264.0, 266.0, 268.0, 270.0, 272.0, 274.0, 276.0, 278.0, 280.0, 282.0, 284.0, 286.0, 288.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, 420.0, 440.0, 450.0, 460.0, 480.0, 500.0, 520.0, 540.0, 550.0, 560.0, 580.0, 600.0, 620.0, 640.0, 650.0, 660.0, 680.0, 700.0, 720.0, 740.0, 750.0, 760.0, 780.0, 800.0, 820.0, 840.0, 850.0, 860.0, 880.0, 900.0, 920.0, 940.0, 950.0, 960.0, 980.0, 1000.0};
    double ggH8_xsec [] = {46.12, 45.04, 43.99, 42.99, 42.01, 41.07, 40.17, 39.29, 38.44, 37.62, 36.80, 36.05, 35.30, 34.58, 33.87, 33.19, 32.53, 31.89, 31.27, 30.66, 30.12, 29.55, 28.99, 28.44, 27.92, 27.39, 26.89, 26.42, 25.95, 25.49, 25.04, 24.82, 24.60, 24.39, 24.18, 24.05, 23.76, 23.56, 23.36, 23.16, 22.96, 22.84, 22.58, 22.39, 22.20, 22.09, 21.90, 21.72, 21.55, 21.37, 21.20, 20.96, 20.86, 20.69, 20.53, 20.37, 20.21, 20.04, 19.89, 19.73, 19.57, 19.42, 19.27, 19.07, 18.97, 18.78, 18.67, 18.53, 18.39, 18.21, 18.11, 17.98, 17.84, 17.68, 17.55, 17.42, 17.29, 17.19, 17.07, 16.92, 16.82, 16.70, 16.60, 16.46, 16.36, 16.24, 16.13, 16.01, 15.89, 15.78, 15.67, 15.43, 15.20, 14.98, 14.79, 14.59, 14.40, 14.31, 14.11, 13.92, 13.65, 13.33, 13.12, 12.97, 12.86, 12.84, 12.92, 12.84, 12.52, 12.17, 11.98, 11.96, 11.27, 10.96, 10.69, 10.40, 10.15, 9.895, 9.637, 9.531, 9.429, 9.203, 8.923, 8.634, 8.410, 8.315, 8.219, 7.994, 7.984, 7.724, 7.443, 7.431, 7.430, 7.237, 7.127, 7.002, 6.885, 6.767, 6.644, 6.534, 6.421, 6.327, 6.249, 6.138, 6.038, 5.938, 5.851, 5.776, 5.680, 5.593, 5.513, 5.438, 5.352, 5.272, 5.183, 5.101, 5.020, 4.945, 4.874, 4.802, 4.734, 4.667, 4.602, 4.539, 4.479, 4.421, 4.366, 4.308, 4.253, 4.198, 4.149, 4.100, 4.053, 4.008, 3.964, 3.921, 3.880, 3.841, 3.803, 3.767, 3.683, 3.606, 3.539, 3.482, 3.434, 3.392, 3.364, 3.349, 3.349, 3.367, 3.405, 3.406, 3.390, 3.336, 3.235, 3.093, 2.924, 2.552, 2.180, 2.003, 1.837, 1.538, 1.283, 1.069, 0.8912, 0.8141, 0.7442, 0.6229, 0.5230, 0.4403, 0.3718, 0.3423, 0.3152, 0.2680, 0.2288, 0.1962, 0.1687, 0.1566, 0.1455, 0.1260, 0.1095, 0.09547, 0.08346, 0.07811, 0.07321, 0.06443, 0.05684, 0.05030, 0.04463, 0.04206, 0.03969, 0.03539, 0.03163};
    double ggH8_errp [] = {+16.7, +16.7, +16.6, +16.6, +16.5, +16.5, +16.4, +16.4, +16.3, +16.2, +16.1, +16.1, +16.0, +16.0, +15.9, +15.9, +15.9, +15.8, +15.8, +15.8, +15.7, +15.7, +15.7, +15.6, +15.6, +15.5, +15.4, +15.4, +15.4, +15.3, +15.3, +15.3, +15.2, +15.1, +15.1, +15.1, +15.1, +15.1, +15.1, +15.1, +15.0, +15.0, +14.9, +14.9, +14.9, +14.9, +14.9, +14.8, +14.8, +14.8, +14.8, +14.8, +14.8, +14.8, +14.8, +14.7, +14.7, +14.7, +14.7, +14.7, +14.7, +14.7, +14.7, +14.7, +14.6, +14.6, +14.6, +14.6, +14.6, +14.6, +14.6, +14.6, +14.6, +14.5, +14.5, +14.5, +14.4, +14.4, +14.4, +14.4, +14.4, +14.4, +14.3, +14.3, +14.3, +14.3, +14.3, +14.3, +14.3, +14.3, +14.3, +14.1, +14.1, +14.1, +14.1, +14.1, +14.0, +14.0, +14.0, +14.1, +14.1, +14.1, +14.0, +14.1, +14.1, +14.1, +14.1, +14.1, +14.0, +14.0, +14.0, +14.0, +13.9, +13.9, +13.9, +13.9, +13.9, +13.8, +13.7, +13.7, +13.7, +13.7, +13.6, +13.6, +13.6, +13.6, +13.5, +13.5, +13.5, +13.5, +13.5, +13.5, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.3, +13.3, +13.2, +13.2, +13.3, +13.3, +13.3, +13.3, +13.3, +13.2, +13.2, +13.2, +13.2, +13.2, +13.2, +13.2, +13.2, +13.2, +13.2, +13.3, +13.3, +13.3, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.3, +13.3, +13.3, +13.3, +13.3, +13.3, +13.3, +13.3, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.5, +13.6, +13.6, +13.6, +13.7, +13.8, +13.9, +13.9, +14.0, +14.0, +14.1, +14.3, +14.4, +14.5, +14.7, +14.9, +15.0, +15.2, +15.2, +15.3, +15.4, +15.4, +15.5, +15.6, +15.6, +15.7, +15.9, +16.1, +16.2, +16.5, +16.5, +16.6, +16.6, +16.7, +16.9, +17.0, +17.1, +17.1, +17.3, +17.4, +17.5, +17.8, +17.9, +18.0, +18.3, +18.6};
    double ggH8_errm [] = {-15.9, -15.8, -15.7, -15.6, -15.6, -15.6, -15.5, -15.5, -15.4, -15.4, -15.4, -15.3, -15.2, -15.2, -15.1, -15.1, -15.1, -15.0, -15.0, -15.0, -14.9, -14.9, -14.9, -14.9, -14.8, -14.8, -14.9, -14.9, -14.8, -14.9, -14.9, -14.9, -14.9, -14.8, -14.8, -14.8, -14.8, -14.9, -14.9, -14.9, -14.9, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.6, -14.6, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.6, -14.6, -14.6, -14.6, -14.6, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.5, -14.5, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.4, -14.4, -14.4, -14.5, -14.5, -14.5, -14.4, -14.5, -14.5, -14.5, -14.5, -14.4, -14.5, -14.5, -14.4, -14.4, -14.3, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -13.9, -13.8, -13.7, -13.7, -13.8, -13.9, -14.0, -14.0, -14.1, -14.1, -14.1, -14.2, -14.2, -14.2, -14.2, -14.1, -14.1, -14.1, -14.1, -14.1, -14.0, -13.9, -14.0, -14.0, -14.0, -14.1, -14.1, -14.0, -14.0, -14.1, -14.1, -14.1, -14.0, -13.8, -13.7, -13.6, -13.6, -13.7, -13.6, -13.6, -13.7, -13.6, -13.7, -13.7, -13.8, -13.8, -13.8, -13.8, -13.9, -14.0, -14.0, -14.2, -14.3, -14.4, -14.6, -14.7, -14.8, -14.9, -15.0, -15.0, -15.0, -15.1, -15.2, -15.3, -15.4, -15.5, -15.7, -16.0, -16.1, -16.2, -16.4, -16.7};
    double ggH8_scap [] = {+8.8, +8.8, +8.7, +8.7, +8.6, +8.6, +8.5, +8.5, +8.4, +8.4, +8.3, +8.3, +8.2, +8.2, +8.1, +8.1, +8.1, +8.0, +8.0, +8.0, +7.9, +7.9, +7.9, +7.8, +7.8, +7.8, +7.7, +7.7, +7.7, +7.6, +7.6, +7.6, +7.6, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.2, +7.2, +7.2, +7.2, +7.2, +7.2, +7.2, +7.2, +7.2, +7.1, +7.1, +7.1, +7.1, +7.1, +7.1, +7.1, +7.1, +7.1, +7.0, +7.0, +7.0, +7.0, +7.0, +7.0, +7.0, +7.0, +7.0, +6.9, +6.9, +6.9, +6.9, +6.9, +6.9, +6.9, +6.9, +6.9, +6.8, +6.8, +6.8, +6.8, +6.8, +6.7, +6.7, +6.7, +6.7, +6.7, +6.7, +6.6, +6.6, +6.6, +6.6, +6.6, +6.6, +6.5, +6.5, +6.5, +6.5, +6.4, +6.4, +6.4, +6.4, +6.4, +6.3, +6.3, +6.3, +6.3, +6.3, +6.2, +6.2, +6.2, +6.2, +6.1, +6.1, +6.1, +6.1, +6.1, +6.1, +6.0, +6.0, +6.0, +6.0, +6.0, +6.0, +6.0, +6.0, +6.0, +6.0, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +6.0, +6.0, +6.0, +6.1, +6.1, +6.1, +6.1, +6.1, +6.2, +6.2, +6.2, +6.2, +6.3, +6.3, +6.3, +6.4, +6.4, +6.4, +6.5, +6.5};
    double ggH8_scam [] = {-9.2, -9.1, -9.1, -9.0, -9.0, -9.0, -8.9, -8.9, -8.8, -8.8, -8.8, -8.7, -8.7, -8.7, -8.6, -8.6, -8.6, -8.5, -8.5, -8.5, -8.4, -8.4, -8.4, -8.4, -8.3, -8.3, -8.3, -8.3, -8.2, -8.2, -8.2, -8.2, -8.2, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.2, -7.2, -7.2, -7.2, -7.2, -7.1, -7.1, -7.1, -7.1, -7.1, -7.0, -7.0, -7.0, -7.0, -7.0, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.8, -6.8, -6.8, -6.8, -6.8, -6.7, -6.7, -6.7, -6.7, -6.7, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.4, -6.4, -6.4, -6.4, -6.4, -6.4, -6.4, -6.3, -6.3, -6.3, -6.3, -6.3, -6.3, -6.3, -6.2, -6.2, -6.2, -6.2, -6.2, -6.2, -6.2, -6.2, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -5.9, -5.9, -5.9, -5.9, -5.9, -5.8, -5.6, -5.5, -5.4, -5.3, -5.3, -5.2, -5.2, -5.2, -5.1, -5.1, -5.1, -5.1, -5.1, -5.1, -5.0, -5.0, -5.0, -5.0, -5.1, -5.1, -5.1, -5.1, -5.1, -5.1, -5.2, -5.2, -5.2, -5.2, -5.2, -5.3, -5.3, -5.3, -5.3, -5.3, -5.4, -5.4, -5.4, -5.4, -5.4};
    double ggH8_pdfp [] = {+7.9, +7.9, +7.9, +7.9, +7.9, +7.9, +7.9, +7.9, +7.9, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.7, +7.7, +7.7, +7.7, +7.7, +7.7, +7.7, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.4, +7.4, +7.4, +7.4, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.3, +7.3, +7.4, +7.4, +7.4, +7.4, +7.4, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.4, +7.4, +7.4, +7.4, +7.5, +7.5, +7.5, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.7, +7.7, +7.7, +7.7, +7.7, +7.7, +7.8, +7.9, +7.9, +7.9, +8.0, +8.0, +8.1, +8.1, +8.2, +8.2, +8.3, +8.5, +8.6, +8.7, +8.9, +9.1, +9.2, +9.4, +9.4, +9.4, +9.5, +9.5, +9.6, +9.7, +9.7, +9.8, +9.9, +10.1, +10.2, +10.4, +10.4, +10.5, +10.5, +10.6, +10.7, +10.8, +10.9, +10.9, +11.0, +11.1, +11.2, +11.4, +11.5, +11.6, +11.8, +12.1};
    double ggH8_pdfm [] = {-6.7, -6.7, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.6, -6.6, -6.6, -6.7, -6.7, -6.7, -6.7, -6.7, -6.7, -6.7, -6.7, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -7.0, -7.0, -7.0, -7.0, -7.1, -7.1, -7.1, -7.1, -7.1, -7.1, -7.1, -7.1, -7.1, -7.2, -7.2, -7.3, -7.3, -7.3, -7.4, -7.4, -7.4, -7.4, -7.4, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.6, -7.6, -7.6, -7.6, -7.7, -7.7, -7.7, -7.7, -7.7, -7.8, -7.8, -7.7, -7.7, -7.7, -7.6, -7.6, -7.6, -7.6, -7.6, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.6, -7.5, -7.4, -7.4, -7.5, -7.6, -7.7, -7.8, -7.9, -7.9, -7.9, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -7.9, -7.9, -8.0, -8.0, -8.0, -8.1, -8.1, -8.1, -8.1, -8.2, -8.2, -8.2, -8.2, -8.2, -8.2, -8.2, -8.3, -8.4, -8.4, -8.4, -8.5, -8.5, -8.6, -8.6, -8.7, -8.7, -8.7, -8.8, -8.9, -9.0, -9.0, -9.1, -9.2, -9.3, -9.5, -9.6, -9.7, -9.7, -9.8, -9.8, -9.8, -9.9, -9.9, -10.0, -10.1, -10.2, -10.4, -10.6, -10.7, -10.8, -11.0, -11.3};

    double qqH8_mass [] = {80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 110.5, 111.0, 111.5, 112.0, 112.5, 113.0, 113.5, 114.0, 114.5, 115.0, 115.5, 116.0, 116.5, 117.0, 117.5, 118.0, 118.5, 119.0, 119.5, 120.0, 120.5, 121.0, 121.5, 122.0, 122.5, 123.0, 123.5, 124.0, 124.5, 125.0, 125.5, 126.0, 126.5, 127.0, 127.5, 128.0, 128.5, 129.0, 129.5, 130.0, 130.5, 131.0, 131.5, 132.0, 132.5, 133.0, 133.5, 134.0, 134.5, 135.0, 135.5, 136.0, 136.5, 137.0, 137.5, 138.0, 138.5, 139.0, 139.5, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 162.0, 164.0, 165.0, 166.0, 168.0, 170.0, 172.0, 174.0, 175.0, 176.0, 178.0, 180.0, 182.0, 184.0, 185.0, 186.0, 188.0, 190.0, 192.0, 194.0, 195.0, 196.0, 198.0, 200.0, 202.0, 204.0, 206.0, 208.0, 210.0, 212.0, 214.0, 216.0, 218.0, 220.0, 222.0, 224.0, 226.0, 228.0, 230.0, 232.0, 234.0, 236.0, 238.0, 240.0, 242.0, 244.0, 246.0, 248.0, 250.0, 252.0, 254.0, 256.0, 258.0, 260.0, 262.0, 264.0, 266.0, 268.0, 270.0, 272.0, 274.0, 276.0, 278.0, 280.0, 282.0, 284.0, 286.0, 288.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, 420.0, 440.0, 450.0, 460.0, 480.0, 500.0, 520.0, 540.0, 550.0, 560.0, 580.0, 600.0, 620.0, 640.0, 650.0, 660.0, 680.0, 700.0, 720.0, 740.0, 750.0, 760.0, 780.0, 800.0, 820.0, 840.0, 850.0, 860.0, 880.0, 900.0, 920.0, 940.0, 950.0, 960.0, 980.0, 1000.0};
    double qqH8_xsec [] = {2.410, 2.384, 2.360, 2.336, 2.311, 2.289, 2.265, 2.243, 2.221, 2.199, 2.176, 2.154, 2.133, 2.112, 2.090, 2.071, 2.050, 2.030, 2.010, 1.991, 1.971, 1.952, 1.934, 1.915, 1.897, 1.878, 1.860, 1.843, 1.826, 1.808, 1.791, 1.783, 1.775, 1.766, 1.758, 1.750, 1.742, 1.733, 1.725, 1.717, 1.709, 1.701, 1.693, 1.686, 1.678, 1.670, 1.661, 1.654, 1.647, 1.639, 1.632, 1.624, 1.617, 1.609, 1.602, 1.595, 1.588, 1.580, 1.573, 1.566, 1.559, 1.552, 1.544, 1.539, 1.531, 1.524, 1.517, 1.511, 1.504, 1.497, 1.490, 1.483, 1.477, 1.470, 1.463, 1.458, 1.451, 1.444, 1.439, 1.432, 1.425, 1.419, 1.413, 1.407, 1.401, 1.395, 1.388, 1.382, 1.376, 1.370, 1.365, 1.352, 1.341, 1.329, 1.317, 1.306, 1.295, 1.284, 1.272, 1.261, 1.251, 1.240, 1.229, 1.218, 1.208, 1.197, 1.187, 1.176, 1.166, 1.155, 1.146, 1.136, 1.123, 1.115, 1.106, 1.088, 1.070, 1.052, 1.035, 1.026, 1.017, 1.000, 0.9820, 0.9670, 0.9558, 0.9496, 0.9429, 0.9286, 0.9139, 0.8998, 0.8854, 0.8783, 0.8714, 0.8574, 0.8441, 0.8309, 0.8178, 0.8051, 0.7927, 0.7805, 0.7687, 0.7568, 0.7452, 0.7340, 0.7229, 0.7120, 0.7016, 0.6913, 0.6808, 0.6707, 0.6610, 0.6513, 0.6418, 0.6326, 0.6234, 0.6144, 0.6056, 0.5969, 0.5885, 0.5802, 0.5720, 0.5640, 0.5562, 0.5484, 0.5408, 0.5333, 0.5259, 0.5187, 0.5116, 0.5047, 0.4978, 0.4910, 0.4845, 0.4780, 0.4715, 0.4652, 0.4590, 0.4530, 0.4470, 0.4716, 0.4562, 0.4408, 0.4266, 0.4131, 0.3999, 0.3875, 0.3753, 0.3637, 0.3526, 0.3422, 0.3303, 0.3200, 0.3028, 0.2896, 0.2776, 0.2660, 0.2543, 0.2317, 0.2103, 0.2002, 0.1905, 0.1724, 0.1561, 0.1414, 0.1283, 0.1223, 0.1166, 0.1062, 0.09688, 0.08861, 0.08121, 0.07784, 0.07459, 0.06865, 0.06330, 0.05853, 0.05420, 0.05235, 0.05032, 0.04682, 0.04365, 0.04078, 0.03815, 0.03706, 0.03579, 0.03363, 0.03164, 0.02986, 0.02820, 0.02745, 0.02669, 0.02524, 0.02399};
    double qqH8_errp [] = {+2.9, +3.0, +2.9, +3.0, +2.9, +2.9, +2.9, +2.9, +2.8, +2.9, +2.9, +2.9, +2.8, +2.9, +2.9, +2.8, +2.9, +2.9, +2.9, +2.8, +2.8, +2.9, +2.7, +2.8, +2.8, +2.8, +2.8, +2.7, +2.7, +2.8, +2.7, +2.7, +2.7, +2.7, +2.7, +2.8, +2.8, +2.8, +2.7, +2.7, +2.7, +2.8, +2.8, +2.8, +2.7, +2.8, +2.9, +2.8, +2.9, +2.8, +2.8, +2.8, +2.8, +2.8, +2.9, +2.8, +2.8, +2.8, +2.9, +2.8, +2.8, +2.8, +2.9, +2.8, +2.9, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.7, +2.7, +2.8, +2.8, +2.7, +2.7, +2.7, +2.7, +2.7, +2.7, +2.7, +2.8, +2.7, +2.7, +2.8, +2.7, +2.8, +2.7, +2.7, +2.8, +2.7, +2.7, +2.7, +2.7, +2.7, +2.8, +2.7, +2.7, +2.7, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.9, +2.8, +2.7, +2.7, +2.7, +2.8, +2.8, +2.7, +2.8, +2.7, +2.8, +2.8, +2.8, +2.7, +2.8, +2.8, +2.8, +2.7, +2.7, +2.7, +2.7, +2.6, +2.7, +2.9, +2.8, +2.8, +2.8, +2.8, +2.8, +2.7, +2.8, +2.7, +2.8, +2.8, +2.7, +2.8, +2.7, +2.8, +2.8, +2.8, +2.7, +2.8, +2.8, +2.8, +2.8, +2.8, +2.7, +2.8, +2.8, +2.8, +2.7, +2.6, +2.7, +2.7, +2.8, +2.7, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.7, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.8, +2.9, +2.9, +2.9, +3.1, +3.1, +3.1, +3.2, +3.2, +3.4, +3.4, +3.5, +3.5, +3.5, +3.7, +3.9, +3.9, +4.0, +4.0, +4.2, +4.3, +4.4, +4.4, +4.5, +4.7, +4.8, +4.9, +5.0, +5.0, +5.2, +5.4, +5.5, +5.5, +5.6, +5.7, +5.8, +6.0, +6.1, +6.2, +6.3};
    double qqH8_errm [] = {-3.3, -3.3, -3.1, -3.0, -3.0, -3.0, -3.0, -3.2, -3.2, -3.2, -2.9, -2.9, -3.2, -2.9, -2.9, -3.1, -3.1, -3.0, -3.1, -3.1, -3.1, -3.1, -3.0, -3.0, -3.1, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -2.8, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -2.9, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.8, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.9, -2.8, -2.8, -2.8, -2.8, -2.8, -2.9, -2.9, -2.8, -2.8, -2.8, -2.9, -2.7, -2.8, -2.8, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.6, -2.6, -2.6, -2.6, -2.8, -2.6, -2.8, -2.8, -2.8, -2.8, -2.8, -2.7, -2.8, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.7, -2.7, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.9, -2.8, -2.8, -2.8, -2.9, -3.0, -3.0, -3.1, -3.2, -3.2, -3.2, -3.3, -3.6, -3.6, -3.6, -3.7, -3.9, -4.1, -4.1, -4.1, -4.3, -4.3, -4.4, -4.6, -4.8, -4.8, -4.9, -5.0, -5.3, -5.2, -5.2, -5.4, -5.6, -5.7, -6.0, -6.1, -6.0, -6.2, -6.3, -6.5, -6.8, -6.9, -6.8, -7.0, -7.1, -7.2};
    double qqH8_scap [] = {+0.2, +0.4, +0.3, +0.4, +0.3, +0.3, +0.3, +0.3, +0.2, +0.3, +0.3, +0.3, +0.2, +0.3, +0.3, +0.2, +0.3, +0.3, +0.3, +0.2, +0.2, +0.3, +0.2, +0.2, +0.2, +0.3, +0.3, +0.2, +0.2, +0.3, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.3, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.3, +0.2, +0.3, +0.2, +0.2, +0.2, +0.2, +0.2, +0.3, +0.2, +0.2, +0.2, +0.3, +0.2, +0.2, +0.2, +0.3, +0.2, +0.3, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.3, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.3, +0.2, +0.2, +0.3, +0.2, +0.3, +0.2, +0.2, +0.3, +0.2, +0.2, +0.2, +0.2, +0.2, +0.3, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.3, +0.3, +0.2, +0.2, +0.2, +0.3, +0.3, +0.2, +0.3, +0.2, +0.3, +0.3, +0.3, +0.2, +0.3, +0.3, +0.3, +0.2, +0.2, +0.2, +0.2, +0.2, +0.2, +0.3, +0.2, +0.2, +0.3, +0.3, +0.3, +0.2, +0.3, +0.2, +0.3, +0.3, +0.2, +0.3, +0.2, +0.3, +0.3, +0.3, +0.2, +0.3, +0.3, +0.3, +0.3, +0.3, +0.2, +0.3, +0.3, +0.3, +0.3, +0.2, +0.3, +0.3, +0.3, +0.2, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.3, +0.4, +0.4, +0.4, +0.4, +0.4, +0.4, +0.4, +0.4, +0.4, +0.4, +0.5, +0.5, +0.5, +0.5, +0.6, +0.6, +0.6, +0.6, +0.6, +0.6, +0.6, +0.6, +0.7, +0.7, +0.7};
    double qqH8_scam [] = {-0.3, -0.3, -0.3, -0.2, -0.2, -0.2, -0.2, -0.3, -0.3, -0.3, -0.2, -0.2, -0.3, -0.2, -0.2, -0.2, -0.2, -0.3, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.1, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.1, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.1, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.1, -0.1, -0.1, -0.1, -0.1, -0.2, -0.2, -0.1, -0.1, -0.1, -0.2, -0.1, -0.2, -0.2, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.4, -0.4, -0.4, -0.4, -0.5, -0.5, -0.5, -0.5, -0.6, -0.6, -0.6, -0.6, -0.7, -0.7, -0.7, -0.8, -0.8, -0.8, -0.8, -0.8, -0.9, -0.9, -0.9, -0.9, -1.0, -1.0, -1.1, -1.0, -1.1, -1.1, -1.1, -1.2, -1.2, -1.2, -1.3, -1.2, -1.3, -1.3};
    double qqH8_pdfp [] = {+2.7, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.5, +2.6, +2.6, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.6, +2.5, +2.6, +2.5, +2.5, +2.5, +2.6, +2.6, +2.6, +2.5, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.5, +2.6, +2.5, +2.5, +2.6, +2.6, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.6, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.4, +2.5, +2.6, +2.6, +2.6, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.4, +2.4, +2.4, +2.4, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.4, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.5, +2.6, +2.6, +2.6, +2.8, +2.8, +2.8, +2.9, +2.9, +3.1, +3.1, +3.2, +3.2, +3.2, +3.4, +3.5, +3.5, +3.6, +3.6, +3.8, +3.9, +4.0, +4.0, +4.1, +4.3, +4.3, +4.4, +4.5, +4.5, +4.6, +4.8, +4.9, +4.9, +5.0, +5.1, +5.2, +5.4, +5.4, +5.5, +5.6};
    double qqH8_pdfm [] = {-3.0, -3.0, -2.8, -2.8, -2.8, -2.8, -2.8, -2.9, -2.9, -2.9, -2.7, -2.7, -2.9, -2.7, -2.7, -2.9, -2.9, -2.7, -2.9, -2.9, -2.9, -2.9, -2.8, -2.8, -2.9, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.6, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.8, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.7, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.5, -2.5, -2.5, -2.5, -2.7, -2.5, -2.7, -2.7, -2.7, -2.7, -2.7, -2.6, -2.7, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.6, -2.5, -2.5, -2.5, -2.6, -2.7, -2.7, -2.8, -2.8, -2.8, -2.8, -2.9, -3.1, -3.1, -3.1, -3.2, -3.3, -3.5, -3.5, -3.5, -3.6, -3.6, -3.7, -3.8, -4.0, -4.0, -4.1, -4.2, -4.4, -4.3, -4.3, -4.5, -4.6, -4.7, -4.9, -5.1, -4.9, -5.1, -5.2, -5.3, -5.6, -5.7, -5.5, -5.8, -5.8, -5.9};


    ggH7TG_xsec = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_xsec);
    ggH7TG_errp = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_errp);
    ggH7TG_errm = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_errm);
    ggH7TG_scap = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_scap);
    ggH7TG_scam = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_scam);
    ggH7TG_pdfp = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_pdfp);
    ggH7TG_pdfm = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_pdfm);

    qqH7TG_xsec = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_xsec);
    qqH7TG_errp = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_errp);
    qqH7TG_errm = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_errm);
    qqH7TG_scap = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_scap);
    qqH7TG_scam = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_scam);
    qqH7TG_pdfp = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_pdfp);
    qqH7TG_pdfm = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_pdfm);

    ggH8TG_xsec = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_xsec);
    ggH8TG_errp = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_errp);
    ggH8TG_errm = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_errm);
    ggH8TG_scap = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_scap);
    ggH8TG_scam = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_scam);
    ggH8TG_pdfp = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_pdfp);
    ggH8TG_pdfm = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_pdfm);

    qqH8TG_xsec = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_xsec);
    qqH8TG_errp = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_errp);
    qqH8TG_errm = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_errm);
    qqH8TG_scap = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_scap);
    qqH8TG_scam = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_scam);
    qqH8TG_pdfp = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_pdfp);
    qqH8TG_pdfm = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_pdfm);

    //#FIXME: extrapolated from 600 to 1TeVmissing points from 650GeV to 1TeV
    double QCDScaleMass   [] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000};
    double QCDScaleK0ggH0 [] = {1.15 , 1.16, 1.17, 1.20, 1.17, 1.19, 1.22, 1.24, 1.25, 1.25, 1.25, 1.25, 1.25};
    double QCDScaleK0ggH1 [] = {0.88, 0.86, 0.84, 0.83, 0.82, 0.81, 0.80, 0.78, 0.78, 0.78, 0.78, 0.78, 0.78};
    double QCDScaleK1ggH1 [] = {1.27, 1.27, 1.27, 1.27, 1.26, 1.26, 1.25, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26};
    double QCDScaleK1ggH2 [] = {0.96, 0.96, 0.95, 0.95, 0.95, 0.95, 0.95,  0.95, 0.94, 0.94, 0.94, 0.94, 0.94};
    double QCDScaleK2ggH2 [] = { 1.20, 1.17, 1.20, 1.21, 1.20, 1.20, 1.17, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19};

    double UEPSf0 []         = {0.952, 0.955, 0.958, 0.964, 0.966, 0.954, 0.946, 0.931, 0.920, 0.920, 0.920, 0.920, 0.920};
    double UEPSf1 []         = {1.055, 1.058, 1.061, 1.068, 1.078, 1.092, 1.102, 1.117, 1.121, 1.121, 1.121, 1.121, 1.121};
    double UEPSf2 []         = {0.059, 0.990, 0.942, 0.889, 0.856, 0.864, 0.868, 0.861, 0.872, 0.872, 0.872, 0.872, 0.872};

    TG_QCDScaleK0ggH0 = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleK0ggH0);
    TG_QCDScaleK0ggH1 = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleK0ggH1);
    TG_QCDScaleK1ggH1 = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleK1ggH1);
    TG_QCDScaleK1ggH2 = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleK1ggH2);
    TG_QCDScaleK2ggH2 = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleK2ggH2);

    TG_UEPSf0         = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, UEPSf0);
    TG_UEPSf1         = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, UEPSf1);
    TG_UEPSf2         = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, UEPSf2);



}



void initializeAtgcMaps()
{

    // Z pt bins -- 8 TeV only, full 2021 dataset
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.112318;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.878952;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 1.24356; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.470437;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.307615;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 1.83817;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 1.45612; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.641511;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.610903;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 3.12009;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 1.71705; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.83278;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 1.02218;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 4.72473;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 2.0137; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 1.04065;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 1.54145;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 6.65207;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 2.33927; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 1.26489;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 2.16871;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 8.90213;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 2.69043; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 1.50637;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 2.90396;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 11.4749;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 3.06586; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 1.76626;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 3.7472;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 14.3704;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 3.46527; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 2.04575;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 4.69843;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 17.5885;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 3.88897; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 2.34594;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 6.92486;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 24.993;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 4.81169; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 3.01225;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 9.58326;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 33.6883;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 5.83979; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 3.77152;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 12.6736;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 43.6745;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 6.97911; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 4.6283;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 16.1959;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 54.9515;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 8.23473; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 5.58582;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 20.1502;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 67.5193;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 9.61082; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 6.64635;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.594503;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 1.23644; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.456885;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.173411;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 1.44044;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 1.45455; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.620225;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.496276;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 2.6233;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 1.72136; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.803765;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.949591;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 4.14307;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 2.024; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 1.00393;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 1.53336;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 5.99976;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 2.35568; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 1.22046;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 2.24758;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 8.19337;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 2.71318; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 1.45419;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 3.09225;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 10.7239;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 3.09526; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 1.70624;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 4.06737;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 13.5913;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 3.50172; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 1.97775;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 5.17294;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 16.7957;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 3.93294; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 2.2698;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 7.77544;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 24.2152;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 4.87234; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 2.91913;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 10.8997;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 32.9823;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 5.91955; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 3.6603;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 14.5459;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 43.0972;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 7.08062; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 4.49761;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 18.7138;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 54.5597;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 8.36075; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 5.4341;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 23.4035;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 67.3698;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 9.7642; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 6.47191;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.0353875;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.260686;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 1.29725; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.454409;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.275214;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.965611;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 1.53447; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.616209;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.67853;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 2.02718;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 1.82364; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.798335;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 1.24533;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 3.4454;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 2.15097; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.997232;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 1.97563;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 5.22026;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 2.5094; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 1.21266;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 2.86941;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 7.35177;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 2.89569; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 1.44543;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 3.92668;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 9.83993;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 3.30869; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 1.69669;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 5.14744;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 12.6847;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 3.74832; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 1.96757;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 6.53169;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 15.8862;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 4.21508; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 2.25914;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 9.79065;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 23.359;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 5.23316; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 2.90796;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 13.7036;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 32.2584;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 6.36985; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 3.64914;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 18.2705;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 42.5844;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 7.63181; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 4.48693;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 23.4913;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 54.337;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 9.0247; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 5.42432;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 29.3661;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 67.5162;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 10.5531; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 6.46339;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.0285903;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.173812;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 1.29656; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.417746;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.207082;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.793067;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 1.52228; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.556224;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.584414;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 1.73577;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 1.79937; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.714352;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 1.16058;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 3.00192;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 2.11468; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.888449;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 1.93559;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 4.59153;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 2.46128; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 1.07799;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 2.90944;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 6.50458;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 2.8359; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 1.28353;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 4.08213;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 8.74108;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 3.23731; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 1.50598;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 5.45365;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 11.301;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 3.66538; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 1.7463;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 7.02402;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 14.1844;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 4.12053; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 2.00538;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 10.7613;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 20.9216;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 5.11499; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 2.58291;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 15.2939;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 28.9525;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 6.22721; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 3.24371;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 20.6218;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 38.2773;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 7.46355; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 3.99143;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 26.7451;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 48.8958;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 8.82938; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 4.82861;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 33.6638;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 60.8081;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 10.3291; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 5.75704;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.705302;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 1.26306; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.456471;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 1.54134;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 1.47882; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.616886;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.169235;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 2.70227;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 1.74422; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.798229;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.527363;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 4.18808;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 2.04663; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.996895;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 1.05117;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 5.99877;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 2.3793; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 1.21262;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 1.74065;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 8.13434;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 2.73901; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 1.44623;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 2.59581;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 10.5948;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 3.12449; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 1.69886;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 3.61664;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 13.3801;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 3.53556; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 1.97165;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 4.80315;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 16.4903;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 3.97259; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 2.26567;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 7.67321;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 23.6854;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 4.92722; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 2.92098;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 11.206;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 32.18;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 5.99452; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 3.6707;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 15.4014;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 41.9742;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 7.18051; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 4.51904;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 20.2596;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 53.0678;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 8.49035; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 5.46889;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 25.7805;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 65.461;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 9.92825; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 6.52228;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.721736;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 1.24599; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.453895;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.0580536;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 1.57893;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 1.45494; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.614676;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.303323;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 2.76261;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 1.71224; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.795639;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.672913;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 4.27279;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 2.0055; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.993152;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 1.16682;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 6.10947;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 2.32798; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 1.20691;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 1.78506;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 8.27265;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 2.67641; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 1.4377;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 2.52761;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 10.7623;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 3.04944; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 1.68661;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 3.39449;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 13.5785;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 3.44679; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 1.95477;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 4.38569;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 16.7212;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 3.86875; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 2.24323;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 6.74104;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 23.986;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 4.78897; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 2.88461;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 9.59369;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 32.5568;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 5.81579; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 3.61672;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 12.9436;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 42.4336;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 6.95498; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 4.44382;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 16.7908;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 53.6164;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 8.21154; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 5.3689;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 21.1353;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 66.1051;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 9.58957; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 6.39405;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.417665;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.554021;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 1.31556; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.4517;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.694369;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 1.25369;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 1.54895; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.608608;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 1.13386;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 2.2714;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 1.83458; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.786017;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 1.73615;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 3.60714;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 2.15881; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.980199;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 2.50123;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 5.26092;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 2.51452; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 1.19078;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 3.42909;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 7.23273;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 2.89836; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 1.41848;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 4.51975;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 9.52258;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 3.30911; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 1.66437;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 5.7732;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 12.1305;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 3.74662; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 1.92955;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 7.18944;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 15.0564;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 4.21135; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 2.21504;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 10.5103;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 21.8623;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 5.22544; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 2.85044;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 14.4823;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 29.9404;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 6.35811; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 3.57639;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 19.1055;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 39.2907;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 7.6159; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 4.39704;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 24.3798;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 49.9131;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 9.00437; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 5.3153;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 30.3053;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 61.8076;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 10.5281; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 6.33318;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.223032;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.516882;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 1.28333; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.427962;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.55704;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 1.17178;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 1.51856; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.568555;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 1.03541;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 2.11937;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 1.80512; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.729184;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 1.65815;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 3.35966;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 2.12932; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.905991;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 2.42525;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 4.89263;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 2.48411; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 1.09837;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 3.33672;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 6.71829;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 2.86629; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 1.30684;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 4.39255;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 8.83665;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 3.2747; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 1.5323;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 5.59275;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 11.2477;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 3.70925; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 1.7757;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 6.93731;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 13.9514;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 4.17043; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 2.03794;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 10.0595;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 20.237;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 5.17579; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 2.62211;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 13.7592;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 27.6933;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 6.29762; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 3.29;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 18.0363;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 36.3204;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 7.54248; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 4.04536;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 22.8909;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 46.1182;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 8.91601; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 4.8908;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 28.323;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 57.0868;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 10.4227; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 5.82815;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.202599;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 1.22436; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.408191;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.639873;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 1.43329; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.534316;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.122733;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 1.35625;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 1.68967; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.680402;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.378289;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 2.35174;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 1.98099; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.84262;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.722256;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 3.62634;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 2.3005; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 1.02023;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 1.15463;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 5.18005;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 2.64489; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 1.2136;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 1.67542;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 7.01286;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 3.01281; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 1.42352;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 2.28462;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 9.12478;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 3.40397; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 1.65083;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 2.98223;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 11.5158;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 3.81861; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 1.89634;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 4.64268;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 17.1352;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 4.72076; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 2.44475;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 6.65677;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 23.871;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 5.72483; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 3.07341;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 9.0245;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 31.7233;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 6.83653; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 3.78563;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 11.7459;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 40.6919;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 8.06085; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 4.5837;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 14.8209;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 50.777;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 9.4019; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 5.46924;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.612365;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 1.25564; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.431486;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 1.35577;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 1.45989; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.577503;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 2.36931;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 1.71229; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.742988;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.0442844;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 3.653;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 2.0007; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.924123;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.380655;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 5.20682;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 2.31838; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 1.12035;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.852039;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 7.03079;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 2.66199; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 1.33223;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 1.45844;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 9.1249;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 3.03011; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 1.56069;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 2.19985;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 11.4892;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 3.4224; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 1.8067;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 3.07628;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 14.1235;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 3.83908; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 2.0712;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 5.23417;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 20.2028;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 4.74794; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 2.65889;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 7.93213;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 27.3625;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 5.76213; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 3.32919;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 11.1701;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 35.6029;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 6.88724; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 4.08602;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 14.9482;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 44.9238;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 8.12816; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 4.93215;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 19.2663;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 55.3253;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 9.48894; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 5.86955;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.256991;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 1.29312; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.410784;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.239211;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.793289;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 1.5168; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.543826;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.631054;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 1.54664;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 1.79083; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.695362;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 1.11675;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 2.51704;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 2.10185; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.86143;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 1.69631;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 3.70449;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 2.44271; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 1.04124;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 2.36972;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 5.10899;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 2.80994; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 1.23513;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 3.13698;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 6.73055;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 3.20214; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 1.44385;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 3.9981;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 8.56915;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 3.61902; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 1.66824;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 4.95308;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 10.6248;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 4.06089; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 1.90912;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 7.1446;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 15.3873;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 5.0222; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 2.44325;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 9.71155;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 21.0179;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 6.0921; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 3.05119;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 12.6539;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 27.5168;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 7.27669; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 3.73659;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 15.9717;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 34.8839;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 8.58133; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 4.50204;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 19.6649;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 43.1192;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 10.0104; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 5.34943;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 0.936655;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.28168;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 1.32026; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.427146;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 1.37099;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.76905;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 1.55117; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.56857;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 1.88619;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 1.52176;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 1.83366; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.729763;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 2.48227;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 2.53981;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 2.15394; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.906883;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 3.15922;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 3.8232;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 2.50469; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 1.09933;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 3.91705;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 5.37194;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 2.88238; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 1.30763;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 4.75574;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 7.18601;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 3.28562; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 1.53266;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 5.67531;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 9.26542;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 3.71411; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 1.7754;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 6.67575;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 11.6102;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 4.16821; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 2.03675;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 8.91926;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 17.0957;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 5.15588; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 2.61841;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 11.4863;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 23.6426;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 6.25492; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 3.28288;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 14.3767;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 31.2509;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 7.47163; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 4.03395;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 17.5907;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 39.9205;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 8.81152; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 4.87427;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 21.1282;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 49.6515;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 10.2791; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 5.8057;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.155385;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.602321;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 1.25494; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.43716;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.346607;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 1.23805;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 1.46193; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.585515;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.614193;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 2.10346;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 1.71685; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.753225;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.958142;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 3.19855;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 2.00722; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.936367;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 1.37845;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 4.52332;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 2.32613; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 1.13433;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 1.87513;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 6.07777;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 2.67013; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 1.34765;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 2.44817;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 7.8619;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 3.03775; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 1.57722;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 3.09757;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 9.87571;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 3.42858; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 1.82403;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 3.82334;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 12.1192;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 3.84281; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 2.089;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 5.50396;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 17.2952;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 4.74367; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 2.67672;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 7.49004;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 23.39;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 5.74566; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 3.34587;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 9.78157;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 30.4034;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 6.85434; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 4.10049;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 12.3786;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 38.3356;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 8.07465; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 4.94343;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 15.281;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 47.1865;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 9.41072; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 5.87675;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.355385;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 1.23934; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.414287;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.859004;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 1.43655; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.549069;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 1.60716;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 1.68088; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.702916;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 2.59986;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 1.96058; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.87196;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.0211238;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 3.83709;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 2.26907; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 1.0555;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.410785;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 5.31886;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 2.60302; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 1.25396;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 0.930173;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 7.04517;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 2.96099; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 1.46814;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 1.57929;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 9.01602;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 3.3426; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 1.69893;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 2.35813;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 11.2314;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 3.74804; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 1.94718;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 4.30499;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 16.3958;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 4.63254; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 2.49904;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 6.77077;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 22.5383;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 5.61967; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 3.12873;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 9.75545;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 29.659;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 6.71481; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 3.83989;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 13.259;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 37.7579;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 7.92267; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 4.63509;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 17.2815;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 46.8349;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 9.24716; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 5.51616;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.0231774;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.336845;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 1.27444; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.40945;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.118325;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.850446;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 1.48785; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.539933;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.334315;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 1.59806;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 1.75057; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.689553;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.671147;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 2.57969;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 2.04989; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.854425;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 1.12882;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 3.79534;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 2.3789; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 1.0338;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 1.70733;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 5.24499;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 2.7342; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 1.22805;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 2.40669;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 6.92867;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 3.1144; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 1.43795;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 3.22689;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 8.84636;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 3.51919; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 1.66436;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 4.16793;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 10.9981;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 3.94883; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 1.90809;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 6.41254;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 16.0035;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 4.88516; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 2.45043;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 9.14051;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 21.945;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 5.92915; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 3.0698;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 12.3518;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 28.8226;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 7.08665; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 3.76973;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 16.0466;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 36.6362;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 8.36278; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 4.55269;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 20.2246;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 45.3859;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 9.76174; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 5.42044;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.263553;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.310829;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 1.27673; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.426782;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.520741;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.853214;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 1.49289; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.569913;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 0.869875;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 1.61759;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 1.7585; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.731586;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 1.31095;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 2.60395;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 2.06065; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.907793;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 1.84397;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 3.81229;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 2.39234; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 1.0978;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 2.46894;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 5.24263;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 2.75016; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 1.30202;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 3.18585;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 6.89495;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 3.1327; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 1.52125;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 3.99471;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 8.76927;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 3.53964; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 1.7564;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 4.89551;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 10.8656;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 3.97126; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 2.00835;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 6.97295;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 15.7241;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 4.91101; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 2.56572;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 9.41816;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 21.4706;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 5.95776; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 3.19868;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 12.2312;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 28.1051;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 7.11741; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 3.91114;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 15.4119;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 35.6275;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 8.39514; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 4.70597;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 18.9605;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 44.0379;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 9.79521; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 5.5852;


    // Z pt bins
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.131416;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.0312657;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.213479; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.0337431;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.193208;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.0594768;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.245928; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.0487644;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.271561;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.0961302;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.286418; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.0649553;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.366476;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.141226;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.333047; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.0822702;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 0.477953;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.194764;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 0.384728; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.100807;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 0.605991;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.256744;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 0.440902; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.120693;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 0.750591;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 0.327166;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 0.501327; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 0.142055;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 0.911752;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 0.406031;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 0.565935; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 0.16501;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 1.08947;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 0.493338;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 0.634757; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 0.189658;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 1.49461;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 0.693279;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 0.785395; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 0.244369;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 1.96598;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 0.926989;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 0.954102; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 0.306734;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 2.5036;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 1.19447;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 1.14177; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 0.377137;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 3.10747;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 1.49572;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 1.34917; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 0.455844;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 3.77759;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 1.83073;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 1.57694; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 0.54304;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.214305;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.0385704;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.215557; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.0390813;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.31544;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.0684303;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.249662; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.0571816;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.425227;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.105675;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.291968; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.076468;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.543665;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.150304;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.340454; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.096946;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 0.670755;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.202318;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 0.393985; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.118748;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 0.806497;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.261716;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 0.451988; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.142029;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 0.950891;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 0.328499;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 0.514217; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 0.16694;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 1.10394;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 0.402667;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 0.580608; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 0.193616;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 1.26563;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 0.484219;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 0.651196; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 0.222179;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 1.61498;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 0.669478;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 0.805335; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 0.285359;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 1.99894;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 0.884275;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 0.97754; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 0.357134;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 2.4175;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 1.12861;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 1.16874; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 0.43797;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 2.87067;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 1.40248;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 1.37976; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 0.528196;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 3.35845;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 1.7059;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 1.61127; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 0.628042;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.0691866;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.0205375;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.238042; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.0536722;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.130584;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.060439;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.273494; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.0768592;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.213621;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.121073;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.317805; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.101982;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.318299;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.20244;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.36889; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.1289;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 0.444617;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.30454;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 0.425539; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.157737;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 0.592576;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.427372;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 0.48712; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.188678;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 0.762174;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 0.570938;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 0.553351; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 0.221913;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 0.953413;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 0.735236;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 0.624146; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 0.257616;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 1.16629;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 0.920266;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 0.699529; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 0.295945;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 1.65697;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 1.35253;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 0.864422; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 0.380994;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 2.23421;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 1.86772;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 1.04894; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 0.477903;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 2.89801;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 2.46584;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 1.25405; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 0.587271;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 3.64838;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 3.14689;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 1.48059; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 0.709514;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 4.4853;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 3.91088;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 1.72927; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 0.844923;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.0500644;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.0359716;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.237463; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.0495965;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.103945;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.0788138;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.272569; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.0702536;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.180154;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.138632;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.316468; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.0927924;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.278692;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.215427;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.367088; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.117018;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 0.399559;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.309197;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 0.423221; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.143012;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 0.542755;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.419944;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 0.484233; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.170927;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 0.708279;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 0.547668;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 0.549835; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 0.200928;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 0.896133;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 0.692367;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 0.619938; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 0.233169;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 1.10631;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 0.854042;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 0.69456; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 0.267788;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 1.59367;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 1.22832;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 0.857711; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 0.344621;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 2.17033;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 1.67051;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 1.04018; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 0.432181;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 2.83631;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 2.1806;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 1.2429; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 0.531005;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 3.59161;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 2.75859;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 1.46673; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 0.641468;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 4.43622;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 3.40449;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 1.71236; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 0.763831;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.);
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.210272; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.0305842;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.00912376;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.240892; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.0438231;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.00979503;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.0407356;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.279239; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.0583541;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.0643334;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.0865887;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.323503; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.0741299;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 0.138336;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.146683;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 0.37262; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.0912491;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 0.231804;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.221019;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 0.426026; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.109839;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 0.344737;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 0.309596;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 0.483461; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 0.130024;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 0.477134;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 0.412415;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 0.544839; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 0.151913;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 0.628995;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 0.529475;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 0.610172; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 0.175599;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 0.991113;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 0.806319;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 0.752995; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 0.228661;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 1.43109;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 1.14013;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 0.912691; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 0.289681;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 1.94892;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 1.5309;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 1.09008; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 0.358972;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 2.54462;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 1.97864;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 1.2859; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 0.436744;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 3.21817;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 2.48335;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 1.50075; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 0.523139;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.0728393;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.21273; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.0275605;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.112127;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.0118336;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.244754; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.0388371;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.167147;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.0336879;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.284711; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.0512897;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.2379;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.064676;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.330708; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.0648119;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 0.324384;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.104798;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 0.381649; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.0794571;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 0.426602;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.154053;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 0.436966; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.0953181;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 0.544551;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 0.212442;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 0.496404; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 0.112492;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 0.678233;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 0.279965;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 0.559885; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 0.131069;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 0.827647;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 0.356621;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 0.627432; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 0.151126;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 1.17367;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 0.537335;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 0.775045; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 0.195936;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 1.58263;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 0.754583;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 0.940064; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 0.247331;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 2.05451;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 1.00837;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 1.12335; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 0.30559;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 2.58933;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 1.29868;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 1.32568; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 0.370902;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 3.18707;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 1.62554;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 1.54769; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 0.443397;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.164083;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.0649986;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.238695; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.0507416;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.253631;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.122659;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.274448; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.0720588;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.357366;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.19722;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.319059; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.0952041;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.475289;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.288681;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.370401; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.11998;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 0.6074;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.397041;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 0.427241; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.146466;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 0.753698;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.522302;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 0.488934; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.174812;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 0.914183;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 0.664462;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 0.555186; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 0.205182;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 1.08886;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 0.823522;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 0.625904; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 0.237731;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 1.27772;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 0.999483;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 0.701108; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 0.272598;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 1.698;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 1.4021;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 0.865323; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 0.349758;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 2.17503;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 1.87232;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 1.04873; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 0.437439;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 2.70882;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 2.41014;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 1.25228; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 0.536202;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 3.29935;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 3.01556;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 1.47684; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 0.646443;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 3.94664;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 3.68858;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 1.72312; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 0.768444;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.0755061;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.237816; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.0397883;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.145384;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.00425104;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.273114; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.0541243;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.232395;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.0547986;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.317196; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.0704281;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.336537;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.12951;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.367964; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.0884359;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 0.457812;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.228385;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 0.424191; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.108163;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 0.59622;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.351424;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 0.485234; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.129709;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 0.751759;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 0.498626;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 0.550796; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 0.153189;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 0.924431;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 0.669993;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 0.620782; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 0.178716;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 1.11424;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 0.865523;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 0.695206; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 0.206387;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 1.54524;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 1.32907;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 0.857709; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 0.268487;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 2.04478;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 1.88928;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 1.03918; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 0.340001;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 2.61284;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 2.54614;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 1.24056; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 0.421277;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 3.24943;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 3.29966;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 1.46271; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 0.512548;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 3.95456;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 4.14983;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 1.70632; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 0.61397;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.0704047;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.0187349;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.212635; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.0323029;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.128747;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.0361378;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.244697; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.0464244;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.197419;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.0593982;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.28466; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.0616554;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.276422;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.0885163;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.330614; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.077921;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.365754;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.123492;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.381451; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.0952968;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 0.465416;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.164325;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 0.436594; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.113894;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 0.575408;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.211016;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 0.49578; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.133825;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 0.695731;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 0.263564;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 0.558928; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 0.155196;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 0.826383;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 0.32197;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 0.626054; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 0.1781;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 1.11868;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 0.456355;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 0.772553; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 0.228821;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 1.45229;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 0.614169;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 0.936079; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 0.286501;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 1.82723;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 0.795414;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 1.11749; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 0.351507;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 2.24348;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 1.00009;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 1.31757; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 0.424097;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 2.70105;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 1.2282;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 1.53693; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 0.504452;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 0.0809639;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.0504655;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 0.209985; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.0352282;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 0.113445;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.082978;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 0.239826; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.0510316;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 0.155347;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.120408;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 0.277253; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.0679051;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.206671;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.162756;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.32049; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.0857908;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.267415;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.210022;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.368472; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.104774;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.33758;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.262206;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.420624; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.124974;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 0.417167;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.319308;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 0.476671; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.14651;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 0.506175;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 0.381327;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 0.536513; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 0.169496;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 0.604603;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 0.448264;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 0.60015; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 0.194032;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 0.829724;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 0.596891;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 0.739063; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 0.248096;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 1.09253;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 0.76519;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 0.894107; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 0.309271;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 1.39302;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 0.95316;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 1.06607; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 0.377973;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 1.73119;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 1.1608;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 1.25566; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 0.4545;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 2.10705;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 1.38811;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 1.46347; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 0.539067;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.0640145;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.0230408;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.237101; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.041279;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.11916;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.0552653;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.271824; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.0564886;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.188629;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.100784;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.315222; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.0734742;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.272422;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.159596;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.365216; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.0919355;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.370538;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.231702;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.420581; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.111864;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 0.482979;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.317103;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 0.480665; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.133343;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 0.609743;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.415797;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 0.545161; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.156479;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 0.75083;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 0.527785;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 0.613965; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 0.181382;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 0.906242;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 0.653067;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 0.687083; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 0.208149;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 1.26004;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 0.943513;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 0.846571; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 0.267619;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 1.67113;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 1.28714;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 1.02446; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 0.335452;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 2.13951;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 1.68393;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 1.22165; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 0.412049;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 2.66519;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 2.13391;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 1.439; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 0.497693;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 3.24816;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 2.63706;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 1.67719; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 0.592583;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.0179958;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 0.233301; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.0475204;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.0531722;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 0.265468; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.0668626;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 0.0325123;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.106297;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 0.30593; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.0880645;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.0882914;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.177371;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.352771; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.110904;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.160689;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.266394;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.404823; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.135443;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.249706;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.373365;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.461445; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.161818;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 0.355341;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.498285;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 0.522321; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.190179;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 0.477595;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 0.641153;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 0.587329; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 0.220671;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 0.616468;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 0.80197;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 0.656459; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 0.253423;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 0.94407;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 1.17745;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 0.80733; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 0.326135;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 1.33815;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 1.62473;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 0.975655; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 0.409026;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 1.7987;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 2.1438;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 1.16227; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 0.502598;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 2.32573;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 2.73466;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 1.36793; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 0.607204;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 2.91923;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 3.39732;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 1.59329; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 0.723089;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.012709;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.0639058;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.207256; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.0331005;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.0123315;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.100103;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.235151; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.0476135;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.0257846;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.138914;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.270363; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.0631394;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.0530683;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.180339;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.311251; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.0795853;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 0.0941827;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.224377;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 0.356808; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.0970107;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 0.149128;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.271029;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 0.406473; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.115514;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 0.217903;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.320295;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 0.459971; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.135198;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 0.30051;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 0.372175;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 0.517194; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 0.156163;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 0.396947;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 0.426668;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 0.578134; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 0.178499;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 0.631313;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 0.543495;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 0.711378; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 0.227593;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 0.921001;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 0.670778;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 0.860329; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 0.283;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 1.26601;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 0.808515;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 1.02571; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 0.345104;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 1.66635;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 0.956707;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 1.20819; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 0.414188;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 2.122;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 1.11535;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 1.40834; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 0.490453;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.100076;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.0260957;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.211794; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.0289609;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.139376;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.0470015;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.242814; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.0409962;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 0.191274;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.0729548;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 0.28163; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.0540664;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 0.255771;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.103956;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 0.326408; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.0680437;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.332865;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.140004;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.376071; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.0829637;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.422559;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.1811;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.430049; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.0989062;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 0.52485;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.227243;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 0.48808; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.11596;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 0.63974;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 0.278434;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 0.550078; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 0.134211;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 0.767229;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 0.334673;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 0.616057; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 0.153735;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 1.06;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 0.462292;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 0.760253; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 0.19687;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 1.40317;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 0.610102;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 0.921438; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 0.245804;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 1.79673;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 0.778101;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 1.10045; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 0.300856;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 2.24068;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 0.966291;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 1.29802; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 0.362253;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 2.73502;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 1.17467;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 1.51477; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 0.430159;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.0571621;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.0444818;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.233102; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.0479436;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.0758366;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.0841477;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.264524; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.0674993;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.10885;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.136266;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.304157; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.0888473;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.156203;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.200838;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.350139; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.111753;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 0.217896;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.277863;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 0.401325; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.136266;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 0.293927;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.36734;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 0.457079; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.162516;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 0.384298;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.46927;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 0.517084; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.190647;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 0.489008;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 0.583653;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 0.581215; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 0.220801;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 0.608058;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 0.71049;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 0.649459; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 0.253103;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 0.889174;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 1.00152;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 0.798515; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 0.324585;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 1.22765;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 1.34236;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 0.964943; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 0.40581;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 1.62348;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 1.73302;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 1.14955; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 0.497295;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 2.07667;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 2.17348;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 1.3531; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 0.599406;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 2.58721;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 2.66376;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 1.5762; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 0.712405;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.0721352;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.0583883;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.237046; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.0462938;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.145436;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.102263;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.271733; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.0647455;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 0.227922;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.156201;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 0.315067; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.0849138;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 0.319594;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.220202;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 0.364962; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.106528;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.420451;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.294265;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.420186; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.129608;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.530494;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.378391;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.480081; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.15426;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 0.649723;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.472579;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 0.544335; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.180611;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 0.778138;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 0.576831;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 0.612838; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 0.208787;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 0.915738;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 0.691144;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 0.685594; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 0.238905;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 1.2185;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 0.94996;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 0.844161; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 0.305366;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 1.558;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 1.24903;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 1.02085; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 0.38067;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 1.93424;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 1.58834;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 1.21658; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 0.46531;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 2.34722;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 1.96791;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 1.43217; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 0.559645;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.);
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 2.79695;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 2.38773;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 1.66833; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 0.663928;

    /*
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.0699713; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.547566;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.774709; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.293071;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.191637; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 1.14514;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.907126; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.399646;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.380578; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 1.94374;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 1.06968; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.518802;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.636794; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 2.94339;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 1.25449; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.648297;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 0.960286; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 4.14408;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 1.45731; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.787995;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 1.35105; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 5.54581;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 1.67607; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.938433;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 1.8091; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 7.14858;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 1.90996; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 1.10034;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 2.33441; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 8.95239;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 2.15878; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 1.27445;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 2.92701; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 10.9572;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 2.42273; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 1.46147;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 4.31402; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 15.5701;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 2.99756; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 1.87656;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 5.97014; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 20.987;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 3.63805; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 2.34956;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 7.89535; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 27.2082;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 4.34782; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 2.88332;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 10.0897; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 34.2335;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 5.13004; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 3.47983;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 12.5531; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 42.0629;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 5.98731; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 4.14051;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.370361;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.770271; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.284628;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.108031; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.897359;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.906149; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.386385;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.309168; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 1.63425;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 1.07237; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.500726;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.591572; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 2.58103;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 1.26091; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.625423;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 0.955245; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 3.73771;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 1.46753; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.76032;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 1.40019; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 5.10427;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 1.69025; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.905926;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 1.9264; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 6.68073;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 1.92827; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 1.06294;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 2.53387; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 8.46708;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 2.18149; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 1.23209;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 3.22262; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 10.4633;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 2.45013; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 1.41403;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 4.84391; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 15.0855;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 3.03535; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 1.81855;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 6.79028; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 20.5472;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 3.68774; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 2.28028;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 9.06172; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 26.8485;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 4.41105; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 2.8019;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 11.6582; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 33.9894;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 5.20854; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 3.38532;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 14.5798; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 41.9698;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 6.08286; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 4.03184;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.0220454; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.162401;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.808156; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.283086;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.171452; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.601552;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.955936; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.383883;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.422707; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 1.26289;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 1.13608; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.497343;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.775813; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 2.1464;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 1.34; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.621251;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 1.23077; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 3.2521;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 1.5633; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.755456;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 1.78757; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 4.57998;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 1.80395; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.90047;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 2.44623; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 6.13004;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 2.06123; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 1.057;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 3.20673; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 7.90228;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 2.33511; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 1.22575;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 4.06909; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 9.89671;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 2.62589; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 1.40739;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 6.09934; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 14.5521;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 3.26013; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 1.81159;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 8.537; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 20.0962;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 3.96826; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 2.27333;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 11.3821; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 26.5291;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 4.75443; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 2.79525;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 14.6345; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 33.8507;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 5.62217; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 3.37922;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 18.2944; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 42.061;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 6.57432; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 4.02654;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.0178111; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.108281;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.807727; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.260246;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.129007; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.494062;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.948342; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.346514;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.364076; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 1.08134;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 1.12096; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.445024;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.723016; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 1.87013;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 1.31739; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.553483;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 1.20583; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 2.86041;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 1.53332; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.671562;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 1.81251; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 4.0522;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 1.7667; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.799611;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 2.54307; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 5.44548;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 2.01677; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 0.938192;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 3.3975; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 7.04027;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 2.28344; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 1.0879;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 4.3758; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 8.83656;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 2.56699; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 1.2493;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 6.70401; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 13.0336;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 3.18651; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 1.60909;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 9.52771; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 18.0367;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 3.8794; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 2.02075;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 12.8469; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 23.8458;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 4.64961; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 2.48656;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 16.6616; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 30.4609;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 5.50049; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 3.00811;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 20.9717; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 37.882;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 6.43477; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 3.5865;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.439386;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.786858; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.284371;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.960221;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.921268; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.384305;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.10543; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 1.68345;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 1.08661; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.497278;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.328534; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 2.60907;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 1.275; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.621042;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 0.654852; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 3.73709;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 1.48225; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.755434;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 1.08438; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 5.0675;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 1.70634; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.900968;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 1.61713; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 6.6003;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 1.94648; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 1.05835;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 2.25308; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 8.3355;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 2.20257; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 1.22829;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 2.99225; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 10.2731;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 2.47483; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 1.41146;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 4.78022; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 14.7554;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 3.06954; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 1.8197;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 6.98105; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 20.0474;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 3.73444; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 2.28676;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 9.59473; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 26.1489;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 4.47328; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 2.81525;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 12.6213; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 33.06;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 5.28929; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 3.40699;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 16.0606; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 40.7806;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 6.18506; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 4.06322;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.449624;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.776224; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.282766;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.0361658; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.983633;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.906395; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.382928;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.188963; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 1.72104;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 1.06668; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.495664;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.419208; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 2.66185;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 1.24938; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.61871;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 0.726904; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 3.80606;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 1.45027; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.751877;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 1.11205; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 5.15366;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 1.66734; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.895652;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 1.57464; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 6.70467;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 1.89973; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 1.05072;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 2.11468; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 8.45907;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 2.14727; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 1.21777;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 2.73218; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 10.4169;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 2.41014; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 1.39748;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 4.19951; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 14.9427;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 2.98341; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 1.79704;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 5.97664; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 20.2821;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 3.6231; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 2.25313;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 8.06356; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 26.4351;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 4.33278; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 2.76839;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 10.4603; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 33.4017;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 5.11559; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 3.34469;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 13.1668; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 41.1819;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 5.97407; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 3.98334;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.260195; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.345142;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.819563; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.281398;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.432575; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.78102;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.96496; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.379148;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.706369; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 1.41503;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 1.1429; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.489669;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 1.08158; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 2.24716;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 1.34489; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.61064;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 1.5582; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 3.27743;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 1.56648; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.741826;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 2.13624; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 4.50582;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 1.80561; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.88368;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 2.81569; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 5.93234;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 2.0615; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 1.03686;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 3.59656; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 7.55699;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 2.33406; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 1.20206;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 4.47885; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 9.37977;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 2.62357; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 1.37992;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 6.54766; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 13.6197;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 3.25532; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 1.77576;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 9.02213; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 18.6522;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 3.96095; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 2.22801;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 11.9023; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 24.4771;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 4.74452; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 2.73925;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 15.188; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 31.0946;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 5.60951; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 3.3113;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 18.8795; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 38.5046;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 6.55872; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 3.94542;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.138944; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.322005;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.799485; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.26661;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.347022; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.729992;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.946027; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.354196;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.645037; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 1.32032;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 1.12455; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.454264;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 1.03299; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 2.09299;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 1.32651; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.564411;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 1.51087; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 3.04799;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 1.54754; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.684258;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 2.0787; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 4.18534;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 1.78563; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.814131;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 2.73645; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 5.50502;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 2.04006; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 0.954586;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 3.48415; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 7.00704;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 2.31077; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 1.10622;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 4.32178; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 8.69141;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 2.59808; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 1.26959;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 6.26684; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 12.6072;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 3.2244; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 1.63351;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 8.57165; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 17.2523;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 3.92326; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 2.04959;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 11.2362; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 22.6267;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 4.69879; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 2.52016;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 14.2605; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 28.7305;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 5.55446; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 3.04685;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 17.6445; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 35.5637;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 6.49311; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 3.6308;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.126215;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.762748; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.254293;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.398625;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.892908; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.332866;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.0764599; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.844914;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 1.05263; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.423874;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.235665; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 1.46508;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 1.23411; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.524932;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.449948; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 2.25912;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 1.43315; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.635576;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 0.719308; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 3.22704;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 1.6477; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.756046;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 1.04375; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 4.36884;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 1.87691; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.886821;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 1.42326; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 5.68452;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 2.12059; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 1.02843;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 1.85786; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 7.17407;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 2.3789; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 1.18137;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 2.89228; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 10.6748;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 2.94092; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 1.52302;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 4.14701; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 14.8711;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 3.56643; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 1.91466;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 5.62205; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 19.7628;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 4.259; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 2.35835;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 7.3174; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 25.3501;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 5.02171; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 2.85554;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 9.23306; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 31.6329;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 5.85716; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 3.4072;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.381489;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 0.782232; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.268805;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.844611;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 0.909475; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.359771;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 1.47602;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 1.06672; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.462864;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.027588; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 2.27573;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 1.24639; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.575706;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.237139; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 3.24373;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 1.44429; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.697952;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.5308; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 4.38002;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 1.65835; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.829949;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 0.908572; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 5.6846;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 1.88768; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.972271;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 1.37045; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 7.15747;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 2.13207; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 1.12553;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 1.91645; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 8.79863;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 2.39166; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 1.29031;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 3.26076; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 12.5858;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 2.95785; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 1.65643;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 4.94152; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 17.0462;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 3.58967; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 2.07401;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 6.95873; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 22.1797;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 4.29058; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 2.54549;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 9.31237; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 27.9864;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 5.06365; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 3.07261;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 12.0025; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 34.4663;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 5.91138; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 3.65659;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.160099;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.80558; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.255909;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.149023; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.4942;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.944929; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.338791;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.393132; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.963519;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 1.11564; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.433194;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.69571; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 1.56805;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 1.3094; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.53665;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 1.05676; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 2.30781;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 1.52175; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.648666;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 1.47628; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 3.18278;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 1.75052; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.769457;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 1.95426; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 4.19297;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 1.99485; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.899484;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 2.49072; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 5.33837;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 2.25456; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 1.03927;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 3.08565; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 6.619;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 2.52984; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 1.18933;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 4.45092; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 9.5859;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 3.12871; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 1.52209;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 6.05006; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 13.0937;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 3.79523; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 1.90082;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 7.88308; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 17.1423;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 4.53321; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 2.3278;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 9.94999; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 21.7318;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 5.34596; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 2.80466;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 12.2508; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 26.8622;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 6.23625; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 3.33257;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 0.583514; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.17548;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 0.822491; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.266102;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 0.854092; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.4791;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 0.966345; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.354206;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 1.17505; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.94802;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 1.14232; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.454625;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 1.54639; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 1.58224;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 1.34185; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.564966;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 1.96812; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 2.38176;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 1.56036; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.684855;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 2.44022; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 3.34659;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 1.79566; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.814619;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 2.96271; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 4.47671;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 2.04686; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.954813;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 3.53558; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 5.77214;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 2.3138; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 1.10603;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 4.15883; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 7.23286;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 2.59669; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 1.26885;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 5.55648; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 10.6502;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 3.21199; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 1.63121;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 7.15566; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 14.7288;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 3.89666; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 2.04516;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 8.95636; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 19.4685;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 4.65464; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 2.51306;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 10.9586; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 24.8695;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 5.48936; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 3.03655;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 13.1624; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 30.9317;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 6.40365; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 3.61681;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.0968013; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.375232;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.781797; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.27234;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.215928; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.771276;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.910747; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.364762;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.382627; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 1.31041;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 1.06956; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.469241;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.596899; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 1.99262;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 1.25045; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.583334;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 0.858744; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 2.81792;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 1.44912; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.706661;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 1.16816; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 3.78631;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 1.66343; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.839553;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 1.52515; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 4.89778;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 1.89245; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.982572;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 1.92971; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 6.15233;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 2.13592; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 1.13633;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 2.38185; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 7.54997;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 2.39398; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 1.3014;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 3.42884; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 10.7745;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 2.95519; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 1.66753;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 4.66612; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 14.5714;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 3.57941; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 2.0844;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 6.09369; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 18.9406;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 4.27009; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 2.5545;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 7.71154; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 23.8822;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 5.03032; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 3.07964;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 9.51969; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 29.396;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 5.86265; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 3.66108;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.221396;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.772081; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.258091;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.535139;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.894938; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.342056;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 1.00122;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 1.04715; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.4379;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 1.61965;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 1.2214; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.54321;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.0131595; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 2.39041;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 1.41358; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.65755;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.255909; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 3.31352;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 1.62162; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.781184;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 0.579475; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 4.38897;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 1.84462; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.914613;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 0.983858; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 5.61676;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 2.08236; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 1.05839;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 1.46906; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 6.99689;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 2.33494; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 1.21304;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 2.68191; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 10.2142;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 2.88596; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 1.55684;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 4.21802; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 14.0408;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 3.50092; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 1.94913;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 6.07741; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 18.4769;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 4.18316; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 2.39216;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 8.26006; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 23.5222;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 4.93563; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 2.88755;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 10.766; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 29.177;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 5.76076; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 3.43643;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.014439; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.209847;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.793943; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.255077;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.0737139; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.529808;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.926894; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.336365;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.20827; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.995554;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 1.09056; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.429575;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.418108; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 1.60709;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 1.27703; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.532286;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 0.703228; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 2.3644;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 1.482; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.644031;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 1.06363; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 3.26751;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 1.70334; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.765047;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 1.49931; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 4.31639;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 1.9402; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.895811;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 2.01027; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 5.51107;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 2.19237; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 1.03686;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 2.59652; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 6.85152;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 2.46003; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 1.1887;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 3.99485; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 9.9698;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 3.04334; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 1.52656;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 5.69432; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 13.6712;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 3.69372; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 1.91241;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 7.6949; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 17.9558;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 4.41481; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 2.34845;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 9.99662; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 22.8235;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 5.20981; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 2.83622;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 12.5995; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 28.2743;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 6.08133; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 3.3768;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.164187; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.193639;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.795372; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.265875;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.324409; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.531532;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.930034; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.355042;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 0.541911; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 1.00772;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 1.0955; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.45576;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 0.816692; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 1.6222;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 1.28373; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.565533;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 1.14875; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 2.37497;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 1.49037; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.683902;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 1.53809; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 3.26603;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 1.71328; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.811125;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 1.98471; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 4.29539;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 1.95159; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.947702;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 2.48861; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 5.46304;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 2.20511; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 1.0942;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 3.04979; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 6.76899;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 2.474; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 1.25115;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 4.34398; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 9.79575;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 3.05944; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 1.59838;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 5.86729; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 13.3757;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 3.71154; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 1.9927;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 7.61971; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 17.5088;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 4.43398; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 2.43655;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 9.60126; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 22.1951;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 5.22997; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 2.93171;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 11.8119; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 27.4346;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 6.10218; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 3.47944;
    */


    // RedMET bins
    /*
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.221054; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.00419153;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.324234; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.0257435;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.336798; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.0214484;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.370325; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.0358989;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.467111; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.0486465;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.428102; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.0472892;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.611993; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.085786;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.494805; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.0598036;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 0.771444; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.132867;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 0.568781; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.0734901;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 0.945464; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.189889;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 0.649132; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.0884366;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 1.13405; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 0.256852;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 0.735427; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 0.104734;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 1.33721; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 0.333757;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 0.827505; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 0.122466;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 1.55494; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 0.420603;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 0.925363; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 0.141703;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 2.0341; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 0.624119;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 1.13879; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 0.18492;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 2.57154; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 0.8674;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 1.37678; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 0.234745;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 3.16725; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 1.15045;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 1.64052; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 0.291414;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 3.82124; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 1.47326;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 1.93112; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 0.355085;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 4.53351; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 1.83584;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 2.24951; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 0.425862;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.274545; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.0447112;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.324398; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.0389861;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.38731; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.0766664;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.370562; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.0570161;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.513896; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.115327;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.428438; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.0762128;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.654301; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.160694;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.495267; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.0965756;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 0.808527; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.212766;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 0.569399; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.118232;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 0.976572; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.271544;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 0.649939; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.141334;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 1.15844; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 0.337028;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 0.73646; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 0.166029;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 1.35412; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 0.409217;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 0.828804; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 0.192451;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 1.56363; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 0.488112;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 0.92697; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 0.220719;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 2.0241; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 0.66602;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 1.14115; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 0.283186;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 2.53985; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 0.87075;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 1.38007; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 0.354082;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 3.11089; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 1.1023;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 1.64495; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 0.433873;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 3.7372; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 1.36068;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 1.93688; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 0.522888;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 4.41879; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 1.64588;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 2.25679; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 0.621361;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.0198148;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.406405; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.0468929;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.0606406;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.460919; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.0710185;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.026816; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.122133;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.529644; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.0963145;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.125867; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.204294;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.609325; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.123077;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 0.257278; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.307121;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 0.697946; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.15159;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 0.421047; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.430615;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 0.794379; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.182112;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 0.617175; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 0.574777;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 0.898053; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 0.214869;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 0.845662; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 0.739606;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 1.00874; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 0.250055;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 1.10651; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 0.925102;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 1.12639; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 0.287835;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 1.72528; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 1.3581;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 1.38298; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 0.371712;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 2.47348; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 1.87376;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 1.66897; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 0.46736;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 3.35112; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 2.47209;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 1.98574; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 0.575371;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 4.3582; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 3.15309;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 2.33459; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 0.696156;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 5.49471; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 3.91676;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 2.71662; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 0.829996;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.0368384; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.0413244;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.405651; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.0422085;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.055663; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.087394;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.459212; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.0638643;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.105851; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.150425;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.526842; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.0865833;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.187402; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.230417;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.60535; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.110624;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 0.300316; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.327371;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 0.692744; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.13624;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 0.444593; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.441286;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 0.787898; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.163663;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 0.620233; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 0.572163;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 0.890237; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 0.193094;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 0.827237; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 0.72;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 0.999522; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 0.224708;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 1.0656; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 0.884799;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 1.11571; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 0.258653;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 1.63643; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 1.26528;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 1.36912; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 0.334017;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 2.3327; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 1.71361;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 1.65159; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 0.419956;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 3.15443; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 2.22978;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 1.96445; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 0.517004;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 4.10161; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 2.8138;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 2.30898; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 0.625528;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 5.17424; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 3.46566;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 2.68626; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 0.745784;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.112816;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.322547; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.0301451;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.183336; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.00209901;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.367605; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.0431161;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.272074; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.0309295;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.424185; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.0573777;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.379029; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.0737962;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.489591; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.0728778;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 0.504203; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.130699;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 0.56219; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.0897113;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 0.647594; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.201639;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 0.641089; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.108003;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 0.809203; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 0.286614;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 0.72585; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 0.127874;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 0.989031; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 0.385626;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 0.816307; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 0.149431;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 1.18708; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 0.498674;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 0.912447; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 0.172766;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 1.63782; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 0.766879;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 1.12213; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 0.225062;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 2.16143; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 1.09123;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 1.3559; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 0.285222;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 2.75792; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 1.47172;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 1.61494; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 0.353552;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 3.42728; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 1.90836;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 1.9003; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 0.430258;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 4.16951; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 2.40115;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 2.21292; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 0.515476;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.244467;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.324902; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.0243833;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.34385; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.00113861;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.371401; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.0336171;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.459948; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.0193211;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.429656; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.0440529;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.592763; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.0466003;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.496886; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.0555602;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 0.742294; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.0829763;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 0.57143; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.0681695;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 0.908541; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.128449;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 0.652393; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.0819547;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 1.0915; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 0.183019;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 0.739345; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 0.096996;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 1.29118; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 0.246685;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 0.832132; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 0.113368;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 1.50758; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 0.319448;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 0.930753; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 0.131134;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 1.99052; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 0.492265;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 1.14589; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 0.171058;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 2.54032; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 0.701469;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 1.38583; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 0.217094;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 3.15699; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 0.94706;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 1.6518; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 0.26946;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 3.84052; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 1.22904;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 1.94492; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 0.328299;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 4.59092; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 1.5474;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 2.2661; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 0.393708;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.211687; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.0776597;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.407078; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.043325;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.30164; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.139619;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.461383; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.0654981;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.407643; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.217002;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.52984; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.0886772;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.529697; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.30981;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.60919; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.113108;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 0.667802; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.418041;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 0.697403; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.139036;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 0.821958; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.541697;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 0.793338; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.166683;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 0.992165; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 0.680777;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 0.89641; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 0.196248;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 1.17842; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 0.835282;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 1.00637; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 0.227904;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 1.38073; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 1.00521;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 1.12319; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 0.2618;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 1.8335; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 1.39134;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 1.37768; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 0.336794;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 2.35048; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 1.83917;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 1.66101; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 0.422022;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 2.93166; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 2.34869;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 1.97453; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 0.51804;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 3.57704; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 2.91991;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 2.31952; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 0.625237;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 4.28663; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 3.55283;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 2.6971; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 0.743888;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.0483471;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.406726; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.0295539;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.105582; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.0220613;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.461224; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.0447234;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.187801; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.0760283;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.529904; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.0609428;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.295003; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.15328;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.609499; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.0784306;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 0.427188; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.253818;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 0.697981; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.0974039;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 0.584356; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.37764;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 0.794209; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.118055;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 0.766508; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 0.524748;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 0.897607; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 0.140544;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 0.973643; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 0.695141;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 1.00793; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 0.165007;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 1.20576; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 0.88882;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 1.12515; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 0.191551;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 1.74495; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 1.34603;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 1.38058; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 0.251218;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 2.38407; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 1.89639;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 1.66503; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 0.32006;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 3.12312; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 2.53988;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 1.97987; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 0.398407;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 3.96211; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 3.27652;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 2.32639; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 0.486472;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 4.90103; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 4.10629;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 2.7057; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 0.584398;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.229934; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.0164997;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.325603; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.0319965;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.383057; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.0323673;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.372877; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.0459351;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.538455; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.0538823;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.431941; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.0609829;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.696128; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.0810448;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.499935; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.0770618;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.856077; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.113855;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.575155; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.0942453;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 1.0183; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.152312;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 0.656686; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.112642;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 1.1828; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.196417;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 0.744086; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.132365;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 1.34957; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 0.24617;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 0.837191; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 0.153518;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 1.51862; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 0.30157;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 0.935997; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 0.176193;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 1.86355; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 0.429312;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 1.15109; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 0.226419;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 2.21758; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 0.579644;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 1.39043; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 0.283549;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 2.58071; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 0.752567;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 1.65526; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 0.347945;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 2.95294; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 0.948079;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 1.94668; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 0.419863;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 3.33427; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 1.16618;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 2.26567; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 0.49948;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 0.256506; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.0501682;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 0.321544; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.0352407;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 0.346951; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.0826095;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 0.365293; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.0510525;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 0.441454; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.120022;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 0.420328; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.0679359;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.540015; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.162405;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.484014; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.0858335;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.642633; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.20976;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.55473; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.104831;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.749309; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.262085;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.631569; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.125048;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 0.860043; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.319382;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 0.714075; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.146604;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 0.974834; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 0.381649;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 0.802061; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 0.169613;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 1.09368; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 0.448888;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 0.895498; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 0.194175;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 1.34355; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 0.598278;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 1.09901; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 0.248301;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 1.60965; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 0.767553;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 1.32552; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 0.309553;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 1.89198; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 0.956711;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 1.57615; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 0.378345;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 2.19054; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 1.16575;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 1.85191; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 0.454977;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 2.50533; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 1.39468;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 2.15371; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 0.539661;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.0946488; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.0334306;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.406109; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.0312987;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.145201; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.0681808;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.459909; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.0471851;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.211776; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.114928;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.52775; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.063912;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.294374; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.173671;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.60639; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.0816536;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.392995; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.244411;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.693805; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.100594;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 0.50764; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.327148;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 0.788845; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.120902;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 0.638307; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.421881;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 0.89092; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.142729;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 0.784997; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 0.528611;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 0.999776; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 0.166202;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 0.947711; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 0.647337;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 1.11536; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 0.191431;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 1.32121; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 0.92078;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 1.36702; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 0.247512;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 1.7588; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 1.24221;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 1.64696; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 0.311537;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 2.26048; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 1.61162;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 1.95652; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 0.383896;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 2.82625; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 2.02903;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 2.29697; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 0.464856;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 3.45611; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 2.49441;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 2.6694; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 0.5546;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.0242669;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 0.402606; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.0392019;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.0620159;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 0.454221; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.0592738;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.116862;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 0.519517; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.0803463;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.0470483; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.188804;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.595402; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.102656;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.116051; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.277844;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.67991; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.126436;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.205505; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.38398;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.771906; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.151903;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 0.31541; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.507213;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 0.870796; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.179242;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 0.445767; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 0.647542;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 0.976313; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 0.208617;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 0.596575; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 0.804969;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 1.08839; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 0.240165;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 0.959546; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 1.17111;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 1.33248; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 0.310223;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 1.40432; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 1.60564;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 1.60404; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 0.390133;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 1.9309; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 2.10856;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 1.90432; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 0.480389;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 2.53929; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 2.67986;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 2.23454; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 0.581329;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 3.22949; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 3.31956;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 2.59574; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 0.69319;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.222953; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.0551377;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.320424; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.0325754;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.286911; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.0879254;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.363272; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.0467907;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.358888; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.124002;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.417286; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.0620339;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.438884; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.163369;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.479898; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.0782114;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 0.526901; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.206025;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 0.549515; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.0953825;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 0.622937; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.25197;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 0.62524; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.113645;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 0.726992; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.301205;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 0.706618; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.133104;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 0.839067; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 0.353729;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 0.79346; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 0.153858;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 0.959162; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 0.409542;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 0.885733; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 0.175996;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 1.22341; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 0.531037;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 1.08684; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 0.22473;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 1.51974; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 0.665689;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 1.31083; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 0.279818;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 1.84814; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 0.813499;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 1.55878; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 0.341638;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 2.20863; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 0.974467;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 1.8317; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 0.410462;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 2.60119; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 1.14859;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 2.13046; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 0.486486;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.315364; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.0044267;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.324297; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.024719;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.434351; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.0165361;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.370073; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.034081;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 0.560883; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.0350455;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 0.427472; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.0445292;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 0.694959; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.0599549;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 0.493741; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.0559077;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.83658; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.0912643;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.567223; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.0682284;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.985744; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.128974;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.647012; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.0815521;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 1.14245; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.173083;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 0.732668; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.0959513;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 1.30671; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 0.223593;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 0.824022; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 0.111496;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 1.4785; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 0.280502;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 0.921065; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 0.128249;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 1.84473; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 0.413521;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 1.13257; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 0.165592;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 2.24113; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 0.572139;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 1.3682; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 0.208323;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 2.66771; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 0.756358;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 1.62915; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 0.256683;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 3.12447; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 0.966177;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 1.9165; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 0.310839;
    extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 3.6114; extrvalpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 1.2016;
    //extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 2.23119; extrerrpoints["7TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 0.370904;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.0518076;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.402986; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.0400233;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.000114047; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.0946174;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.454709; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.0604533;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.0252127; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.148741;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.520155; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.0818147;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.0724259; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.214179;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.596232; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.104327;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 0.141754; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.290931;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 0.680978; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.128214;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 0.233196; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.378998;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 0.773265; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.153677;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 0.346753; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.478378;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 0.872502; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.180901;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 0.482424; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 0.589072;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 0.978427; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 0.210043;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 0.64021; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 0.711081;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 1.09098; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 0.24124;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 1.02213; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 0.98904;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 1.33622; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 0.310244;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 1.4925; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 1.31226;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 1.60922; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 0.388643;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 2.05133; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 1.68073;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 1.91123; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 0.476951;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 2.69862; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 2.09446;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 2.24349; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 0.575528;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 3.43437; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 2.55344;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 2.60705; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 0.684626;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.0562472; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.0604374;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.405418; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.0378468;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.0879915; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.105112;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.458537; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.0570986;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 0.136598; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.159445;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 0.525596; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.0772022;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 0.202066; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.223438;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 0.603401; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.0983518;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.284395; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.297089;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.689945; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.120748;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.383587; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.3804;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.784084; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.144578;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 0.49964; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.47337;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 0.885226; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.170008;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 0.632555; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 0.575999;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 0.993113; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 0.197185;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 0.782332; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 0.688287;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 1.10769; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 0.226235;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 1.13247; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 0.941839;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 1.35719; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 0.290376;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 1.55006; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 1.23403;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 1.63477; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 0.363115;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 2.03509; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 1.56485;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 1.94173; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 0.444942;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 2.58757; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 1.93432;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 2.27934; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 0.536201;
    extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 3.2075; extrvalpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 2.34241;
    //extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 2.64867; extrerrpoints["7TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 0.637136;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.0690766; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.505606;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 1.09078; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.296728;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.198031; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 1.08287;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 1.26971; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.403411;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.416818; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 1.85828;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 1.49027; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.522853;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 0.725438; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 2.83184;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 1.74171; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.652701;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 1.12389; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 4.00354;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 2.01805; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.792759;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 1.61217; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 5.37338;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 2.3163; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.943533;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 2.19029; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 6.94137;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 2.63518; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 1.10573;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 2.85824; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 8.7075;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 2.97435; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 1.28009;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 3.61602; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 10.6718;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 3.33397; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 1.46729;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 5.40107; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 15.1947;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 4.11646; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 1.88259;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 7.54546; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 20.5103;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 4.98728; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 2.35558;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 10.0492; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 26.6184;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 5.95125; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 2.88912;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 12.9122; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 33.5191;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 7.01266; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 3.48523;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 16.1346; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 41.2124;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 8.17505; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 4.14535;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.0570039; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.33299;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 1.08884; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.289027;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.242376; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.845793;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 1.27127; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.391454;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.523834; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 1.56756;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 1.49559; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.506666;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 0.901378; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 2.49829;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 1.75086; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.632339;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 1.37501; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 3.63799;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 2.03107; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.768272;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 1.94472; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 4.98664;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 2.33328; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.914953;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 2.61053; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 6.54427;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 2.65623; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 1.07308;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 3.37242; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 8.31085;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 2.99962; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 1.24336;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 4.23039; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 10.2864;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 3.36366; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 1.42646;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 6.2346; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 14.8644;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 4.15565; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 1.8334;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 8.62315; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 20.2782;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 5.03702; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 2.29769;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 11.396; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 26.5279;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 6.0127; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 2.82207;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 14.5533; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 33.6135;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 7.08705; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 3.40844;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 18.0949; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 41.5349;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 8.26368; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 4.05815;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 0.185863; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.276236;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][6] = 1.21654; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.002"][7] = 0.208007;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 0.381045; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.742693;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][6] = 1.42057; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.003"][7] = 0.296227;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 0.696814; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 1.42276;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][6] = 1.67148; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.004"][7] = 0.393353;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 1.13317; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 2.31645;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][6] = 1.95703; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.005"][7] = 0.498895;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 1.69011; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 3.42375;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][6] = 2.27055; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.006"][7] = 0.613434;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 2.36764; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 4.74467;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][6] = 2.60874; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.007"][7] = 0.737785;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 3.16576; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 6.2792;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][6] = 2.97023; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.008"][7] = 0.872755;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 4.08446; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 8.02734;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][6] = 3.35468; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.009"][7] = 1.01907;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 5.12375; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 9.9891;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][6] = 3.76233; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.01"][7] = 1.17734;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 7.56409; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 14.5535;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][6] = 4.64948; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.012"][7] = 1.53175;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 10.4868; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 19.9723;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][6] = 5.63709; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.014"][7] = 1.93913;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 13.8918; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 26.2456;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][6] = 6.73069; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.016"][7] = 2.40161;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 17.7792; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 33.3733;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][6] = 7.93517; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.018"][7] = 2.9206;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 22.1489; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 41.3555;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][6] = 9.25455; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=-0.02"][7] = 3.49705;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 0.274531; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.180393;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][6] = 1.20549; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.002"][7] = 0.184575;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 0.512531; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.569346;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][6] = 1.40425; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.003"][7] = 0.258332;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 0.870813; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 1.15162;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][6] = 1.64923; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.004"][7] = 0.340778;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 1.34938; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 1.92722;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][6] = 1.92853; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.005"][7] = 0.431173;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 1.94822; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 2.89615;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][6] = 2.23558; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.006"][7] = 0.52989;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 2.66735; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 4.05839;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][6] = 2.56711; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.007"][7] = 0.637573;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 3.50675; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 5.41396;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][6] = 2.92174; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.008"][7] = 0.754888;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 4.46644; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 6.96285;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][6] = 3.29912; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.009"][7] = 0.882436;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 5.54641; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 8.70507;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][6] = 3.69946; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.01"][7] = 1.02074;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 8.0672; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 12.7695;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][6] = 4.57119; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.012"][7] = 1.33125;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 11.0691; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 17.6072;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][6] = 5.54216; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.014"][7] = 1.68903;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 14.5521; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 23.2182;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][6] = 6.61778; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.016"][7] = 2.09581;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 18.5163; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 29.6024;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][6] = 7.80278; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.018"][7] = 2.55275;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 22.9616; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 36.76;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][6] = 9.10112; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4z=0.02"][7] = 3.06062;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.180537; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.432102;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 1.10236; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.28864;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.366099; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.947076;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 1.28697; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.38896;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 0.639016; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 1.66215;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 1.51398; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.502528;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 0.999287; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 2.57732;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 1.77228; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.626998;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 1.44691; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 3.69258;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 2.05583; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.762154;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 1.98189; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 5.00794;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 2.36161; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.908483;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 2.60423; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 6.5234;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 2.68837; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 1.06668;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 3.31392; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 8.23896;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 3.03578; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 1.23745;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 4.11096; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 10.1546;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 3.40406; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 1.42146;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 5.96711; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 14.5862;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 4.20522; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 1.83144;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 8.17268; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 19.8182;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 5.09669; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 2.30031;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 10.7277; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 25.8506;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 6.08347; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 2.83072;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 13.6321; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 32.6833;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 7.16997; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 3.4245;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 16.8859; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 40.3165;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 8.35985; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 4.08291;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.044853; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.487716;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 1.09557; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.292044;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.272082; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 1.02138;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 1.27622; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.395591;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 0.582308; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 1.75049;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 1.49881; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.511923;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 0.975531; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 2.67503;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 1.75248; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.63867;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 1.45175; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 3.79501;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 2.03125; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.775608;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 2.01097; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 5.11044;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 2.33212; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.923215;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 2.65318; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 6.6213;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 2.65382; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 1.08218;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 3.3784; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 8.3276;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 2.99603; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 1.25322;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 4.1866; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 10.2293;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 3.35892; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 1.43699;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 6.05201; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 14.6191;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 4.1487; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 1.84506;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 8.24941; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 19.7907;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 5.02786; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 2.31021;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 10.7788; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 25.744;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 6.0013; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 2.83522;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 13.6402; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 32.4791;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 7.07333; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 3.42205;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 16.8335; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 39.9959;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 8.24754; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 4.07207;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 0.554305; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.407637;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][6] = 1.22201; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.002"][7] = 0.219581;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 0.847244; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.879359;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][6] = 1.42536; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.003"][7] = 0.314406;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 1.25145; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 1.53803;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][6] = 1.6757; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.004"][7] = 0.417969;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 1.76692; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 2.38364;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][6] = 1.96086; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.005"][7] = 0.529823;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 2.39366; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 3.4162;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][6] = 2.27414; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.006"][7] = 0.650592;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 3.13166; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 4.6357;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][6] = 2.61224; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.007"][7] = 0.781127;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 3.98093; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 6.04215;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][6] = 2.97378; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.008"][7] = 0.922273;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 4.94147; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 7.63554;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][6] = 3.35841; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.009"][7] = 1.07479;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 6.01327; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 9.41588;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][6] = 3.76635; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.01"][7] = 1.23934;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 8.49068; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 13.5374;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][6] = 4.65445; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.012"][7] = 1.60667;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 11.4132; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 18.4067;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][6] = 5.64344; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.014"][7] = 2.02767;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 14.7807; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 24.0238;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][6] = 6.73885; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.016"][7] = 2.5047;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 18.5933; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 30.3886;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][6] = 7.94554; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.018"][7] = 3.03933;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 22.851; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 37.5013;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][6] = 9.26754; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=-0.02"][7] = 3.63264;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 0.42946; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.361958;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][6] = 1.21381; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.002"][7] = 0.202292;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 0.85057; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.79567;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][6] = 1.42059; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.003"][7] = 0.286526;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 1.35682; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 1.40533;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][6] = 1.67432; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.004"][7] = 0.379397;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 1.9482; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 2.19095;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][6] = 1.96259; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.005"][7] = 0.480268;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 2.62472; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 3.15252;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][6] = 2.27869; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.006"][7] = 0.589612;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 3.38638; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 4.29004;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][6] = 2.61933; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.007"][7] = 0.708165;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 4.23317; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 5.60351;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][6] = 2.98315; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.008"][7] = 0.836674;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 5.1651; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 7.09294;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][6] = 3.36983; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.009"][7] = 0.975818;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 6.18217; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 8.75831;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][6] = 3.77964; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.01"][7] = 1.12618;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 8.47171; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 12.6169;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][6] = 4.67091; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.012"][7] = 1.46247;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 11.1018; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 17.1793;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][6] = 5.66245; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.014"][7] = 1.84858;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 14.0724; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 22.4456;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][6] = 6.75987; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.016"][7] = 2.28656;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 17.3836; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 28.4156;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][6] = 7.9681; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.018"][7] = 2.7778;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 21.0354; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 35.0895;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][6] = 9.29124; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5z=0.02"][7] = 3.32323;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.128316;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 1.08642; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.258389;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.395344;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 1.26426; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.337096;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 0.0532857; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.831824;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 1.48338; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.428372;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 0.25098; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 1.43776;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 1.73301; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.529738;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 0.522583; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 2.21314;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 2.00716; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.640672;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 0.868096; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 3.15798;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 2.30281; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.761378;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 1.28752; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 4.27226;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 2.61864; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.892316;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 1.78085; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 5.556;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 2.95427; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 1.034;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 2.34809; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 7.0092;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 3.30985; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 1.18694;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 3.7043; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 10.4239;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 4.08266; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 1.52831;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 5.35614; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 14.5165;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 4.94155; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 1.91934;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 7.30362; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 19.2868;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 5.89128; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 2.3621;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 9.54674; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 24.735;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 6.93607; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 2.85805;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 12.0855; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 30.861;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 8.07949; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 3.40822;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.36424;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 1.09145; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.273203;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.823103;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 1.26655; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.364731;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 1.44898;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 1.48295; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.468567;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 0.113163; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 2.24186;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 1.73006; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.582238;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 0.383686; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 3.20175;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 2.00196; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.705348;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 0.742148; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 4.32866;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 2.29563; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.838221;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 1.18855; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 5.62257;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 2.60975; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.981415;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 1.7229; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 7.08349;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 2.94393; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 1.13554;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 2.34518; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 8.71143;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 3.29829; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 1.30117;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 3.85357; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 12.4683;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 4.06937; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 1.66899;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 5.71372; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 16.8933;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 4.92741; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 2.08825;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 7.92563; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 21.9862;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 5.87711; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 2.56145;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 10.4893; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 27.7472;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 6.92264; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 3.09033;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 13.4047; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 34.1763;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 8.0675; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 3.67613;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 0.422381; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.205388;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][6] = 1.21749; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.002"][7] = 0.181532;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 0.788028; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.52409;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][6] = 1.41924; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.003"][7] = 0.252827;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 1.20051; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.977237;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][6] = 1.66734; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.004"][7] = 0.332125;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 1.65983; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 1.56483;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][6] = 1.94951; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.005"][7] = 0.418499;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 2.16598; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 2.28687;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][6] = 2.25897; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.006"][7] = 0.512172;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 2.71896; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 3.14335;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][6] = 2.59232; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.007"][7] = 0.613684;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 3.31878; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 4.13428;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][6] = 2.94808; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.008"][7] = 0.723624;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 3.96544; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 5.25965;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][6] = 3.32581; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.009"][7] = 0.842547;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 4.65893; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 6.51947;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][6] = 3.72569; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.01"][7] = 0.970944;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 6.18641; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 9.44244;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][6] = 4.59389; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.012"][7] = 1.25776;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 7.90123; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 12.9032;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][6] = 5.55775; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.014"][7] = 1.58668;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 9.80339; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 16.9017;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][6] = 6.62263; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.016"][7] = 1.95947;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 11.8929; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 21.438;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][6] = 7.79332; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.018"][7] = 2.37735;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 14.1697; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 26.5121;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][6] = 9.07385; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=-0.02"][7] = 2.84115;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 1.42315; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.241305;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][6] = 1.23483; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.002"][7] = 0.202281;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 1.76332; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.58287;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][6] = 1.43756; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.003"][7] = 0.286506;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 2.15836; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 1.08039;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][6] = 1.68718; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.004"][7] = 0.379178;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 2.60828; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 1.73387;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][6] = 1.97136; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.005"][7] = 0.479626;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 3.11307; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 2.54331;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][6] = 2.28327; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.006"][7] = 0.588299;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 3.67272; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 3.50871;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][6] = 2.61945; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.007"][7] = 0.705911;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 4.28725; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 4.63006;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][6] = 2.97838; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.008"][7] = 0.833195;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 4.95665; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 5.90738;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][6] = 3.35963; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.009"][7] = 0.970822;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 5.68092; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 7.34065;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][6] = 3.76335; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.01"][7] = 1.11938;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 7.29406; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 10.6751;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][6] = 4.64023; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.012"][7] = 1.45116;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 9.1267; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 14.6333;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][6] = 5.6141; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.014"][7] = 1.83159;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 11.1788; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 19.2154;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][6] = 6.69034; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.016"][7] = 2.26276;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 13.4504; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 24.4213;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][6] = 7.87378; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.018"][7] = 2.74608;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 15.9415; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 30.251;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][6] = 9.16847; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f4g=0.02"][7] = 3.2825;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.151923; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.380363;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 1.09808; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.277078;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.363752; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.768617;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 1.27983; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.370093;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 0.627371; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 1.29787;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 1.50342; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.475421;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 0.94278; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 1.96813;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 1.75782; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.590528;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 1.30998; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 2.77939;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 2.03694; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.714989;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 1.72897; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 3.73165;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 2.33772; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.849115;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 2.19975; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 4.82492;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 2.65884; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.993461;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 2.72232; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 6.05918;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 2.99992; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 1.14863;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 3.29668; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 7.43445;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 3.36111; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 1.31522;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 4.60077; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 10.608;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 4.1457; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 1.68465;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 6.11202; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 14.3455;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 5.01719; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 2.10522;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 7.83043; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 18.6471;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 5.98043; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 2.57945;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 9.756; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 23.5127;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 7.03975; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 3.10915;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 11.8887; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 28.9422;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 8.19878; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 3.69561;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.000481091; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.22117;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 1.08904; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.262757;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.519389;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 1.26027; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.347034;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.967502;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 1.47241; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.443433;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 0.0739256; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 1.56551;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 1.71513; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.549449;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 0.282866; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 2.31341;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 1.98259; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.664592;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 0.584037; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 3.2112;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 2.27177; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.789101;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 0.977437; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 4.25889;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 2.58132; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.923464;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 1.46307; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 5.45647;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 2.91083; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 1.06823;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 2.04093; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 6.80394;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 3.26042; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 1.22392;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 3.47333; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 9.94857;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 4.02151; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 1.56994;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 5.27466; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 13.6928;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 4.86887; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 1.96467;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 7.4449; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 18.0365;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 5.80707; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 2.41039;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 9.98407; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 22.9799;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 6.84021; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 2.90871;
    extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 12.8921; extrvalpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 28.5228;
    //extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 7.97171; extrerrpoints["8TeV_ee_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 3.46079;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 0.703883; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.263283;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][6] = 1.20778; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.002"][7] = 0.185129;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 0.906554; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.597904;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][6] = 1.4021; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.003"][7] = 0.258681;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 1.18248; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 1.07183;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][6] = 1.6421; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.004"][7] = 0.34053;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 1.53166; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 1.68506;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][6] = 1.91601; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.005"][7] = 0.429829;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 1.95409; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 2.43759;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][6] = 2.21725; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.006"][7] = 0.526871;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 2.44977; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 3.32942;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][6] = 2.54246; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.007"][7] = 0.632248;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 3.0187; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 4.36055;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][6] = 2.89015; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.008"][7] = 0.746591;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 3.66089; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 5.53099;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][6] = 3.2599; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.009"][7] = 0.870484;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 4.37633; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 6.84073;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][6] = 3.65185; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.01"][7] = 1.00444;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 6.02697; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 9.87813;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][6] = 4.50428; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.012"][7] = 1.30418;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 7.97062; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 13.4727;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][6] = 5.45233; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.014"][7] = 1.64848;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 10.2073; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 17.6245;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][6] = 6.50119; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.016"][7] = 2.03914;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 12.737; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 22.3336;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][6] = 7.65549; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.018"][7] = 2.47738;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 15.5596; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 27.5998;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][6] = 8.91912; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=-0.02"][7] = 2.96401;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 0.514953; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.133128;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][6] = 1.20613; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.002"][7] = 0.176796;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 0.901126; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.433211;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][6] = 1.40298; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.003"][7] = 0.245059;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 1.33206; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.876443;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][6] = 1.64553; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.004"][7] = 0.321198;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 1.80775; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 1.46282;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][6] = 1.92182; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.005"][7] = 0.404225;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 2.32821; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 2.19235;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][6] = 2.22516; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.006"][7] = 0.494308;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 2.89342; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 3.06503;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][6] = 2.55219; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.007"][7] = 0.591941;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 3.50339; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 4.08085;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][6] = 2.90142; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.008"][7] = 0.697679;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 4.15812; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 5.23983;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][6] = 3.27241; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.009"][7] = 0.812048;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 4.85762; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 6.54195;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][6] = 3.6653; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.01"][7] = 0.935516;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 6.39088; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 9.57564;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][6] = 4.51872; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.012"][7] = 1.21128;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 8.10319; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 13.1819;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][6] = 5.4666; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.014"][7] = 1.52747;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 9.99454; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 17.3608;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][6] = 6.51415; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.016"][7] = 1.88579;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 12.0649; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 22.1123;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][6] = 7.66608; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.018"][7] = 2.28742;
    extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 14.3144; extrvalpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 27.4363;
    //extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"] = std::vector<double>(8, -1.); extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][6] = 8.9263; extrerrpoints["8TeV_mumu_ZZ#rightarrow 2l2#nu f5g=0.02"][7] = 2.73316;
    */
}
