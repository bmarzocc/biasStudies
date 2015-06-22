#include <cstring>
#include <cerrno>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <unistd.h>
#include <errno.h>
#include <iomanip>
// ROOT headers
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TIterator.h"

#include "TLatex.h"
#include "TString.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TFitResult.h" 

// RooFit headers
#include "RooAbsPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TStyle.h"
// RooStats headers
#include "RooStats/HLFactory.h"
#include "RooStats/RooStatsUtils.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsData.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooProduct.h"
#include "RooExtendPdf.h"
#include "RooBernstein.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooGenericPdf.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooCmdArg.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooChi2MCSModule.h"

#include "HiggsCSandWidth.h"

using namespace RooFit;
using namespace RooStats ;

Int_t NCAT;
TString inDir;
Int_t resMass;
bool withCorr;
RooAbsPdf* BkgMggModelFitNoCorr;
int NDOFMgg_noCorr;
int NDOFMgg;
int NDOFMjj;



void AddBkgData(RooWorkspace*, int);
void AddSigData(RooWorkspace*, int);
void SigModelFit(RooWorkspace*, int);
RooAbsPdf* BkgMggModelFit(RooWorkspace*, int, int);
RooAbsPdf* BkgMjjModelFit(RooWorkspace*, int, int);
void BkgModelBias(RooWorkspace*,int,RooAbsPdf*,RooAbsPdf*,FILE*,FILE*, FILE*);
void SetParamNames(RooWorkspace*);
void SetConstantParams(const RooArgSet* params);
Double_t effSigma(TH1 *hist);
void style();
 
RooArgSet* defineVariables()
{
  // define variables of the input ntuple
  RooRealVar* mJJ  = new RooRealVar("mjj","M(jj)",60,180,"GeV");
  RooRealVar* mGG  = new RooRealVar("mgg","M(#gamma#gamma)",100,180,"GeV");
  if(resMass == 270) mGG->setRange(100,155);
  else if (resMass == 300) mGG->setRange(100,175);
  RooRealVar* mRad = new RooRealVar("mtot","M(#gamma#gamma jj)",0,1500,"GeV");
  RooRealVar* wei  = new RooRealVar("evWeight","event weight",0,100,"");
  RooCategory* bJetTagCategory = new RooCategory("cut_based_ct","event category 4") ;
  bJetTagCategory->defineType("cat0",0);
  bJetTagCategory->defineType("cat1",1);
  bJetTagCategory->defineType("cat2",2);
  bJetTagCategory->defineType("cat3",3);

  RooArgSet* ntplVars = new RooArgSet(*mGG,*mJJ,*mRad,*bJetTagCategory, *wei);
   
 
  return ntplVars;
}


void runfits(int cat=0, int modelNumMgg=0, int modelNumMjj=0, int inDirNum=0)
{

  //create truth models
  RooAbsPdf *MggBkgTruth=0;
  RooAbsPdf *MjjBkgTruth=0;

  style();

  //TString card_name("hgghbb_models_Pol_8TeV.rs");
  TString card_name("/afs/cern.ch/work/b/bmarzocc/BiasStudies/CMSSW_7_4_0/src/biasStudies/models_2D_correlation.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  AddBkgData(w,cat);
  AddSigData(w,cat);
  SigModelFit(w, cat);

  FILE *fout = fopen("resultsBias2D.txt","a");
  FILE *foutTotal = fopen("resultsBias2D_total.txt","a");
  FILE *foutCorr = fopen("resultsCorr2D.txt","a");
  //if(modelNumMgg==0 && modelNumMjj ==0) fprintf(fout,"%s\n\n",inDir.Data());

  if(modelNumMgg==2 || modelNumMjj==2) return;//skip Landau, it sucks.
  if(modelNumMgg==3 || modelNumMjj==3) return;//skip Laurent, it sucks.
  if(withCorr && modelNumMgg==0 ) return; //skip Ber for mgg correlation

  MggBkgTruth = BkgMggModelFit(w,cat,modelNumMgg); //Ber, Exp, Lan, Lau, Pow
  MjjBkgTruth = BkgMjjModelFit(w,cat,modelNumMjj); //Ber, Exp, Lan, Lau, Pow
  BkgModelBias(w,cat,MggBkgTruth,MjjBkgTruth,fout,foutTotal,foutCorr);


  //if(modelNumMgg==4 && modelNumMjj==4) fprintf(fout,"\n\n");
  fclose(fout);
  fclose(foutTotal);
  fclose(foutCorr);
  return;
}


void AddBkgData(RooWorkspace* w, int cat) {

  Int_t ncat = NCAT;

// common preselection cut
  TString mainCut("1");
  //TString mainCut("mRad>200 && mRad<700 && mJJ>90 && mJJ<170");
  //if(cat==0) mainCut = "60<mJJ && mJJ<180";
  //else if(cat==1) mainCut = "85<mJJ && mJJ<155 && mRad>225 && mRad<290";
  //else if(cat==2) mainCut = "85<mJJ && mJJ<155 && mRad>225 && mRad<290";

//****************************//
// Signal Data Set
//****************************//

  // Variables
  RooArgSet* ntplVars = defineVariables();

//****************************//
// CMS Data Set
//****************************//
// retrieve the data tree;
// no common preselection cut applied yet; 

  TFile dataFile(TString::Format("%sDataCS_m%d.root",inDir.Data(),resMass));
  //TFile dataFile(TString::Format("%sData_m%d.root",inDir.Data(),resMass));   
  //TFile dataFile(TString::Format("%ssum_bkg_m%d.root",inDir.Data(),resMass));   
  //TFile dataFile(TString::Format("mcSum_m%d.root",resMass));   
  TTree* dataTree     = (TTree*) dataFile.Get("TCVARS");

  RooDataSet Data("Data","dataset",dataTree,*ntplVars,"","evWeight");

// apply a common preselection cut;
// split into NCAT  categories;

  
  RooDataSet* dataToFit[9];
  for (int c = 0; c < ncat; ++c) {
// Real data
    dataToFit[c]   = (RooDataSet*) Data.reduce(RooArgList(*w->var("mgg"),*w->var("mjj")),mainCut+TString::Format(" && cut_based_ct==%d",c));
    w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
  }

// Create full data set without categorization
  RooDataSet* data    = (RooDataSet*) Data.reduce(RooArgList(*w->var("mgg"),*w->var("mjj")),mainCut);
  w->import(*data, Rename("Data"));
  data->Print("v");

}

void AddSigData(RooWorkspace* w, int cat) {
  const Int_t ncat = NCAT;

  RooArgSet* ntplVars = defineVariables();
  TFile* sigFile;
  if (resMass==0)
    sigFile = new TFile(TString::Format("%sggHH_Lam_1d0_Yt_1d0_c2_0d0_8TeV_m0.root",inDir.Data()));
  else 
    sigFile = new TFile(TString::Format("%sRadion_m%d_8TeV_m%d.root",inDir.Data(),resMass,resMass));
  TTree* sigTree = (TTree*) sigFile->Get("TCVARS");
  // common preselection cut
  TString mainCut("1");
  RooDataSet sigScaled(
		       "sigScaled",
		       "dataset",
		       sigTree,
		       *ntplVars,
		       mainCut,
		       "evWeight");

  RooDataSet* sigToFit[ncat];
  TString cut0 = " && 1>0";
  //
  // we take only mtot to fit to the workspace, we include the cuts
  for ( int i=0; i<ncat; ++i){
    sigToFit[i] = (RooDataSet*) sigScaled.reduce(
						 RooArgList(*w->var("mgg"),*w->var("mjj")),
						 mainCut+TString::Format(" && cut_based_ct==%d ",i)+cut0);
    w->import(*sigToFit[i],Rename(TString::Format("Sig_cat%d",i)));
  }
  // Create full signal data set without categorization
  RooDataSet* sigToFitAll = (RooDataSet*) sigScaled.reduce(
							   RooArgList(*w->var("mgg"),*w->var("mjj")),
							   mainCut);

  w->import(*sigToFitAll,Rename("Sig"));
  sigFile->Close();
  return;
} // end add signal function

void SigModelFit(RooWorkspace* w, int cat) {
  const Int_t ncat = NCAT;

  float MASS=125.03;
  //******************************************//
  // Fit signal with model pdfs
  //******************************************//
  // four categories to fit
  RooDataSet* sigToFit[ncat];
  RooAbsPdf* mggSig[ncat];
  RooAbsPdf* mjjSig[ncat];
  RooProdPdf* SigPdf[ncat];
  // fit range
  Float_t minSigFitMgg(115),maxSigFitMgg(135);
  Float_t minSigFitMjj(75),maxSigFitMjj(180);
  RooRealVar* mgg = w->var("mgg");
  RooRealVar* mjj = w->var("mjj");
  mgg->setRange("SigFitRange",minSigFitMgg,maxSigFitMgg);
  mjj->setRange("SigFitRange",minSigFitMjj,maxSigFitMjj);

  for (int c = 0; c < ncat; ++c) {
    // import sig and data from workspace
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    mggSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggSig_cat%d",c));
    mjjSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mjjSig_cat%d",c));
    SigPdf[c] = new RooProdPdf(TString::Format("SigPdf_cat%d",c),"",RooArgSet(*mggSig[c], *mjjSig[c]));

    ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(MASS);
    //RooRealVar* peak = w->var(TString::Format("mgg_sig_m0_cat%d",c));
    //peak->setVal(MASS);

    SigPdf[c]->fitTo(*sigToFit[c],Range("SigFitRange"),SumW2Error(kTRUE));

    double mPeak = ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal()+(MASS-125.0); // shift the peak
    ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(mPeak); // shift the peak


    // IMPORTANT: fix all pdf parameters to constant, why?
    RooArgSet sigParams( *w->var(TString::Format("mgg_sig_m0_cat%d",c)),
			 *w->var(TString::Format("mgg_sig_sigma_cat%d",c)),
			 *w->var(TString::Format("mgg_sig_alpha_cat%d",c)),
			 *w->var(TString::Format("mgg_sig_n_cat%d",c)),
			 *w->var(TString::Format("mgg_sig_gsigma_cat%d",c)),
			 *w->var(TString::Format("mgg_sig_frac_cat%d",c)));
    sigParams.add(RooArgSet(
			   *w->var(TString::Format("mjj_sig_m0_cat%d",c)),
		           *w->var(TString::Format("mjj_sig_sigma_cat%d",c)),
		           *w->var(TString::Format("mjj_sig_alpha_cat%d",c)),
		           *w->var(TString::Format("mjj_sig_n_cat%d",c)),
		           *w->var(TString::Format("mjj_sig_gsigma_cat%d",c)),
		           *w->var(TString::Format("mjj_sig_frac_cat%d",c))) );

    w->defineSet(TString::Format("SigPdfParam_cat%d",c), sigParams);
    SetConstantParams(w->set(TString::Format("SigPdfParam_cat%d",c)));

    w->import(*SigPdf[c]);

  } // close for ncat
} // close signal model fit



RooAbsPdf *BkgMggModelFit(RooWorkspace* w, int c, int modelNum) {

  std::vector<TString> catdesc;

  catdesc.push_back("#scale[0.8]{cat0}");
  catdesc.push_back("#scale[0.8]{cat1}");
  catdesc.push_back("#scale[0.8]{cat2}");
  catdesc.push_back("#scale[0.8]{cat3}");

  RooDataSet* data[9];
  RooFitResult* fitresult[9];;
  RooPlot* plotMggBkg[9];

  Float_t minMassFit(100),maxMassFit(180); 
  if(resMass == 270) maxMassFit = 155;
  if(resMass == 300) maxMassFit = 175;

  RooRealVar* mGG     = w->var("mgg");
  mGG->setUnit("GeV");
  RooRealVar* mJJ     = w->var("mjj");
  mJJ->setUnit("GeV");

  mGG->setRange("fitrange",minMassFit,maxMassFit);
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);

  data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
  
  RooFormulaVar *p1modMgg = new RooFormulaVar(TString::Format("mggp1modMgg_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
  RooFormulaVar *p2modMgg = new RooFormulaVar(TString::Format("mggp2modMgg_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c)));
  RooFormulaVar *p3modMgg = new RooFormulaVar(TString::Format("mggp3modMgg_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c)));
  RooFormulaVar *p4modMgg = new RooFormulaVar(TString::Format("mggp4modMgg_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope4_cat%d",c)));
  RooFormulaVar *p5modMgg = new RooFormulaVar(TString::Format("mggp5modMgg_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope5_cat%d",c)));
  RooFormulaVar *p1argMgg = new RooFormulaVar(TString::Format("mggp1argMgg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_arg1_cat%d",c)));
  RooFormulaVar *p2argMgg = new RooFormulaVar(TString::Format("mggp2argMgg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_arg2_cat%d",c)));
  RooFormulaVar *p3argMgg = new RooFormulaVar(TString::Format("mggp3argMgg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_arg3_cat%d",c)));
  RooFormulaVar *p4argMgg = new RooFormulaVar(TString::Format("mggp4argMgg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_arg4_cat%d",c)));
  RooFormulaVar *p5argMgg = new RooFormulaVar(TString::Format("mggp5argMgg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_arg5_cat%d",c)));
  RooFormulaVar *p1widMgg = new RooFormulaVar(TString::Format("mggp1widMgg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_wid1_cat%d",c)));
  RooFormulaVar *p2widMgg = new RooFormulaVar(TString::Format("mggp2widMgg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_wid2_cat%d",c)));
  RooFormulaVar *p3widMgg = new RooFormulaVar(TString::Format("mggp3widMgg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_wid3_cat%d",c)));
  RooFormulaVar *p4widMgg = new RooFormulaVar(TString::Format("mggp4widMgg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_wid4_cat%d",c)));
  RooFormulaVar *p5widMgg = new RooFormulaVar(TString::Format("mggp5widMgg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_wid5_cat%d",c)));

  RooExponential* expo1Mgg = new RooExponential("exp1Mgg","",*mGG,*p1argMgg);
  RooExponential* expo2Mgg = new RooExponential("exp2Mgg","",*mGG,*p2argMgg);
  RooExponential* expo3Mgg = new RooExponential("exp3Mgg","",*mGG,*p3argMgg);
  RooExponential* expo4Mgg = new RooExponential("exp4Mgg","",*mGG,*p4argMgg);
  RooExponential* expo5Mgg = new RooExponential("exp5Mgg","",*mGG,*p5argMgg);
  RooLandau* lan1Mgg = new RooLandau("lan1Mgg","",*mGG,*p1argMgg,*p1widMgg);
  RooLandau* lan2Mgg = new RooLandau("lan2Mgg","",*mGG,*p2argMgg,*p2widMgg);
  RooLandau* lan3Mgg = new RooLandau("lan3Mgg","",*mGG,*p3argMgg,*p3widMgg);
  RooLandau* lan4Mgg = new RooLandau("lan4Mgg","",*mGG,*p4argMgg,*p4widMgg);
  RooLandau* lan5Mgg = new RooLandau("lan5Mgg","",*mGG,*p5argMgg,*p5widMgg);

  const int totalNDOF=5;
  int NDOF[totalNDOF]={0};
  float minNLL[totalNDOF]={0.0}, chi2prob[totalNDOF]={-1.0};
  int bestN = 0;

  RooAbsPdf* MggBkgTmp[totalNDOF] = {0};
  char fitName[10];

  switch (modelNum){

  case 0: //Bernstein
    MggBkgTmp[0] = new RooBernstein("BerN0Mgg", "", *mGG, RooArgList(*p1modMgg));
    MggBkgTmp[1] = new RooBernstein("BerN1Mgg", "", *mGG, RooArgList(*p1modMgg,*p2modMgg));
    MggBkgTmp[2] = new RooBernstein("BerN2Mgg", "", *mGG, RooArgList(*p1modMgg,*p2modMgg,*p3modMgg));
    MggBkgTmp[3] = new RooBernstein("BerN3Mgg", "", *mGG, RooArgList(*p1modMgg,*p2modMgg,*p3modMgg,*p4modMgg));
    MggBkgTmp[4] = new RooBernstein("BerN4Mgg", "", *mGG, RooArgList(*p1modMgg,*p2modMgg,*p3modMgg,*p4modMgg,*p5modMgg));
    sprintf(fitName,"Ber");
    break;

  case 1: //Exponential
    w->factory(TString::Format("mgg_bkg_8TeV_norm_cat%d[50.0,0.0,100000]",c));
    MggBkgTmp[0] = new RooExtendPdf("ExpN1Mgg","",*expo1Mgg,*w->var(TString::Format("mgg_bkg_8TeV_norm_cat%d",c)));
    MggBkgTmp[1] = new RooAddPdf("ExpN2Mgg", "", RooArgList(*expo1Mgg,*expo2Mgg), RooArgList(*p1modMgg,*p2modMgg));
    MggBkgTmp[2] = new RooAddPdf("ExpN3Mgg", "", RooArgList(*expo1Mgg,*expo2Mgg,*expo3Mgg), RooArgList(*p1modMgg,*p2modMgg,*p3modMgg));
    MggBkgTmp[3] = new RooAddPdf("ExpN4Mgg", "", RooArgList(*expo1Mgg,*expo2Mgg,*expo3Mgg,*expo4Mgg), RooArgList(*p1modMgg,*p2modMgg,*p3modMgg,*p4modMgg));
    MggBkgTmp[4] = new RooAddPdf("ExpN5Mgg", "", RooArgList(*expo1Mgg,*expo2Mgg,*expo3Mgg,*expo4Mgg,*expo5Mgg), RooArgList(*p1modMgg,*p2modMgg,*p3modMgg,*p4modMgg,*p5modMgg));
    sprintf(fitName,"Exp");
    sprintf(fitName,"Exp");
    break;

  case 2: //Landau
    w->factory(TString::Format("mgg_bkg_8TeV_norm_cat%d[800.0,0.0,100000]",c));
    MggBkgTmp[0] = new RooExtendPdf("LanN1Mgg","",*lan1Mgg,*w->var(TString::Format("mgg_bkg_8TeV_norm_cat%d",c)));
    MggBkgTmp[1] = new RooAddPdf("LanN2Mgg", "", RooArgList(*lan1Mgg,*lan2Mgg), RooArgList(*p1modMgg,*p2modMgg));
    MggBkgTmp[2] = new RooAddPdf("LanN3Mgg", "", RooArgList(*lan1Mgg,*lan2Mgg,*lan3Mgg), RooArgList(*p1modMgg,*p2modMgg,*p3modMgg));
    MggBkgTmp[3] = new RooAddPdf("LanN4Mgg", "", RooArgList(*lan1Mgg,*lan2Mgg,*lan3Mgg,*lan4Mgg), RooArgList(*p1modMgg,*p2modMgg,*p3modMgg,*p4modMgg));
    MggBkgTmp[4] = new RooAddPdf("LanN5Mgg", "", RooArgList(*lan1Mgg,*lan2Mgg,*lan3Mgg,*lan4Mgg,*lan5Mgg), RooArgList(*p1modMgg,*p2modMgg,*p3modMgg,*p4modMgg,*p5modMgg));
    sprintf(fitName,"Lan");
    break;

  case 3: //Laurent
    MggBkgTmp[0] = new RooGenericPdf("LauN1Mgg","@1*pow(@0,-4)",RooArgList(*mGG,*p1modMgg));
    MggBkgTmp[1] = new RooGenericPdf("LauN2Mgg","@1*pow(@0,-4)+@2*pow(@0,-3)",RooArgList(*mGG,*p1modMgg,*p2modMgg));
    MggBkgTmp[2] = new RooGenericPdf("LauN3Mgg","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)",RooArgList(*mGG,*p1modMgg,*p2modMgg,*p3modMgg));
    MggBkgTmp[3] = new RooGenericPdf("LauN4Mgg","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)+@4*pow(@0,-2)",RooArgList(*mGG,*p1modMgg,*p2modMgg,*p3modMgg,*p4modMgg));
    MggBkgTmp[4] = new RooGenericPdf("LauN5Mgg","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)+@4*pow(@0,-2)+@5*pow(@0,-6)",RooArgList(*mGG,*p1modMgg,*p2modMgg,*p3modMgg,*p4modMgg,*p5modMgg));
    sprintf(fitName,"Lau");
    sprintf(fitName,"Lau");
    break;

  case 4: //Power
    MggBkgTmp[0] = new RooGenericPdf("PowN1Mgg","@1*pow(@0,@2)",RooArgList(*mGG,*p1modMgg,*p1argMgg));
    MggBkgTmp[1] = new RooGenericPdf("PowN2Mgg","@1*pow(@0,@2)+@3*pow(@0,@4)",RooArgList(*mGG,*p1modMgg,*p1argMgg,*p2modMgg,*p2argMgg));
    MggBkgTmp[2] = new RooGenericPdf("PowN3Mgg","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)",RooArgList(*mGG,*p1modMgg,*p1argMgg,*p2modMgg,*p2argMgg,*p3modMgg,*p3argMgg));
    MggBkgTmp[3] = new RooGenericPdf("PowN4Mgg","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)",RooArgList(*mGG,*p1modMgg,*p1argMgg,*p2modMgg,*p2argMgg,*p3modMgg,*p3argMgg,*p4modMgg,*p4argMgg));
    MggBkgTmp[4] = new RooGenericPdf("PowN5Mgg","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)",RooArgList(*mGG,*p1modMgg,*p1argMgg,*p2modMgg,*p2argMgg,*p3modMgg,*p3argMgg,*p4modMgg,*p4argMgg));
    //MggBkgTmp[4] = new RooGenericPdf(TString::Format("MggBkg_cat%d",c),"@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)+@9*pow(@0,@10)",RooArgList(*mGG,*p1modMgg,*p1argMgg,*p2modMgg,*p2argMgg,*p3modMgg,*p3argMgg,*p4modMgg,*p4argMgg,*p5modMgg,*p5argMgg));
    sprintf(fitName,"Pow");
    break;
  }

  bool breakLoop=false;
  for(int i=0; i<totalNDOF; ++i){
    
    fitresult[c] = MggBkgTmp[i]->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE));//strategy 1 or 2?
    //w->import(*MggBkgTmp[i]);
    
    minNLL[i] = fitresult[c]->minNll();
    switch(modelNum){
    case 0: NDOF[i] = i+1; break;
    case 1: NDOF[i] = 2*(i+1); break;
    case 2: NDOF[i] = 3*(i+1); break;
    case 3: NDOF[i] = i+1; break;
    case 4: NDOF[i] = 2*(i+1); break;
    }
    //************************************************//
    // Plot Mgg background fit results per categories 
    //************************************************//
    // Plot Background Categories 
    //****************************//
    
    TCanvas* ctmp = new TCanvas("ctmp","Mgg Background Categories",0,0,500,500);
    Int_t nBinsMass(80);
    plotMggBkg[c] = mGG->frame();
    data[c]->plotOn(plotMggBkg[c],LineColor(kWhite),MarkerColor(kWhite));    
    MggBkgTmp[i]->plotOn(plotMggBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 
    data[c]->plotOn(plotMggBkg[c]);    

    plotMggBkg[c]->SetTitle("");      
    plotMggBkg[c]->SetMinimum(1e-5);
    plotMggBkg[c]->SetMaximum(1.40*plotMggBkg[c]->GetMaximum());
    plotMggBkg[c]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    plotMggBkg[c]->Draw();  

    TLegend *legmc = new TLegend(0.55,0.75,0.98,0.9);
    legmc->AddEntry(plotMggBkg[c]->getObject(2),"Data CS","LPE");
    legmc->AddEntry(plotMggBkg[c]->getObject(1),TString::Format("%.5s Truth",MggBkgTmp[i]->GetName()),"L");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    

    TLatex *lat  = new TLatex(minMassFit+3.0,0.85*plotMggBkg[c]->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
    lat->Draw();
    TLatex *lat2 = new TLatex(minMassFit+3.0,0.7*plotMggBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
 
    ctmp->SaveAs(TString::Format("plots/dataBkgMgg_%.5s_cat%d_mass%d.png",MggBkgTmp[i]->GetName(),c,resMass));
    ctmp->SaveAs(TString::Format("plots/dataBkgMgg_%.5s_cat%d_mass%d.pdf",MggBkgTmp[i]->GetName(),c,resMass));
    delete ctmp;

    if(i>0){
      float chi2 = 2*(minNLL[i-1]-minNLL[i]);
      int chi2dof = NDOF[i]-NDOF[i-1];
      chi2prob[i-1] = chi2<0 ? 1.0 : TMath::Prob(chi2,chi2dof);
      if(chi2prob[i-1]>0.05 ) breakLoop=true;
    }
    if(data[c]->sumEntries()-1<=NDOF[i]) breakLoop=true;
    if(!breakLoop)
      bestN=i;
  }

  FILE *results = fopen("resultsModel2DMgg.txt","a");
  fprintf(results,"### Mgg spectrum, %s fit for cat%i. The best order is %i.\n",fitName,c,bestN+(modelNum!=0));//for Bernstein, N is the order of the polynomial fit.

  for(int i=0; i<totalNDOF; ++i){
    if(NDOF[i]==0) break;
    fprintf(results,"N=%i NLL=%.3f NDOF=%i\n",i+(modelNum!=0),minNLL[i],NDOF[i]);
  }
  for(int i=1; i<totalNDOF; ++i){
    if(NDOF[i]==0) break;
    float chi2 = 2*(minNLL[i-1]-minNLL[i]);
    int chi2dof = NDOF[i]-NDOF[i-1];
    fprintf(results,"N=%i chi2=%.3f chi2dof=%i chi2prob=%.3f\n",i-1+(modelNum!=0),chi2,chi2dof,chi2prob[i-1]);
  }
  fprintf(results,"\n\n");
  fclose(results);

  //redo the fit since parameters are used across multiple truth models
  MggBkgTmp[bestN]->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE));//strategy 1 or 2?
  BkgMggModelFitNoCorr = MggBkgTmp[bestN];
  NDOFMgg_noCorr = NDOF[bestN];
  
  RooRealVar* parCondBern[10];
  for(int ii = 0; ii < 10; ii++){
    parCondBern[ii] = new RooRealVar(TString::Format("parCondBern_%d",ii),"",0.1,-100.,100.);
  }
  RooRealVar* parCond[5];
  for(int ii = 0; ii < 5; ii++){
    parCond[ii] = new RooRealVar(TString::Format("parCond_%d",ii),"",0.1,-10.,10.);
  }
  RooRealVar* argCond[5];
  for(int ii = 0; ii < 5; ii++){
    argCond[ii] = new RooRealVar(TString::Format("argCond_%d",ii),"",-0.1,-10.,0.);
  }
  RooRealVar* corrCond[5];
  for(int ii = 0; ii < 5; ii++){
    corrCond[ii] = new RooRealVar(TString::Format("corrCond_%d",ii),"",0.001,-1.,1.);
  }
  
  RooArgList arglist(*mGG,*mJJ);
  switch (modelNum){

  case 0: //Bernstein
      arglist.add(*parCondBern[0]);
      arglist.add(*parCondBern[1]);
      MggBkgTmp[0] = new RooGenericPdf("BerN0Mgg", "(@3+@1*@2)", arglist);
      arglist.add(*parCondBern[2]);
      arglist.add(*parCondBern[3]);
      MggBkgTmp[1] = new RooGenericPdf("BerN1Mgg", "(@3+@1*@2)*(1-@0)+(@5+@1*@4)*(@0)",arglist);
      arglist.add(*parCondBern[4]);
      arglist.add(*parCondBern[5]);
      MggBkgTmp[2] = new RooGenericPdf("BerN2Mgg", "(@3+@1*@2)*(1-@0)**2+(@5+@1*@4)*2*(@0)*(1-@0)+(@7+@1*@6)*(@0)**2",arglist); 
      arglist.add(*parCondBern[6]);
      arglist.add(*parCondBern[7]);
      MggBkgTmp[3] = new RooGenericPdf("BerN3Mgg", "(@3+@1*@2)*(1-@0)**3+(@5+@1*@4)*3*(@0)*(1-@0)**2+(@7+@1*@6)*3*(@0)**2*(1-@0)+(@9+@1*@8)*(@0)**3",arglist);
      arglist.add(*parCondBern[8]);
      arglist.add(*parCondBern[9]);
      MggBkgTmp[4] = new RooGenericPdf("BerN4Mgg", "(@3+@1*@2)*(1-@0)**4+(@5+@1*@4)*4*(@0)*(1-@0)**3+(@7+@1*@6)*6*(@0)**2*(1-@0)**2+(@9+@1*@8)*4*(@0)**3*(1-@0)+(@11+@1*@10)*@0**4",arglist);
  break;

  arglist.Clear();
  arglist.add(*mGG);
  arglist.add(*mJJ);
  case 1: //Exponential
      arglist.add(*corrCond[0]);
      arglist.add(*argCond[0]);
      arglist.add(*parCond[0]);
      MggBkgTmp[0] = new RooGenericPdf("ExpN1Mgg", "@4*exp((@3+@1*@2)*@0*0.0001)",arglist);
      arglist.add(*corrCond[1]);
      arglist.add(*argCond[1]);
      arglist.add(*parCond[1]);
      MggBkgTmp[1] = new RooGenericPdf("ExpN2Mgg", "@4*exp((@3+@1*@2)*@0*0.0001)+@7*exp((@6+@1*@5)*@0*0.0001)", arglist);
      arglist.add(*corrCond[2]);
      arglist.add(*argCond[2]);
      arglist.add(*parCond[2]);
      MggBkgTmp[2] = new RooGenericPdf("ExpN3Mgg", "@4*exp((@3+@1*@2)*@0*0.0001)+@7*exp((@6+@1*@5)*@0*0.0001)+@10*exp((@9+@1*@8)*@0*0.0001)", arglist);
      arglist.add(*corrCond[3]);
      arglist.add(*argCond[3]);
      arglist.add(*parCond[3]);
      MggBkgTmp[3] = new RooGenericPdf("ExpN4Mgg", "@4*exp((@3+@1*@2)*@0*0.0001)+@7*exp((@6+@1*@5)*@0*0.0001)+@10*exp((@9+@1*@8)*@0*0.0001)+@13*exp((@12+@1*@11)*@0*0.0001)", arglist);
      arglist.add(*corrCond[4]);
      arglist.add(*argCond[4]);
      arglist.add(*parCond[4]);
      MggBkgTmp[4] = new RooGenericPdf("ExpN5Mgg", "@4*exp((@3+@1*@2)*@0*0.0001)+@7*exp((@6+@1*@5)*@0*0.0001)+@10*exp((@9+@1*@8)*@0*0.0001)+@13*exp((@12+@1*@11)*@0*0.0001)+@16*exp((@15+@1*@14)*@0*0.0001)", arglist);
  break;

  arglist.Clear();
  arglist.add(*mGG);
  arglist.add(*mJJ);
  case 4: //Power
      arglist.add(*corrCond[0]);
      arglist.add(*argCond[0]);
      arglist.add(*parCond[0]);
      MggBkgTmp[0] = new RooGenericPdf("PowN1Mgg", "@4*exp((@3+@1*@2)*log(@)*0.0001)",arglist);
      arglist.add(*corrCond[1]);
      arglist.add(*argCond[1]);
      arglist.add(*parCond[1]);
      MggBkgTmp[1] = new RooGenericPdf("PowN2Mgg", "@4*exp((@3+@1*@2)*log(@)*0.0001)+@7*exp((@6+@1*@5)*log(@)*0.0001)", arglist);
      arglist.add(*corrCond[2]);
      arglist.add(*argCond[2]);
      arglist.add(*parCond[2]);
      MggBkgTmp[2] = new RooGenericPdf("PowN3Mgg", "@4*exp((@3+@1*@2)*log(@)*0.0001)+@7*exp((@6+@1*@5)*log(@)*0.0001)+@10*exp((@9+@1*@8)*log(@)*0.0001)", arglist);
      arglist.add(*corrCond[3]);
      arglist.add(*argCond[3]);
      arglist.add(*parCond[3]);
      MggBkgTmp[3] = new RooGenericPdf("PowN4Mgg", "@4*exp((@3+@1*@2)*log(@)*0.0001)+@7*exp((@6+@1*@5)*log(@)*0.0001)+@10*exp((@9+@1*@8)*log(@)*0.0001)+@13*exp((@12+@1*@11)*log(@)*0.0001)", arglist);
      arglist.add(*corrCond[4]);
      arglist.add(*argCond[4]);
      arglist.add(*parCond[4]);
      MggBkgTmp[4] = new RooGenericPdf("PowN5Mgg", "@4*exp((@3+@1*@2)*log(@)*0.0001)+@7*exp((@6+@1*@5)*log(@)*0.0001)+@10*exp((@9+@1*@8)*log(@)*0.0001)+@13*exp((@12+@1*@11)*log(@)*0.0001)+@16*exp((@15+@1*@14)*log(@)*0.0001)", arglist);
  break;
  }
  
  for(int i = 0; i < totalNDOF; i++)
  {
    switch(modelNum){
    case 0: NDOF[i] = 2*(i+1); break;
    case 1: NDOF[i] = 3*(i+1); break;
    case 4: NDOF[i] = 3*(i+1); break;
    }
  }
  NDOFMgg = NDOF[bestN];
  //w->import(*MggBkgTmp[bestN]);

  return MggBkgTmp[bestN];
}


RooAbsPdf *BkgMjjModelFit(RooWorkspace* w, int c, int modelNum) {

  std::vector<TString> catdesc;

  catdesc.push_back("#scale[0.8]{cat0}");
  catdesc.push_back("#scale[0.8]{cat1}");
  catdesc.push_back("#scale[0.8]{cat2}");
  catdesc.push_back("#scale[0.8]{cat3}");

  RooDataSet* data[9];
  RooFitResult* fitresult[9];;
  RooPlot* plotMjjBkg[9];

  Float_t minMassFit(60),maxMassFit(180); 
  if(resMass == 0) minMassFit = 75;

  RooRealVar* mJJ     = w->var("mjj");
  mJJ->setUnit("GeV");
  mJJ->setRange("fitrange",minMassFit,maxMassFit); 
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);

  data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));

  RooFormulaVar *p1modMjj = new RooFormulaVar(TString::Format("mjjp1modMjj_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope1_cat%d",c)));
  RooFormulaVar *p2modMjj = new RooFormulaVar(TString::Format("mjjp2modMjj_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope2_cat%d",c)));
  RooFormulaVar *p3modMjj = new RooFormulaVar(TString::Format("mjjp3modMjj_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope3_cat%d",c)));
  RooFormulaVar *p4modMjj = new RooFormulaVar(TString::Format("mjjp4modMjj_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope4_cat%d",c)));
  RooFormulaVar *p5modMjj = new RooFormulaVar(TString::Format("mjjp5modMjj_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope5_cat%d",c)));
  RooFormulaVar *p1argMjj = new RooFormulaVar(TString::Format("mjjp1argMjj_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_arg1_cat%d",c)));
  RooFormulaVar *p2argMjj = new RooFormulaVar(TString::Format("mjjp2argMjj_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_arg2_cat%d",c)));
  RooFormulaVar *p3argMjj = new RooFormulaVar(TString::Format("mjjp3argMjj_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_arg3_cat%d",c)));
  RooFormulaVar *p4argMjj = new RooFormulaVar(TString::Format("mjjp4argMjj_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_arg4_cat%d",c)));
  RooFormulaVar *p5argMjj = new RooFormulaVar(TString::Format("mjjp5argMjj_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_arg5_cat%d",c)));
  RooFormulaVar *p1widMjj = new RooFormulaVar(TString::Format("mjjp1widMjj_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_wid1_cat%d",c)));
  RooFormulaVar *p2widMjj = new RooFormulaVar(TString::Format("mjjp2widMjj_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_wid2_cat%d",c)));
  RooFormulaVar *p3widMjj = new RooFormulaVar(TString::Format("mjjp3widMjj_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_wid3_cat%d",c)));
  RooFormulaVar *p4widMjj = new RooFormulaVar(TString::Format("mjjp4widMjj_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_wid4_cat%d",c)));
  RooFormulaVar *p5widMjj = new RooFormulaVar(TString::Format("mjjp5widMjj_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_wid5_cat%d",c)));

  RooExponential* expo1Mjj = new RooExponential("exp1Mjj","",*mJJ,*p1argMjj);
  RooExponential* expo2Mjj = new RooExponential("exp2Mjj","",*mJJ,*p2argMjj);
  RooExponential* expo3Mjj = new RooExponential("exp3Mjj","",*mJJ,*p3argMjj);
  RooExponential* expo4Mjj = new RooExponential("exp4Mjj","",*mJJ,*p4argMjj);
  RooExponential* expo5Mjj = new RooExponential("exp5Mjj","",*mJJ,*p5argMjj);
  RooLandau* lan1Mjj = new RooLandau("lan1Mjj","",*mJJ,*p1argMjj,*p1widMjj);
  RooLandau* lan2Mjj = new RooLandau("lan2Mjj","",*mJJ,*p2argMjj,*p2widMjj);
  RooLandau* lan3Mjj = new RooLandau("lan3Mjj","",*mJJ,*p3argMjj,*p3widMjj);
  RooLandau* lan4Mjj = new RooLandau("lan4Mjj","",*mJJ,*p4argMjj,*p4widMjj);
  RooLandau* lan5Mjj = new RooLandau("lan5Mjj","",*mJJ,*p5argMjj,*p5widMjj);

  const int totalNDOF=5;
  int NDOF[totalNDOF]={0};
  float minNLL[totalNDOF]={0.0}, chi2prob[totalNDOF]={-1.0};
  int bestN = 0;

  RooAbsPdf* MjjBkgTmp[totalNDOF] = {0};
  char fitName[10];

  switch (modelNum){

  case 0: //Bernstein
    MjjBkgTmp[0] = new RooBernstein("BerN0Mjj", "", *mJJ, RooArgList(*p1modMjj));
    MjjBkgTmp[1] = new RooBernstein("BerN1Mjj", "", *mJJ, RooArgList(*p1modMjj,*p2modMjj));
    MjjBkgTmp[2] = new RooBernstein("BerN2Mjj", "", *mJJ, RooArgList(*p1modMjj,*p2modMjj,*p3modMjj));
    MjjBkgTmp[3] = new RooBernstein("BerN3Mjj", "", *mJJ, RooArgList(*p1modMjj,*p2modMjj,*p3modMjj,*p4modMjj));
    MjjBkgTmp[4] = new RooBernstein("BerN4Mjj", "", *mJJ, RooArgList(*p1modMjj,*p2modMjj,*p3modMjj,*p4modMjj,*p5modMjj));
    sprintf(fitName,"Ber");
    break;

  case 1: //Exponential
    w->factory(TString::Format("mjj_bkg_8TeV_norm_cat%d[800.0,0.0,100000]",c));
    MjjBkgTmp[0] = new RooExtendPdf("ExpN1Mjj","",*expo1Mjj,*w->var(TString::Format("mjj_bkg_8TeV_norm_cat%d",c)));
    MjjBkgTmp[1] = new RooAddPdf("ExpN2Mjj", "", RooArgList(*expo1Mjj,*expo2Mjj), RooArgList(*p1modMjj,*p2modMjj));
    MjjBkgTmp[2] = new RooAddPdf("ExpN3Mjj", "", RooArgList(*expo1Mjj,*expo2Mjj,*expo3Mjj), RooArgList(*p1modMjj,*p2modMjj,*p3modMjj));
    MjjBkgTmp[3] = new RooAddPdf("ExpN4Mjj", "", RooArgList(*expo1Mjj,*expo2Mjj,*expo3Mjj,*expo4Mjj), RooArgList(*p1modMjj,*p2modMjj,*p3modMjj,*p4modMjj));
    MjjBkgTmp[4] = new RooAddPdf("ExpN5Mjj", "", RooArgList(*expo1Mjj,*expo2Mjj,*expo3Mjj,*expo4Mjj,*expo5Mjj), RooArgList(*p1modMjj,*p2modMjj,*p3modMjj,*p4modMjj,*p5modMjj));
    sprintf(fitName,"Exp");
    break;

  case 2: //Landau
    w->factory(TString::Format("mjj_bkg_8TeV_norm_cat%d[800.0,0.0,100000]",c));
    MjjBkgTmp[0] = new RooExtendPdf("LanN1Mjj","",*lan1Mjj,*w->var(TString::Format("mjj_bkg_8TeV_norm_cat%d",c)));
    MjjBkgTmp[1] = new RooAddPdf("LanN2Mjj", "", RooArgList(*lan1Mjj,*lan2Mjj), RooArgList(*p1modMjj,*p2modMjj));
    MjjBkgTmp[2] = new RooAddPdf("LanN3Mjj", "", RooArgList(*lan1Mjj,*lan2Mjj,*lan3Mjj), RooArgList(*p1modMjj,*p2modMjj,*p3modMjj));
    MjjBkgTmp[3] = new RooAddPdf("LanN4Mjj", "", RooArgList(*lan1Mjj,*lan2Mjj,*lan3Mjj,*lan4Mjj), RooArgList(*p1modMjj,*p2modMjj,*p3modMjj,*p4modMjj));
    MjjBkgTmp[4] = new RooAddPdf("LanN5Mjj", "", RooArgList(*lan1Mjj,*lan2Mjj,*lan3Mjj,*lan4Mjj,*lan5Mjj), RooArgList(*p1modMjj,*p2modMjj,*p3modMjj,*p4modMjj,*p5modMjj));
    sprintf(fitName,"Lan");
    break;

  case 3: //Laurent
    MjjBkgTmp[0] = new RooGenericPdf("LauN1Mjj","@1*pow(@0,-4)",RooArgList(*mJJ,*p1modMjj));
    MjjBkgTmp[1] = new RooGenericPdf("LauN2Mjj","@1*pow(@0,-4)+@2*pow(@0,-3)",RooArgList(*mJJ,*p1modMjj,*p2modMjj));
    MjjBkgTmp[2] = new RooGenericPdf("LauN3Mjj","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)",RooArgList(*mJJ,*p1modMjj,*p2modMjj,*p3modMjj));
    MjjBkgTmp[3] = new RooGenericPdf("LauN4Mjj","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)+@4*pow(@0,-2)",RooArgList(*mJJ,*p1modMjj,*p2modMjj,*p3modMjj,*p4modMjj));
    MjjBkgTmp[4] = new RooGenericPdf("LauN5Mjj","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)+@4*pow(@0,-2)+@5*pow(@0,-6)",RooArgList(*mJJ,*p1modMjj,*p2modMjj,*p3modMjj,*p4modMjj,*p5modMjj));
    sprintf(fitName,"Lau");
    break;

  case 4: //Power
    MjjBkgTmp[0] = new RooGenericPdf("PowN1Mjj","@1*pow(@0,@2)",RooArgList(*mJJ,*p1modMjj,*p1argMjj));
    MjjBkgTmp[1] = new RooGenericPdf("PowN2Mjj","@1*pow(@0,@2)+@3*pow(@0,@4)",RooArgList(*mJJ,*p1modMjj,*p1argMjj,*p2modMjj,*p2argMjj));
    MjjBkgTmp[2] = new RooGenericPdf("PowN3Mjj","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)",RooArgList(*mJJ,*p1modMjj,*p1argMjj,*p2modMjj,*p2argMjj,*p3modMjj,*p3argMjj));
    MjjBkgTmp[3] = new RooGenericPdf("PowN4Mjj","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)",RooArgList(*mJJ,*p1modMjj,*p1argMjj,*p2modMjj,*p2argMjj,*p3modMjj,*p3argMjj,*p4modMjj,*p4argMjj));
    MjjBkgTmp[4] = new RooGenericPdf("PowN5Mjj","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)",RooArgList(*mJJ,*p1modMjj,*p1argMjj,*p2modMjj,*p2argMjj,*p3modMjj,*p3argMjj,*p4modMjj,*p4argMjj));
    //MjjBkgTmp[4] = new RooGenericPdf(TString::Format("MjjBkg_cat%d",c),"@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)+@9*pow(@0,@10)",RooArgList(*mJJ,*p1modMjj,*p1argMjj,*p2modMjj,*p2argMjj,*p3modMjj,*p3argMjj,*p4modMjj,*p4argMjj,*p5modMjj,*p5argMjj));
    sprintf(fitName,"Pow");
    break;
  }

  bool breakLoop=false;
  for(int i=0; i<totalNDOF; ++i){
    
    fitresult[c] = MjjBkgTmp[i]->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE));//strategy 1 or 2?
    //w->import(*MjjBkgTmp[i]);
    
    minNLL[i] = fitresult[c]->minNll();
    switch(modelNum){
    case 0: NDOF[i] = i+1; break;
    case 1: NDOF[i] = 2*(i+1); break;
    case 2: NDOF[i] = 3*(i+1); break;
    case 3: NDOF[i] = i+1; break;
    case 4: NDOF[i] = 2*(i+1); break;
    }
    //************************************************//
    // Plot Mjj background fit results per categories 
    //************************************************//
    // Plot Background Categories 
    //****************************//
    
    TCanvas* ctmp = new TCanvas("ctmp","Mjj Background Categories",0,0,500,500);
    Int_t nBinsMass(80);
    plotMjjBkg[c] = mJJ->frame();
    data[c]->plotOn(plotMjjBkg[c],LineColor(kWhite),MarkerColor(kWhite));    
    MjjBkgTmp[i]->plotOn(plotMjjBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 
    data[c]->plotOn(plotMjjBkg[c]);    

    plotMjjBkg[c]->SetTitle("");      
    plotMjjBkg[c]->SetMinimum(1e-5);
    plotMjjBkg[c]->SetMaximum(1.40*plotMjjBkg[c]->GetMaximum());
    plotMjjBkg[c]->GetXaxis()->SetTitle("m_{jj} (GeV)");
    plotMjjBkg[c]->Draw();  

    TLegend *legmc = new TLegend(0.55,0.75,0.98,0.9);
    legmc->AddEntry(plotMjjBkg[c]->getObject(2),"Data CS","LPE");
    legmc->AddEntry(plotMjjBkg[c]->getObject(1),TString::Format("%.5s Truth",MjjBkgTmp[i]->GetName()),"L");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    

    TLatex *lat  = new TLatex(minMassFit+3.0,0.85*plotMjjBkg[c]->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
    lat->Draw();
    TLatex *lat2 = new TLatex(minMassFit+3.0,0.7*plotMjjBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
 
    ctmp->SaveAs(TString::Format("plots/dataBkgMjj_%.5s_cat%d_mass%d.png",MjjBkgTmp[i]->GetName(),c,resMass));
    ctmp->SaveAs(TString::Format("plots/dataBkgMjj_%.5s_cat%d_mass%d.pdf",MjjBkgTmp[i]->GetName(),c,resMass));
    delete ctmp;

    if(i>0){
      float chi2 = 2*(minNLL[i-1]-minNLL[i]);
      int chi2dof = NDOF[i]-NDOF[i-1];
      chi2prob[i-1] = chi2<0 ? 1.0 : TMath::Prob(chi2,chi2dof);
      if(chi2prob[i-1]>0.05 ) breakLoop=true;
    }
    if(data[c]->sumEntries()-1<=NDOF[i]) breakLoop=true;
    if(!breakLoop)
      bestN=i;
  }

  FILE *results = fopen("resultsModel2DMjj.txt","a");
  fprintf(results,"### Mjj spectrum, %s fit for cat%i. The best order is %i.\n",fitName,c,bestN+(modelNum!=0));//for Bernstein, N is the order of the polynomial fit.

  for(int i=0; i<totalNDOF; ++i){
    if(NDOF[i]==0) break;
    fprintf(results,"N=%i NLL=%.3f NDOF=%i\n",i+(modelNum!=0),minNLL[i],NDOF[i]);
  }
  for(int i=1; i<totalNDOF; ++i){
    if(NDOF[i]==0) break;
    float chi2 = 2*(minNLL[i-1]-minNLL[i]);
    int chi2dof = NDOF[i]-NDOF[i-1];
    fprintf(results,"N=%i chi2=%.3f chi2dof=%i chi2prob=%.3f\n",i-1+(modelNum!=0),chi2,chi2dof,chi2prob[i-1]);
  }
  fprintf(results,"\n\n");
  fclose(results);

  MjjBkgTmp[bestN]->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE));//strategy 1 or 2?
  w->import(*MjjBkgTmp[bestN]);
  NDOFMjj = NDOF[bestN];

  return MjjBkgTmp[bestN];

}


void BkgModelBias(RooWorkspace* w,int c,RooAbsPdf* MggBkgTruth, RooAbsPdf* MjjBkgTruth, FILE *fout, FILE *foutTotal, FILE *foutCorr){

  std::vector<TString> catdesc;

  catdesc.push_back("#scale[0.8]{cat0}");
  catdesc.push_back("#scale[0.8]{cat1}");
  catdesc.push_back("#scale[0.8]{cat2}");
  catdesc.push_back("#scale[0.8]{cat3}");

  Float_t minMggMassFit(100),maxMggMassFit(180); 
  Float_t minMjjMassFit(60),maxMjjMassFit(180); 
  if(resMass == 270) maxMggMassFit = 155;
  if(resMass == 300) maxMggMassFit = 175;
  if(resMass == 0) minMjjMassFit = 75; 
  float sigMeanMgg=0, sigFWHM_Mgg=0, sigMeanMjj=0, sigFWHM_Mjj=0;
  switch(c){
  case 0: case 2: sigMeanMgg=124.81; sigFWHM_Mgg=1.76; sigMeanMjj=122.20; sigFWHM_Mjj=22.43; break;
  case 1: case 3: sigMeanMgg=124.78; sigFWHM_Mgg=1.92; sigMeanMjj=120.41; sigFWHM_Mjj=31.04; break;
  }

  RooRealVar* mGG     = w->var("mgg");  
  mGG->setUnit("GeV");
  mGG->setRange("sigRegion",sigMeanMgg-sigFWHM_Mgg,sigMeanMgg+sigFWHM_Mgg);
  RooRealVar* mJJ     = w->var("mjj");  
  mJJ->setUnit("GeV");
  mJJ->setRange("sigRegion",sigMeanMjj-sigFWHM_Mjj,sigMeanMjj+sigFWHM_Mjj);
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);

  RooDataSet* data = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
  
  RooFormulaVar *p0 = new RooFormulaVar("p0","","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
  RooFormulaVar *p1 = new RooFormulaVar("p1","","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c)));
  RooFormulaVar *p2 = new RooFormulaVar("p2","","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c)));
  RooFormulaVar *p3 = new RooFormulaVar("p3","","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope1_cat%d",c)));
  RooFormulaVar *p4 = new RooFormulaVar("p4","","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope2_cat%d",c)));
  RooFormulaVar *p5 = new RooFormulaVar("p5","","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope3_cat%d",c)));
  RooRealVar *p6 = new RooRealVar("p6","",-0.1,-10.,0.);
  RooRealVar *p7 = new RooRealVar("p7","",-0.1,-10.,0.);
  RooRealVar *p8 = new RooRealVar("p8","",-0.1,-10.,0.);
  RooRealVar *p9 = new RooRealVar("p9","",-0.1,-10.,0.);

  RooExponential* expo1FitMgg = new RooExponential("exp1FitMgg","",*mGG,*p6);
  RooExponential* expo2FitMgg = new RooExponential("exp2FitMgg","",*mGG,*p7);
  RooExponential* expo1FitMjj = new RooExponential("exp1FitMjj","",*mJJ,*p8);
  RooExponential* expo2FitMjj = new RooExponential("exp2FitMjj","",*mJJ,*p9);
  
  RooGenericPdf* pow1FitMgg = new RooGenericPdf("pow1FitMgg","pow(@0,@1)",RooArgList(*mGG,*p6));
  RooGenericPdf* pow2FitMgg = new RooGenericPdf("pow2FitMgg","pow(@0,@1)",RooArgList(*mGG,*p7));
  RooGenericPdf* pow1FitMjj = new RooGenericPdf("pow1FitMjj","pow(@0,@1)",RooArgList(*mJJ,*p8));
  RooGenericPdf* pow2FitMjj = new RooGenericPdf("pow2FitMjj","pow(@0,@1)",RooArgList(*mJJ,*p9));

  const int totalNDOF=7;
  RooAbsPdf* MggBkgTmp[totalNDOF] = {0};
  RooAbsPdf* MjjBkgTmp[totalNDOF] = {0};
  MggBkgTmp[0] = new RooExponential(TString::Format("ExpN%dMggCand",1), "", *mGG,*p6);
  MggBkgTmp[1] = new RooAddPdf(TString::Format("ExpN%dMggCand",2), "", RooArgList(*expo1FitMgg,*expo2FitMgg), RooArgList(*p6,*p7));
  MggBkgTmp[2] = new RooGenericPdf(TString::Format("PowN%dMggCand",1),"pow(@0,@1)",RooArgList(*mGG,*p6));
  MggBkgTmp[3] = new RooAddPdf(TString::Format("PowN%dMggCand",2), "", RooArgList(*pow1FitMgg,*pow2FitMgg), RooArgList(*p6,*p8));
  MggBkgTmp[4] = new RooBernstein(TString::Format("BerN%dMggCand",1), "", *mGG,RooArgList(*p0,*p1));
  MggBkgTmp[5] = new RooBernstein(TString::Format("BerN%dMggCand",1), "", *mGG,RooArgList(*p0,*p1));
  MggBkgTmp[6] = new RooBernstein(TString::Format("BerN%dMggCand",2), "", *mGG,RooArgList(*p0,*p1,*p2));
  MjjBkgTmp[0] = new RooExponential(TString::Format("ExpN%dMjjCand",1), "", *mJJ,*p8);
  MjjBkgTmp[1] = new RooAddPdf(TString::Format("ExpN%dMjjCand",2), "", RooArgList(*expo1FitMjj,*expo2FitMjj), RooArgList(*p8,*p9));
  MjjBkgTmp[2] = new RooGenericPdf(TString::Format("PowN%dMjjCand",1),"pow(@0,@1)",RooArgList(*mJJ,*p8));
  MjjBkgTmp[3] = new RooAddPdf(TString::Format("PowN%dMjjCand",2), "", RooArgList(*pow1FitMjj,*pow2FitMjj), RooArgList(*p8,*p9));
  MjjBkgTmp[4] = new RooBernstein(TString::Format("BerN%dMjjCand",1), "", *mJJ,RooArgList(*p3,*p4));
  MjjBkgTmp[5] = new RooBernstein(TString::Format("BerN%dMjjCand",2), "", *mJJ,RooArgList(*p3,*p4,*p5));
  MjjBkgTmp[6] = new RooBernstein(TString::Format("BerN%dMjjCand",1), "", *mJJ,RooArgList(*p3,*p4));

  float n_gen_sr=0;

  const int Npse = 20;
  float results[totalNDOF];
  float resultsTotal[totalNDOF][Npse];
  float resultsBkg[totalNDOF];
  TH1F *corrHist = new TH1F("corrHist","corrHist",200,-1,1);

  mGG->setRange("massFit",100,180);
  mJJ->setRange("massFit",60,180); 
    
  if(resMass == 270){
     mGG->setRange("massFit",100,155);
  }
  if(resMass == 300){
     mGG->setRange("massFit",100,175);
  }
  if(resMass == 0){
     mJJ->setRange("massFit",75,180);
  }

  std::cout << "Fit without correlation!" << std::endl;
  RooProdPdf* BkgTruthTmp = new RooProdPdf("BkgTruthTmp","",RooArgList(*BkgMggModelFitNoCorr,*MjjBkgTruth));
  RooRealVar* nbkgTruth = new RooRealVar("nbkgTruth","",data->sumEntries());
  RooExtendPdf* BkgTruth = new RooExtendPdf("BkgTruth","",*BkgTruthTmp,*nbkgTruth);
  RooFitResult* fitResultTruth = BkgTruth->fitTo(*data, Strategy(1),Minos(kFALSE), Range("massFit"),SumW2Error(kTRUE), Save(kTRUE));//strategy 1 or 2?
  w->import(*BkgTruth);
  float minNLL_noCorr = fitResultTruth->minNll();

  std::cout << "Fit with correlation!" << std::endl;
  RooProdPdf* BkgTruthCondtmp = new RooProdPdf("MggBkgTruthCondtmp","",*MjjBkgTruth, Conditional(*MggBkgTruth,*mGG));
  RooRealVar* nbkgTruthCond = new RooRealVar("nbkgTruthCond","",data->sumEntries());
  RooExtendPdf* BkgTruthCond = new RooExtendPdf("BkgTruthCond","",*BkgTruthCondtmp,*nbkgTruthCond);
  RooFitResult* fitResultTruthCorr = BkgTruthCond->fitTo(*data, Strategy(1),Minos(kFALSE), Range("massFit"),SumW2Error(kTRUE), Save(kTRUE));//strategy 1 or 2?
  w->import(*BkgTruthCond);
  float minNLL_Corr = fitResultTruthCorr->minNll();

  bool noCorr = false;
  float chi2 = 2*(minNLL_noCorr-minNLL_Corr);
  std::cout << "chi2 = " << chi2 << std::endl;
  int chi2dof = NDOFMgg-NDOFMgg_noCorr;
  std::cout << "chi2dof = " << chi2dof << std::endl;
  float chi2prob = chi2<0 ? 1.0 : TMath::Prob(chi2,chi2dof);
  if(chi2prob>0.05 ) noCorr=true;
  if(data->sumEntries()-1<=(NDOFMjj+NDOFMgg+1)) noCorr=true;

  std::cout << "noCorr = " << noCorr << " , chi2/chi2dof = " << chi2/chi2dof << std::endl;

  for(int k=0; k<totalNDOF; ++k){
      //for(int k=0; k<1; ++k){

    if(noCorr == 1) continue;

    if(resMass == 0 && c == 0 && k != 0) continue; 
    if(resMass == 0 && c == 1 && k != 2) continue; 
    if(resMass == 0 && c == 2 && k != 0) continue; 
    if(resMass == 0 && c == 3 && k != 2) continue; 

    TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
    c1->Divide(3,3);
    TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
    c2->Divide(3,3);

    int iCount = 0;

    std::vector<double> pulls,pullsBkg,minNLogL, powMgg, powMjj;
    pulls.clear();
    pullsBkg.clear();
    minNLogL.clear();
    powMgg.clear();
    powMjj.clear();
    for(int i=0; i<Npse; ++i){
      
      resultsTotal[k][i] = -9999;

      TRandom rndm(0);
      float nEvents = rndm.Poisson(data->sumEntries());
      if(nEvents == 0) continue;

      /*RooProdPdf BkgTruthTmp("BkgTruthTmp","",*BkgTruthCondtmp);
      RooRealVar nbkgTruth("nbkgTruth","",nEvents);
      RooExtendPdf BkgTruth("BkgTruth","",BkgTruthTmp,nbkgTruth);
      //RooRealVar nsigTruth("nsigTruth","",0.);
      //RooAddPdf BkgTruth("BkgTruth", "", RooArgList(BkgTruthTmp,*w->pdf(TString::Format("SigPdf_cat%d",c))), RooArgList(nbkgTruth,nsigTruth));*/

      RooProdPdf BkgFitTmp("BkgFitTmp","",RooArgList(*MggBkgTmp[k],*MjjBkgTmp[k]));
      RooRealVar nbkg("nbkg","",nEvents,0.5*nEvents,1.5*nEvents);
      RooRealVar nsig("nsig","",0.,-1.*nEvents,nEvents);
      //RooRealVar nsig("nsig","",nsigTruth.getVal(),-3.*nsigTruth.getVal(),3*nsigTruth.getVal());
      RooAddPdf BkgFit(TString::Format("BkgFit_cat%d",c), "", RooArgList(BkgFitTmp,*w->pdf(TString::Format("SigPdf_cat%d",c))), RooArgList(nbkg,nsig));

      //n_gen_sr=BkgTruth.createIntegral(RooArgSet(*mGG,*mJJ),NormSet(RooArgSet(*mGG,*mJJ)),Range("sigRegion"))->getVal()*nbkgTruth.getVal();
      
      //set parameters to initial values      
      RooArgSet* paramsMgg = (RooArgSet*) MggBkgTmp[k]->getParameters(*mGG);
      TIterator* iterMgg = paramsMgg->createIterator();
      TObject* tempObjMgg=0;
      while((tempObjMgg=iterMgg->Next())){
                RooRealVar* var = (RooRealVar*)tempObjMgg;
                var->setVal(0.1);
                //std::cout << "k = " << k << " , varNameMgg = " << tempObjMgg->GetName() << " , val = " << var->getVal() << std::endl;
      }
         
      RooArgSet* paramsMjj = (RooArgSet*) MjjBkgTmp[k]->getParameters(*mJJ);
      TIterator* iterMjj = paramsMjj->createIterator();
      TObject* tempObjMjj=0;
      while((tempObjMjj=iterMjj->Next())){
                 RooRealVar* var = (RooRealVar*)tempObjMjj;
                 var->setVal(0.1);
                 //std::cout << "k = " << k << " , varNameMjj = " << tempObjMjj->GetName() << " , val = " << var->getVal() << std::endl;
      }
 
      p6->setVal(-0.1);
      p7->setVal(-0.1);
      p8->setVal(-0.1);
      p9->setVal(-0.1);


      //nEvents = nEvents + nsigTruth.getVal();
      //RooRandom::randomGenerator()->SetSeed(nEvent); 
      //RooDataSet *dataGen = BkgTruth.generate(RooArgSet(*mGG,*mJJ),nEvents,Range("massFit"),Extended(kTRUE));
      RooDataSet *dataGen = BkgTruthCond->generate(RooArgSet(*mGG,*mJJ),nEvents,Range("massFit"),Extended(kTRUE));
      RooFitResult* fitResults = BkgFit.fitTo(*dataGen,Range("massFit"),Save(),SumW2Error(kTRUE)) ; 

     /*if(k == 1){
       std::cout << "TOY = " << i << " , PowerMgg = " << p8->getVal()<< "+/-" << p8->getError() << std::endl;
       std::cout << "TOY = " << i << " , PowerMjj = " << p9->getVal()<< "+/-" << p9->getError() << std::endl;
     }*/
       
      //corrHist->Fill(genDataset->correlation(*mGG,*mJJ));
      int fitStatus = fitResults->status();
      float corr = fitResults->correlation("nsig","nbkg");
      float fitN = nsig.getVal();
      float fitNerr = nsig.getError();
      float fitBkg = nbkg.getVal();
      float fitBkgerr = nbkg.getError();
      //float pull = (nsigTruth.getVal()-fitN)/(fitNerr);
      float pull = (0-fitN)/(fitNerr);
      float pullBkg = (nEvents-fitBkg)/(fitBkgerr);
      //std::cout <<  " k = " << k << " , TOY = " << i << " , status = " << fitStatus << " , nbkg = " << fitBkg << " , nbkgErr = " << fitBkgerr << " , nsig = " << fitN << " , nsigErr = " << fitNerr << " , corr = " << corr << " , pullSig = "   << pull << " , pullBkg = "   << pullBkg << " , -LogLike = " << fitResults->minNll() << " , nEvents = " << nEvents << " " << data->sumEntries() << std::endl;
      if(fitNerr == 0. || fitBkgerr== 0.) continue;
      if(fabs(fitResults->minNll()) > 10e6) continue;
      if(fabs(nsig.getVal())> 0.85*nEvents) continue;
      //if(fabs(nsig.getVal())> 0.85*3*nsigTruth.getVal()) continue;
      if(nbkg.getVal() < 0.5*1.15*nEvents || nbkg.getVal() > 1.5*0.85*nEvents) continue;
      if(fabs(corr)>0.85) continue;
      iCount++;
      pulls.push_back(pull);
      pullsBkg.push_back(pullBkg);
      resultsTotal[k][i]=pull;
      minNLogL.push_back(fitResults->minNll());
      if(k == 1)powMgg.push_back(p8->getVal());
      if(k == 1)powMjj.push_back(p9->getVal());
      if(iCount > 0 && iCount <= 9){
      c1->cd(iCount);
      RooPlot *frame = mGG->frame();
      if(iCount==1){
	data->plotOn(frame);
	MggBkgTruth->plotOn(frame,LineColor(kBlue),Range("fitrange"),NormRange("fitrange"));
	frame->SetTitle("");
	frame->SetMinimum(1e-5);
	frame->SetMaximum(1.40*frame->GetMaximum());
	frame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
	frame->Draw();

	TLegend *legmc = new TLegend(0.55,0.75,0.98,0.9);
	legmc->AddEntry(frame->getObject(0),"Data CS","LPE");
	legmc->AddEntry(frame->getObject(1),TString::Format("%.5s Truth",MggBkgTruth->GetName()),"L");
	legmc->SetBorderSize(0);
	legmc->SetFillStyle(0);
	legmc->Draw();    
	TLatex *lat  = new TLatex(minMggMassFit+3.0,0.85*frame->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
	lat->Draw();
	TLatex *lat2 = new TLatex(minMggMassFit+3.0,0.7*frame->GetMaximum(),catdesc.at(c));
	lat2->Draw();

      }else{
      dataGen->plotOn(frame);
      RooRealVar *q0 = new RooRealVar("q0","",p0->getVal());
      RooRealVar *q1 = new RooRealVar("q1","",p1->getVal());
      RooRealVar *q2 = new RooRealVar("q2","",p2->getVal());
      RooRealVar *q3 = new RooRealVar("q3","",p3->getVal());
      RooRealVar *q4 = new RooRealVar("q4","",p4->getVal());
      RooRealVar *q5 = new RooRealVar("q5","",p5->getVal());
      RooRealVar *q6 = new RooRealVar("q6","",p6->getVal());
      RooRealVar *q7 = new RooRealVar("q7","",p7->getVal());
      RooRealVar *q8 = new RooRealVar("q8","",p8->getVal());
      RooRealVar *q9 = new RooRealVar("q9","",p9->getVal());
      RooExponential* expo1Mgg_tmp = new RooExponential("exp1FitMgg_tmp","",*mGG,*q6);
      RooExponential* expo2Mgg_tmp = new RooExponential("exp2FitMgg_tmp","",*mGG,*q7);
      RooExponential* expo1Mjj_tmp = new RooExponential("exp1FitMjj_tmp","",*mJJ,*q8);
      RooExponential* expo2Mjj_tmp = new RooExponential("exp2FitMjj_tmp","",*mJJ,*q9);
      RooGenericPdf* pow1Mgg_tmp = new RooGenericPdf("pow1FitMgg_tmp","pow(@0,@1)",RooArgList(*mGG,*q6));
      RooGenericPdf* pow2Mgg_tmp = new RooGenericPdf("pow2FitMgg_tmp","pow(@0,@1)",RooArgList(*mGG,*q7));
      RooGenericPdf* pow1Mjj_tmp = new RooGenericPdf("pow1FitMjj_tmp","pow(@0,@1)",RooArgList(*mJJ,*q8));
      RooGenericPdf* pow2Mjj_tmp = new RooGenericPdf("pow2FitMjj_tmp","pow(@0,@1)",RooArgList(*mJJ,*q9));
      RooAbsPdf *fitFuncMgg, *fitFuncMjj;
      switch(k){
      case 0: fitFuncMgg = new RooExponential("fitFuncMgg","",*mGG,*q6); break;
      case 1: fitFuncMgg = new RooAddPdf("fitFuncMgg","",RooArgList(*expo1Mgg_tmp,*expo2Mgg_tmp),RooArgList(*q6,*q7)); break;
      case 2: fitFuncMgg = new RooGenericPdf("fitFuncMgg","pow(@0,@1)",RooArgList(*mGG,*q6)); break;
      case 3: fitFuncMgg = new RooAddPdf("fitFuncMgg","",RooArgList(*pow1Mgg_tmp,*pow2Mgg_tmp),RooArgList(*q6,*q7)); break;
      case 4: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1)); break;
      case 5: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1)); break;
      case 6: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1,*q2)); break;
      }
      switch(k){
      case 0: fitFuncMjj = new RooExponential("fitFuncMjj","",*mJJ,*q8); break;
      case 1: fitFuncMjj = new RooAddPdf("fitFuncMjj","",RooArgList(*expo1Mjj_tmp,*expo2Mjj_tmp),RooArgList(*q8,*q9)); break;
      case 2: fitFuncMjj = new RooGenericPdf("fitFuncMjj","pow(@0,@1)",RooArgList(*mJJ,*q8)); break;
      case 3: fitFuncMjj = new RooAddPdf("fitFuncMjj","",RooArgList(*pow1Mjj_tmp,*pow2Mjj_tmp),RooArgList(*q8,*q9)); break;
      case 4: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q3,*q4)); break;
      case 5: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q3,*q4,*q5)); break;
      case 6: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q3,*q4)); break;
      }
      RooProdPdf *fitFunc= new RooProdPdf("fitFunc","",RooArgSet(*fitFuncMgg,*fitFuncMjj));
      fitFunc->plotOn(frame,LineColor(kGreen),Range("fitrange"),NormRange("fitrange"));
      frame->SetTitle("");                                                                             
      frame->SetMinimum(1e-5);                                                                          
      frame->SetMaximum(1.40*frame->GetMaximum());                                                     
      frame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");                                           
      frame->Draw();

      TLegend *legmc = new TLegend(0.52,0.75,0.99,0.9);
      legmc->AddEntry(frame->getObject(0),TString::Format("Gen PSE #%i",i-1),"LPE");
      legmc->AddEntry(frame->getObject(1),TString::Format("%.5s Candidate",MggBkgTmp[k]->GetName()),"L");
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();    
      TLatex *lat  = new TLatex(minMggMassFit+3.0,0.85*frame->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
      lat->Draw();
      TLatex *lat2 = new TLatex(minMggMassFit+3.0,0.7*frame->GetMaximum(),catdesc.at(c));
      lat2->Draw();
      //TLatex *lat3 = new TLatex();
      //lat3->SetNDC();
      //lat3->DrawLatex(0.62,0.68,TString::Format("#chi^{2}/N = %.2f",mcs->fitParams(i)->getRealValue("chi2")/mcs->fitParams(i)->getRealValue("ndof")));
      //lat3->DrawLatex(0.62,0.62,TString::Format("p(#chi^{2},N) = %.2f",mcs->fitParams(i)->getRealValue("prob")));
      }
      c2->cd(iCount);
      frame = mJJ->frame();
      if(iCount==1){
	data->plotOn(frame);
	MjjBkgTruth->plotOn(frame,LineColor(kBlue),Range("fitrange"),NormRange("fitrange"));
	frame->SetTitle("");
	frame->SetMinimum(1e-5);
	frame->SetMaximum(1.40*frame->GetMaximum());
	frame->GetXaxis()->SetTitle("m_{jj} (GeV)");
	frame->Draw();

	TLegend *legmc = new TLegend(0.55,0.75,0.98,0.9);
	legmc->AddEntry(frame->getObject(0),"Data CS","LPE");
	legmc->AddEntry(frame->getObject(1),TString::Format("%.5s Truth",MjjBkgTruth->GetName()),"L");
	legmc->SetBorderSize(0);
	legmc->SetFillStyle(0);
	legmc->Draw();    
	TLatex *lat  = new TLatex(minMjjMassFit+3.0,0.85*frame->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
	lat->Draw();
	TLatex *lat2 = new TLatex(minMjjMassFit+3.0,0.7*frame->GetMaximum(),catdesc.at(c));
	lat2->Draw();

      }else{
      dataGen->plotOn(frame);
      RooRealVar *q0 = new RooRealVar("q0","",p0->getVal());
      RooRealVar *q1 = new RooRealVar("q1","",p1->getVal());
      RooRealVar *q2 = new RooRealVar("q2","",p2->getVal());
      RooRealVar *q3 = new RooRealVar("q3","",p3->getVal());
      RooRealVar *q4 = new RooRealVar("q4","",p4->getVal());
      RooRealVar *q5 = new RooRealVar("q5","",p5->getVal());
      RooRealVar *q6 = new RooRealVar("q6","",p6->getVal());
      RooRealVar *q7 = new RooRealVar("q7","",p7->getVal());
      RooRealVar *q8 = new RooRealVar("q8","",p8->getVal());
      RooRealVar *q9 = new RooRealVar("q9","",p9->getVal());
      RooExponential* expo1Mgg_tmp = new RooExponential("exp1FitMgg_tmp","",*mGG,*q6);
      RooExponential* expo2Mgg_tmp = new RooExponential("exp2FitMgg_tmp","",*mGG,*q7);
      RooExponential* expo1Mjj_tmp = new RooExponential("exp1FitMjj_tmp","",*mJJ,*q8);
      RooExponential* expo2Mjj_tmp = new RooExponential("exp2FitMjj_tmp","",*mJJ,*q9);
      RooGenericPdf* pow1Mgg_tmp = new RooGenericPdf("pow1FitMgg_tmp","pow(@0,@1)",RooArgList(*mGG,*q6));
      RooGenericPdf* pow2Mgg_tmp = new RooGenericPdf("pow2FitMgg_tmp","pow(@0,@1)",RooArgList(*mGG,*q7));
      RooGenericPdf* pow1Mjj_tmp = new RooGenericPdf("pow1FitMjj_tmp","pow(@0,@1)",RooArgList(*mJJ,*q8));
      RooGenericPdf* pow2Mjj_tmp = new RooGenericPdf("pow2FitMjj_tmp","pow(@0,@1)",RooArgList(*mJJ,*q9));
      RooAbsPdf *fitFuncMgg, *fitFuncMjj;
      switch(k){
      case 0: fitFuncMgg = new RooExponential("fitFuncMgg","",*mGG,*q6); break;
      case 1: fitFuncMgg = new RooAddPdf("fitFuncMgg","",RooArgList(*expo1Mgg_tmp,*expo2Mgg_tmp),RooArgList(*q6,*q7)); break;
      case 2: fitFuncMgg = new RooGenericPdf("fitFuncMgg","pow(@0,@1)",RooArgList(*mGG,*q6)); break;
      case 3: fitFuncMgg = new RooAddPdf("fitFuncMgg","",RooArgList(*pow1Mgg_tmp,*pow2Mgg_tmp),RooArgList(*q6,*q7)); break;
      case 4: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1)); break;
      case 5: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1)); break;
      case 6: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1,*q2)); break;
      }
      switch(k){
      case 0: fitFuncMjj = new RooExponential("fitFuncMjj","",*mJJ,*q8); break;
      case 1: fitFuncMjj = new RooAddPdf("fitFuncMjj","",RooArgList(*expo1Mjj_tmp,*expo2Mjj_tmp),RooArgList(*q8,*q9)); break;
      case 2: fitFuncMjj = new RooGenericPdf("fitFuncMjj","pow(@0,@1)",RooArgList(*mJJ,*q8)); break;
      case 3: fitFuncMjj = new RooAddPdf("fitFuncMjj","",RooArgList(*pow1Mjj_tmp,*pow2Mjj_tmp),RooArgList(*q8,*q9)); break;
      case 4: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q3,*q4)); break;
      case 5: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q3,*q4,*q5)); break;
      case 6: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q3,*q4)); break;
      }
      RooProdPdf*fitFunc= new RooProdPdf("fitFunc","",RooArgSet(*fitFuncMgg,*fitFuncMjj));
      fitFunc->plotOn(frame,LineColor(kGreen),Range("fitrange"),NormRange("fitrange"));
      frame->SetTitle("");                                                                             
      frame->SetMinimum(1e-5);                                                                          
      frame->SetMaximum(1.40*frame->GetMaximum());                                                     
      frame->GetXaxis()->SetTitle("m_{jj} (GeV)");                                           
      frame->Draw();

      TLegend* legmc = new TLegend(0.52,0.75,0.99,0.9);
      legmc->AddEntry(frame->getObject(0),TString::Format("Gen PSE #%i",i-1),"LPE");
      legmc->AddEntry(frame->getObject(1),TString::Format("%.5s Candidate",MjjBkgTmp[k]->GetName()),"L");
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();    
      TLatex *lat  = new TLatex(minMjjMassFit+3.0,0.85*frame->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
      lat->Draw();
      TLatex *lat2 = new TLatex(minMjjMassFit+3.0,0.7*frame->GetMaximum(),catdesc.at(c));
      lat2->Draw();
      //TLatex *lat3 = new TLatex();
      //lat3->SetNDC();
      //lat3->DrawLatex(0.62,0.68,TString::Format("#chi^{2}/N = %.2f",mcs->fitParams(i)->getRealValue("chi2")/mcs->fitParams(i)->getRealValue("ndof")));
      //lat3->DrawLatex(0.62,0.62,TString::Format("p(#chi^{2},N) = %.2f",mcs->fitParams(i)->getRealValue("prob")));
      }
      }

    }

    float h1_low=0, h1_high=0;
    if (c==0 || c==2){
      h1_low=-1;
      h1_high=1;
    }
    else{
      h1_low=-2;
      h1_high=2;
    }
    
    for(int i=0; i<pulls.size(); ++i){
      for(int j=i; j<pulls.size(); ++j){
	if(pulls[i]>pulls[j]){
	  float tmp=pulls[i];
	  pulls[i]=pulls[j];
	  pulls[j]=tmp;
	}}}
    if(pulls.size()==0)
      results[k] = -9999;
    else if(pulls.size()%2==0)
      results[k] = 0.5*(pulls[pulls.size()/2]+pulls[pulls.size()/2-1]);
    else
      results[k] = pulls[pulls.size()/2];

    if(pulls.size()==0){ pulls.push_back(-9999);  pulls.push_back(-9998);}

    TH1F *h1 = new TH1F("h1","",50,-5,5);
    h1->GetXaxis()->SetTitle("pull(nSig)");
    for(int i=0; i<pulls.size(); i++) h1->Fill(pulls[i]);
    TCanvas *c0 = new TCanvas("c0","c0",700,500);
    h1->SetMaximum(h1->GetMaximum()*1.50);
    h1->Draw();
    if(h1->Integral()>4){
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(h1_low+0.15,0.91*h1->GetMaximum(),TString::Format("Num PSE: %d",pulls.size()) );
      lat13->DrawLatex(h1_low+0.15,0.85*h1->GetMaximum(),TString::Format("Gen function: %.5s,%.5s",MggBkgTruth->GetName(),MjjBkgTruth->GetName()) );
      lat13->DrawLatex(h1_low+0.15,0.79*h1->GetMaximum(),TString::Format("Fit function: %.5s,%.5s",MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName()) );
      lat13->DrawLatex(h1_low+0.15,0.73*h1->GetMaximum(),TString::Format("Median = %.2f ",results[k]) );
      lat13->DrawLatex(h1_low+0.15,0.67*h1->GetMaximum(),TString::Format("N_{gen}^{SR} = %.2f ",n_gen_sr) );
    }
    else{
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(0.0,0.5,"Not enough entries." );
      lat13->DrawLatex(0.0,0.4,TString::Format("Median = %.2f ",results[k])  );
    }

    for(int i=0; i<pullsBkg.size(); ++i){
      for(int j=i; j<pullsBkg.size(); ++j){
	if(pullsBkg[i]>pullsBkg[j]){
	  float tmp=pullsBkg[i];
	  pullsBkg[i]=pullsBkg[j];
	  pullsBkg[j]=tmp;
	}}}
    if(pullsBkg.size()==0)
      resultsBkg[k] = -9999;
    else if(pullsBkg.size()%2==0)
      resultsBkg[k] = 0.5*(pullsBkg[pullsBkg.size()/2]+pullsBkg[pullsBkg.size()/2-1]);
    else
      resultsBkg[k] = pullsBkg[pullsBkg.size()/2];

    if(pullsBkg.size()==0){ pullsBkg.push_back(-9999);  pullsBkg.push_back(-9998);}

    TH1F *h2 = new TH1F("h2","",50,-5,5);
    h2->GetXaxis()->SetTitle("pull(nBkg)");
    for(int i=0; i<pullsBkg.size(); i++) h2->Fill(pullsBkg[i]);
    TCanvas *c3 = new TCanvas("c3","c3",700,500);
    h2->SetMaximum(h2->GetMaximum()*1.50);
    h2->Draw();
    if(h2->Integral()>4){
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(h1_low+0.15,0.91*h2->GetMaximum(),TString::Format("Num PSE: %d",pullsBkg.size()) );
      lat13->DrawLatex(h1_low+0.15,0.85*h2->GetMaximum(),TString::Format("Gen function: %.5s,%.5s",MggBkgTruth->GetName(),MjjBkgTruth->GetName()) );
      lat13->DrawLatex(h1_low+0.15,0.79*h2->GetMaximum(),TString::Format("Fit function: %.5s,%.5s",MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName()) );
      lat13->DrawLatex(h1_low+0.15,0.73*h2->GetMaximum(),TString::Format("Median = %.2f ",resultsBkg[k]) );
      lat13->DrawLatex(h1_low+0.15,0.67*h2->GetMaximum(),TString::Format("N_{gen}^{SR} = %.2f ",n_gen_sr) );
    }
    else{
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(0.0,0.5,"Not enough entries." );
      lat13->DrawLatex(0.0,0.4,TString::Format("Median = %.2f ",resultsBkg[k])  );
    }

    TH1F *h3 = new TH1F("h3","",40,0.,400.);
    h3->GetXaxis()->SetTitle("-minLogLikelihood");
    for(int i=0; i<minNLogL.size(); i++) h3->Fill(minNLogL[i]);
    TCanvas *c4 = new TCanvas("c4","c4",700,500);
    h3->SetMaximum(h3->GetMaximum()*1.50);
    h3->Draw();
    if(h3->Integral()>4){
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(h1_low+0.15,0.91*h3->GetMaximum(),TString::Format("Num PSE: %d",pullsBkg.size()) );
      lat13->DrawLatex(h1_low+0.15,0.85*h3->GetMaximum(),TString::Format("Gen function: %.5s,%.5s",MggBkgTruth->GetName(),MjjBkgTruth->GetName()) );
      lat13->DrawLatex(h1_low+0.15,0.79*h3->GetMaximum(),TString::Format("Fit function: %.5s,%.5s",MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName()) );
      lat13->DrawLatex(h1_low+0.15,0.67*h3->GetMaximum(),TString::Format("N_{gen}^{SR} = %.2f ",n_gen_sr) );
    }
    else{
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(0.0,0.5,"Not enough entries." );
      //lat13->DrawLatex(0.0,0.4,TString::Format("Median = %.2f ",results[k])  );
    }
    
    TH1F *h4 = new TH1F("h4","",100,-10,0.);
    h4->GetXaxis()->SetTitle("powMgg");
    for(int i=0; i<powMgg.size(); i++) h4->Fill(powMgg.at(i));
    TCanvas *c5 = new TCanvas("c5","c5",700,500);
    h4->SetMaximum(h4->GetMaximum()*1.50);
    h4->Draw();
    if(h4->Integral()>4){
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(h1_low+0.15,0.91*h3->GetMaximum(),TString::Format("Num PSE: %d",powMgg.size()) );
      lat13->DrawLatex(h1_low+0.15,0.85*h3->GetMaximum(),TString::Format("Gen function: %.5s,%.5s",MggBkgTruth->GetName(),MjjBkgTruth->GetName()) );
      lat13->DrawLatex(h1_low+0.15,0.79*h3->GetMaximum(),TString::Format("Fit function: %.5s,%.5s",MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName()) );
      lat13->DrawLatex(h1_low+0.15,0.67*h3->GetMaximum(),TString::Format("N_{gen}^{SR} = %.2f ",n_gen_sr) );
    }
    else{
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(0.0,0.5,"Not enough entries." );
      //lat13->DrawLatex(0.0,0.4,TString::Format("Median = %.2f ",results[k])  );
    }

    TH1F *h5 = new TH1F("h5","",100,-10.,0.);
    h5->GetXaxis()->SetTitle("powMjj");
    for(int i=0; i<powMjj.size(); i++) h5->Fill(powMjj.at(i));
    TCanvas *c6 = new TCanvas("c6","c6",700,500);
    h5->SetMaximum(h5->GetMaximum()*1.50);
    h5->Draw();
    if(h5->Integral()>4){
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(h1_low+0.15,0.91*h3->GetMaximum(),TString::Format("Num PSE: %d",powMjj.size()) );
      lat13->DrawLatex(h1_low+0.15,0.85*h3->GetMaximum(),TString::Format("Gen function: %.5s,%.5s",MggBkgTruth->GetName(),MjjBkgTruth->GetName()) );
      lat13->DrawLatex(h1_low+0.15,0.79*h3->GetMaximum(),TString::Format("Fit function: %.5s,%.5s",MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName()) );
      lat13->DrawLatex(h1_low+0.15,0.67*h3->GetMaximum(),TString::Format("N_{gen}^{SR} = %.2f ",n_gen_sr) );
    }
    else{
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(0.0,0.5,"Not enough entries." );
      //lat13->DrawLatex(0.0,0.4,TString::Format("Median = %.2f ",results[k])  );
    }

    c0->SaveAs(TString::Format("plots/pulls2D_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.png",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));  
    c0->SaveAs(TString::Format("plots/pulls2D_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.pdf",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));  
    
    c3->SaveAs(TString::Format("plots/pulls2D_bkg_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.png",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));  
    c3->SaveAs(TString::Format("plots/pulls2D_bkg_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.pdf",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));

    c4->SaveAs(TString::Format("plots/minNLL_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.png",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));  
    c4->SaveAs(TString::Format("plots/minNLL_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.pdf",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));  

    if(k==1){
    c5->SaveAs(TString::Format("plots/powMgg_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.png",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));  
    c5->SaveAs(TString::Format("plots/powMgg_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.pdf",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));  

    c6->SaveAs(TString::Format("plots/powMjj_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.png",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));  
    c6->SaveAs(TString::Format("plots/powMjj_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.pdf",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));  }

    c1->SaveAs(TString::Format("plots/toymc2D_Mgg_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.png",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));
    c1->SaveAs(TString::Format("plots/toymc2D_Mgg_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.pdf",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));

    c2->SaveAs(TString::Format("plots/toymc2D_Mjj_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.png",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));
    c2->SaveAs(TString::Format("plots/toymc2D_Mjj_gen%.5s%.5s_fit%.5s%.5s_cat%d_mass%d.pdf",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),MggBkgTmp[k]->GetName(),MjjBkgTmp[k]->GetName(),c,resMass));

    delete c0;
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete c5;
    delete c6;
    delete h1;
    delete h2;
    delete h3;
    delete h4;
    delete h5;
  }

  fprintf(fout,"resMass=%d,cat%d,withCorr=%d\t%.5s,%.5s\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n",resMass,c,withCorr,MggBkgTruth->GetName(),MjjBkgTruth->GetName(),results[0],results[1],results[2],results[3],results[4],results[5],results[6]);
  fprintf(foutCorr,"resMass=%d,cat%d,withCorr=%d\t%.5s,%.5s\t%.4f +/- %.4f\n",resMass,c,withCorr,MggBkgTruth->GetName(),MjjBkgTruth->GetName(),corrHist->GetMean(),effSigma(corrHist));
  
  for(int i = 0; i<Npse;i++)
    fprintf(foutTotal,"resMass=%d,cat%d,withCorr=%d\t%.5s,%.5s\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n",resMass,c,withCorr,MggBkgTruth->GetName(),MjjBkgTruth->GetName(),resultsTotal[0][i],resultsTotal[1][i],resultsTotal[2][i],resultsTotal[3][i],resultsTotal[4][i],resultsTotal[5][i],resultsTotal[6][i]);

  return;
}


void SetConstantParams(const RooArgSet* params) {

  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  

}


Double_t effSigma(TH1 *hist) {

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }

  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

int main(int argc, const char* argv[]){

  int cat=0, modelNum1=0, modelNum2=0, searchMass=0;

  if (argc>1)
    cat=atoi(argv[1]);
  if (argc>2)
    modelNum1=atoi(argv[2]);
  if (argc>3)
    modelNum2=atoi(argv[3]);
  if (argc>4)
    searchMass=atoi(argv[4]);
  if (argc>5)
    withCorr=atoi(argv[5]);

  resMass=searchMass;

  if (searchMass==0){
    inDir  = "/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_nonresSearch_withKinFit/";
    NCAT=4;
  }
  else{
    NCAT = 2;
    inDir  = "/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v44/v44_fitTo2D_resSearch_withRegKinFit/";
  }

  runfits(cat,modelNum1,modelNum2);

  return 0;
}

void style(){
  TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
  defaultStyle->SetOptStat(0000);
  defaultStyle->SetOptFit(000);
  defaultStyle->SetPalette(1);
  /////// pad ////////////
  defaultStyle->SetPadBorderMode(1);
  defaultStyle->SetPadBorderSize(1);
  defaultStyle->SetPadColor(0);
  defaultStyle->SetPadTopMargin(0.05);
  defaultStyle->SetPadBottomMargin(0.13);
  defaultStyle->SetPadLeftMargin(0.13);
  defaultStyle->SetPadRightMargin(0.02);
  /////// canvas /////////
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetCanvasColor(0);
  defaultStyle->SetCanvasDefH(600);
  defaultStyle->SetCanvasDefW(600);
  /////// frame //////////
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetFrameBorderSize(1);
  defaultStyle->SetFrameFillColor(0);
  defaultStyle->SetFrameLineColor(1);
  /////// label //////////
  defaultStyle->SetLabelOffset(0.005,"XY");
  defaultStyle->SetLabelSize(0.05,"XY");
  defaultStyle->SetLabelFont(42,"XY");
  /////// title //////////
  defaultStyle->SetTitleOffset(1.1,"X");
  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.25,"Y");
  defaultStyle->SetTitleSize(0.05,"Y");
  defaultStyle->SetTitleFont(42, "XYZ");
  /////// various ////////
  defaultStyle->SetNdivisions(505,"Y");
  defaultStyle->SetLegendBorderSize(0); // For the axis titles:

  defaultStyle->SetTitleColor(1, "XYZ");
  defaultStyle->SetTitleFont(42, "XYZ");
  defaultStyle->SetTitleSize(0.06, "XYZ");
 
  // defaultStyle->SetTitleYSize(Float_t size = 0.02);
  defaultStyle->SetTitleXOffset(0.9);
  defaultStyle->SetTitleYOffset(1.05);
  // defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:
  defaultStyle->SetLabelColor(1, "XYZ");
  defaultStyle->SetLabelFont(42, "XYZ");
  defaultStyle->SetLabelOffset(0.007, "XYZ");
  defaultStyle->SetLabelSize(0.04, "XYZ");

  // For the axis:
  defaultStyle->SetAxisColor(1, "XYZ");
  defaultStyle->SetStripDecimals(kTRUE);
  defaultStyle->SetTickLength(0.03, "XYZ");
  defaultStyle->SetNdivisions(510, "XYZ");
  defaultStyle->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  defaultStyle->SetPadTickY(1);
  defaultStyle->cd();
  return;
}
