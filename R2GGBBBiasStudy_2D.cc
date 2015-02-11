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

Int_t NCAT = 2;
TString inDir   = "/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v42/v42_fitTo2D_nonresSearch_withKinFit/";

void AddBkgData(RooWorkspace*, int);
RooAbsPdf* BkgMggModelFit(RooWorkspace*, int, int);
RooAbsPdf* BkgMjjModelFit(RooWorkspace*, int, int);
void BkgModelBias(RooWorkspace*,int,RooAbsPdf*,RooAbsPdf*,FILE*);
void SetParamNames(RooWorkspace*);
void SetConstantParams(const RooArgSet* params);
Double_t effSigma(TH1 *hist);
 
RooArgSet* defineVariables()
{
  // define variables of the input ntuple
  RooRealVar* mJJ  = new RooRealVar("mjj","M(jj)",60,180,"GeV");
  RooRealVar* mGG  = new RooRealVar("mgg","M(#gamma#gamma)",100,180,"GeV");
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
  const int nTruthModels=5;
  RooAbsPdf *MggBkgTruth=0;
  RooAbsPdf *MjjBkgTruth=0;

  TString card_name("hgghbb_models_Pol_8TeV.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  AddBkgData(w,cat);

  FILE *fout = fopen("resultsBias2D.txt","a");
  if(modelNumMgg==0 && modelNumMjj ==0) fprintf(fout,"%s\n\n",inDir.Data());

  if(modelNumMgg==2 || modelNumMjj==2) return;//skip Landau, it sucks.
  MggBkgTruth = BkgMggModelFit(w,cat,modelNumMgg); //Ber, Exp, Lan, Lau, Pow
  MjjBkgTruth = BkgMjjModelFit(w,cat,modelNumMjj); //Ber, Exp, Lan, Lau, Pow
  BkgModelBias(w,cat,MggBkgTruth,MjjBkgTruth,fout);


  if(modelNumMgg==4 && modelNumMjj==4) fprintf(fout,"\n\n");
  fclose(fout);
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
  RooRealVar weightVar("evWeight","",1,0,1000);

//****************************//
// CMS Data Set
//****************************//
// retrieve the data tree;
// no common preselection cut applied yet; 

  TFile dataFile(inDir+"DataCS_m0.root");   
  TTree* dataTree     = (TTree*) dataFile.Get("TCVARS");
  weightVar.setVal(1.);
  ntplVars->add(RooArgList(weightVar));

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

  RooRealVar* mGG     = w->var("mgg");
  mGG->setUnit("GeV");
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);

  data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));

  RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("p1mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
  RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("p2mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c)));
  RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("p3mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c)));
  RooFormulaVar *p4mod = new RooFormulaVar(TString::Format("p4mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope4_cat%d",c)));
  RooFormulaVar *p5mod = new RooFormulaVar(TString::Format("p5mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope5_cat%d",c)));
  RooFormulaVar *p1arg = new RooFormulaVar(TString::Format("p1arg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_arg1_cat%d",c)));
  RooFormulaVar *p2arg = new RooFormulaVar(TString::Format("p2arg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_arg2_cat%d",c)));
  RooFormulaVar *p3arg = new RooFormulaVar(TString::Format("p3arg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_arg3_cat%d",c)));
  RooFormulaVar *p4arg = new RooFormulaVar(TString::Format("p4arg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_arg4_cat%d",c)));
  RooFormulaVar *p5arg = new RooFormulaVar(TString::Format("p5arg_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_arg5_cat%d",c)));
  RooFormulaVar *p1wid = new RooFormulaVar(TString::Format("p1wid_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_wid1_cat%d",c)));
  RooFormulaVar *p2wid = new RooFormulaVar(TString::Format("p2wid_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_wid2_cat%d",c)));
  RooFormulaVar *p3wid = new RooFormulaVar(TString::Format("p3wid_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_wid3_cat%d",c)));
  RooFormulaVar *p4wid = new RooFormulaVar(TString::Format("p4wid_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_wid4_cat%d",c)));
  RooFormulaVar *p5wid = new RooFormulaVar(TString::Format("p5wid_cat%d",c),"","@0",*w->var(TString::Format("mgg_bkg_8TeV_wid5_cat%d",c)));

  RooExponential* expo1 = new RooExponential("exp1","",*mGG,*p1arg);
  RooExponential* expo2 = new RooExponential("exp2","",*mGG,*p2arg);
  RooExponential* expo3 = new RooExponential("exp3","",*mGG,*p3arg);
  RooExponential* expo4 = new RooExponential("exp4","",*mGG,*p4arg);
  RooExponential* expo5 = new RooExponential("exp5","",*mGG,*p5arg);
  RooLandau* lan1 = new RooLandau("lan1","",*mGG,*p1arg,*p1wid);
  RooLandau* lan2 = new RooLandau("lan2","",*mGG,*p2arg,*p2wid);
  RooLandau* lan3 = new RooLandau("lan3","",*mGG,*p3arg,*p3wid);
  RooLandau* lan4 = new RooLandau("lan4","",*mGG,*p4arg,*p4wid);
  RooLandau* lan5 = new RooLandau("lan5","",*mGG,*p5arg,*p5wid);

  const int totalNDOF=5;
  int NDOF[totalNDOF]={0};
  float minNLL[totalNDOF]={0.0}, chi2prob[totalNDOF]={-1.0};
  int bestN = 0;

  RooAbsPdf* MggBkgTmp[totalNDOF] = {0};
  char fitName[10];

  switch (modelNum){

  case 0: //Bernstein
    MggBkgTmp[0] = new RooBernstein("BerN0Mgg", "", *mGG, RooArgList(*p1mod));
    MggBkgTmp[1] = new RooBernstein("BerN1Mgg", "", *mGG, RooArgList(*p1mod,*p2mod));
    MggBkgTmp[2] = new RooBernstein("BerN2Mgg", "", *mGG, RooArgList(*p1mod,*p2mod,*p3mod));
    MggBkgTmp[3] = new RooBernstein("BerN3Mgg", "", *mGG, RooArgList(*p1mod,*p2mod,*p3mod,*p4mod));
    MggBkgTmp[4] = new RooBernstein("BerN4Mgg", "", *mGG, RooArgList(*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Ber");
    break;

  case 1: //Exponential
    w->factory(TString::Format("mgg_bkg_8TeV_norm_cat%d[800.0,0.0,100000]",c));
    MggBkgTmp[0] = new RooExtendPdf("ExpN1Mgg","",*expo1,*w->var(TString::Format("mgg_bkg_8TeV_norm_cat%d",c)));
    MggBkgTmp[1] = new RooAddPdf("ExpN2Mgg", "", RooArgList(*expo1,*expo2), RooArgList(*p1mod,*p2mod));
    MggBkgTmp[2] = new RooAddPdf("ExpN3Mgg", "", RooArgList(*expo1,*expo2,*expo3), RooArgList(*p1mod,*p2mod,*p3mod));
    MggBkgTmp[3] = new RooAddPdf("ExpN4Mgg", "", RooArgList(*expo1,*expo2,*expo3,*expo4), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod));
    MggBkgTmp[4] = new RooAddPdf("ExpN5Mgg", "", RooArgList(*expo1,*expo2,*expo3,*expo4,*expo5), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Exp");
    break;

  case 2: //Landau
    w->factory(TString::Format("mgg_bkg_8TeV_norm_cat%d[800.0,0.0,100000]",c));
    MggBkgTmp[0] = new RooExtendPdf("LanN1Mgg","",*lan1,*w->var(TString::Format("mgg_bkg_8TeV_norm_cat%d",c)));
    MggBkgTmp[1] = new RooAddPdf("LanN2Mgg", "", RooArgList(*lan1,*lan2), RooArgList(*p1mod,*p2mod));
    MggBkgTmp[2] = new RooAddPdf("LanN3Mgg", "", RooArgList(*lan1,*lan2,*lan3), RooArgList(*p1mod,*p2mod,*p3mod));
    MggBkgTmp[3] = new RooAddPdf("LanN4Mgg", "", RooArgList(*lan1,*lan2,*lan3,*lan4), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod));
    MggBkgTmp[4] = new RooAddPdf("LanN5Mgg", "", RooArgList(*lan1,*lan2,*lan3,*lan4,*lan5), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Lan");
    break;

  case 3: //Laurent
    MggBkgTmp[0] = new RooGenericPdf("LauN1Mgg","@1*pow(@0,-4)",RooArgList(*mGG,*p1mod));
    MggBkgTmp[1] = new RooGenericPdf("LauN2Mgg","@1*pow(@0,-4)+@2*pow(@0,-3)",RooArgList(*mGG,*p1mod,*p2mod));
    MggBkgTmp[2] = new RooGenericPdf("LauN3Mgg","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)",RooArgList(*mGG,*p1mod,*p2mod,*p3mod));
    MggBkgTmp[3] = new RooGenericPdf("LauN4Mgg","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)+@4*pow(@0,-2)",RooArgList(*mGG,*p1mod,*p2mod,*p3mod,*p4mod));
    MggBkgTmp[4] = new RooGenericPdf("LauN5Mgg","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)+@4*pow(@0,-2)+@5*pow(@0,-6)",RooArgList(*mGG,*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Lau");
    break;

  case 4: //Power
    MggBkgTmp[0] = new RooGenericPdf("PowN1Mgg","@1*pow(@0,@2)",RooArgList(*mGG,*p1mod,*p1arg));
    MggBkgTmp[1] = new RooGenericPdf("PowN2Mgg","@1*pow(@0,@2)+@3*pow(@0,@4)",RooArgList(*mGG,*p1mod,*p1arg,*p2mod,*p2arg));
    MggBkgTmp[2] = new RooGenericPdf("PowN3Mgg","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)",RooArgList(*mGG,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg));
    MggBkgTmp[3] = new RooGenericPdf("PowN4Mgg","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)",RooArgList(*mGG,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg,*p4mod,*p4arg));
    MggBkgTmp[4] = new RooGenericPdf("PowN5Mgg","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)",RooArgList(*mGG,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg,*p4mod,*p4arg));
    //MggBkgTmp[4] = new RooGenericPdf(TString::Format("MggBkg_cat%d",c),"@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)+@9*pow(@0,@10)",RooArgList(*mGG,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg,*p4mod,*p4arg,*p5mod,*p5arg));
    sprintf(fitName,"Pow");
    break;
  }


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
    plotMggBkg[c]->SetMinimum(0.0);
    plotMggBkg[c]->SetMaximum(1.40*plotMggBkg[c]->GetMaximum());
    plotMggBkg[c]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    plotMggBkg[c]->Draw();  

    TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
    legmc->AddEntry(plotMggBkg[c]->getObject(2),"Weighted CS","LPE");
    legmc->AddEntry(plotMggBkg[c]->getObject(1),"Bkg Model","L");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    

    TLatex *lat  = new TLatex(minMassFit+3.0,0.85*plotMggBkg[c]->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
    lat->Draw();
    TLatex *lat2 = new TLatex(minMassFit+3.0,0.7*plotMggBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
 
    ctmp->SaveAs(TString::Format("plots/dataBkg2DMgg_%sN%i_cat%d.png",fitName,i+(modelNum!=0),c));

    if(i>0){
      float chi2 = 2*(minNLL[i-1]-minNLL[i]);
      int chi2dof = NDOF[i]-NDOF[i-1];
      chi2prob[i-1] = chi2<0 ? 1.0 : TMath::Prob(chi2,chi2dof);
      if(chi2prob[i-1]>0.05 && (modelNum!=0 || i!=1) ) break;
    }
    if(data[c]->sumEntries()-1<=NDOF[i]) break;
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

  MggBkgTmp[bestN]->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE));//strategy 1 or 2?
  w->import(*MggBkgTmp[bestN]);
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

  RooRealVar* mJJ     = w->var("mjj");
  mJJ->setUnit("GeV");
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);

  data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));

  RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("p1mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope1_cat%d",c)));
  RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("p2mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope2_cat%d",c)));
  RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("p3mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope3_cat%d",c)));
  RooFormulaVar *p4mod = new RooFormulaVar(TString::Format("p4mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope4_cat%d",c)));
  RooFormulaVar *p5mod = new RooFormulaVar(TString::Format("p5mod_cat%d",c),"","@0*@0",*w->var(TString::Format("mjj_bkg_8TeV_slope5_cat%d",c)));
  RooFormulaVar *p1arg = new RooFormulaVar(TString::Format("p1arg_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_arg1_cat%d",c)));
  RooFormulaVar *p2arg = new RooFormulaVar(TString::Format("p2arg_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_arg2_cat%d",c)));
  RooFormulaVar *p3arg = new RooFormulaVar(TString::Format("p3arg_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_arg3_cat%d",c)));
  RooFormulaVar *p4arg = new RooFormulaVar(TString::Format("p4arg_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_arg4_cat%d",c)));
  RooFormulaVar *p5arg = new RooFormulaVar(TString::Format("p5arg_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_arg5_cat%d",c)));
  RooFormulaVar *p1wid = new RooFormulaVar(TString::Format("p1wid_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_wid1_cat%d",c)));
  RooFormulaVar *p2wid = new RooFormulaVar(TString::Format("p2wid_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_wid2_cat%d",c)));
  RooFormulaVar *p3wid = new RooFormulaVar(TString::Format("p3wid_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_wid3_cat%d",c)));
  RooFormulaVar *p4wid = new RooFormulaVar(TString::Format("p4wid_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_wid4_cat%d",c)));
  RooFormulaVar *p5wid = new RooFormulaVar(TString::Format("p5wid_cat%d",c),"","@0",*w->var(TString::Format("mjj_bkg_8TeV_wid5_cat%d",c)));

  RooExponential* expo1 = new RooExponential("exp1","",*mJJ,*p1arg);
  RooExponential* expo2 = new RooExponential("exp2","",*mJJ,*p2arg);
  RooExponential* expo3 = new RooExponential("exp3","",*mJJ,*p3arg);
  RooExponential* expo4 = new RooExponential("exp4","",*mJJ,*p4arg);
  RooExponential* expo5 = new RooExponential("exp5","",*mJJ,*p5arg);
  RooLandau* lan1 = new RooLandau("lan1","",*mJJ,*p1arg,*p1wid);
  RooLandau* lan2 = new RooLandau("lan2","",*mJJ,*p2arg,*p2wid);
  RooLandau* lan3 = new RooLandau("lan3","",*mJJ,*p3arg,*p3wid);
  RooLandau* lan4 = new RooLandau("lan4","",*mJJ,*p4arg,*p4wid);
  RooLandau* lan5 = new RooLandau("lan5","",*mJJ,*p5arg,*p5wid);

  const int totalNDOF=5;
  int NDOF[totalNDOF]={0};
  float minNLL[totalNDOF]={0.0}, chi2prob[totalNDOF]={-1.0};
  int bestN = 0;

  RooAbsPdf* MjjBkgTmp[totalNDOF] = {0};
  char fitName[10];

  switch (modelNum){

  case 0: //Bernstein
    MjjBkgTmp[0] = new RooBernstein("BerN0Mjj", "", *mJJ, RooArgList(*p1mod));
    MjjBkgTmp[1] = new RooBernstein("BerN1Mjj", "", *mJJ, RooArgList(*p1mod,*p2mod));
    MjjBkgTmp[2] = new RooBernstein("BerN2Mjj", "", *mJJ, RooArgList(*p1mod,*p2mod,*p3mod));
    MjjBkgTmp[3] = new RooBernstein("BerN3Mjj", "", *mJJ, RooArgList(*p1mod,*p2mod,*p3mod,*p4mod));
    MjjBkgTmp[4] = new RooBernstein("BerN4Mjj", "", *mJJ, RooArgList(*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Ber");
    break;

  case 1: //Exponential
    w->factory(TString::Format("mjj_bkg_8TeV_norm_cat%d[800.0,0.0,100000]",c));
    MjjBkgTmp[0] = new RooExtendPdf("ExpN1Mjj","",*expo1,*w->var(TString::Format("mjj_bkg_8TeV_norm_cat%d",c)));
    MjjBkgTmp[1] = new RooAddPdf("ExpN2Mjj", "", RooArgList(*expo1,*expo2), RooArgList(*p1mod,*p2mod));
    MjjBkgTmp[2] = new RooAddPdf("ExpN3Mjj", "", RooArgList(*expo1,*expo2,*expo3), RooArgList(*p1mod,*p2mod,*p3mod));
    MjjBkgTmp[3] = new RooAddPdf("ExpN4Mjj", "", RooArgList(*expo1,*expo2,*expo3,*expo4), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod));
    MjjBkgTmp[4] = new RooAddPdf("ExpN5Mjj", "", RooArgList(*expo1,*expo2,*expo3,*expo4,*expo5), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Exp");
    break;

  case 2: //Landau
    w->factory(TString::Format("mjj_bkg_8TeV_norm_cat%d[800.0,0.0,100000]",c));
    MjjBkgTmp[0] = new RooExtendPdf("LanN1Mjj","",*lan1,*w->var(TString::Format("mjj_bkg_8TeV_norm_cat%d",c)));
    MjjBkgTmp[1] = new RooAddPdf("LanN2Mjj", "", RooArgList(*lan1,*lan2), RooArgList(*p1mod,*p2mod));
    MjjBkgTmp[2] = new RooAddPdf("LanN3Mjj", "", RooArgList(*lan1,*lan2,*lan3), RooArgList(*p1mod,*p2mod,*p3mod));
    MjjBkgTmp[3] = new RooAddPdf("LanN4Mjj", "", RooArgList(*lan1,*lan2,*lan3,*lan4), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod));
    MjjBkgTmp[4] = new RooAddPdf("LanN5Mjj", "", RooArgList(*lan1,*lan2,*lan3,*lan4,*lan5), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Lan");
    break;

  case 3: //Laurent
    MjjBkgTmp[0] = new RooGenericPdf("LauN1Mjj","@1*pow(@0,-4)",RooArgList(*mJJ,*p1mod));
    MjjBkgTmp[1] = new RooGenericPdf("LauN2Mjj","@1*pow(@0,-4)+@2*pow(@0,-3)",RooArgList(*mJJ,*p1mod,*p2mod));
    MjjBkgTmp[2] = new RooGenericPdf("LauN3Mjj","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)",RooArgList(*mJJ,*p1mod,*p2mod,*p3mod));
    MjjBkgTmp[3] = new RooGenericPdf("LauN4Mjj","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)+@4*pow(@0,-2)",RooArgList(*mJJ,*p1mod,*p2mod,*p3mod,*p4mod));
    MjjBkgTmp[4] = new RooGenericPdf("LauN5Mjj","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)+@4*pow(@0,-2)+@5*pow(@0,-6)",RooArgList(*mJJ,*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Lau");
    break;

  case 4: //Power
    MjjBkgTmp[0] = new RooGenericPdf("PowN1Mjj","@1*pow(@0,@2)",RooArgList(*mJJ,*p1mod,*p1arg));
    MjjBkgTmp[1] = new RooGenericPdf("PowN2Mjj","@1*pow(@0,@2)+@3*pow(@0,@4)",RooArgList(*mJJ,*p1mod,*p1arg,*p2mod,*p2arg));
    MjjBkgTmp[2] = new RooGenericPdf("PowN3Mjj","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)",RooArgList(*mJJ,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg));
    MjjBkgTmp[3] = new RooGenericPdf("PowN4Mjj","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)",RooArgList(*mJJ,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg,*p4mod,*p4arg));
    MjjBkgTmp[4] = new RooGenericPdf("PowN5Mjj","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)",RooArgList(*mJJ,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg,*p4mod,*p4arg));
    //MjjBkgTmp[4] = new RooGenericPdf(TString::Format("MjjBkg_cat%d",c),"@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)+@9*pow(@0,@10)",RooArgList(*mJJ,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg,*p4mod,*p4arg,*p5mod,*p5arg));
    sprintf(fitName,"Pow");
    break;
  }


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
    plotMjjBkg[c]->SetMinimum(0.0);
    plotMjjBkg[c]->SetMaximum(1.40*plotMjjBkg[c]->GetMaximum());
    plotMjjBkg[c]->GetXaxis()->SetTitle("m_{jj} (GeV)");
    plotMjjBkg[c]->Draw();  

    TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
    legmc->AddEntry(plotMjjBkg[c]->getObject(2),"Weighted CS","LPE");
    legmc->AddEntry(plotMjjBkg[c]->getObject(1),"Bkg Model","L");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    

    TLatex *lat  = new TLatex(minMassFit+3.0,0.85*plotMjjBkg[c]->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
    lat->Draw();
    TLatex *lat2 = new TLatex(minMassFit+3.0,0.7*plotMjjBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
 
    ctmp->SaveAs(TString::Format("plots/dataBkg2DMjj_%sN%i_cat%d.png",fitName,i+(modelNum!=0),c));

    if(i>0){
      float chi2 = 2*(minNLL[i-1]-minNLL[i]);
      int chi2dof = NDOF[i]-NDOF[i-1];
      chi2prob[i-1] = chi2<0 ? 1.0 : TMath::Prob(chi2,chi2dof);
      if(chi2prob[i-1]>0.05 && (modelNum!=0 || i!=1) ) break;
    }
    if(data[c]->sumEntries()-1<=NDOF[i]) break;
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
  return MjjBkgTmp[bestN];

}


void BkgModelBias(RooWorkspace* w,int c,RooAbsPdf* MggBkgTruth, RooAbsPdf* MjjBkgTruth, FILE *fout){

  std::vector<TString> catdesc;

  catdesc.push_back("#scale[0.8]{cat0}");
  catdesc.push_back("#scale[0.8]{cat1}");
  catdesc.push_back("#scale[0.8]{cat2}");
  catdesc.push_back("#scale[0.8]{cat3}");

  Float_t minMggMassFit(100),maxMggMassFit(180); 
  Float_t minMjjMassFit( 60),maxMjjMassFit(180); 
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

  RooRealVar *p0 = new RooRealVar("p0","",0.1,0,10000);
  RooRealVar *p1 = new RooRealVar("p1","",0.1,0,10000);
  RooRealVar *p2 = new RooRealVar("p2","",0.1,0,10000);
  RooRealVar *p3 = new RooRealVar("p3","",0.1,0,10000);
  RooRealVar *p4 = new RooRealVar("p4","",0.1,0,10000);
  RooRealVar *p5 = new RooRealVar("p5","",0.1,0,10000);
  RooRealVar *p6 = new RooRealVar("p6","",0.1,0,10000);
  RooRealVar *p7 = new RooRealVar("p7","",0.1,0,10000);
  RooRealVar *p8 = new RooRealVar("p8","",-0.1,-100000.0,0.0);
  RooRealVar *p9 = new RooRealVar("p9","",-0.1,-100000.0,0.0);

  const int totalNDOF=5;
  RooAbsPdf* MggBkgTmp[totalNDOF] = {0};
  RooAbsPdf* MjjBkgTmp[totalNDOF] = {0};
  MggBkgTmp[0] = new RooExponential(TString::Format("MggExp%d",1), "", *mGG,*p8);
  MggBkgTmp[1] = new RooGenericPdf("MggPow1","pow(@0,@1)",RooArgList(*mGG,*p8));
  MggBkgTmp[2] = new RooBernstein(TString::Format("MggPol%d",1), "", *mGG,RooArgList(*p0,*p1));
  MggBkgTmp[3] = new RooBernstein(TString::Format("MggPol%d",2), "", *mGG,RooArgList(*p0,*p1,*p2));
  MggBkgTmp[4] = new RooBernstein(TString::Format("MggPol%d",3), "", *mGG,RooArgList(*p0,*p1,*p2,*p3));
  MjjBkgTmp[0] = new RooExponential(TString::Format("MjjExp%d",1), "", *mGG,*p9);
  MjjBkgTmp[1] = new RooGenericPdf("MjjPow1","pow(@0,@1)",RooArgList(*mGG,*p9));
  MjjBkgTmp[2] = new RooBernstein(TString::Format("MjjPol%d",1), "", *mGG,RooArgList(*p4,*p5));
  MjjBkgTmp[3] = new RooBernstein(TString::Format("MjjPol%d",2), "", *mGG,RooArgList(*p4,*p5,*p6));
  MjjBkgTmp[4] = new RooBernstein(TString::Format("MjjPol%d",3), "", *mGG,RooArgList(*p4,*p5,*p6,*p7));


  if(MggBkgTruth->GetName()[0]=='B'){
    fprintf(fout,"Mgg x Mjj spectrum, bias results for cat%d\n",c);
    fprintf(fout,"Model\tExp\tPow\tBer1\tBer2\tBer3\n");
  }

  RooProdPdf *BkgTruth = new RooProdPdf("BkgTruth","",RooArgList(*MggBkgTruth,*MjjBkgTruth));

  const int Npse = 100;
  float results[totalNDOF];
  for(int k=0; k<totalNDOF; ++k){

    RooProdPdf *BkgFitTmp = new RooProdPdf("BkgFitTmp","",RooArgList(*MggBkgTmp[k],*MjjBkgTmp[k]));
    RooRealVar *nbkg = new RooRealVar("nbkg","",1,0,100000);
    RooExtendPdf *BkgFit = new RooExtendPdf(TString::Format("BkgFit_cat%d",c),"",*BkgFitTmp,*nbkg);

    RooMCStudy * mcs = new RooMCStudy(*BkgTruth, RooArgSet(*mGG,*mJJ), FitModel(*BkgFit),Silence(), Extended(kFALSE), Binned(kFALSE),//13.81 corresponds to a mu that gives p(N=0)=0.000001
				      FitOptions(Extended(kTRUE),PrintEvalErrors(0), Save()));
    RooChi2MCSModule chi2mod;
    mcs->addModule(chi2mod);

    if(c==0) mcs->generateAndFit(Npse*2,data->sumEntries(),kTRUE);
    else mcs->generateAndFit(Npse,data->sumEntries(),kTRUE);

    std::vector<double> pulls;
    //float genFraction = MggBkgTruth->createIntegral(*mGG,*mGG,"sigRegion")->getVal();
    for(int i=0; i<Npse; ++i){

      RooRealVar *q0 = new RooRealVar("q0","",mcs->fitParams(i)->getRealValue("p0"));
      RooRealVar *q1 = new RooRealVar("q1","",mcs->fitParams(i)->getRealValue("p1"));
      RooRealVar *q2 = new RooRealVar("q2","",mcs->fitParams(i)->getRealValue("p2"));
      RooRealVar *q3 = new RooRealVar("q3","",mcs->fitParams(i)->getRealValue("p3"));
      RooRealVar *q4 = new RooRealVar("q4","",mcs->fitParams(i)->getRealValue("p4"));
      RooRealVar *q5 = new RooRealVar("q5","",mcs->fitParams(i)->getRealValue("p5"));
      RooRealVar *q6 = new RooRealVar("q6","",mcs->fitParams(i)->getRealValue("p6"));
      RooRealVar *q7 = new RooRealVar("q7","",mcs->fitParams(i)->getRealValue("p7"));
      RooRealVar *q8 = new RooRealVar("q8","",mcs->fitParams(i)->getRealValue("p8"));
      RooRealVar *q9 = new RooRealVar("q9","",mcs->fitParams(i)->getRealValue("p9"));
      RooAbsPdf *fitFuncMgg,*fitFuncMjj;
      switch(k){
      case 0: fitFuncMgg = new RooExponential("fitFuncMgg","",*mGG,*q8); break;
      case 1: fitFuncMgg = new RooGenericPdf("fitFuncMgg","pow(@0,@1)",RooArgList(*mGG,*q8)); break;
      case 2: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1)); break;
      case 3: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1,*q2)); break;
      case 4: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1,*q2,*q3)); break;
      }
      switch(k){
      case 0: fitFuncMjj = new RooExponential("fitFuncMjj","",*mJJ,*q9); break;
      case 1: fitFuncMjj = new RooGenericPdf("fitFuncMjj","pow(@0,@1)",RooArgList(*mJJ,*q9)); break;
      case 2: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q4,*q5)); break;
      case 3: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q4,*q5,*q6)); break;
      case 4: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q4,*q5,*q6,*q7)); break;
      }

      RooProdPdf * fitFunc = new RooProdPdf("fitFunc","",RooArgSet(*fitFuncMgg,*fitFuncMjj));
      float fitFraction = fitFunc->createIntegral(RooArgSet(*mGG,*mJJ),RooArgSet(*mGG,*mJJ),"sigRegion")->getVal();
      //float genN = mcs->fitParams(i)->getRealValue("ngen");
      const RooAbsData* genDataset = mcs->genData(i);
      char fitRangeString[80]; sprintf(fitRangeString,"mgg>%f && mgg<%f && mjj>%f && mjj<%f",sigMeanMgg-sigFWHM_Mgg,sigMeanMgg+sigFWHM_Mgg,sigMeanMjj-sigFWHM_Mjj,sigMeanMjj+sigFWHM_Mjj);
      float genN = genDataset->sumEntries(fitRangeString);

      float fitN2 = mcs->fitParams(i)->getRealValue("nbkg");
      float a=0.318/2.0, nl=fitN2*fitFraction-0.5*TMath::ChisquareQuantile(a,2*fitN2*fitFraction), nh=0.5*TMath::ChisquareQuantile(1-a,2*(fitN2*fitFraction+1))-fitN2*fitFraction;
      float fitNerr2 = 0.5*(nl+nh);
      //float fitNerr2 = mcs->fitParams(i)->getRealValue("nbkgerr")*sqrt(sigFWHM*2/80.0);
      RooAbsReal *intRange = fitFunc->createIntegral(RooArgSet(*mGG,*mJJ),RooArgSet(*mGG,*mJJ),"sigRegion");
      nbkg = new RooRealVar("nbkg","",mcs->fitParams(i)->getRealValue("nbkg"),0,100000);
      RooFormulaVar *normIntRange= new RooFormulaVar(Form("normIntRange_m125_cat%d",c),Form("normIntRange_m125_cat%d",c),"@0*@1",RooArgList(*intRange,*nbkg));

      float fitN = normIntRange->getVal();
      float fitNerr = normIntRange->getPropagatedError(*mcs->fitResult(i));

      //cout<<"Pull check (NgenTot, Ngen, Nfit, fitErr, NfitOld, fitErrOld): "<<genDataset->sumEntries()<<' '<<genN<<' '<<fitN<<' '<<fitNerr<<' '<<fitN2*fitFraction<<' '<<fitNerr2<<' '<<endl;
      if(fitNerr>0)// && mcs->fitParams(i)->getRealValue("chi2")/mcs->fitParams(i)->getRealValue("ndof")<10)
	pulls.push_back((genN-fitN)/(fitNerr2));
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
    //TH1F *h1 = new TH1F("h1","",50,pulls[0],pulls[pulls.size()-1]);
    TH1F *h1 = new TH1F("h1","",25,-3,3);
    h1->GetXaxis()->SetTitle("pull(N_{bgd}^{SR})");
    for(int i=0; i<pulls.size(); i++) h1->Fill(pulls[i]);
    TCanvas *c0 = new TCanvas("c0","c0",700,500);
    h1->SetMaximum(h1->GetMaximum()*1.40);
    h1->Draw();
    TFitResultPtr pullFit;
    if(h1->Integral()>4){
      pullFit = h1->Fit("gaus","s");
      TLatex *lat1  = new TLatex(-2.6,0.85*h1->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
      lat1->Draw();
      TLatex *lat12 = new TLatex(0.8,0.60*h1->GetMaximum(),TString::Format("#splitline{#scale[1.0]{#mu = %.2f #pm %.2f}}{#scale[1.0]{#sigma = %.2f #pm %.2f}}",pullFit->GetParams()[1],pullFit->GetErrors()[1],pullFit->GetParams()[2],pullFit->GetErrors()[2]));
      lat12->SetTextColor(kRed);
      lat12->Draw();
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(0.8,0.91*h1->GetMaximum(),TString::Format("Num PSE: %d",pulls.size()) );
      lat13->DrawLatex(0.8,0.85*h1->GetMaximum(),TString::Format("Gen function: %s",MggBkgTruth->GetName()) );
      lat13->DrawLatex(0.8,0.79*h1->GetMaximum(),TString::Format("Fit function: %s",MggBkgTmp[k]->GetName()) );
      lat13->DrawLatex(0.8,0.73*h1->GetMaximum(),TString::Format("Median = %.2f ",results[k]) );
    }
    else{
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(0.0,0.5,"Not enough entries." );
      lat13->DrawLatex(0.0,0.4,TString::Format("Median = %.2f ",results[k])  );
    }

    c0->SaveAs(TString::Format("plots/pulls2D_gen%s_fitN%d_cat%d.png",MggBkgTruth->GetName(),k+1,c));  

    TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
    c1->Divide(3,3);
    for(int i=1;i<=9; ++i){
      c1->cd(i);
      RooPlot *frame = mGG->frame();
      if(i==1){
	data->plotOn(frame);
	MggBkgTruth->plotOn(frame,LineColor(kBlue),Range("fitrange"),NormRange("fitrange"));
	frame->SetTitle("");
	frame->SetMinimum(0.0);
	frame->SetMaximum(1.40*frame->GetMaximum());
	frame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
	frame->Draw();

	TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
	legmc->AddEntry(frame->getObject(0),"Weighted CS","LPE");
	legmc->AddEntry(frame->getObject(1),"Bkg Model","L");
	legmc->SetBorderSize(0);
	legmc->SetFillStyle(0);
	legmc->Draw();    
	TLatex *lat  = new TLatex(minMggMassFit+3.0,0.85*frame->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
	lat->Draw();
	TLatex *lat2 = new TLatex(minMggMassFit+3.0,0.7*frame->GetMaximum(),catdesc.at(c));
	lat2->Draw();

	continue;
      }
      mcs->genData(i)->plotOn(frame);
      RooRealVar *q0 = new RooRealVar("q0","",mcs->fitParams(i)->getRealValue("p0"));
      RooRealVar *q1 = new RooRealVar("q1","",mcs->fitParams(i)->getRealValue("p1"));
      RooRealVar *q2 = new RooRealVar("q2","",mcs->fitParams(i)->getRealValue("p2"));
      RooRealVar *q3 = new RooRealVar("q3","",mcs->fitParams(i)->getRealValue("p3"));
      RooRealVar *q4 = new RooRealVar("q4","",mcs->fitParams(i)->getRealValue("p4"));
      RooRealVar *q5 = new RooRealVar("q5","",mcs->fitParams(i)->getRealValue("p5"));
      RooRealVar *q6 = new RooRealVar("q6","",mcs->fitParams(i)->getRealValue("p6"));
      RooRealVar *q7 = new RooRealVar("q7","",mcs->fitParams(i)->getRealValue("p7"));
      RooRealVar *q8 = new RooRealVar("q8","",mcs->fitParams(i)->getRealValue("p8"));
      RooRealVar *q9 = new RooRealVar("q9","",mcs->fitParams(i)->getRealValue("p9"));
      RooAbsPdf *fitFuncMgg, *fitFuncMjj;
      switch(k){
      case 0: fitFuncMgg = new RooExponential("fitFuncMgg","",*mGG,*q8); break;
      case 1: fitFuncMgg = new RooGenericPdf("fitFuncMgg","pow(@0,@1)",RooArgList(*mGG,*q8)); break;
      case 2: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1)); break;
      case 3: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1,*q2)); break;
      case 4: fitFuncMgg = new RooBernstein("fitFuncMgg","",*mGG,RooArgList(*q0,*q1,*q2,*q3)); break;
      }
      switch(k){
      case 0: fitFuncMjj = new RooExponential("fitFuncMjj","",*mJJ,*q9); break;
      case 1: fitFuncMjj = new RooGenericPdf("fitFuncMjj","pow(@0,@1)",RooArgList(*mJJ,*q9)); break;
      case 2: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q4,*q5)); break;
      case 3: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q4,*q5,*q6)); break;
      case 4: fitFuncMjj = new RooBernstein("fitFuncMjj","",*mJJ,RooArgList(*q4,*q5,*q6,*q7)); break;
      }
      RooProdPdf *fitFunc= new RooProdPdf("fitFunc","",RooArgSet(*fitFuncMgg,*fitFuncMjj));
      fitFunc->plotOn(frame,LineColor(kGreen),Range("fitrange"),NormRange("fitrange"));
      frame->SetTitle("");                                                                             
      frame->SetMinimum(0.0);                                                                          
      frame->SetMaximum(1.40*frame->GetMaximum());                                                     
      frame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");                                           
      frame->Draw();

      TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
      legmc->AddEntry(frame->getObject(0),TString::Format("Gen PSE #%i",i-1),"LPE");
      legmc->AddEntry(frame->getObject(1),"Bkg Model","L");
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();    
      TLatex *lat  = new TLatex(minMggMassFit+3.0,0.85*frame->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
      lat->Draw();
      TLatex *lat2 = new TLatex(minMggMassFit+3.0,0.7*frame->GetMaximum(),catdesc.at(c));
      lat2->Draw();
      TLatex *lat3 = new TLatex();
      lat3->SetNDC();
      lat3->DrawLatex(0.62,0.68,TString::Format("#chi^{2}/N = %.2f",mcs->fitParams(i)->getRealValue("chi2")/mcs->fitParams(i)->getRealValue("ndof")));
      lat3->DrawLatex(0.62,0.62,TString::Format("p(#chi^{2},N) = %.2f",mcs->fitParams(i)->getRealValue("prob")));
    }
    c1->SaveAs(TString::Format("plots/toymc2D_gen%s_fitN%d_cat%d.png",MggBkgTruth->GetName(),k+1,c));

    delete mcs;
  }

  fprintf(fout,"%s,%s\t%8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",MggBkgTruth->GetName(),MjjBkgTruth->GetName(),results[0],results[1],results[2],results[3],results[4]);

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
