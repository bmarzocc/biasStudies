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

Int_t NCAT = 3;
TString inDir   = "./Trees/";

void AddBkgData(RooWorkspace*, int);
RooAbsPdf* BkgModelFit(RooWorkspace*, int, int);
void BkgModelBias(RooWorkspace*,int,RooAbsPdf*,FILE*);
void SetParamNames(RooWorkspace*);
void SetConstantParams(const RooArgSet* params);
Double_t effSigma(TH1 *hist);
 
RooArgSet* defineVariables()
{
  // define variables of the input ntuple
  RooRealVar* mJJ  = new RooRealVar("mJJ","M(jj)",60,180,"GeV");
  RooRealVar* mGG  = new RooRealVar("mGG","M(#gamma#gamma)",0,3000,"GeV");
  RooRealVar* mRad = new RooRealVar("mRad","M(#gamma#gamma jj)",0,1500,"GeV");
  RooRealVar* pt    = new RooRealVar("pt","p_{T}(jj)",0,3000,"GeV");
  RooRealVar* pt1   = new RooRealVar("pt1","p_{T}^{1}(jj)",0.0,3000,"GeV");
  RooRealVar* pt2   = new RooRealVar("pt2","p_{T}^{2}(jj)",0.0,3000,"GeV");
  RooRealVar* wei   = new RooRealVar("wei","HqT x PUwei",0,100,"");
  RooRealVar* weiS  = new RooRealVar("weiS","HqT x PUwei",0,1000,"");
  RooCategory* bJetTagCategory = new RooCategory("bJetTagCategory","event category 4") ;
  bJetTagCategory->defineType("cat4_0",0);
  bJetTagCategory->defineType("cat4_1",1);
  bJetTagCategory->defineType("cat4_2",2);
  bJetTagCategory->defineType("cat4_3",3);

  RooCategory* catMva = new RooCategory("catMva","MVA event category") ;
  catMva->defineType("MVAcat_0",0);
  catMva->defineType("MVAcat_1",1);
  catMva->defineType("MVAcat_2",2);
  catMva->defineType("MVAcat_3",3);
  catMva->defineType("MVAcat_4",4);
  catMva->defineType("MVAcat_5",5);
  catMva->defineType("MVAcat_6",6);
  catMva->defineType("MVAcat_7",7);
  catMva->defineType("MVAcat_8",8);



  RooCategory* catjet = new RooCategory("CATJET","event category par Njets") ;
  catjet->defineType("catjet_0",0);
  catjet->defineType("catjet_1",1);

  RooCategory* vbftag = new RooCategory("vbfTag","Passed VBF selection") ;
  vbftag->defineType("vbftag_0",0);
  vbftag->defineType("vbftag_1",1);

  RooCategory* vhtag = new RooCategory("metTag","Passed VH hadronic selection") ;
  vhtag->defineType("vhtag_0",0);
  vhtag->defineType("vhtag_1",1);

  RooCategory* leptag = new RooCategory("lepTag","Passed lepton selection") ;
  leptag->defineType("leptag_0",0);
  leptag->defineType("leptag_1",1);

  RooCategory* excltag = new RooCategory("NOEXCLTAG","Passed exclusive selection") ;
  excltag->defineType("excltag_0",0);
  excltag->defineType("excltag_1",1);

  RooRealVar* MITBDTG   = new RooRealVar("MITBDTG","MIT BDTG",0.05,1,"");
  RooRealVar* phomvaid1  = new RooRealVar("phomvaid1","phomvaid 1",-0.3,1,"");
  RooRealVar* phomvaid2  = new RooRealVar("phomvaid2","phomvaid 2",-0.3,1,"");

  RooArgSet* ntplVars = new RooArgSet(*mJJ,*mGG,*pt,*pt1, *pt2, *bJetTagCategory, *wei, *weiS, *catjet);
  ntplVars->add(*mRad);
  ntplVars->add(*vbftag);
  ntplVars->add(*vhtag);
  ntplVars->add(*leptag);
  ntplVars->add(*catMva);
  ntplVars->add(*MITBDTG);
  ntplVars->add(*phomvaid1);   
  ntplVars->add(*phomvaid2);   
   
 
  return ntplVars;
}


void runfits(int cat=0, int modelNum=-1, int inDirNum=0)
{

  //create truth models
  const int nTruthModels=5;
  RooAbsPdf *MjjBkgTruth[nTruthModels] = {0};

  switch(inDirNum){
  case 0: inDir=" /afs/cern.ch/work/h/hebda/HggHbb/ntuples/TreesCS/"; break;
  case 1: inDir=" /afs/cern.ch/work/h/hebda/HggHbb/ntuples/Trees_reg/"; break;
  case 2: inDir=" /afs/cern.ch/work/h/hebda/HggHbb/ntuples/Trees_reg_kinfit/"; break;
  case 3: inDir=" /afs/cern.ch/work/h/hebda/HggHbb/ntuples/Trees_kinfit/"; break;
  }

  TString card_name("hbbSM_models_Pol_8TeV.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  AddBkgData(w,cat);

  FILE *fout = fopen("resultsBiasMjj.txt","a");
  if(modelNum<=0) fprintf(fout,"%s\n",inDir.Data());

  for(int m=0; m<nTruthModels; ++m) {
    if(modelNum>=0 && m!=modelNum) continue;
    if(m==2) continue;//skip Landau, it sucks.
    MjjBkgTruth[m] = BkgModelFit(w,cat,m); //Ber, Exp, Lan, Lau, Pow
    BkgModelBias(w,cat,MjjBkgTruth[m],fout);
  }

  if(modelNum<0 || modelNum+1==nTruthModels) fprintf(fout,"\n\n");
  fclose(fout);
  return;
}


void AddBkgData(RooWorkspace* w, int cat) {

  Int_t ncat = NCAT;

// common preselection cut
  //TString mainCut("1");
  TString mainCut("mGG>122 && mGG<128");// && mRad>225 && mRad<295");
  //if(cat==0) mainCut = "122<mGG && mGG<128";
  //else if(cat==1) mainCut = "122<mGG && mGG<129";
  //else if(cat==2) mainCut = "122<mGG && mGG<126";

//****************************//
// Signal Data Set
//****************************//

  // Variables
  RooArgSet* ntplVars = defineVariables();
  RooRealVar weightVar("weightVar","",1,0,1000);

//****************************//
// CMS Data Set
//****************************//
// retrieve the data tree;
// no common preselection cut applied yet; 

  TFile dataFile(inDir+"tree_data_weiMjj.root");   
  TTree* dataTree     = (TTree*) dataFile.Get("Events");
  weightVar.setVal(1.);
  ntplVars->add(RooArgList(weightVar));

  RooDataSet Data("Data","dataset",dataTree,*ntplVars,"","wei");

// apply a common preselection cut;
// split into NCAT  categories;

  
  RooDataSet* dataToFit[9];
  for (int c = 0; c < ncat; ++c) {
// Real data
    dataToFit[c]   = (RooDataSet*) Data.reduce(*w->var("mJJ"),mainCut+TString::Format(" && bJetTagCategory==%d",c));
    w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
  }

// Create full data set without categorization
  RooDataSet* data    = (RooDataSet*) Data.reduce(*w->var("mJJ"),mainCut);
  w->import(*data, Rename("Data"));
  data->Print("v");

}


RooAbsPdf *BkgModelFit(RooWorkspace* w, int c, int modelNum) {

  std::vector<TString> catdesc;

  catdesc.push_back("#scale[0.8]{Untagged}");
  catdesc.push_back("#scale[0.8]{1-bjet tagged}");
  catdesc.push_back("#scale[0.8]{2-bjet tagged}");

  RooDataSet* data[9];
  RooFitResult* fitresult[9];;
  RooPlot* plotMjjBkg[9];

  Float_t minMassFit(60),maxMassFit(180); 

  RooRealVar* mJJ     = w->var("mJJ");  
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
    MjjBkgTmp[0] = new RooBernstein("BerN0", "", *mJJ, RooArgList(*p1mod));
    MjjBkgTmp[1] = new RooBernstein("BerN1", "", *mJJ, RooArgList(*p1mod,*p2mod));
    MjjBkgTmp[2] = new RooBernstein("BerN2", "", *mJJ, RooArgList(*p1mod,*p2mod,*p3mod));
    MjjBkgTmp[3] = new RooBernstein("BerN3", "", *mJJ, RooArgList(*p1mod,*p2mod,*p3mod,*p4mod));
    MjjBkgTmp[4] = new RooBernstein("BerN4", "", *mJJ, RooArgList(*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Ber");
    break;

  case 1: //Exponential
    w->factory(TString::Format("mjj_bkg_8TeV_norm_cat%d[800.0,0.0,100000]",c));
    MjjBkgTmp[0] = new RooExtendPdf("ExpN1","",*expo1,*w->var(TString::Format("mjj_bkg_8TeV_norm_cat%d",c)));
    MjjBkgTmp[1] = new RooAddPdf("ExpN2", "", RooArgList(*expo1,*expo2), RooArgList(*p1mod,*p2mod));
    MjjBkgTmp[2] = new RooAddPdf("ExpN3", "", RooArgList(*expo1,*expo2,*expo3), RooArgList(*p1mod,*p2mod,*p3mod));
    MjjBkgTmp[3] = new RooAddPdf("ExpN4", "", RooArgList(*expo1,*expo2,*expo3,*expo4), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod));
    MjjBkgTmp[4] = new RooAddPdf("ExpN5", "", RooArgList(*expo1,*expo2,*expo3,*expo4,*expo5), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Exp");
    break;

  case 2: //Landau
    w->factory(TString::Format("mjj_bkg_8TeV_norm_cat%d[800.0,0.0,100000]",c));
    MjjBkgTmp[0] = new RooExtendPdf("LanN1","",*lan1,*w->var(TString::Format("mjj_bkg_8TeV_norm_cat%d",c)));
    MjjBkgTmp[1] = new RooAddPdf("LanN2", "", RooArgList(*lan1,*lan2), RooArgList(*p1mod,*p2mod));
    MjjBkgTmp[2] = new RooAddPdf("LanN3", "", RooArgList(*lan1,*lan2,*lan3), RooArgList(*p1mod,*p2mod,*p3mod));
    MjjBkgTmp[3] = new RooAddPdf("LanN4", "", RooArgList(*lan1,*lan2,*lan3,*lan4), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod));
    MjjBkgTmp[4] = new RooAddPdf("LanN5", "", RooArgList(*lan1,*lan2,*lan3,*lan4,*lan5), RooArgList(*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Lan");
    break;

  case 3: //Laurent
    MjjBkgTmp[0] = new RooGenericPdf("LauN1","@1*pow(@0,-4)",RooArgList(*mJJ,*p1mod));
    MjjBkgTmp[1] = new RooGenericPdf("LauN2","@1*pow(@0,-4)+@2*pow(@0,-3)",RooArgList(*mJJ,*p1mod,*p2mod));
    MjjBkgTmp[2] = new RooGenericPdf("LauN3","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)",RooArgList(*mJJ,*p1mod,*p2mod,*p3mod));
    MjjBkgTmp[3] = new RooGenericPdf("LauN4","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)+@4*pow(@0,-2)",RooArgList(*mJJ,*p1mod,*p2mod,*p3mod,*p4mod));
    MjjBkgTmp[4] = new RooGenericPdf("LauN5","@1*pow(@0,-4)+@2*pow(@0,-3)+@3*pow(@0,-5)+@4*pow(@0,-2)+@5*pow(@0,-6)",RooArgList(*mJJ,*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    sprintf(fitName,"Lau");
    break;

  case 4: //Power
    MjjBkgTmp[0] = new RooGenericPdf("PowN1","@1*pow(@0,@2)",RooArgList(*mJJ,*p1mod,*p1arg));
    MjjBkgTmp[1] = new RooGenericPdf("PowN2","@1*pow(@0,@2)+@3*pow(@0,@4)",RooArgList(*mJJ,*p1mod,*p1arg,*p2mod,*p2arg));
    MjjBkgTmp[2] = new RooGenericPdf("PowN3","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)",RooArgList(*mJJ,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg));
    MjjBkgTmp[3] = new RooGenericPdf("PowN4","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)",RooArgList(*mJJ,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg,*p4mod,*p4arg));
    MjjBkgTmp[4] = new RooGenericPdf("PowN5","@1*pow(@0,@2)+@3*pow(@0,@4)+@5*pow(@0,@6)+@7*pow(@0,@8)",RooArgList(*mJJ,*p1mod,*p1arg,*p2mod,*p2arg,*p3mod,*p3arg,*p4mod,*p4arg));
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
    plotMjjBkg[c] = mJJ->frame(nBinsMass);
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
 
    //ctmp->SaveAs(TString::Format("databkgoversig_cat%d.pdf",c));
    //ctmp->SaveAs(TString::Format("databkgoversig_cat%d.eps",c));
    ctmp->SaveAs(TString::Format("plots/dataBkgMjj_%sN%i_cat%d.png",fitName,i+(modelNum!=0),c));
    //ctmp->SaveAs(TString::Format("databkgoversig_cat%d.C",c));

    if(i>0){
      float chi2 = 2*(minNLL[i-1]-minNLL[i]);
      int chi2dof = NDOF[i]-NDOF[i-1];
      chi2prob[i-1] = chi2<0 ? 1.0 : TMath::Prob(chi2,chi2dof);
      if(chi2prob[i-1]>0.05 && (modelNum!=0 || i!=1) ) break;
    }
    if(data[c]->sumEntries()-1<=NDOF[i]) break;
    bestN=i;
  }

  FILE *results = fopen("resultsModelMjj.txt","a");
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


void BkgModelBias(RooWorkspace* w,int c,RooAbsPdf* MjjBkgTruth, FILE *fout){

  std::vector<TString> catdesc;

  catdesc.push_back("#scale[0.8]{Untagged}");
  catdesc.push_back("#scale[0.8]{1-bjet tagged}");
  catdesc.push_back("#scale[0.8]{2-bjet tagged}");

  Float_t minMassFit(60),maxMassFit(180); 
  float sigMean=0, sigFWHM=0;
  switch(c){
  case 0: sigMean=118.43; sigFWHM=30.25; break;
  case 1: sigMean=117.23; sigFWHM=25.52; break;
  case 2: sigMean=121.89; sigFWHM=17.75; break;
  }

  RooRealVar* mJJ     = w->var("mJJ");  
  mJJ->setUnit("GeV");
  mJJ->setRange("sigRegion",sigMean-sigFWHM,sigMean+sigFWHM);
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);

  RooDataSet* data = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));

  RooRealVar *p0 = new RooRealVar("p0","",0.1,0.0,10000);
  RooRealVar *p1 = new RooRealVar("p1","",0.1,0.0,10000);
  RooRealVar *p2 = new RooRealVar("p2","",0.1,0.0,10000);
  RooRealVar *p3 = new RooRealVar("p3","",0.1,0.0,10000);
  RooRealVar *p4 = new RooRealVar("p4","",0.1,0.0,10000);
  RooRealVar *p5 = new RooRealVar("p5","",0.1,0.0,10000);

  const int totalNDOF=5;
  RooAbsPdf* MjjBkgTmp[totalNDOF] = {0};
  MjjBkgTmp[0] = new RooBernstein(TString::Format("Pol%d",1), "", *mJJ,RooArgList(*p0,*p1));
  MjjBkgTmp[1] = new RooBernstein(TString::Format("Pol%d",2), "", *mJJ,RooArgList(*p0,*p1,*p2));
  MjjBkgTmp[2] = new RooBernstein(TString::Format("Pol%d",3), "", *mJJ,RooArgList(*p0,*p1,*p2,*p3));
  MjjBkgTmp[3] = new RooBernstein(TString::Format("Pol%d",4), "", *mJJ,RooArgList(*p0,*p1,*p2,*p3,*p4));
  MjjBkgTmp[4] = new RooBernstein(TString::Format("Pol%d",5), "", *mJJ,RooArgList(*p0,*p1,*p2,*p3,*p4,*p5));

  if(MjjBkgTruth->GetName()[0]=='B'){
    fprintf(fout,"Mjj spectrum, bias results for cat%d\n",c);
    fprintf(fout,"Model\t    %s      %s      %s      %s      %s\n",
	    MjjBkgTmp[0]->GetName(),MjjBkgTmp[1]->GetName(),MjjBkgTmp[2]->GetName(),MjjBkgTmp[3]->GetName(),MjjBkgTmp[4]->GetName());
  }

  const int Npse = 1000;
  float results[totalNDOF];
  for(int k=0; k<totalNDOF; ++k){
    RooRealVar *nbkg = new RooRealVar("nbkg","",1,0,100000);
    RooExtendPdf *MjjBkgFit = new RooExtendPdf(TString::Format("MjjBkgFit_cat%d",c),"",*MjjBkgTmp[k],*nbkg);

    RooMCStudy * mcs = new RooMCStudy(*MjjBkgTruth, *mJJ, FitModel(*MjjBkgFit),Silence(), Extended(kTRUE), Binned(kFALSE),
				      FitOptions(Range(minMassFit,maxMassFit),Extended(kTRUE),PrintEvalErrors(0)));
    RooChi2MCSModule chi2mod;
    mcs->addModule(chi2mod);

    if(c==2) mcs->generateAndFit(Npse*2,data->sumEntries(),kTRUE);
    else mcs->generateAndFit(Npse,data->sumEntries(),kTRUE);

    std::vector<double> pulls;
    //float genFraction = MjjBkgTruth->createIntegral(*mJJ,*mJJ,"sigRegion")->getVal();
    for(int i=0; i<Npse; ++i){

      RooRealVar *q0 = new RooRealVar("q0","",mcs->fitParams(i)->getRealValue("p0"));
      RooRealVar *q1 = new RooRealVar("q1","",mcs->fitParams(i)->getRealValue("p1"));
      RooRealVar *q2 = new RooRealVar("q2","",mcs->fitParams(i)->getRealValue("p2"));
      RooRealVar *q3 = new RooRealVar("q3","",mcs->fitParams(i)->getRealValue("p3"));
      RooRealVar *q4 = new RooRealVar("q4","",mcs->fitParams(i)->getRealValue("p4"));
      RooRealVar *q5 = new RooRealVar("q5","",mcs->fitParams(i)->getRealValue("p5"));
      RooBernstein *fitFunc;
      switch(k){
      case 0: fitFunc = new RooBernstein("fitFunc","",*mJJ,RooArgList(*q0,*q1)); break;
      case 1: fitFunc = new RooBernstein("fitFunc","",*mJJ,RooArgList(*q0,*q1,*q2)); break;
      case 2: fitFunc = new RooBernstein("fitFunc","",*mJJ,RooArgList(*q0,*q1,*q2,*q3)); break;
      case 3: fitFunc = new RooBernstein("fitFunc","",*mJJ,RooArgList(*q0,*q1,*q2,*q3,*q4)); break;
      case 4: fitFunc = new RooBernstein("fitFunc","",*mJJ,RooArgList(*q0,*q1,*q2,*q3,*q4,*q5)); break;
      }

      float fitFraction = fitFunc->createIntegral(*mJJ,*mJJ,"sigRegion")->getVal();
      //float genN = mcs->fitParams(i)->getRealValue("ngen");
      //float genN = mcs->fitParams(i)->getRealValue("ngen");
      const RooAbsData* genDataset = mcs->genData(i);
      char fitRangeString[80]; sprintf(fitRangeString,"mJJ>%f && mJJ<%f",sigMean-sigFWHM,sigMean+sigFWHM);
      float genN = genDataset->sumEntries(fitRangeString);

      float fitN = mcs->fitParams(i)->getRealValue("nbkg");
      float a=0.318/2.0, nl=fitN*fitFraction-0.5*TMath::ChisquareQuantile(a,2*fitN*fitFraction), nh=0.5*TMath::ChisquareQuantile(1-a,2*(fitN*fitFraction+1))-fitN*fitFraction;
      float fitNerr = 0.5*(nl+nh);
      if(fitNerr*fitFraction>0)// && mcs->fitParams(i)->getRealValue("chi2")/mcs->fitParams(i)->getRealValue("ndof")<10)
	pulls.push_back((genN-fitN*fitFraction)/(fitNerr));
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
      lat13->DrawLatex(0.8,0.85*h1->GetMaximum(),TString::Format("Gen function: %s",MjjBkgTruth->GetName()) );
      lat13->DrawLatex(0.8,0.79*h1->GetMaximum(),TString::Format("Fit function: %s",MjjBkgTmp[k]->GetName()) );
      lat13->DrawLatex(0.8,0.73*h1->GetMaximum(),TString::Format("Median = %.2f ",results[k]) );
    }
    else{
      TLatex *lat13 = new TLatex();
      lat13->DrawLatex(0.0,0.5,"Not enough entries." );
      lat13->DrawLatex(0.0,0.4,TString::Format("Median = %.2f ",results[k])  );
    }

    c0->SaveAs(TString::Format("plots/pullsMjj_gen%s_fitN%d_cat%d.png",MjjBkgTruth->GetName(),k+1,c));  

    TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
    c1->Divide(3,3);
    for(int i=1;i<=9; ++i){
      c1->cd(i);
      RooPlot *frame = mJJ->frame();
      if(i==1){
	data->plotOn(frame);
	MjjBkgTruth->plotOn(frame,LineColor(kBlue),Range("fitrange"),NormRange("fitrange"));
	frame->SetTitle("");
	frame->SetMinimum(0.0);
	frame->SetMaximum(1.40*frame->GetMaximum());
	frame->GetXaxis()->SetTitle("m_{jj} (GeV)");
	frame->Draw();

	TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
	legmc->AddEntry(frame->getObject(0),"Weighted CS","LPE");
	legmc->AddEntry(frame->getObject(1),"Bkg Model","L");
	legmc->SetBorderSize(0);
	legmc->SetFillStyle(0);
	legmc->Draw();    
	TLatex *lat  = new TLatex(minMassFit+3.0,0.85*frame->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
	lat->Draw();
	TLatex *lat2 = new TLatex(minMassFit+3.0,0.7*frame->GetMaximum(),catdesc.at(c));
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
      RooBernstein *fitFunc;
      switch(k){
      case 0: fitFunc = new RooBernstein("fitFunc","",*mJJ,RooArgList(*q0,*q1)); break;
      case 1: fitFunc = new RooBernstein("fitFunc","",*mJJ,RooArgList(*q0,*q1,*q2)); break;
      case 2: fitFunc = new RooBernstein("fitFunc","",*mJJ,RooArgList(*q0,*q1,*q2,*q3)); break;
      case 3: fitFunc = new RooBernstein("fitFunc","",*mJJ,RooArgList(*q0,*q1,*q2,*q3,*q4)); break;
      case 4: fitFunc = new RooBernstein("fitFunc","",*mJJ,RooArgList(*q0,*q1,*q2,*q3,*q4,*q5)); break;
      }
      fitFunc->plotOn(frame,LineColor(kGreen),Range("fitrange"),NormRange("fitrange"));
      frame->SetTitle("");                                                                             
      frame->SetMinimum(0.0);                                                                          
      frame->SetMaximum(1.40*frame->GetMaximum());                                                     
      frame->GetXaxis()->SetTitle("m_{jj} (GeV)");
      frame->Draw();

      TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
      legmc->AddEntry(frame->getObject(0),TString::Format("Gen PSE #%i",i-1),"LPE");
      legmc->AddEntry(frame->getObject(1),"Bkg Model","L");
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();    
      TLatex *lat  = new TLatex(minMassFit+3.0,0.85*frame->GetMaximum(),"#splitline{#scale[1.0]{CMS Preliminary}}{#scale[0.8]{#sqrt{s} = 8 TeV}}");
      lat->Draw();
      TLatex *lat2 = new TLatex(minMassFit+3.0,0.7*frame->GetMaximum(),catdesc.at(c));
      lat2->Draw();
      TLatex *lat3 = new TLatex();
      lat3->SetNDC();
      lat3->DrawLatex(0.62,0.68,TString::Format("#chi^{2}/N = %.2f",mcs->fitParams(i)->getRealValue("chi2")/mcs->fitParams(i)->getRealValue("ndof")));
      lat3->DrawLatex(0.62,0.62,TString::Format("p(#chi^{2},N) = %.2f",mcs->fitParams(i)->getRealValue("prob")));
    }
    c1->SaveAs(TString::Format("plots/toymcMjj_gen%s_fitN%d_cat%d.png",MjjBkgTruth->GetName(),k+1,c));

    delete mcs;
  }

  fprintf(fout,"%s\t%8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",MjjBkgTruth->GetName(),results[0],results[1],results[2],results[3],results[4]);

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
