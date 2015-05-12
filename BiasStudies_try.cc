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
#include <fstream>
#include <iostream>
#include <memory>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
 #include <fstream>
#include <iterator>
#include <algorithm>

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
#include "RooArgusBG.h"
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

#include "RooGlobalFunc.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "RooFitResult.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TRandom.h"

using namespace RooFit;
using namespace RooStats ;

int main(int argc, const char* argv[])
{
  int nEvents = atoi(argv[1]);
  int Npse = atoi(argv[2]);

  TFile* f = new TFile((std::string("output_")+std::string(argv[1])+".root").c_str(),"RECREATE");
  
  TH1F* h_pulls = new TH1F("h_pulls","h_pulls",6000,-3,3);

  TString card_name("models_2D.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  
  // --- Observable --- 
  RooRealVar mes("mgg","mgg (GeV)",100.,180.) ; 
  // --- Build Gaussian signal PDF --- 
  RooRealVar sigmean("sigmean","mgg mass",125.,123.,127.) ; 
  RooRealVar sigwidth("sigwidth","mgg width",1.7,0.7,2.7) ; 
  RooGaussian gauss("gauss","gaussian PDF",mes,sigmean,sigwidth) ; 
  // --- Build Expnential background PDF --- 
  RooRealVar a("a","-a",-0.1,-100000.0,0.) ; 
  RooExponential expo("exp","",mes,a);

  // -- Build Berstein background PDF ---
  RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("mggp1mod_cat%d",0),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",0)));
  RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("mggp2mod_cat%d",0),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",0)));
  RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("mggp3mod_cat%d",0),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",0)));
  RooFormulaVar *p4mod = new RooFormulaVar(TString::Format("mggp4mod_cat%d",0),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope4_cat%d",0)));
  RooFormulaVar *p5mod = new RooFormulaVar(TString::Format("mggp5mod_cat%d",0),"","@0*@0",*w->var(TString::Format("mgg_bkg_8TeV_slope5_cat%d",0)));
  RooBernstein bern1("BerN1Mgg", "",mes, RooArgList(*p1mod,*p2mod));
  RooBernstein bern4("BerN4Mgg", "",mes, RooArgList(*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));

  // --- Construct signal+background PDF --- 
  RooRealVar nsig("nsig","#signal events",0.,-1000.,1000.) ; 
  RooRealVar nbkg("nbkg","#background events",nEvents,0.5*nEvents,1.5*nEvents) ; 
  //RooAddPdf sum("sum","g+a",RooArgList(gauss,expo),RooArgList(nsig,nbkg)) ; 
  //RooAddPdf sum("sum","g+a",RooArgList(gauss,bern1),RooArgList(nsig,nbkg)) ; 
  RooAddPdf sum("sum","g+a",RooArgList(gauss,bern4),RooArgList(nsig,nbkg)) ; 
  
  // --- Constraints ---
  float tmp_sigma_bkg = sqrt(nEvents);
  float tmp_sigma_sig = 0.1*6./80.*nEvents;//sqrt(sigFrac*data->sumEntries());
  RooRealVar mean_sig("mean_sig","",0.0);
  RooRealVar sigma_sig("sigma_sig","",tmp_sigma_sig);
  RooRealVar mean_bkg("mean_bkg","",nEvents);
  RooRealVar sigma_bkg("sigma_bkg","",tmp_sigma_bkg);
  RooGaussian nsig_constraint("mu_constaint","mu_constraint",nsig,mean_sig,sigma_sig);
  RooGaussian nbkg_constraint("nbkg_constaint","nbkg_constraint",nbkg,mean_bkg,sigma_bkg);
  
  // --- Truth model ---
  RooRealVar b("b","-b",-0.1,-100000.0,0.) ; 
  RooExponential BkgTruthTmp("exp","",mes,b);
  RooRealVar nbkgTruth("nbkgTruth","",nEvents);
  RooExtendPdf BkgTruth("BkgTruth","",BkgTruthTmp,nbkgTruth);

  RooMCStudy * mcs = new RooMCStudy(BkgTruth, RooArgSet(mes), FitModel(sum),Silence(), Extended(kTRUE), Binned(kFALSE),
    				      FitOptions(Range(100.,180.), Save(kFALSE), SumW2Error(kTRUE)));

  RooChi2MCSModule chi2mod;
  mcs->addModule(chi2mod);

  mcs->generateAndFit(Npse,nEvents,kTRUE);

  std::vector<float> pulls;
  pulls.clear();
  
  for(int ii = 0; ii < Npse; ii++){
      if(ii > 995) continue;
      float pull;
      const RooAbsData* genDataset = mcs->genData(ii);
      float fitN = mcs->fitParams(ii)->getRealValue("nsig");
      float fitNerr = mcs->fitParams(ii)->getRealValue("nsigerr");
      pull = (0-fitN)/(fitNerr);
      h_pulls->Fill(pull);
      pulls.push_back(pull);

      std::cout << "TOY = " << ii << " , nsig = " << fitN << " , nsigErr = " << fitNerr << " , pull = "   << pull << " , nEvents = " << genDataset->sumEntries() << std::endl;

      /*RooPlot* xframe = x.frame(Title("Mgg"));
      data->plotOn(xframe); 
      fitFunk.plotOn(xframe);
      
      char name[100]; // enough to hold all numbers up to 64-bits
      sprintf(name, "plot_%d_%d",ii,nEvents);
      TCanvas c1("c1","c1");
      c1.cd();
      xframe->Draw();
      c1.Print((std::string(name)+".png").c_str(),"png");
      c1.Print((std::string(name)+".pdf").c_str(),"pdf");*/
  }
  
  float median =-9999.;
  float mean = 0.; 
  float Mean = -9999.; 

   for(unsigned int i=0; i<pulls.size(); ++i)
      if(fabs(pulls.at(i)))mean = mean + pulls.at(i);
  if(pulls.size()!= 0) Mean = mean/pulls.size();

  for(unsigned int i=0; i<pulls.size(); ++i){
      for(unsigned int j=i; j<pulls.size(); ++j){
	if(pulls[i]>pulls[j]){
	  float tmp=pulls[i];
	  pulls[i]=pulls[j];
	  pulls[j]=tmp;
	}}}
    if(pulls.size()==0)
      median = -9999;
    else if(pulls.size()%2==0)
      median = 0.5*(pulls[pulls.size()/2]+pulls[pulls.size()/2-1]);
    else
      median = pulls[pulls.size()/2];


  std::cout << "MEDIAN = " << median << std::endl;
  //std::cout << "MEAN = " << Mean << std::endl;

  TCanvas* c = new TCanvas("c","c");
  c->cd();
  h_pulls->Draw();
  c -> Print("pulls_1.png","png");
  c -> Print("pulls_1.pdf","pdf");
  f->cd();
  h_pulls->Write();
  f->Close();
}
