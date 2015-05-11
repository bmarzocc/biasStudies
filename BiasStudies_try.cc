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
  
  TH1F* h_pulls = new TH1F("h_pulls","h_pulls",1000,-1,1);
  
  // --- Observable --- 
  RooRealVar mes("mes","m_{ES} (GeV)",100.,180.) ; 
  // --- Build Gaussian signal PDF --- 
  RooRealVar sigmean("sigmean","mgg mass",125.,123.,127.) ; 
  RooRealVar sigwidth("sigwidth","mgg width",1.7,0.7,2.7) ; 
  RooGaussian gauss("gauss","gaussian PDF",mes,sigmean,sigwidth) ; 
  // --- Build Expnential background PDF --- 
  RooRealVar a("a","-a",-0.1,-100000.0,0.) ; 
  RooExponential expo("exp","",mes,a);
  RooRealVar par1("par1","par1",0.1,0.,10000.) ; 
  RooRealVar par2("par2","par2",0.1,0.,10000.) ; 
  RooRealVar par3("par3","par3",0.1,0.,10000.) ; 
  RooRealVar par4("par4","par4",0.1,0.,10000.) ; 
  RooRealVar par5("par5","par5",0.1,0.,10000.) ; 
  RooBernstein bern1("BerN1Mgg", "",mes, RooArgList(par1,par2));
  RooBernstein bern4("BerN4Mgg", "",mes, RooArgList(par1,par2,par3,par4,par5));
  // --- Construct signal+background PDF --- 
  RooRealVar nsig("nsig","#signal events",0.,-1.*nEvents,1.*nEvents) ; 
  RooRealVar nbkg("nbkg","#background events",nEvents,0.5*nEvents,1.5*nEvents) ; 
  //RooAddPdf sum("sum","g+a",RooArgList(gauss,expo),RooArgList(nsig,nbkg)) ; 
  RooAddPdf sum("sum","g+a",RooArgList(bern1,expo),RooArgList(nsig,nbkg)) ; 
  //RooAddPdf sum("sum","g+a",RooArgList(bern4,expo),RooArgList(nsig,nbkg)) ; 
  
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
  /*RooRealVar nsig_truth("nsig_truth","#signal events",0.1*6./80.*nEvents,-1.*nEvents,1.*nEvents) ; 
  RooRealVar nbkg_truth("nbkg_truth","#background events",nEvents,0.5*nEvents,1.5*nEvents) ; 
  RooAddPdf BkgTruth("BkgTruth","",RooArgList(gauss,BkgTruthTmp),RooArgList(nsig_truth,nbkg_truth)) ; */

  /*RooMCStudy * mcs = new RooMCStudy(BkgTruth, RooArgSet(mes), FitModel(sum),Silence(), Extended(kTRUE), Binned(kFALSE),
    				      FitOptions(Range(100.,180.), Save(kFALSE), SumW2Error(kTRUE),
    						 ExternalConstraints(RooArgSet(nsig_constraint,nbkg_constraint )) ));

  RooChi2MCSModule chi2mod;
  mcs->addModule(chi2mod);

  mcs->generateAndFit(Npse,nEvents,kTRUE);*/

  std::vector<float> pulls;
  pulls.clear();
  
  for(int ii = 0; ii < Npse; ii++){
      /*float pull;
      const RooAbsData* genDataset = mcs->genData(ii);
      float fitN = mcs->fitParams(ii)->getRealValue("nsig");
      float fitNerr = mcs->fitParams(ii)->getRealValue("nsigerr");
      pull = (0-fitN)/(fitNerr);
      h_pulls->Fill(pull);
      pulls.push_back(pull);

      std::cout << "TOY = " << ii << " , nsig = " << fitN << " , nsigErr = " << fitNerr << " , pull = "   << pull << " , nEvents = " << genDataset->sumEntries() << std::endl;*/
      TRandom Poisson(0);
      int nEventsSmear = Poisson.Poisson(nEvents);
      float pull;
      // --- Generate a toyMC sample from composite PDF --- 
      RooDataSet *data = BkgTruth.generate(mes,nEventsSmear,Silence(), Extended(kTRUE), Binned(kFALSE)) ; 
      // --- Perform extended ML fit of composite PDF to toy data --- 
      sum.fitTo(*data,Verbose(kFALSE),Range(100.,180.), Silence(),Save(kFALSE), SumW2Error(kTRUE),ExternalConstraints(RooArgSet(nsig_constraint,nbkg_constraint ))) ; 
      pull = nsig.getVal()/nsig.getError();
      h_pulls->Fill(pull);
      pulls.push_back(pull);

      std::cout << "TOY = " << ii << " , nsig = " << nsig.getVal() << " , nsigErr = " << nsig.getError() << " , pull = "   << pull << " , nEvents = " << nEventsSmear << std::endl;

      /*RooPlot* xframe = x.frame(Title("Mgg"));
      data->plotOn(xframe); 
      fitFunk.plotOn(xframe);*/

      /*RooPlot* xframe = mes.frame() ; 
      data->plotOn(xframe) ; 
      sum.plotOn(xframe) ; 
      //sum.plotOn(xframe,Components(expo),LineStyle(kDashed)) ; 
      
      char name[100]; // enough to hold all numbers up to 64-bits
      sprintf(name, "plot_%d_%d",ii,nEvents);
      TCanvas c1("c1","c1");
      c1.cd();
      xframe->Draw();
      c1.Print((std::string(name)+".png").c_str(),"png");
      c1.Print((std::string(name)+".pdf").c_str(),"pdf");*/
  }
  
  float median =-9999.;

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

  TCanvas* c = new TCanvas("c","c");
  c->cd();
  h_pulls->Draw();
  c -> Print("pulls_1.png","png");
  c -> Print("pulls_1.pdf","pdf");
  f->cd();
  h_pulls->Write();
  f->Close();
}
