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
#include "RooStats/HLFactory.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsData.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooProduct.h"
#include "RooExtendPdf.h"
#include "RooBernstein.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooCmdArg.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooExponential.h"

using namespace RooFit;
using namespace RooStats ;

int main(int argc, const char* argv[])
{
  int nEvents = atoi(argv[1]);
  int Npse = atoi(argv[2]);

  TFile* f = new TFile((std::string("output_")+std::string(argv[1])+".root").c_str(),"RECREATE");
  
  TH1F* h_pulls = new TH1F("h_pulls","h_pulls",20,-5,5);
  TH1F* h_Num = new TH1F("h_Num","h_Num",150,-3,3);
  TH1F* h_DeNum = new TH1F("h_DeNum","h_DeNum",200000,-100,100);
  TH1F* h_LogLike = new TH1F("h_LogLike","h_LogLike",100,0.,100.);

  TString card_name("models_2D.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  
  // --- Observable --- 
  RooRealVar mgg("mgg","mgg (GeV)",100.,180.) ; 

  // --- Truth model ---
  RooRealVar b("b","-b",-0.1,-100000.0,0.) ; 
  RooExponential BkgTruthTmp("exp","",mgg,b);
  RooRealVar nbkgTruth("nbkgTruth","",nEvents);
  RooExtendPdf BkgTruth("BkgTruth","",BkgTruthTmp,nbkgTruth);

  // --- Build Gaussian signal PDF --- 
  RooRealVar sigmean("sigmean","mgg mass",125.,123.,127.) ; 
  RooRealVar sigwidth("sigwidth","mgg width",1.7,0.7,2.7) ; 
  RooGaussian gauss("gauss","gaussian PDF",mgg,sigmean,sigwidth) ; 
  sigmean.setConstant(kTRUE);
  sigwidth.setConstant(kTRUE);

  // --- Build Expnential background PDF --- 
  RooRealVar a("a","-a",-0.1,-100000.0,0.) ; 
  RooExponential expo("exp","",mgg,a);

  // --- Construct signal+background PDF --- 
   RooRealVar nsig("nsig","#signal events",0.,-1*nEvents,1*nEvents) ; 
   RooRealVar nbkg("nbkg","#background events",nEvents,0.5*nEvents,1.5*nEvents) ; 
   RooAddPdf sum("sum","g+a",RooArgList(gauss,expo),RooArgList(nsig,nbkg)) ; 

  float aa[Npse] = {0};
  float bb[Npse] = {0};

  std::vector<float> pulls;
  pulls.clear();

  for(int i=0;i<Npse;i++){
   aa[i] = 0.0;
   bb[i] = 0.0;
   nsig.setVal(0);
   nbkg.setVal(nEvents);
   a.setVal(-0.1);

   RooArgSet* params = (RooArgSet*) sum.getParameters(mgg);
          std::cout << params->getSize() << std::endl;
          params->Print("v");
          TIterator* iter = params->createIterator();
          TObject* tempObj=0;
          while((tempObj=iter->Next())){
                 RooRealVar* var = (RooRealVar*)tempObj;
                 //var->setVal(0.1);
                 std::cout << "VALNAMEMgg = " << var->GetName() << " " << var->getVal() << std::endl;
          }
   
   TRandom rndm(0);
   int nEvent = rndm.Poisson(nEvents);
   nEvent = 0;
   //RooRandom::randomGenerator()->SetSeed(nEvent); 

   double expEvents = BkgTruth.expectedEvents(RooArgSet(mgg));
   std::cout << "expEvents = " << expEvents << std::endl;
   
   RooDataSet *data;
   if(nEvent !=0) data= BkgTruth.generate(mgg,nEvent,Silence(), Extended(kTRUE), Binned(kFALSE));
   else data= BkgTruth.generate(mgg,expEvents,Silence(), Extended(kTRUE), Binned(kFALSE));
   RooFitResult* fitResults = sum.fitTo(*data,Save(),SumW2Error(kTRUE)) ; 
   aa[i] = nsig.getVal();
   bb[i] = nsig.getError();
   if(bb[i] == 0) continue;
   if(fabs(fitResults->minNll()) > 10e6) continue;
   if(fabs(nsig.getVal())> 0.85*nEvents) continue;
   h_Num->Fill(aa[i]);
   h_DeNum->Fill(bb[i]);
   h_pulls->Fill(aa[i]/bb[i]);
   h_LogLike->Fill(fitResults->minNll());
   pulls.push_back(aa[i]/bb[i]);
   std::cout << "TOY = " << i << " , nsig = " << aa[i] << " , nsigError = " << bb[i] << " , -LogLike = " << fitResults->minNll() << " , sigmean = " << sigmean.getVal() << " , sigwidth = " << sigwidth.getVal() << std::endl;
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


  f->cd();
  h_pulls->Write();
  h_Num->Write();
  h_DeNum->Write();
  h_LogLike->Write();
  f->Close();
}
