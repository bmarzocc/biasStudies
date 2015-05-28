/*
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
#include "RooPower.h"
*/
using namespace std;
using namespace RooFit;
void plot_data_pull_v6_band()
{
  //  gROOT->ForceStyle();
  TFile* f = new TFile("data_had_v7.root","READ");
  TTree* tree = (TTree*)f->Get("HToGG_v1");

  RooRealVar mgg_cut_v1("mgg_cut_v1","mgg_cut_v1",100,180);
  RooDataSet ds("ds","ds",RooArgSet(mgg_cut_v1),Import(*tree) );
  ds.Print();
  int n_entry = ds.numEntries();
  cout<<n_entry<<endl;

  //signal shape

  RooRealVar signal_nuisance_NonLinearity_scale("signal_nuisance_NonLinearity_scale","signal_nuisance_NonLinearity_scale",-5,5);
  RooConstVar signal_const_mean_scale("signal_const_mean_scale","signal_const_mean_scale",0.001);
  RooRealVar signal_nuisance_mean_scale("signal_nuisance_mean_scale","signal_nuisance_mean_scale",-5,5);
  RooConstVar signal_const_sigma_scale("signal_const_sigma_scale","signal_const_sigma_scale",0.001);
  RooRealVar signal_nuisance_sigma_scale("signal_nuisance_sigma_scale","signal_nuisance_sigma_scale",-5,5);
  RooConstVar signal_const_sigma_smear("signal_const_sigma_smear","signal_const_sigma_smear",0.001);
  RooRealVar signal_nuisance_sigma_smear("signal_nuisance_sigma_smear","signal_nuisance_sigma_smear",-5,5);
  RooConstVar signal_const_mean_MaterialEBCentral_scale("signal_const_mean_MaterialEBCentral_scale","signal_const_mean_MaterialEBCentral_scale",0.0034);
  RooRealVar signal_nuisance_MaterialEBCentral_scale("signal_nuisance_MaterialEBCentral_scale","signal_nuisance_MaterialEBCentral_scale",-5,5);
  RooConstVar signal_const_mean_MaterialEBOuterEE_scale("signal_const_mean_MaterialEBOuterEE_scale","signal_const_mean_MaterialEBOuterEE_scale",0.0034);
  RooRealVar signal_nuisance_MaterialEBOuterEE_scale("signal_nuisance_MaterialEBOuterEE_scale","signal_nuisance_MaterialEBOuterEE_scale",-5,5);

  RooRealVar signal_gaussian_1_mean("signal_gaussian_1_mean","signal_gaussian_1_mean",1.24906e+02,1.24906e+02,1.24906e+02);
  RooRealVar signal_gaussian_2_mean("signal_gaussian_2_mean","signal_gaussian_2_mean",1.23877e+02,1.23877e+02,1.23877e+02);
  RooRealVar signal_gaussian_3_mean("signal_gaussian_3_mean","signal_gaussian_3_mean",1.19919e+02,1.19919e+02,1.19919e+02);
  RooRealVar signal_gaussian_1_sigma("signal_gaussian_1_sigma","signal_gaussian_1_sigma",1.03959e+00,1.03959e+00,1.03959e+00);
  RooRealVar signal_gaussian_2_sigma("signal_gaussian_2_sigma","signal_gaussian_2_sigma",2.41697e+00,2.41697e+00,2.41697e+00);
  RooRealVar signal_gaussian_3_sigma("signal_gaussian_3_sigma","signal_gaussian_3_sigma",7.60366e+00,7.60366e+00,7.60366e+00);

  RooFormulaVar signal_gauss_1_mean("signal_gauss_1_mean","signal_gauss_1_mean","@0*(1.+@1+@2*@3+@4*@5+@6*@7)",RooArgList(signal_gaussian_1_mean,signal_nuisance_NonLinearity_scale,signal_const_mean_scale,signal_nuisance_mean_scale,signal_const_mean_MaterialEBCentral_scale,signal_nuisance_MaterialEBCentral_scale,signal_const_mean_MaterialEBOuterEE_scale,signal_nuisance_MaterialEBOuterEE_scale));
  RooFormulaVar signal_gauss_2_mean("signal_gauss_2_mean","signal_gauss_2_mean","@0*(1.+@1+@2*@3+@4*@5+@6*@7)",RooArgList(signal_gaussian_2_mean,signal_nuisance_NonLinearity_scale,signal_const_mean_scale,signal_nuisance_mean_scale,signal_const_mean_MaterialEBCentral_scale,signal_nuisance_MaterialEBCentral_scale,signal_const_mean_MaterialEBOuterEE_scale,signal_nuisance_MaterialEBOuterEE_scale));
  RooFormulaVar signal_gauss_3_mean("signal_gauss_3_mean","signal_gauss_3_mean","@0*(1.+@1+@2*@3+@4*@5+@6*@7)",RooArgList(signal_gaussian_3_mean,signal_nuisance_NonLinearity_scale,signal_const_mean_scale,signal_nuisance_mean_scale,signal_const_mean_MaterialEBCentral_scale,signal_nuisance_MaterialEBCentral_scale,signal_const_mean_MaterialEBOuterEE_scale,signal_nuisance_MaterialEBOuterEE_scale));

  RooFormulaVar signal_gauss_1_sigma("signal_gauss_1_sigma","signal_gauss_1_sigma","TMath::Max(@0*(1.+TMath::Sqrt(0.+@1*@2*@1*@2+@3*@4*@3*@4+@5*@6*@5*@6+@7*@8*@7*@8)),1.e-6)",RooArgList(signal_gaussian_1_sigma,signal_const_sigma_scale,signal_nuisance_sigma_scale,signal_const_sigma_smear,signal_nuisance_sigma_smear,signal_const_mean_MaterialEBCentral_scale,signal_nuisance_MaterialEBCentral_scale,signal_const_mean_MaterialEBOuterEE_scale,signal_nuisance_MaterialEBOuterEE_scale));
  RooFormulaVar signal_gauss_2_sigma("signal_gauss_2_sigma","signal_gauss_2_sigma","TMath::Max(@0*(1.+TMath::Sqrt(0.+@1*@2*@1*@2+@3*@4*@3*@4+@5*@6*@5*@6+@7*@8*@7*@8)),1.e-6)",RooArgList(signal_gaussian_2_sigma,signal_const_sigma_scale,signal_nuisance_sigma_scale,signal_const_sigma_smear,signal_nuisance_sigma_smear,signal_const_mean_MaterialEBCentral_scale,signal_nuisance_MaterialEBCentral_scale,signal_const_mean_MaterialEBOuterEE_scale,signal_nuisance_MaterialEBOuterEE_scale));
  RooFormulaVar signal_gauss_3_sigma("signal_gauss_3_sigma","signal_gauss_3_sigma","TMath::Max(@0*(1.+TMath::Sqrt(0.+@1*@2*@1*@2+@3*@4*@3*@4+@5*@6*@5*@6+@7*@8*@7*@8)),1.e-6)",RooArgList(signal_gaussian_3_sigma,signal_const_sigma_scale,signal_nuisance_sigma_scale,signal_const_sigma_smear,signal_nuisance_sigma_smear,signal_const_mean_MaterialEBCentral_scale,signal_nuisance_MaterialEBCentral_scale,signal_const_mean_MaterialEBOuterEE_scale,signal_nuisance_MaterialEBOuterEE_scale));


  RooRealVar signal_ngs1("signal gs events 1","signal gs events 1",1.5,0,6.26);
  RooRealVar signal_ngs2("signal gs events 2","signal gs events 2",1.0,0,6.26);
  RooRealVar signal_ngs3("signal gs events 3","signal gs events 3",0.2,0,6.26);

  
     RooRealVar fraction1_sig("fraction1_sig","fraction1_sig",6.84794e-01,6.84794e-01,6.84794e-01);
     RooRealVar fraction2_sig("fraction2_sig","fraction2_sig",2.86962e-01,2.86962e-01,2.86962e-01);
//     RooRealVar fraction3_sig("fraction3_sig","fraction3_sig",0.1,0.0,1);
     

  RooGaussian MggGauss_Sig_1("MggGauss_Sig_1","MggGauss_Sig_1",mgg_cut_v1, signal_gaussian_1_mean, signal_gaussian_1_sigma); 
  RooGaussian MggGauss_Sig_2("MggGauss_Sig_2","MggGauss_Sig_2",mgg_cut_v1, signal_gaussian_2_mean, signal_gaussian_2_sigma); 
  RooGaussian MggGauss_Sig_3("MggGauss_Sig_3","MggGauss_Sig_3",mgg_cut_v1, signal_gaussian_3_mean, signal_gaussian_3_sigma); 

/*
   RooGaussian MggGauss_Sig_1("MggGauss_Sig_1","MggGauss_Sig_1",mgg_cut_v1, signal_gauss_1_mean, signal_gauss_1_sigma); 
   RooGaussian MggGauss_Sig_2("MggGauss_Sig_2","MggGauss_Sig_2",mgg_cut_v1, signal_gauss_2_mean, signal_gauss_2_sigma); 
   RooGaussian MggGauss_Sig_3("MggGauss_Sig_3","MggGauss_Sig_3",mgg_cut_v1, signal_gauss_3_mean, signal_gauss_3_sigma); 
*/
//  RooAddPdf MggSig("MggSig","MggSig",RooArgList(MggGauss_Sig_1,MggGauss_Sig_2,MggGauss_Sig_3),RooArgList(signal_ngs1,signal_ngs2,signal_ngs3));
    RooAddPdf MggSig("MggSig","MggSig",RooArgList(MggGauss_Sig_1,MggGauss_Sig_2,MggGauss_Sig_3),RooArgList(fraction1_sig,fraction2_sig));

  RooRealVar MggSig_norm("MggSig_norm","MggSig_norm",0.00005,0,6.26);
//  RooRealVar MggSig_norm("MggSig_norm","MggSig_norm",0.05,0,6.26);
//  RooRealVar MggSig_norm("MggSig_norm","MggSig_norm",0.6,0.6,0.6);
//  RooRealVar MggSig_norm_1("MggSig_norm_1","MggSig_norm_1",1.006,-2.02,16.26);
  RooRealVar MggSig_norm_1("MggSig_norm_1","MggSig_norm_1",0.0006,-2.3,9);
//   RooRealVar MggSig_norm("MggSig_norm","MggSig_norm",0,0,0);
//   RooRealVar MggSig_norm("MggSig_norm","MggSig_norm",1,0,3);

  //Higgs bkg shape

  RooRealVar higgs_bkg_nuisance_NonLinearity_scale("higgs_bkg_nuisance_NonLinearity_scale","higgs_bkg_nuisance_NonLinearity_scale",-5,5);
  RooConstVar higgs_bkg_const_mean_scale("higgs_bkg_const_mean_scale","higgs_bkg_const_mean_scale",0.001);
  RooRealVar higgs_bkg_nuisance_mean_scale("higgs_bkg_nuisance_mean_scale","higgs_bkg_nuisance_mean_scale",-5,5);
  RooConstVar higgs_bkg_const_sigma_scale("higgs_bkg_const_sigma_scale","higgs_bkg_const_sigma_scale",0.001);
  RooRealVar higgs_bkg_nuisance_sigma_scale("higgs_bkg_nuisance_sigma_scale","higgs_bkg_nuisance_sigma_scale",-5,5);
  RooConstVar higgs_bkg_const_sigma_smear("higgs_bkg_const_sigma_smear","higgs_bkg_const_sigma_smear",0.001);
  RooRealVar higgs_bkg_nuisance_sigma_smear("higgs_bkg_nuisance_sigma_smear","higgs_bkg_nuisance_sigma_smear",-5,5);
  RooConstVar higgs_bkg_const_mean_MaterialEBCentral_scale("higgs_bkg_const_mean_MaterialEBCentral_scale","higgs_bkg_const_mean_MaterialEBCentral_scale",0.0034);
  RooRealVar higgs_bkg_nuisance_MaterialEBCentral_scale("higgs_bkg_nuisance_MaterialEBCentral_scale","higgs_bkg_nuisance_MaterialEBCentral_scale",-5,5);
  RooConstVar higgs_bkg_const_mean_MaterialEBOuterEE_scale("higgs_bkg_const_mean_MaterialEBOuterEE_scale","higgs_bkg_const_mean_MaterialEBOuterEE_scale",0.0034);
  RooRealVar higgs_bkg_nuisance_MaterialEBOuterEE_scale("higgs_bkg_nuisance_MaterialEBOuterEE_scale","higgs_bkg_nuisance_MaterialEBOuterEE_scale",-5,5);

  RooRealVar higgs_bkg_gaussian_1_mean("higgs_bkg_gaussian_1_mean","higgs_bkg_gaussian_1_mean",1.24906e+02,1.24906e+02,1.24906e+02);
  RooRealVar higgs_bkg_gaussian_2_mean("higgs_bkg_gaussian_2_mean","higgs_bkg_gaussian_2_mean",1.23877e+02,1.23877e+02,1.23877e+02);
  RooRealVar higgs_bkg_gaussian_3_mean("higgs_bkg_gaussian_3_mean","higgs_bkg_gaussian_3_mean",1.19919e+02,1.19919e+02,1.19919e+02);
  RooRealVar higgs_bkg_gaussian_1_sigma("higgs_bkg_gaussian_1_sigma","higgs_bkg_gaussian_1_sigma",1.03959e+00,1.03959e+00,1.03959e+00);
  RooRealVar higgs_bkg_gaussian_2_sigma("higgs_bkg_gaussian_2_sigma","higgs_bkg_gaussian_2_sigma",2.41697e+00,2.41697e+00,2.41697e+00);
  RooRealVar higgs_bkg_gaussian_3_sigma("higgs_bkg_gaussian_3_sigma","higgs_bkg_gaussian_3_sigma",7.60366e+00,7.60366e+00,7.60366e+00);


  RooFormulaVar higgs_bkg_gauss_1_mean("higgs_bkg_gauss_1_mean","higgs_bkg_gauss_1_mean","@0*(1.+@1+@2*@3+@4*@5+@6*@7)",RooArgList(higgs_bkg_gaussian_1_mean,higgs_bkg_nuisance_NonLinearity_scale,higgs_bkg_const_mean_scale,higgs_bkg_nuisance_mean_scale,higgs_bkg_const_mean_MaterialEBCentral_scale,higgs_bkg_nuisance_MaterialEBCentral_scale,higgs_bkg_const_mean_MaterialEBOuterEE_scale,higgs_bkg_nuisance_MaterialEBOuterEE_scale));
  RooFormulaVar higgs_bkg_gauss_2_mean("higgs_bkg_gauss_2_mean","higgs_bkg_gauss_2_mean","@0*(1.+@1+@2*@3+@4*@5+@6*@7)",RooArgList(higgs_bkg_gaussian_2_mean,higgs_bkg_nuisance_NonLinearity_scale,higgs_bkg_const_mean_scale,higgs_bkg_nuisance_mean_scale,higgs_bkg_const_mean_MaterialEBCentral_scale,higgs_bkg_nuisance_MaterialEBCentral_scale,higgs_bkg_const_mean_MaterialEBOuterEE_scale,higgs_bkg_nuisance_MaterialEBOuterEE_scale));
  RooFormulaVar higgs_bkg_gauss_3_mean("higgs_bkg_gauss_3_mean","higgs_bkg_gauss_3_mean","@0*(1.+@1+@2*@3+@4*@5+@6*@7)",RooArgList(higgs_bkg_gaussian_3_mean,higgs_bkg_nuisance_NonLinearity_scale,higgs_bkg_const_mean_scale,higgs_bkg_nuisance_mean_scale,higgs_bkg_const_mean_MaterialEBCentral_scale,higgs_bkg_nuisance_MaterialEBCentral_scale,higgs_bkg_const_mean_MaterialEBOuterEE_scale,higgs_bkg_nuisance_MaterialEBOuterEE_scale));

  RooFormulaVar higgs_bkg_gauss_1_sigma("higgs_bkg_gauss_1_sigma","higgs_bkg_gauss_1_sigma","TMath::Max(@0*(1.+TMath::Sqrt(0.+@1*@2*@1*@2+@3*@4*@3*@4+@5*@6*@5*@6+@7*@8*@7*@8)),1.e-6)",RooArgList(higgs_bkg_gaussian_1_sigma,higgs_bkg_const_sigma_scale,higgs_bkg_nuisance_sigma_scale,higgs_bkg_const_sigma_smear,higgs_bkg_nuisance_sigma_smear,higgs_bkg_const_mean_MaterialEBCentral_scale,higgs_bkg_nuisance_MaterialEBCentral_scale,higgs_bkg_const_mean_MaterialEBOuterEE_scale,higgs_bkg_nuisance_MaterialEBOuterEE_scale));
  RooFormulaVar higgs_bkg_gauss_2_sigma("higgs_bkg_gauss_2_sigma","higgs_bkg_gauss_2_sigma","TMath::Max(@0*(1.+TMath::Sqrt(0.+@1*@2*@1*@2+@3*@4*@3*@4+@5*@6*@5*@6+@7*@8*@7*@8)),1.e-6)",RooArgList(higgs_bkg_gaussian_2_sigma,higgs_bkg_const_sigma_scale,higgs_bkg_nuisance_sigma_scale,higgs_bkg_const_sigma_smear,higgs_bkg_nuisance_sigma_smear,higgs_bkg_const_mean_MaterialEBCentral_scale,higgs_bkg_nuisance_MaterialEBCentral_scale,higgs_bkg_const_mean_MaterialEBOuterEE_scale,higgs_bkg_nuisance_MaterialEBOuterEE_scale));
  RooFormulaVar higgs_bkg_gauss_3_sigma("higgs_bkg_gauss_3_sigma","higgs_bkg_gauss_3_sigma","TMath::Max(@0*(1.+TMath::Sqrt(0.+@1*@2*@1*@2+@3*@4*@3*@4+@5*@6*@5*@6+@7*@8*@7*@8)),1.e-6)",RooArgList(higgs_bkg_gaussian_3_sigma,higgs_bkg_const_sigma_scale,higgs_bkg_nuisance_sigma_scale,higgs_bkg_const_sigma_smear,higgs_bkg_nuisance_sigma_smear,higgs_bkg_const_mean_MaterialEBCentral_scale,higgs_bkg_nuisance_MaterialEBCentral_scale,higgs_bkg_const_mean_MaterialEBOuterEE_scale,higgs_bkg_nuisance_MaterialEBOuterEE_scale));


  RooRealVar higgs_bkg_ngs1("higgs bkg gs events 1","higgs bkg gs events 1",0.08,0.01,0.15166);
  RooRealVar higgs_bkg_ngs2("higgs bkg gs events 2","higgs bkg gs events 2",0.06,0.01,0.15166);
  RooRealVar higgs_bkg_ngs3("higgs bkg gs events 3","higgs bkg gs events 3",0.04,0.01,0.15166);
  
     RooRealVar fraction1_bkg("fraction1_bkg","fraction1_bkg",6.84794e-01,6.84794e-01,6.84794e-01);
     RooRealVar fraction2_bkg("fraction2_bkg","fraction2_bkg",2.86962e-01,2.86962e-01,2.86962e-01);
//     RooRealVar fraction3_bkg("fraction3_bkg","fraction3_bkg",0.1,0.0,1);
     

  RooGaussian MggGauss_higgs_bkg_1("MggGauss_higgs_bkg_1","MggGauss_higgs_bkg_1",mgg_cut_v1, higgs_bkg_gaussian_1_mean, higgs_bkg_gaussian_1_sigma); 
  RooGaussian MggGauss_higgs_bkg_2("MggGauss_higgs_bkg_2","MggGauss_higgs_bkg_2",mgg_cut_v1, higgs_bkg_gaussian_2_mean, higgs_bkg_gaussian_2_sigma); 
  RooGaussian MggGauss_higgs_bkg_3("MggGauss_higgs_bkg_3","MggGauss_higgs_bkg_3",mgg_cut_v1, higgs_bkg_gaussian_3_mean, higgs_bkg_gaussian_3_sigma); 

/*
   RooGaussian MggGauss_higgs_bkg_1("MggGauss_higgs_bkg_1","MggGauss_higgs_bkg_1",mgg_cut_v1, higgs_bkg_gauss_1_mean, higgs_bkg_gauss_1_sigma); 
   RooGaussian MggGauss_higgs_bkg_2("MggGauss_higgs_bkg_2","MggGauss_higgs_bkg_2",mgg_cut_v1, higgs_bkg_gauss_2_mean, higgs_bkg_gauss_2_sigma); 
   RooGaussian MggGauss_higgs_bkg_3("MggGauss_higgs_bkg_3","MggGauss_higgs_bkg_3",mgg_cut_v1, higgs_bkg_gauss_3_mean, higgs_bkg_gauss_3_sigma); 
*/
//  RooAddPdf MggHiggsbkg("MggHiggsbkg","MggHiggsbkg",RooArgList(MggGauss_higgs_bkg_1,MggGauss_higgs_bkg_2,MggGauss_higgs_bkg_3),RooArgList(higgs_bkg_ngs1,higgs_bkg_ngs2,higgs_bkg_ngs3));
    RooAddPdf MggHiggsbkg("MggHiggsbkg","MggHiggsbkg",RooArgList(MggGauss_higgs_bkg_1,MggGauss_higgs_bkg_2,MggGauss_higgs_bkg_3),RooArgList(fraction1_bkg,fraction2_bkg));

  RooRealVar MggHiggsbkg_norm("MggHiggsbkg_norm","MggHiggsbkg_norm",0.15166,0.13166,0.17166);
//  RooRealVar MggHiggsbkg_norm("MggHiggsbkg_norm","MggHiggsbkg_norm",0.15166,0.15166,0.15166);
//  RooRealVar MggHiggsbkg_norm_1("MggHiggsbkg_norm_1","MggHiggsbkg_norm_1",0.15166,0.13166,0.17166);
  RooRealVar MggHiggsbkg_norm_1("MggHiggsbkg_norm_1","MggHiggsbkg_norm_1",0.15166,0.15166,0.15166);
//  RooRealVar MggHiggsbkg_norm_1("MggHiggsbkg_norm_1","MggHiggsbkg_norm_1",0.,0.,0.);



//non-resonant bkg shape


  
/*
     RooRealVar a0("a0","a0",0.1,-10,10) ;
     RooRealVar a1("a1","a1",0.1,-10,10) ;
     RooRealVar a2("a2","a2",0.1,-10,10) ;
     RooRealVar a3("a3","a3",0.1,-10,10) ;
  RooBernstein Mgg_Nonres_bkg("Mgg_Nonres_bkg", "Mgg_Nonres_bkg", mgg_cut_v1, RooArgList(RooConst(1.0),a0,a1,a2,a3));

     RooRealVar a0("a0","a0",9.5956e-01,9.5956e-01,9.5956e-01) ;
     RooRealVar a1("a1","a1",-5.9474e-02,-5.9474e-02,-5.9474e-02) ;
     RooRealVar a2("a2","a2",1.8817e-01,1.8817e-01,1.8817e-01) ;
     RooRealVar a3("a3","a3",9.0115e-02,9.0115e-02,9.0115e-02) ;


  RooRealVar a0("a0","a0",9.5956e-01,0.1,1.5) ;
  RooRealVar a1("a1","a1",-5.9474e-02,-0.01,0.1) ;
  RooRealVar a2("a2","a2",1.8817e-01,0.01,0.8) ;
  RooRealVar a3("a3","a3",9.0115e-02,0.01,0.1) ;
  RooBernstein Mgg_Nonres_bkg("Mgg_Nonres_bkg", "Mgg_Nonres_bkg", mgg_cut_v1, RooArgList(RooConst(1.0),a0,a1,a2,a3));


  RooRealVar a0_v1("a0_v1","a0_v1",9.5956e-01,0.1,1.5) ;
  RooRealVar a1_v1("a1_v1","a1_v1",-5.9474e-02,-0.01,0.1) ;
  RooRealVar a2_v1("a2_v1","a2_v1",1.8817e-01,0.01,0.8) ;
  RooRealVar a3_v1("a3_v1","a3_v1",9.0115e-02,0.01,0.1) ;
  RooBernstein Mgg_Nonres_bkg_v1("Mgg_Nonres_bkg_v1", "Mgg_Nonres_bkg_v1", mgg_cut_v1, RooArgList(RooConst(1.0),a0_v1,a1_v1,a2_v1,a3_v1));

*/
/*
  //fourth order Bernstein Polynomial
  RooRealVar a0("a0","a0",9.5956e-01,9.5956e-01 - 6.18e-02,9.5956e-01 + 6.18e-02) ;
  RooRealVar a1("a1","a1",-5.9474e-02,-5.9474e-02 - 4.78e-02,-5.9474e-02 + 4.78e-02) ;
  RooRealVar a2("a2","a2",1.8817e-01,1.8817e-01 - 2.82e-02,1.8817e-01 + 2.82e-02) ;
  RooRealVar a3("a3","a3",9.0115e-02,9.0115e-02 - 7.72e-03,9.0115e-02 + 7.72e-03) ;
  RooBernstein Mgg_Nonres_bkg("Mgg_Nonres_bkg", "Mgg_Nonres_bkg", mgg_cut_v1, RooArgList(RooConst(1.0),a0,a1,a2,a3));


  RooRealVar a0_v1("a0_v1","a0_v1",9.5956e-01,9.5956e-01 - 6.18e-02,9.5956e-01 + 6.18e-02) ;
  RooRealVar a1_v1("a1_v1","a1_v1",-5.9474e-02,-5.9474e-02 - 4.78e-02,-5.9474e-02 + 4.78e-02) ;
  RooRealVar a2_v1("a2_v1","a2_v1",1.8817e-01,1.8817e-01 - 2.82e-02,1.8817e-01 + 2.82e-02) ;
  RooRealVar a3_v1("a3_v1","a3_v1",9.0115e-02,9.0115e-02 - 7.72e-03,9.0115e-02 + 7.72e-03) ;
  RooBernstein Mgg_Nonres_bkg_v1("Mgg_Nonres_bkg_v1", "Mgg_Nonres_bkg_v1", mgg_cut_v1, RooArgList(RooConst(1.0),a0_v1,a1_v1,a2_v1,a3_v1));
*/



  //third order Bernstein Polynomial
  RooRealVar a0("a0","a0",4.4718e-01,4.4718e-01 - 2.86e-02,4.4718e-01 + 2.86e-02) ;
  RooRealVar a1("a1","a1",2.9732e-02,2.9732e-02 - 1.63e-02,2.9732e-02 + 1.63e-02) ;
  RooRealVar a2("a2","a2",9.6354e-02,9.6354e-02 - 6.72e-03,9.6354e-02 + 6.72e-03) ;
  RooBernstein Mgg_Nonres_bkg("Mgg_Nonres_bkg", "Mgg_Nonres_bkg", mgg_cut_v1, RooArgList(RooConst(1.0),a0,a1,a2));

  RooRealVar a0_v1("a0_v1","a0_v1",4.4718e-01,4.4718e-01 - 2.86e-02,4.4718e-01 + 2.86e-02) ;
  RooRealVar a1_v1("a1_v1","a1_v1",2.9732e-02,2.9732e-02 - 1.63e-02,2.9732e-02 + 1.63e-02) ;
  RooRealVar a2_v1("a2_v1","a2_v1",9.6354e-02,9.6354e-02 - 6.72e-03,9.6354e-02 + 6.72e-03) ;
  RooBernstein Mgg_Nonres_bkg_v1("Mgg_Nonres_bkg_v1", "Mgg_Nonres_bkg_v1", mgg_cut_v1, RooArgList(RooConst(1.0),a0_v1,a1_v1,a2_v1));

/*
  //second order Bernstein Polynomial
  RooRealVar a0("a0","a0",6.2484e-02,6.2484e-02 - 1.08e-02,6.2484e-02 + 1.08e-02) ;
  RooRealVar a1("a1","a1",8.6070e-02,8.6070e-02 - 1.08e-02,8.6070e-02 + 1.08e-02) ;
  RooBernstein Mgg_Nonres_bkg("Mgg_Nonres_bkg", "Mgg_Nonres_bkg", mgg_cut_v1, RooArgList(RooConst(1.0),a0,a1));

  RooRealVar a0_v1("a0_v1","a0_v1",6.2484e-02,6.2484e-02 - 1.08e-02,6.2484e-02 + 1.08e-02) ;
  RooRealVar a1_v1("a1_v1","a1_v1",8.6070e-02,8.6070e-02 - 1.08e-02,8.6070e-02 + 1.08e-02) ;
 RooBernstein Mgg_Nonres_bkg_v1("Mgg_Nonres_bkg_v1", "Mgg_Nonres_bkg_v1", mgg_cut_v1, RooArgList(RooConst(1.0),a0_v1,a1_v1));
  //first order exponential
  RooRealVar exp_lamda("exp_lamda","exp_lamda",-3.2609e-02,-3.2609e-02 - 3.57e-04,-3.2609e-02 + 3.57e-04) ;
  RooExponential Mgg_Nonres_bkg("Mgg_Nonres_bkg","",mgg_cut_v1,exp_lamda);

  RooRealVar exp_lamda_v1("exp_lamda_v1","exp_lamda_v1",-3.2609e-02,-3.2609e-02 - 3.57e-04,-3.2609e-02 + 3.57e-04) ;
  RooExponential Mgg_Nonres_bkg_v1("Mgg_Nonres_bkg_v1","",mgg_cut_v1,exp_lamda_v1);



  //Landau
  RooRealVar lan_mean("lan_mean","lan_mean",1.0778e+02,1.0778e+02 - 3.55e-01,1.0778e+02 + 3.55e-01) ;
  RooRealVar lan_sigma("lan_sigma","lan_sigma",8.5910e+00,8.5910e+00 - 1.19e-01,8.5910e+00 + 1.19e-01) ;
  RooLandau Mgg_Nonres_bkg("Mgg_Nonres_bkg","",mgg_cut_v1,lan_mean,lan_sigma);

  RooRealVar lan_mean_v1("lan_mean_v1","lan_mean_v1",1.0778e+02,1.0778e+02 - 3.55e-01,1.0778e+02 + 3.55e-01) ;
  RooRealVar lan_sigma_v1("lan_sigma_v1","lan_sigma_v1",8.5910e+00,8.5910e+00 - 1.19e-01,8.5910e+00 + 1.19e-01) ;
  RooLandau Mgg_Nonres_bkg_v1("Mgg_Nonres_bkg_v1","",mgg_cut_v1,lan_mean_v1,lan_sigma_v1);


   

  //first order power law
  RooRealVar pw("pw", "",-4.3903e+00 ,-4.3903e+00 - 4.58e-02,-4.3903e+00 + 4.58e-02 );
  RooGenericPdf Mgg_Nonres_bkg("Mgg_Nonres_bkg", "pow(@0,@1)", RooArgList(mgg_cut_v1, pw));

  RooRealVar pw_v1("pw_v1", "",-4.3903e+00 ,-4.3903e+00 - 4.58e-02,-4.3903e+00 + 4.58e-02 );
  RooGenericPdf Mgg_Nonres_bkg_v1("Mgg_Nonres_bkg_v1", "pow(@0,@1)", RooArgList(mgg_cut_v1, pw_v1));

*/

  RooRealVar Mgg_Nonres_bkg_norm("Mgg_Nonres_bkg_norm", "Mgg_Nonres_bkg_norm",20,-15,70);
  RooRealVar Mgg_Nonres_bkg_norm_1("Mgg_Nonres_bkg_norm_1", "Mgg_Nonres_bkg_norm_1",20,-15,70);


//signal + Higgs bkg + non-resonant bkg

//  RooAddPdf Mggall("Mggall","Mggall",RooArgList(Mgg_Nonres_bkg),RooArgList(Mgg_Nonres_bkg_norm));
//  RooAddPdf Mggall_bkg("Mggall_bkg","Mggall_bkg",RooArgList(MggSig,Mgg_Nonres_bkg_v1),RooArgList(MggSig_norm_1,Mgg_Nonres_bkg_norm_1));
//  RooAddPdf Mggall("Mggall","Mggall",RooArgList(MggHiggsbkg,Mgg_Nonres_bkg),RooArgList(MggHiggsbkg_norm,Mgg_Nonres_bkg_norm));
  RooAddPdf Mggall("Mggall","Mggall",RooArgList(MggSig,MggHiggsbkg,Mgg_Nonres_bkg),RooArgList(MggSig_norm,MggHiggsbkg_norm,Mgg_Nonres_bkg_norm));
//  RooAddPdf Mggall_bkg("Mggall_bkg","Mggall_bkg",RooArgList(MggSig,MggHiggsbkg,Mgg_Nonres_bkg),RooArgList(MggSig_norm,MggHiggsbkg_norm,Mgg_Nonres_bkg_norm));
  RooAddPdf Mggall_bkg("Mggall_bkg","Mggall_bkg",RooArgList(MggSig,MggHiggsbkg,Mgg_Nonres_bkg_v1),RooArgList(MggSig_norm_1,MggHiggsbkg_norm_1,Mgg_Nonres_bkg_norm_1));
//  RooAddPdf Mggall_bkg("Mggall_bkg","Mggall_bkg",RooArgList(MggHiggsbkg,Mgg_Nonres_bkg),RooArgList(MggHiggsbkg_norm,Mgg_Nonres_bkg_norm));

  RooFitResult* r_ml_wgt = Mggall.fitTo(ds,Save() ) ;

  TCanvas* c = new TCanvas("c1","c1",800,800) ;
  c->Divide(1,1) ;

  RooPlot* frame = mgg_cut_v1.frame(Title("Data"),Bins(20)) ;
  ds.plotOn(frame,DataError(RooAbsData::SumW2)) ;
  r_ml_wgt->Print();

    Mggall.plotOn(frame,Components(Mgg_Nonres_bkg),LineStyle(kDashed),LineColor(kBlue),LineWidth(4));


  RooAbsPdf *cpdf; cpdf = Mgg_Nonres_bkg;
  TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
  TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
 
  RooRealVar *nlim = new RooRealVar("nlim","",0.0,0.0,10.0);
  nlim->removeRange();

  RooCurve *nomcurve = dynamic_cast<RooCurve*>(frame->getObject(1));

      for (int i=1; i<(frame->GetXaxis()->GetNbins()+2); ++i) {
         double lowedge = frame->GetXaxis()->GetBinLowEdge(i);
         double upedge  = frame->GetXaxis()->GetBinUpEdge(i);
         double center  = frame->GetXaxis()->GetBinCenter(i);
//         double nombkg = nomcurve->interpolate(center);
         double nombkg = nomcurve->interpolate(lowedge);
         nlim->setVal(nombkg);

//         cout<<lowedge<<", "<<upedge<<", "<<center<<", "<<nombkg<<endl;

         mgg_cut_v1.setRange("errRange",lowedge,upedge);
         RooAbsPdf *epdf = 0;
         epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
 
         RooAbsReal *nll = epdf->createNLL(ds,Extended());
         RooMinimizer minim(*nll);
         minim.setStrategy(0);
         double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
         double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
 
         minim.migrad();
         minim.minos(*nlim);
         // printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
 
//         onesigma->SetPoint(i-1,center,nombkg);
         onesigma->SetPoint(i-1,lowedge,nombkg);
         onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
 
         minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
         // eventually if cl = 0.95 this is the usual 1.92!      
 
 
         minim.migrad();
         minim.minos(*nlim);
 
//         twosigma->SetPoint(i-1,center,nombkg);
         twosigma->SetPoint(i-1,lowedge,nombkg);
         twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
 
 
         delete nll;
         delete epdf;


         }



  //  Mggall.plotOn(frame,FillColor(kOrange),VisualizeError(*r_ml_wgt,1)) ;
  //  Mggall.plotOn(frame,VisualizeError(*r_ml_wgt,1,kTRUE),DrawOption("L"),LineWidth(2),LineColor(kBlue)) ;
//    Mggall.plotOn(frame,Components(Mgg_Nonres_bkg),LineStyle(kDashed),LineColor(kBlue),LineWidth(4));
//    Mggall.plotOn(frame,Components(Mgg_Nonres_bkg),FillColor(kOrange),VisualizeError(*r_ml_wgt,1)) ;
    Mggall.plotOn(frame,Components("MggHiggsbkg,Mgg_Nonres_bkg"),LineStyle(1),LineColor(kBlue),LineWidth(4));
  Mggall.plotOn(frame,LineColor(kRed),LineWidth(1));
//  ds.plotOn(frame) ;
  //  Mggall.plotOn(frame,Components(MggBkgAll),LineStyle(kDashed),LineColor(kRed));
  //  Mggall.plotOn(frame,Components(MggBkgAll),LineColor(kRed));
  //  Mggall.plotOn(frame,Components(MggHiggsbkg),LineColor(kRed));
  //  Mggall.plotOn(frame,Components(MggSig),LineColor(kGreen));
  //  Mggall.plotOn(frame,Components(MggSig),LineColor(kGreen));

       onesigma->SetLineColor(kGreen);
       onesigma->SetFillColor(kGreen);
       onesigma->SetMarkerColor(kGreen);
 
       twosigma->SetLineColor(kYellow);
       twosigma->SetFillColor(kYellow);
       twosigma->SetMarkerColor(kYellow);


  c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.2) ;frame->Draw() ;

       twosigma->Draw("L3 ");
       onesigma->Draw("L3 ");
       frame->Draw("SAME") ;


  TLegend *legmc = new TLegend(0.32,0.6,0.85,0.8);
  legmc->AddEntry(frame->getObject(7),"hadronic channel","");
  legmc->AddEntry(frame->getObject(0),"Data","LPE");
  legmc->AddEntry(frame->getObject(3),"Signal + total background fit","L");
  //  legmc->AddEntry(frame->getObject(3),"Bkg. + SM Higgs Bkg. + Sig.","L");
  legmc->AddEntry(frame->getObject(2),"Total background","L");
  legmc->AddEntry(frame->getObject(1),"Non-resonant diphoton background","L");

  legmc->SetBorderSize(0);
  legmc->SetTextSize(0.03);
  legmc->SetFillStyle(0);
  legmc->Draw("same");

  TLegend* leg1 = new TLegend(0.2,0.82,0.85,0.89);
  leg1->SetHeader("CMS Preliminary,   #sqrt{s} = 8 TeV,   #int Ldt = 19.7 fb^{-1}");
  leg1->SetFillColor(0);
  leg1->SetLineColor(0);
  leg1->SetBorderSize( 0);
  leg1->SetTextSize(0.03);
  leg1->Draw("same");

  frame->SetTitle("");
  frame->SetMinimum(0.0);
  frame->SetMaximum(1.10*frame->GetMaximum());
  frame->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
  frame->GetYaxis()->SetTitle("Events / 4 GeV");


     mgg_cut_v1.setRange("ttt",123,127);
     RooAbsReal* igx_sig = Mgg_Nonres_bkg.createIntegral(mgg_cut_v1,NormSet(mgg_cut_v1),Range("ttt"));
     cout<<igx_sig->getVal()*28.9 <<endl ;

/*  

     higgs_bkg_gaussian_1_mean.setConstant(kTRUE);
     higgs_bkg_gaussian_2_mean.setConstant(kTRUE);
     higgs_bkg_gaussian_3_mean.setConstant(kTRUE);
     higgs_bkg_gaussian_1_sigma.setConstant(kTRUE);
     higgs_bkg_gaussian_2_sigma.setConstant(kTRUE);
     higgs_bkg_gaussian_3_sigma.setConstant(kTRUE);
     fraction1_bkg.setConstant(kTRUE);
     fraction2_bkg.setConstant(kTRUE);
//     higgs_bkg_ngs1.setConstant(kTRUE);
//     higgs_bkg_ngs2.setConstant(kTRUE);
//     higgs_bkg_ngs3.setConstant(kTRUE);
     a0.setConstant(kTRUE);
     a1.setConstant(kTRUE);
     a2.setConstant(kTRUE);
//     a3.setConstant(kTRUE);
     MggHiggsbkg_norm.setConstant(kTRUE);
//     Mgg_Nonres_bkg_norm.setConstant(kTRUE);


     RooWorkspace *workspace_bkg = new RooWorkspace("workspace_bkg","workspace_bkg");
  //  workspace_bkg->import(MggNonresbkgpdf);
  workspace_bkg->import(MggHiggsbkg);
  workspace_bkg->import(Mgg_Nonres_bkg);
  workspace_bkg->import(MggHiggsbkg_norm);
  workspace_bkg->import(Mgg_Nonres_bkg_norm);
  workspace_bkg->writeToFile("input_bkg_had_v8.root");

*/

/*
  float aa[1000] = {0};
  float bb[1000] = {0};

  for(int i=0;i<1000;i++){
   aa[i] = 0.0;
   bb[i] = 0.0;

   RooDataSet *data = Mggall.generate(mgg_cut_v1,29);
   Mggall_bkg.fitTo(*data) ; 

  RooPlot* frame222 = mgg_cut_v1.frame(Title(""),Bins(20)) ;
  data->plotOn(frame222);
  Mggall_bkg.plotOn(frame222);

//  TCanvas* c4 = new TCanvas("c4","c4",800,600) ;
//  c4->cd() ; gPad->SetLeftMargin(0.15) ; frame222->GetYaxis()->SetTitleOffset(1.1) ; frame222->Draw() ;

//     RooAbsReal* test = MggSig_norm_1.createIntegral(mgg_cut_v1,NormSet(mgg_cut_v1),Range("ttt"));
//     cout<<MggSig_norm_1.getVal() <<","<<MggSig_norm_1.getError()<<endl ;
     aa[i] = MggSig_norm_1.getVal();
     bb[i] = MggSig_norm_1.getError();
//     aa[i] = Mgg_Nonres_bkg_norm_1.getVal();
//     bb[i] = Mgg_Nonres_bkg_norm_1.getError();
  }

  TH1F* h = new TH1F("h",";pull;events",24,-3,3);
  
  for(int i=0;i<1000;i++){
  h->Fill((aa[i]-0.6)/bb[i]);
//  cout<<aa[i]<<","<<bb[i]<<endl;
//  cout<<"----------"<<endl;
  }
  
  TCanvas* c2 = new TCanvas("c2","c2",800,600) ;
  c2->cd();
  h->Draw("PE");
*/
 
/*
//  RooMCStudy* mcstudy = new RooMCStudy(Mggall,Mggall_bkg,mgg_cut_v1,"","mhv");

//  RooMCStudy* mcstudy = new RooMCStudy(Mggall,mgg_cut_v1,Binned(kTRUE),Silence(),Extended(),FitModel(Mggall_bkg),
  RooMCStudy* mcstudy = new RooMCStudy(Mggall,mgg_cut_v1,Binned(kTRUE),Silence(),Extended(),
      FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;

  // G e n e r a t e   a n d   f i t   e v e n t s
  // ---------------------------------------------

  // Generate and fit 1000 samples of Poisson(nExpected) events
  mcstudy->generateAndFit(1000,30) ;
//  mcstudy->generate(1000,28) ;
//  mcstudy->fit(3,);



  // E x p l o r e   r e s u l t s   o f   s t u d y 
  // ------------------------------------------------

  // Make plots of the distributions of mean, the error on mean and the pull of mean
  RooPlot* frame1 = mcstudy->plotParam(Mgg_Nonres_bkg_norm,Bins(25)) ;
  RooPlot* frame2 = mcstudy->plotError(Mgg_Nonres_bkg_norm,Bins(25)) ;
//  RooPlot* frame3 = mcstudy->plotPull(MggSig_norm,Bins(25)) ;
  RooPlot* frame3 = mcstudy->plotPull(Mgg_Nonres_bkg_norm,Bins(25),FitGauss(kTRUE)) ;
//  RooPlot* frame1 = mcstudy->plotParam(MggSig_norm,Bins(25)) ;
//  RooPlot* frame2 = mcstudy->plotError(MggSig_norm,Bins(25)) ;
  frame1->SetTitle("");
  frame1->GetXaxis()->SetTitle("Number of non-resonant diphoton background Pull");
  frame2->SetTitle("");
  frame2->GetXaxis()->SetTitle("Number of non-resonant diphoton background Pull");
  frame3->SetTitle("");
  frame3->GetXaxis()->SetTitle("Number of non-resonant diphoton background Pull");

  // Plot distribution of minimized likelihood
  RooPlot* frame4 = mcstudy->plotNLL(Bins(30)) ;
  frame4->SetTitle("");

  // Make some histograms from the parameter dataset
  //  TH1* hh_cor_a0_s1f = mcstudy->fitParDataSet().createHistogram("hh",a1,YVar(sig1frac)) ;
  //  TH1* hh_cor_a0_a1  = mcstudy->fitParDataSet().createHistogram("hh",a0,YVar(a1)) ;

  // Access some of the saved fit results from individual toys
  //  TH2* corrHist000 = mcstudy->fitResult(0)->correlationHist("c000") ;
  //  TH2* corrHist127 = mcstudy->fitResult(127)->correlationHist("c127") ;
  //  TH2* corrHist953 = mcstudy->fitResult(953)->correlationHist("c953") ;



  // Draw all plots on a canvas
  //  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  TCanvas* c2 = new TCanvas("c2","c2",800,600) ;
  TCanvas* c3 = new TCanvas("c3","c3",800,600) ;
  TCanvas* c4 = new TCanvas("c4","c4",800,600) ;
  TCanvas* c5 = new TCanvas("c5","c5",800,600) ;
  c2->cd() ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.1) ; frame1->Draw() ;
  c3->cd() ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.1) ; frame2->Draw() ;
  c4->cd() ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.1) ; frame3->Draw() ;
  c5->cd() ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.1) ; frame4->Draw() ;

  TLine* l1 = new TLine(53.7197,0,53.7197,140);
  l1->SetLineWidth(2);
  l1->Draw("same");

*/


}

