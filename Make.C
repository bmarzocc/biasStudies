int Make(int mode, int cat)
{
  /* find location of include file by issuing
     scramv1 tool info roofitcore
   */

  TString incpath = gSystem->GetIncludePath();
  incpath.Append(" -I/afs/cern.ch/cms/slc5_amd64_gcc472/lcg/roofit/5.34.04-cms2/include");
  incpath.Append(" -L/afs/cern.ch/cms/slc5_amd64_gcc472/lcg/roofit/5.34.04-cms2/lib");
  gSystem->SetIncludePath(incpath.Data());

  //gSystem->Load("libRooFit") ;
  gSystem->Load("libRooFitCore") ;
  gSystem->Load("libRooStats") ;
  using namespace RooFit ;

  //gROOT->ProcessLine(".L HiggsCSandWidth.cc++");
  gSystem->Load("HiggsCSandWidth_cc.so");
  gStyle->SetStripDecimals(kFALSE);

  gROOT->ProcessLine(".L GaussExp.cxx+");
  gROOT->ProcessLine(".L ExpGaussExp.cxx+");

  if(mode==2) gROOT->ProcessLine(".L R2GGBBBiasStudy_2D.cc+");
  if(mode==1) gROOT->ProcessLine(".L R2GGBBBiasStudy_mgg.cc+");
  if(mode==0) gROOT->ProcessLine(".L R2GGBBBiasStudy_mggjj.cc+");

  if(mode!=2)
    runfits(cat);
  else{
    for( int i=0; i<1; i++)
      for ( int j=0; j<1; j++)
	runfits(cat,i,j);
  }

  return 0;

}
