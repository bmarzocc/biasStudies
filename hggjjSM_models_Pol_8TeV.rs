mRad[320,1200];
mggjj_sig_m0[110.0, 70, 160];
mggjj_sig_sigma[1.5, 0.5, 2.0];
mggjj_sig_alpha[1.0, 0.5, 3]; 
mggjj_sig_n[4.0, 0.5, 10]; 
mggjj_sig_gsigma[2.5, 2.0, 5];
mggjj_sig_frac[0.3, 0, 0.5];

MggjjGaussSig = Gaussian(mRad, mggjj_sig_m0, mggjj_sig_gsigma);
MggjjCBSig    = CBShape(mRad, mggjj_sig_m0, mggjj_sig_sigma, mggjj_sig_alpha, mggjj_sig_n);
MggjjSig      = AddPdf(MggjjGaussSig, MggjjCBSig, mggjj_sig_frac);

mggjj_sig_m0_cat0[110.0, 105, 155];
mggjj_sig_sigma_cat0[1.2, 0.8, 2.5];
mggjj_sig_alpha_cat0[2.0, 1.0, 2.5]; 
mggjj_sig_n_cat0[2.0, 1.0, 5.0]; 
mggjj_sig_gsigma_cat0[5.0, 3.0, 8.0];
mggjj_sig_frac_cat0[0.1, 0, 0.4];


MggjjGaussSig_cat0 = Gaussian(mRad, mggjj_sig_m0_cat0, mggjj_sig_gsigma_cat0);
MggjjCBSig_cat0    = CBShape(mRad, mggjj_sig_m0_cat0, mggjj_sig_sigma_cat0, mggjj_sig_alpha_cat0, mggjj_sig_n_cat0);
MggjjSig_cat0      = AddPdf(MggjjGaussSig_cat0, MggjjCBSig_cat0, mggjj_sig_frac_cat0);

mggjj_sig_m0_cat1[110.0, 70, 160];
mggjj_sig_sigma_cat1[1.8, 1.0, 2.5];
mggjj_sig_alpha_cat1[2.0, 1.2, 5]; 
mggjj_sig_n_cat1[2.0, 1.5, 10]; 
mggjj_sig_gsigma_cat1[5.0, 3.0, 8.0];
mggjj_sig_frac_cat1[0.1, 0, 0.4];

MggjjGaussSig_cat1 = Gaussian(mRad, mggjj_sig_m0_cat1, mggjj_sig_gsigma_cat1);
MggjjCBSig_cat1    = CBShape(mRad, mggjj_sig_m0_cat1, mggjj_sig_sigma_cat1, mggjj_sig_alpha_cat1, mggjj_sig_n_cat1);
MggjjSig_cat1      = AddPdf(MggjjGaussSig_cat1, MggjjCBSig_cat1, mggjj_sig_frac_cat1);

mggjj_sig_m0_cat2[110.0, 70, 160];
mggjj_sig_sigma_cat2[2.0, 1.0, 3.5];
mggjj_sig_alpha_cat2[1.0, 1.0, 3.5]; 
mggjj_sig_n_cat2[2.0, 1.5, 10]; 
mggjj_sig_gsigma_cat2[5.0, 2.0, 8.0];
mggjj_sig_frac_cat2[0.1, 0, 0.3];

MggjjGaussSig_cat2 = Gaussian(mRad, mggjj_sig_m0_cat2, mggjj_sig_gsigma_cat2);
MggjjCBSig_cat2    = CBShape(mRad, mggjj_sig_m0_cat2, mggjj_sig_sigma_cat2, mggjj_sig_alpha_cat2, mggjj_sig_n_cat2);
MggjjSig_cat2      = AddPdf(MggjjGaussSig_cat2, MggjjCBSig_cat2, mggjj_sig_frac_cat2);

mggjj_sig_m0_cat3[110.0, 70, 160];
mggjj_sig_sigma_cat3[2.5, 2.0, 3.5];
mggjj_sig_alpha_cat3[1.0, 1.0, 3.5]; 
mggjj_sig_n_cat3[2.0, 1.5, 10]; 
mggjj_sig_gsigma_cat3[5.0, 3.5, 10.0];
mggjj_sig_frac_cat3[0.1, 0, 0.3];

MggjjGaussSig_cat3 = Gaussian(mRad, mggjj_sig_m0_cat3, mggjj_sig_gsigma_cat3);
MggjjCBSig_cat3    = CBShape(mRad, mggjj_sig_m0_cat3, mggjj_sig_sigma_cat3, mggjj_sig_alpha_cat3, mggjj_sig_n_cat3);
MggjjSig_cat3      = AddPdf(MggjjGaussSig_cat3, MggjjCBSig_cat3, mggjj_sig_frac_cat3);

mggjj_sig_m0_cat4[110.0, 70, 160];
mggjj_sig_sigma_cat4[2.5, 1.2, 3.0];
mggjj_sig_alpha_cat4[1.0, 1.0, 3.5]; 
mggjj_sig_n_cat4[2.0, 1.5, 10]; 
mggjj_sig_gsigma_cat4[5.0, 3.0, 10.0];
mggjj_sig_frac_cat4[0.1, 0, 0.3];

MggjjGaussSig_cat4 = Gaussian(mRad, mggjj_sig_m0_cat4, mggjj_sig_gsigma_cat4);
MggjjCBSig_cat4    = CBShape(mRad, mggjj_sig_m0_cat4, mggjj_sig_sigma_cat4, mggjj_sig_alpha_cat4, mggjj_sig_n_cat4);
MggjjSig_cat4      = AddPdf(MggjjGaussSig_cat4, MggjjCBSig_cat4, mggjj_sig_frac_cat4);


mggjj_sig_m0_cat5[110.0, 70, 160];
mggjj_sig_sigma_cat5[2.5, 1.2, 3.0];
mggjj_sig_alpha_cat5[1.0, 1.0, 3.5]; 
mggjj_sig_n_cat5[2.0, 1.5, 10]; 
mggjj_sig_gsigma_cat5[5.0, 3.0, 10.0];
mggjj_sig_frac_cat5[0.1, 0, 0.3];

MggjjGaussSig_cat5 = Gaussian(mRad, mggjj_sig_m0_cat5, mggjj_sig_gsigma_cat5);
MggjjCBSig_cat5    = CBShape(mRad, mggjj_sig_m0_cat5, mggjj_sig_sigma_cat5, mggjj_sig_alpha_cat5, mggjj_sig_n_cat5);
MggjjSig_cat5      = AddPdf(MggjjGaussSig_cat5, MggjjCBSig_cat5, mggjj_sig_frac_cat5);


mggjj_sig_m0_cat6[110.0, 70, 160];
mggjj_sig_sigma_cat6[2.5, 1.2, 3.0];
mggjj_sig_alpha_cat6[1.0, 1.0, 3.5]; 
mggjj_sig_n_cat6[2.0, 1.5, 10]; 
mggjj_sig_gsigma_cat6[5.0, 3.0, 10.0];
mggjj_sig_frac_cat6[0.1, 0, 0.3];

MggjjGaussSig_cat6 = Gaussian(mRad, mggjj_sig_m0_cat6, mggjj_sig_gsigma_cat6);
MggjjCBSig_cat6    = CBShape(mRad, mggjj_sig_m0_cat6, mggjj_sig_sigma_cat6, mggjj_sig_alpha_cat6, mggjj_sig_n_cat6);
MggjjSig_cat6      = AddPdf(MggjjGaussSig_cat6, MggjjCBSig_cat6, mggjj_sig_frac_cat6);


mggjj_sig_m0_cat7[110.0, 70, 160];
mggjj_sig_sigma_cat7[2.5, 1.2, 3.0];
mggjj_sig_alpha_cat7[1.0, 1.0, 3.5]; 
mggjj_sig_n_cat7[2.0, 1.5, 10]; 
mggjj_sig_gsigma_cat7[5.0, 3.0, 10.0];
mggjj_sig_frac_cat7[0.1, 0, 0.3];

MggjjGaussSig_cat7 = Gaussian(mRad, mggjj_sig_m0_cat7, mggjj_sig_gsigma_cat7);
MggjjCBSig_cat7    = CBShape(mRad, mggjj_sig_m0_cat7, mggjj_sig_sigma_cat7, mggjj_sig_alpha_cat7, mggjj_sig_n_cat7);
MggjjSig_cat7      = AddPdf(MggjjGaussSig_cat7, MggjjCBSig_cat7, mggjj_sig_frac_cat7);


mggjj_sig_m0_cat8[110.0, 70, 160];
mggjj_sig_sigma_cat8[2.5, 1.2, 3.0];
mggjj_sig_alpha_cat8[1.0, 1.0, 3.5]; 
mggjj_sig_n_cat8[2.0, 1.5, 10]; 
mggjj_sig_gsigma_cat8[5.0, 3.0, 10.0];
mggjj_sig_frac_cat8[0.1, 0, 0.3];

MggjjGaussSig_cat8 = Gaussian(mRad, mggjj_sig_m0_cat8, mggjj_sig_gsigma_cat8);
MggjjCBSig_cat8    = CBShape(mRad, mggjj_sig_m0_cat8, mggjj_sig_sigma_cat8, mggjj_sig_alpha_cat8, mggjj_sig_n_cat8);
MggjjSig_cat8      = AddPdf(MggjjGaussSig_cat8, MggjjCBSig_cat8, mggjj_sig_frac_cat8);



mggjj_bkg_8TeV_slope1[0.1,-100.0, 100.0];
mggjj_bkg_8TeV_slope2[0.1,-100.0, 100.0];
mggjj_bkg_8TeV_slope3[0.1,-100.0, 100.0];
mggjj_bkg_8TeV_slope4[0.1,-100.0, 100.0];
mggjj_bkg_8TeV_slope5[0.1,-100.0, 100.0];
mggjj_bkg_8TeV_arg1[0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg2[0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg3[0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg4[0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg5[0.1,-1000.0, 1000.0];

mggjj_bkg_8TeV_slope1_cat0[0.1,-100.0, 100.0];
mggjj_bkg_8TeV_slope2_cat0[0.1,-100.0, 100.0];
mggjj_bkg_8TeV_slope3_cat0[0.1,-100.0, 100.0];
mggjj_bkg_8TeV_slope4_cat0[0.1,-100.0, 100.0];
mggjj_bkg_8TeV_slope5_cat0[0.1,-100.0, 100.0];
mggjj_bkg_8TeV_arg1_cat0[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg2_cat0[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg3_cat0[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg4_cat0[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg5_cat0[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_wid1_cat0[250,0.0, 1000.0];
mggjj_bkg_8TeV_wid2_cat0[40,0.0, 1000.0];
mggjj_bkg_8TeV_wid3_cat0[1.0,0.0, 1000.0];
mggjj_bkg_8TeV_wid4_cat0[1.0,0.0, 1000.0];
mggjj_bkg_8TeV_wid5_cat0[0.1,0.0, 1000.0];

mggjj_bkg_8TeV_slope1_cat1[1.0,0, 10000.0];
mggjj_bkg_8TeV_slope2_cat1[1.0,0, 10000.0];
mggjj_bkg_8TeV_slope3_cat1[1.0,0, 10000.0];
mggjj_bkg_8TeV_slope4_cat1[1.0,0, 10000.0];
mggjj_bkg_8TeV_slope5_cat1[1.0,0, 10000.0];
mggjj_bkg_8TeV_arg1_cat1[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg2_cat1[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg3_cat1[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg4_cat1[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg5_cat1[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_wid1_cat1[250,0.0, 1000.0];
mggjj_bkg_8TeV_wid2_cat1[40,0.0, 1000.0];
mggjj_bkg_8TeV_wid3_cat1[1.0,0.0, 1000.0];
mggjj_bkg_8TeV_wid4_cat1[1.0,0.0, 1000.0];
mggjj_bkg_8TeV_wid5_cat1[0.1,0.0, 1000.0];

mggjj_bkg_8TeV_slope1_cat2[1.0,0, 1000.0];
mggjj_bkg_8TeV_slope2_cat2[1.0,0, 1000.0];
mggjj_bkg_8TeV_slope3_cat2[1.0,0, 1000.0];
mggjj_bkg_8TeV_slope4_cat2[1.0,0, 1000.0];
mggjj_bkg_8TeV_slope5_cat2[1.0,0, 1000.0];
mggjj_bkg_8TeV_arg1_cat2[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg2_cat2[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg3_cat2[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg4_cat2[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_arg5_cat2[-0.1,-1000.0, 1000.0];
mggjj_bkg_8TeV_wid1_cat2[250,0.0, 1000.0];
mggjj_bkg_8TeV_wid2_cat2[40,0.0, 1000.0];
mggjj_bkg_8TeV_wid3_cat2[1.0,0.0, 1000.0];
mggjj_bkg_8TeV_wid4_cat2[1.0,0.0, 1000.0];
mggjj_bkg_8TeV_wid5_cat2[0.1,0.0, 1000.0];

wei[1,0,10];
