mJJ[60,180];
mjj_sig_m0[110.0, 70, 160];
mjj_sig_sigma[1.5, 0.5, 2.0];
mjj_sig_alpha[1.0, 0.5, 3]; 
mjj_sig_n[4.0, 0.5, 10]; 
mjj_sig_gsigma[2.5, 2.0, 5];
mjj_sig_frac[0.3, 0, 0.5];

MjjGaussSig = Gaussian(mJJ, mjj_sig_m0, mjj_sig_gsigma);
MjjCBSig    = CBShape(mJJ, mjj_sig_m0, mjj_sig_sigma, mjj_sig_alpha, mjj_sig_n);
MjjSig      = AddPdf(MjjGaussSig, MjjCBSig, mjj_sig_frac);

mjj_sig_m0_cat0[110.0, 105, 155];
mjj_sig_sigma_cat0[11.0, 10.0, 30.0];
mjj_sig_alpha_cat0[0.5, 0.1, 5.0]; 
mjj_sig_n_cat0[20.0, 1.0, 145.0]; 
mjj_sig_gsigma_cat0[40.0, 30.0, 2000.0];
mjj_sig_frac_cat0[0.4, 0, 0.8];


MjjGaussSig_cat0 = Gaussian(mJJ, mjj_sig_m0_cat0, mjj_sig_gsigma_cat0);
MjjCBSig_cat0    = CBShape(mJJ, mjj_sig_m0_cat0, mjj_sig_sigma_cat0, mjj_sig_alpha_cat0, mjj_sig_n_cat0);
MjjSig_cat0      = AddPdf(MjjGaussSig_cat0, MjjCBSig_cat0, mjj_sig_frac_cat0);

mjj_sig_m0_cat1[110.0, 70, 160];
mjj_sig_sigma_cat1[10.0, 5.0, 40.0];
mjj_sig_alpha_cat1[2.0, 0.1, 10.0]; 
mjj_sig_n_cat1[2.0, 1.0, 20.0]; 
mjj_sig_gsigma_cat1[70.0, 20.0, 120.0];
mjj_sig_frac_cat1[0.1, 0, 0.6];

MjjGaussSig_cat1 = Gaussian(mJJ, mjj_sig_m0_cat1, mjj_sig_gsigma_cat1);
MjjCBSig_cat1    = CBShape(mJJ, mjj_sig_m0_cat1, mjj_sig_sigma_cat1, mjj_sig_alpha_cat1, mjj_sig_n_cat1);
MjjSig_cat1      = AddPdf(MjjGaussSig_cat1, MjjCBSig_cat1, mjj_sig_frac_cat1);

mjj_sig_m0_cat2[110.0, 70, 160];
mjj_sig_sigma_cat2[2.0, 1.5, 40.0];
mjj_sig_alpha_cat2[2.0, 0.1, 10.0]; 
mjj_sig_n_cat2[2.0, 1.0, 20]; 
mjj_sig_gsigma_cat2[70.0, 30.0, 120.0];
mjj_sig_frac_cat2[0.1, 0, 0.4];

MjjGaussSig_cat2 = Gaussian(mJJ, mjj_sig_m0_cat2, mjj_sig_gsigma_cat2);
MjjCBSig_cat2    = CBShape(mJJ, mjj_sig_m0_cat2, mjj_sig_sigma_cat2, mjj_sig_alpha_cat2, mjj_sig_n_cat2);
MjjSig_cat2      = AddPdf(MjjGaussSig_cat2, MjjCBSig_cat2, mjj_sig_frac_cat2);

mjj_sig_m0_cat3[110.0, 70, 160];
mjj_sig_sigma_cat3[2.5, 2.0, 3.5];
mjj_sig_alpha_cat3[1.0, 1.0, 3.5]; 
mjj_sig_n_cat3[2.0, 1.5, 10]; 
mjj_sig_gsigma_cat3[5.0, 3.5, 10.0];
mjj_sig_frac_cat3[0.1, 0, 0.3];

MjjGaussSig_cat3 = Gaussian(mJJ, mjj_sig_m0_cat3, mjj_sig_gsigma_cat3);
MjjCBSig_cat3    = CBShape(mJJ, mjj_sig_m0_cat3, mjj_sig_sigma_cat3, mjj_sig_alpha_cat3, mjj_sig_n_cat3);
MjjSig_cat3      = AddPdf(MjjGaussSig_cat3, MjjCBSig_cat3, mjj_sig_frac_cat3);

mjj_sig_m0_cat4[110.0, 70, 160];
mjj_sig_sigma_cat4[2.5, 1.2, 3.0];
mjj_sig_alpha_cat4[1.0, 1.0, 3.5]; 
mjj_sig_n_cat4[2.0, 1.5, 10]; 
mjj_sig_gsigma_cat4[5.0, 3.0, 10.0];
mjj_sig_frac_cat4[0.1, 0, 0.3];

MjjGaussSig_cat4 = Gaussian(mJJ, mjj_sig_m0_cat4, mjj_sig_gsigma_cat4);
MjjCBSig_cat4    = CBShape(mJJ, mjj_sig_m0_cat4, mjj_sig_sigma_cat4, mjj_sig_alpha_cat4, mjj_sig_n_cat4);
MjjSig_cat4      = AddPdf(MjjGaussSig_cat4, MjjCBSig_cat4, mjj_sig_frac_cat4);


mjj_sig_m0_cat5[110.0, 70, 160];
mjj_sig_sigma_cat5[2.5, 1.2, 3.0];
mjj_sig_alpha_cat5[1.0, 1.0, 3.5]; 
mjj_sig_n_cat5[2.0, 1.5, 10]; 
mjj_sig_gsigma_cat5[5.0, 3.0, 10.0];
mjj_sig_frac_cat5[0.1, 0, 0.3];

MjjGaussSig_cat5 = Gaussian(mJJ, mjj_sig_m0_cat5, mjj_sig_gsigma_cat5);
MjjCBSig_cat5    = CBShape(mJJ, mjj_sig_m0_cat5, mjj_sig_sigma_cat5, mjj_sig_alpha_cat5, mjj_sig_n_cat5);
MjjSig_cat5      = AddPdf(MjjGaussSig_cat5, MjjCBSig_cat5, mjj_sig_frac_cat5);


mjj_sig_m0_cat6[110.0, 70, 160];
mjj_sig_sigma_cat6[2.5, 1.2, 3.0];
mjj_sig_alpha_cat6[1.0, 1.0, 3.5]; 
mjj_sig_n_cat6[2.0, 1.5, 10]; 
mjj_sig_gsigma_cat6[5.0, 3.0, 10.0];
mjj_sig_frac_cat6[0.1, 0, 0.3];

MjjGaussSig_cat6 = Gaussian(mJJ, mjj_sig_m0_cat6, mjj_sig_gsigma_cat6);
MjjCBSig_cat6    = CBShape(mJJ, mjj_sig_m0_cat6, mjj_sig_sigma_cat6, mjj_sig_alpha_cat6, mjj_sig_n_cat6);
MjjSig_cat6      = AddPdf(MjjGaussSig_cat6, MjjCBSig_cat6, mjj_sig_frac_cat6);


mjj_sig_m0_cat7[110.0, 70, 160];
mjj_sig_sigma_cat7[2.5, 1.2, 3.0];
mjj_sig_alpha_cat7[1.0, 1.0, 3.5]; 
mjj_sig_n_cat7[2.0, 1.5, 10]; 
mjj_sig_gsigma_cat7[5.0, 3.0, 10.0];
mjj_sig_frac_cat7[0.1, 0, 0.3];

MjjGaussSig_cat7 = Gaussian(mJJ, mjj_sig_m0_cat7, mjj_sig_gsigma_cat7);
MjjCBSig_cat7    = CBShape(mJJ, mjj_sig_m0_cat7, mjj_sig_sigma_cat7, mjj_sig_alpha_cat7, mjj_sig_n_cat7);
MjjSig_cat7      = AddPdf(MjjGaussSig_cat7, MjjCBSig_cat7, mjj_sig_frac_cat7);


mjj_sig_m0_cat8[110.0, 70, 160];
mjj_sig_sigma_cat8[2.5, 1.2, 3.0];
mjj_sig_alpha_cat8[1.0, 1.0, 3.5]; 
mjj_sig_n_cat8[2.0, 1.5, 10]; 
mjj_sig_gsigma_cat8[5.0, 3.0, 10.0];
mjj_sig_frac_cat8[0.1, 0, 0.3];

MjjGaussSig_cat8 = Gaussian(mJJ, mjj_sig_m0_cat8, mjj_sig_gsigma_cat8);
MjjCBSig_cat8    = CBShape(mJJ, mjj_sig_m0_cat8, mjj_sig_sigma_cat8, mjj_sig_alpha_cat8, mjj_sig_n_cat8);
MjjSig_cat8      = AddPdf(MjjGaussSig_cat8, MjjCBSig_cat8, mjj_sig_frac_cat8);



mjj_bkg_8TeV_slope1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope3[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope4[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope5[0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg3[0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg4[0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg5[0.1,-100.0, 100.0];

mjj_bkg_8TeV_slope1_cat0[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope2_cat0[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope3_cat0[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope4_cat0[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope5_cat0[0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg1_cat0[-0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg2_cat0[-0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg3_cat0[-0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg4_cat0[-0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg5_cat0[-0.1,-100.0, 100.0];
mjj_bkg_8TeV_wid1_cat0[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid2_cat0[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid3_cat0[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid4_cat0[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid5_cat0[0.1,0.0, 100.0];

mjj_bkg_8TeV_slope1_cat1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope2_cat1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope3_cat1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope4_cat1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope5_cat1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg1_cat1[-0.1,-10.0, 0.0];
mjj_bkg_8TeV_arg2_cat1[-0.1,-10.0, 0.0];
mjj_bkg_8TeV_arg3_cat1[-0.1,-10.0, 0.0];
mjj_bkg_8TeV_arg4_cat1[-0.1,-10.0, 0.0];
mjj_bkg_8TeV_arg5_cat1[-0.1,-10.0, 0.0];
mjj_bkg_8TeV_wid1_cat1[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid2_cat1[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid3_cat1[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid4_cat1[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid5_cat1[0.1,0.0, 100.0];

mjj_bkg_8TeV_slope1_cat2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope2_cat2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope3_cat2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope4_cat2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope5_cat2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_arg1_cat2[-0.1,-10.0, 0.0];
mjj_bkg_8TeV_arg2_cat2[-0.1,-10.0, 0.0];
mjj_bkg_8TeV_arg3_cat2[-0.1,-10.0, 0.0];
mjj_bkg_8TeV_arg4_cat2[-0.1,-10.0, 0.0];
mjj_bkg_8TeV_arg5_cat2[-0.1,-10.0, 0.0];
mjj_bkg_8TeV_wid1_cat2[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid2_cat2[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid3_cat2[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid4_cat2[0.1,0.0, 100.0];
mjj_bkg_8TeV_wid5_cat2[0.1,0.0, 100.0];

wei[1,0,10];
