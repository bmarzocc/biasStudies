mgg[100,180];
mgg_sig_m0[110.0, 70, 160];
mgg_sig_sigma[1.5, 0.5, 2.0];
mgg_sig_alpha[1.0, 0.5, 3]; 
mgg_sig_n[4.0, 0.5, 10]; 
mgg_sig_gsigma[2.5, 2.0, 5];
mgg_sig_frac[0.3, 0, 0.5];

MggGaussSig = Gaussian(mgg, mgg_sig_m0, mgg_sig_gsigma);
MggCBSig    = CBShape(mgg, mgg_sig_m0, mgg_sig_sigma, mgg_sig_alpha, mgg_sig_n);
MggSig      = AddPdf(MggGaussSig, MggCBSig, mgg_sig_frac);

mgg_sig_m0_cat0[110.0, 105, 155];
mgg_sig_sigma_cat0[1.2, 0.8, 2.5];
mgg_sig_alpha_cat0[2.0, 1.0, 2.5]; 
mgg_sig_n_cat0[2.0, 1.0, 5.0]; 
mgg_sig_gsigma_cat0[5.0, 3.0, 8.0];
mgg_sig_frac_cat0[0.1, 0, 0.4];


MggGaussSig_cat0 = Gaussian(mgg, mgg_sig_m0_cat0, mgg_sig_gsigma_cat0);
MggCBSig_cat0    = CBShape(mgg, mgg_sig_m0_cat0, mgg_sig_sigma_cat0, mgg_sig_alpha_cat0, mgg_sig_n_cat0);
MggSig_cat0      = AddPdf(MggGaussSig_cat0, MggCBSig_cat0, mgg_sig_frac_cat0);

mgg_sig_m0_cat1[110.0, 70, 160];
mgg_sig_sigma_cat1[1.8, 1.0, 2.5];
mgg_sig_alpha_cat1[2.0, 1.2, 5]; 
mgg_sig_n_cat1[2.0, 1.5, 10]; 
mgg_sig_gsigma_cat1[5.0, 3.0, 8.0];
mgg_sig_frac_cat1[0.1, 0, 0.4];

MggGaussSig_cat1 = Gaussian(mgg, mgg_sig_m0_cat1, mgg_sig_gsigma_cat1);
MggCBSig_cat1    = CBShape(mgg, mgg_sig_m0_cat1, mgg_sig_sigma_cat1, mgg_sig_alpha_cat1, mgg_sig_n_cat1);
MggSig_cat1      = AddPdf(MggGaussSig_cat1, MggCBSig_cat1, mgg_sig_frac_cat1);

mgg_sig_m0_cat2[110.0, 70, 160];
mgg_sig_sigma_cat2[2.0, 1.0, 3.5];
mgg_sig_alpha_cat2[1.0, 1.0, 3.5]; 
mgg_sig_n_cat2[2.0, 1.5, 10]; 
mgg_sig_gsigma_cat2[5.0, 2.0, 8.0];
mgg_sig_frac_cat2[0.1, 0, 0.3];

MggGaussSig_cat2 = Gaussian(mgg, mgg_sig_m0_cat2, mgg_sig_gsigma_cat2);
MggCBSig_cat2    = CBShape(mgg, mgg_sig_m0_cat2, mgg_sig_sigma_cat2, mgg_sig_alpha_cat2, mgg_sig_n_cat2);
MggSig_cat2      = AddPdf(MggGaussSig_cat2, MggCBSig_cat2, mgg_sig_frac_cat2);

mgg_sig_m0_cat3[110.0, 70, 160];
mgg_sig_sigma_cat3[2.5, 2.0, 3.5];
mgg_sig_alpha_cat3[1.0, 1.0, 3.5]; 
mgg_sig_n_cat3[2.0, 1.5, 10]; 
mgg_sig_gsigma_cat3[5.0, 3.5, 10.0];
mgg_sig_frac_cat3[0.1, 0, 0.3];

MggGaussSig_cat3 = Gaussian(mgg, mgg_sig_m0_cat3, mgg_sig_gsigma_cat3);
MggCBSig_cat3    = CBShape(mgg, mgg_sig_m0_cat3, mgg_sig_sigma_cat3, mgg_sig_alpha_cat3, mgg_sig_n_cat3);
MggSig_cat3      = AddPdf(MggGaussSig_cat3, MggCBSig_cat3, mgg_sig_frac_cat3);

mgg_sig_m0_cat4[110.0, 70, 160];
mgg_sig_sigma_cat4[2.5, 1.2, 3.0];
mgg_sig_alpha_cat4[1.0, 1.0, 3.5]; 
mgg_sig_n_cat4[2.0, 1.5, 10]; 
mgg_sig_gsigma_cat4[5.0, 3.0, 10.0];
mgg_sig_frac_cat4[0.1, 0, 0.3];

MggGaussSig_cat4 = Gaussian(mgg, mgg_sig_m0_cat4, mgg_sig_gsigma_cat4);
MggCBSig_cat4    = CBShape(mgg, mgg_sig_m0_cat4, mgg_sig_sigma_cat4, mgg_sig_alpha_cat4, mgg_sig_n_cat4);
MggSig_cat4      = AddPdf(MggGaussSig_cat4, MggCBSig_cat4, mgg_sig_frac_cat4);


mgg_sig_m0_cat5[110.0, 70, 160];
mgg_sig_sigma_cat5[2.5, 1.2, 3.0];
mgg_sig_alpha_cat5[1.0, 1.0, 3.5]; 
mgg_sig_n_cat5[2.0, 1.5, 10]; 
mgg_sig_gsigma_cat5[5.0, 3.0, 10.0];
mgg_sig_frac_cat5[0.1, 0, 0.3];

MggGaussSig_cat5 = Gaussian(mgg, mgg_sig_m0_cat5, mgg_sig_gsigma_cat5);
MggCBSig_cat5    = CBShape(mgg, mgg_sig_m0_cat5, mgg_sig_sigma_cat5, mgg_sig_alpha_cat5, mgg_sig_n_cat5);
MggSig_cat5      = AddPdf(MggGaussSig_cat5, MggCBSig_cat5, mgg_sig_frac_cat5);


mgg_sig_m0_cat6[110.0, 70, 160];
mgg_sig_sigma_cat6[2.5, 1.2, 3.0];
mgg_sig_alpha_cat6[1.0, 1.0, 3.5]; 
mgg_sig_n_cat6[2.0, 1.5, 10]; 
mgg_sig_gsigma_cat6[5.0, 3.0, 10.0];
mgg_sig_frac_cat6[0.1, 0, 0.3];

MggGaussSig_cat6 = Gaussian(mgg, mgg_sig_m0_cat6, mgg_sig_gsigma_cat6);
MggCBSig_cat6    = CBShape(mgg, mgg_sig_m0_cat6, mgg_sig_sigma_cat6, mgg_sig_alpha_cat6, mgg_sig_n_cat6);
MggSig_cat6      = AddPdf(MggGaussSig_cat6, MggCBSig_cat6, mgg_sig_frac_cat6);


mgg_sig_m0_cat7[110.0, 70, 160];
mgg_sig_sigma_cat7[2.5, 1.2, 3.0];
mgg_sig_alpha_cat7[1.0, 1.0, 3.5]; 
mgg_sig_n_cat7[2.0, 1.5, 10]; 
mgg_sig_gsigma_cat7[5.0, 3.0, 10.0];
mgg_sig_frac_cat7[0.1, 0, 0.3];

MggGaussSig_cat7 = Gaussian(mgg, mgg_sig_m0_cat7, mgg_sig_gsigma_cat7);
MggCBSig_cat7    = CBShape(mgg, mgg_sig_m0_cat7, mgg_sig_sigma_cat7, mgg_sig_alpha_cat7, mgg_sig_n_cat7);
MggSig_cat7      = AddPdf(MggGaussSig_cat7, MggCBSig_cat7, mgg_sig_frac_cat7);


mgg_sig_m0_cat8[110.0, 70, 160];
mgg_sig_sigma_cat8[2.5, 1.2, 3.0];
mgg_sig_alpha_cat8[1.0, 1.0, 3.5]; 
mgg_sig_n_cat8[2.0, 1.5, 10]; 
mgg_sig_gsigma_cat8[5.0, 3.0, 10.0];
mgg_sig_frac_cat8[0.1, 0, 0.3];

MggGaussSig_cat8 = Gaussian(mgg, mgg_sig_m0_cat8, mgg_sig_gsigma_cat8);
MggCBSig_cat8    = CBShape(mgg, mgg_sig_m0_cat8, mgg_sig_sigma_cat8, mgg_sig_alpha_cat8, mgg_sig_n_cat8);
MggSig_cat8      = AddPdf(MggGaussSig_cat8, MggCBSig_cat8, mgg_sig_frac_cat8);



mgg_bkg_8TeV_slope1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope3[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope4[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope5[0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg3[0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg4[0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg5[0.1,-100.0, 100.0];

mgg_bkg_8TeV_slope1_cat0[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope2_cat0[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope3_cat0[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope4_cat0[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope5_cat0[0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg1_cat0[-0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg2_cat0[-0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg3_cat0[-0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg4_cat0[-0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg5_cat0[-0.1,-100.0, 100.0];
mgg_bkg_8TeV_wid1_cat0[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid2_cat0[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid3_cat0[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid4_cat0[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid5_cat0[0.1,0.0, 100.0];

mgg_bkg_8TeV_slope1_cat1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope2_cat1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope3_cat1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope4_cat1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope5_cat1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg1_cat1[-0.1,-10,0.0];
mgg_bkg_8TeV_arg2_cat1[-0.1,-10,0.0];
mgg_bkg_8TeV_arg3_cat1[-0.1,-10,0.0];
mgg_bkg_8TeV_arg4_cat1[-0.1,-10,0.0];
mgg_bkg_8TeV_arg5_cat1[-0.1,-10,0.0];
mgg_bkg_8TeV_wid1_cat1[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid2_cat1[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid3_cat1[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid4_cat1[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid5_cat1[0.1,0.0, 100.0];

mgg_bkg_8TeV_slope1_cat2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope2_cat2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope3_cat2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope4_cat2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope5_cat2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_arg1_cat2[-0.1,-10,0.0];
mgg_bkg_8TeV_arg2_cat2[-0.1,-10,0.0];
mgg_bkg_8TeV_arg3_cat2[-0.1,-10,0.0];
mgg_bkg_8TeV_arg4_cat2[-0.1,-10,0.0];
mgg_bkg_8TeV_arg5_cat2[-0.1,-10,0.0];
mgg_bkg_8TeV_wid1_cat2[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid2_cat2[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid3_cat2[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid4_cat2[0.1,0.0, 100.0];
mgg_bkg_8TeV_wid5_cat2[0.1,0.0, 100.0];

wei[1,0,10];
