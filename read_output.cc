#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TString.h"

#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TCanvas.h"

#include "TPad.h"

//#include "WinPairFilling.h"
#include "CalcRoutineFB.h"

//#include "/opt/mygit/alice_NN_PtPt_2016/utils.C"
#include "/Users/macbookpro/work/alice_NN_PtPt_2016/utils.C" //"/opt/mygit/alice_NN_PtPt_2016/utils.C"


//const int nWinPairWithSubsamples = 20;//15; //10;


// #############
//int read_output_FB_EPOS_May2020() - BASED ON THIS!
int read_output()
{
    TLatex *   tex;

    //    TFile *fIn = new TFile( "output_with_FB.root" ); //"MergedOutput.root" );
    //    TFile *fIn = new TFile( "output_May2020_DIR_0_9.root" ); //"MergedOutput.root" );
    //    TFile *fIn = new TFile( "output_May2020_DIR_0_4.root" ); //"MergedOutput.root" );
    //    TFile *fIn = new TFile( "MergedOutput.root" ); //"MergedOutput.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_50_1Mln_events.root" ); //"MergedOutput.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_50_500k_events_SRC_only_8etaWins.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_25_1mln_events_SRC_only_8etaWins.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_25_1mln_events_SRC_AND_ebye_Kpi_8etaWins.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_25_1mln_events_SRCKpKm_NoEbyeFluctKpi_8etaWins.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_25_800k_events_NoSRC_NoEbyeFluctKpi_AND_meanPt_8etaWins_MORE_VARS.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_1000_1_1M_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPt_8etaWins_EVEN_MORE_VARS.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_1000_1_1M_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPt_8etaWins_EVEN_MORE_VARS_eff08.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_1000_40_1M_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPt_8etaWins_EVEN_MORE_VARS.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_1000_1_1M_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPt_8etaWins_EVEN_MORE_VARS_effPtDep3ranges.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_1000_1_1M_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPt_8etaWins_EVEN_MORE_VARS_effPtDep4ranges_WithCORRECTIONS.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_1000_1_1M_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPt_8etaWins_effPtDep4ranges_GEN_REC_CORRECTED.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_1000_1_30k_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPtBW_TkinFluct_8etaWins_effPtDep4ranges_GEN_REC_CORRECTED.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_1000_1_40k_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPtBW_betaTFluct_8etaWins_effPtDep4ranges_GEN_REC_CORRECTED.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_400_10_1M_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPtExp_8etaWins_effPtDep4ranges_GEN_REC_CORRECTED.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_10_1Mln_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPtExp_8etaWins_effPtDep4ranges_GEN_REC_CORRECTED.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_10_10k_test_withFullEtaInDenom.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_800_10_200k_test_withFullEtaInDenom.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_10_WITH_SIGNAL.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_10.root" );
    //    TFile *fIn = new TFile( "toy_output_gaus_500_10_WITH_SIGNAL_500k.root" );
    //    TFile *fIn = new TFile( "FIST_CE_50k_p_pi_only_plus.root" );
    //    TFile *fIn = new TFile( "FIST_CE_50k_K_pi_only_plus.root" );
    //    TFile *fIn = new TFile( "FIST_GCE_50k_p_pi_only_plus.root" );
    //    TFile *fIn = new TFile( "FIST_GCE_PointParticle_100k_T_145_155_165.root" );
    //    TFile *fIn = new TFile( "toy_output_fraction_K_fixed0.2_nEv1Mln.root" );
    //    TFile *fIn = new TFile( "toy_output_fraction_K_fluctGaus0.7_0.02_nEv1Mln.root" );
    //    TFile *fIn = new TFile( "toy_output_fraction_K_fluctGaus0.8_0.02_nEv1Mln.root" );
    //    TFile *fIn = new TFile( "toy_output_fraction_K_fluctGaus0.8_0.02_SRC_nEv1Mln.root" );

    //    TFile *fIn = new TFile( "toy_output_fraction_K_fixed0.2_SRC_SS_OS_nEv1Mln.root" );
    //        TFile *fIn = new TFile( "toy_output_gaus_500_10.root" );



    // NUCLEUS-2021: For GCE, CE, SF plot, with No SRC:
    //    TString dirName = "output_for_NUCLEUS2021";
    TString dirName = "./";
    //    TString dirName = "/Users/macbookpro/work/FIST2/grid_run_IA_18_10_2021_FIST/out_downloaded/out_from_grid_2021_11_13_first_big_FIST_run_500k_events/";
    //    TString dirName = "/Users/macbookpro/work/FIST2/grid_run_IA_18_10_2021_FIST/out_downloaded/out_from_grid_2021_11_16_first_big_FIST_run_several_Mln_events/";
    //    TString dirName = "/Users/macbookpro/work/FIST2/grid_run_IA_18_10_2021_FIST/out_downloaded/out_from_grid_2021_11_16_big_FIST_run_several_Mln_events_CE_Volume3/";
    //    TString dirName = "/Users/macbookpro/work/FIST2/grid_run_IA_18_10_2021_FIST/out_downloaded/out_from_grid_2021_11_21_big_FIST_run_several_Mln_events_GCE_Volume3/";
    //    TString dirName = "/Users/macbookpro/work/ALICE_analysis/2019_Particle_ratios_NEW_CODE_Pythia_pp_PbPb_analysis/grid_run_IA_25_01_2022_Ang_nuFB_vs_centrality/";
    //    TString dirName = "/Users/macbookpro/work/ALICE_analysis/2019_Particle_ratios_NEW_CODE_Pythia_pp_PbPb_analysis/grid_run_IA_27_01_2022_Ang_nuFB_vs_centrality_OneWideEtaWinPair/";
    //    TString dirName = "/Users/macbookpro/work/ALICE_analysis/2019_Particle_ratios_NEW_CODE_Pythia_pp_PbPb_analysis/grid_run_IA_07_02_2022_Ang_nuFB_vs_centrality_etaBins01/";
    //    TString fileName = "toy_output_fraction_K_fraction_binomial_0.2_NoSRC_nPartGaus80_4_nEv500k.root";
    //    TString fileName = "toy_output_fraction_K_fraction_fixed_0.2_NoSRC_nPartGaus80_4_nEv500k.root";
    //    TString fileName = "toy_output_fraction_K_fraction_binomial_0.2_FLUCT_0.04_NoSRC_nPartGaus80_4_nEv500k.root";

    // NUCLEUS-2021: For GCE, CE, SF plot, WITH SRC:
    //    TString fileName = "toy_output_fraction_K_fraction_binomial_0.2_SRC_OSonly_nPartGaus80_4_nEv500k.root";
    //    TString fileName = "toy_output_fraction_K_fraction_fixed_0.2_SRC_OSonly_nPartGaus80_4_nEv500k.root";
    //    TString fileName = "toy_output_fraction_K_fraction_binomial_0.2_SRC_SS_OS_nPartGaus80_4_nEv500k.root";
    //    TString fileName = "toy_output_fraction_K_fraction_fixed_0.2_SRC_SS_OS_nPartGaus80_4_nEv500k.root";
//    TString fileName = "toy_output_gaus_500_10.root";
    TString fileName = "output_toy.root";
//    TString fileName = "before_phi_output_toy.root";
    //    TString fileName = "toy_output_gaus_500_10_BEFORE_CHANGES.root";
    //    TString fileName = "toy_output_FIST.root";
    //    TString fileName = "output_heavy_ions_test_10k.root";
    //    TString fileName = "MergedOutput_output_heavy_ions_700k.root";
    //    TString fileName = "MergedOutput_output_heavy_ions_4.5mln_events_pt0_100_nSub10.root";
    //    TString fileName = "MergedOutput_output_heavy_ions_8.8mln_events_pt02_20_nSub15.root";
    //    TString fileName = "MergedOutput_10mln_pt_02_20_manyEtaBins.root";


    TFile *fIn = new TFile( dirName+"/"+fileName );





    //    TFile *fIn = new TFile( "toy_output_gaus_500_15_1M_events_NoSRC_WithEbyeFluctKpiGAUS_AND_meanPt_8etaWins_EVEN_MORE_VARS.root" );
    //    TFile *fIn = new TFile( "toy_output_more_npart_fluct.root" ); //"MergedOutput.root" );
    //    TFile *fIn = new TFile( "toy_output_no_collectivity.root" ); //"MergedOutput.root" );
    //    TFile *fIn = new TFile( "toy_output_collectivity_prop_to_Nch.root" ); //"MergedOutput.root" );


    //TFile *fIn = new TFile( "../out_ROPE_HADR/MergedOutput.root" );
    //    TFile *fIn = new TFile( "../out_ROPE_HADR_try2/MergedOutput.root" );

    // #############
    const int nCW = 1;//2;//1;
    const int nCentrBins[] = { 1 };//10, 20 };
    const double cBinWidths[] = { 20 };//10, 5 }; // in %

    const double eRange = 0.8; //1.0;//1.4;//1.6;//0.8;
    //const double eSize = 0.2;//4;
    //const double eStep = 0.2;
    const double eSize = 0.2;//1;//1;//4;
    const double eStep = 0.2;//1;
    //const double eSize = 0.2;//4;
    //const double eStep = 0.2;



    //const int nPhiWins = 16;
    //const int nEtaBins = /*2**/nPhiWins;//23; //(eRange-eSize) / eStep + 2;//1;
//    const int nEtaBins = (eRange-eSize) / eStep + 1 + 0.0000001;//2;//1;
    //const int nEtaBins = 1;
    //const int nEtaBins = /* +1 - for full eta! */ 1 + (eRange-eSize) / eStep + 1 + 0.0000001;//2;//1;


    const int nPtBins = 1;
    //    double ptmin[nPtBins] = { 0 }; //0.2 };//0.1,  };//0.2 }; //0.3, };//0.2 };
    //    double ptmax[nPtBins] = { 100 }; //2.0,};// 2.0 }; //1.5,};// 5.0 };
    double ptmin[nPtBins] = { 0. };
    double ptmax[nPtBins] = { 100.0 };


    int whichInputHist = 1;

    //
    const int nPartTypes = 1;//5;//10;//6;//3;//4;
    int arrPartTypes[nPartTypes][4] =
    { // F,B, X,Y
      //  { /*421*/0, 0, 0, 0 },         // 0
      //  { 0, 0, 0, 0 },         // 1
      //  { 0, 0, 0, 0 },         // 2
      //  { 321, 321, 211, 211 }, // 6

      //      { 2212, 2212, 211, 211 }, // 7
      //      { 2212, 2212, 211, 211 }, // 7
      //      { 2212, 2212, 211, 211 }, // 7
      //      { 2212, 2212, 211, 211 }, // 7
      //      { 2212, 2212, 211, 211 }, // 7


      { 321, 321, 211, 211 }, // 7
//      { 321, 321, 211, 211 }, // 8
//      { 321, 321, 211, 211 }, // 9
//      { 321, 321, 211, 211 }, // 8
//      { 321, 321, 211, 211 }, // 9


      //  { 0, 321, 0, 211 }, // 7
      //  { 0, 321, 0, 211 }, // 8
      //  { 0, 321, 0, 211 }, // 9
    };

    int arrCharges[/*nPartTypes*/][4] =
    { // F,B, X,Y
//      {  0,  0,  0,  0 },
//      { +1, +1, +1, +1 },
//      { -1, -1, -1, -1 },
      { +1, -1, +1, -1 },
//      { -1, +1, -1, +1 },
    };

    const double pSize = TMath::TwoPi()/8;


//    TString strAnLevels[] = { "GEN", "REC", "CORRECTED" };
    TString strAnLevels[] = { "REC" };
    const int N_AN_LEVELS = sizeof(strAnLevels)/sizeof(*strAnLevels);

    TStopwatch timer;
    timer.Start();



    // NEW WRAPPER:
//    TList* fOutputWinPairsLists[nPartTypes];

    CalcWithSubsamples observables[N_AN_LEVELS][nPartTypes][nCW][ nCentrBins[nCW-1] ][nPtBins];
    CalcWithSubsamples observables_FULL_ETA_DENOM[N_AN_LEVELS][nPartTypes][nCW][ nCentrBins[nCW-1] ][nPtBins];
    CalcWithSubsamples observables_FULL_ETA_NUM_AND_DENOM[N_AN_LEVELS][nPartTypes][nCW][ nCentrBins[nCW-1] ][nPtBins];
    //    return 0;


    for ( int iAnLevel = 0; iAnLevel < N_AN_LEVELS; ++iAnLevel )
    {
        cout << ">>>>>> DO " << strAnLevels[iAnLevel] << "..." << endl;
        for ( int iType = 0; iType < nPartTypes; iType++)
            for ( int iCW = 0; iCW < nCW; iCW++)
            {
                //                fIn->cd( Form( "cBinWidth_%.1f", cBinWidths[iCW] ) );
                TString strPID = Form( "pid_%d_%d_%d_%d_charge_%d_%d_%d_%d",
                                       arrPartTypes[iType][0], arrPartTypes[iType][1], arrPartTypes[iType][2], arrPartTypes[iType][3],
                        arrCharges[iType][0], arrCharges[iType][1], arrCharges[iType][2], arrCharges[iType][3] );

//                fOutputWinPairsLists[iType] = new TList();
//                fOutputWinPairsLists[iType] = (TList*)gDirectory->Get("list_PIDcorr_"+strPID );
                TList *inList = (TList*)gDirectory->Get("list_PIDcorr_"+strPID );
                cout << "taking list " << strPID << ": " << inList << endl;
//                fOutputWinPairsLists[iType]->ls();

                bool ifIdent = ifIdenticalParticlesFB( arrPartTypes[iType], arrCharges[iType] );
                cout << "ifIdentical = " << ifIdent << endl;

                for ( int cBin = 0; cBin < nCentrBins[iCW]; ++cBin )
                    for ( int iPt = 0; iPt < nPtBins; ++iPt )
                    {
                        CalcWithSubsamples *calcObj = &observables[iAnLevel][iType][iCW][ cBin ][iPt];
                        TString strPostfix = Form("%s_cBin%d", strAnLevels[iAnLevel].Data(), cBin);
                        calcObj->calc( inList, strPostfix, ifIdent, whichInputHist ); // eSize, eSize, ifIdent, pSize );

                        // FULL ETA DENOM:
                        calcObj = &observables_FULL_ETA_DENOM[iAnLevel][iType][iCW][ cBin ][iPt];
                        strPostfix = Form("%s_FULL_ETA_DENOM_cBin%d", strAnLevels[iAnLevel].Data(), cBin);
                        calcObj->calc( inList, strPostfix, ifIdent, whichInputHist );


                        // FULL ETA FOR both Num and Denom: to calc nu_dyn!
                        calcObj = &observables_FULL_ETA_NUM_AND_DENOM[iAnLevel][iType][iCW][ cBin ][iPt];
                        strPostfix = Form("%s_FULL_ETA_NUM_AND_DENOM_cBin%d", strAnLevels[iAnLevel].Data(), cBin);
                        calcObj->calc( inList, strPostfix, ifIdent, 0 );
                    }
            }
    }
    cout << "##### End of reading win info from histograms." << endl;

//    const int nEtaBins = observables[0][0][0][ 0 ][0].h3D->GetYaxis()->GetNbins();
//    cout << "nEtaBins = " << nEtaBins << endl;
//    return 0;


    // ### End of file reading. Time:
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();

    printf("End of calc: RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

    //    return 0;



    // ##### DRAWING:
    TStopwatch timer2;
    timer2.Start();

    //    int pTypeForAn = 1;
//    int ptBinForAn = 0;//1;//0;//1;


    // QA: nF vs centrality:
    TCanvas *canv_nF = new TCanvas("canv_nF","canv_nF",0,0,800,600);
    tuneCanvas( canv_nF );
//    drawGraph( gr_nF[ 0 ][0][0][0], 20, kRed, "APz" );










    // #################################
    // #### Several basic observables
    TCanvas *canv_basic_1D_obs_VS_ETA = new TCanvas("canv_basic_1D_obs_VS_ETA","canv_basic_obs_VS_ETA",220,55,800,600);
    tuneCanvas( canv_basic_1D_obs_VS_ETA );
    gPad->SetGridy();


    //    TGraphErrors *gr3;
    //    gr0 = gr_SIGMA_VS_ETA[ 0 ][0][0];
    //    gr1 = gr_SIGMA_VS_ETA[ 1 ][0][0];
    //    gr2 = gr_SIGMA_VS_ETA[ 2 ][0][0];
    //    gr0 = gr_SIGMA_VS_ETA   [ 0 ][0][0];
    //    gr1 = gr_SIGMA_XY_VS_ETA[ 0 ][0][0];
    //    gr2 = gr_SIGMA_FY_VS_ETA[ 0 ][0][0];
    //    gr3 = gr_SIGMA_XB_VS_ETA[ 0 ][0][0];

    TGraphErrors *gr_F = observables[0][0][0][0][0].getGraphVsWinIdQA( "avF" );
    TGraphErrors *gr_B = observables[0][0][0][0][0].getGraphVsWinIdQA( "avB" );
    TGraphErrors *gr_X = observables[0][0][0][0][0].getGraphVsWinIdQA( "avX" );
    TGraphErrors *gr_Y = observables[0][0][0][0][0].getGraphVsWinIdQA( "avY" );


    gr_F->GetYaxis()->SetRangeUser(-1, 2);

    drawGraph( gr_F, 20, kRed, "APz" );
    drawGraph( gr_B, 24, kBlue, "Pz" );
    drawGraph( gr_X, 25, kMagenta, "Pz" );
    drawGraph( gr_Y, 5, kGray+2, "Pz" );



    TLegend *leg_1D_obs_vs_eta = new TLegend(0.64,0.55,0.94,0.82);
    tuneLegend( leg_1D_obs_vs_eta );
    leg_1D_obs_vs_eta->AddEntry( gr_F, "F", "p");
    leg_1D_obs_vs_eta->AddEntry( gr_B, "B", "p");
    leg_1D_obs_vs_eta->AddEntry( gr_X, "X", "p");
    leg_1D_obs_vs_eta->AddEntry( gr_Y, "Y", "p");
    leg_1D_obs_vs_eta->Draw();


    // ###
    TCanvas *canv_basic_2D_obs_VS_ETA = new TCanvas("canv_basic_2D_obs_VS_ETA","canv_basic_2D_obs_VS_ETA",220,55,800,600);
    tuneCanvas( canv_basic_2D_obs_VS_ETA );
    gPad->SetGridy();

    TGraphErrors *gr_FB = observables[0][0][0][0][0].getGraphVsWinIdQA( "FB" );
    TGraphErrors *gr_XY = observables[0][0][0][0][0].getGraphVsWinIdQA( "XY" );
    TGraphErrors *gr_FY = observables[0][0][0][0][0].getGraphVsWinIdQA( "FY" );
    TGraphErrors *gr_XB = observables[0][0][0][0][0].getGraphVsWinIdQA( "XB" );


    drawGraph( gr_FB, 20, kRed, "APz" );
    drawGraph( gr_XY, 24, kBlue, "Pz" );
    drawGraph( gr_FY, 25, kMagenta, "Pz" );
    drawGraph( gr_XB, 5, kGray+2, "Pz" );









    // #######################################
    // ### Sigma correlations VS winId
    TCanvas *canv_sigma_VS_winId = new TCanvas("canv_sigma_VS_winId","canv_sigma_VS_winId",220,55,800,600);
    tuneCanvas( canv_sigma_VS_winId );

    TGraphErrors *gr_sigma_FB_vs_winId = observables[0][0][0][0][0].getGraphVsWinIdQA( "sigma_FB" );
    TGraphErrors *gr_sigma_XY_vs_winId = observables[0][0][0][0][0].getGraphVsWinIdQA( "sigma_XY" );
    TGraphErrors *gr_sigma_FY_vs_winId = observables[0][0][0][0][0].getGraphVsWinIdQA( "sigma_FY" );
    TGraphErrors *gr_sigma_XB_vs_winId = observables[0][0][0][0][0].getGraphVsWinIdQA( "sigma_XB" );




    tuneGraphAxisLabels( gr_sigma_FB_vs_winId );
    gr_sigma_FB_vs_winId->GetYaxis()->SetRangeUser(0.88, 1.2);
    gr_sigma_FB_vs_winId->SetTitle( ";winId;#Sigma_{FB}" );

    drawGraph( gr_sigma_FB_vs_winId, 20, kRed, "APz" );
    drawGraph( gr_sigma_XY_vs_winId, 24, kBlue, "Pz" );
    drawGraph( gr_sigma_FY_vs_winId, 25, kMagenta, "Pz" );
    drawGraph( gr_sigma_XB_vs_winId, 5, kGray+2, "Pz" );


    TLegend *leg_SIGMAFB_vs_winId = new TLegend(0.64,0.55,0.94,0.82);
    tuneLegend( leg_SIGMAFB_vs_winId );
    //    leg_SIGMAFB_vs_eta->AddEntry( gr0, "K/#pi, K/#pi, formula", "p");
    //    leg_SIGMAFB_vs_eta->AddEntry( gr1, "K+/#pi+, K+/#pi+, formula", "p");
    //    leg_SIGMAFB_vs_eta->AddEntry( gr2, "K+/#pi+, K#minus/#pi#minus, formula", "p");
    leg_SIGMAFB_vs_winId->AddEntry( gr_sigma_FB_vs_winId, "FB", "p");
    leg_SIGMAFB_vs_winId->AddEntry( gr_sigma_XY_vs_winId, "XY", "p");
    leg_SIGMAFB_vs_winId->AddEntry( gr_sigma_FY_vs_winId, "FY", "p");
    leg_SIGMAFB_vs_winId->AddEntry( gr_sigma_XB_vs_winId, "XB", "p");
    leg_SIGMAFB_vs_winId->Draw();


    gPad->SetGrid();


    // #######################################
    // ### Sigma correlations VS dEta:
    TCanvas *canv_sigma_VS_ETA = new TCanvas("canv_sigma_VS_ETA","canv_sigma_vs_eta",220,55,800,600);
    tuneCanvas( canv_sigma_VS_ETA );


    //    TGraphErrors *gr3;
    //    gr0 = gr_SIGMA_VS_ETA[ 0 ][0][0];
    //    gr1 = gr_SIGMA_VS_ETA[ 1 ][0][0];
    //    gr2 = gr_SIGMA_VS_ETA[ 2 ][0][0];
    //    gr0 = gr_SIGMA_VS_ETA   [ 0 ][0][0];
    //    gr1 = gr_SIGMA_XY_VS_ETA[ 0 ][0][0];
    //    gr2 = gr_SIGMA_FY_VS_ETA[ 0 ][0][0];
    //    gr3 = gr_SIGMA_XB_VS_ETA[ 0 ][0][0];

    TGraphErrors *gr_sigma_FB = observables[0][0][0][0][0].getGraph( "sigma_FB" );
    TGraphErrors *gr_sigma_XY = observables[0][0][0][0][0].getGraph( "sigma_XY" );
    TGraphErrors *gr_sigma_FY = observables[0][0][0][0][0].getGraph( "sigma_FY" );
    TGraphErrors *gr_sigma_XB = observables[0][0][0][0][0].getGraph( "sigma_XB" );




    tuneGraphAxisLabels( gr_sigma_FB );
    gr_sigma_FB->GetYaxis()->SetRangeUser(0.88, 1.2);
    gr_sigma_FB->SetTitle( ";#Delta#eta;#Sigma_{FB}" );

    drawGraph( gr_sigma_FB, 20, kRed, "APz" );
    drawGraph( gr_sigma_XY, 24, kBlue, "Pz" );
    drawGraph( gr_sigma_FY, 25, kMagenta, "Pz" );
    drawGraph( gr_sigma_XB, 5, kGray+2, "Pz" );


    TLegend *leg_SIGMAFB_vs_eta = new TLegend(0.64,0.55,0.94,0.82);
    tuneLegend( leg_SIGMAFB_vs_eta );
    //    leg_SIGMAFB_vs_eta->AddEntry( gr0, "K/#pi, K/#pi, formula", "p");
    //    leg_SIGMAFB_vs_eta->AddEntry( gr1, "K+/#pi+, K+/#pi+, formula", "p");
    //    leg_SIGMAFB_vs_eta->AddEntry( gr2, "K+/#pi+, K#minus/#pi#minus, formula", "p");
    leg_SIGMAFB_vs_eta->AddEntry( gr_sigma_FB, "FB", "p");
    leg_SIGMAFB_vs_eta->AddEntry( gr_sigma_XY, "XY", "p");
    leg_SIGMAFB_vs_eta->AddEntry( gr_sigma_FY, "FY", "p");
    leg_SIGMAFB_vs_eta->AddEntry( gr_sigma_XB, "XB", "p");
    leg_SIGMAFB_vs_eta->Draw();


    gPad->SetGrid();






    // ###### Draw dEta-dPhi 2D
    for ( int iType = 0; iType < 1/*nPartTypes*/; iType++)
    {
        TString strPID = Form( "pid_%d_%d_%d_%d_charge_%d_%d_%d_%d",
                               arrPartTypes[iType][0], arrPartTypes[iType][1], arrPartTypes[iType][2], arrPartTypes[iType][3],
                arrCharges[iType][0], arrCharges[iType][1], arrCharges[iType][2], arrCharges[iType][3] );

        TString strCanvName = Form( "canv_Graph2D_type_%s", strPID.Data() );
        TCanvas *canv_ratio_ratio_VS_ETA_Graph2D = new TCanvas( strCanvName, strCanvName, 10+30*iType,10 +30*iType,800,800);
        canv_ratio_ratio_VS_ETA_Graph2D->Divide(2,2);

        canv_ratio_ratio_VS_ETA_Graph2D->cd(1);
//        TGraph2D *gr2D_0_rr_GEN = observables[0][iType][0][0][0].getGraph2D( "corr_rr_formula" );
        TH2D *gr2D_0_sigmaXY_GEN = observables[0][iType][0][0][0].getHist2D( "sigma_XY" );
        gr2D_0_sigmaXY_GEN->Draw("surf1");


        canv_ratio_ratio_VS_ETA_Graph2D->cd(2);
        TH2D *gr2D_0_rr_GEN = observables[0][iType][0][0][0].getHist2D( "corr_rr_formula" );
        gr2D_0_rr_GEN->Draw("surf1");

        canv_ratio_ratio_VS_ETA_Graph2D->cd(3);
        TH2D *gr2D_0_rPt_GEN = observables[0][iType][0][0][0].getHist2D( "corr_rPt_formula" );
        gr2D_0_rPt_GEN->Draw("surf1");

        canv_ratio_ratio_VS_ETA_Graph2D->cd(4);
        TH2D *gr2D_0_PtPt_GEN = observables[0][iType][0][0][0].getHist2D( "avPtF_avPtB_formula" );
        gr2D_0_PtPt_GEN->Draw("surf1");

    }






return 0;












    // #####################################################
    // ### Ratio - Ratio correlations vs centrality:
    // #####################################################
    if (0)
    {
//        TCanvas *canv_ratio_ratio = new TCanvas("canv_ratio_ratio","canv_ratio_ratio",120,55,800,600);
//        tuneCanvas( canv_ratio_ratio );
//        gPad->SetGridy();

//        ////    TGraphErrors *grRatioPt = gr_ptpt_bcorr[ 1 ][0][1];
//        tuneGraphAxisLabels( gr_corr_rr_formula[ 0 ][0][0][0] );
//        gr_corr_rr_formula[ 0 ][0][0][0]->GetYaxis()->SetRangeUser(-0.03, 0.03 );
//        drawGraph( gr_corr_rr_formula[ 0 ][0][0][0], 20, kRed, "APz" );
//        drawGraph( gr_corr_rr_formula[ 0 ][1][0][0], 24, kBlue, "Pz" );
//        drawGraph( gr_corr_rr_formula[ 0 ][3][0][0], 24, kMagenta, "Pz" );

//        drawGraph( gr_corr_rr_direct[ 0 ][0][0][0], 20, kRed, "Lz" );
//        drawGraph( gr_corr_rr_direct[ 0 ][1][0][0], 24, kBlue, "Lz" );
//        drawGraph( gr_corr_rr_direct[ 0 ][3][0][0], 24, kMagenta, "Lz" );
    }


    // write to file (VS CENTR)
    TString strGrNames[] = { "gr_all", "graph_PP", "graph_MM", "graph_PM", "graph_MP" };
    if(0)
    {


        TFile *fileWithGraphsVsCentr = new TFile( Form( "graphs_oneWideEtaWinPair_VS_CENTR_cW%.0f_p_pi_from_%s", cBinWidths[0], fileName.Data() ), "RECREATE" );
        fileWithGraphsVsCentr->mkdir("formula")->cd();
        for ( int iType = 0; iType < nPartTypes; iType++)
        {
//            gr_corr_rr_formula[ 0 ][iType][0][0]->SetName( strGrNames[iType] );
//            gr_corr_rr_formula[ 0 ][iType][0][0]->Write();
            observables[0][iType][0][0][0].getGraph( "corr_rr_formula" )->SetName( strGrNames[iType] );  // [iAnLevel][iType][iCW][cBin][ptBinForAn]
            observables[0][iType][0][0][0].getGraph( "corr_rr_formula" )->Write();
        }
        fileWithGraphsVsCentr->mkdir("direct")->cd();
        for ( int iType = 0; iType < nPartTypes; iType++)
        {
//            gr_corr_rr_direct[ 0 ][iType][0][0]->SetName( strGrNames[iType] );
//            gr_corr_rr_direct[ 0 ][iType][0][0]->Write();
            observables[0][iType][0][0][0].getGraph( "corr_rr_direct" )->SetName( strGrNames[iType] );
            observables[0][iType][0][0][0].getGraph( "corr_rr_direct" )->Write();
        }
        fileWithGraphsVsCentr->Close();

    }






    // ##### VS ETA:
    TCanvas *canv_ratio_ratio_VS_ETA = new TCanvas("canv_ratio_ratio_VS_ETA","canv_ratio_ratio_VS_ETA",140,75,800,600);
    tuneCanvas( canv_ratio_ratio_VS_ETA );

    //    TGraphErrors *grOnlySim = gr_corr_rr_formula_VS_ETA[ 0 ][0][0][0];
    //    multiplyPointsByFactor( grOnlySim, 1./0.2, true );
    //    grOnlySim->GetYaxis()->SetRangeUser(0,2);
    //    grOnlySim->SetTitle( ";#Delta#eta;#nu_{FB} #times #LTdN_{#pi^{+}}/d#eta#GT" );
    //    tuneGraphAxisLabels( grOnlySim );
    //    drawGraph( grOnlySim, 20, kRed, "APz" );

    //    return 0;




    // formula:
//    TGraphErrors *gr0 = gr_corr_rr_formula_VS_ETA[ 0 ][0][0][0];
//    //    TGraphErrors *gr1 = gr_corr_rr_formula_VS_ETA[ 1 ][0][0][0];
//    TGraphErrors *gr1 = gr_corr_rr_formula_VS_ETA[ 0 ][1][0][0];
//    TGraphErrors *gr2 = gr_corr_rr_formula_VS_ETA[ /*2*/0 ][2][0][0];
//    TGraphErrors *gr3 = gr_corr_rr_formula_VS_ETA[ /*2*/0 ][3][0][0];
//    TGraphErrors *gr4 = gr_corr_rr_formula_VS_ETA[ /*2*/0 ][4][0][0];


    TGraphErrors *gr0_corr_rr_formula = observables[0][0][0][0][0].getGraph( "corr_rr_formula" ); // [iAnLevel][iType][iCW][cBin][ptBinForAn]
    TGraphErrors *gr1_corr_rr_formula = observables[0][1][0][0][0].getGraph( "corr_rr_formula" );
    TGraphErrors *gr2_corr_rr_formula = observables[0][2][0][0][0].getGraph( "corr_rr_formula" );
    TGraphErrors *gr3_corr_rr_formula = observables[0][3][0][0][0].getGraph( "corr_rr_formula" );
    TGraphErrors *gr4_corr_rr_formula = observables[0][4][0][0][0].getGraph( "corr_rr_formula" );


    double factor_n_to_N = (2*eRange)/eSize;
    cout << ">>>>>> Integral rr = " << observables[0][0][0][0][0].getIntegral( "corr_rr_formula" ) /factor_n_to_N /factor_n_to_N << endl;
    cout << ">>>>>> observables_FULL_ETA_NUM_AND_DENOM: Integral rr = " << observables_FULL_ETA_NUM_AND_DENOM[0][0][0][0][0].getIntegral( "corr_rr_formula" ) << endl;



    //    TGraphErrors *gr0 = gr_corr_rr_formula_FullEtaInDenom_VS_ETA[ 0 ][0][0][0];
    //    TGraphErrors *gr1 = gr_corr_rr_formula_FullEtaInDenom_VS_ETA[ 1 ][0][0][0];
    //    TGraphErrors *gr2 = gr_corr_rr_formula_FullEtaInDenom_VS_ETA[ 2 ][0][0][0];

    // if scale factor is *<dN/deta>
    //    multiplyPointsByFactor( gr0, 1./2, true ); // from all pi to pi+
    //    multiplyPointsByFactor( gr0, 1./0.1, true );
    //    multiplyPointsByFactor( gr1, 1./0.1, true );
    //    multiplyPointsByFactor( gr2, 1./0.1, true );

    // if scale factor is <na><nb>/(<na> + <nb>), we don't need to scale


    //    multiplyPointsByFactor( gr0_FullEtaInDenom, 1./0.2, true );

    gr0_corr_rr_formula->GetYaxis()->SetRangeUser(0,2);
    gr0_corr_rr_formula->SetTitle( ";#Delta#eta;#nu_{FB} #times #LTdN_{#pi^{+}}/d#eta#GT" );
    tuneGraphAxisLabels( gr0_corr_rr_formula );

    drawGraph( gr0_corr_rr_formula, 27, kGreen+1,   "APz" );
    drawGraph( gr1_corr_rr_formula, 20, kRed,       "Pz" );
    drawGraph( gr2_corr_rr_formula, 24, kRed+2,     "Pz" );
    drawGraph( gr3_corr_rr_formula, 21, kBlue,      "Pz" );
    drawGraph( gr4_corr_rr_formula, 25, kBlue+2,    "Pz" );


    TLegend *leg_nuFB_vs_eta = new TLegend(0.34,0.55,0.94,0.82);
    tuneLegend( leg_nuFB_vs_eta );
    leg_nuFB_vs_eta->SetNColumns(2);
    leg_nuFB_vs_eta->AddEntry( gr0_corr_rr_formula, "K/#pi, K/#pi, formula", "p");
    leg_nuFB_vs_eta->AddEntry( gr1_corr_rr_formula, "K+/#pi+, K+/#pi+, formula", "p");
    leg_nuFB_vs_eta->AddEntry( gr2_corr_rr_formula, "K#minus/#pi#minus, K#minus/#pi#minus, formula", "p");
    leg_nuFB_vs_eta->AddEntry( gr3_corr_rr_formula, "K+/#pi+, K#minus/#pi#minus, formula", "p");
    leg_nuFB_vs_eta->AddEntry( gr4_corr_rr_formula, "K#minus/#pi#minus, K+/#pi+, formula", "p");
    //    leg_nuFB_vs_eta->AddEntry( gr0, "gen K/#pi, K/#pi, formula", "p");
    //    leg_nuFB_vs_eta->AddEntry( gr1, "rec K/#pi, K/#pi, formula", "p");
    //    leg_nuFB_vs_eta->AddEntry( gr2, "corrected K/#pi, K/#pi, formula", "p");
    //    leg_nuFB_vs_eta->AddEntry( gr0_FullEtaInDenom, "gen FULL ETA K/#pi, K/#pi, formula", "p");


    // ##### SAVING TO FILE:
    //    {  0,  0,  0,  0 },
    //    { +1, +1, +1, +1 },
    //    { -1, -1, -1, -1 },

    //    { +1, -1, +1, -1 },
    //    { -1, +1, -1, +1 },


    //    gr1->SetName("graph_PP");
    //    gr2->SetName("graph_MM");
    //    gr3->SetName("graph_PM");
    //    gr4->SetName("graph_MP");
    //    //    gr0->SaveAs("graph_FIST_500k_GCE_DiagonalEV_500k_kaons_plus_plus.root");
    //    //    cout << fIn->GetName() << endl;
    //    //    TFile *fileWithGraphs = new TFile( "graphs_from.root", "RECREATE" );
    //    gr0->Write();
    //    gr1->Write();
    //    gr2->Write();
    //    gr3->Write();
    //    gr4->Write();

    // write to file (VS ETA)
    if(0)
    {
        //        TFile *fileWithGraphsVsEta = new TFile( Form( "graphs_oneWideEtaWinPair_VS_ETA_cW%.0f_p_pi_from_%s", cBinWidths[0], fileName.Data() ), "RECREATE" );
        TFile *fileWithGraphsVsEta = new TFile( Form( "graphs_VS_ETA_cW%.0f_p_pi_from_%s", cBinWidths[0], fileName.Data() ), "RECREATE" );
        for ( int iBin = 0; iBin < nCentrBins[0]; iBin++)
        {
            //            fileWithGraphsVsEta->mkdir("formula")->cd();
            fileWithGraphsVsEta->mkdir( Form("formula_cBin%d", iBin) )->cd();
            for ( int iType = 0; iType < nPartTypes; iType++)
            {
//                gr_corr_rr_formula_VS_ETA[ 0 ][iType][0][iBin]->SetName( strGrNames[iType] );
//                gr_corr_rr_formula_VS_ETA[ 0 ][iType][0][iBin]->Write();
                observables[0][iType][0][iBin][0].getGraph( "corr_rr_formula" )->SetName( strGrNames[iType] ); // [iAnLevel][iType][iCW][cBin][ptBinForAn]
                observables[0][iType][0][iBin][0].getGraph( "corr_rr_formula" )->Write();
            }

        }
        for ( int iBin = 0; iBin < nCentrBins[0]; iBin++)
        {
            //        fileWithGraphsVsEta->mkdir("direct")->cd();
            fileWithGraphsVsEta->mkdir( Form("direct_cBin%d", iBin) )->cd();
            for ( int iType = 0; iType < nPartTypes; iType++)
            {
//                gr_corr_rr_direct_VS_ETA[ 0 ][iType][0][iBin]->SetName( strGrNames[iType] );
//                gr_corr_rr_direct_VS_ETA[ 0 ][iType][0][iBin]->Write();
                observables[0][iType][0][iBin][0].getGraph( "corr_rr_direct" )->SetName( strGrNames[iType] ); // [iAnLevel][iType][iCW][cBin][ptBinForAn]
                observables[0][iType][0][iBin][0].getGraph( "corr_rr_direct" )->Write();
            }
        }
        fileWithGraphsVsEta->Close();
    }

    // #### write graphs vs ETA to file
    //    if(0)
    //    {
    //        TFile *file_output_graphs = new TFile( Form( "graphs_Toy_cW%.0f.root", cBinWidths[0] ), "RECREATE" );
    //        for ( int iType = 0; iType < nPartTypes; iType++)
    //        {
    //            TString strPID = Form( "pid_%d_%d_%d_%d_charge_%d_%d_%d_%d",
    //                                   arrPartTypes[iType][0], arrPartTypes[iType][1], arrPartTypes[iType][2], arrPartTypes[iType][3],
    //                    arrCharges[iType][0], arrCharges[iType][1], arrCharges[iType][2], arrCharges[iType][3] );

    //            for ( int cBin = 0; cBin < nCentrBins[0]; ++cBin )
    //            {
    //                file_output_graphs->WriteObject( gr_corr_rr_direct_VS_ETA[ /*N_AN_LEVELS*/0 ][iType][/*nCW*/0][cBin]
    //                        , Form( "gr_corr_rr_direct_VS_ETA_cW%.0f_cBin%d_%s", cBinWidths[0], cBin, strPID.Data() ) );

    //                file_output_graphs->WriteObject( gr_corr_rr_formula_VS_ETA[ /*N_AN_LEVELS*/0 ][iType][/*nCW*/0][cBin]
    //                        , Form( "gr_corr_rr_formula_VS_ETA_cW%.0f_cBin%d_%s", cBinWidths[0], cBin, strPID.Data() ) );

    //            }
    //        }
    //    }



    gr0_corr_rr_formula->Fit("pol0");
    gr0_corr_rr_formula->GetFunction("pol0")->Draw("same");


    // direct:
    //    gr0 = gr_corr_rr_direct_VS_ETA[ 0 ][0][0];
    //    gr1 = gr_corr_rr_direct_VS_ETA[ 1 ][0][0];
    //    gr2 = gr_corr_rr_direct_VS_ETA[ 2 ][0][0];
//    gr0 = gr_corr_rr_direct_VS_ETA[ 0 ][0][0][0];
//    gr1 = gr_corr_rr_direct_VS_ETA[ 0 ][1][0][0];
//    gr2 = gr_corr_rr_direct_VS_ETA[ 0 ][2][0][0];

    TGraphErrors *gr0_corr_rr_direct = observables[0][0][0][0][0].getGraph( "corr_rr_direct" );
    TGraphErrors *gr1_corr_rr_direct = observables[0][1][0][0][0].getGraph( "corr_rr_direct" );
    TGraphErrors *gr2_corr_rr_direct = observables[0][2][0][0][0].getGraph( "corr_rr_direct" );
    TGraphErrors *gr3_corr_rr_direct = observables[0][3][0][0][0].getGraph( "corr_rr_direct" );
    TGraphErrors *gr4_corr_rr_direct = observables[0][4][0][0][0].getGraph( "corr_rr_direct" );

    //    gr1 = gr_corr_rr_direct_VS_ETA[ 1 ][0][0][0];
    //    gr2 = gr_corr_rr_direct_VS_ETA[ 2 ][0][0][0];
    //     gr0 = gr_corr_rr_formula_FullEtaInDenom_VS_ETA[ 0 ][2][0][0];
    //     gr1 = gr_corr_rr_formula_FullEtaInDenom_VS_ETA[ 1 ][2][0][0];
    //     gr2 = gr_corr_rr_formula_FullEtaInDenom_VS_ETA[ 2 ][2][0][0];

    // if scale factor is *<dN/deta>
    //    multiplyPointsByFactor( gr0, 1./2, true ); // from all pi to pi+
    //    multiplyPointsByFactor( gr0, 1./0.1, true );
    //    multiplyPointsByFactor( gr1, 1./0.1, true );
    //    multiplyPointsByFactor( gr2, 1./0.1, true );


    gr0_corr_rr_direct->SetFillColor( kGreen+1 );
    gr1_corr_rr_direct->SetFillColor( kRed );
    gr2_corr_rr_direct->SetFillColor( kRed+2 );
    gr3_corr_rr_direct->SetFillColor( kBlue );
    gr4_corr_rr_direct->SetFillColor( kBlue+2 );

    gr0_corr_rr_direct->SetFillStyle(3005);
    gr1_corr_rr_direct->SetFillStyle(3005);
    gr2_corr_rr_direct->SetFillStyle(3005);
    gr3_corr_rr_direct->SetFillStyle(3005);
    gr4_corr_rr_direct->SetFillStyle(3005);


    drawGraph( gr0_corr_rr_direct, 1, kGreen+1, "le3" );
    drawGraph( gr1_corr_rr_direct, 1, kRed, "le3" );
    drawGraph( gr2_corr_rr_direct, 1, kRed+2, "le3" );
    drawGraph( gr3_corr_rr_direct, 1, kBlue, "le3" );
    drawGraph( gr4_corr_rr_direct, 1, kBlue+2, "le3" );

    leg_nuFB_vs_eta->AddEntry( gr0_corr_rr_direct, "K/#pi, K/#pi, direct", "lf");
    leg_nuFB_vs_eta->AddEntry( gr1_corr_rr_direct, "K+/#pi+, K+/#pi+, direct", "lf");
    leg_nuFB_vs_eta->AddEntry( gr2_corr_rr_direct, "K+/#pi+, K#minus/#pi#minus, direct", "lf");
    leg_nuFB_vs_eta->AddEntry( gr3_corr_rr_direct, "K+/#pi+, K#minus/#pi#minus, direct", "lf");
    leg_nuFB_vs_eta->AddEntry( gr4_corr_rr_direct, "K#minus/#pi#minus, K+/#pi+, direct", "lf");





    //    leg_nuFB_vs_eta->AddEntry( gr0, "gen K/#pi, K/#pi, direct", "lf");
    //    leg_nuFB_vs_eta->AddEntry( gr0, "gen K/#pi, K/#pi, direct", "lf");
    //    leg_nuFB_vs_eta->AddEntry( gr1, "gen K/#pi, K/#pi, direct", "lf");
    //    leg_nuFB_vs_eta->AddEntry( gr1, "rec K/#pi, K/#pi, direct", "lf");
    //    leg_nuFB_vs_eta->AddEntry( gr2, "corrected K/#pi, K/#pi, direct", "lf");


    //    return 0;






    // r-r RECONSTRUCTED:
    TGraphErrors *gr0_rr_REC = observables[1][0][0][0][0].getGraph( "corr_rr_formula" );  // gr_corr_rr_formula_VS_ETA[ 1 ][0][0][0];
    TGraphErrors *gr1_rr_REC = observables[1][1][0][0][0].getGraph( "corr_rr_formula" );  // gr_corr_rr_formula_VS_ETA[ 1 ][1][0][0];
    TGraphErrors *gr2_rr_REC = observables[1][2][0][0][0].getGraph( "corr_rr_formula" );  // gr_corr_rr_formula_VS_ETA[ 1 ][2][0][0];
    TGraphErrors *gr3_rr_REC = observables[1][3][0][0][0].getGraph( "corr_rr_formula" );  // gr_corr_rr_formula_VS_ETA[ 1 ][3][0][0];
    TGraphErrors *gr4_rr_REC = observables[1][4][0][0][0].getGraph( "corr_rr_formula" );  // gr_corr_rr_formula_VS_ETA[ 1 ][4][0][0];

    gr0_rr_REC->GetYaxis()->SetRangeUser(-2, 2);

    shiftPointX( gr0_rr_REC, 0.04 );
    shiftPointX( gr1_rr_REC, 0.04 );
    shiftPointX( gr2_rr_REC, 0.04 );
    shiftPointX( gr3_rr_REC, 0.04 );
    shiftPointX( gr4_rr_REC, 0.04 );

    drawGraph( gr0_rr_REC, 24, kGreen+1,   "Pz", 1.5 );
    drawGraph( gr1_rr_REC, 24, kRed,       "Pz", 1.5 );
    drawGraph( gr2_rr_REC, 24, kRed+1,     "Pz", 1.5 );
    drawGraph( gr3_rr_REC, 24, kBlue,      "Pz", 1.5 );
    drawGraph( gr4_rr_REC, 24, kBlue+1,    "Pz", 1.5 );


    // r-r CORRECTED:
    TGraphErrors *gr0_rr_CORR = observables[2][0][0][0][0].getGraph( "corr_rr_formula" );  //gr_corr_rr_formula_VS_ETA[ 2 ][0][0][0];
    TGraphErrors *gr1_rr_CORR = observables[2][1][0][0][0].getGraph( "corr_rr_formula" );  //gr_corr_rr_formula_VS_ETA[ 2 ][1][0][0];
    TGraphErrors *gr2_rr_CORR = observables[2][2][0][0][0].getGraph( "corr_rr_formula" );  //gr_corr_rr_formula_VS_ETA[ 2 ][2][0][0];
    TGraphErrors *gr3_rr_CORR = observables[2][3][0][0][0].getGraph( "corr_rr_formula" );  //gr_corr_rr_formula_VS_ETA[ 2 ][3][0][0];
    TGraphErrors *gr4_rr_CORR = observables[2][4][0][0][0].getGraph( "corr_rr_formula" );  //gr_corr_rr_formula_VS_ETA[ 2 ][4][0][0];

    shiftPointX( gr0_rr_CORR, 0.04 );
    shiftPointX( gr1_rr_CORR, 0.04 );
    shiftPointX( gr2_rr_CORR, 0.04 );
    shiftPointX( gr3_rr_CORR, 0.04 );
    shiftPointX( gr4_rr_CORR, 0.04 );

    drawGraph( gr0_rr_CORR, 25, kGreen+1,   "Pz", 1.5 );
    drawGraph( gr1_rr_CORR, 25, kRed,       "Pz", 1.5 );
    drawGraph( gr2_rr_CORR, 25, kRed+1,     "Pz", 1.5 );
    drawGraph( gr3_rr_CORR, 25, kBlue,      "Pz", 1.5 );
    drawGraph( gr4_rr_CORR, 25, kBlue+1,    "Pz", 1.5 );




    // #### now FULL DENOM formula:
//    TGraphErrors *gr0_fullDenom = gr_corr_rr_formula_FullEtaInDenom_VS_ETA_SS[ 0 ][0][0][0];
//    TGraphErrors *gr1_fullDenom = gr_corr_rr_formula_FullEtaInDenom_VS_ETA_SS[ 0 ][1][0][0];
//    TGraphErrors *gr2_fullDenom = gr_corr_rr_formula_FullEtaInDenom_VS_ETA_SS[ 0 ][2][0][0];
//    TGraphErrors *gr3_fullDenom = gr_corr_rr_formula_FullEtaInDenom_VS_ETA_OS[ 0 ][3][0][0]; // !!! Note OS!
//    TGraphErrors *gr4_fullDenom = gr_corr_rr_formula_FullEtaInDenom_VS_ETA_OS[ 0 ][4][0][0]; // !!! Note OS!

    TGraphErrors *gr0_rr_formula_fullDenom = observables_FULL_ETA_DENOM[0][0][0][0][0].getGraph( "corr_rr_FULL_ETA_DENOM_formula" );
    TGraphErrors *gr1_rr_formula_fullDenom = observables_FULL_ETA_DENOM[0][1][0][0][0].getGraph( "corr_rr_FULL_ETA_DENOM_formula" );
    TGraphErrors *gr2_rr_formula_fullDenom = observables_FULL_ETA_DENOM[0][2][0][0][0].getGraph( "corr_rr_FULL_ETA_DENOM_formula" );
    TGraphErrors *gr3_rr_formula_fullDenom = observables_FULL_ETA_DENOM[0][3][0][0][0].getGraph( "corr_rr_FULL_ETA_DENOM_formula" ); // !!! Note OS!
    TGraphErrors *gr4_rr_formula_fullDenom = observables_FULL_ETA_DENOM[0][4][0][0][0].getGraph( "corr_rr_FULL_ETA_DENOM_formula" ); // !!! Note OS!


    shiftPointX( gr0_rr_formula_fullDenom, 0.04 );
    shiftPointX( gr1_rr_formula_fullDenom, 0.04 );
    shiftPointX( gr2_rr_formula_fullDenom, 0.04 );
    shiftPointX( gr3_rr_formula_fullDenom, 0.04 );
    shiftPointX( gr4_rr_formula_fullDenom, 0.04 );

    drawGraph( gr0_rr_formula_fullDenom, 30, kGreen+1, "Pz", 1.5 );
    drawGraph( gr1_rr_formula_fullDenom, 30, kRed, "Pz", 1.5 );
    drawGraph( gr2_rr_formula_fullDenom, 30, kRed+2, "Pz", 1.5 );
    drawGraph( gr3_rr_formula_fullDenom, 30, kBlue, "Pz", 1.5 );
    drawGraph( gr4_rr_formula_fullDenom, 30, kBlue+2, "Pz", 1.5 );


    // #### now FULL DENOM direct:
    TGraphErrors *gr0_rr_direct_fullDenom = observables_FULL_ETA_DENOM[0][0][0][0][0].getGraph( "corr_rr_FULL_ETA_DENOM_direct" );
    TGraphErrors *gr1_rr_direct_fullDenom = observables_FULL_ETA_DENOM[0][1][0][0][0].getGraph( "corr_rr_FULL_ETA_DENOM_direct" );
    TGraphErrors *gr2_rr_direct_fullDenom = observables_FULL_ETA_DENOM[0][2][0][0][0].getGraph( "corr_rr_FULL_ETA_DENOM_direct" );
    TGraphErrors *gr3_rr_direct_fullDenom = observables_FULL_ETA_DENOM[0][3][0][0][0].getGraph( "corr_rr_FULL_ETA_DENOM_direct" );  // !!! Note OS!
    TGraphErrors *gr4_rr_direct_fullDenom = observables_FULL_ETA_DENOM[0][4][0][0][0].getGraph( "corr_rr_FULL_ETA_DENOM_direct" );  // !!! Note OS!

//    shiftPointX( gr0_rr_direct_fullDenom, 0.04 );
//    shiftPointX( gr1_rr_direct_fullDenom, 0.04 );
//    shiftPointX( gr2_rr_direct_fullDenom, 0.04 );
//    shiftPointX( gr3_rr_direct_fullDenom, 0.04 );
//    shiftPointX( gr4_rr_direct_fullDenom, 0.04 );

//    drawGraph( gr0_rr_direct_fullDenom, 21, kGreen+1, "L", 0.8 );
//    drawGraph( gr1_rr_direct_fullDenom, 21, kRed, "L", 0.8 );
//    drawGraph( gr2_rr_direct_fullDenom, 21, kRed+1, "L", 0.8 );
//    drawGraph( gr3_rr_direct_fullDenom, 21, kBlue, "L", 0.8 );
//    drawGraph( gr4_rr_direct_fullDenom, 21, kBlue+1, "L", 0.8 );

    gr0_rr_direct_fullDenom->SetFillColor( kGreen+1 );
    gr1_rr_direct_fullDenom->SetFillColor( kRed );
    gr2_rr_direct_fullDenom->SetFillColor( kRed+2 );
    gr3_rr_direct_fullDenom->SetFillColor( kBlue );
    gr4_rr_direct_fullDenom->SetFillColor( kBlue+2 );

    gr0_rr_direct_fullDenom->SetLineStyle(2);
    gr1_rr_direct_fullDenom->SetLineStyle(2);
    gr2_rr_direct_fullDenom->SetLineStyle(2);
    gr3_rr_direct_fullDenom->SetLineStyle(2);
    gr4_rr_direct_fullDenom->SetLineStyle(2);


    drawGraph( gr0_rr_direct_fullDenom, 1, kGreen+1, "l" );
    drawGraph( gr1_rr_direct_fullDenom, 1, kRed, "l" );
    drawGraph( gr2_rr_direct_fullDenom, 1, kRed+2, "l" );
    drawGraph( gr3_rr_direct_fullDenom, 1, kBlue, "l" );
    drawGraph( gr4_rr_direct_fullDenom, 1, kBlue+2, "l" );


    gPad->SetGrid();
    leg_nuFB_vs_eta->Draw();






    // ############## ratio_ratio BIS
    TCanvas *canv_ratio_ratio_VS_ETA_bis = new TCanvas("canv_ratio_ratio_VS_ETA_bis","canv_ratio_ratio_VS_ETA_bis",170,95,800,600);
    tuneCanvas( canv_ratio_ratio_VS_ETA );

    drawGraph( gr0_corr_rr_formula, 27, kGreen+1,   "APz" );
    drawGraph( gr1_corr_rr_formula, 20, kRed,       "Pz" );
//    drawGraph( gr2_corr_rr_formula, 24, kRed+2,     "Pz" );
    drawGraph( gr3_corr_rr_formula, 21, kBlue,      "Pz" );
//    drawGraph( gr4_corr_rr_formula, 25, kBlue+2,    "Pz" );

    drawGraph( gr0_corr_rr_direct, 1, kGreen+1, "le3" );
    drawGraph( gr1_corr_rr_direct, 1, kRed, "le3" );
//    drawGraph( gr2_corr_rr_direct, 1, kRed+2, "le3" );
    drawGraph( gr3_corr_rr_direct, 1, kBlue, "le3" );
//    drawGraph( gr4_corr_rr_direct, 1, kBlue+2, "le3" );


    drawGraph( gr0_rr_REC, 24, kGreen+1,   "Pz", 1.5 );
    drawGraph( gr1_rr_REC, 24, kRed,       "Pz", 1.5 );
//    drawGraph( gr2_rr_REC, 24, kRed+1,     "Pz", 1.5 );
    drawGraph( gr3_rr_REC, 24, kBlue,      "Pz", 1.5 );
//    drawGraph( gr4_rr_REC, 24, kBlue+1,    "Pz", 1.5 );


    drawGraph( gr0_rr_CORR, 25, kGreen+1,   "Pz", 1.5 );
    drawGraph( gr1_rr_CORR, 25, kRed,       "Pz", 1.5 );
//    drawGraph( gr2_rr_CORR, 25, kRed+1,     "Pz", 1.5 );
    drawGraph( gr3_rr_CORR, 25, kBlue,      "Pz", 1.5 );
//    drawGraph( gr4_rr_CORR, 25, kBlue+1,    "Pz", 1.5 );


    drawGraph( gr0_rr_formula_fullDenom, 30, kGreen+1, "Pz", 1.5 );
    drawGraph( gr1_rr_formula_fullDenom, 30, kRed, "Pz", 1.5 );
//    drawGraph( gr2_rr_formula_fullDenom, 30, kRed+2, "Pz", 1.5 );
    drawGraph( gr3_rr_formula_fullDenom, 30, kBlue, "Pz", 1.5 );
//    drawGraph( gr4_rr_formula_fullDenom, 30, kBlue+2, "Pz", 1.5 );





    //    return 0;







    // ### R2:
    // VS ETA:
    TCanvas *canv_R2_VS_ETA = new TCanvas("canv_R2_VS_ETA","canv_R2_VS_ETA",220,55,800,600);
    tuneCanvas( canv_sigma_VS_ETA );

    TGraphErrors *grRaa = observables[0][1][0][0][0].getGraph( "corr_R2_aa" );  //  gr_corr_R2_aa_VS_ETA   [0][ 1 ][0][0];
    TGraphErrors *grRbb = observables[0][1][0][0][0].getGraph( "corr_R2_bb" );  //  gr_corr_R2_bb_VS_ETA   [0][ 1 ][0][0];
    TGraphErrors *grRab = observables[0][1][0][0][0].getGraph( "corr_R2_ab" );  //  gr_corr_R2_ab_VS_ETA   [0][ 1 ][0][0];
    TGraphErrors *grRba = observables[0][1][0][0][0].getGraph( "corr_R2_ba" );  //  gr_corr_R2_ba_VS_ETA   [0][ 1 ][0][0];

    tuneGraphAxisLabels( grRaa );
    grRaa->GetYaxis()->SetRangeUser(0, 2);
    grRaa->SetTitle( ";#Delta#eta;#R_{2}" );

    drawGraph( grRaa, 20, kRed, "APz" );
    drawGraph( grRbb, 24, kBlue, "Pz" );
    drawGraph( grRab, 25, kMagenta, "Pz" );
    drawGraph( grRba, 5, kGray+2, "Pz" );


    TLegend *leg_R2_vs_eta = new TLegend(0.64,0.55,0.94,0.82);
    tuneLegend( leg_R2_vs_eta );
    //    leg_SIGMAFB_vs_eta->AddEntry( gr0, "K/#pi, K/#pi, formula", "p");
    //    leg_SIGMAFB_vs_eta->AddEntry( gr1, "K+/#pi+, K+/#pi+, formula", "p");
    //    leg_SIGMAFB_vs_eta->AddEntry( gr2, "K+/#pi+, K#minus/#pi#minus, formula", "p");
    leg_R2_vs_eta->AddEntry( grRaa, "Raa", "p");
    leg_R2_vs_eta->AddEntry( grRbb, "Rbb", "p");
    leg_R2_vs_eta->AddEntry( grRab, "Rab", "p");
    leg_R2_vs_eta->AddEntry( grRba, "Rba", "p");
    leg_R2_vs_eta->Draw();

    gPad->SetGrid();





    // RECONSTRUCTED R2:
    TGraphErrors *gr0_R2_aa_REC = observables[1][1][0][0][0].getGraph( "corr_R2_aa" );  // gr_corr_R2_aa_VS_ETA[ 1 ][1][0][0];
    TGraphErrors *gr1_R2_bb_REC = observables[1][1][0][0][0].getGraph( "corr_R2_bb" );  // gr_corr_R2_bb_VS_ETA[ 1 ][1][0][0];
    TGraphErrors *gr2_R2_ab_REC = observables[1][1][0][0][0].getGraph( "corr_R2_ab" );  // gr_corr_R2_ab_VS_ETA[ 1 ][1][0][0];
    TGraphErrors *gr3_R2_ba_REC = observables[1][1][0][0][0].getGraph( "corr_R2_ba" );  // gr_corr_R2_ba_VS_ETA[ 1 ][1][0][0];

    gr0_R2_aa_REC->GetYaxis()->SetRangeUser(-2, 2);

    shiftPointX( gr0_R2_aa_REC, 0.04 );
    shiftPointX( gr1_R2_bb_REC, 0.04 );
    shiftPointX( gr2_R2_ab_REC, 0.04 );
    shiftPointX( gr3_R2_ba_REC, 0.04 );

    drawGraph( gr0_R2_aa_REC, 24, kGreen+1,   "Pz", 1.5 );
    drawGraph( gr1_R2_bb_REC, 24, kRed,       "Pz", 1.5 );
    drawGraph( gr2_R2_ab_REC, 24, kRed+1,     "Pz", 1.5 );
    drawGraph( gr3_R2_ba_REC, 24, kBlue,      "Pz", 1.5 );


    // CORRECTED:
    TGraphErrors *gr0_R2_aa_CORR = observables[2][1][0][0][0].getGraph( "corr_R2_aa" );  // gr_corr_R2_aa_VS_ETA[ 2 ][1][0][0];
    TGraphErrors *gr1_R2_bb_CORR = observables[2][1][0][0][0].getGraph( "corr_R2_bb" );  // gr_corr_R2_bb_VS_ETA[ 2 ][1][0][0];
    TGraphErrors *gr2_R2_ab_CORR = observables[2][1][0][0][0].getGraph( "corr_R2_ab" );  // gr_corr_R2_ab_VS_ETA[ 2 ][1][0][0];
    TGraphErrors *gr3_R2_ba_CORR = observables[2][1][0][0][0].getGraph( "corr_R2_ba" );  // gr_corr_R2_ba_VS_ETA[ 2 ][1][0][0];

    shiftPointX( gr0_R2_aa_CORR, 0.04 );
    shiftPointX( gr1_R2_bb_CORR, 0.04 );
    shiftPointX( gr2_R2_ab_CORR, 0.04 );
    shiftPointX( gr3_R2_ba_CORR, 0.04 );

    drawGraph( gr0_R2_aa_CORR, 25, kGreen+1,   "Pz", 1.5 );
    drawGraph( gr1_R2_bb_CORR, 25, kRed,       "Pz", 1.5 );
    drawGraph( gr2_R2_ab_CORR, 25, kRed+1,     "Pz", 1.5 );
    drawGraph( gr3_R2_ba_CORR, 25, kBlue,      "Pz", 1.5 );





    // ############################################################
    // ##### new ratio-avPt VS ETA (with new expansion, July 2022):
    if(1)
    {
        TCanvas *canv_rPt_NEW_EXPANSION_VS_ETA = new TCanvas("canv_rPt_NEW_EXPANSION_VS_ETA","canv_rPt_NEW_EXPANSION_VS_ETA",170,95,800,600);
        tuneCanvas( canv_rPt_NEW_EXPANSION_VS_ETA );

        TGraphErrors *gr0_rAvPt_FullFormula = observables[0][0][0][0][0].getGraph( "corr_rPt_formula" );  //  gr_corr_rPt_full_formula_VS_ETA[ 0 ][0][0][0];
        TGraphErrors *gr1_rAvPt_FullFormula = observables[0][1][0][0][0].getGraph( "corr_rPt_formula" );  //  gr_corr_rPt_full_formula_VS_ETA[ 0 ][1][0][0];
        TGraphErrors *gr2_rAvPt_FullFormula = observables[0][2][0][0][0].getGraph( "corr_rPt_formula" );  //  gr_corr_rPt_full_formula_VS_ETA[ 0 ][2][0][0];
        TGraphErrors *gr3_rAvPt_FullFormula = observables[0][3][0][0][0].getGraph( "corr_rPt_formula" );  //  gr_corr_rPt_full_formula_VS_ETA[ 0 ][3][0][0];
        TGraphErrors *gr4_rAvPt_FullFormula = observables[0][4][0][0][0].getGraph( "corr_rPt_formula" );  //  gr_corr_rPt_full_formula_VS_ETA[ 0 ][4][0][0];


        gr0_rAvPt_FullFormula->GetYaxis()->SetRangeUser(-2, 2);

        // formula:
        drawGraph( gr0_rAvPt_FullFormula, 27, kGreen+1,   "APz" );
        drawGraph( gr1_rAvPt_FullFormula, 20, kRed,       "Pz" );
        drawGraph( gr2_rAvPt_FullFormula, 24, kRed+1,     "Pz" );
        drawGraph( gr3_rAvPt_FullFormula, 21, kBlue,      "Pz" );
        drawGraph( gr4_rAvPt_FullFormula, 25, kBlue+1,    "Pz" );

        //    drawGraph( gr0_rPt, 24, kRed, "APz" );
        //    drawGraph( gr_corr_rPt_new_formula_VS_ETA[ 1 ][0][0][0], 24, kBlue, "Pz" );
        //    drawGraph( gr_corr_rPt_new_formula_VS_ETA[ 2 ][0][0][0], 24, kMagenta, "Pz" );

        // direct:
        TGraphErrors *gr0_rAvPt_direct = observables[0][0][0][0][0].getGraph( "corr_rPt_direct" );  //gr_corr_rPt_direct_VS_ETA[ 0 ][0][0][0];
        TGraphErrors *gr1_rAvPt_direct = observables[0][1][0][0][0].getGraph( "corr_rPt_direct" );  //gr_corr_rPt_direct_VS_ETA[ 0 ][1][0][0];
        TGraphErrors *gr2_rAvPt_direct = observables[0][2][0][0][0].getGraph( "corr_rPt_direct" );  //gr_corr_rPt_direct_VS_ETA[ 0 ][2][0][0];
        TGraphErrors *gr3_rAvPt_direct = observables[0][3][0][0][0].getGraph( "corr_rPt_direct" );  //gr_corr_rPt_direct_VS_ETA[ 0 ][3][0][0];
        TGraphErrors *gr4_rAvPt_direct = observables[0][4][0][0][0].getGraph( "corr_rPt_direct" );  //gr_corr_rPt_direct_VS_ETA[ 0 ][4][0][0];

        drawGraph( gr0_rAvPt_direct, 20, kGreen+1, "Lz" );
        drawGraph( gr1_rAvPt_direct, 20, kRed,     "Lz" );
        drawGraph( gr2_rAvPt_direct, 20, kRed+1,   "Lz" );
        drawGraph( gr3_rAvPt_direct, 20, kBlue,    "Lz" );
        drawGraph( gr4_rAvPt_direct, 20, kBlue+1,  "Lz" );

        gPad->SetGrid();


        // RECONSTRUCTED:
        TGraphErrors *gr0_rPt_REC = observables[1][0][0][0][0].getGraph( "corr_rPt_formula" );  // gr_corr_rPt_full_formula_VS_ETA[ 1 ][0][0][0];
        TGraphErrors *gr1_rPt_REC = observables[1][1][0][0][0].getGraph( "corr_rPt_formula" );  // gr_corr_rPt_full_formula_VS_ETA[ 1 ][1][0][0];
        TGraphErrors *gr2_rPt_REC = observables[1][2][0][0][0].getGraph( "corr_rPt_formula" );  // gr_corr_rPt_full_formula_VS_ETA[ 1 ][2][0][0];
        TGraphErrors *gr3_rPt_REC = observables[1][3][0][0][0].getGraph( "corr_rPt_formula" );  // gr_corr_rPt_full_formula_VS_ETA[ 1 ][3][0][0];
        TGraphErrors *gr4_rPt_REC = observables[1][4][0][0][0].getGraph( "corr_rPt_formula" );  // gr_corr_rPt_full_formula_VS_ETA[ 1 ][4][0][0];

        gr0_rPt_REC->GetYaxis()->SetRangeUser(-2, 2);

        shiftPointX( gr0_rPt_REC, 0.04 );
        shiftPointX( gr1_rPt_REC, 0.04 );
        shiftPointX( gr2_rPt_REC, 0.04 );
        shiftPointX( gr3_rPt_REC, 0.04 );
        shiftPointX( gr4_rPt_REC, 0.04 );

        drawGraph( gr0_rPt_REC, 24, kGreen+1,   "Pz", 1.5 );
        drawGraph( gr1_rPt_REC, 24, kRed,       "Pz", 1.5 );
        drawGraph( gr2_rPt_REC, 24, kRed+1,     "Pz", 1.5 );
        drawGraph( gr3_rPt_REC, 24, kBlue,      "Pz", 1.5 );
        drawGraph( gr4_rPt_REC, 24, kBlue+1,    "Pz", 1.5 );


        // CORRECTED:
        TGraphErrors *gr0_rPt_CORR = observables[2][0][0][0][0].getGraph( "corr_rPt_formula" );  // gr_corr_rPt_full_formula_VS_ETA[ 2 ][0][0][0];
        TGraphErrors *gr1_rPt_CORR = observables[2][1][0][0][0].getGraph( "corr_rPt_formula" );  // gr_corr_rPt_full_formula_VS_ETA[ 2 ][1][0][0];
        TGraphErrors *gr2_rPt_CORR = observables[2][2][0][0][0].getGraph( "corr_rPt_formula" );  // gr_corr_rPt_full_formula_VS_ETA[ 2 ][2][0][0];
        TGraphErrors *gr3_rPt_CORR = observables[2][3][0][0][0].getGraph( "corr_rPt_formula" );  // gr_corr_rPt_full_formula_VS_ETA[ 2 ][3][0][0];
        TGraphErrors *gr4_rPt_CORR = observables[2][4][0][0][0].getGraph( "corr_rPt_formula" );  // gr_corr_rPt_full_formula_VS_ETA[ 2 ][4][0][0];

        shiftPointX( gr0_rPt_CORR, 0.04 );
        shiftPointX( gr1_rPt_CORR, 0.04 );
        shiftPointX( gr2_rPt_CORR, 0.04 );
        shiftPointX( gr3_rPt_CORR, 0.04 );
        shiftPointX( gr4_rPt_CORR, 0.04 );

        drawGraph( gr0_rPt_CORR, 25, kGreen+1,   "Pz", 1.5 );
        drawGraph( gr1_rPt_CORR, 25, kRed,       "Pz", 1.5 );
        drawGraph( gr2_rPt_CORR, 25, kRed+1,     "Pz", 1.5 );
        drawGraph( gr3_rPt_CORR, 25, kBlue,      "Pz", 1.5 );
        drawGraph( gr4_rPt_CORR, 25, kBlue+1,    "Pz", 1.5 );


        //    TLegend *leg_rPt_vs_eta = new TLegend(0.64,0.55,0.94,0.82);
        //    tuneLegend( leg_rPt_vs_eta );
        //    leg_rPt_vs_eta->AddEntry( gr0, "K/#pi, K/#pi, formula", "p");
        //    leg_rPt_vs_eta->AddEntry( gr1, "K+/#pi+, K+/#pi+, formula", "p");
        //    leg_rPt_vs_eta->AddEntry( gr2, "K+/#pi+, K#minus/#pi#minus, formula", "p");
        //    leg_rPt_vs_eta->AddEntry( gr_corr_rPt_formula_VS_ETA[ 0 ][0][0][0], "gen pt-K/#pi, formula", "p");
        //    leg_rPt_vs_eta->AddEntry( gr_corr_rPt_formula_VS_ETA[ 1 ][0][0][0], "rec pt-K/#pi, formula", "p");
        //    leg_rPt_vs_eta->AddEntry( gr_corr_rPt_formula_VS_ETA[ 2 ][0][0][0], "corrected pt-K/#pi, formula", "p");

    }



    // ############################################################
    // ##### R2 avPt-avPt VS ETA (new expansion, July 2022):
    if(1)
    {
        TCanvas *canv_avPtavPt_formula_VS_ETA = new TCanvas("canv_avPtavPt_formula_VS_ETA","canv_avPtavPt_formula_VS_ETA",170,95,800,600);
        tuneCanvas( canv_avPtavPt_formula_VS_ETA );


        TGraphErrors *gr0_avPtavPt_formula = observables[0][0][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 0 ][0][0][0];
        TGraphErrors *gr1_avPtavPt_formula = observables[0][1][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 0 ][1][0][0];
        TGraphErrors *gr2_avPtavPt_formula = observables[0][2][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 0 ][2][0][0];
        TGraphErrors *gr3_avPtavPt_formula = observables[0][3][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 0 ][3][0][0];
        TGraphErrors *gr4_avPtavPt_formula = observables[0][4][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 0 ][4][0][0];


        gr0_avPtavPt_formula->GetYaxis()->SetRangeUser(-2, 2);

        // formula:
        drawGraph( gr0_avPtavPt_formula, 27, kGreen+1,   "APz" );
        drawGraph( gr1_avPtavPt_formula, 20, kRed,       "Pz" );
        drawGraph( gr2_avPtavPt_formula, 24, kRed+1,     "Pz" );
        drawGraph( gr3_avPtavPt_formula, 21, kBlue,      "Pz" );
        drawGraph( gr4_avPtavPt_formula, 25, kBlue+1,    "Pz" );

        //    drawGraph( gr0_rPt, 24, kRed, "APz" );
        //    drawGraph( gr_corr_rPt_new_formula_VS_ETA[ 1 ][0][0][0], 24, kBlue, "Pz" );
        //    drawGraph( gr_corr_rPt_new_formula_VS_ETA[ 2 ][0][0][0], 24, kMagenta, "Pz" );

        // direct:
        TGraphErrors *gr0_avPtavPt_direct = observables[0][0][0][0][0].getGraph( "avPtF_avPtB_direct" );  // gr_avPtavPt_direct_VS_ETA[ 0 ][0][0][0];
        TGraphErrors *gr1_avPtavPt_direct = observables[0][1][0][0][0].getGraph( "avPtF_avPtB_direct" );  // gr_avPtavPt_direct_VS_ETA[ 0 ][1][0][0];
        TGraphErrors *gr2_avPtavPt_direct = observables[0][2][0][0][0].getGraph( "avPtF_avPtB_direct" );  // gr_avPtavPt_direct_VS_ETA[ 0 ][2][0][0];
        TGraphErrors *gr3_avPtavPt_direct = observables[0][3][0][0][0].getGraph( "avPtF_avPtB_direct" );  // gr_avPtavPt_direct_VS_ETA[ 0 ][3][0][0];
        TGraphErrors *gr4_avPtavPt_direct = observables[0][4][0][0][0].getGraph( "avPtF_avPtB_direct" );  // gr_avPtavPt_direct_VS_ETA[ 0 ][4][0][0];

        drawGraph( gr0_avPtavPt_direct, 20, kGreen+1, "Lz" );
        drawGraph( gr1_avPtavPt_direct, 20, kRed,     "Lz" );
        drawGraph( gr2_avPtavPt_direct, 20, kRed+3,   "Lz" );
        drawGraph( gr3_avPtavPt_direct, 20, kBlue,    "Lz" );
        drawGraph( gr4_avPtavPt_direct, 20, kBlue+1,  "Lz" );

        gPad->SetGrid();


        // RECONSTRUCTED:
        TGraphErrors *gr0_avPtavPt_REC = observables[1][0][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 1 ][0][0][0];
        TGraphErrors *gr1_avPtavPt_REC = observables[1][1][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 1 ][1][0][0];
        TGraphErrors *gr2_avPtavPt_REC = observables[1][2][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 1 ][2][0][0];
        TGraphErrors *gr3_avPtavPt_REC = observables[1][3][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 1 ][3][0][0];
        TGraphErrors *gr4_avPtavPt_REC = observables[1][4][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 1 ][4][0][0];

        gr0_avPtavPt_REC->GetYaxis()->SetRangeUser(-2, 2);

        shiftPointX( gr0_avPtavPt_REC, 0.04 );
        shiftPointX( gr1_avPtavPt_REC, 0.04 );
        shiftPointX( gr2_avPtavPt_REC, 0.04 );
        shiftPointX( gr3_avPtavPt_REC, 0.04 );
        shiftPointX( gr4_avPtavPt_REC, 0.04 );

        drawGraph( gr0_avPtavPt_REC, 24, kGreen+1,   "Pz", 1.5 );
        drawGraph( gr1_avPtavPt_REC, 24, kRed,       "Pz", 1.5 );
        drawGraph( gr2_avPtavPt_REC, 24, kRed+3,     "Pz", 1.5 );
        drawGraph( gr3_avPtavPt_REC, 24, kBlue,      "Pz", 1.5 );
        drawGraph( gr4_avPtavPt_REC, 24, kBlue+1,    "Pz", 1.5 );


        // CORRECTED:
        TGraphErrors *gr0_avPtavPt_CORR = observables[2][0][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 2 ][0][0][0];
        TGraphErrors *gr1_avPtavPt_CORR = observables[2][1][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 2 ][1][0][0];
        TGraphErrors *gr2_avPtavPt_CORR = observables[2][2][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 2 ][2][0][0];
        TGraphErrors *gr3_avPtavPt_CORR = observables[2][3][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 2 ][3][0][0];
        TGraphErrors *gr4_avPtavPt_CORR = observables[2][4][0][0][0].getGraph( "avPtF_avPtB_formula" );  // gr_avPtavPt_formula_VS_ETA[ 2 ][4][0][0];

        gr0_avPtavPt_REC->GetYaxis()->SetRangeUser(-2, 2);

        shiftPointX( gr0_avPtavPt_CORR, 0.04 );
        shiftPointX( gr1_avPtavPt_CORR, 0.04 );
        shiftPointX( gr2_avPtavPt_CORR, 0.04 );
        shiftPointX( gr3_avPtavPt_CORR, 0.04 );
        shiftPointX( gr4_avPtavPt_CORR, 0.04 );

//        gr1_avPtavPt_CORR->SetLineStyle(2);
//        gr2_avPtavPt_CORR->SetLineStyle(9);
        drawGraph( gr0_avPtavPt_CORR, 25, kGreen+1,   "Pz", 1.5 );
        drawGraph( gr1_avPtavPt_CORR, 25, kRed,       "PLz", 1.5 );
        drawGraph( gr2_avPtavPt_CORR, 25, kRed+3,     "PLz", 1.5 );
        drawGraph( gr3_avPtavPt_CORR, 25, kBlue,      "Pz", 1.5 );
        drawGraph( gr4_avPtavPt_CORR, 25, kBlue+1,    "Pz", 1.5 );
    }




    return 0;





    timer2.Stop();
    Double_t rtime2 = timer2.RealTime();
    Double_t ctime2 = timer2.CpuTime();

    printf("End of drawing: RealTime=%f seconds, CpuTime=%f seconds\n",rtime2,ctime2);





//    // ###### Draw Graph 2D
//    TCanvas *canv_ratio_ratio_VS_ETA_Graph2D = new TCanvas("canv_ratio_ratio_VS_ETA_Graph2D","canv_ratio_ratio_VS_ETA_Graph2D",140,75,800,600);
//    TGraph2D *gr2D_0_rr_GEN = observables[0][0][0][0][0].getGraph2D( "corr_rr_formula" );
//    gr2D_0_rr_GEN->Draw("surf1");
//    gr2D_0_rr_GEN->SetMarkerStyle(24);
//    gr2D_0_rr_GEN->SetMarkerColor( kRed );
//    gr2D_0_rr_GEN->Draw("same P");






    return 0;











}






