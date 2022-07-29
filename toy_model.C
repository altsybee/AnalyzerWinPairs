#include "TTree.h"
#include "TFile.h"

#include "TRandom.h"
#include "TString.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"

//#include "WinPairFilling.h"
#include "/Users/macbookpro/work/AnalyzerWinPairs/WinPairFilling.h"

#include <iostream>
#include <fstream>
#include <string>






// #################
// function blastwave
Double_t fitBW(Double_t *x, Double_t *par)
{
    // taken from https://arxiv.org/pdf/1303.0737.pdf (ALICE)
    double pT = x[0];

    double norm = par[0];

    double Tkin = par[1]; //0.15;
    double beta = par[2]; //0.5;
    double n = par[3];
    double m = par[4];
    int use_transf_to_betaS = par[5];


    double mT = sqrt(m*m+pT*pT);

    double R = 1.;


    // loop for r:
    int nSteps = 200;
    double dr = R/nSteps;
    // if <betaT> is provided, calc <betaMax> (=<betaSurface>)
    if(use_transf_to_betaS)
    {
        //        double transform_coeff = 0;
        //        for ( int i = 0; i < nSteps; i++ )
        //        {
        //            double r = dr*i;
        //            transform_coeff += TMath::Power(r/R, n);
        //        }
        //        transform_coeff /= nSteps;
        //        cout << "n=" << n << ", transform_coeff = " << transform_coeff << ", betaBefore=" << beta;
        //        beta /= transform_coeff;
        //        cout << ", betaAfter=" << beta << endl;
        //        int aa;
        //        cin >> aa;
        beta = beta*(2+n)/2;
    }
    double res = 0;
    for ( int i = 0; i < nSteps; i++ )
    {
        double r = dr*i;
        double rho = TMath::ATanH( TMath::Power(r/R, n) * beta);

        double x_I0 = pT * TMath::SinH( rho ) / Tkin;
        double I0 = TMath::BesselI0( x_I0 );

        if(0)if(i==50)cout << "rho = " << rho << ", i=" << i << ", r/R=" << r/R << ", pT=" << pT
                           << ", norm = " << par[0]
                           << ", T = " << par[1]
                           << ", beta = " << par[2]
                           << endl;
        double x_K1 = mT * TMath::CosH( rho ) / Tkin;
        double K1 = TMath::BesselK1( x_K1 );

        res += r * dr * mT * I0 * K1;
    }

    res *= pT;

    return /*10000 **/ norm * res;
}

// from Table 5 in 2013_pi_K_p_ALICE_PbPb_1303.0737.pdf:

//        <betaT>                 <Tkin>                  n
//0–5%   0.651 ± 0.004 ± 0.020  0.095 ± 0.004 ± 0.010   0.712 ± 0.019 ± 0.086
//5–10%  0.646 ± 0.004 ± 0.023  0.097 ± 0.003 ± 0.011	0.723 ± 0.019 ± 0.116
//10–20% 0.639 ± 0.004 ± 0.022	0.099 ± 0.004 ± 0.011   0.738 ± 0.020 ± 0.118
//20–30% 0.625 ± 0.004 ± 0.025	0.101 ± 0.004 ± 0.012   0.779 ± 0.022 ± 0.133
//30–40% 0.604 ± 0.005 ± 0.022	0.106 ± 0.004 ± 0.012	0.841 ± 0.025 ± 0.168
//40–50% 0.574 ± 0.005 ± 0.016	0.112 ± 0.004 ± 0.013	0.944 ± 0.029 ± 0.142
//50–60% 0.535 ± 0.007 ± 0.018	0.118 ± 0.004 ± 0.014	1.099 ± 0.038 ± 0.187
//60–70% 0.489 ± 0.008 ± 0.024	0.129 ± 0.005 ± 0.017	1.292 ± 0.052 ± 0.194
//70–80% 0.438 ± 0.011 ± 0.039	0.139 ± 0.005 ± 0.027	1.578 ± 0.081 ± 0.205
//80–90% 0.357 ± 0.016 ± 0.084	0.151 ± 0.006 ± 0.044	2.262 ± 0.191 ± 0.498

double table5_arr_avBetaT[] = {
    0.651,
    0.646,
    0.639,
    0.625,
    0.604,
    0.574,
    0.535,
    0.489,
    0.438,
    0.357,
};

double table5_arr_Tkin[] = {
    0.095,
    0.097,
    0.099,
    0.101,
    0.106,
    0.112,
    0.118,
    0.129,
    0.139,
    0.151,
};

double table5_par_n[] = {
    0.712,
    0.723,
    0.738,
    0.779,
    0.841,
    0.944,
    1.099,
    1.292,
    1.578,
    2.262,
};


const int nSubsamples = 3;//10;//3;//15;//3;//10;//8;//20;//4;//20;//15;//20;//15;//3;//20;//5;//3;//15;//8;//30;


const int nEtaBins = 6;//8;
const int nPhiBins = 6;//8;//8;//8;

double eRange[] = { -0.8, 0.8 };


//const double eRange = 0.8; //1.0;//1.4;//1.6;//0.8;

const int nPtBins = 1;
double ptRange[] = { 0., 100 };




//double ptmin[nPtBins] = { 0 }; //0.2 };//0.1,  };//0.2 }; //0.3, };//0.2 };
//double ptmax[nPtBins] = { 100 }; //2.0,};// 2.0 }; //1.5,};// 5.0 };


//
const int nPartTypes = 1;//5;//3;//4;
int arrPartTypes[nPartTypes][4] =
{ // F,B, X,Y
  //  { /*421*/0, 0, 0, 0 },         // 0
  //  { 0, 0, 0, 0 },         // 1
  //  { 0, 0, 0, 0 },         // 2
  //  { 321, 321, 211, 211 }, // 6

  { 321, 321, 211, 211 }, // 7
//  { 321, 321, 211, 211 }, // 8
//  { 321, 321, 211, 211 }, // 9
//  { 321, 321, 211, 211 }, // 8
//  { 321, 321, 211, 211 }, // 9

  //  { 0, 321, 0, 211 }, // 7
  //  { 0, 321, 0, 211 }, // 8
  //  { 0, 321, 0, 211 }, // 9
};

int arrCharges[nPartTypes][4] =
{ // F,B, X,Y
//  {  0,  0,  0,  0 },
//  { +1, +1, +1, +1 },
//  { -1, -1, -1, -1 },
  { +1, -1, +1, -1 },
//  { -1, +1, -1, +1 },
};


const int nCW = 1;
//const int nCentrBins[] = { 1 };//10 };//9 };
const int nCentrBins = 1;//10 };//9 };
const double cBinWidths[] = { 100 }; // in %



TString strAnLevels[] = { "GEN", "REC", "CORRECTED" };
//TString strAnLevels[] = { "REC" };
const int N_AN_LEVELS = sizeof(strAnLevels)/sizeof(*strAnLevels);


//bool whichHistos[] = { true, true, true };   // flags for hAllWins, hDetaDphi, hAllEtaDphi
bool whichHistos[] = { false, true, true };   // flags for hAllWins, hDetaDphi, hAllEtaDphi



// analyzer objects:
WinPairWrapper winPairWrapper[N_AN_LEVELS][nPartTypes][nCW][ nCentrBins ][nPtBins];
WinPairWrapper winPairWrapper_FULL_ETA_DENOM[N_AN_LEVELS][nPartTypes][nCW][ nCentrBins ][nPtBins];
WinPairWrapper winPairWrapper_FULL_ETA_NUM_AND_DENOM[N_AN_LEVELS][nPartTypes][nCW][ nCentrBins ][nPtBins];


// ########################
void toy_model( int _nEv = 2e3 //25e3 //15000 //25000 //1e3 //800 //50000 //25000 //25000 //3e2 //250000
        )
{


    TFile *fileEff = new TFile("file_EFF_in_cBins_from_HIJING_FOR_NN_CORRECTIONS_9ptBins_recPrim_to_kine.root");
    TH1D *histEff = 0x0;
    if( fileEff )
    {
        histEff = (TH1D*) fileEff->Get( "cBin0" );
        if( histEff )
        {
            cout << "histEff->GetNbinsX() = " << histEff->GetNbinsX() << endl;
            if(0)
            {
                histEff->SetBinContent( 1, 0.32 );
                histEff->SetBinContent( 2, 0.38 );
                histEff->SetBinContent( 3, 0.45 );
                histEff->SetBinContent( 4, 0.57 );
                histEff->SetBinContent( 5, 0.32 );
                histEff->SetBinContent( 6, 0.48 );
                histEff->SetBinContent( 7, 0.60 );
                histEff->SetBinContent( 8, 0.75 );
                histEff->SetBinContent( 9, 0.68 );
            }

        }

    }

    //    histEff->DrawCopy();
    //    return;


    TF1 *fSlopeEta = new TF1( "slopeEta", "[0]+x*[1]", -1.5, 1.5 );
    fSlopeEta->SetParameters( 3, -1 );

//    fSlopeEta->Draw();
//    return;


    TList* fOutputWinPairsLists[nPartTypes];       //! output list
    for ( int iType = 0; iType < nPartTypes; iType++)
    {
        TString strPID = Form( "pid_%d_%d_%d_%d_charge_%d_%d_%d_%d",
                               arrPartTypes[iType][0], arrPartTypes[iType][1], arrPartTypes[iType][2], arrPartTypes[iType][3],
                arrCharges[iType][0], arrCharges[iType][1], arrCharges[iType][2], arrCharges[iType][3] );

        fOutputWinPairsLists[iType] = new TList();
        //        fOutputWinPairsLists[iType]->SetName( "histForAnalysis_"+strPID );
        fOutputWinPairsLists[iType]->SetName( "list_PIDcorr_"+strPID );
        //        fOutputList->Add( fOutputWinPairsLists[iType] );
    }




    //    cout << "### nVars = " << nVars << endl;






    TH2I *histAccMap = new TH2I( "histAccMap", "histAccMap;#eta;#varphi", nEtaBins, eRange[0], eRange[1], nPhiBins, 0, TMath::TwoPi() );
    for( int e1 = 0; e1 < nEtaBins; e1++ )
        for( int p1 = 0; p1 < nPhiBins; p1++ )
//            histAccMap->SetBinContent( e1+1, p1+1, 1 );
            histAccMap->SetBinContent( e1+1, p1+1, gRandom->Uniform() > 0.4 ? 1 : 0 );

//    histAccMap->Fill( -0.6, 4.8, -1 );
//    histAccMap->Fill( 0.2, 2.5, -1 );



    // ### FB:
    for ( int iAnLevel = 0; iAnLevel < N_AN_LEVELS; ++iAnLevel )
        for ( int iType = 0; iType < nPartTypes; iType++)
            for ( int iCW = 0; iCW < nCW; iCW++)
                for ( int cBin = 0; cBin < /*nCentrBins[iCW]*/nCentrBins; ++cBin )
                    for ( int iPt = 0; iPt < nPtBins; ++iPt )
                    {
                        winPairWrapper[iAnLevel][iType][iCW][cBin][iPt].setup( strAnLevels[iAnLevel], cBin, nSubsamples, arrPartTypes[iType], arrCharges[iType], nEtaBins, nPhiBins, eRange, ptRange, whichHistos, fOutputWinPairsLists[iType] );
                        if ( strAnLevels[iAnLevel].Contains( "REC" ) || strAnLevels[iAnLevel].Contains( "CORRECTED" ) )   //iAnLevel == 1 ) // TEST ACC MAP FOR RECO!!
                        {
                            winPairWrapper[iAnLevel][iType][iCW][cBin][iPt].hAccMap = histAccMap;
                            fOutputWinPairsLists[iType]->Add( winPairWrapper[iAnLevel][iType][iCW][cBin][iPt].hAccMap );
                            fOutputWinPairsLists[iType]->Add( winPairWrapper[iAnLevel][iType][iCW][cBin][iPt].QA_hist_n[0] );
                            fOutputWinPairsLists[iType]->Add( winPairWrapper[iAnLevel][iType][iCW][cBin][iPt].QA_hist_pt[0] );

                        }

                        // ##### FULL ETA FOR DENOM:
                        winPairWrapper_FULL_ETA_DENOM[iAnLevel][iType][iCW][cBin][iPt].setup( strAnLevels[iAnLevel]+"_FULL_ETA_DENOM", cBin, nSubsamples, arrPartTypes[iType], arrCharges[iType], nEtaBins, nPhiBins, eRange, ptRange, whichHistos, fOutputWinPairsLists[iType], true );

                        // ##### FULL ETA FOR both Num and Denom: to calc nu_dyn!
                        bool whichHistosFullAcceptance[] = { true, false, false };   // flags for hAllWins, hDetaDphi, hAllEtaDphi
                        winPairWrapper_FULL_ETA_NUM_AND_DENOM[iAnLevel][iType][iCW][cBin][iPt].setup( strAnLevels[iAnLevel]+"_FULL_ETA_NUM_AND_DENOM", cBin, nSubsamples, arrPartTypes[iType], arrCharges[iType], 1, 1, eRange, ptRange, whichHistosFullAcceptance, fOutputWinPairsLists[iType], true );

                    } // end of pT loop


    //    double multRanges[] = { 3199.5, 283.382, 195.397, 133.02, 88.058, 55.5156, 32.8752, 18.0094, 8.92059, 3.94449,  };
    double multRanges[] = { 3199.5, -1  }; // for pp!

    //    int nCbins = sizeof(multRanges)/sizeof(*multRanges);
    int nCbins = nCentrBins; // for pp!

    //        return;

    //
    TDatabasePDG *db_PDG = new TDatabasePDG;

    TH1D *h_QA_event = new TH1D("h_QA_event",";n tracks;Entries", 5, -0.5, 4.5 );

    h_QA_event->GetXaxis()->SetBinLabel(1,"noTracksInEta1");
    h_QA_event->GetXaxis()->SetBinLabel(2,"Analyzed");

    // basic plots for tree branches:
    //    TH1D *h_bImp = new TH1D("h_bImp", ";b, fm;entries", 250, 0, 25 );
    //    TH1D *h_bImpTriggered = new TH1D("h_bImpTriggered", ";b, fm;entries", 250, 0, 25 );
    TH1D *h_nparticles_gen = new TH1D("h_nparticles_gen", ";npart;entries", 2501, -0.5, 2500.5);
    TH1D *h_nparticles_gen_in01win = new TH1D("h_nparticles_gen_in01win", ";npart;entries", 2501, -0.5, 2500.5);

    //    TH1D *h_px = new TH1D("h_px", ";px;entries", 400, -10, 10);
    //    TH1D *h_py = new TH1D("h_py", ";py;entries", 400, -10, 10);
    //    TH1D *h_pz = new TH1D("h_pz", ";pz;entries", 400, -10, 10);
    TH1D *h_pid = new TH1D("h_pid", ";pid;entries", 10001, -5000.5, 5000.5);
    TH1D *h_pid_REC = new TH1D("h_pid_REC", ";pid_REC;entries", 10001, -5000.5, 5000.5);



    // ### SIM:
    int _nYBins = 800;
    double etaMin = -8;
    double etaMax = 8;
    int _nPtBins = 500;
    double ptMin = 0;
    double ptMax = 5;

    TH1D *h_Pt = new TH1D("h_Pt", "h_Pt", _nPtBins, ptMin, ptMax);
    TH1D *h_Eta = new TH1D("h_Eta", "h_Eta", _nYBins, etaMin, etaMax );
    TH1D *h_Phi = new TH1D("h_Phi", "h_Phi", 400, -2*TMath::TwoPi(), 2*TMath::TwoPi());
    TH1D *h_Rapidity = new TH1D("h_Rapidity", "h_Rapidity", _nYBins, etaMin, etaMax );

    TH1D *h_Pt_REC = new TH1D("h_Pt_REC", "h_Pt_REC", _nPtBins, ptMin, ptMax);
    TH1D *h_Eta_REC = new TH1D("h_Eta_REC", "h_Eta_REC", _nYBins, etaMin, etaMax );
    TH1D *h_Phi_REC = new TH1D("h_Phi_REC", "h_Phi_REC", 400, -2*TMath::TwoPi(), 2*TMath::TwoPi());
    TH1D *h_Rapidity_REC = new TH1D("h_Rapidity_REC", "h_Rapidity_REC", _nYBins, etaMin, etaMax );

    TH1D *h_Pt_REC_pi = new TH1D("h_Pt_REC_pi", "h_Pt_REC_pi", _nPtBins, ptMin, ptMax);
    TH1D *h_Pt_REC_K = new TH1D("h_Pt_REC_K", "h_Pt_REC_K", _nPtBins, ptMin, ptMax);

    TH1D *h_weights = new TH1D("h_weights", "h_weights", _nPtBins, ptMin, ptMax);


    TH1D *h_Pt_in_eta16 = new TH1D("h_Pt_in_eta16", "h_Pt_in_eta16", _nPtBins, ptMin, ptMax);       // create your histogra
    TH1D *h_Eta_in_pt01_10 = new TH1D("h_Eta_in_pt01_10", "h_Eta_in_pt01_10", _nYBins, etaMin, etaMax );

    TH1D *h_nTracksForMultCentrality = new TH1D("h_nTracksForMultCentrality",";n tracks;Entries", 3201, -0.5, 3200.5 ); //302,-1.5,300.5);

    TH1D *h_nTracksInEtaCuts = new TH1D("h_nTracksInEtaCuts",";n tracks;Entries", 3201, -0.5, 3200.5 ); //302,-1.5,300.5);
    TH1D *h_nTracksInEtaCuts_REC = new TH1D("h_nTracksInEtaCuts_REC",";n tracks REC;Entries", 3201, -0.5, 3200.5 ); //302,-1.5,300.5);

    TH1D *h_QA_centr_bins = new TH1D("h_QA_centr_bins",";centrality bin;Entries", 101, -0.5, 100.5 );

    TH1D *h_QA_collectivity_degree = new TH1D("h_QA_collectivity_degree",";collectivity degree;Entries", 150, 0, 1.5 );

    const int MAX_N_PART = 10000;
    //    Double_t impact_b;
    Int_t npart;
    Bool_t emptyEv;
    //    Double_t tree_px[MAX_N_PART], tree_py[MAX_N_PART], tree_pz[MAX_N_PART];
    Double_t arr_pt[MAX_N_PART], arr_eta[MAX_N_PART], arr_phi[MAX_N_PART];
    //    Int_t pdgcode[MAX_N_PART];
    //    Int_t track_charge[MAX_N_PART];
    Int_t arr_charge[MAX_N_PART];
    Int_t arr_pdg[MAX_N_PART];

    //    double eff = 0.8;



    // ##### BW for spectra
    double min_for_BW = 0.2;//0.1;
    double max_for_BW = 2.0;//5.0;
    TF1 *my_func_BW = new TF1( "my_BW", fitBW, min_for_BW, max_for_BW, /*5*/6 );
    my_func_BW->SetNpx(12);
    my_func_BW->SetParameter( 0, 1.);
    my_func_BW->FixParameter( 5, 1); // use or not transformation from <betaT> to betaS within function

    TH1D *hist_fluctT_mean_pt_pions     = new TH1D( "hist_fluctT_mean_pt_pions", ";#LTp_{T}#GT;entries", 400, 0.2, 1.8 );
    TH1D *hist_fluctT_mean_pt_kaons     = new TH1D( "hist_fluctT_mean_pt_kaons", ";#LTp_{T}#GT;entries", 400, 0.2, 1.8 );
    TH1D *hist_fluctT_mean_pt_protons   = new TH1D( "hist_fluctT_mean_pt_protons", ";#LTp_{T}#GT;entries", 400, 0.2, 1.8 );
    TH2D *hist_fluctT_mean_pt_pions_vs_Tkin      = new TH2D( "hist_fluctT_mean_pt_pions_vs_Tkin",  ";T_{kin}, GeV;#LTp_{T}#GT;entries", 100, 0.0, 0.2, 200, 0.2, 1.8 );
    TH2D *hist_fluctT_mean_pt_pions_vs_betaT     = new TH2D( "hist_fluctT_mean_pt_pions_vs_betaT", ";#beta, GeV;#LTp_{T}#GT;entries",   80, 0.0, 0.8, 200, 0.2, 1.8 );


    double ptMinForInt = 0;
    double ptMaxForInt = 5;
    //    bool doCheckWithHist = false;
    //    bool doDrawBW = false;
    bool doCalcMeanTF1 = true;//false;


    TStopwatch timer;
    timer.Start();

    // ### event loop
    int nEvents = _nEv;//4e3;//30000;//3e4;//1e6; //chain->GetEntries();
    for (int iEv = 0; iEv < nEvents; iEv++)
    {
        if( iEv%10 == 0 )
            cout << "generating " << (int)iEv << "/" << nEvents << " event... \r"; cout.flush();

        //        cout << "event " << iEv << ": npart = " << npart << ", impact_b = " << impact_b << ", isEmpty = " << emptyEv << endl;

        //        h_bImp->Fill(impact_b);

        // number of particles in this event
        //        int nPart = gRandom->Gaus(500,25);
        //        int nPart = gRandom->Gaus(500,20);//0.00001);
        //        int nPart = gRandom->Gaus(50,2);//0.00001);
        //        int nPart = 80;//gRandom->Gaus(50/2,2/2) * 2; // to have even an number
        //        int nPart = gRandom->Gaus(80/2,4/2) * 2; // to have even an number
        int nPart = gRandom->Gaus( 450, 20 );//5 );
        //        int nPart = gRandom->Gaus( 15, 1 );//5 );
        //        const int nPions = gRandom->Poisson(400);
        //        const int nKaons = 200;//6;
        //        int nPart = nPions+nKaons;
        //        int nPart = gRandom->Gaus(750,20);
        //        int nPart = gRandom->Gaus(1000,20);
        //        int nPart = gRandom->Gaus(1000,1);
        //        int nPart = gRandom->Gaus(1000,1);
        //        int nPart = gRandom->Gaus(500,15);
        //        int nPart = gRandom->Gaus(500,50);
        if ( nPart < 1 )
            continue;

        //        double collectivity_degree = 0.8;
        //        double collectivity_degree = gRandom->Gaus(0.8, 0.02);
        //        const double av_collectivity = 0.7; //0.65;
        const double av_collectivity = 0.65;//0.8; //0.65;
        //        double collectivity_degree = av_collectivity;
        //        double collectivity_degree = gRandom->Gaus( av_collectivity, 0.02 );
        double collectivity_degree = gRandom->Gaus( av_collectivity, 0.03 ); //0.03 );//0.04 );
        //        double collectivity_degree = gRandom->Uniform( 0.75, 0.85 );
        //        double collectivity_degree = 0.8;
        //        double collectivity_degree = gRandom->Gaus(0.8, -(nPart-500.)/500. );

        if( collectivity_degree<0 || collectivity_degree>1 )
            continue;

        h_QA_collectivity_degree->Fill( collectivity_degree );

        //                double meanPt = 0.5 * (1 - 2*(collectivity_degree-0.65) ) ;
        //        double meanPt = 0.6;
        double meanPt = 0.6 * (1 - 4*( collectivity_degree - av_collectivity) );
        //        double meanPt = gRandom->Gaus( 0.6, 0.08 );
        //        while ( meanPt < 0.2 )
        //            meanPt = gRandom->Gaus( 0.6, 0.08 );

        int counterNpartIn01Win = 0;

        if(0) // use BW or not
        {
            // set Tkin and betaT (Table 5 in pi,K,p in Pb-Pb paper!)
            const int centrBinBW = 9;//0;//9;
            double Tkin = table5_arr_Tkin[centrBinBW] ;// - 0.4*(collectivity_degree - av_collectivity); //+ gRandom->Gaus(0,0.005)
            double betaT = table5_arr_avBetaT[centrBinBW] - 0.5*(collectivity_degree - av_collectivity); //+ gRandom->Gaus(0,0.02);
            double nBW = table5_par_n[centrBinBW];
            if ( iEv == 0 )
                cout << "Tkin = " << Tkin << ", betaT = " << betaT << ", nBW = " << nBW << endl;

            my_func_BW->SetParameter( 1, Tkin );//+ gRandom->Gaus(0,0.005) ); // Tkin
            my_func_BW->SetParameter( 2, betaT ); // beta // 1./(2./3) is needed to convert from <betaT> back to betaSurface
            my_func_BW->SetParameter( 3, nBW ); // n
            my_func_BW->SetParameter( 0, 1.);

            //pions
            my_func_BW->SetParameter( 4, 0.140);

            if( doCalcMeanTF1 && iEv<200)
            {
                double meanPt_pions_from_TF1 = my_func_BW->Mean( ptMinForInt, ptMaxForInt );
                //        cout << "Check: MEAN pT FROM BW for pions from TF1 = " << meanPt_pions_from_TF1 << endl;
                //            gr_mean_pt_pions->SetPoint( i, cBinMiddle, meanPt_pions_from_TF1 );
                hist_fluctT_mean_pt_pions->Fill( meanPt_pions_from_TF1 );
                hist_fluctT_mean_pt_pions_vs_Tkin ->Fill( Tkin, meanPt_pions_from_TF1 );
                hist_fluctT_mean_pt_pions_vs_betaT->Fill( betaT, meanPt_pions_from_TF1 );
            }
        }


        double prevEta = 0;
        double prevPhi = 0;
        double prevPt = 0;
        double prevPdg = 0;
        double prevCharge = 0;

        if(0) // NEW MODEL Fall 2021 for NUCLEUS 2021:
        {
            // ### track pre-loop:
            //        int nParticlesForTrigger = 0;
            //        int nTracksForMultCentrality = 0;
            //        cout << nPart << endl;


            // pre-loop to define pdgs
            double arr_eta_K_plus[MAX_N_PART];
            double arr_eta_pi_plus[MAX_N_PART];
            double arr_eta_K_minus[MAX_N_PART];
            double arr_eta_pi_minus[MAX_N_PART];

            int nK_plus = 0;
            int npi_plus = 0;

            for (int i = 0; i < nPart/2; i++) // create only positive particles
            {
                int pdg = gRandom->Uniform() > collectivity_degree ? 321 : 211;
                //                        int pdg = i < 0.2*(nPart/2) ? 321 : 211;
                //            int pdg = i < (1-collectivity_degree)*(nPart/2) ? 321 : 211;
                //            arr_pdg[i] = pdg;
                if (pdg == 321)
                    nK_plus++;
                else
                    npi_plus++;
            }

            // loop to define positions of KAONS in eta
            for (int i = 0; i < nK_plus; i++)
            {
                //                double eta = gRandom->Uniform( -2., 2. );
                double eta = gRandom->Uniform( -1, 1 );
                if(0) // if we want repulsion for pairs of KAONS
                {
                    if ( i%2 == 0 )
                        prevEta = eta;
                    else
                    {
                        double sepEta = 0.2 + fabs( gRandom->Gaus(0,0.25) );  // REPULSION of a pair of same-sign charges
                        while( fabs(eta-prevEta) < sepEta )
                        {
                            //                            cout << "check dEta = " << fabs(eta-prevEta) << endl;
                            eta = gRandom->Uniform( -2., 2. );
                        }
                    }
                }
                arr_eta_K_plus[i] = eta;

                // now opp. sign pairing:
                double randWhere = gRandom->Uniform();
                double etaPartner = 0;
                if(0) // if we want opp sign SR correlations for some fraction of pairs:
                {
                    if ( randWhere < 0.25 ) // opposite-sign SR correlation (25% prob.)
                        etaPartner = gRandom->Gaus(eta,0.25);
                    else // somewhere else
                        etaPartner = gRandom->Uniform( -2., 2. );
                }
                else // always random eta
                    etaPartner = gRandom->Uniform( -2., 2. );


                arr_eta_K_minus[i] = etaPartner;
            }

            // loop to define positions of PIONS in eta
            for (int i = 0; i < npi_plus; i++)
            {
                double eta = gRandom->Uniform( -2., 2. );
                if(0) // if we want repulsion for pairs of PIONS
                {
                    if ( i%2 == 0 )
                        prevEta = eta;
                    else
                    {
                        double sepEta = 0.2 + fabs( gRandom->Gaus(0,0.3) );  // REPULSION of a pair of same-sign charges
                        while( fabs(eta-prevEta) < sepEta )
                        {
                            //                            cout << "check dEta = " << fabs(eta-prevEta) << endl;
                            eta = gRandom->Uniform( -2., 2. );
                        }
                    }
                }
                arr_eta_pi_plus[i] = eta;

                // now opp. sign pairing:
                double randWhere = gRandom->Uniform();
                double etaPartner = 0;
                if(0) // if we want opp sign SR correlations for some fraction of pairs:
                {
                    if ( randWhere < 0.25 ) // opposite-sign SR correlation (25% prob.)
                        etaPartner = gRandom->Gaus(eta,0.4);
                    else // somewhere else
                        etaPartner = gRandom->Uniform( -2., 2. );
                }
                else // always random eta
                    etaPartner = gRandom->Uniform( -2., 2. );


                arr_eta_pi_minus[i] = etaPartner;
            }

            // now combine all K and pi
            int counter = 0;
            // K:
            for (int i = 0; i < nK_plus; i++)
            {
                arr_eta[counter] = arr_eta_K_plus[i];
                arr_pdg[counter] = 321;
                arr_charge[counter] = +1;
                arr_pt[counter] = 0.6;
                arr_phi[counter] = 0.2;
                counter++;

                arr_eta[counter] = arr_eta_K_minus[i];
                arr_pdg[counter] = -321;
                arr_charge[counter] = -1;
                arr_pt[counter] = 0.6;
                arr_phi[counter] = 0.2;
                counter++;
            }
            // pi:
            for (int i = 0; i < npi_plus; i++)
            {
                arr_eta[counter] = arr_eta_pi_plus[i];
                arr_pdg[counter] = 211;
                arr_charge[counter] = +1;
                arr_pt[counter] = 0.6;
                arr_phi[counter] = 0.2;
                counter++;

                arr_eta[counter] = arr_eta_pi_minus[i];
                arr_pdg[counter] = -211;
                arr_charge[counter] = -1;
                arr_pt[counter] = 0.6;
                arr_phi[counter] = 0.2;
                counter++;
            }
        } // end of new model Fall 2021


        // ########
        else if(1) // PREVIOUS MODEL
        {
            for (int i = 0; i < nPart; i++)
            {
                //            if ( track_charge[itrack] == 0 )
                //                continue;

                int charge = gRandom->Uniform() > 0.5 ? +1 : -1;
                //                int charge = i%2 ? +1 : -1;
                //            double pt = gRandom->Exp(0.5); // TMath::Sqrt(px*px + py*py);
                double pt = gRandom->Exp( meanPt ); // TMath::Sqrt(px*px + py*py);
                //            double pt = my_func_BW->GetRandom();
                //            double eta = gRandom->Uniform(-2, 2); // 0.5 * TMath::Log( ( p + pz )/( p - pz ) );
                //            double eta = gRandom->Uniform(-1.2, 1.2); // 0.5 * TMath::Log( ( p + pz )/( p - pz ) );
//                double eta = gRandom->Uniform( -0.8, 0.8 ); // 0.5 * TMath::Log( ( p + pz )/( p - pz ) );
                double eta = gRandom->Uniform( -1.5, 1.5 ); // 0.5 * TMath::Log( ( p + pz )/( p - pz ) );
//                double eta = fSlopeEta->GetRandom();
                double phi = gRandom->Uniform( 0, TMath::TwoPi() );

                // !!!
                double pdg = gRandom->Uniform() > collectivity_degree ? 321 : 211;
                //                double pdg = ( i < nPions ) ? 211 : 321;
                //                double pdg = i < 0.2*nPart ? 321 : 211;
                if ( charge < 0 )
                    pdg *= -1;

                // ### "SRC":
                if(1)
                {
                    //            if ( i%2 == 0 )
                    if ( i%2 == 0 )
                    {
                        prevEta = eta;
                        prevPhi = phi;
                        prevPt = pt;
                        prevPdg = pdg;
                        prevCharge = charge;
                    }
                    //            else
                    //                else if ( abs(prevPdg) == 321 )
                    //                {
                    //                    eta = gRandom->Gaus(prevEta,0.2);
                    //                    pdg = -prevPdg;
                    //                    charge = -prevCharge;
                    //                    //                charge *= -1;
                    //                }
                    else
                    {
                        //                        double randWhat = gRandom->Uniform();
                        //                        if ( randWhat < 0.25 ) // opposite-sign correlation (25% prob.)
                        {
                            eta = gRandom->Gaus(prevEta,0.32);
                            phi = gRandom->Gaus(prevPhi,0.32);
                            if( phi < 0 )
                                phi += TMath::TwoPi();
                            else if( phi > TMath::TwoPi() )
                                phi -= TMath::TwoPi();
                            pt = gRandom->Gaus(prevPt, 0.02);
                            pdg = -prevPdg;
                            //charge = -prevCharge;
                            //                charge *= -1;
                        }
                        if(0)//else if ( randWhat < 0.8 && abs(prevPdg) == 321 ) // same-sign correlation (25% prob.)
                        {
                            double sepEta = 0.2 + fabs( gRandom->Gaus(0,0.4) );
                            while( fabs(eta-prevEta) < sepEta )
                            {
                                //                            cout << "check dEta = " << fabs(eta-prevEta) << endl;
                                eta = gRandom->Uniform( -2., 2. );
                            }
                            pdg = prevPdg;
                            charge = prevCharge;
                            //                charge *= -1;
                        }
                    }
                } // end of if SRC


                //            double pdg = gRandom->Uniform() > collectivity_degree ? 321 : 211;
                //            if ( charge < 0 )
                //                pdg *= -1;


                arr_charge[i] = charge;
                arr_pt[i] = pt;
                arr_eta[i] = eta;
                arr_phi[i] = phi;
                arr_pdg[i] = pdg;

                if( fabs(eta) < 0.05 )
                    counterNpartIn01Win++;

                //            if ( fabs(eta) < 1 )
                //                nParticlesForTrigger++;

                //            if( fabs( eta ) < 1.0 && ( pt > 0.15) ) //&& ( fnTPChits[iP] > 30 ) )
                //                nTracksForMultCentrality++;




                //            // ALTERNATIVE SRC:
                //            if( 1 && i<nPart-1 )
                //            {

                //                double randWhat = gRandom->Uniform();
                //                if ( randWhat < 0.25 ) // opposite-sign correlation (25% prob.)
                //                {
                //                    prevEta = eta;
                //                    prevPdg = pdg;
                //                    prevCharge = charge;

                //                    i++;

                //                    eta = gRandom->Gaus(prevEta,0.25);
                //                    pdg = -prevPdg;
                //                    charge = -prevCharge;
                //                    //                charge *= -1;
                //                }
                //                else if ( randWhat < 0.8 && abs(prevPdg) == 321 ) // same-sign correlation (25% prob.)
                //                {
                //                    prevEta = eta;
                //                    prevPdg = pdg;
                //                    prevCharge = charge;

                //                    i++;

                //                    double sepEta = 0.2 + fabs( gRandom->Gaus(0,0.4) );
                //                    while( fabs(eta-prevEta) < sepEta )
                //                    {
                ////                            cout << "check dEta = " << fabs(eta-prevEta) << endl;
                //                        eta = gRandom->Uniform( -2., 2. );
                //                    }
                //                    pdg = prevPdg;
                //                    charge = prevCharge;
                //                    //                charge *= -1;
                //                }



                //                int charge = gRandom->Uniform() > 0.5 ? +1 : -1;
                //                //            double pt = gRandom->Exp(0.5); // TMath::Sqrt(px*px + py*py);
                //                double pt = gRandom->Exp( meanPt ); // TMath::Sqrt(px*px + py*py);
                //                //            double pt = my_func_BW->GetRandom();
                //                //            double eta = gRandom->Uniform(-2, 2); // 0.5 * TMath::Log( ( p + pz )/( p - pz ) );
                //                //            double eta = gRandom->Uniform(-1.2, 1.2); // 0.5 * TMath::Log( ( p + pz )/( p - pz ) );
                //                double eta = gRandom->Uniform( -2., 2. ); // 0.5 * TMath::Log( ( p + pz )/( p - pz ) );
                //                double phi = gRandom->Uniform( 0, TMath::TwoPi() );

                //                arr_charge[i] = charge;
                //                arr_pt[i] = pt;
                //                arr_eta[i] = eta;
                //                arr_phi[i] = phi;
                //                arr_pdg[i] = pdg;
                //            }


            } // end of particle loop
        } // end of "old" model




        // "trigger"!
        //        if( nParticlesForTrigger == 0 )
        //        {
        //            h_QA_event->Fill(0);
        //            continue;
        //        }
        //        h_QA_event->Fill(1);


        // fill some histos:
        //        h_nTracksForMultCentrality->Fill( nTracksForMultCentrality );
        //        h_bImpTriggered->Fill(impact_b);
        h_nparticles_gen->Fill(nPart);
        h_nparticles_gen_in01win->Fill(counterNpartIn01Win);




        //        // determine centrality bin!
        //        int cBin = -1;
        //        for( int j = 0; j < nCbins; j++ )
        //            //            if ( bImp < bImpRanges[j] )
        //            //            if ( f_zdc_total_energy < energyRanges[j] )
        //            if ( nTracksForMultCentrality < multRanges[j] && nTracksForMultCentrality > multRanges[j+1] )
        //            {
        //                cBin = j;
        //                break;
        //            }
        //        if ( cBin >= nCbins ) // just outside of max b_imp thresh.
        //            continue;
        //        if ( cBin == -1 )
        //        {
        //            //            cout << "AHTUNG!!! cBin == -1" << endl;
        //            continue;
        //        }
        int cBinId[10];
        cBinId[0] = 0;//cBin; // tmp
        //        h_QA_centr_bins->Fill( cBin );


        //        continue;

        // assign subsample id
        int subsampleId = gRandom->Integer(nSubsamples);


        // ### main loop:
        int nSimTracksInCuts = 0;
        int nRecTracksInCuts = 0;
        for (int i = 0; i < /*nPart*//*counter*/nPart; i++)
        {
            //Double_t pt = TMath::Sqrt( px*px + py*py);
            //            cout << "pdgcode = " << pdgcode[itrack] << ", px,py,pz = " << px[itrack] << " " << py[itrack] << " " << pz[itrack] << endl;

            int charge = arr_charge[i];

            double pt = arr_pt[i];
            double eta = arr_eta[i];
            double phi = arr_phi[i];
            int pdg = arr_pdg[i];

            // VERY NAIVE corrections (Wrong):
            //            if ( pt < 0.3 )
            //                pt *= 1./0.5;
            //            else if ( pt > 0.3 && pt < 0.5 )
            //                pt *= 1./0.68;
            //            else if ( pt > 0.5 )
            //                pt *= 1./0.85;


            // FOR EFFICIENCY PLOTS:
            if( fabs(eta)<1.6 )
                h_Pt_in_eta16->Fill( pt );
            if( pt > 0.1 && pt < 10 )
                h_Eta_in_pt01_10->Fill( eta );

            //            if( fabs(eta)<1.0 && (pt > 0.15 ) ) // && pt < 10) )
            //                nSimTracksInCuts++;

            if( fabs(eta)>0.8 )
                continue;

            if( pt < 0.2 || pt > 2.0 )
                continue;

            nSimTracksInCuts++;

            h_Pt->Fill( pt );
            h_Eta->Fill( eta );
            //            h_Rapidity->Fill( rapidity );
            h_Phi->Fill( phi );

            h_pid->Fill( pdg );


            // add track to winPairs GEN
            for ( int iType = 0; iType < nPartTypes; iType++)
                for ( int iCW = 0; iCW < nCW; iCW++)
                    for ( int iPt = 0; iPt < nPtBins; ++iPt )
                    {
                        winPairWrapper[0][iType][iCW][cBinId[iCW]][iPt].addTrack( pdg, eta, phi, pt, charge, /*weight*/1.0 );
                        winPairWrapper_FULL_ETA_DENOM[0][iType][iCW][cBinId[iCW]][iPt].addTrack( pdg, eta, phi, pt, charge, /*weight*/1.0 );
                        winPairWrapper_FULL_ETA_NUM_AND_DENOM[0][iType][iCW][cBinId[iCW]][iPt].addTrack( pdg, eta, phi, pt, charge, /*weight*/1.0 );
                    }



            // artificial efficiency:
            //            if ( gRandom->Uniform() > eff )
            double randForEff = gRandom->Uniform();
            //            if (    (pt < 0.3 && randForEff > 0.5)
            //                 || ( (pt > 0.3 && pt < 0.5) && randForEff > 0.68 )
            //                 || ( pt > 0.5 && randForEff > 0.85 )
            //                 )
            if(1) // "toy" efficiency
            {
                if( fabs(pdg) == 211 ) // pi
                {
                    if (    (pt < 0.3 && randForEff > 0.43)
                            || ( (pt > 0.3 && pt < 0.5) && randForEff > 0.81 )
                            || ( (pt > 0.5 && pt < 0.7) && randForEff > 0.72 )
                            || ( pt > 0.7 && randForEff > 0.57 )
                            )
                        //                    if( randForEff > 0.65)
                        continue;
                }
                else // K or p
                {
                    if (    (pt < 0.3 && randForEff > 0.39)
                            || ( (pt > 0.3 && pt < 0.5) && randForEff > 0.45 )
                            || ( (pt > 0.5 && pt < 0.7) && randForEff > 0.58 )
                            || ( pt > 0.7 && randForEff > 0.32 )
                            )
                        continue;
                }
            }

            // realistic efficiency
            double eff = 1;
            if(0)
            {
                double eff = histEff->GetBinContent( histEff->FindBin( pt ) );
                //            cout << "pt = " << pt << ", eff = " << eff << ", bin " << histEff->FindBin( pt ) << endl;

                if ( gRandom->Uniform() > eff )
                    continue;
            }




            h_Pt_REC->Fill( pt );
            if ( abs(pdg) == 211)
                h_Pt_REC_pi->Fill( pt );
            else if ( abs(pdg) == 321)
                h_Pt_REC_K->Fill( pt );
            h_Eta_REC->Fill( eta );
            //            h_Rapidity->Fill( rapidity );
            h_Phi_REC->Fill( phi );
            h_pid_REC->Fill( pdg );

            nRecTracksInCuts++;


            // add track to winPairs REC
            for ( int iType = 0; iType < nPartTypes; iType++)
                for ( int iCW = 0; iCW < nCW; iCW++)
                    for ( int iPt = 0; iPt < nPtBins; ++iPt )
                    {
                        winPairWrapper[1][iType][iCW][cBinId[iCW]][iPt].addTrack( pdg, eta, phi, pt, charge, /*weight*/1.0 );
                        winPairWrapper_FULL_ETA_DENOM[1][iType][iCW][cBinId[iCW]][iPt].addTrack( pdg, eta, phi, pt, charge, /*weight*/1.0 );
                        winPairWrapper_FULL_ETA_NUM_AND_DENOM[1][iType][iCW][cBinId[iCW]][iPt].addTrack( pdg, eta, phi, pt, charge, /*weight*/1.0 );
                    }



            double weight = 1.;
            //            if ( pt < 0.3 )
            //                weight = 1./0.5;
            //            else if ( pt > 0.3 && pt < 0.5 )
            //                weight = 1./0.68;
            //            else if ( pt > 0.5 )
            //                weight = 1./0.85;
            if(1) // "toy" efficiency correction
            {
                if( fabs(pdg) == 211 ) // pi
                {
                    if ( pt < 0.3 )
                        weight = 1./0.43;
                    else if ( pt > 0.3 && pt < 0.5 )
                        weight = 1./0.81;
                    else if ( pt > 0.5 && pt < 0.7 )
                        weight = 1./0.72;
                    else if ( pt > 0.7 )
                        weight = 1./0.57;
                }
                else // K or p
                {
                    if ( pt < 0.3 )
                        weight = 1./0.39;
                    else if ( pt > 0.3 && pt < 0.5 )
                        weight = 1./0.45;
                    else if ( pt > 0.5 && pt < 0.7 )
                        weight = 1./0.58;
                    else if ( pt > 0.7 )
                        weight = 1./0.32;
                }
                //                weight = 1./0.65;
            }
            // realistic weight
            //            weight = 1./eff;


            h_weights->Fill(weight);

            //            TParticlePDG *partPDG = db_PDG->GetParticle( pdg );
            //            double m = partPDG->Mass();
            //            double E = TMath::Sqrt(m*m + px*px + py*py + pz*pz );
            //            double rapidity = 0.5 * TMath::Log( (E+pz) / (E-pz) );



            // add track to winPairs CORRECTED
            for ( int iType = 0; iType < nPartTypes; iType++)
                for ( int iCW = 0; iCW < nCW; iCW++)
                    for ( int iPt = 0; iPt < nPtBins; ++iPt )
                    {
                        winPairWrapper[2][iType][iCW][cBinId[iCW]][iPt].addTrack( pdg, eta, phi, pt, charge, weight );
                        winPairWrapper_FULL_ETA_DENOM[2][iType][iCW][cBinId[iCW]][iPt].addTrack( pdg, eta, phi, pt, charge, weight );
                        winPairWrapper_FULL_ETA_NUM_AND_DENOM[2][iType][iCW][cBinId[iCW]][iPt].addTrack( pdg, eta, phi, pt, charge, weight );
                    }


            //            primGen->AddTrack(pdgID, px, py, pz, 0., 0., 0.);

        } // end of main track loop

        h_nTracksInEtaCuts->Fill( nSimTracksInCuts );
        h_nTracksInEtaCuts_REC->Fill( nRecTracksInCuts );

        // finish event for winPairs:
        for ( int iAnLevel = 0; iAnLevel < N_AN_LEVELS; ++iAnLevel )
            for ( int iType = 0; iType < nPartTypes; iType++)
                for ( int iCW = 0; iCW < nCW; iCW++)
                    for ( int iPt = 0; iPt < nPtBins; ++iPt )
                    {
                        winPairWrapper[iAnLevel][iType][iCW][cBinId[iCW]][iPt].finishEvent(subsampleId);
                        winPairWrapper_FULL_ETA_DENOM[iAnLevel][iType][iCW][cBinId[iCW]][iPt].finishEvent(subsampleId);
                        winPairWrapper_FULL_ETA_NUM_AND_DENOM[iAnLevel][iType][iCW][cBinId[iCW]][iPt].finishEvent(subsampleId);
                    }


    } // end of event loop



    // ### End of file reading. Time:
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();

    printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);





    // canv_mean_pt_vs_T
    TCanvas *canv_mean_pt_vs_T = new TCanvas( "canv_mean_pt_vs_T", "canv_mean_pt_vs_T", 20,30,800,600 );
    canv_mean_pt_vs_T->Divide(4,2);
    canv_mean_pt_vs_T->cd(1);
    hist_fluctT_mean_pt_pions->DrawCopy();
    //    canv_mean_pt_vs_T->cd(2);
    //    hist_fluctT_mean_pt_kaons->DrawCopy();
    //    canv_mean_pt_vs_T->cd(3);
    //    hist_fluctT_mean_pt_protons->DrawCopy();
    canv_mean_pt_vs_T->cd(4);
    hist_fluctT_mean_pt_pions_vs_Tkin->DrawCopy("colz");
    canv_mean_pt_vs_T->cd(5);
    hist_fluctT_mean_pt_pions_vs_betaT->DrawCopy("colz");



    TCanvas *canv_QA = new TCanvas("canv_QA","canv_QA",20,20,1200,600 );
    canv_QA->Divide(2,2);

    //    canv_QA->cd(1);
    //    h_bImp->Draw("same");
    canv_QA->cd(1);
    h_nparticles_gen->DrawCopy();
    h_nparticles_gen_in01win->SetLineColor(kRed);
    h_nparticles_gen_in01win->Draw("same");
    //
    //    canv_QA->cd(3);
    //    h_px->Draw("same");
    //    canv_QA->cd(4);
    //    h_py->Draw("same");
    //    canv_QA->cd(5);
    //    h_pz->Draw("same");
    canv_QA->cd(2);
    h_pid->Draw("same");



    // ### write output
    //    TString strOutputDir = "output";
    TString strOutputDir = ".";
    TString strOutFileName = "output_toy.root";
    //    TString strOutFileName = "toy_output_gaus_500_50.root";
    //    TString strOutFileName = "toy_output_no_collectivity.root";
    //    TString strOutFileName = "toy_output_collectivity_prop_to_Nch.root";
    TString fullOutputName = Form("%s/%s", strOutputDir.Data(), strOutFileName.Data() );
    cout << "writing the output to " << fullOutputName << "..." << endl;
    TFile *fileOutput = new TFile( fullOutputName, "recreate");

    h_QA_event->Write();

    h_Pt_in_eta16->Write();
    h_Eta_in_pt01_10->Write();


    //    h_bImp->Write();
    //    h_bImpTriggered->Write();
    h_nparticles_gen->Write();
    h_Pt->Write();
    h_Eta->Write();
    h_Phi->Write();
    h_pid->Write();
    //    h_Rapidity->Write();

    h_Pt_REC->Write();
    h_Pt_REC_pi->Write();
    h_Pt_REC_K->Write();
    h_Eta_REC->Write();
    h_Phi_REC->Write();
    h_pid_REC->Write();

    h_weights->Write();

    h_nTracksForMultCentrality->Write();
    h_nTracksInEtaCuts->Write();
    h_nTracksInEtaCuts_REC->Write();
    h_QA_centr_bins->Write();
    h_QA_collectivity_degree->Write();

    // ### WRITE FB:
    if(1)
        for ( int iType = 0; iType < nPartTypes; iType++)
        {
            fOutputWinPairsLists[iType]->Write( fOutputWinPairsLists[iType]->GetName(), TObject::kSingleKey );
        }

    //histEff->Write();

    fileOutput->Close();


}
