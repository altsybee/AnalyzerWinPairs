#ifndef SimpleCalculations_cxx
#define SimpleCalculations_cxx


#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THn.h"
#include "TDirectory.h"
#include "TString.h"
#include "TGraphErrors.h"
//#include "TGraph2D.h"
#include "TGraph2DErrors.h"

#include <iostream>
#include <string>
using namespace std;



//int nEtaWins = 8;
//int nPhiWins = 1;

bool DO_GLOBAL_QA = false;//true;
bool FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_1D = true;
bool FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_2D = true;
bool FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_3D = true;
bool FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_3D_bis = true;
bool FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_Obs = true;
bool FLAG_QA_TRIPLETS = true;


// names of final quantities
const char *obsNames[] =
{
    //    "bcorr_FB",
    //    "bcorr_PtN",
    //    "bcorr_PtPt",

    "nu_dyn_FB",
    //    "nu_dyn_PtN",
    //    "nu_dyn_PtPt",

    "sigma_FB",
    //    "sigma_PtN",
    //    "sigma_PtPt",


    "bcorr_XY",
    "nu_dyn_XY",
    "sigma_XY",

    "nu_dyn_FY",
    "sigma_FY",

    "nu_dyn_XB",
    "sigma_XB",


    "avF",
    "avB",
    "avX",
    "avY",

    "FB",
    "XY",
    "FY",
    "XB",

    "avPtF",
    "avPtB",
    "avPtX",
    "avPtY",

    "nB*PX",
    "nY*PX",

    "PF*PB",
    "nF*PB",
    "nB*PF",



    "corr_rr_formula",
    "corr_rr_direct",

    "corr_R2_aa",
    "corr_R2_bb",
    "corr_R2_ab",
    "corr_R2_ba",


    // ##### meanPt-meanPt VS ETA (July 2022):
    "avPtF_avPtB_direct",
    "avPtF_avPtB_formula",

    "avPtX_avPtY_direct",
    "avPtX_avPtY_formula",


    // new ratio-ratio - June 2021: when denominator is in full acceptance
    "corr_rr_FULL_ETA_DENOM_formula",
//    "corr_rr_FULL_ETA_DENOM_formula_with_minusC2X",
    "corr_rr_FULL_ETA_DENOM_direct",

    "corr_rPt_formula",   // new ratio-meanPt - July 2022: when expansion for <pT> is done as <nB*PX>/.. + <nY*nX>/.. - <nB*nX>/.. - <nY*PX>/..
    "corr_rPt_direct",



    // ### March 2023: C3
//    "c3",



    // new ratio-pt - April 2021: when pt is not averaged in each event
//    "corr_rSumPt_formula",
//    "corr_rSumPt_direct",


    //
//    "coeff_dptdpt",
//    "coeff_ptpt_CHECK_VV_formula",
//    "coeff_ptpt_BOZEK_2017",
//    "C_BOZEK_F",
//    "C_BOZEK_B",
};
const int nObs = sizeof(obsNames)/sizeof(*obsNames); // number of observables



// names of final quantities
const char *obsTripletNames[] =
{
    // ### March 2023: C3
    "c3",
    "c3_direct",
    "c3_direct_Poisson",
};
const int nTripletObs = sizeof(obsTripletNames)/sizeof(*obsTripletNames); // number of observables





// ################################################
// check if particles are identical: F==B, X==Y
bool ifIdenticalParticlesFB( int *_pTypes, int *_pCharges )
{
    if( _pTypes[0] != _pTypes[1] ) return false;
    if( _pTypes[2] != _pTypes[3] ) return false;
    if( _pCharges[0] != _pCharges[1] ) return false;
    if( _pCharges[2] != _pCharges[3] ) return false;

//    if( h3D_winTripletsInfo )   // March 15: something should be done here for triplets! more complex combinatorics...


    return true;
}





// ################################################
//const int MAX_N_WIN_PAIRS = 500;//150;
struct CalcWithSubsamples
{
    int whichInputHist;

    TH3D *h3D_singleWinInfo; // incoming hist with variables
    TH3D *h3D_winPairInfo; // incoming hist with variables
    TH3D *h3D_winTripletsInfo; //
    map<string, int> mapVarIdByNameSingleWin;
    map<string, int> mapVarIdByNameWinPairs;
    map<string, int> mapVarIdByNameWinTriplets;

    // current status vars:
    int currentWinPairId;
    int currentDetaDphiPairId;
    int currentEtaWinsDphiPairId;

    int currentWinTripletId;


    int currentSubId;
    double currenEtaSep;
    double currenPhiSep;


    double currenEtaSepFB;
    double currenEtaSepFC;
    double currenEtaSepBC;


    int nEtaWins;
    int nPhiWins;
    double eRangeMin;
    double eRangeMax;

    //    int current_eWinId;
    //    int current_pWinId;
    int current_winIdF;
    int current_winIdB;
    int current_winIdC;

    int current_wp_for_triplets_FB;
    int current_wp_for_triplets_FC;
    int current_wp_for_triplets_BC;

    int current_wp_for_triplets_BF;
    int current_wp_for_triplets_CF;
    int current_wp_for_triplets_CB;


    TGraphErrors *grVsWinId[nObs];
    TGraphErrors *grVsDeltaEta[nObs];
    TGraph2DErrors *gr2D_dEta_dPhi[nObs];
    TH2D *h2D_dEta_dPhi[nObs];
    TH2D *h2D_dEta_dPhi_binError2[nObs];
    TH2D *h2D_dEta_dPhi_binCounts[nObs];
    //    TH2D *hist2D_dEta_dPhi[nObs];
    double integrals[nObs];

    TGraphErrors *grVsTripletId[nObs];
    TGraphErrors *grTriplets_TwoGap0_vs_Third[nObs];
    TGraphErrors *grTriplets_TwoGap1_vs_Third[nObs];
    TGraphErrors *grTriplets_TwoGap2_vs_Third[nObs];

    //    int nVars;
    int nWinPairs;
    int nWinTriplets;

    int nSubs;

    TH3D *histCalcObs;   // hist with calculated observables
    map<string, int> mapObsIdByName;   // observable name-binId correspondance


    TH3D *histCalcTripletObs;   // hist with calculated observables
    map<string, int> mapTripletObsIdByName;   // observable name-binId correspondance



    double eSizeNum;
    double eSizeDenom;
    bool if_Identical_FB_XY;
    double pSizeNum;

    bool QA_FLAG;

    // possibility to set bins zero by hand (artificial acceptance map) - Sept 12, 2022
    vector<int> arrZeroWinIdsByHand;

    CalcWithSubsamples()
    {
        h3D_singleWinInfo = 0x0;
        h3D_winPairInfo = 0x0;
        h3D_winTripletsInfo = 0x0;

        QA_FLAG = false;
//        nArrSizeZeroWinsByHand = 0;
    }

    TGraphErrors *getGraph( const char *strName )
    {
        int idObs = mapObsIdByName[strName];
        return grVsDeltaEta[idObs];
    }

    TGraphErrors *getGraphVsWinIdQA( const char *strName )
    {
        int idObs = mapObsIdByName[strName];
        return grVsWinId[idObs];
    }

    TGraphErrors *getGraphVsTripletIdQA( const char *strName )
    {
        int idObs = mapTripletObsIdByName[strName];
        return grVsTripletId[idObs];
    }

    TGraphErrors *getGraphTriplets_TwoGap0_vs_Third( const char *strName )
    {
        int idObs = mapTripletObsIdByName[strName];
        return grTriplets_TwoGap0_vs_Third[idObs];
    }


    TGraph2D *getGraph2D( const char *strName )
    {
        int idObs = mapObsIdByName[strName];
        return gr2D_dEta_dPhi[idObs];
    }

    TH2D *getHist2D( const char *strName )
    {
        int idObs = mapObsIdByName[strName];
        return h2D_dEta_dPhi[idObs];
    }


    double getIntegral( const char *strName )
    {
        int idObs = mapObsIdByName[strName];
        return integrals[idObs];
    }


    // ###########
    void takeEtaSep( int _iWinPair )
    {
        double _etaSep = -1000;
        double _phiSep = -1000;

        TString strEtaBin;
        if ( h3D_winPairInfo )
            strEtaBin = h3D_winPairInfo->GetYaxis()->GetBinLabel( _iWinPair + 1 );

        if(0) cout << "strEtaBin: " << strEtaBin << endl;

        if ( whichInputHist == 0 )   //strEtaBin.Contains("etaB") ) // i.e. we do pair-by-pair hist3D
        {
            float eBounds[4];
            float phiBounds[4];
            //            sscanf( strEtaBin.Data(), "eB_%f_%f_eF_%f_%f", &eBounds[0], &eBounds[1], &eBounds[2], &eBounds[3] );
            sscanf( strEtaBin.Data(), "etaB_%f_%f_phiB_%f_%f_etaF_%f_%f_phiF_%f_%f"
                    , &eBounds[0], &eBounds[1], &phiBounds[0], &phiBounds[1]
                    , &eBounds[2], &eBounds[3], &phiBounds[2], &phiBounds[3] );


            if(0)for ( int k = 0; k < 4; k++ )
            {
                //                    eBounds[k] = round( eBounds[k]*100 ) / 100;
                cout << "eBounds[" << k << "] = " << eBounds[k]  << endl;
            }
            //            eSize = eBounds[1] - eBounds[0];
            double eBpos = ( eBounds[0] + eBounds[1] )/2;
            double eFpos = ( eBounds[2] + eBounds[3] )/2;
            _etaSep = round( (eFpos - eBpos)*100 ) / 100;  // ="1 - 2"

            double phiBpos = ( phiBounds[0] + phiBounds[1] )/2;
            double phiFpos = ( phiBounds[2] + phiBounds[3] )/2;
            _phiSep = round( (phiFpos - phiBpos)*100 ) / 100;  // ="1 - 2"



            // single wins:
            TString strBinF = Form( "eta_%.2f_%.2f_phi_%.2f_%.2f", eBounds[2], eBounds[3], phiBounds[2], phiBounds[3] );
            TString strBinB = Form( "eta_%.2f_%.2f_phi_%.2f_%.2f", eBounds[0], eBounds[1], phiBounds[0], phiBounds[1] );
            //            cout << "strBinF = " << strBinF << ", strBinB = " << strBinB << endl;
            current_winIdF = h3D_singleWinInfo->GetYaxis()->FindBin( strBinF.Data() ) - 1;
            current_winIdB = h3D_singleWinInfo->GetYaxis()->FindBin( strBinB.Data() ) - 1;

        }
        else if ( whichInputHist == 1 )   // if dEta-dPhi
        {
            float dEta, dPhi;
            sscanf( strEtaBin.Data(), "dEta_%f_dPhi_%f", &dEta, &dPhi );
            //            dEta = h3D->GetYaxis()->GetBinCenter();
            _etaSep = round( dEta *100 ) / 100;
            _phiSep = round( dPhi *100 ) / 100;
            if(0)
                cout << "_etaSep = " << _etaSep << ", _phiSep = " << _phiSep << endl;
        }
        else if ( whichInputHist == 2 )   // hist AllEta-Dphi
        {
            float eBounds[4];
            //            float phiBounds[4];
            float dPhi;
            sscanf( strEtaBin.Data(), "etaB_%f_%f_etaF_%f_%f_dPhi_%f", &eBounds[0], &eBounds[1], &eBounds[2], &eBounds[3], &dPhi );


            if(0)for ( int k = 0; k < 4; k++ )
            {
                //                    eBounds[k] = round( eBounds[k]*100 ) / 100;
                cout << "eBounds[" << k << "] = " << eBounds[k]  << endl;
            }
            //            eSize = eBounds[1] - eBounds[0];
            double eBpos = ( eBounds[0] + eBounds[1] )/2;
            double eFpos = ( eBounds[2] + eBounds[3] )/2;
            _etaSep = round( (eFpos - eBpos)*100 ) / 100;  // ="1 - 2"
            _phiSep = dPhi;
        }

        if ( _etaSep < -999 ) { cout << "AHTUNG!!!  _etaSep < -999!" << endl; int aa; cin>>aa; }
        if ( _phiSep < -999 ) { cout << "AHTUNG!!!  _phiSep < -999!" << endl; int aa; cin>>aa; }

        currenEtaSep = _etaSep;
        currenPhiSep = _phiSep;
    }





    // ###########
    void takeWinsForTriplets( int _iTriplet )
    {
//        double _etaSep = -1000;
//        double _phiSep = -1000;

        TString strEtaBin;
        if ( h3D_winTripletsInfo )
            strEtaBin = h3D_winTripletsInfo->GetYaxis()->GetBinLabel( _iTriplet + 1 );

        if(0) cout << "strEtaBin: " << strEtaBin << endl;

//        if ( whichInputHist == 0 )   //strEtaBin.Contains("etaB") ) // i.e. we do pair-by-pair hist3D

            float eBounds[6];
            float phiBounds[6];
            //            sscanf( strEtaBin.Data(), "eB_%f_%f_eF_%f_%f", &eBounds[0], &eBounds[1], &eBounds[2], &eBounds[3] );
            sscanf( strEtaBin.Data(), "etaB_%f_%f_phiB_%f_%f_etaF_%f_%f_phiF_%f_%f_etaC_%f_%f_phiC_%f_%f"
                    , &eBounds[0], &eBounds[1], &phiBounds[0], &phiBounds[1]
                    , &eBounds[2], &eBounds[3], &phiBounds[2], &phiBounds[3]
                    , &eBounds[4], &eBounds[5], &phiBounds[4], &phiBounds[5]
                    );


            if(0)for ( int k = 0; k < 4; k++ )
            {
                //                    eBounds[k] = round( eBounds[k]*100 ) / 100;
                cout << "eBounds[" << k << "] = " << eBounds[k]  << endl;
            }
            //            eSize = eBounds[1] - eBounds[0];
            double eBpos = ( eBounds[0] + eBounds[1] )/2;
            double eFpos = ( eBounds[2] + eBounds[3] )/2;
            double eCpos = ( eBounds[5] + eBounds[4] )/2;
            currenEtaSepFB = round( (eFpos - eBpos)*100 ) / 100;  // ="1 - 2"
            currenEtaSepFC = round( (eFpos - eCpos)*100 ) / 100;  // ="1 - 2"
            currenEtaSepBC = round( (eBpos - eCpos)*100 ) / 100;  // ="1 - 2"


//            double phiBpos = ( phiBounds[0] + phiBounds[1] )/2;
//            double phiFpos = ( phiBounds[2] + phiBounds[3] )/2;
//            _phiSep = round( (phiFpos - phiBpos)*100 ) / 100;  // ="1 - 2"



            // single wins:
            TString strBinF = Form( "eta_%.2f_%.2f_phi_%.2f_%.2f", eBounds[2], eBounds[3], phiBounds[2], phiBounds[3] );
            TString strBinB = Form( "eta_%.2f_%.2f_phi_%.2f_%.2f", eBounds[0], eBounds[1], phiBounds[0], phiBounds[1] );
            TString strBinC = Form( "eta_%.2f_%.2f_phi_%.2f_%.2f", eBounds[4], eBounds[5], phiBounds[4], phiBounds[5] );
            //            cout << "strBinF = " << strBinF << ", strBinB = " << strBinB << endl;
            current_winIdF = h3D_singleWinInfo->GetYaxis()->FindBin( strBinF.Data() ) - 1;
            current_winIdB = h3D_singleWinInfo->GetYaxis()->FindBin( strBinB.Data() ) - 1;
            current_winIdC = h3D_singleWinInfo->GetYaxis()->FindBin( strBinC.Data() ) - 1;



            //
            TString strBinFB = Form( "etaB_%.2f_%.2f_phiB_%.2f_%.2f_etaF_%.2f_%.2f_phiF_%.2f_%.2f"
                                     , eBounds[0], eBounds[1], phiBounds[0], phiBounds[1]
                                     , eBounds[2], eBounds[3], phiBounds[2], phiBounds[3]
                    );
            TString strBinFC = Form( "etaB_%.2f_%.2f_phiB_%.2f_%.2f_etaF_%.2f_%.2f_phiF_%.2f_%.2f"
                                     , eBounds[4], eBounds[5], phiBounds[4], phiBounds[5]
                                     , eBounds[2], eBounds[3], phiBounds[2], phiBounds[3]
                    );
            TString strBinBC = Form( "etaB_%.2f_%.2f_phiB_%.2f_%.2f_etaF_%.2f_%.2f_phiF_%.2f_%.2f"
                                     , eBounds[4], eBounds[5], phiBounds[4], phiBounds[5]
                                     , eBounds[0], eBounds[1], phiBounds[0], phiBounds[1]
                    );
            //
            TString strBinBF = Form( "etaB_%.2f_%.2f_phiB_%.2f_%.2f_etaF_%.2f_%.2f_phiF_%.2f_%.2f"
                                     , eBounds[2], eBounds[3], phiBounds[2], phiBounds[3]
                                     , eBounds[0], eBounds[1], phiBounds[0], phiBounds[1]
                    );
            TString strBinCF = Form( "etaB_%.2f_%.2f_phiB_%.2f_%.2f_etaF_%.2f_%.2f_phiF_%.2f_%.2f"
                                     , eBounds[2], eBounds[3], phiBounds[2], phiBounds[3]
                                     , eBounds[4], eBounds[5], phiBounds[4], phiBounds[5]
                    );
            TString strBinCB = Form( "etaB_%.2f_%.2f_phiB_%.2f_%.2f_etaF_%.2f_%.2f_phiF_%.2f_%.2f"
                                     , eBounds[0], eBounds[1], phiBounds[0], phiBounds[1]
                                     , eBounds[4], eBounds[5], phiBounds[4], phiBounds[5]
                    );

            current_wp_for_triplets_FB = h3D_winPairInfo->GetYaxis()->FindBin( strBinFB.Data() ) - 1;
            current_wp_for_triplets_FC = h3D_winPairInfo->GetYaxis()->FindBin( strBinFC.Data() ) - 1;
            current_wp_for_triplets_BC = h3D_winPairInfo->GetYaxis()->FindBin( strBinBC.Data() ) - 1;

            current_wp_for_triplets_BF = h3D_winPairInfo->GetYaxis()->FindBin( strBinBF.Data() ) - 1;
            current_wp_for_triplets_CF = h3D_winPairInfo->GetYaxis()->FindBin( strBinCF.Data() ) - 1;
            current_wp_for_triplets_CB = h3D_winPairInfo->GetYaxis()->FindBin( strBinCB.Data() ) - 1;

//        if ( _etaSep < -999 ) { cout << "AHTUNG!!!  _etaSep < -999!" << endl; int aa; cin>>aa; }
//        if ( _phiSep < -999 ) { cout << "AHTUNG!!!  _phiSep < -999!" << endl; int aa; cin>>aa; }

//        currenEtaSep = _etaSep;
//        currenPhiSep = _phiSep;

    }






    // #####################
    void calc( TList *inputList, /*const char *strAnLevel, int cBin,*/ TString strPostfix, bool _if_Identical_FB_XY, int _whichInputHist, bool _calcTriplets = false ) //  double _eSizeNum, double _eSizeDenom, bool _if_Identical_FB_XY, double _pSize )
    {
//        cout << "_calcTriplets = " << _calcTriplets << endl;

        // extract meta info:
        TH1D *hMetaInfo = (TH1D*)inputList->FindObject( Form("hMetaInfo_%s", strPostfix.Data() ) );
//        cout << "hMetaInfo = " << hMetaInfo << endl;

        int nInstances = hMetaInfo->GetBinContent( hMetaInfo->GetXaxis()->FindBin( "nInstances" ) );
        if ( nInstances == 0 ) // TMP!!! for compatibility with task_2022_08_01_Data_PbPb_LHC10h_FB768_NewWP_e8p1_all_histos
            nInstances = 8837;
        eRangeMin = hMetaInfo->GetBinContent( hMetaInfo->GetXaxis()->FindBin( "etaRangeMin" ) ) / nInstances;
        eRangeMax = hMetaInfo->GetBinContent( hMetaInfo->GetXaxis()->FindBin( "etaRangeMax" ) ) / nInstances;
        eSizeNum = hMetaInfo->GetBinContent( hMetaInfo->GetXaxis()->FindBin( "etaSize" ) ) / nInstances;
        eSizeDenom = eSizeNum;
        if ( hMetaInfo->GetBinContent( hMetaInfo->GetXaxis()->FindBin( "fullAcceptanceForDenom" ) ) )
            eSizeDenom = eRangeMax - eRangeMin;
        double _pSize = hMetaInfo->GetBinContent( hMetaInfo->GetXaxis()->FindBin( "phiSize" ) ) / nInstances;


        nEtaWins = hMetaInfo->GetBinContent( hMetaInfo->GetXaxis()->FindBin( "nEtaWins" ) ) / nInstances;
        nPhiWins = hMetaInfo->GetBinContent( hMetaInfo->GetXaxis()->FindBin( "nPhiWins" ) ) / nInstances;

//        cout << " eRangeMax = " << eRangeMax << endl;
//        cout << " eSizeNum = " << eSizeNum << endl;
//        cout << " nEntries = " << eSizeNum / 0.2 << endl;
//        cout << " hMetaInfo->GetBinContent( hMetaInfo->GetXaxis()->FindBin( fullAcceptanceForDenom ) ) = " << hMetaInfo->GetBinContent( hMetaInfo->GetXaxis()->FindBin( "fullAcceptanceForDenom" ) ) << endl;
//        int aa;
//        cin >> aa;

        // create graphs per each observable
        for( int iObs = 0; iObs < nObs; iObs++ )
        {
            grVsDeltaEta[iObs] = new TGraphErrors;
            grVsWinId[iObs] = new TGraphErrors;
            integrals[iObs] = 0;

            gr2D_dEta_dPhi[iObs] = new TGraph2DErrors;
        }

        for( int iObs = 0; iObs < nTripletObs; iObs++ )
        {
            grVsTripletId[iObs] = new TGraphErrors;
            grTriplets_TwoGap0_vs_Third[iObs] = new TGraphErrors;
            grTriplets_TwoGap1_vs_Third[iObs] = new TGraphErrors;
            grTriplets_TwoGap2_vs_Third[iObs] = new TGraphErrors;
        }



        int nVars1D = -1;
        h3D_singleWinInfo = (TH3D*) inputList->FindObject( Form("hSingleWin_%s", strPostfix.Data() ) );
        if( h3D_singleWinInfo )
            nVars1D = h3D_singleWinInfo->GetNbinsX();

        // artificial making of zero bins (Sept 2022):
        if( arrZeroWinIdsByHand.size() > 0 )
        {
            for( int i = 0; i < arrZeroWinIdsByHand.size(); i++ )
            {
                int winIdToMakeZero = arrZeroWinIdsByHand[i];
                for( int j = 0; j < h3D_singleWinInfo->GetNbinsX(); j++ ) // nVar
                    for( int k = 0; k < h3D_singleWinInfo->GetNbinsZ(); k++ ) // nSub
                        h3D_singleWinInfo->SetBinContent( j+1, winIdToMakeZero+1, k+1, 0 );
            }
        }

        whichInputHist = _whichInputHist;
        if( whichInputHist == 0 )
            h3D_winPairInfo = (TH3D*) inputList->FindObject( Form("hAllWins_%s", strPostfix.Data() ) );
        else if ( whichInputHist == 1 )
            h3D_winPairInfo = (TH3D*) inputList->FindObject( Form("hDetaDphi_%s", strPostfix.Data() ) );
        else if ( whichInputHist == 2 )
            h3D_winPairInfo = (TH3D*) inputList->FindObject( Form("hAllEtaDphi_%s", strPostfix.Data() ) );


        if( _calcTriplets )
        {
            cout << "### Calculating triplets!" << endl;
            h3D_winTripletsInfo = (TH3D*) inputList->FindObject( Form("hAllWinsTriplets_%s", strPostfix.Data() ) );
        }


        // calc integral for F, B, X, Y
//        {
//            double value = 0;
//            int varId = mapVarIdByNameSingleWin.at( "Nf" );
//            for( int winId = 0; winId < h3D_singleWinInfo->GetNbinsY(); winId++ )
//                for( int subId = 0; subId < h3D_singleWinInfo->GetNbinsZ(); subId++ )
//                    value += h3D_singleWinInfo->GetBinContent( varId+1, winId+1, subId+1 );
//        }


        if_Identical_FB_XY = _if_Identical_FB_XY;
        pSizeNum = _pSize;

        int nVars2D = -1;
        if( h3D_winPairInfo )
        {
            nVars2D = h3D_winPairInfo->GetNbinsX();
            nWinPairs = h3D_winPairInfo->GetNbinsY();
            nSubs = h3D_winPairInfo->GetNbinsZ();
        }




        cout << "##### nObs = " << nObs << ", nVars1D = " << nVars1D << ", nVars2D = " << nVars2D << ", nSubs = " << nSubs << ", nWinPairs = " << nWinPairs << ", eSizeNum = " << eSizeNum << ", eSizeDenom = " << eSizeDenom << ", nSubs = " << nSubs << endl;

        createBinNameBinIdMap();

        // !!! check for empty bins:
        if(0)for( int i = 0; i < nVars2D; i++ )
            for( int j = 0; j < nWinPairs; j++ )
                for( int k = 0; k < nSubs; k++ )
                {
                    double binContent = h3D_winPairInfo->GetBinContent(i+1, j+1, k+1);
                    if ( binContent < 0.00001 )
                    {
                        cout << "AHTUNG!!! binContent < 0.00001 !!" << " i=" << i << " j=" << j << " k=" << k <<
                                "bin label: " << h3D_winPairInfo->GetXaxis()->GetBinLabel(i+1) << endl;
                    }

                }


        // calculated observables
        TString strHistObsTitle = Form( "%s_observables", h3D_winPairInfo->GetName() );
        histCalcObs = new TH3D( strHistObsTitle, strHistObsTitle  //strAccumHistName , strAccumHistName //"histAccumulatedValues"
                                , nObs,-0.5, nObs-0.5
                                , nWinPairs, -0.5, nWinPairs-0.5
                                , nSubs+1, -0.5, nSubs-0.5 + 1  // +1 more for full calc (=sum of subsamples)
                                );

        // set x labels for final hist with observables
        for( int i=0; i < nObs; i++ )
        {
            mapObsIdByName.insert( pair<string, int>( obsNames[i], i ) );
            histCalcObs->GetXaxis()->SetBinLabel( i+1, obsNames[i] );
        }

        // loop over win pairs:
        for ( int iWinPair = 0; iWinPair < nWinPairs; iWinPair++ )
        {
//            cout << "###### STARTING iWinPair = " << iWinPair << endl;

            // find out eSep, phiSep:
            takeEtaSep( iWinPair );
            //            cout << "eSep = " << currenEtaSep << ", currenPhiSep = " << currenPhiSep << endl;

            //            cout << "current_winIdF = " << current_winIdF << endl;
            //            cout << "current_winIdB = " << current_winIdB << endl;

            //            int aa;
            //            cin >> aa;


            //            current_eWinId = iWinPair / nEtaWins;
            //            current_pWinId = iWinPair % nEtaWins;


            currentWinPairId = iWinPair;
            currentSubId = nSubs; // i.e. calc for sum of all subsamples

            if( DO_GLOBAL_QA )
                cout << "##### WIN PAIR " << h3D_winPairInfo->GetYaxis()->GetBinLabel( currentWinPairId+1 ) << endl;
//            if ( valueByNameWP( "Nf*Nb" ) < 0.00001 )
//                continue;
            bool isGoodCalc = finalCalcPerWP(); // in case of holes in the acceptance - win pair will be skipped
            if ( !isGoodCalc )
                continue;

            // subsamples
            for ( int iSub = 0; iSub < nSubs; iSub++ )
            {
                // cout << "iSub=" << iSub << ", iType=" << iType << ", iCW=" << iCW << ", iWinPair=" << iWinPair << endl;
                currentSubId = iSub;
                finalCalcPerWP();
            }

            // calc mean and errors:
            for( int iObs = 0; iObs < nObs; iObs++ )
            {
                // mean
                double mean = histCalcObs->GetBinContent( iObs+1, iWinPair+1, nSubs+1 );

                // stdDev:
                double std_dev = 0;
                for ( int iSub = 0; iSub < nSubs; iSub++)
                {
                    double value = histCalcObs->GetBinContent( iObs+1, iWinPair+1, iSub+1 );
                    float diff = value - mean;
                    std_dev += diff*diff;
                }
                std_dev = sqrt(std_dev / nSubs / (nSubs-1) );

                //                cout << "iObs = " << iObs << ", obs = " << obsNames[iObs] << ": mean +/- std_dev: " << mean << "   " << std_dev << endl;

                histCalcObs->SetBinError( iObs+1, iWinPair+1, nSubs+1, std_dev );

                //                cout << wpObs[iWinPair]->histCalcObs->GetXaxis()->GetBinLabel( bin+1 ) << " = " << mean << ", std_dev = " << std_dev << endl;

                // add points to graphs:
                int nGrP = grVsDeltaEta[iObs]->GetN();
                grVsDeltaEta[iObs]->SetPoint( nGrP, currenEtaSep, mean );
                grVsDeltaEta[iObs]->SetPointError( nGrP, 0, std_dev );

                grVsWinId[iObs]->SetPoint( nGrP, iWinPair, mean );
                grVsWinId[iObs]->SetPointError( nGrP, 0, std_dev );



                gr2D_dEta_dPhi[iObs]->SetPoint( nGrP, currenEtaSep, currenPhiSep, mean );
                gr2D_dEta_dPhi[iObs]->SetPointError( nGrP, 0, 0, std_dev );//currenEtaSep, currenPhiSep, std_dev );
                //                cout << " gr2D_dEta_dPhi[iObs]->GetN() = " << gr2D_dEta_dPhi[iObs]->GetN() << endl;
                //                cout << "iWinPair = " << iWinPair << ", obs = " << obsNames[iObs] << ", currenEtaSep = " << currenEtaSep << ", currenPhiSep = " << currenPhiSep << ", mean = " << mean << endl;

                integrals[iObs] += mean;
            }

//            int aa;
//            cin >> aa;
        }  // end of loop over win pairs




        // ##### calculate TRIPLET observables (March 2023)
        if( h3D_winTripletsInfo )
        {
            cout << "### Calculating triplets - test point here" << endl;


            int nVars3D = h3D_winTripletsInfo->GetNbinsX();
            nWinTriplets = h3D_winTripletsInfo->GetNbinsY();

            cout << "##### nTripletObs = " << nTripletObs << ", nVars3D = " << nVars3D << ", nWinTriplets = " << nWinTriplets << endl;


            TString strHistTripletObsTitle = Form( "%s_observables", h3D_winTripletsInfo->GetName() );
            histCalcTripletObs = new TH3D( strHistTripletObsTitle, strHistTripletObsTitle  //strAccumHistName , strAccumHistName //"histAccumulatedValues"
                                           , nTripletObs,-0.5, nTripletObs-0.5
                                           , nWinTriplets, -0.5, nWinTriplets-0.5
                                           , nSubs+1, -0.5, nSubs-0.5 + 1  // +1 more for full calc (=sum of subsamples)
                                           );

            // set x labels for final hist with observables
            for( int i=0; i < nTripletObs; i++ )
            {
                mapTripletObsIdByName.insert( pair<string, int>( obsTripletNames[i], i ) );
                histCalcTripletObs->GetXaxis()->SetBinLabel( i+1, obsTripletNames[i] );
            }
            cout << "##### nTripletObs = " << nTripletObs << ", nVars1D = " << nVars1D << ", nVars2D = " << nVars2D << ", nSubs = " << nSubs << ", nWinPairs = " << nWinPairs << ", eSizeNum = " << eSizeNum << ", eSizeDenom = " << eSizeDenom << ", nSubs = " << nSubs << endl;


            // loop over TROPLETS:
            for ( int iTriplet = 0; iTriplet < nWinTriplets; iTriplet++ )
            {
                //            cout << "###### STARTING iTriplet = " << iTriplet << endl;

                // find out eSep, phiSep:
                takeWinsForTriplets( iTriplet );
                //            cout << "eSep = " << currenEtaSep << ", currenPhiSep = " << currenPhiSep << endl;

                //            cout << "current_winIdF = " << current_winIdF << endl;
                //            cout << "current_winIdB = " << current_winIdB << endl;

                //            int aa;
                //            cin >> aa;


                //            current_eWinId = iWinPair / nEtaWins;
                //            current_pWinId = iWinPair % nEtaWins;
//                cout << "triplet id = " << iTriplet << ", win ids: F = " << current_winIdF << ", B = " << current_winIdB << ", C = " << current_winIdC << endl;
//                int aa;
//                cin >> aa;


                currentWinTripletId = iTriplet;
                currentSubId = nSubs; // i.e. calc for sum of all subsamples

                if( DO_GLOBAL_QA )
                {
                    cout << "   " << endl;
                    cout << "   " << endl;
                    cout << "##### WIN PAIR " << h3D_winPairInfo->GetYaxis()->GetBinLabel( currentWinPairId+1 ) << endl;
                }
                //            if ( valueByNameWP( "Nf*Nb" ) < 0.00001 )
                //                continue;
                FLAG_QA_TRIPLETS = true;
                bool isGoodCalc = finalCalcPerTriplet(); // in case of holes in the acceptance - win pair will be skipped
                if ( !isGoodCalc )
                    continue;
                FLAG_QA_TRIPLETS = false;



                // subsamples
                for ( int iSub = 0; iSub < nSubs; iSub++ )
                {
                    // cout << "iSub=" << iSub << ", iType=" << iType << ", iCW=" << iCW << ", iWinPair=" << iWinPair << endl;
                    currentSubId = iSub;
                    finalCalcPerTriplet();
                }

                // calc mean and errors:
                for( int iObs = 0; iObs < nTripletObs; iObs++ )
                {
                    // mean
                    double mean = histCalcTripletObs->GetBinContent( iObs+1, iTriplet+1, nSubs+1 );

                    // stdDev:
                    double std_dev = 0;
                    for ( int iSub = 0; iSub < nSubs; iSub++)
                    {
                        double value = histCalcTripletObs->GetBinContent( iObs+1, iTriplet+1, iSub+1 );
                        float diff = value - mean;
                        std_dev += diff*diff;
                    }
                    std_dev = sqrt(std_dev / nSubs / (nSubs-1) );

                    //                cout << "iObs = " << iObs << ", obs = " << obsNames[iObs] << ": mean +/- std_dev: " << mean << "   " << std_dev << endl;

                    histCalcTripletObs->SetBinError( iObs+1, iTriplet+1, nSubs+1, std_dev );

                    //                cout << wpObs[iWinPair]->histCalcObs->GetXaxis()->GetBinLabel( bin+1 ) << " = " << mean << ", std_dev = " << std_dev << endl;

                    cout << ">>> eta separations: FB = " << currenEtaSepFB << ", FC = " << currenEtaSepFC << ", BC = " << currenEtaSepBC
                         << ", c3 = " << mean << " +- " << std_dev //<< endl;
                      << ", win ids: F = " << current_winIdF << ", B = " << current_winIdB << ", C = " << current_winIdC
                        << endl;


                    // add points to graphs:
                    int nGrP = grVsTripletId[iObs]->GetN();
                    grVsTripletId[iObs]->SetPoint( nGrP, iTriplet, mean );
                    grVsTripletId[iObs]->SetPointError( nGrP, 0, std_dev );


                    // TRIPLET WIN SEPARATIONS:
                    int w1 = -1;
                    int w2 = -1;
                    int wMoving = -1;
                    if( fabs( current_winIdF - current_winIdB) == 1 && current_winIdC != current_winIdF && current_winIdC != current_winIdB )
                    {
                        w1 = current_winIdF;
                        w2 = current_winIdB;
                        wMoving = current_winIdC;
                    }
                    if( fabs( current_winIdF - current_winIdC) == 1 && current_winIdB != current_winIdF && current_winIdB != current_winIdC )
                    {
                        w1 = current_winIdF;
                        w2 = current_winIdC;
                        wMoving = current_winIdB;
                    }
                    if( fabs( current_winIdB - current_winIdC) == 1 && current_winIdF != current_winIdB && current_winIdF != current_winIdC )
                    {
                        w1 = current_winIdB;
                        w2 = current_winIdC;
                        wMoving = current_winIdF;
                    }


                    if( fabs( w1 - w2) == 1 )
                    {
                        double middleId = (w1 + w2) / 2;
                        double xPosOfThirdWin = wMoving - middleId;
                        int nGrP = grTriplets_TwoGap0_vs_Third[iObs]->GetN();
                        grTriplets_TwoGap0_vs_Third[iObs]->SetPoint( nGrP, xPosOfThirdWin, mean );
                        grTriplets_TwoGap0_vs_Third[iObs]->SetPointError( nGrP, 0, std_dev );
                    };



                    // add points to graphs:
//                    int nGrP = grVsDeltaEta[iObs]->GetN();
//                    grVsDeltaEta[iObs]->SetPoint( nGrP, currenEtaSep, mean );
//                    grVsDeltaEta[iObs]->SetPointError( nGrP, 0, std_dev );

//                    grVsWinId[iObs]->SetPoint( nGrP, iWinPair, mean );
//                    grVsWinId[iObs]->SetPointError( nGrP, 0, std_dev );



//                    gr2D_dEta_dPhi[iObs]->SetPoint( nGrP, currenEtaSep, currenPhiSep, mean );
//                    gr2D_dEta_dPhi[iObs]->SetPointError( nGrP, 0, 0, std_dev );//currenEtaSep, currenPhiSep, std_dev );
//                    //                cout << " gr2D_dEta_dPhi[iObs]->GetN() = " << gr2D_dEta_dPhi[iObs]->GetN() << endl;
//                    //                cout << "iWinPair = " << iWinPair << ", obs = " << obsNames[iObs] << ", currenEtaSep = " << currenEtaSep << ", currenPhiSep = " << currenPhiSep << ", mean = " << mean << endl;

//                    integrals[iObs] += mean;
                }

                //            int aa;
                //            cin >> aa;
            }  // end of loop over TRIPLETS
        } // END OF TRIPLETS CALC













//        cout << "BEFORE HIST2D" << endl;

        // fill 2D hist
        TH1::AddDirectory(kFALSE); // to suppress Warning in <TFile::Append>: Replacing existing TH1, suggested in https://root-forum.cern.ch/t/switching-off-warning-in-tfile-append-replacing-existing-th1/32041
        if(1)for( int iObs = 0; iObs < nObs; iObs++ )
        {
            double minEdge = -(eRangeMax-eRangeMin) + eSizeNum - eSizeNum/2;
            double maxEdge = (eRangeMax-eRangeMin) - eSizeNum + eSizeNum/2;
//            cout << "minEdge = " << minEdge << ", maxEdge = " << maxEdge << endl;


            TString strHist2D = Form("h2D_dEta_dPhi_%s", obsNames[iObs]);
            h2D_dEta_dPhi[iObs] = new TH2D( strHist2D, ";#Delta#eta;#Delta#varphi", 2*nEtaWins-1, minEdge, maxEdge, 2*nPhiWins-1, -TMath::TwoPi(), TMath::TwoPi() );
            h2D_dEta_dPhi_binError2[iObs] = new TH2D( strHist2D+"_binError2", ";#Delta#eta;#Delta#varphi", 2*nEtaWins-1, minEdge, maxEdge, 2*nPhiWins-1, -TMath::TwoPi(), TMath::TwoPi() );
            h2D_dEta_dPhi_binCounts[iObs] = new TH2D( strHist2D+"_binCounts", ";#Delta#eta;#Delta#varphi", 2*nEtaWins-1, minEdge, maxEdge, 2*nPhiWins-1, -TMath::TwoPi(), TMath::TwoPi() );

            TH2D *hist = h2D_dEta_dPhi[iObs];
            TH2D *histError2 = h2D_dEta_dPhi_binError2[iObs];
            TH2D *histCounts = h2D_dEta_dPhi_binCounts[iObs];

            double eSep, pSep, value;
            for ( int iP = 0; iP < gr2D_dEta_dPhi[iObs]->GetN(); iP++)
            {
                gr2D_dEta_dPhi[iObs]->GetPoint( iP, eSep, pSep, value );
                hist->Fill( eSep, pSep, value );
                histCounts->Fill( eSep, pSep, 1 );

                double err = gr2D_dEta_dPhi[iObs]->GetErrorZ( iP );
//                cout << " err = " << err << endl;
                histError2->Fill( eSep, pSep, err*err );
            }
            for ( int i = 0; i < hist->GetNbinsX(); i++)
                for ( int j = 0; j < hist->GetNbinsY(); j++)
                {
                    double binContent = hist->GetBinContent( i+1, j+1 );
                    double sumBinError2 = histError2->GetBinContent( i+1, j+1 );
                    double nCounts = histCounts->GetBinContent( i+1, j+1 );
                    if( nCounts == 0 )
                    {
                        //                    cout << "AHTUNG!!! Somehow, in h2D_dEta_dPhi_binCounts, nCounts == 0! " << endl;
                        //                    int aa;
                        //                    cin >> aa;
                    }
                    else
                    {
                        hist->SetBinContent( i+1, j+1, binContent / nCounts );
                        hist->SetBinError( i+1, j+1, sqrt( sumBinError2 / nCounts / nCounts ) );
                    }
                }
            hist->SetMarkerStyle(24);

        }


    } // end of calc()





    // ###########
    void createBinNameBinIdMap() // const char* valueBinName, const char* nEvBinName )
    {
        // single win:
        for( int i=0; i < h3D_singleWinInfo->GetNbinsX(); i++ )
        {
            const char *strName = h3D_singleWinInfo->GetXaxis()->GetBinLabel(i+1);
            mapVarIdByNameSingleWin.insert( pair<string,int>( strName, i ) );
        }

        // winPairs:
        for( int i=0; i < h3D_winPairInfo->GetNbinsX(); i++ )
        {
            const char *strName = h3D_winPairInfo->GetXaxis()->GetBinLabel(i+1);
            mapVarIdByNameWinPairs.insert( pair<string,int>( strName, i ) );
        }


        // winTriplets:
        if( h3D_winTripletsInfo )
            for( int i=0; i < h3D_winTripletsInfo->GetNbinsX(); i++ )
            {
                const char *strName = h3D_winTripletsInfo->GetXaxis()->GetBinLabel(i+1);
                mapVarIdByNameWinTriplets.insert( pair<string,int>( strName, i ) );
            }

    }



    // ###########
    void fillHistWithValue( const char *obsName, double value )
    {
        int obsId = -1;
        try
        {
            obsId = mapObsIdByName.at( obsName );
        }
        catch (const std::out_of_range& oor) {
            if(FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_Obs) { std::cerr << " Obs " << obsName << ": out of range error: " << oor.what() << "this message is shown only once.\n"; FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_Obs = false; }
        }
        if( obsId >= 0 )
            histCalcObs->SetBinContent( obsId+1, currentWinPairId+1, currentSubId+1, value );
        //        cout << "obsName = " << obsName << ", mapObsIdByName.at( obsName ) = " << mapObsIdByName.at( obsName ) << ", value = " << value << endl;
    }

    // ###########
    void fillHistWithTripletValue( const char *obsName, double value )
    {
        int obsId = -1;
        try
        {
            obsId = mapTripletObsIdByName.at( obsName );
        }
        catch (const std::out_of_range& oor) {
            if(FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_Obs) { std::cerr << " Obs " << obsName << ": out of range error: " << oor.what() << "this message is shown only once.\n"; FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_Obs = false; }
        }
        if( obsId >= 0 )
            histCalcTripletObs->SetBinContent( obsId+1, currentWinTripletId+1, currentSubId+1, value );
        //        cout << "obsName = " << obsName << ", mapObsIdByName.at( obsName ) = " << mapObsIdByName.at( obsName ) << ", value = " << value << endl;
    }


    double valueByNameSW( int whichWin, const char* binName )
    {
        //        cout << binName << endl;

        double value = -1000;
        int varId = -1;   //varIdByName( binName );

        try
        {
            varId = mapVarIdByNameSingleWin.at( binName );
        }
        catch (const std::out_of_range& oor) {
            if(FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_1D) { std::cerr << " Var " << binName << ": out of range error: " << oor.what() << "this message is shown only once.\n"; FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_1D = false; }
        }


        if ( varId == -1 )
            return -1000;
        //        cout << "getValueByName: binName = " << binName << ", mapVarIdByNameSingleWin[binName]+1 = " << mapVarIdByNameSingleWin[binName]+1 << endl;

        if ( whichInputHist == 0 ) // if we have AllWinPairs histo. For this case, we EXPLICITLY find F and B win id-s (this is fast!)
        {
//            int winId = whichWin==0 ? current_winIdF : current_winIdB;
            int winId = -1;
            if ( whichWin==0 )
                winId = current_winIdF;
            else if ( whichWin==1 )
                winId = current_winIdB;
            else if ( whichWin==2 )
                winId = current_winIdC;

            if ( currentSubId < nSubs )
                value = h3D_singleWinInfo->GetBinContent( varId+1, winId+1, currentSubId+1 );
            else // i.e. we calc the full hist (=sum of all subsamples)
            {
                value = 0;
                for( int subId = 0; subId < h3D_singleWinInfo->GetNbinsZ(); subId++ )
                    value += h3D_singleWinInfo->GetBinContent( varId+1, winId+1, subId+1 );
                //                h2D_QA_FsingleWinInfo[currentWinPairId][0]->Fill(  h3D_singleWinInfo->GetYaxis()->GetBinLabel() );
                if( QA_FLAG )
                    cout << "  >> B: " << h3D_singleWinInfo->GetYaxis()->GetBinLabel( current_winIdB+1 )
                         << "  >> F: " << h3D_singleWinInfo->GetYaxis()->GetBinLabel( current_winIdF+1 ) << endl;
            }


        }
        else if ( whichInputHist == 1 || whichInputHist == 2 ) // if dEta_dPhi histo OR AllEta-Dphi histo
        {
            //            cout << ">>>>>> currenEtaSep = " << currenEtaSep << ", currenPhiSep = " << currenPhiSep << endl;
            //            cout << "nEtaWins = " << nEtaWins << ", nPhiWins = " << nPhiWins << endl;

            value = 0;

            int nDphiWP = 2*nPhiWins-1;

            int counterWinsQA = 0;

//            int wpId = 0;
            for( int e1 = 0; e1 < nEtaWins; e1++ )  // 1==Forward
                for( int p1 = 0; p1 < nPhiWins; p1++ )
                {
                    int thisWinIdF = e1*nPhiWins + p1;

                    for( int e2 = 0; e2 < nEtaWins; e2++ )  // 2==Backward
                        for( int p2 = 0; p2 < nPhiWins; p2++ )
                        {
                            int thisWinIdB = e2*nPhiWins + p2;

                            // currentWinPairId = wpId;

                            // winId for dEta-dPhi:
                            int id_dEta = nEtaWins-1 + e1-e2;
                            int id_dPhi = nPhiWins-1 + p1-p2;

                            bool flagTake = false;

                            if ( whichInputHist == 1 )
                            {
                                int thistDetaDphiPairId = nDphiWP*id_dEta + id_dPhi;
                                if ( thistDetaDphiPairId == currentWinPairId )
                                    flagTake = true;
                            }
                            else if ( whichInputHist == 2 )
                            {
                                int thistEtaWinsDphiPairId = nDphiWP*(nEtaWins*e1+e2) + id_dPhi;
                                if ( thistEtaWinsDphiPairId == currentWinPairId )
                                    flagTake = true;
                            }

                            if ( flagTake )
                            {
                                int winId = whichWin==0 ? thisWinIdF : thisWinIdB;
                                if ( currentSubId < nSubs )
                                    value += h3D_singleWinInfo->GetBinContent( varId+1, winId+1, currentSubId+1 );
                                else // i.e. we calc now the full hist (=sum of all subsamples)
                                {
                                    double sumOverSub = 0;
                                    for( int subId = 0; subId < h3D_singleWinInfo->GetNbinsZ(); subId++ )
                                        sumOverSub += h3D_singleWinInfo->GetBinContent( varId+1, winId+1, subId+1 );
                                    value += sumOverSub;
                                    if( QA_FLAG )
                                    {
                                        cout << (whichWin==0 ? "  >> F: " : "  >> B: ") << h3D_singleWinInfo->GetYaxis()->GetBinLabel( winId+1 )
                                             << ", sumOverSub = " << sumOverSub << endl;
                                    }
                                }
                            }
                        }
                }



            if (0) // value == 0 )
            {
                cout << "valueByNameSW(): value == 0!    counterWinsQA = " << counterWinsQA << ", binName = " << binName << endl;
                cout << ">>>>>> currenEtaSep = " << currenEtaSep << ", currenPhiSep = " << currenPhiSep << endl;
                cout << "nEtaWins = " << nEtaWins << ", nPhiWins = " << nPhiWins << endl;
                int aa;
                cin >> aa;
            }
            //            if( counterWinsToAverage > 0 )
            //                value /= counterWinsToAverage;

        } // end of if dEta_dPhi histo


        return value;
    }



    double valueByNameWP( const char* binName )
    {
        //        cout << binName << endl;

        double value = -1000;
        int varId = -1;   //varIdByName( binName );

        try
        {
            varId = mapVarIdByNameWinPairs.at( binName );
        }
        catch (const std::out_of_range& oor) {
            if(FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_2D) { std::cerr << " Var " << binName << ": out of range error: " << oor.what() << "this message is shown only once.\n"; FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_2D = false; }
        }


        if ( varId == -1 )
            return -1000;
        //        cout << "getValueByName: binName = " << binName << ", mapVarIdByNameWinPairs[binName]+1 = " << mapVarIdByNameWinPairs[binName]+1 << endl;
        if ( currentSubId < nSubs )
            value = h3D_winPairInfo->GetBinContent( varId+1, currentWinPairId+1, currentSubId+1 );
        else // i.e. we calc now the full hist (=sum of all subsamples)
        {
            value = 0;
            for( int subId = 0; subId < h3D_winPairInfo->GetNbinsZ(); subId++ )
                value += h3D_winPairInfo->GetBinContent( varId+1, currentWinPairId+1, subId+1 );
        }

        //            value = currentVarFullWinPairsHist->GetBinContent(  varId+1 );

        return value;
    }



    double valueByNameTriplet( const char* binName )
    {
        //        cout << binName << endl;

        double value = -1000;
        int varId = -1;   //varIdByName( binName );

        try
        {
            varId = mapVarIdByNameWinTriplets.at( binName );
        }
        catch (const std::out_of_range& oor) {
            if(FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_3D) { std::cerr << " Var " << binName << ": out of range error: " << oor.what() << "this message is shown only once. (FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_3D)\n"; FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_3D = false; }
        }


        if ( varId == -1 )
            return -1000;
        //        cout << "getValueByName: binName = " << binName << ", mapVarIdByNameWinPairs[binName]+1 = " << mapVarIdByNameWinPairs[binName]+1 << endl;
        if ( currentSubId < nSubs )
            value = h3D_winTripletsInfo->GetBinContent( varId+1, currentWinTripletId+1, currentSubId+1 );
        else // i.e. we calc now the full hist (=sum of all subsamples)
        {
            value = 0;
            for( int subId = 0; subId < h3D_winTripletsInfo->GetNbinsZ(); subId++ )
                value += h3D_winTripletsInfo->GetBinContent( varId+1, currentWinTripletId+1, subId+1 );
        }

        //            value = currentVarFullWinPairsHist->GetBinContent(  varId+1 );

        return value;
    }



    // #######
    double valueByNameWPforTriplet( int _w1, int _w2, const char* binName )
    {
        //        cout << binName << endl;

        double value = -1000;
        int varId = -1;   //varIdByName( binName );

        try
        {
            varId = mapVarIdByNameWinPairs.at( binName );
        }
        catch (const std::out_of_range& oor) {
            if(FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_3D_bis) { std::cerr << " Var " << binName << ": out of range error: " << oor.what() << "this message is shown only once. (FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_3D_bis)\n"; FLAG_CATCH_VAR_NOT_IN_MAP_SHOW_ONCE_3D_bis = false; }
        }



        if ( varId == -1 )
            return -1000;
        //        cout << "getValueByName: binName = " << binName << ", mapVarIdByNameWinPairs[binName]+1 = " << mapVarIdByNameWinPairs[binName]+1 << endl;

        // F = 0, B = 1, C = 2
        int thisWinPair = -1;
        if      ( _w1 == 0 && _w2 == 1)     thisWinPair = current_wp_for_triplets_FB;
        else if ( _w1 == 0 && _w2 == 2)     thisWinPair = current_wp_for_triplets_FC;
        else if ( _w1 == 1 && _w2 == 2)     thisWinPair = current_wp_for_triplets_BC;

        else if ( _w1 == 1 && _w2 == 0)     thisWinPair = current_wp_for_triplets_BF;
        else if ( _w1 == 2 && _w2 == 0)     thisWinPair = current_wp_for_triplets_CF;
        else if ( _w1 == 2 && _w2 == 1)     thisWinPair = current_wp_for_triplets_CB;




        if ( currentSubId < nSubs )
        {
//            value = h3D_singleWinInfo->GetBinContent( varId+1, winId+1, currentSubId+1 );
            value = h3D_winPairInfo->GetBinContent( varId+1, thisWinPair+1, currentSubId+1 );
        }
        else // i.e. we calc now the full hist (=sum of all subsamples)
        {
            value = 0;
            for( int subId = 0; subId < h3D_winPairInfo->GetNbinsZ(); subId++ )
                value += h3D_winPairInfo->GetBinContent( varId+1, thisWinPair+1, subId+1 );
        }

        //            value = currentVarFullWinPairsHist->GetBinContent(  varId+1 );

        return value;
    }









    // ###########
    bool finalCalcPerWP()
    {
        double eSep = currenEtaSep;
        double phiSep = currenPhiSep;

        // we assume that eSep=0 means windows are completely overlapped!
        bool identical = ( eSep==0 && phiSep==0 && if_Identical_FB_XY );
//        cout << "finalCalcPerWP: " << eSep << endl;


        // ##############################
        // get n Events for F-, B-only, and for both FB - to deal with Acceptance Gaps!
        double _nEventsAccFB = valueByNameWP( "NeventsBothWin" );
        if ( _nEventsAccFB < 0.0001 ) // GAP IN ACCEPTANCE (one of the two wins is zero), skip the rest
            return false;

//        cout << "_nEventsAccFB = " << _nEventsAccFB << endl;
        if( DO_GLOBAL_QA )
            QA_FLAG = true;
        double _nEventsF     = valueByNameSW(  0, "Nevents" ); // take e.g. from F win
        double _nEventsB     = valueByNameSW(  1, "Nevents" );
        QA_FLAG = false;
//        int aaa;
//        cin >> aaa;

        if ( _nEventsF < 0.0001 ) // GAP IN ACCEPTANCE, skip the rest
        {
            cout << "check _nEventsF < 0.0001: should never appear!.." << endl;
            return false;
        }
        if ( _nEventsB < 0.0001 ) // GAP IN ACCEPTANCE, skip the rest
        {
            cout << "check _nEventsB < 0.0001: never should appear!.." << endl;
            return false;
        }

        // ##############################
        // now calc the observables

        double FB = valueByNameWP( "Nf*Nb" ) / _nEventsAccFB;
        double XY = valueByNameWP( "Nx*Ny" ) / _nEventsAccFB;
        double FY = valueByNameWP( "Nf*Ny" ) / _nEventsAccFB;
        double XB = valueByNameWP( "Nb*Nx" ) / _nEventsAccFB;

        double F = valueByNameSW(  0,   "Nf" )  / _nEventsF;
        double B = valueByNameSW(  1,   "Nb" )  / _nEventsB;
        double X = valueByNameSW(  0,   "Nx" )  / _nEventsF;
        double Y = valueByNameSW(  1,   "Ny" )  / _nEventsB;

        double F2 = valueByNameSW(  0,  "Nf2" ) / _nEventsF;
        double B2 = valueByNameSW(  1,  "Nb2" ) / _nEventsB;
        double X2 = valueByNameSW(  0,  "Nx2" ) / _nEventsF;
        double Y2 = valueByNameSW(  1,  "Ny2" ) / _nEventsB;


        double PF      = valueByNameSW(  0,  "PF" )  / _nEventsF;
        double PB      = valueByNameSW(  1,  "PB" )  / _nEventsB;
        double PX      = valueByNameSW(  0,  "PX" )  / _nEventsF;
        double PY      = valueByNameSW(  1,  "PY" )  / _nEventsB;


        double nB_PX = valueByNameWP(  "nB*PX" ) / _nEventsAccFB;
        double nY_PX = valueByNameWP(  "nY*PX" ) / _nEventsAccFB;

        double PF_PB = valueByNameWP(  "PF*PB" ) / _nEventsAccFB;
        double nF_PB = valueByNameWP(  "nF*PB" ) / _nEventsAccFB;
        double nB_PF = valueByNameWP(  "nB*PF" ) / _nEventsAccFB;

        double PX_PY = valueByNameWP(  "PX*PY" ) / _nEventsAccFB;
        double nX_PY = valueByNameWP(  "nX*PY" ) / _nEventsAccFB;
//        double nY_PX = valueByNameWP(  "nY*PX" ) / _nEventsAccFB; // already have

        double nX_PX = valueByNameSW(  0,  "nX*PX" ) / _nEventsF;


        // if( _nEventsAccFB != _nEventsF || _nEventsF != _nEventsB || _nEventsAccFB != _nEventsB )
        if(0)
        {
            cout << "currentWinPairId = " << currentWinPairId << ", F = " << F << ", B = " << B << ", X = " << X << ", Y = " << Y
                 << ", F2 = " << F2 << ", B2 = " << B2 << ", X2 = " << X2 << ", Y2 = " << Y2 << ", FY = " << FY << ", XB = " << XB
                 << ", _nEventsF = " << _nEventsF << ", _nEventsB = " << _nEventsB << ", _nEventsAccFB = " << _nEventsAccFB
                 << endl;
            int aa;
            cin >> aa;
        }



        fillHistWithValue( "avF", F );
        fillHistWithValue( "avB", B );
        fillHistWithValue( "avX", X );
        fillHistWithValue( "avY", Y );

        fillHistWithValue( "FB", FB );
        fillHistWithValue( "XY", XY );
        fillHistWithValue( "FY", FY );
        fillHistWithValue( "XB", XB );

        fillHistWithValue( "avPtF", PF );
        fillHistWithValue( "avPtB", PB );
        fillHistWithValue( "avPtX", PX );
        fillHistWithValue( "avPtY", PY );


        fillHistWithValue( "nB*PX", nB_PX );
        fillHistWithValue( "nY*PX", nY_PX );

        fillHistWithValue( "PF*PB", PF_PB );
        fillHistWithValue( "nF*PB", nF_PB );
        fillHistWithValue( "nB*PF", nB_PF );


        // ###########
        // ##### R2:
        // ###########
        {
            double R2_aa = FB/F/B - 1;
            double R2_bb = XY/X/Y - 1;
            if ( identical )
            {
                R2_aa = F2/F/F -1/F - 1;
                R2_bb = X2/X/X -1/X - 1;
            }
            double R2_ab = FY/F/Y - 1;
            double R2_ba = XB/X/B - 1;

            fillHistWithValue( "corr_R2_aa", R2_aa );
            fillHistWithValue( "corr_R2_bb", R2_bb );
            fillHistWithValue( "corr_R2_ab", R2_ab );
            fillHistWithValue( "corr_R2_ba", R2_ba );
        }


        // ##############################
        // ##### coeff ratio-ratio
        // ##############################
        double rFB      = valueByNameWP(  "Nf_OVER_Nx_vs_Nb_OVER_Ny" ) / valueByNameWP(  "xy_Nevents" );
        double ratioF   = valueByNameWP(  "Nf_OVER_Nx" ) / valueByNameWP(  "xy_Nevents" );
        double ratioB   = valueByNameWP(  "Nb_OVER_Ny" ) / valueByNameWP(  "xy_Nevents" );
        {
//            double Norm = F/eSizeNum * X/eSizeDenom / (F/eSizeNum + X/eSizeDenom);
            double Norm = (X+Y)/2 /eSizeDenom;
            double corr_rr_direct = Norm * ( rFB/ratioF/ratioB    - 1 );
            double corr_rr_formula = Norm * (FB/F/B + XY/X/Y - FY/F/Y - XB/X/B   );

            if (identical)
            {
                corr_rr_direct = Norm * ( rFB/ratioF/ratioB - 1/F - 1/X   - 1 );
                corr_rr_formula = Norm * (F2/F/F + X2/X/X - FY/F/Y - XB/X/B - 1/F - 1/X );
            }
            fillHistWithValue( "corr_rr_direct", corr_rr_direct );
            fillHistWithValue( "corr_rr_formula", corr_rr_formula );
            if(0)
            {
                cout << "corr_rr_formula = " << corr_rr_formula
                      << ", F = " << F << ", X = " << X << ", Y = " << Y
                      << ", F2 = " << F2 << ", X2 = " << X2 << ", FY = " << FY << ", XB = " << XB
                      << ", _nEventsAccFB = " << _nEventsAccFB
                      << endl;
                        int aa;
                        cin >> aa;
            }

            // when denominator is in full acceptance:
            {
                double FX = FY;

                double corr_rr_FULL_ETA_DENOM_formula = Norm * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X  ) ;
//                double corr_rr_FULL_ETA_DENOM_formula = Norm * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -(X2-X*X)/X/X  ) ;
                double corr_rr_FULL_ETA_DENOM_direct =  Norm * ( rFB/ratioF/ratioB - 1   -1/X  ) ;

                // March 2023: what if instead of 1/<X> take (<X2>-<X>^2) / <X^2>
                // poluchaetsya bredovo -> disable
//                double corr_rr_FULL_ETA_DENOM_formula_with_minusC2X = Norm * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -(X2-X*X)/X/X  ) ;

                if (identical)
                {
                    corr_rr_FULL_ETA_DENOM_formula = Norm * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X  - 1/F ) ;
                    corr_rr_FULL_ETA_DENOM_direct =  Norm * ( rFB/ratioF/ratioB - 1   -1/X   - 1/F ) ;
                }

                //            fillHistWithValue( "corr_rr_FULL_ETA_DENOM_formula", (F+B)/2 * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X) );
                //            fillHistWithValue( "corr_rr_FULL_ETA_DENOM_formula", F*B/(F+B) * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X) );
                fillHistWithValue( "corr_rr_FULL_ETA_DENOM_formula", corr_rr_FULL_ETA_DENOM_formula );
//                fillHistWithValue( "corr_rr_FULL_ETA_DENOM_formula_with_minusC2X", corr_rr_FULL_ETA_DENOM_formula_with_minusC2X );
                fillHistWithValue( "corr_rr_FULL_ETA_DENOM_direct", corr_rr_FULL_ETA_DENOM_direct );

            }
        }



        // ##############################
        // ##### coeff ratio-pT
        // ##############################
        {
            // B/Y vs pT_X
            // direct
            //            double Nb_OVER_Ny_PX = getValueByName(  "Nb_OVER_Ny*PX" );
            double PxPy_avPx      = valueByNameWP(  "PxPy_avPx" );
            double Nb_OVER_Ny_vs_Px = valueByNameWP(  "Nb_OVER_Ny_vs_avPx" );

//            double Norm = X/eSizeNum;
            double Norm = (X+Y)/2 /eSizeDenom;

            fillHistWithValue( "corr_rPt_direct", Norm * ( Nb_OVER_Ny_vs_Px/PxPy_avPx/ratioB - 1 ) );

            // formula

            // new ratio-meanPt formula - July 2022: when expansion for <pT> is done as <nB*PX>/.. + <nY*nX>/.. - <nB*nX>/.. - <nY*PX>/..
            //            double _avSumPtInEvX      = getValueByName(  "sumPtAllEvX" )  / _nEvents;
            double rPt_full_formula = Norm * (  nB_PX /B/PX  + XY/X/Y - XB/X/B - nY_PX /Y/PX  ); // here even if wins are identical, +1/X-1/X cancels!
            if (identical)
            {
                //                double subtrX = identical ? 1/X : 0;
                rPt_full_formula = Norm * (  nB_PX /B/PX  + X2/X/X - XB/X/B - nX_PX /X/PX );
            }
            fillHistWithValue( "corr_rPt_formula", rPt_full_formula );
        }



        // ####### nu_dyn, sigma, avX, avY
        if(1)
        {
            //            double numerator = XY - X * Y;
            //            double denominator_bCorr = X2 - X*X;


            if ( F != 0 && B != 0 && X != 0 && Y != 0 )
            {
                // FB
                double nu_dyn_FB = ( F2 - F ) / F / F
                        + ( B2 - B ) / B / B
                        - 2*( FB/F/B);
                if ( identical )
                    nu_dyn_FB = 0;
                fillHistWithValue( "nu_dyn_FB", nu_dyn_FB );
                fillHistWithValue( "sigma_FB",  nu_dyn_FB / (1./F + 1./B) + 1. );
                if(0)cout << "F2 = " << F2
                          << ", F = " << F
                          << ", B2 = " << B2
                          << ", B = " << B
                          << ", FB = " << FB
                          << ", sigmaFB = " << nu_dyn_FB / (1./F + 1./B) + 1.
                          << endl;

                // << "( F2 - F ) / F / F = " << ( F2 - F ) / F / F
                ////                     << ", ( B2 - B ) / B / B = " << ( B2 - B ) / B / B
                ////                     << ", FB/F/B = " << FB/F/B
                //                     << "sigma_FB = " << nu_dyn_FB / (1./F + 1./B) + 1. << endl;

                // XY
                double nu_dyn_XY = ( X2 - X ) / X / X
                        + ( Y2 - Y ) / Y / Y
                        - 2*( XY/X/Y);
                if ( identical )
                    nu_dyn_XY = 0;
                fillHistWithValue( "nu_dyn_XY", nu_dyn_XY );
                fillHistWithValue( "sigma_XY", nu_dyn_XY / (1./X + 1./Y) + 1.  );

                //                cout << "sigma_XY = " << nu_dyn_XY / (1./X + 1./Y) + 1. << endl;




                // FY
                double nu_dyn_FY = ( F2 - F ) / F / F
                        + ( Y2 - Y ) / Y / Y
                        - 2*( FY/F/Y);
                fillHistWithValue( "nu_dyn_FY", nu_dyn_FY );
                fillHistWithValue( "sigma_FY", nu_dyn_FY / (1./F + 1./Y) + 1. );

                // XB
                double nu_dyn_XB = ( X2 - X ) / X / X
                        + ( B2 - B ) / B / B
                        - 2*( XB/X/B);
                fillHistWithValue( "nu_dyn_XB", nu_dyn_XB );
                fillHistWithValue( "sigma_XB", nu_dyn_XB / (1./X + 1./B) + 1. );
            }
        }


        // ##### coeff avPt-avPt FB formula (new expansion, July 2022)
        {
//            double Norm = F/eSizeNum;
            double Norm = (F+B)/2 /eSizeNum;

            double avPtF_avPtB_formula = Norm * ( PF_PB/PF/PB + FB/F/B - nF_PB/F/PB  - nB_PF/B/PF );

            if (identical)
            {
                double PF2 =    valueByNameSW(  0,  "PF2" ) / _nEventsF;
                double nF_PF =  valueByNameSW(  0,  "nF*PF" ) / _nEventsF;
                double piF2 =   valueByNameSW(  0,  "piF2" ) / _nEventsF;

                //                double subtrX = identical ? 1/X : 0;
                avPtF_avPtB_formula = Norm * ( (PF2 - piF2)/PF/PF + F2/F/F  - nF_PF/F/PF  - nF_PF/F/PF     - 1/F + 1/F + 1/F );
            }

            fillHistWithValue( "avPtF_avPtB_formula", avPtF_avPtB_formula );

            // ##### direct
            {
                double _Pf = valueByNameWP(  "PfPb_avPf" ) / valueByNameWP(  "fb_Nevents" );
                double _Pb = valueByNameWP(  "PfPb_avPb" ) / valueByNameWP(  "fb_Nevents" );
                double _Pf_Pb = valueByNameWP(  "PfPb_avPf_avPb" ) / valueByNameWP(  "fb_Nevents" );

                fillHistWithValue( "avPtF_avPtB_direct", Norm * ( _Pf_Pb/(_Pf * _Pb) - 1 ) );
                //            cout << " >> Pf_Pb = " << Pf_Pb << ", Pf = " << Pf << ", Pb = " << Pb << endl;
            }
        }




        // ##### coeff avPt-avPt XY formula
        {
//            double Norm = X/eSizeNum;
            double Norm = (X+Y)/2 /eSizeDenom;
            double avPtX_avPtY_formula = Norm * ( PX_PY/PX/PY + XY/X/Y - nX_PY/X/PY  - nY_PX/Y/PX );

            if (identical)
            {
                double PX2 =    valueByNameSW(  0,  "PX2" ) / _nEventsF;
                double piX2 =   valueByNameSW(  0,  "piX2" ) / _nEventsF;

                //                double subtrX = identical ? 1/X : 0;
                avPtX_avPtY_formula = Norm * ( (PX2 - piX2)/PX/PX + X2/X/X  - nX_PX/X/PX  - nX_PX/X/PX     - 1/X + 1/X + 1/X );
            }

            fillHistWithValue( "avPtX_avPtY_formula", avPtX_avPtY_formula );

            // ##### direct
            {
                double _Px = valueByNameWP(  "PxPy_avPx" ) / valueByNameWP(  "xy_Nevents" );
                double _Py = valueByNameWP(  "PxPy_avPy" ) / valueByNameWP(  "xy_Nevents" );
                double _Px_Py = valueByNameWP(  "PxPy_avPx_avPy" ) / valueByNameWP(  "xy_Nevents" );

                double avPtX_avPtY_direct = Norm * ( _Px_Py/(_Px * _Py) - 1 );
                if (identical)
                {
//                    avPtX_avPtY_direct = Norm * ( _Px_Py/(_Px * _Py) - 1 - 1/X ); // this is most probably wrong, need <<pT>^2>!..
                }

                fillHistWithValue( "avPtX_avPtY_direct", avPtX_avPtY_direct );
                //            cout << " >> Pf_Pb = " << Pf_Pb << ", Pf = " << Pf << ", Pb = " << Pb << endl;
            }
        }








        // ##### dptdpt
        if(0)
        {
            double _avPtAllEvF     = valueByNameSW(  0,  "sumPtAllEvF" )  / ( F*_nEventsF );
            double _avPtAllEvB     = valueByNameSW(  0,  "sumPtAllEvB" )  / ( B*_nEventsB );
            double _piFpjB = valueByNameWP(  "piFpjB" );
            double _nF_sum_pB = valueByNameWP(  "nF*PB" ); // /2; // /2 is TMP!!!
            double _nB_sum_pF = valueByNameWP(  "nB*PF" ); // /2;

            // from /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/results/_common_files/routine.h:169
            fillHistWithValue( "coeff_dptdpt", (
                                   ( _piFpjB - _nF_sum_pB*_avPtAllEvF - _nB_sum_pF*_avPtAllEvB)/ (FB*_nEventsAccFB)
                                   + _avPtAllEvF*_avPtAllEvB
                                   ) / (_avPtAllEvF*_avPtAllEvB) );
        }

        //        }



        // ### coeff ptpt with means in denom
        if(0){

            fillHistWithValue( "coeff_ptpt_CHECK_VV_formula",
                               (valueByNameSW(  0,  "Nf" )
                                *valueByNameSW(  1,  "Nb" )
                                *valueByNameWP(  "PF*PB" )
                                -
                                valueByNameSW(  1,  "Nb" )
                                *valueByNameSW(  0,  "PF" )
                                *valueByNameWP(  "nF*PB" )// /2 // /2 is TMP!!!
                                -
                                valueByNameSW(  0,  "Nf" )
                                *valueByNameSW(  1,  "PB" )
                                *valueByNameWP(  "nB*PF" )// /2 // /2 is TMP!!!
                                +
                                valueByNameSW(  0,  "PF" )
                                *valueByNameSW(  1,  "PB" )
                                *valueByNameWP(  "Nf*Nb" )
                                ) / (
                                   valueByNameSW(  0,  "PF" )
                                   *valueByNameSW(  1,  "PB" )
                                   *valueByNameWP(  "Nf*Nb" )
                                   )
                               );
        }

        // ### coeff ptpt BOZEK_2017_USING_QM17_RESULTS_1704.02777.pdf
        if(0)
        {
            double _avPtAllEvF     = valueByNameSW(  0,  "sumPtAllEvF" )  / ( F*_nEventsF );
            double _avPtAllEvB     = valueByNameSW(  1,  "sumPtAllEvB" )  / ( B*_nEventsB );

            double Pf = valueByNameWP(  "PfPb_Pf" ) ;
            double Pb = valueByNameWP(  "PfPb_Pb" ) ;
            double Pf_Pb = valueByNameWP(  "PfPb_Pf_Pb" );

            double _pipjF = valueByNameSW(  0,  "pipjF" );
            double _pipjB = valueByNameSW(  1,  "pipjB" );

            double secondTermF     = valueByNameSW(  0,  "(nF-1)*PF" );
            double secondTermB     = valueByNameSW(  1,  "(nB-1)*PB" );

            double nn_1F     = valueByNameSW(  0, "nF*(nF-1)" );
            double nn_1B     = valueByNameSW(  1, "nB*(nB-1)" );


            double C_Bozek_F_nPairs_OUTSIDE_sum = ( _pipjF*2 - 2*_avPtAllEvF*secondTermF + _avPtAllEvF*_avPtAllEvF*nn_1F ) / nn_1F;
            double C_Bozek_B_nPairs_OUTSIDE_sum = ( _pipjB*2 - 2*_avPtAllEvB*secondTermB + _avPtAllEvB*_avPtAllEvB*nn_1B ) / nn_1B;

            //            double C_Bozek_F = pipjF_avPerEv/nEventsF*2 - 2*avPtAllEvF*secondTermFnew/nEventsF + avPtAllEvF*avPtAllEvF;
            //            double C_Bozek_B = pipjB_avPerEv/nEventsB*2 - 2*avPtAllEvB*secondTermBnew/nEventsB + avPtAllEvB*avPtAllEvB;

            fillHistWithValue( "coeff_ptpt_BOZEK_2017", (Pf_Pb - Pf*Pb)/sqrt( C_Bozek_F_nPairs_OUTSIDE_sum * C_Bozek_B_nPairs_OUTSIDE_sum ) );

            //            if(0)cout <<  ">>> C_Bozek_B_nPairs_OUTSIDE_sum: " << _pipjB << " "<< _avPtAllEvB << " " << _avPtAllEvB << " " << nn_1B << ", coeff_ptpt_BOZEK_2017=" << coeff_ptpt_BOZEK_2017 << endl;

            fillHistWithValue( "C_BOZEK_F", C_Bozek_F_nPairs_OUTSIDE_sum );
            fillHistWithValue( "C_BOZEK_B", C_Bozek_B_nPairs_OUTSIDE_sum );
        }









        return true;
    } // end of final calc for PAIRS





    // ###########
    bool finalCalcPerTriplet()
    {
//        double eSep = currenEtaSep;
//        double phiSep = currenPhiSep;

        // we assume that eSep=0 means windows are completely overlapped!
//        bool identical = ( eSep==0 && phiSep==0 && if_Identical_FB_XY );
//        cout << "finalCalcPerWP: " << eSep << endl;


        // ##############################
        // get n Events for F-, B-only, and for both FB - to deal with Acceptance Gaps!
        double _nEventsAccFB = valueByNameWPforTriplet( 0, 1, "NeventsBothWin" );
        if ( _nEventsAccFB < 0.0001 ) // GAP IN ACCEPTANCE (one of the two wins is zero), skip the rest
            return false;

//        cout << "_nEventsAccFB = " << _nEventsAccFB << endl;
        if( DO_GLOBAL_QA )
            QA_FLAG = true;
        double _nEventsF     = valueByNameSW(  0, "Nevents" ); // take e.g. from F win
        double _nEventsB     = valueByNameSW(  1, "Nevents" );
        double _nEventsC     = valueByNameSW(  2, "Nevents" );
        QA_FLAG = false;
//        int aaa;
//        cin >> aaa;

        if ( _nEventsF < 0.0001 ) // GAP IN ACCEPTANCE, skip the rest
        {
            cout << "check _nEventsF < 0.0001: should never appear!.." << endl;
            return false;
        }
        if ( _nEventsB < 0.0001 ) // GAP IN ACCEPTANCE, skip the rest
        {
            cout << "check _nEventsB < 0.0001: never should appear!.." << endl;
            return false;
        }
        if ( _nEventsC < 0.0001 ) // GAP IN ACCEPTANCE, skip the rest
        {
            cout << "check _nEventsC < 0.0001: never should appear!.." << endl;
            return false;
        }

        // ##############################
        // now calc the observables
//        current_wp_for_triplets_FB
//        current_wp_for_triplets_FC
//        current_wp_for_triplets_BC

//        current_wp_for_triplets_BF
//        current_wp_for_triplets_CF
//        current_wp_for_triplets_CB
        double FB = valueByNameWPforTriplet( 0, 1, "Nf*Nb" ) / _nEventsAccFB;
        double FC = valueByNameWPforTriplet( 0, 2, "Nf*Nb" ) / _nEventsAccFB;
        double BC = valueByNameWPforTriplet( 1, 2, "Nf*Nb" ) / _nEventsAccFB;

        double BF = valueByNameWPforTriplet( 1, 0, "Nf*Nb" ) / _nEventsAccFB;
        double CF = valueByNameWPforTriplet( 2, 0, "Nf*Nb" ) / _nEventsAccFB;
        double CB = valueByNameWPforTriplet( 2, 1, "Nf*Nb" ) / _nEventsAccFB;



        double FXX = valueByNameSW( 0, "Nf*Nx*Nx" ) / _nEventsAccFB;
        double BXX = valueByNameSW( 1, "Nb*Nx*Nx" ) / _nEventsAccFB;  // temporary, VALID ONLY IF USE FULL-ETA DENOM (X)!
        double CXX = valueByNameSW( 2, "Nc*Nx*Nx" ) / _nEventsAccFB;  // temporary, VALID ONLY IF USE FULL-ETA DENOM (X)!

        //        double XY = valueByNameWP( "Nx*Ny" ) / _nEventsAccFB;
//        double FY = valueByNameWP( "Nf*Ny" ) / _nEventsAccFB;
        double FX = valueByNameSW( 0, "Nf*Nx" ) / _nEventsAccFB;
        double BX = valueByNameWPforTriplet( 0, 1, "Nb*Nx" ) / _nEventsAccFB;
        double CX = valueByNameWPforTriplet( 0, 2, "Nc*Nx" ) / _nEventsAccFB;

        double F = valueByNameSW(  0,   "Nf" )  / _nEventsF;
        double B = valueByNameSW(  1,   "Nb" )  / _nEventsB;
        double C = valueByNameSW(  2,   "Nc" )  / _nEventsC;
        double X = valueByNameSW(  0,   "Nx" )  / _nEventsF;
//        double Y = valueByNameSW(  1,   "Ny" )  / _nEventsB;

//        double F2 = valueByNameSW(  0,  "Nf2" ) / _nEventsF;
//        double B2 = valueByNameSW(  1,  "Nb2" ) / _nEventsB;
        double X2 = valueByNameSW(  0,  "Nx2" ) / _nEventsF;
        double X3 = valueByNameSW(  0,  "Nx3" ) / _nEventsF;
//        double Y2 = valueByNameSW(  1,  "Ny2" ) / _nEventsB;

        double FBC = valueByNameTriplet( "Nf*Nb*Nc" ) / _nEventsAccFB;
        double FBX = valueByNameTriplet( "Nf*Nb*Nx" ) / _nEventsAccFB;
        double FXC = valueByNameTriplet( "Nf*Nx*Nc" ) / _nEventsAccFB;
        double XBC = valueByNameTriplet( "Nx*Nb*Nc" ) / _nEventsAccFB;


        if(FLAG_QA_TRIPLETS) cout << " !!! CX = " << CX << ", F=" << F << ", B=" << B << ", C=" << C
              << ", current_wp_for_triplets_FB = " << current_wp_for_triplets_FB
              << ", current_wp_for_triplets_FC = " << current_wp_for_triplets_FC
              << ", current_wp_for_triplets_BC = " << current_wp_for_triplets_BC

              << ", current_wp_for_triplets_BF = " << current_wp_for_triplets_BF
              << ", current_wp_for_triplets_CF = " << current_wp_for_triplets_CF
              << ", current_wp_for_triplets_CB = " << current_wp_for_triplets_CB
             << endl;




        double Norm = /*(X+Y)/2*/ X /eSizeDenom;
        // ##### c3 (March 2023)
        {
            double c3 = Norm*Norm *
                    (
                        FBC/F/B/C - ( FBX/F/B/X + FXC/F/X/C + XBC/X/B/C )
                        + 2*( FXX/F/X/X + BXX/B/X/X + CXX/C/X/X )
                        - 2*( FX/F/X + BX/B/X + CX/C/X )
                        - X3/X/X/X - 3*X2/X/X + 6
                        + 1/X/X  // this is the CORRECTION for poissonian baseline
                        );

            if(FLAG_QA_TRIPLETS) cout << "Norm = " << Norm << ",  X = " << X << ",  F = " << F << ",  1/X/X = " << +1/X/X << ",  C(X)/X3 = " << X3/X/X/X - 3*X2/X/X +2 << endl;
            if(FLAG_QA_TRIPLETS) cout << "       FBC/F/B/C = " << FBC/F/B/C << endl;
            if(FLAG_QA_TRIPLETS) cout << "       FBX/F/B/X + FXC/F/X/C + XBC/X/B/C = " << FBX/F/B/X + FXC/F/X/C + XBC/X/B/C << endl;
            if(FLAG_QA_TRIPLETS) cout << "       FXX/F/X/X + BXX/B/X/X + CXX/C/X/X = " << FXX/F/X/X + BXX/B/X/X + CXX/C/X/X << endl;
            if(FLAG_QA_TRIPLETS) cout << "       FXX/F/X/X = " << FXX/F/X/X << ", BXX/B/X/X = " << BXX/B/X/X << ", CXX/C/X/X = " << CXX/C/X/X << ", compare with: 1/X = " << 1/X << endl;
            if(FLAG_QA_TRIPLETS) cout << "       FX/F/X + BX/B/X + CX/C/X = " << FX/F/X + BX/B/X + CX/C/X << endl;
            if(FLAG_QA_TRIPLETS) cout << "       FX/F/X = " << FX/F/X << ",  BX/B/X = " << BX/B/X << ",  CX/C/X = " << CX/C/X << endl;
            if(FLAG_QA_TRIPLETS) cout << "       - X3/X/X/X - 3*X2/X/X + 6 = " << - X3/X/X/X - 3*X2/X/X + 6 << endl;



            // ######### !!! TRY NON-RATIO - ordinary cumulant (March 29)
            double c3_notRatio = Norm*Norm *
                    (
                        FBC/F/B/C - ( FB/F/B + FC/F/C + BC/B/C )
                        + 2
                        );



//            if (identical)
//            {
//                double PX2 =    valueByNameSW(  0,  "PX2" ) / _nEventsF;
//                double piX2 =   valueByNameSW(  0,  "piX2" ) / _nEventsF;

//                //                double subtrX = identical ? 1/X : 0;
//                avPtX_avPtY_formula = Norm * ( (PX2 - piX2)/PX/PX + X2/X/X  - nX_PX/X/PX  - nX_PX/X/PX     - 1/X + 1/X + 1/X );
//            }

//            fillHistWithTripletValue( "c3", c3 );
//            fillHistWithTripletValue( "c3", FBC/F/B/C );
            fillHistWithTripletValue( "c3", c3_notRatio );

            // ##### c3 direct:
            {
                double _nEventsAccXY = valueByNameWPforTriplet( 0, 1, "xy_Nevents" );
cout << "        >>> _nEventsAccXY = " << _nEventsAccXY  << ", _nEventsAccFB = " << _nEventsAccFB << endl;
                double rF  = valueByNameWPforTriplet( 0, 1, "Nf_OVER_Nx" ) / _nEventsAccXY;
                double rFB  = valueByNameWPforTriplet( 0, 1, "Nf_OVER_Nx_vs_Nb_OVER_Ny" ) / _nEventsAccXY;
                double rFC  = valueByNameWPforTriplet( 0, 2, "Nf_OVER_Nx_vs_Nc_OVER_Nz" ) / _nEventsAccFB;
                double rBC  = valueByNameWPforTriplet( 1, 2, "Nb_OVER_Ny_vs_Nc_OVER_Nz" ) / _nEventsAccFB;
                double rFBC = valueByNameTriplet( "Nfbc_OVER_NNN" ) / _nEventsAccFB;

                if(FLAG_QA_TRIPLETS) cout << "rF = " << rF << endl;
                if(FLAG_QA_TRIPLETS) cout << "rFBC = " << rFBC << endl;
                if(FLAG_QA_TRIPLETS) cout << "rFB = " << rFB << endl;
                if(FLAG_QA_TRIPLETS) cout << "rFC = " << rFC << endl;
                if(FLAG_QA_TRIPLETS) cout << "rBC = " << rBC << endl;
                if(FLAG_QA_TRIPLETS) cout << "1/X/X = " << 1/X/X << endl;


                double c3_direct = Norm*Norm *
                        (
                            rFBC / rF / rF / rF  // TMP! use only rF here, while should be also rB and rC!
                            - (rFB + rFC + rBC) / rF / rF
                            + 2
//                            - X3/X/X/X - 3*X2/X/X + 6
//                            + 1/X/X  // this is the CORRECTION for poissonian baseline
                            + 1/X/X // this is the CORRECTION for poissonian baseline
                            );
                if(FLAG_QA_TRIPLETS) cout << "######## c3_direct = " << c3_direct << endl;
                fillHistWithTripletValue( "c3_direct", c3_direct );
//                fillHistWithTripletValue( "c3_direct", rFB / rF / rF );
//                fillHistWithTripletValue( "c3_direct", rFBC / rF / rF / rF );
//                fillHistWithTripletValue( "c3_direct", rBC / rF/ rF  );
//                fillHistWithTripletValue( "c3_direct", rFBC   );
//                fillHistWithTripletValue( "c3_direct", rF  );

//            }
            // ##### c3 CHECK POISSONIAN BASELINE:
//            {
//                double _nEventsAccXY = valueByNameWPforTriplet( 0, 1, "xy_Nevents" );

//                double rF  = valueByNameWPforTriplet( 0, 1, "Nf_OVER_Nx" ) / _nEventsAccFB;
                double one_OVER_NN  = valueByNameWPforTriplet( 0, 1, "1_OVER_NN" ) / _nEventsAccFB;
//                double rFC  = valueByNameWPforTriplet( 0, 2, "Nf_OVER_Nx_vs_Nc_OVER_Nz" ) / _nEventsAccFB;
//                double rBC  = valueByNameWPforTriplet( 1, 2, "Nb_OVER_Ny_vs_Nc_OVER_Nz" ) / _nEventsAccFB;
                double one_OVER_NNN = valueByNameTriplet( "1_OVER_NNN" ) / _nEventsAccFB;

//                if(FLAG_QA_TRIPLETS) cout << "rF = " << rF << endl;
                if(FLAG_QA_TRIPLETS) cout << "one_OVER_NNN = " << one_OVER_NNN << endl;
                if(FLAG_QA_TRIPLETS) cout << "one_OVER_NN = " << one_OVER_NN << endl;
                if(FLAG_QA_TRIPLETS) cout << "<F>*<B> = " << F*B << ", <FB> = " << FB << endl;
                if(FLAG_QA_TRIPLETS) cout << "<F>*<B>*<C> = " << F*B*C << ", <FBC> = " << FBC << endl;
//                if(FLAG_QA_TRIPLETS) cout << "rFC = " << rFC << endl;
//                if(FLAG_QA_TRIPLETS) cout << "rBC = " << rBC << endl;
//                if(FLAG_QA_TRIPLETS) cout << "1/X/X = " << 1/X/X << endl;


                double c3_POISSON = Norm*Norm *
                        (
                            one_OVER_NNN / (1/X) / (1/X) / (1/X)  // TMP! use only rF here, while should be also rB and rC!
                            - 3*one_OVER_NN / (1/X) / (1/X)
                            + 2
                            + 1/X/X // this is the CORRECTION for poissonian baseline
                            );
                if(FLAG_QA_TRIPLETS) cout << "######## c3_POISSON = " << c3_POISSON << endl;
                if(FLAG_QA_TRIPLETS) cout << "##### one_OVER_NNN / (1/X) / (1/X) / (1/X)  = " << one_OVER_NNN / (1/X) / (1/X) / (1/X) << ", rFBC / rF / rF / rF = rFBC / rF / rF / rF = " << rFBC / rF / rF / rF << endl;
                fillHistWithTripletValue( "c3_direct_Poisson", c3_POISSON );
//                fillHistWithTripletValue( "c3_direct_Poisson", one_OVER_NN / (1/X) / (1/X)  );
//                fillHistWithTripletValue( "c3_direct_Poisson", one_OVER_NNN / (1/X) / (1/X) / (1/X)  );
//                fillHistWithTripletValue( "c3_direct_Poisson", one_OVER_NN / (1/X) / (1/X) );
//                fillHistWithTripletValue( "c3_direct_Poisson", one_OVER_NNN * F * F * F );
//                fillHistWithTripletValue( "c3_direct_Poisson", F/X );

            }

        }




        return true;

    } // end of final calc for TRIPLETS




}; // end of class



#endif

