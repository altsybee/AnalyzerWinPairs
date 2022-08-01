#ifndef WinPairFilling_cxx
#define WinPairFilling_cxx


#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THn.h"

#include "TDirectory.h"
#include "TString.h"

//#include "WinPairBase.h"

#include <iostream>
//#include <map>
//#include <string>

#define stringify( name ) #name

using namespace std;

// #########################
// ##### available 1D vars:
enum vars_1D
{
    _Nevents_    ,   // 1D var
//    _f_Nevents_,

    // ##### ratio-ratio correlations by approx. formula:
    _Nf_    ,   // 1D var
    _Nb_    ,   // 1D var
    _Nx_    ,   // 1D var
    _Ny_    ,   // 1D var

    _Nf2_   ,   // 1D var
    _Nb2_   ,   // 1D var
    _Nx2_   ,   // 1D var
    _Ny2_   ,   // 1D var


    // ##### for pt-pt FB
    _PF_ ,   // 1D var
    _PB_ ,   // 1D var

    _nF_PF_, // for same-window case    // 1D var

    // ##### for pt-pt XY
    _PX_ ,   // 1D var
    _PY_ ,   // 1D var
    //    nY_PX ,  // added above!

    _nX_PX_, // for same-window case   // 1D var

    _PF2_ ,   // 1D var
    _PX2_ ,   // 1D var

    // for corrections:
    _piF2_ ,   // 1D var
    _piX2_ ,   // 1D var

    // total number of vars:
    _nVars_1D_

};

// #########################
// ##### available vars:
enum vars_2D
{
    _Nevents_bothWinAcceptance_    ,   // 1D var
//    _x_Nevents_  ,
    _fb_Nevents_ ,
    _xy_Nevents_ ,
//    _fy_Nevents_ ,

    // ##### ratio-ratio correlations by approx. formula:
    _Nf_Nb_ ,       // "pure" 2D var
    _Nx_Ny_ ,       // "pure" 2D var
    _Nf_Ny_ ,       // "pure" 2D var
    _Nb_Nx_ ,       // "pure" 2D var

    // ratio-ratio correlations by direct formula:
    _Nf_OVER_Nx_vs_Nb_OVER_Ny_  ,
    _Nf_OVER_Nx_                ,   // 1D var
    _Nb_OVER_Ny_                ,   // 1D var


    // ##### for r-Pt
    _Nb_OVER_Ny_vs_avPx_,
//    _PfNb_Pf_           ,
    //    sumPtAllEvX      , // replaced by PX
    _nY_PX_            ,       // "pure" 2D var
    _nB_PX_            ,       // "pure" 2D var

    // ##### for pt-pt FB
    _PfPb_avPf_avPb_,
    _PfPb_avPf_,
    _PfPb_avPb_,
    _PF_PB_ ,       // "pure" 2D var
    _nF_PB_ ,       // "pure" 2D var
    _nB_PF_ ,       // "pure" 2D var


    // ##### for pt-pt XY
    _PxPy_avPx_avPy_ ,
    _PxPy_avPx_ ,
    _PxPy_avPy_ ,
    _PX_PY_ ,       // "pure" 2D var
    _nX_PY_ ,       // "pure" 2D var
    //    nY_PX ,  // added above!

    // total number of vars:
    _nVars_2D_

};



//// #########################
//// ##### available vars:
//enum vars
//{
//    _Nevents_    ,   // 1D var
//    _f_Nevents_,
//    _x_Nevents_  ,
//    _fb_Nevents_ ,
//    _xy_Nevents_ ,
////    _fy_Nevents_ ,

//    // ##### ratio-ratio correlations by approx. formula:
//    _Nf_Nb_ ,       // "pure" 2D var
//    _Nx_Ny_ ,       // "pure" 2D var
//    _Nf_Ny_ ,       // "pure" 2D var
//    _Nb_Nx_ ,       // "pure" 2D var

//    _Nf_    ,   // 1D var
//    _Nb_    ,   // 1D var
//    _Nx_    ,   // 1D var
//    _Ny_    ,   // 1D var

//    _Nf2_   ,   // 1D var
//    _Nb2_   ,   // 1D var
//    _Nx2_   ,   // 1D var
//    _Ny2_   ,   // 1D var

//    // ratio-ratio correlations by direct formula:
//    _Nf_OVER_Nx_vs_Nb_OVER_Ny_  ,
//    _Nf_OVER_Nx_                ,   // 1D var
//    _Nb_OVER_Ny_                ,   // 1D var


//    // ##### for r-Pt
//    _Nb_OVER_Ny_vs_avPx_,
////    _PfNb_Pf_           ,
//    //    sumPtAllEvX      , // replaced by PX
//    _nY_PX_            ,       // "pure" 2D var
//    _nB_PX_            ,       // "pure" 2D var

//    // ##### for pt-pt FB
//    _PfPb_avPf_avPb_,
//    _PfPb_avPf_,
//    _PfPb_avPb_,
//    _PF_ ,   // 1D var
//    _PB_ ,   // 1D var
//    _PF_PB_ ,       // "pure" 2D var
//    _nF_PB_ ,       // "pure" 2D var
//    _nB_PF_ ,       // "pure" 2D var

//    _nF_PF_, // for same-window case    // 1D var

//    // ##### for pt-pt XY
//    _PxPy_avPx_avPy_ ,
//    _PxPy_avPx_ ,
//    _PxPy_avPy_ ,
//    _PX_ ,   // 1D var
//    _PY_ ,   // 1D var
//    _PX_PY_ ,       // "pure" 2D var
//    _nX_PY_ ,       // "pure" 2D var
//    //    nY_PX ,  // added above!

//    _nX_PX_, // for same-window case   // 1D var

//    _PF2_ ,   // 1D var
//    _PX2_ ,   // 1D var

//    // for corrections:
//    _piF2_ ,   // 1D var
//    _piX2_ ,   // 1D var

//    // total number of vars:
//    _nVars_

//};



const char* enum1D_toStr(int e) //throw()
{
    switch (e)
    {
    case _Nevents_      :              return "Nevents"    ;
    case _Nf_     :                    return "Nf"    ;
    case _Nb_     :                    return "Nb"    ;
    case _Nx_     :                    return "Nx"    ;
    case _Ny_     :                    return "Ny"    ;
    case _Nf2_    :                    return "Nf2";
    case _Nb2_    :                    return "Nb2";
    case _Nx2_    :                    return "Nx2";
    case _Ny2_    :                    return "Ny2";
    case _PF_  :                       return "PF" ;
    case _PB_  :                       return "PB" ;
    case _nF_PF_ :                     return "nF*PF";
    case _PX_  :                       return "PX" ;
    case _PY_  :                       return "PY" ;
    case _nX_PX_ :                     return "nX*PX";
    case _PF2_  :                      return "PF2" ;
    case _PX2_  :                      return "PX2" ;
    case _piF2_  :                     return "piF2" ;
    case _piX2_  :                     return "piX2" ;
    default: //throw std::invalid_argument("Unimplemented item");
    {
        cout << "AHTUNG! Unimplemented var name in enum1D_toStr!" << endl;
        int aa;
        cin >> aa;
    }

    }

    return "";
}



//constexpr
const char* enum2D_toStr(int e) //throw()
{
    switch (e)
    {
//    case _f_Nevents_    :              return "f_Nevents";
//    case _x_Nevents_    :              return "x_Nevents"  ;
    case _Nevents_bothWinAcceptance_:  return "NeventsBothWin" ;
    case _fb_Nevents_   :              return "fb_Nevents" ;
    case _xy_Nevents_   :              return "xy_Nevents" ;
//    case _fy_Nevents_   :              return "fy_Nevents" ;
    case _Nf_Nb_  :                    return "Nf*Nb" ;
    case _Nx_Ny_  :                    return "Nx*Ny" ;
    case _Nf_Ny_  :                    return "Nf*Ny" ;
    case _Nb_Nx_  :                    return "Nb*Nx" ;
    case _Nf_OVER_Nx_vs_Nb_OVER_Ny_ :  return "Nf_OVER_Nx_vs_Nb_OVER_Ny"  ;
    case _Nb_OVER_Ny_vs_avPx_ :        return "Nb_OVER_Ny_vs_avPx";
    case _Nf_OVER_Nx_               :  return "Nf_OVER_Nx"                ;
    case _Nb_OVER_Ny_               :  return "Nb_OVER_Ny"                ;
//    case _PfNb_Pf_            :        return "PfNb_Pf"           ;
    case _nY_PX_             :         return "nY*PX"            ;
    case _nB_PX_             :         return "nB*PX"            ;
    case _PfPb_avPf_avPb_ :            return "PfPb_avPf_avPb";
    case _PfPb_avPf_ :                 return "PfPb_avPf";
    case _PfPb_avPb_ :                 return "PfPb_avPb";
    case _PF_PB_  :                    return "PF*PB" ;
    case _nF_PB_  :                    return "nF*PB" ;
    case _nB_PF_  :                    return "nB*PF" ;
    case _PxPy_avPx_avPy_  :           return "PxPy_avPx_avPy" ;
    case _PxPy_avPx_  :                return "PxPy_avPx" ;
    case _PxPy_avPy_  :                return "PxPy_avPy" ;
    case _PX_PY_  :                    return "PX*PY" ;
    case _nX_PY_  :                    return "nX*PY" ;
    default: //throw std::invalid_argument("Unimplemented item");
    {
        cout << "AHTUNG! Unimplemented var name in enum2D_toStr!" << endl;
        int aa;
        cin >> aa;
    }

    }

    return "";
}


////constexpr
//const char* enumToStr(int e) //throw()
//{
//    switch (e)
//    {
//    case _Nevents_      :              return "Nevents"    ;
//    case _f_Nevents_    :              return "f_Nevents";
//    case _x_Nevents_    :              return "x_Nevents"  ;
//    case _fb_Nevents_   :              return "fb_Nevents" ;
//    case _xy_Nevents_   :              return "xy_Nevents" ;
////    case _fy_Nevents_   :              return "fy_Nevents" ;
//    case _Nf_Nb_  :                    return "Nf*Nb" ;
//    case _Nx_Ny_  :                    return "Nx*Ny" ;
//    case _Nf_Ny_  :                    return "Nf*Ny" ;
//    case _Nb_Nx_  :                    return "Nb*Nx" ;
//    case _Nf_     :                    return "Nf"    ;
//    case _Nb_     :                    return "Nb"    ;
//    case _Nx_     :                    return "Nx"    ;
//    case _Ny_     :                    return "Ny"    ;
//    case _Nf2_    :                    return "Nf2";
//    case _Nb2_    :                    return "Nb2";
//    case _Nx2_    :                    return "Nx2";
//    case _Ny2_    :                    return "Ny2";
//    case _Nf_OVER_Nx_vs_Nb_OVER_Ny_ :  return "Nf_OVER_Nx_vs_Nb_OVER_Ny"  ;
//    case _Nf_OVER_Nx_               :  return "Nf_OVER_Nx"                ;
//    case _Nb_OVER_Ny_               :  return "Nb_OVER_Ny"                ;
//    case _Nb_OVER_Ny_vs_avPx_ :        return "Nb_OVER_Ny_vs_avPx";
////    case _PfNb_Pf_            :        return "PfNb_Pf"           ;
//    case _nY_PX_             :         return "nY*PX"            ;
//    case _nB_PX_             :         return "nB*PX"            ;
//    case _PfPb_avPf_avPb_ :            return "PfPb_avPf_avPb";
//    case _PfPb_avPf_ :                 return "PfPb_avPf";
//    case _PfPb_avPb_ :                 return "PfPb_avPb";
//    case _PF_  :                       return "PF" ;
//    case _PB_  :                       return "PB" ;
//    case _PF_PB_  :                    return "PF*PB" ;
//    case _nF_PB_  :                    return "nF*PB" ;
//    case _nB_PF_  :                    return "nB*PF" ;
//    case _nF_PF_ :                     return "nF*PF";
//    case _PxPy_avPx_avPy_  :           return "PxPy_avPx_avPy" ;
//    case _PxPy_avPx_  :                return "PxPy_avPx" ;
//    case _PxPy_avPy_  :                return "PxPy_avPy" ;
//    case _PX_  :                       return "PX" ;
//    case _PY_  :                       return "PY" ;
//    case _PX_PY_  :                    return "PX*PY" ;
//    case _nX_PY_  :                    return "nX*PY" ;
//    case _nX_PX_ :                     return "nX*PX";
//    case _PF2_  :                      return "PF2" ;
//    case _PX2_  :                      return "PX2" ;
//    case _piF2_  :                     return "piF2" ;
//    case _piX2_  :                     return "piX2" ;
//    default: //throw std::invalid_argument("Unimplemented item");
//    {
//        cout << "AHTUNG! Unimplemented var name!" << endl;
//        int aa;
//        cin >> aa;
//    }

//    }

//    return "";
//}




// ####################################################################################################################
// class which collects all window pairs and performs filling of TH3D histogram with e-by-e calculated values
// ####################################################################################################################

class WinPairWrapper
{
public:
    // windows setup:
    int nEtaWins;
    int nPhiWins;
    int partTypes[4];   // F,B, X,Y  // 0 - NO PID SELECTION
    int partCharges[4]; // F,B, X,Y  // 0 - NO CHARGE SELECTION
    double ptWin[2];

    // June 2021: new for ratio-ratio when denominator is in full acceptance
    bool fullAcceptanceForDenom;
    //    double etaForDenomMin;
    //    double etaForDenomMax;
    double etaMin;
    double etaMax;



    // e-by-e histograms:
    TH2D *hist_n[4]; //!
    TH2D *hist_pt[4];  //!
    TH2D *hist_w_wMinus1[4];  //!
    TH2D *hist_w_wMinus1_pt[4];  //!
    TH2D *hist_pt2[4];  //!

    TH2D *QA_hist_n[4];  //!
    TH2D *QA_hist_pt[4];  //!


    TH1D *hMetaInfo; //!

    // histos with SINGLE win info:
    TH3D *hSingleWin; //!

    // histos with win PAIR info:
    TH3D *hAllWins; //!
    TH3D *hDetaDphi; //!
    TH3D *hAllEtaDphi; //!


    bool flagHistAllWins;
    bool flagHistDetaDphi;
    bool flagHistAllEtaDphi;


    TH2I *hAccMap; //!
//    TH3D *hSPEC_DetaDphi; //!
    //    THnD *hDeltaEta; //!
//    TString strAnLevel;

    // some vars to make flexible hist filling (i.e. don't fill if bin name is not in varNames array):
    int currentSubsampleId;
    int currentSingleWinId;
    int currentWinPairId;
    int currentDetaDphiPairId;
    int currentEtaWinsDphiPairId;


//    TString strListName;

    // ############
    WinPairWrapper()
    {
        hMetaInfo = 0x0;
        hSingleWin = 0x0;

        hAllWins = 0x0;
        hDetaDphi = 0x0;
        hAllEtaDphi = 0x0;

        flagHistAllWins = true;
        flagHistDetaDphi = true;
        flagHistAllEtaDphi = true;


        hAccMap = 0x0;

        nEtaWins = 0;
        nPhiWins = 0;
    }

    ~WinPairWrapper() {}

    void setup( const char* strPrefix, int cBin, int nSub,
                int *_pTypes, int *_pCharges, int _nEtaWins, int _nPhiWins, double *_etaRange, double *_ptRange, bool *whichHistosToTake, TList *_outputList
                , bool _fullAcceptanceForDenom = false //, double _etaForDenomMin = -0.8, double _etaForDenomMax = 0.8
            )
    {
//        strListName = _outputList->GetName();
//        cout << "TList name: " << strListName << endl;

        nEtaWins = _nEtaWins;
        nPhiWins = _nPhiWins;

        //        nWPs = (nEtaWins*nPhiWins) * (nEtaWins*nPhiWins);
        //        if ( nWPs >= MAX_N_WIN_PAIRS )
        //        {
        //            cout << "AHTUNG!!! nWPs >= MAX_N_WIN_PAIRS" << endl;
        //            int aa;
        //            cin >> aa;
        //        }
        ptWin[0] = _ptRange[0];
        ptWin[1] = _ptRange[1];

        for( int i = 0; i < 4; i++ )
        {
            partTypes[i] = _pTypes[i];
            partCharges[i] = _pCharges[i];
        }

        fullAcceptanceForDenom = _fullAcceptanceForDenom;
        //        etaForDenomMin = _etaForDenomMin;
        //        etaForDenomMax = _etaForDenomMax;
        etaMin = _etaRange[0];
        etaMax = _etaRange[1];

        flagHistAllWins     = whichHistosToTake[0];
        flagHistDetaDphi    = whichHistosToTake[1];
        flagHistAllEtaDphi  = whichHistosToTake[2];


        //        double eBpos = ( _eWins[1] + _eWins[0] )/2;
        //        double eFpos = ( _eWins[3] + _eWins[2] )/2;

        //        eSep[iWin] = round( (eFpos - eBpos)*100 ) / 100;
        //        eFsize[iWin] = round( (_eWins[3] - _eWins[2])*100 ) / 100;
        //        eBsize[iWin] = round( (_eWins[1] - _eWins[0])*100 ) / 100;

        //        double phiBpos = ( _phiWins[1] + _phiWins[0] )/2;
        //        double phiFpos = ( _phiWins[3] + _phiWins[2] )/2;
        //        phiSep[iWin] = round( (phiFpos - phiBpos)*100 ) / 100;
        //        phiFsize[iWin] = round( (_phiWins[3] - _phiWins[2])*100 ) / 100;
        //        phiBsize[iWin] = round( (_phiWins[3] - _phiWins[2])*100 ) / 100;


        // ############
//        int _nVars2D = _nVars_2D_;
        cout << "_nVars_1D_ = " << _nVars_1D_ << ", _nVars_2D_ = " << _nVars_2D_ << endl;


        TString strAnLevel = Form("%s", strPrefix);

        TString strPID = Form( "pid_%d_%d_%d_%d_charge_%d_%d_%d_%d",
                               partTypes[0], partTypes[1], partTypes[2], partTypes[3],
                partCharges[0], partCharges[1], partCharges[2], partCharges[3] );

        // e-by-e histos:
        TString strWinLabels[] = { "F", "B", "X", "Y" };
        for( int i = 0; i < 4; i++ )
        {
            int _nEtaBins = nEtaWins;
            int _nPhiBins = nPhiWins;
            if ( fullAcceptanceForDenom && (i==2 || i==3) ) // use full acceptance for X and Y
            {
                _nEtaBins = 1;
                _nPhiBins = 1;
            }
            TString strHistName = Form( "%s_hist_n_%s_cBin%d_", strPID.Data(), strAnLevel.Data(), cBin) + strWinLabels[i];
            hist_n[i] = new TH2D( strHistName, strHistName,     _nEtaBins, etaMin, etaMax,   _nPhiBins, 0, TMath::TwoPi() );
            QA_hist_n[i] = new TH2D( strHistName+"QA", strHistName+"QA",     _nEtaBins, etaMin, etaMax,   _nPhiBins, 0, TMath::TwoPi() );

            strHistName = Form( "%s_hist_pt_%s_cBin%d_", strPID.Data(), strAnLevel.Data(), cBin) + strWinLabels[i];
            hist_pt[i] = new TH2D( strHistName, strHistName,    _nEtaBins, etaMin, etaMax,   _nPhiBins, 0, TMath::TwoPi() );
            QA_hist_pt[i] = new TH2D( strHistName+"QA", strHistName+"QA",    _nEtaBins, etaMin, etaMax,   _nPhiBins, 0, TMath::TwoPi() );


            strHistName = Form( "%s_hist_w_wMinus1_%s_cBin%d_", strPID.Data(), strAnLevel.Data(), cBin) + strWinLabels[i];
            hist_w_wMinus1[i] = new TH2D( strHistName, strHistName,     _nEtaBins, etaMin, etaMax,   _nPhiBins, 0, TMath::TwoPi() );

            strHistName = Form( "%s_hist_w_wMinus1_pt_%s_cBin%d_", strPID.Data(), strAnLevel.Data(), cBin) + strWinLabels[i];
            hist_w_wMinus1_pt[i] = new TH2D( strHistName, strHistName,  _nEtaBins, etaMin, etaMax,   _nPhiBins, 0, TMath::TwoPi() );

            strHistName = Form( "%s_hist_pt2_pt_%s_cBin%d_", strPID.Data(), strAnLevel.Data(), cBin) + strWinLabels[i];
            hist_pt2[i] = new TH2D( strHistName, strHistName,           _nEtaBins, etaMin, etaMax,   _nPhiBins, 0, TMath::TwoPi() );
        }

        cout << " >>> etaMin = " << etaMin << ", etaMax = " << etaMax << endl;

        if(0)for( int e1 = 0; e1 < nEtaWins; e1++ )
            for( int p1 = 0; p1 < nPhiWins; p1++ )
            {
                cout <<  "   >> F: eBin=" << hist_n[0]->GetXaxis()->GetBinCenter(e1+1)
                     <<  "   pBin=" << hist_n[0]->GetYaxis()->GetBinCenter(p1+1) << endl;
            }


        double etaSize  = ( etaMax - etaMin ) / nEtaWins;
        double min_dEta = -( etaMax - etaMin ) + etaSize;
        double max_dEta = ( etaMax - etaMin ) - etaSize;

        // find max dPhi:
        double phiSize  = TMath::TwoPi() / nPhiWins;
        double min_dPhi = -TMath::TwoPi() + phiSize;
        double max_dPhi = TMath::TwoPi() - phiSize;


        // hist with meta info:
        TString strHistMetaInfo = Form("hMetaInfo_%s_cBin%d", strAnLevel.Data(), cBin);
        hMetaInfo = new TH1D( strHistMetaInfo, strHistMetaInfo
                              , 10, -0.5, 10-0.5 );
        int metaBinId = 1;
        hMetaInfo->GetXaxis()->SetBinLabel( metaBinId++, "nEtaWins" );
        hMetaInfo->GetXaxis()->SetBinLabel( metaBinId++, "nPhiWins" );
        hMetaInfo->GetXaxis()->SetBinLabel( metaBinId++, "etaRangeMin" );
        hMetaInfo->GetXaxis()->SetBinLabel( metaBinId++, "etaRangeMax" );
        hMetaInfo->GetXaxis()->SetBinLabel( metaBinId++, "etaSize" );
        hMetaInfo->GetXaxis()->SetBinLabel( metaBinId++, "phiSize" );
        hMetaInfo->GetXaxis()->SetBinLabel( metaBinId++, "ptMin" );
        hMetaInfo->GetXaxis()->SetBinLabel( metaBinId++, "ptMax" );
        hMetaInfo->GetXaxis()->SetBinLabel( metaBinId++, "fullAcceptanceForDenom" );
        hMetaInfo->Fill( "nEtaWins", nEtaWins );
        hMetaInfo->Fill( "nPhiWins", nPhiWins );
        hMetaInfo->Fill( "etaRangeMin",  _etaRange[0] );
        hMetaInfo->Fill( "etaRangeMax",  _etaRange[1] );
        hMetaInfo->Fill( "etaSize",   etaSize );
        hMetaInfo->Fill( "phiSize", phiSize );
        hMetaInfo->Fill( "ptMin", _ptRange[0] );
        hMetaInfo->Fill( "ptMax", _ptRange[1] );
        hMetaInfo->Fill( "fullAcceptanceForDenom", _fullAcceptanceForDenom );


//        void setup( const char* strPrefix, int cBin, int nSub,
//                    int *_pTypes, int *_pCharges, int _nEtaWins, int _nPhiWins, double *_etaRange, double *_ptRange, bool *whichHistosToTake, TList *_outputList
//                    , bool _fullAcceptanceForDenom = false //, double _etaForDenomMin = -0.8, double _etaForDenomMax = 0.8
//                )





        // single win histos:
        TString strHist1DName = Form("hSingleWin_%s_cBin%d", strAnLevel.Data(), cBin);
        hSingleWin = new TH3D( strHist1DName, strHist1DName
                              , _nVars_1D_, -0.5, _nVars_1D_-0.5
                              , nEtaWins*nPhiWins, -0.5, nEtaWins*nPhiWins-0.5
                              , nSub, -0.5, nSub-0.5
                              );


        // win pair info histos:
        TString strHistName = Form("hAllWins_%s_cBin%d", strAnLevel.Data(), cBin);

        int nAllEtaWPs = nEtaWins*nEtaWins;
        int nAllPhiWPs = nPhiWins*nPhiWins;
        int nAllWPs = nAllEtaWPs * nAllPhiWPs;

        // all window pairs separately:
        if( flagHistAllWins )
            hAllWins = new TH3D( strHistName, strHistName
                                  , _nVars_2D_, -0.5, _nVars_2D_-0.5
                                  , nAllWPs, -0.5, nAllWPs-0.5
                                  , nSub, -0.5, nSub-0.5
                                  );

        int nDetaWP = 2*nEtaWins-1;
        int nDphiWP = 2*nPhiWins-1;

        TString str_dEta_dPhi_HistName = Form("hDetaDphi_%s_cBin%d", strAnLevel.Data(), cBin);
        if( flagHistDetaDphi )
            hDetaDphi = new TH3D( str_dEta_dPhi_HistName, str_dEta_dPhi_HistName
                              , _nVars_2D_, -0.5, _nVars_2D_-0.5
                              , nDetaWP*nDphiWP, -0.5, nDetaWP*nDphiWP-0.5
                              , nSub, -0.5, nSub-0.5
                              );
        cout << "nDetaWP = " << nDetaWP << ", min_dEta = " << min_dEta << ", max_dEta = " << max_dEta << endl;
        cout << "nDphiWP = " << nDphiWP << endl;


        int nAllEtaDphiWP = nEtaWins*nEtaWins * nDphiWP;
        TString str_AllEta_dPhi_HistName = Form("hAllEtaDphi_%s_cBin%d", strAnLevel.Data(), cBin);
        if( flagHistAllEtaDphi )
            hAllEtaDphi = new TH3D( str_AllEta_dPhi_HistName, str_AllEta_dPhi_HistName
                                          , _nVars_2D_, -0.5, _nVars_2D_-0.5
                                          , nAllEtaDphiWP, -0.5, nAllEtaDphiWP-0.5
                                          , nSub, -0.5, nSub-0.5
                                          );




        // set labels x and prepare a map with vars
//        cout << "_nVars = " << _nVars << endl;
        for( int i = 0; i < _nVars_1D_; i++ )
        {
            const char* _varName = enum1D_toStr(i);
            if(0)cout << "i = " << i << ", _varName = " << _varName << endl;

            if( hSingleWin )
                hSingleWin->GetXaxis()->SetBinLabel( i+1, _varName ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
        }

        for( int i = 0; i < _nVars_2D_; i++ )
        {
            const char* _varName = enum2D_toStr(i);
            if(0)cout << "i = " << i << ", _varName = " << _varName << endl;

            if( flagHistAllWins )
                hAllWins->GetXaxis()->SetBinLabel( i+1, _varName ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
            if( flagHistDetaDphi )
                hDetaDphi->GetXaxis()->SetBinLabel( i+1, _varName ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
            if( flagHistAllEtaDphi )
                hAllEtaDphi->GetXaxis()->SetBinLabel( i+1, _varName ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
        }



        // set labels y for hSingleWin
        if( 1 )
        {
            int wpId = 0;
            for( int e1 = 0; e1 < nEtaWins; e1++ )
                for( int p1 = 0; p1 < nPhiWins; p1++ )
                {
                    double eMin = hist_n[0]->GetXaxis()->GetBinLowEdge(e1+1);
                    double pMin = hist_n[0]->GetYaxis()->GetBinLowEdge(p1+1);
                    double eMax = hist_n[0]->GetXaxis()->GetBinUpEdge(e1+1);
                    double pMax = hist_n[0]->GetYaxis()->GetBinUpEdge(p1+1);

                    hSingleWin->GetYaxis()->SetBinLabel( wpId+1, Form("eta_%.2f_%.2f_phi_%.2f_%.2f", eMin, eMax, pMin, pMax ) );
                    wpId++;
                }
        }

        // set labels y for hAllWins
        if( flagHistAllWins )
        {
            int wpId = 0;
            for( int e1 = 0; e1 < nEtaWins; e1++ )
                for( int p1 = 0; p1 < nPhiWins; p1++ )
                {
                    double eFmin = hist_n[0]->GetXaxis()->GetBinLowEdge(e1+1);
                    double pFmin = hist_n[0]->GetYaxis()->GetBinLowEdge(p1+1);
                    double eFmax = hist_n[0]->GetXaxis()->GetBinUpEdge(e1+1);
                    double pFmax = hist_n[0]->GetYaxis()->GetBinUpEdge(p1+1);

//                    int winId1 = nPhiWins*e1 + p1;

                    for( int e2 = 0; e2 < nEtaWins; e2++ )
                        for( int p2 = 0; p2 < nPhiWins; p2++ )
                        {
                            double eBmin = hist_n[0]->GetXaxis()->GetBinLowEdge(e2+1);
                            double pBmin = hist_n[0]->GetYaxis()->GetBinLowEdge(p2+1);
                            double eBmax = hist_n[0]->GetXaxis()->GetBinUpEdge(e2+1);
                            double pBmax = hist_n[0]->GetYaxis()->GetBinUpEdge(p2+1);

//                            int winId2 = nPhiWins*e2 + p2;

                            hAllWins->GetYaxis()->SetBinLabel( wpId+1, Form("etaB_%.2f_%.2f_phiB_%.2f_%.2f_etaF_%.2f_%.2f_phiF_%.2f_%.2f"
                                                                            , eBmin, eBmax,    pBmin, pBmax
                                                                            , eFmin, eFmax,    pFmin, pFmax
                                                                            ) );

                            // for QA (alternative binLabelling for hDetaDphi):
                            if(0)
                            {
                                int id_dEta = nEtaWins-1 + e1-e2;
                                int id_dPhi = nPhiWins-1 + p1-p2;
                                int sepId = nDphiWP*id_dEta + id_dPhi;
    //                            cout << "sepId = " << sepId << endl;
                                hDetaDphi->GetYaxis()->SetBinLabel( sepId+1, Form("dEta_%.2f_dPhi_%.2f", round( (e1-e2)*etaSize *100 )/100, round( (p1-p2)*phiSize *100 ) / 100 ) );
                            }

                            wpId++;
                        }
                }
        }


        // QA:
        if(0)
        {
            cout << "hDetaDphi->GetYaxis()->GetBinLabel( j+1 ):" << endl;
            for( int j = 0; j < hDetaDphi->GetNbinsY(); j++ )
                cout << hDetaDphi->GetYaxis()->GetBinLabel( j+1 ) << " ";
            cout << endl;
        }


        // set labels y for hDetaDphi
//        cout << "was like this:" << endl;
        if( flagHistDetaDphi )
        {
            for( int i = 0; i < nDetaWP; i++ )
                for( int j = 0; j < nDphiWP; j++ )
                {
                    int sepId = nDphiWP*i + j;
                    hDetaDphi->GetYaxis()->SetBinLabel( sepId+1, Form("dEta_%.2f_dPhi_%.2f", min_dEta + etaSize*i, min_dPhi + phiSize*j ) );
//                    cout << Form("dEta_%.2f_dPhi_%.2f", min_dEta + etaSize*i, min_dPhi + phiSize*j ) << " ";

                    if(0) cout << "i = " << i << ", j = " << j << ", sepId = " << sepId
                         << ", min_dEta + etaStep*i = " << min_dEta + etaSize*i
                         << ", min_dPhi + phiSize*j = " << min_dPhi + phiSize*j  << endl;
                }
        }
//        cout << endl;
//        int aa;
//        cin >> aa;


        // set labels y for hAllEtaDphi
        if( flagHistAllEtaDphi )
        {
            for( int e1 = 0; e1 < nEtaWins; e1++ )
            {
                double eFmin = hist_n[0]->GetXaxis()->GetBinLowEdge(e1+1);
                double eFmax = hist_n[0]->GetXaxis()->GetBinUpEdge(e1+1);

                for( int e2 = 0; e2 < nEtaWins; e2++ )
                {
                    double eBmin = hist_n[0]->GetXaxis()->GetBinLowEdge(e2+1);
                    double eBmax = hist_n[0]->GetXaxis()->GetBinUpEdge(e2+1);

                    for( int j = 0; j < nDphiWP; j++ )
                    {
                        int sepId = nDphiWP*(nEtaWins*e1+e2) + j;
                        hAllEtaDphi->GetYaxis()->SetBinLabel( sepId+1, Form("etaB_%.2f_%.2f_etaF_%.2f_%.2f_dPhi_%.2f"
                                                                        , eBmin, eBmax
                                                                        , eFmin, eFmax,    min_dPhi + phiSize*j
                                                                        ) );
                    }
                }
            }
        }

        _outputList->Add( hMetaInfo );
        _outputList->Add( hSingleWin );
        if ( flagHistAllWins    )  _outputList->Add( hAllWins );
        if ( flagHistDetaDphi   )  _outputList->Add( hDetaDphi );
        if ( flagHistAllEtaDphi )  _outputList->Add( hAllEtaDphi );


//        cout << "end of setup: nEtaWins = " << nEtaWins << ", nPhiWins = " << nPhiWins << endl;
    }
    

    // ############
    void addTrack( int pid, double eta, double phi, double pt, int charge, double weight = 1.0 )
    {
        // cout << "in addTrack: nEtaWins = " << nEtaWins << ", nPhiWins = " << nPhiWins << endl;
        // cout << "TList name: " << strListName << endl;
//        cout << "pid = " << pid << ", eta = " << eta << ", phi = " << phi << ", pt = " << pt << ", charge = " << charge << endl;
        if ( pt < ptWin[0] || pt > ptWin[1] )
            return;

        int pidAbs = abs(pid);

//        cout << "pidAbs = " << pidAbs << ", eta = " << eta << ", phi = " << phi << ", pt = " << pt << ", charge = " << charge << endl;


        // ##### backward window:
        // B
        if ( partTypes[1]==0 || pidAbs == partTypes[1] )
            if ( partCharges[1]==0 || charge == partCharges[1] )
            {
//                cout << ">>>>> B: pidAbs = " << pidAbs << ", eta = " << eta << ", phi = " << phi << ", pt = " << pt << ", charge = " << charge << endl;
                hist_n[1]->Fill( eta, phi, weight );
                hist_pt[1]->Fill( eta, phi, pt*weight );

                hist_w_wMinus1[1]->Fill( eta, phi, weight*(weight-1) );
            }


        // Y
        if ( partTypes[3]==0 || pidAbs == partTypes[3] )
            if ( partCharges[3]==0 || charge == partCharges[3] )
            {
//                cout << ">>>>> Y: pidAbs = " << pidAbs << ", eta = " << eta << ", phi = " << phi << ", pt = " << pt << ", charge = " << charge << endl;
                hist_n[3]->Fill( eta, phi, weight );
//                cout << "hist_n[3]->GetEntries() = " << hist_n[3]->GetEntries() << endl;
                hist_pt[3]->Fill( eta, phi, pt*weight );

                hist_w_wMinus1[3]->Fill( eta, phi, weight*(weight-1) );
            }

//        int aa;
//        cin >> aa;

        // ##### forward window:
        // F
        if ( partTypes[0]==0 || pidAbs == partTypes[0] )
            if ( partCharges[0]==0 || charge == partCharges[0] )
            {
//                cout << ">>>>> F: pidAbs = " << pidAbs << ", eta = " << eta << ", phi = " << phi << ", pt = " << pt << ", charge = " << charge << endl;
                hist_n[0]->Fill( eta, phi, weight );
                hist_pt[0]->Fill( eta, phi, pt*weight );

                hist_w_wMinus1[0]->Fill( eta, phi, weight*(weight-1) );
                hist_w_wMinus1_pt[0]->Fill( eta, phi, weight*(weight-1)*pt );
                hist_pt2[0]->Fill( eta, phi, pt*pt *weight*weight );
            }
        // X
        if ( partTypes[2]==0 || pidAbs == partTypes[2] )
            if ( partCharges[2]==0 || charge == partCharges[2] )
            {
//                cout << ">>>>> X: pidAbs = " << pidAbs << ", eta = " << eta << ", phi = " << phi << ", pt = " << pt << ", charge = " << charge << endl;
                hist_n[2]->Fill( eta, phi, weight );
                hist_pt[2]->Fill( eta, phi, pt*weight );

                hist_w_wMinus1[2]->Fill( eta, phi, weight*(weight-1) );
                hist_w_wMinus1_pt[2]->Fill( eta, phi, weight*(weight-1)*pt );
                hist_pt2[2]->Fill( eta, phi, pt*pt *weight*weight );
            }

//        int aa;
//        cin >> aa;
    }

    void fillWithValueHist1D( /*const char *varName,*/ int varId, double value ) //, bool forceFilling = false )
    {
        if( hSingleWin )
            hSingleWin->Fill(  varId, currentSingleWinId, currentSubsampleId, value  );
    }

    void fillWithValueHist2D( /*const char *varName,*/ int varId, double value ) //, bool forceFilling = false )
    {
        //                return;
        if( flagHistAllWins )
            hAllWins->Fill(  varId, currentWinPairId, currentSubsampleId, value  );

        if( flagHistDetaDphi )
            hDetaDphi->Fill(  varId, currentDetaDphiPairId, currentSubsampleId, value  );

        if( flagHistAllEtaDphi )
            hAllEtaDphi->Fill(  varId, currentEtaWinsDphiPairId, currentSubsampleId, value  );
    }




    // ############
    void finishEvent( int subId )
    {
        currentSubsampleId = subId;

        int nDphiWP = 2*nPhiWins-1;

//         cout << "in finishEvent(): nEtaWins = " << nEtaWins << ", nPhiWins = " << nPhiWins << endl;
        // cout << "TList name: " << strListName << endl;
        //
        // cout << "finishEvent():   hist_n[0]->GetEntries() = " << hist_n[0]->GetEntries() << ", integral: " << hist_n[0]->Integral() << endl;
        // cout << "finishEvent():   hist_n[1]->GetEntries() = " << hist_n[1]->GetEntries() << ", integral: " << hist_n[1]->Integral()<< endl;
        // cout << "finishEvent():   hist_n[2]->GetEntries() = " << hist_n[2]->GetEntries() << ", integral: " << hist_n[2]->Integral()<< endl;
        // cout << "finishEvent():   hist_n[3]->GetEntries() = " << hist_n[3]->GetEntries() << ", integral: " << hist_n[3]->Integral()<< endl;

//        int aa;
//        cin >> aa;

        int singleWinId = 0;
        for( int e1 = 0; e1 < nEtaWins; e1++ )
        {
            double center_e1 = hist_n[0]->GetXaxis()->GetBinCenter(e1+1);
            int e1Den = !fullAcceptanceForDenom ? e1 : 0;

            for( int p1 = 0; p1 < nPhiWins; p1++ )
            {

                if( hAccMap && ( !hAccMap->GetBinContent( e1+1, p1+1 )  ) )
                {
//                    cout << "hAccMap CHECK!!!" << endl;
                    singleWinId++;   // IMPORTANT!!!
                    continue;
                }

                double center_p1 = hist_n[0]->GetYaxis()->GetBinCenter(p1+1);
                int p1Den = !fullAcceptanceForDenom ? p1 : 0;

                currentSingleWinId = singleWinId;

                // F, X:
                double nF = hist_n[0]->GetBinContent(e1+1, p1+1);
                double nX = hist_n[2]->GetBinContent(e1Den+1, p1Den+1);

                double ptF = hist_pt[0]->GetBinContent(e1+1, p1+1);
                double ptX = hist_pt[2]->GetBinContent(e1Den+1, p1Den+1);


                double w_wMinus1_F = hist_w_wMinus1[0]->GetBinContent(e1+1, p1+1);
                double w_wMinus1_X = hist_w_wMinus1[2]->GetBinContent(e1Den+1, p1Den+1);

                double w_wMinus1_PF = hist_w_wMinus1_pt[0]->GetBinContent(e1+1, p1+1);
                double w_wMinus1_PX = hist_w_wMinus1_pt[2]->GetBinContent(e1Den+1, p1Den+1);

                double ptF2 = hist_pt2[0]->GetBinContent(e1+1, p1+1);
                double ptX2 = hist_pt2[2]->GetBinContent(e1Den+1, p1Den+1);


                // B, Y:
                double nB = hist_n[1]->GetBinContent(e1+1, p1+1);
                double nY = hist_n[3]->GetBinContent(e1Den+1, p1Den+1);

                double ptB = hist_pt[1]->GetBinContent(e1+1, p1+1);
                double ptY = hist_pt[3]->GetBinContent(e1Den+1, p1Den+1);

                double w_wMinus1_B = hist_w_wMinus1[1]->GetBinContent(e1+1, p1+1);
                double w_wMinus1_Y = hist_w_wMinus1[3]->GetBinContent(e1Den+1, p1Den+1);

                double w_wMinus1_PB = hist_w_wMinus1_pt[1]->GetBinContent(e1+1, p1+1);
                double w_wMinus1_PY = hist_w_wMinus1_pt[3]->GetBinContent(e1Den+1, p1Den+1);

                double ptB2 = hist_pt2[1]->GetBinContent(e1+1, p1+1);
                double ptY2 = hist_pt2[3]->GetBinContent(e1Den+1, p1Den+1);

//                double meanPtF = -1;
//                double meanPtB = -1;
//                double meanPtX = -1;
//                double meanPtY = -1;

//                if ( nF > 0 ) { meanPtF = ptF / nF; }
//                if ( nB > 0 ) { meanPtB = ptB / nB; }
//                if ( nX > 0 ) { meanPtX = ptX / nX; }
//                if ( nY > 0 ) { meanPtY = ptY / nY; }



                // QA check
                if( hAccMap && hAccMap->GetBinContent( e1+1, p1+1 ) )
                {
                    QA_hist_n[0]->Fill( center_e1, center_p1, nF );
                    QA_hist_pt[0]->Fill( center_e1, center_p1, ptF );
                }


                fillWithValueHist1D( _Nevents_, 1 );//, true ); // last argument - forcing filling the 0th bin of the hist

                // NfNb:
                fillWithValueHist1D( _Nf_    ,   nF             );
                fillWithValueHist1D( _Nb_    ,   nB             );
                fillWithValueHist1D( _Nf2_  ,   nF * nF - w_wMinus1_F     );
                fillWithValueHist1D( _Nb2_  ,   nB * nB - w_wMinus1_B     );
//                fillHistWithValue( _Nf_Nb_ ,   nF * nB      );

                // NxNy:
                fillWithValueHist1D( _Nx_      ,   nX             );
                fillWithValueHist1D( _Ny_      ,   nY             );
                fillWithValueHist1D( _Nx2_     ,   nX * nX - w_wMinus1_X   ); // !fullEtaForDenom ? _nX*_nX  : _nX*_nX - _w_wMinus1_X   );
                fillWithValueHist1D( _Ny2_     ,   nY * nY - w_wMinus1_Y   ); // !fullEtaForDenom ? _nY*_nY  : _nY*_nY - _w_wMinus1_Y   );
//                fillHistWithValue( _Nx_Ny_   ,   nX * nY      );

//                fillHistWithValue( _PF_PB_ ,      ptF*ptB );
                fillWithValueHist1D( _PF_ ,              ptF );
                fillWithValueHist1D( _PB_ ,              ptB );

//                fillHistWithValue( _PX_PY_ ,      ptX*ptY );
                fillWithValueHist1D( _PX_ ,              ptX );
                fillWithValueHist1D( _PY_ ,              ptY );


                // for new ratio-pt
                fillWithValueHist1D( _nF_PF_ , nF*ptF - w_wMinus1_PF );
                fillWithValueHist1D( _nX_PX_ , nX*ptX - w_wMinus1_PX );

                fillWithValueHist1D( _PF2_ , ptF*ptF ); //- w_wMinus1_PF );
                fillWithValueHist1D( _PX2_ , ptX*ptX ); //- w_wMinus1_PX );

                fillWithValueHist1D( _piF2_ , ptF2 );
                fillWithValueHist1D( _piX2_ , ptX2 );


                singleWinId++;
            }
        }



        if(0)for( int e1 = 0; e1 < nEtaWins; e1++ )
            for( int p1 = 0; p1 < nPhiWins; p1++ )
            {
                cout <<  "   >> F: eBin=" << hist_n[0]->GetXaxis()->GetBinCenter(e1+1)
                     <<  "   pBin=" << hist_n[0]->GetYaxis()->GetBinCenter(p1+1)
                     << "  bin content: " << hist_n[0]->GetBinContent(e1+1, p1+1) << "    ";
                cout <<  "   >> B: eBin=" << hist_n[1]->GetXaxis()->GetBinCenter(e1+1)
                     <<  "   pBin=" << hist_n[1]->GetYaxis()->GetBinCenter(p1+1)
                     << "  bin content: " << hist_n[1]->GetBinContent(e1+1, p1+1) << "    ";
                cout <<  "   >> X: eBin=" << hist_n[2]->GetXaxis()->GetBinCenter(e1+1)
                     <<  "   pBin=" << hist_n[2]->GetYaxis()->GetBinCenter(p1+1)
                     << "  bin content: " << hist_n[2]->GetBinContent(e1+1, p1+1) << "    ";
                cout <<  "   >> Y: eBin=" << hist_n[3]->GetXaxis()->GetBinCenter(e1+1)
                     <<  "   pBin=" << hist_n[3]->GetYaxis()->GetBinCenter(p1+1)
                     << "  bin content: " << hist_n[3]->GetBinContent(e1+1, p1+1) << endl;

            }


//        cout << endl;

//        int aa;
//        cin >> aa;


        // fill TH3D with win pairs info
        int wpId = 0;
        for( int e1 = 0; e1 < nEtaWins; e1++ )  // 1==Forward
        {
            double center_e1 = hist_n[0]->GetXaxis()->GetBinCenter(e1+1);

            int e1Den = !fullAcceptanceForDenom ? e1 : 0;

            for( int p1 = 0; p1 < nPhiWins; p1++ )
            {
                double center_p1 = hist_n[0]->GetYaxis()->GetBinCenter(p1+1);

                int p1Den = !fullAcceptanceForDenom ? p1 : 0;


                double nF = hist_n[0]->GetBinContent(e1+1, p1+1);
                double nX = hist_n[2]->GetBinContent(e1Den+1, p1Den+1);

                double ptF = hist_pt[0]->GetBinContent(e1+1, p1+1);
                double ptX = hist_pt[2]->GetBinContent(e1Den+1, p1Den+1);


                double w_wMinus1_F = hist_w_wMinus1[0]->GetBinContent(e1+1, p1+1);
                double w_wMinus1_X = hist_w_wMinus1[2]->GetBinContent(e1Den+1, p1Den+1);

                double w_wMinus1_PF = hist_w_wMinus1_pt[0]->GetBinContent(e1+1, p1+1);
                double w_wMinus1_PX = hist_w_wMinus1_pt[2]->GetBinContent(e1Den+1, p1Den+1);

                double ptF2 = hist_pt2[0]->GetBinContent(e1+1, p1+1);
                double ptX2 = hist_pt2[2]->GetBinContent(e1Den+1, p1Den+1);


                // QA check
//                if( hAccMap && hAccMap->GetBinContent( e1+1, p1+1 ) )
//                {
//                    QA_hist_n[0]->Fill( center_e1, center_p1, nF );
//                    QA_hist_pt[0]->Fill( center_e1, center_p1, ptF );
//                }


                // loop over pairing windows:
                for( int e2 = 0; e2 < nEtaWins; e2++ )  // 2==Backward
                {
//                    double center_e2 = hist_n[0]->GetXaxis()->GetBinCenter(e2+1);
//                    double etaSep = center_e2 - center_e1;

                    for( int p2 = 0; p2 < nPhiWins; p2++ )
                    {
                        // check if we take this bin or not!
                        if( hAccMap && ( !hAccMap->GetBinContent( e1+1, p1+1 ) || !hAccMap->GetBinContent( e2+1, p2+1 ) ) )
                        {
//                            cout << "hAccMap CHECK!!!" << endl;
                            wpId++;   // IMPORTANT!!!
                            continue;
                        }

//                        double center_p2 = hist_n[0]->GetYaxis()->GetBinCenter(p2+1);
//                        double phiSep = center_p2 - center_p1;
                        //                        continue;

                        //                        cout << "e1, p1, e2, p2 = " << e1 << ", " << p1 << ", " << e2 << ", " << p2 << endl;
                        currentWinPairId = wpId;

                        // winId for dEta-dPhi:
                        int id_dEta = nEtaWins-1 + e1-e2;
                        int id_dPhi = nPhiWins-1 + p1-p2;
                        currentDetaDphiPairId = nDphiWP*id_dEta + id_dPhi;

//                        int sepId = nDphiWP*(nEtaWins*e1+e2) + j;
                        currentEtaWinsDphiPairId = nDphiWP*(nEtaWins*e1+e2) + id_dPhi;


                        if(0)
                        {
                            double etaStep  = ( etaMax - etaMin ) / nEtaWins;
                            double min_dEta = -( etaMax - etaMin ) + etaStep;
                            double phiStep  = TMath::TwoPi() / nPhiWins;
                            double min_dPhi = -TMath::TwoPi() + phiStep;

                            cout << "CHECK: e1=" << e1 << ", e2=" << e2 << ", p1=" << p1 << ", p2=" << p2 << endl;
                            cout << "i = " << id_dEta << ", j = " << id_dPhi << ", SPEC_PairId = " << currentDetaDphiPairId
                                 << ", min_dEta + etaStep*i = " << min_dEta + etaStep*id_dEta
                                 << ", min_dPhi + phiStep*j = " << min_dPhi + phiStep*id_dPhi << endl;

                            int aa;
                            cin >> aa;
                        }


//                        if ( currentDetaDphiPairId >= sizeOfMap )
//                            cout << "AHTUNG!!! currentDetaDphiPairId > map_DetaDphi_bin.size() !" << endl;

                        int e2Den = !fullAcceptanceForDenom ? e2 : 0;
                        int p2Den = !fullAcceptanceForDenom ? p2 : 0;

                        double nB = hist_n[1]->GetBinContent(e2+1, p2+1);
                        double nY = hist_n[3]->GetBinContent(e2Den+1, p2Den+1);

                        double ptB = hist_pt[1]->GetBinContent(e2+1, p2+1);
                        double ptY = hist_pt[3]->GetBinContent(e2Den+1, p2Den+1);

                        double w_wMinus1_B = hist_w_wMinus1[1]->GetBinContent(e2+1, p2+1);
                        double w_wMinus1_Y = hist_w_wMinus1[3]->GetBinContent(e2Den+1, p2Den+1);

                        double w_wMinus1_PB = hist_w_wMinus1_pt[1]->GetBinContent(e2+1, p2+1);
                        double w_wMinus1_PY = hist_w_wMinus1_pt[3]->GetBinContent(e2Den+1, p2Den+1);

                        double ptB2 = hist_pt2[1]->GetBinContent(e2+1, p2+1);
                        double ptY2 = hist_pt2[3]->GetBinContent(e2Den+1, p2Den+1);



//                        continue;

                        double meanPtF = -1;
                        double meanPtB = -1;
                        double meanPtX = -1;
                        double meanPtY = -1;

                        if ( nF > 0 ) { meanPtF = ptF / nF; }
                        if ( nB > 0 ) { meanPtB = ptB / nB; }
                        if ( nX > 0 ) { meanPtX = ptX / nX; }
                        if ( nY > 0 ) { meanPtY = ptY / nY; }


//                        cout << " >>> nF = " << nF << ", nB = " << nB << ", ptF = " << ptF << ", ptB = " << ptB << endl;

//                        fillHistWithValue( _Nevents_, 1 );//, true ); // last argument - forcing filling the 0th bin of the hist
                        fillWithValueHist2D( _Nevents_bothWinAcceptance_ ,   1      );

                        // NfNb:
//                        fillHistWithValue( _Nf_    ,   nF             );
//                        fillHistWithValue( _Nb_    ,   nB             );
//                        fillHistWithValue( _Nf2_  ,   nF * nF - w_wMinus1_F     );
//                        fillHistWithValue( _Nb2_  ,   nB * nB - w_wMinus1_B     );
                        fillWithValueHist2D( _Nf_Nb_ ,   nF * nB      );

                        // NxNy:
//                        fillHistWithValue( _Nx_      ,   nX             );
//                        fillHistWithValue( _Ny_      ,   nY             );
//                        fillHistWithValue( _Nx2_     ,   nX * nX - w_wMinus1_X   ); // !fullEtaForDenom ? _nX*_nX  : _nX*_nX - _w_wMinus1_X   );
//                        fillHistWithValue( _Ny2_     ,   nY * nY - w_wMinus1_Y   ); // !fullEtaForDenom ? _nY*_nY  : _nY*_nY - _w_wMinus1_Y   );
                        fillWithValueHist2D( _Nx_Ny_   ,   nX * nY      );


                        // Nf_Nx Nf_Ny Nb_Nx Nb_Ny:
                        //                        if(0) fillHistWithValue( _Nf_Nx_  ,   nF * nX       );
                        fillWithValueHist2D( _Nf_Ny_  ,   nF * nY       );
                        fillWithValueHist2D( _Nb_Nx_  ,   nB * nX       );
                        //                        if(0) fillHistWithValue( _Nb_Ny_  ,   nB * nY       );


                        // FROM /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/AliForwardBackwardAnalysis.cxx:1583
                        //                        fillHistWithValue( _sumPtAllEvF_ , ptF );
                        //                        if(1) fillHistWithValue( _sumPtAllEvB_ , ptB );
                        //                        if(0) fillHistWithValue( _sumPtAllEvX_ , ptX );
                        //                        if(0) fillHistWithValue( _sumPtAllEvY_ , ptY );

                        // for dptdpt:
                        //                if(0) fillHistWithValue( _piFpjB_ , w->piFpjB );
                        fillWithValueHist2D( _nF_PB_ ,  nF * ptB );
                        fillWithValueHist2D( _nB_PF_ ,  nB * ptF );

                        // for C when av over pairs is OUTSIDE sum:
                        //                if(0) fillHistWithValue( _pipjF_ ,          w->pipjF );
                        //                        if(0) fillHistWithValue( _(nF-1)_PF_ ,      (nF-1)*ptF ); // (n-1)*sumPt
                        //                        if(0) fillHistWithValue( _nF_(nF-1)_ ,          (nF-1)*nF ); // (n-1)*n

                        //                if(0) fillHistWithValue( _pipjB_ ,         w->pipjB );
                        //                        if(0) fillHistWithValue( _(nB-1)_PB_ , (nB-1)*ptB ); // (n-1)*sumPt
                        //                        if(0) fillHistWithValue( _nB_(nB-1)_ ,     (nB-1)*nB ); // (n-1)*n



                        // to check VV comparison ptpt vs dptdpt: special terms:
                        fillWithValueHist2D( _PF_PB_ ,      ptF*ptB );
//                        fillHistWithValue( _PF_ ,              ptF );
//                        fillHistWithValue( _PB_ ,              ptB );

                        fillWithValueHist2D( _PX_PY_ ,      ptX*ptY );
//                        fillHistWithValue( _PX_ ,              ptX );
//                        fillHistWithValue( _PY_ ,              ptY );


                        // for new ratio-pt (April 2021)
                        //                        fillHistWithValue( _nY_PF_ , nY*ptF );
                        fillWithValueHist2D( _nB_PX_ , nB*ptX );
                        fillWithValueHist2D( _nY_PX_ , nY*ptX );

                        //                        fillHistWithValue( _nX_PB_ , nX*ptB );
                        //                        fillHistWithValue( _nF_PY_ , nF*ptY );
                        fillWithValueHist2D( _nX_PY_ , nX*ptY );

//                        fillHistWithValue( _nF_PF_ , nF*ptF - w_wMinus1_PF );
//                        fillHistWithValue( _nX_PX_ , nX*ptX - w_wMinus1_PX );

//                        fillHistWithValue( _PF2_ , ptF*ptF ); //- w_wMinus1_PF );
//                        fillHistWithValue( _PX2_ , ptX*ptX ); //- w_wMinus1_PX );

//                        fillHistWithValue( _piF2_ , ptF2 );
//                        fillHistWithValue( _piX2_ , ptX2 );

                        if ( nY > 0 )
                        {
                            //                            if(0) fillHistWithValue( _Nb_OVER_Ny*PF_ , nB/(double)nY * ptF );
                            //                            if(0) fillHistWithValue( _Nb_OVER_Ny*PX_ , nB/(double)nY * ptX );
                        }
                        if ( nX > 0 )
                        {
                            //                            if(0) fillHistWithValue( _Nf_OVER_Nx*PB_ , nF/(double)nX * ptB );
                            //                            if(0) fillHistWithValue( _Nf_OVER_Nx*PY_ , nF/(double)nX * ptY );
                        }



                        // NfPb:
                        if ( nB > 0 )
                        {
                            //                            if(0) fillHistWithValue( _b_Nevents_ ,    1                  );

                            //                            if(0) fillHistWithValue( _NfPb_Nf_ ,       nF             );
                            //                            if(0) fillHistWithValue( _NfPb_avPb_ ,       meanPtB             );
                            //                            if(0) fillHistWithValue( _NfPb_Nf2_ ,      nF*nF      );
                            //                            if(0) fillHistWithValue( _NfPb_avPb2_ ,      meanPtB*meanPtB      );
                            //                            if(0) fillHistWithValue( _NfPb_Nf_avPb_ ,    nF*meanPtB      );

                            //                            if(0) fillHistWithValue( _Nf_OVER_Nb_ ,    nF/(double)nB      );
                            //                            if(0) fillHistWithValue( _Nx_OVER_Nb_ ,    nX/(double)nB      );
                            //                            if(0) fillHistWithValue( _Ny_OVER_Nb_ ,    nY/(double)nB      );

                            //                            if(0) fillHistWithValue( _avPb_Nx_ ,   meanPtB*nX      );
                            //                            if(0) fillHistWithValue( _avPb_Ny_ ,   meanPtB*nY      );


                            //            // for C when av over pairs is OUTSIDE sum:
                            //            fillHistWithValue( _pipjB_ ,         pipjB );
                            //            fillHistWithValue( _(nB-1)*sum_pB_ , (_nB-1)*(_nB*_ptB) ); // (n-1)*sumPt
                            //            fillHistWithValue( _nB*(nB-1)_ ,     (_nB-1)*_nB ); // (n-1)*n

                            // for C when av over pairs is INSIDE sum: (like in GOOD_Ck_definition_STAR_2005_0504031.pdf)
                            //                            int nPairsB = (nB-1)*nB;
                            //                            if ( nPairsB > 0 )
                            //                            {
                            //                                //                fillHistWithValue( _pipjB_avPerEv_ ,         pipjB / nPairsB );
                            //                                //                fillHistWithValue( _(nB-1)*PB_avPerEv_ , (_nB-1)*_ptB / nPairsB );
                            //                            }

                        }

                        // PfNb:
                        if ( nF > 0 )
                        {
//                            fillWithValueHist2D( _f_Nevents_ ,    1                  );

                            //                            fillHistWithValue( _PfNb_avPf_ ,         meanPtF           );
                            //                            if(0) fillHistWithValue( _PfNb_Nb_ ,        nB            );
                            //                            if(0) fillHistWithValue( _PfNb_avPf2_ ,     meanPtF*meanPtF      );
                            //                            if(0) fillHistWithValue( _PfNb_Nb2_ ,      nB*nB      );
                            //                            if(0) fillHistWithValue( _PfNb_avPf_Nb_ ,    meanPtF*nB      );

                            //                            if(0) fillHistWithValue( _Nb_OVER_Nf_ ,    nB/(double)nF      );
                            //                            if(0) fillHistWithValue( _Nx_OVER_Nf_ ,    nX/(double)nF      );
                            //                            if(0) fillHistWithValue( _Ny_OVER_Nf_ ,    nY/(double)nF      );

                            //                            if(0) fillHistWithValue( _avPf_Nx_ ,   meanPtF*nX      );
                            //                            if(0) fillHistWithValue( _avPf_Ny_ ,   meanPtF*nY      );

                            //            // for C when av over pairs is OUTSIDE sum:
                            //            fillHistWithValue( _pipjF_ ,              pipjF );
                            //            fillHistWithValue( _(nF-1)*sum_pF_ ,      (_nF-1)*(_nF*_ptF) ); // (n-1)*sumPt
                            //            fillHistWithValue( _nF*(nF-1)_ ,          (_nF-1)*_nF ); // (n-1)*n

                            // for C when av over pairs is INSIDE sum:
                            //                            int nPairsF = (nF-1)*nF;
                            //                            if ( nPairsF > 0 )
                            //                            {
                            //                                //                fillHistWithValue( _pipjF_avPerEv_ ,          pipjF / nPairsF );
                            //                                //                fillHistWithValue( _(nF-1)*PF_avPerEv_ ,  (_nF-1)*_ptF / nPairsF );
                            //                            }

                        }

                        // fb:
                        if ( nF > 0 && nB > 0 )
                        {
                            if(1) fillWithValueHist2D( _fb_Nevents_ ,   1                  );

                            if(1) fillWithValueHist2D( _PfPb_avPf_ ,       meanPtF             );
                            if(1) fillWithValueHist2D( _PfPb_avPb_ ,       meanPtB             );
                            //                            if(0) fillHistWithValue( _PfPb_avPf2_ ,     meanPtF*meanPtF      );
                            //                            if(0) fillHistWithValue( _PfPb_avPb2_ ,     meanPtB*meanPtB      );
                            if(1) fillWithValueHist2D( _PfPb_avPf_avPb_ ,   meanPtF*meanPtB      );

                            //                            if(0) fillHistWithValue( _Nx_OVER_Nf_vs_avPb_ ,   nX/(double)nF *meanPtB      );
                            //                            if(0) fillHistWithValue( _Ny_OVER_Nb_vs_avPf_ ,   nY/(double)nB *meanPtF      );

                            //                            if(0) fillHistWithValue( _Nx_OVER_Nf_vs_Ny_OVER_Nb_ ,   nX/(double)nF * nY/(double)nB      );

                            //            // for dptdpt:
                            //            fillHistWithValue( _piFpjB_ ,     piFpjB );
                            //            fillHistWithValue( _nF*sum_pB_ ,  _nF*(_nB*_ptB) );
                            //            fillHistWithValue( _nB*sum_pF_ ,  _nB*(_nF*_ptF) );

                        }


                        // y:
                        if ( nY > 0 )
                        {
                            //                            if(0) fillHistWithValue( _y_Nevents_ ,    1                  );

                            //                            if(0) fillHistWithValue( _NxPy_Nx_ ,        nX             );
                            //                            if(0) fillHistWithValue( _NxPy_avPy_ ,       meanPtY             );
                            //                            if(0) fillHistWithValue( _NxPy_Nx2_ ,      nX*nX      );
                            //                            if(0) fillHistWithValue( _NxPy_avPy2_ ,      meanPtY*meanPtY      );
                            //                            if(0) fillHistWithValue( _NxPy_Nx_avPy_ ,    nX*meanPtY      );


                            //                            if(0) fillHistWithValue( _Nf_OVER_Ny_ ,    nF/(double)nY      );
                            //                            if(0) fillHistWithValue( _Nx_OVER_Ny_ ,    nX/(double)nY      );

                            //                            if(0) fillHistWithValue( _Nf_avPy_ ,   nF*meanPtY       );
                            //                            if(0) fillHistWithValue( _Nb_avPy_ ,   nB*meanPtY       );
                        }

                        // x:
                        if ( nX > 0 )
                        {
//                            fillWithValueHist2D( _x_Nevents_ ,    1                  );

                            //                            if(0) fillHistWithValue( _PxNy_avPx_ ,         meanPtX           );
                            //                            if(0) fillHistWithValue( _PxNy_Ny_ ,        nY            );
                            //                            if(0) fillHistWithValue( _PxNy_avPx2_ ,     meanPtX*meanPtX      );
                            //                            if(0) fillHistWithValue( _PxNy_Ny2_ ,      nY*nY      );
                            //                            if(0) fillHistWithValue( _PxNy_avPx_Ny_ ,    meanPtX*nY      );

                            //                            if(0) fillHistWithValue( _Nb_OVER_Nx_ ,    nB/(double)nX      );
                            //                            if(0) fillHistWithValue( _Ny_OVER_Nx_ ,    nY/(double)nX      );

                            //                            if(0) fillHistWithValue( _Nf_avPx_ ,   nF*meanPtX       );
                            //                            if(0) fillHistWithValue( _Nb_avPx_ ,   nB*meanPtX       );
                        }


                        // xy:
                        if ( nX > 0 && nY > 0 )
                        {
                            fillWithValueHist2D( _xy_Nevents_ ,   1                  );

                            fillWithValueHist2D( _Nb_OVER_Ny_ ,    nB/(double)nY      );
                            fillWithValueHist2D( _Nf_OVER_Nx_ ,    nF/(double)nX      );

                            if(1) fillWithValueHist2D( _PxPy_avPx_ ,       meanPtX             );
                            if(1) fillWithValueHist2D( _PxPy_avPy_ ,       meanPtY             );
                            //                            if(0) fillHistWithValue( _PxPy_avPx2_ ,     meanPtX*meanPtX      );
                            //                            if(0) fillHistWithValue( _PxPy_avPy2_ ,     meanPtY*meanPtY      );
                            if(1) fillWithValueHist2D( _PxPy_avPx_avPy_ ,   meanPtX*meanPtY      );

                            //                            if(0) fillHistWithValue( _Nf_OVER_Nx_vs_avPy_ ,   nF/(double)nX*meanPtY      );
                            if(1) fillWithValueHist2D( _Nb_OVER_Ny_vs_avPx_ ,   nB/(double)nY*meanPtX      );

                            fillWithValueHist2D( _Nf_OVER_Nx_vs_Nb_OVER_Ny_ ,   nF/(double)nX * nB/(double)nY      );
                        }

                        // fx:
                        if ( nF > 0 && nX > 0 )
                        {
                            //                            if(0) fillHistWithValue( _fx_Nevents_ ,   1                  );

                            //                            if(0) fillHistWithValue( _avPf_avPx_ ,   meanPtF*meanPtX      );

                            //                            if(0) fillHistWithValue( _Nb_OVER_Nf_vs_Ny_OVER_Nx_ ,   nB/(double)nF * nY/(double)nX      );
                        }
                        // fy:
                        if ( nF > 0 && nY > 0 )
                        {
//                            fillHistWithValue( _fy_Nevents_ ,   1                  );

                            //                            if(0) fillHistWithValue( _avPf_avPy_ ,   meanPtF*meanPtY      );
                            //                            if(0) fillHistWithValue( _Nx_OVER_Nf_vs_avPy_ ,   nX/(double)nF *meanPtY      );
                            //                            if(0) fillHistWithValue( _Nb_OVER_Ny_vs_avPf_ ,   nB/(double)nY *meanPtF      );
                        }
                        // bx:
                        if ( nB > 0 && nX > 0 )
                        {
                            //                            if(0) fillHistWithValue( _bx_Nevents_ ,   1                  );

                            //                            if(0) fillHistWithValue( _Pb_Px_ ,   meanPtB*meanPtX      );
                            //                            if(0) fillHistWithValue( _Nf_OVER_Nx_vs_avPb_ ,   nF/(double)nX *meanPtB      );
                            //                            if(0) fillHistWithValue( _Ny_OVER_Nb_vs_avPx_ ,   nY/(double)nB *meanPtX      );
                        }
                        // by:
                        if ( nB > 0 && nY > 0 )
                        {
                            //                            if(0) fillHistWithValue( _by_Nevents_ ,   1                  );
                            //                            if(0) fillHistWithValue( _avPb_avPy_ ,   meanPtB*meanPtY      );

                            //                            if(0) fillHistWithValue( _Nf_OVER_Nb_vs_Nx_OVER_Ny_ ,   nF/(double)nB * nX/(double)nY      );
                        }



                        wpId++;
                    } // end of p2 loop
                } // end of e2 loop
            } // end of p1 loop
        } // end of e1 loop



        // reset e-by-e histos
        for( int i = 0; i < 4; i++ )
        {
            hist_n[i]->Reset();
            hist_pt[i]->Reset();
            hist_w_wMinus1[i]->Reset();
            hist_w_wMinus1_pt[i]->Reset();
            hist_pt2[i]->Reset();
        }

    }   // end of finishEvent


//    ClassDef(WinPairWrapper, 1);

};


#endif












