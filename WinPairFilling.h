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
// ##### available vars:
enum vars
{
    _Nevents_    ,
    _f_Nevents_,
    _x_Nevents_  ,
    _fb_Nevents_ ,
    _xy_Nevents_ ,
    _fy_Nevents_ ,

    // ##### ratio-ratio correlations by approx. formula:
    _Nf_Nb_ ,
    _Nx_Ny_ ,
    _Nf_Ny_ ,
    _Nb_Nx_ ,

    _Nf_    ,
    _Nb_    ,
    _Nx_    ,
    _Ny_    ,

    _Nf2_   ,
    _Nb2_   ,
    _Nx2_   ,
    _Ny2_   ,

    // ratio-ratio correlations by direct formula:
    _Nf_OVER_Nx_vs_Nb_OVER_Ny_  ,
    _Nf_OVER_Nx_                ,
    _Nb_OVER_Ny_                ,


    // ##### for r-Pt
    _Nb_OVER_Ny_vs_avPx_,
//    _PfNb_Pf_           ,
    //    sumPtAllEvX      , // replaced by PX
    _nY_PX_            ,
    _nB_PX_            ,

    // ##### for pt-pt FB
    _PfPb_avPf_avPb_,
    _PfPb_avPf_,
    _PfPb_avPb_,
    _PF_ ,
    _PB_ ,
    _PF_PB_ ,
    _nF_PB_ ,
    _nB_PF_ ,

    _nF_PF_, // for same-window case

    // ##### for pt-pt XY
    _PxPy_avPx_avPy_ ,
    _PxPy_avPx_ ,
    _PxPy_avPy_ ,
    _PX_ ,
    _PY_ ,
    _PX_PY_ ,
    _nX_PY_ ,
    //    nY_PX ,  // added above!

    _nX_PX_, // for same-window case

    _PF2_ ,
    _PX2_ ,

    // for corrections:
    _piF2_ ,
    _piX2_ ,

    // total number of vars:
    _nVars_

};




//constexpr
const char* enumToStr(int e) //throw()
{
    switch (e)
    {
    case _Nevents_      :              return "Nevents"    ;
    case _f_Nevents_    :              return "f_Nevents";
    case _x_Nevents_    :              return "x_Nevents"  ;
    case _fb_Nevents_   :              return "fb_Nevents" ;
    case _xy_Nevents_   :              return "xy_Nevents" ;
    case _fy_Nevents_   :              return "fy_Nevents" ;
    case _Nf_Nb_  :                    return "Nf*Nb" ;
    case _Nx_Ny_  :                    return "Nx*Ny" ;
    case _Nf_Ny_  :                    return "Nf*Ny" ;
    case _Nb_Nx_  :                    return "Nb*Nx" ;
    case _Nf_     :                    return "Nf"    ;
    case _Nb_     :                    return "Nb"    ;
    case _Nx_     :                    return "Nx"    ;
    case _Ny_     :                    return "Ny"    ;
    case _Nf2_    :                    return "Nf2";
    case _Nb2_    :                    return "Nb2";
    case _Nx2_    :                    return "Nx2";
    case _Ny2_    :                    return "Ny2";
    case _Nf_OVER_Nx_vs_Nb_OVER_Ny_ :  return "Nf_OVER_Nx_vs_Nb_OVER_Ny"  ;
    case _Nf_OVER_Nx_               :  return "Nf_OVER_Nx"                ;
    case _Nb_OVER_Ny_               :  return "Nb_OVER_Ny"                ;
    case _Nb_OVER_Ny_vs_avPx_ :        return "Nb_OVER_Ny_vs_avPx";
//    case _PfNb_Pf_            :        return "PfNb_Pf"           ;
    case _nY_PX_             :         return "nY*PX"            ;
    case _nB_PX_             :         return "nB*PX"            ;
    case _PfPb_avPf_avPb_ :            return "PfPb_avPf_avPb";
    case _PfPb_avPf_ :                 return "PfPb_avPf";
    case _PfPb_avPb_ :                 return "PfPb_avPb";
    case _PF_  :                       return "PF" ;
    case _PB_  :                       return "PB" ;
    case _PF_PB_  :                    return "PF*PB" ;
    case _nF_PB_  :                    return "nF*PB" ;
    case _nB_PF_  :                    return "nB*PF" ;
    case _nF_PF_ :                     return "nF*PF";
    case _PxPy_avPx_avPy_  :           return "PxPy_avPx_avPy" ;
    case _PxPy_avPx_  :                return "PxPy_avPx" ;
    case _PxPy_avPy_  :                return "PxPy_avPy" ;
    case _PX_  :                       return "PX" ;
    case _PY_  :                       return "PY" ;
    case _PX_PY_  :                    return "PX*PY" ;
    case _nX_PY_  :                    return "nX*PY" ;
    case _nX_PX_ :                     return "nX*PX";
    case _PF2_  :                      return "PF2" ;
    case _PX2_  :                      return "PX2" ;
    case _piF2_  :                     return "piF2" ;
    case _piX2_  :                     return "piX2" ;
    default: //throw std::invalid_argument("Unimplemented item");
    {
        cout << "AHTUNG! Unimplemented var name!" << endl;
        int aa;
        cin >> aa;
    }

    }

    return "";
}




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
    TH2D *hist_n[4];
    TH2D *hist_pt[4];
    TH2D *hist_w_wMinus1[4];
    TH2D *hist_w_wMinus1_pt[4];
    TH2D *hist_pt2[4];

    TH2D *QA_hist_n[4];
    TH2D *QA_hist_pt[4];


    //    int nWPs;
    TH3D *hAllWins; //!
    TH3D *hDetaDphi; //!
    TH3D *hAllEtaDphi; //!

    bool flagHistAllWins;
    bool flagHistDetaDphi;
    bool flagHistAllEtaDphi;


    TH2I *hAccMap; //!
//    TH3D *hSPEC_DetaDphi; //!
    //    THnD *hDeltaEta; //!
    TString strAnLevel;

    // some vars to make flexible hist filling (i.e. don't fill if bin name is not in varNames array):
    int currentSubsampleId;
    int currentWinId;
    int currentDetaDphiPairId;
    int currentEtaWinsDphiPairId;


    // ############
    WinPairWrapper()
    {
        hAllWins = 0x0;
        hDetaDphi = 0x0;
        hAllEtaDphi = 0x0;

        flagHistAllWins = true;
        flagHistDetaDphi = true;
        flagHistAllEtaDphi = true;


        hAccMap = 0x0;
    }

    void setup( const char* strPrefix, int cBin, int nSub,
                int *_pTypes, int *_pCharges, int _nEtaWins, int _nPhiWins, double *_etaRange, double *_ptRange, bool *whichHistosToTake
                , bool _fullAcceptanceForDenom = false //, double _etaForDenomMin = -0.8, double _etaForDenomMax = 0.8
            )
    {
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

        flagHistAllWins = whichHistosToTake[0];
        flagHistDetaDphi = whichHistosToTake[1];
        flagHistAllEtaDphi = whichHistosToTake[2];


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
//    }

    // ############
//    void setHistAllWins( const char* strPrefix, int cBin, int nSub ) //, const char* _varNames[] )//, int _nVars )
//    {
        int _nVars = _nVars_;
        cout << "_nVars = " << _nVars << endl;
        //        if ( nWPs == 0 )
        //        {
        //            cout << "AHTUNG!!! No win pairs!" << endl;
        //            int tmpA;
        //            cin >> tmpA;
        //        }

        strAnLevel = Form("%s", strPrefix);

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

        double etaStep  = ( etaMax - etaMin ) / nEtaWins;
        double min_dEta = -( etaMax - etaMin ) + etaStep;
        double max_dEta = ( etaMax - etaMin ) - etaStep;

        // find max dPhi:
        double phiStep  = TMath::TwoPi() / nPhiWins;
        double min_dPhi = -TMath::TwoPi() + phiStep;
        double max_dPhi = TMath::TwoPi() - phiStep;


        TString strHistName = Form("hAllWins_%s_cBin%d", strAnLevel.Data(), cBin);

        int nAllEtaWPs = nEtaWins*nEtaWins;
        int nAllPhiWPs = nPhiWins*nPhiWins;
        int nAllWPs = nAllEtaWPs * nAllPhiWPs;

        // all window pairs separately:
        if( flagHistAllWins )
            hAllWins = new TH3D( strHistName, strHistName
                                  , _nVars, -0.5, _nVars-0.5
                                  , nAllWPs, -0.5, nAllWPs-0.5
                                  , nSub, -0.5, nSub-0.5
                                  );

        int nDetaWP = 2*nEtaWins-1;
        int nDphiWP = 2*nPhiWins-1;

        TString str_dEta_dPhi_HistName = Form("hDetaDphi_%s_cBin%d", strAnLevel.Data(), cBin);
        if( flagHistDetaDphi )
            hDetaDphi = new TH3D( str_dEta_dPhi_HistName, str_dEta_dPhi_HistName
                              , _nVars, -0.5, _nVars-0.5
                              , nDetaWP*nDphiWP, -0.5, nDetaWP*nDphiWP-0.5
                              , nSub, -0.5, nSub-0.5
                              );
        cout << "nDetaWP = " << nDetaWP << ", min_dEta = " << min_dEta << ", max_dEta = " << max_dEta << endl;
        cout << "nDphiWP = " << nDphiWP << endl;


        int nAllEtaDphiWP = nEtaWins*nEtaWins * nDphiWP;
        TString str_AllEta_dPhi_HistName = Form("hAllEtaDphi_%s_cBin%d", strAnLevel.Data(), cBin);
        if( flagHistAllEtaDphi )
            hAllEtaDphi = new TH3D( str_AllEta_dPhi_HistName, str_AllEta_dPhi_HistName
                                          , _nVars, -0.5, _nVars-0.5
                                          , nAllEtaDphiWP, -0.5, nAllEtaDphiWP-0.5
                                          , nSub, -0.5, nSub-0.5
                                          );




        // set labels x and prepare a map with vars
        cout << "_nVars = " << _nVars << endl;
        for( int i = 0; i < _nVars; i++ )
        {
            const char* _varName = enumToStr(i);
            if(0)cout << "i = " << i << ", _varName = " << _varName << endl;

            if( flagHistAllWins )
                hAllWins->GetXaxis()->SetBinLabel( i+1, _varName ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
            if( flagHistDetaDphi )
                hDetaDphi->GetXaxis()->SetBinLabel( i+1, _varName ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
            if( flagHistAllEtaDphi )
                hAllEtaDphi->GetXaxis()->SetBinLabel( i+1, _varName ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
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
                            wpId++;
                        }
                }
        }


        // set labels y for hDetaDphi
        if( flagHistDetaDphi )
        {
            for( int i = 0; i < nDetaWP; i++ )
                for( int j = 0; j < nDphiWP; j++ )
                {
                    int sepId = nDphiWP*i + j;
                    hDetaDphi->GetYaxis()->SetBinLabel( sepId+1, Form("dEta_%.2f_dPhi_%.2f", min_dEta + etaStep*i, min_dPhi + phiStep*j ) );

                    if(0) cout << "i = " << i << ", j = " << j << ", sepId = " << sepId
                         << ", min_dEta + etaStep*i = " << min_dEta + etaStep*i
                         << ", min_dPhi + phiStep*j = " << min_dPhi + phiStep*j  << endl;
                }
        }


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
                                                                        , eFmin, eFmax,    min_dPhi + phiStep*j
                                                                        ) );
                    }
                }
            }
        }

    }
    

    // ############
    void addTrack( int pid, double eta, double phi, double pt, int charge, double weight = 1.0 )
    {
        if ( pt < ptWin[0] || pt > ptWin[1] )
            return;

        int pidAbs = abs(pid);

        // ##### backward window:
        // B
        if ( partTypes[1]==0 || pidAbs == partTypes[1] )
            if ( partCharges[1]==0 || charge == partCharges[1] )
            {
                hist_n[1]->Fill( eta, phi, weight );
                hist_pt[1]->Fill( eta, phi, pt*weight );

                hist_w_wMinus1[1]->Fill( eta, phi, weight*(weight-1) );
            }


        // Y
        if ( partTypes[3]==0 || pidAbs == partTypes[3] )
            if ( partCharges[3]==0 || charge == partCharges[3] )
            {
                hist_n[3]->Fill( eta, phi, weight );
                hist_pt[3]->Fill( eta, phi, pt*weight );

                hist_w_wMinus1[3]->Fill( eta, phi, weight*(weight-1) );
            }

        // ##### forward window:
        // F
        if ( partTypes[0]==0 || pidAbs == partTypes[0] )
            if ( partCharges[0]==0 || charge == partCharges[0] )
            {
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
                hist_n[2]->Fill( eta, phi, weight );
                hist_pt[2]->Fill( eta, phi, pt*weight );

                hist_w_wMinus1[2]->Fill( eta, phi, weight*(weight-1) );
                hist_w_wMinus1_pt[2]->Fill( eta, phi, weight*(weight-1)*pt );
                hist_pt2[2]->Fill( eta, phi, pt*pt *weight*weight );
            }
    }



    void fillHistWithValue( /*const char *varName,*/ int varId, double value ) //, bool forceFilling = false )
    {
        //                return;
        if( flagHistAllWins )
            hAllWins->Fill(  varId, currentWinId, currentSubsampleId, value  );

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

//        int sizeOfMap = map_DetaDphi_bin.size();
        //        cout << "in finishEvent: " << nEta << endl;



        //        // fill TH3D with single eta-phi win info
        //        for( int e1 = 0; e1 < nEtaWins; e1++ )
        //        {
        //            double center_e1 = hist_n[0]->GetXaxis()->GetBinCenter(e1+1);
        //            for( int p1 = 0; p1 < nPhiWins; p1++ )
        //            {
        //                int e1Den = !fullAcceptanceForDenom ? e1 : 0;
        //                int p1Den = !fullAcceptanceForDenom ? p1 : 0;


        //                double nF = hist_n[0]->GetBinContent(e1+1, p1+1);
        //                double nX = hist_n[2]->GetBinContent(e1Den+1, p1Den+1);

        //                double ptF = hist_pt[0]->GetBinContent(e1+1, p1+1);
        //                double ptX = hist_pt[2]->GetBinContent(e1Den+1, p1Den+1);


        //                double w_wMinus1_F = hist_w_wMinus1[0]->GetBinContent(e1+1, p1+1);
        //                double w_wMinus1_X = hist_w_wMinus1[2]->GetBinContent(e1Den+1, p1Den+1);

        //                double w_wMinus1_PF = hist_w_wMinus1_pt[0]->GetBinContent(e1+1, p1+1);
        //                double w_wMinus1_PX = hist_w_wMinus1_pt[2]->GetBinContent(e1Den+1, p1Den+1);

        //                double ptF2 = hist_pt2[0]->GetBinContent(e1+1, p1+1);
        //                double ptX2 = hist_pt2[2]->GetBinContent(e1Den+1, p1Den+1);
        //            }
        //        }

        // fill TH3D with win pairs info
        int wpId = 0;
        for( int e1 = 0; e1 < nEtaWins; e1++ )
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
                if( hAccMap && hAccMap->GetBinContent( e1+1, p1+1 ) )
                {
                    QA_hist_n[0]->Fill( center_e1, center_p1, nF );
                    QA_hist_pt[0]->Fill( center_e1, center_p1, ptF );
                }


                // loop over pairing windows:
                for( int e2 = 0; e2 < nEtaWins; e2++ )
                {
                    double center_e2 = hist_n[0]->GetXaxis()->GetBinCenter(e2+1);
                    double etaSep = center_e2 - center_e1;

                    for( int p2 = 0; p2 < nPhiWins; p2++ )
                    {
                        // check if we take this bin or not!
                        if( hAccMap && ( !hAccMap->GetBinContent( e1+1, p1+1 ) || !hAccMap->GetBinContent( e2+1, p2+1 ) ) )
                        {
                            wpId++;   // IMPORTANT!!!
                            continue;
                        }

                        double center_p2 = hist_n[0]->GetYaxis()->GetBinCenter(p2+1);
                        double phiSep = center_p2 - center_p1;
                        //                        continue;

                        //                        cout << "e1, p1, e2, p2 = " << e1 << ", " << p1 << ", " << e2 << ", " << p2 << endl;
                        currentWinId = wpId;

                        // winId for dEta-dPhi:
                        int id_dEta = nEtaWins-1 + e2-e1;
                        int id_dPhi = nPhiWins-1 + p2-p1;
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



                        fillHistWithValue( _Nevents_, 1 );//, true ); // last argument - forcing filling the 0th bin of the hist

                        // NfNb:
                        fillHistWithValue( _Nf_    ,   nF             );
                        fillHistWithValue( _Nb_    ,   nB             );
                        fillHistWithValue( _Nf2_  ,   nF * nF - w_wMinus1_F     );
                        fillHistWithValue( _Nb2_  ,   nB * nB - w_wMinus1_B     );
                        fillHistWithValue( _Nf_Nb_ ,   nF * nB      );

                        // NxNy:
                        fillHistWithValue( _Nx_      ,   nX             );
                        fillHistWithValue( _Ny_      ,   nY             );
                        fillHistWithValue( _Nx2_     ,   nX * nX - w_wMinus1_X   ); // !fullEtaForDenom ? _nX*_nX  : _nX*_nX - _w_wMinus1_X   );
                        fillHistWithValue( _Ny2_     ,   nY * nY - w_wMinus1_Y   ); // !fullEtaForDenom ? _nY*_nY  : _nY*_nY - _w_wMinus1_Y   );
                        fillHistWithValue( _Nx_Ny_   ,   nX * nY      );


                        // Nf_Nx Nf_Ny Nb_Nx Nb_Ny:
                        //                        if(0) fillHistWithValue( _Nf_Nx_  ,   nF * nX       );
                        fillHistWithValue( _Nf_Ny_  ,   nF * nY       );
                        fillHistWithValue( _Nb_Nx_  ,   nB * nX       );
                        //                        if(0) fillHistWithValue( _Nb_Ny_  ,   nB * nY       );


                        // FROM /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/AliForwardBackwardAnalysis.cxx:1583
                        //                        fillHistWithValue( _sumPtAllEvF_ , ptF );
                        //                        if(1) fillHistWithValue( _sumPtAllEvB_ , ptB );
                        //                        if(0) fillHistWithValue( _sumPtAllEvX_ , ptX );
                        //                        if(0) fillHistWithValue( _sumPtAllEvY_ , ptY );

                        // for dptdpt:
                        //                if(0) fillHistWithValue( _piFpjB_ , w->piFpjB );
                        fillHistWithValue( _nF_PB_ ,  nF * ptB );
                        fillHistWithValue( _nB_PF_ ,  nB * ptF );

                        // for C when av over pairs is OUTSIDE sum:
                        //                if(0) fillHistWithValue( _pipjF_ ,          w->pipjF );
                        //                        if(0) fillHistWithValue( _(nF-1)_PF_ ,      (nF-1)*ptF ); // (n-1)*sumPt
                        //                        if(0) fillHistWithValue( _nF_(nF-1)_ ,          (nF-1)*nF ); // (n-1)*n

                        //                if(0) fillHistWithValue( _pipjB_ ,         w->pipjB );
                        //                        if(0) fillHistWithValue( _(nB-1)_PB_ , (nB-1)*ptB ); // (n-1)*sumPt
                        //                        if(0) fillHistWithValue( _nB_(nB-1)_ ,     (nB-1)*nB ); // (n-1)*n



                        // to check VV comparison ptpt vs dptdpt: special terms:
                        fillHistWithValue( _PF_PB_ ,      ptF*ptB );
                        fillHistWithValue( _PF_ ,              ptF );
                        fillHistWithValue( _PB_ ,              ptB );

                        fillHistWithValue( _PX_PY_ ,      ptX*ptY );
                        fillHistWithValue( _PX_ ,              ptX );
                        fillHistWithValue( _PY_ ,              ptY );


                        // for new ratio-pt (April 2021)
                        //                        fillHistWithValue( _nY_PF_ , nY*ptF );
                        fillHistWithValue( _nB_PX_ , nB*ptX );
                        fillHistWithValue( _nY_PX_ , nY*ptX );

                        //                        fillHistWithValue( _nX_PB_ , nX*ptB );
                        //                        fillHistWithValue( _nF_PY_ , nF*ptY );
                        fillHistWithValue( _nX_PY_ , nX*ptY );

                        fillHistWithValue( _nF_PF_ , nF*ptF - w_wMinus1_PF );
                        fillHistWithValue( _nX_PX_ , nX*ptX - w_wMinus1_PX );

                        fillHistWithValue( _PF2_ , ptF*ptF ); //- w_wMinus1_PF );
                        fillHistWithValue( _PX2_ , ptX*ptX ); //- w_wMinus1_PX );

                        fillHistWithValue( _piF2_ , ptF2 );
                        fillHistWithValue( _piX2_ , ptX2 );

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
                            fillHistWithValue( _f_Nevents_ ,    1                  );

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
                            if(1) fillHistWithValue( _fb_Nevents_ ,   1                  );

                            if(1) fillHistWithValue( _PfPb_avPf_ ,       meanPtF             );
                            if(1) fillHistWithValue( _PfPb_avPb_ ,       meanPtB             );
                            //                            if(0) fillHistWithValue( _PfPb_avPf2_ ,     meanPtF*meanPtF      );
                            //                            if(0) fillHistWithValue( _PfPb_avPb2_ ,     meanPtB*meanPtB      );
                            if(1) fillHistWithValue( _PfPb_avPf_avPb_ ,   meanPtF*meanPtB      );

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
                            fillHistWithValue( _Nb_OVER_Ny_ ,    nB/(double)nY      );
                            //                            if(0) fillHistWithValue( _Nx_OVER_Ny_ ,    nX/(double)nY      );

                            //                            if(0) fillHistWithValue( _Nf_avPy_ ,   nF*meanPtY       );
                            //                            if(0) fillHistWithValue( _Nb_avPy_ ,   nB*meanPtY       );
                        }

                        // x:
                        if ( nX > 0 )
                        {
                            fillHistWithValue( _x_Nevents_ ,    1                  );

                            //                            if(0) fillHistWithValue( _PxNy_avPx_ ,         meanPtX           );
                            //                            if(0) fillHistWithValue( _PxNy_Ny_ ,        nY            );
                            //                            if(0) fillHistWithValue( _PxNy_avPx2_ ,     meanPtX*meanPtX      );
                            //                            if(0) fillHistWithValue( _PxNy_Ny2_ ,      nY*nY      );
                            //                            if(0) fillHistWithValue( _PxNy_avPx_Ny_ ,    meanPtX*nY      );

                            fillHistWithValue( _Nf_OVER_Nx_ ,    nF/(double)nX      );
                            //                            if(0) fillHistWithValue( _Nb_OVER_Nx_ ,    nB/(double)nX      );
                            //                            if(0) fillHistWithValue( _Ny_OVER_Nx_ ,    nY/(double)nX      );

                            //                            if(0) fillHistWithValue( _Nf_avPx_ ,   nF*meanPtX       );
                            //                            if(0) fillHistWithValue( _Nb_avPx_ ,   nB*meanPtX       );
                        }


                        // xy:
                        if ( nX > 0 && nY > 0 )
                        {
                            fillHistWithValue( _xy_Nevents_ ,   1                  );
                            if(1) fillHistWithValue( _PxPy_avPx_ ,       meanPtX             );
                            if(1) fillHistWithValue( _PxPy_avPy_ ,       meanPtY             );
                            //                            if(0) fillHistWithValue( _PxPy_avPx2_ ,     meanPtX*meanPtX      );
                            //                            if(0) fillHistWithValue( _PxPy_avPy2_ ,     meanPtY*meanPtY      );
                            if(1) fillHistWithValue( _PxPy_avPx_avPy_ ,   meanPtX*meanPtY      );

                            //                            if(0) fillHistWithValue( _Nf_OVER_Nx_vs_avPy_ ,   nF/(double)nX*meanPtY      );
                            if(1) fillHistWithValue( _Nb_OVER_Ny_vs_avPx_ ,   nB/(double)nY*meanPtX      );

                            fillHistWithValue( _Nf_OVER_Nx_vs_Nb_OVER_Ny_ ,   nF/(double)nX * nB/(double)nY      );
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
                            fillHistWithValue( _fy_Nevents_ ,   1                  );

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



};


#endif












