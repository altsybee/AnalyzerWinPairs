#ifndef WinPairFilling_cxx
#define WinPairFilling_cxx


#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include "TDirectory.h"
#include "TString.h"

#include "WinPairBase.h"

#include <iostream>
#include <map>
//#include <unordered_map>
#include <string>
//#include <bits/stdc++.h>


using namespace std;


//bool MAKE_TEST_FIRST_FILL_FOR_NON_EXISTING_VALUES = true;



// ######################################################################
// class which performs accumulation of quantities for a given event
// ######################################################################

class WinPairFilling : public WinPairBase //WinPair2019
{
public:
    // e-by-e data:
    //    double _nF; double _nB;
    int _int_counter_nF; int _int_counter_nB;
    double _nF; double _nB;
    double _ptF; double _ptB;

    //    double _nX; double _nY;
    double _nX; double _nY;
    double _ptX; double _ptY;

    double _ptF2; double _ptX2; // for avPt-avPt in same window


    double _w_wMinus1_F; double _w_wMinus1_B; // for corrections in one window!
    double _w_wMinus1_X; double _w_wMinus1_Y; // for corrections in one window!
    double _w_wMinus1_PF; // for corrections in one window! July 2022
    double _w_wMinus1_PX; // for corrections in one window! July 2022
//    double _w_wMinus1_PF2; // for corrections in one window! July 2022
//    double _w_wMinus1_PX2; // for corrections in one window! July 2022

    // extra info:
    Double_t pipjF; //, sum_piF;
    Double_t pipjB; //, sum_piB;
    // to calc C:
    double *arr_pF; //[FB_max_tracks_in_win];
    double *arr_pB; //[FB_max_tracks_in_win];
    int maxNtracks;


    // for dptdpt:
    Double_t piFpjB;


    WinPairFilling():
        _int_counter_nF(0), _int_counter_nB(0),
        _nF(0), _nB(0),
        _ptF(0), _ptB(0),
        _nX(0), _nY(0),
        _w_wMinus1_F(0), _w_wMinus1_B(0),
        _w_wMinus1_X(0), _w_wMinus1_Y(0),
        _w_wMinus1_PF(0), _w_wMinus1_PX(0), //_w_wMinus1_PF2(0), _w_wMinus1_PX2(0),
        _ptX(0), _ptY(0),
        _ptF2(0), _ptX2(0)
    {}

    void init(int _nMaxTracks) //setMaxNumberOfTracks(int _n)
    {
        arr_pF = new double[_nMaxTracks];
        arr_pB = new double[_nMaxTracks];
    }


    void addTrack( int pid, double eta, double pt, int charge, double weight = 1.0 )
    {
        int pidAbs = abs(pid);

        if ( pt < ptWin[0] || pt > ptWin[1] )
            return;

        // ##### backward window:
        // B
        if ( eta > eWin[0] && eta < eWin[1] )
            if ( partTypes[1]==0 || pidAbs == partTypes[1] )
                if ( partCharges[1]==0 || charge == partCharges[1] )
                {
                    //                    arr_pB[_nB] = pt;
                    //                    _nB++;
                    //                    _ptB += pt;
                    arr_pB[_int_counter_nB++] = pt *weight;
                    _nB += 1*weight;
                    _ptB += pt *weight;

                    _w_wMinus1_B += weight*(weight-1);
                }

        // Y
        if ( (!fullEtaForDenom && eta > eWin[0] && eta < eWin[1])
             || (fullEtaForDenom && eta > etaForDenomMin && eta < etaForDenomMax) )
            if ( partTypes[3]==0 || pidAbs == partTypes[3] )
                if ( partCharges[3]==0 || charge == partCharges[3] )
                {
                    //                    _nY++;
                    //                    _ptY += pt;
                    _nY += weight;
                    _ptY += pt *weight;

                    _w_wMinus1_Y += weight*(weight-1);
                }
        //        cout << "test12" << endl;

        // ##### forward window:
        // F
        if ( eta > eWin[2] && eta < eWin[3] )
            if ( partTypes[0]==0 || pidAbs == partTypes[0] )
                if ( partCharges[0]==0 || charge == partCharges[0] )
                {
                    //                    arr_pF[_nF] = pt;
                    //                    _nF++;
                    //                    _ptF += pt;
                    arr_pF[_int_counter_nF++] = pt *weight;
                    _nF += 1*weight ;
                    _ptF += pt *weight;

                    _w_wMinus1_F += weight*(weight-1);
                    _w_wMinus1_PF += weight*(weight-1)*pt;
//                    _w_wMinus1_PF2 += weight*(weight-1)*pt*pt;

                    _ptF2  += pt*pt *weight*weight;
                }
        // X
        if ( (!fullEtaForDenom && eta > eWin[2] && eta < eWin[3] )
             || (fullEtaForDenom && eta > etaForDenomMin && eta < etaForDenomMax) )
            if ( partTypes[2]==0 || pidAbs == partTypes[2] )
                if ( partCharges[2]==0 || charge == partCharges[2] )
                {
                    //                    _nX++;
                    //                    _ptX += pt;
                    _nX += 1*weight;
                    _ptX += pt*weight;

                    _w_wMinus1_X += weight*(weight-1);
                    _w_wMinus1_PX += weight*(weight-1)*pt;
//                    _w_wMinus1_PX2 += weight*(weight-1)*pt*pt;

                    _ptX2  += pt*pt *weight*weight;
                }

        //        cout << "test2" << endl;

    }

    void finishEvent()
    {
        // to calc C:
        if(0)  // 0 to speed up
        {
            // F
            for( Int_t i = 0; i < _int_counter_nF; i++ )
                for( Int_t j = i+1; j < _int_counter_nF; j++ )
                    pipjF += arr_pF[i]*arr_pF[j];

            // B
            for( Int_t i = 0; i < _int_counter_nB; i++ )
                for( Int_t j = i+1; j < _int_counter_nB; j++ )
                    pipjB += arr_pB[i]*arr_pB[j];

            // for dpt-dpt
            for( Int_t i = 0; i < _int_counter_nF; i++ )
                for( Int_t j = 0; j < _int_counter_nB; j++ )
                    piFpjB += arr_pF[i]*arr_pB[j];
        }
    }

    void resetEvent()
    {
        // reset values:
        _int_counter_nF = 0; _int_counter_nB = 0;
        _nF = 0; _nB = 0;
        _ptF = 0; _ptB = 0;

        _nX = 0; _nY = 0;
        _ptX = 0; _ptY = 0;

        _ptF2 = 0; _ptX2 = 0;

        _w_wMinus1_F = 0;
        _w_wMinus1_B = 0;
        _w_wMinus1_X = 0;
        _w_wMinus1_Y = 0;
        _w_wMinus1_PF = 0;
        _w_wMinus1_PX = 0;
//        _w_wMinus1_PF2 = 0;
//        _w_wMinus1_PX2 = 0;


        // extra info:
        pipjF = 0; //sum_piF = 0;
        pipjB = 0; //sum_piB = 0;

        piFpjB = 0;

    }

    //    ClassDef(WinPairFilling, 1);


};








const int MAX_N_WIN_PAIRS = 500;//150;

// ####################################################################################################################
// class which collects all window pairs and performs filling of TH3D histogram with e-by-e calculated values
// ####################################################################################################################

class WinPairWrapper
{
public:
    WinPairFilling *wp[MAX_N_WIN_PAIRS]; //! // win pairs
    double eSep[MAX_N_WIN_PAIRS];
    double eFsize[MAX_N_WIN_PAIRS];
    double eBsize[MAX_N_WIN_PAIRS];
    int nEta;
    TH3D *hAllWins; //!
    TH3D *hDeltaEta; //!
    TString strAnLevel;

    // some vars to make flexible hist filling (i.e. don't fill if bin name is not in varNames array):
    int currentSubsampleId;
    int currentWinId;
    double currentDeltaEta;

    map<const char*,int> mapVar;   // helps checking which vars are not requested in varNames array

    // ############
    WinPairWrapper()
    {
        nEta = 0;
    }

    // ############
    void addWinPair( //const char* strPrefix,
                     int *_pTypes, int *_pCharges
                     , double _eMinB, double _eMaxB, double _eMinF, double _eMaxF, double _ptMin, double _ptMax
                     , bool _fullEtaForDenom = false, double _etaForDenomMin = -0.8, double _etaForDenomMax = 0.8 )
    {
        if ( nEta >= MAX_N_WIN_PAIRS )
        {
            cout << "AHTUNG!!! nEta >= MAX_N_WIN_PAIRS" << endl;
            int aa;
            cin >> aa;
        }
        int iEta = nEta;
        WinPairFilling *_wp = new WinPairFilling;
        wp[iEta] = _wp;
        _wp->setParticleTypes( _pTypes, _pCharges );
        _wp->setWindows( _eMinB, _eMaxB, _eMinF, _eMaxF, _ptMin, _ptMax, _fullEtaForDenom, _etaForDenomMin, _etaForDenomMax );
        _wp->init( 1500 );

//        double eBsize = _eMaxB-_eMinB;
//        double eFsize = _eMaxF-_eMinF;
        double eBpos = ( _eMaxB + _eMinB )/2;
        double eFpos = ( _eMaxF + _eMinF )/2;

        eSep[nEta] = round( (eFpos - eBpos)*100 ) / 100;
        eFsize[nEta] = round( (_eMaxF - _eMinF)*100 ) / 100;
        eBsize[nEta] = round( (_eMaxB - _eMinB)*100 ) / 100;


        nEta++;
    }



    // ############
    void setHistAllWins( const char* strPrefix, int cBin, int nSub, const char* _varNames[], int _nVars )
    {
        if ( nEta == 0 )
        {
            cout << "AHTUNG!!! No win pairs!" << endl;
            int tmpA;
            cin >> tmpA;
        }

        strAnLevel = Form("%s", strPrefix);

        int nWins = nEta;

        TString strHistName = Form("hAllWins_%s_cBin%d", strAnLevel.Data(), cBin);


        // all window pairs separately:
        hAllWins = new TH3D( strHistName, strHistName
                             //, nBinsX, -0.5, nBinsX-0.5
                             , _nVars, -0.5, _nVars-0.5
                             , nWins, -0.5, nWins-0.5
                             , nSub, -0.5, nSub-0.5
                             );

        //
        // count unique dEta-s:
        map <double, int> _mapUnique_dEta;
        for( int i=0; i < nEta; i++)
            _mapUnique_dEta[ eSep[i] ]++;
        int nUnique_dEta = _mapUnique_dEta.size();

        // QA:
        if(0)
        {
            for( int i=0; i<nEta; i++)
                cout << "i = " << i << ", etaSep = " << eSep[i] << ": nWins = " << _mapUnique_dEta[ eSep[i] ] << endl;

            int aa;
            cin >> aa;
        }

        // find max dEta:
        double min_dEta = *min_element(eSep, eSep+nEta);
        double max_dEta = *max_element(eSep, eSep+nEta);
        TString str_dEta_HistName = Form("hDeltaEta_%s_cBin%d", strAnLevel.Data(), cBin);
        hDeltaEta = new TH3D( str_dEta_HistName, str_dEta_HistName
                             , _nVars, -0.5, _nVars-0.5
                             , nUnique_dEta, min_dEta - eFsize[0]/2, max_dEta + eFsize[0]/2
                             , nSub, -0.5, nSub-0.5
                             );
        cout << "nUnique_dEta = " << nUnique_dEta << ", min_dEta = " << min_dEta << ", max_dEta = " << max_dEta << endl;


        // set labels x and prepare a map with vars
        for( int i = 0; i < _nVars; i++ )
        {
            mapVar.insert( pair<const char*,int>( _varNames[i], i ) );
            hAllWins->GetXaxis()->SetBinLabel( i+1, _varNames[i] ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
            hDeltaEta->GetXaxis()->SetBinLabel( i+1, _varNames[i] ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
        }
//        cout << "check mapVar size: " << mapVar.size() << endl;
        // check for duplicate vars:
        if ( (int)mapVar.size() < _nVars )
        {
            cout << "AHTUNG!!! there are duplicates in the varNames array!!!" << endl;
            int aa;
            cin >> aa;
        }

        // set labels y for hAllWins
        for( int winId = 0; winId < nWins; winId++ )
            hAllWins->GetYaxis()->SetBinLabel( winId+1, Form("eB_%.2f_%.2f_eF_%.2f_%.2f", wp[winId]->eWin[0], wp[winId]->eWin[1], wp[winId]->eWin[2], wp[winId]->eWin[3] ) );
    }


    // ############
    void addTrack( int pid, double eta, double pt, int charge, double weight = 1.0 )
    {
        //        return;
        for( int eWin = 0; eWin < nEta; eWin++ )
            wp[eWin]->addTrack( pid, eta, /*phi,*/ pt, charge, weight );
    }



    void fillHistWithValue( const char *varName, double value, bool forceFilling = false )
    {
//        return;
        int varId = mapVar[varName];
        if( varId != 0 || forceFilling ) // because if varName is not in the varNames array, map will add a new key-value pair with value = 0.
        {
            hAllWins->Fill(  varId, currentWinId, currentSubsampleId, value  );
            hDeltaEta->Fill(  varId, currentDeltaEta, currentSubsampleId, value  );
        }
        // here we are forcing filling "Nevents" (varId==0) by hand
    }


    // ############
    void finishEvent( int subId )
    {
        currentSubsampleId = subId;

//        cout << "in finishEvent: " << nEta << endl;
        for( int wpId = 0; wpId < nEta; wpId++ )
        {
            currentWinId = wpId;
            currentDeltaEta = eSep[wpId];
//            varCounter = 0;

            //
            WinPairFilling *w = wp[wpId];
            w->finishEvent();


            double meanPtF = -1;
            double meanPtB = -1;
            double meanPtX = -1;
            double meanPtY = -1;

            if ( w->_nF > 0 ) { meanPtF = w->_ptF / w->_nF; }
            if ( w->_nB > 0 ) { meanPtB = w->_ptB / w->_nB; }
            if ( w->_nX > 0 ) { meanPtX = w->_ptX / w->_nX; }
            if ( w->_nY > 0 ) { meanPtY = w->_ptY / w->_nY; }



            fillHistWithValue( "Nevents", 1, true ); // last argument - forcing filling the 0th bin of the hist

            // NfNb:
            fillHistWithValue( "Nf"    ,   w->_nF             );
            fillHistWithValue( "Nb"    ,   w->_nB             );
            fillHistWithValue( "Nf2"  ,   w->_nF * w->_nF - w->_w_wMinus1_F     );
            fillHistWithValue( "Nb2"  ,   w->_nB * w->_nB - w->_w_wMinus1_B     );
            fillHistWithValue( "Nf*Nb" ,   w->_nF * w->_nB      );

            // NxNy:
            fillHistWithValue( "Nx"      ,   w->_nX             );
            fillHistWithValue( "Ny"      ,   w->_nY             );
            fillHistWithValue( "Nx2"     ,   w->_nX * w->_nX - w->_w_wMinus1_X   ); // !fullEtaForDenom ? _nX*_nX  : _nX*_nX - _w_wMinus1_X   );
            fillHistWithValue( "Ny2"     ,   w->_nY * w->_nY - w->_w_wMinus1_Y   ); // !fullEtaForDenom ? _nY*_nY  : _nY*_nY - _w_wMinus1_Y   );
            fillHistWithValue( "Nx*Ny"   ,   w->_nX * w->_nY      );


            // Nf_Nx Nf_Ny Nb_Nx Nb_Ny:
            if(0) fillHistWithValue( "Nf*Nx"  ,   w->_nF * w->_nX       );
            fillHistWithValue( "Nf*Ny"  ,   w->_nF * w->_nY       );
            fillHistWithValue( "Nb*Nx"  ,   w->_nB * w->_nX       );
            if(0) fillHistWithValue( "Nb*Ny"  ,   w->_nB * w->_nY       );


            // FROM /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/AliForwardBackwardAnalysis.cxx:1583
            fillHistWithValue( "sumPtAllEvF" , w->_ptF );
            if(1) fillHistWithValue( "sumPtAllEvB" , w->_ptB );
            if(0) fillHistWithValue( "sumPtAllEvX" , w->_ptX );
            if(0) fillHistWithValue( "sumPtAllEvY" , w->_ptY );

            // for dptdpt:
            if(0) fillHistWithValue( "piFpjB" , w->piFpjB );
            fillHistWithValue( "nF*PB" ,  w->_nF * w->_ptB );
            fillHistWithValue( "nB*PF" ,  w->_nB * w->_ptF );

            // for C when av over pairs is OUTSIDE sum:
            if(0) fillHistWithValue( "pipjF" ,          w->pipjF );
            if(0) fillHistWithValue( "(nF-1)*PF" ,      (w->_nF-1)*w->_ptF ); // (n-1)*sumPt
            if(0) fillHistWithValue( "nF*(nF-1)" ,          (w->_nF-1)*w->_nF ); // (n-1)*n

            if(0) fillHistWithValue( "pipjB" ,         w->pipjB );
            if(0) fillHistWithValue( "(nB-1)*PB" , (w->_nB-1)*w->_ptB ); // (n-1)*sumPt
            if(0) fillHistWithValue( "nB*(nB-1)" ,     (w->_nB-1)*w->_nB ); // (n-1)*n



            // to check VV comparison ptpt vs dptdpt: special terms:
            fillHistWithValue( "PF*PB" ,      w->_ptF*w->_ptB );
            fillHistWithValue( "PF" ,              w->_ptF );
            fillHistWithValue( "PB" ,              w->_ptB );

            fillHistWithValue( "PX*PY" ,      w->_ptX*w->_ptY );
            fillHistWithValue( "PX" ,              w->_ptX );
            fillHistWithValue( "PY" ,              w->_ptY );


            // for new ratio-pt (April 2021)
            fillHistWithValue( "nY*PF" , w->_nY*w->_ptF );
            fillHistWithValue( "nB*PX" , w->_nB*w->_ptX );
            fillHistWithValue( "nY*PX" , w->_nY*w->_ptX );

            fillHistWithValue( "nX*PB" , w->_nX*w->_ptB );
            fillHistWithValue( "nF*PY" , w->_nF*w->_ptY );
            fillHistWithValue( "nX*PY" , w->_nX*w->_ptY );

            fillHistWithValue( "nF*PF" , w->_nF*w->_ptF - w->_w_wMinus1_PF );
            fillHistWithValue( "nX*PX" , w->_nX*w->_ptX - w->_w_wMinus1_PX );

            fillHistWithValue( "PF2" , w->_ptF*w->_ptF ); //- w->_w_wMinus1_PF );
            fillHistWithValue( "PX2" , w->_ptX*w->_ptX ); //- w->_w_wMinus1_PX );

            fillHistWithValue( "piF2" , w->_ptF2 );
            fillHistWithValue( "piX2" , w->_ptX2 );

            if ( w->_nY > 0 )
            {
                if(0) fillHistWithValue( "Nb_OVER_Ny*PF" , w->_nB/(double)w->_nY * w->_ptF );
                if(0) fillHistWithValue( "Nb_OVER_Ny*PX" , w->_nB/(double)w->_nY * w->_ptX );
            }
            if ( w->_nX > 0 )
            {
                if(0) fillHistWithValue( "Nf_OVER_Nx*PB" , w->_nF/(double)w->_nX * w->_ptB );
                if(0) fillHistWithValue( "Nf_OVER_Nx*PY" , w->_nF/(double)w->_nX * w->_ptY );
            }



            // NfPb:
            if ( w->_nB > 0 )
            {
                if(0) fillHistWithValue( "b_Nevents" ,    1                  );

                if(0) fillHistWithValue( "NfPb_Nf" ,       w->_nF             );
                if(0) fillHistWithValue( "NfPb_avPb" ,       meanPtB             );
                if(0) fillHistWithValue( "NfPb_Nf2" ,      w->_nF*w->_nF      );
                if(0) fillHistWithValue( "NfPb_avPb2" ,      meanPtB*meanPtB      );
                if(0) fillHistWithValue( "NfPb_Nf_avPb" ,    w->_nF*meanPtB      );

                if(0) fillHistWithValue( "Nf_OVER_Nb" ,    w->_nF/(double)w->_nB      );
                if(0) fillHistWithValue( "Nx_OVER_Nb" ,    w->_nX/(double)w->_nB      );
                if(0) fillHistWithValue( "Ny_OVER_Nb" ,    w->_nY/(double)w->_nB      );

                if(0) fillHistWithValue( "avPb_Nx" ,   meanPtB*w->_nX      );
                if(0) fillHistWithValue( "avPb_Ny" ,   meanPtB*w->_nY      );


                //            // for C when av over pairs is OUTSIDE sum:
                //            fillHistWithValue( "pipjB" ,         pipjB );
                //            fillHistWithValue( "(nB-1)*sum_pB" , (_nB-1)*(_nB*_ptB) ); // (n-1)*sumPt
                //            fillHistWithValue( "nB*(nB-1)" ,     (_nB-1)*_nB ); // (n-1)*n

                // for C when av over pairs is INSIDE sum: (like in GOOD_Ck_definition_STAR_2005_0504031.pdf)
                int nPairsB = (w->_nB-1)*w->_nB;
                if ( nPairsB > 0 )
                {
                    //                fillHistWithValue( "pipjB_avPerEv" ,         pipjB / nPairsB );
                    //                fillHistWithValue( "(nB-1)*PB_avPerEv" , (_nB-1)*_ptB / nPairsB );
                }

            }

            // PfNb:
            if ( w->_nF > 0 )
            {
                fillHistWithValue( "f_Nevents" ,    1                  );

                fillHistWithValue( "PfNb_avPf" ,         meanPtF           );
                if(0) fillHistWithValue( "PfNb_Nb" ,        w->_nB            );
                if(0) fillHistWithValue( "PfNb_avPf2" ,     meanPtF*meanPtF      );
                if(0) fillHistWithValue( "PfNb_Nb2" ,      w->_nB*w->_nB      );
                if(0) fillHistWithValue( "PfNb_avPf_Nb" ,    meanPtF*w->_nB      );

                if(0) fillHistWithValue( "Nb_OVER_Nf" ,    w->_nB/(double)w->_nF      );
                if(0) fillHistWithValue( "Nx_OVER_Nf" ,    w->_nX/(double)w->_nF      );
                if(0) fillHistWithValue( "Ny_OVER_Nf" ,    w->_nY/(double)w->_nF      );

                if(0) fillHistWithValue( "avPf_Nx" ,   meanPtF*w->_nX      );
                if(0) fillHistWithValue( "avPf_Ny" ,   meanPtF*w->_nY      );

                //            // for C when av over pairs is OUTSIDE sum:
                //            fillHistWithValue( "pipjF" ,              pipjF );
                //            fillHistWithValue( "(nF-1)*sum_pF" ,      (_nF-1)*(_nF*_ptF) ); // (n-1)*sumPt
                //            fillHistWithValue( "nF*(nF-1)" ,          (_nF-1)*_nF ); // (n-1)*n

                // for C when av over pairs is INSIDE sum:
                int nPairsF = (w->_nF-1)*w->_nF;
                if ( nPairsF > 0 )
                {
                    //                fillHistWithValue( "pipjF_avPerEv" ,          pipjF / nPairsF );
                    //                fillHistWithValue( "(nF-1)*PF_avPerEv" ,  (_nF-1)*_ptF / nPairsF );
                }

            }

            // fb:
            if ( w->_nF > 0 && w->_nB > 0 )
            {
                if(1) fillHistWithValue( "fb_Nevents" ,   1                  );

                if(1) fillHistWithValue( "PfPb_avPf" ,       meanPtF             );
                if(1) fillHistWithValue( "PfPb_avPb" ,       meanPtB             );
                if(0) fillHistWithValue( "PfPb_avPf2" ,     meanPtF*meanPtF      );
                if(0) fillHistWithValue( "PfPb_avPb2" ,     meanPtB*meanPtB      );
                if(1) fillHistWithValue( "PfPb_avPf_avPb" ,   meanPtF*meanPtB      );

                if(0) fillHistWithValue( "Nx_OVER_Nf_vs_avPb" ,   w->_nX/(double)w->_nF *meanPtB      );
                if(0) fillHistWithValue( "Ny_OVER_Nb_vs_avPf" ,   w->_nY/(double)w->_nB *meanPtF      );

                if(0) fillHistWithValue( "Nx_OVER_Nf_vs_Ny_OVER_Nb" ,   w->_nX/(double)w->_nF * w->_nY/(double)w->_nB      );

                //            // for dptdpt:
                //            fillHistWithValue( "piFpjB" ,     piFpjB );
                //            fillHistWithValue( "nF*sum_pB" ,  _nF*(_nB*_ptB) );
                //            fillHistWithValue( "nB*sum_pF" ,  _nB*(_nF*_ptF) );

            }


            // y:
            if ( w->_nY > 0 )
            {
                if(0) fillHistWithValue( "y_Nevents" ,    1                  );

                if(0) fillHistWithValue( "NxPy_Nx" ,        w->_nX             );
                if(0) fillHistWithValue( "NxPy_avPy" ,       meanPtY             );
                if(0) fillHistWithValue( "NxPy_Nx2" ,      w->_nX*w->_nX      );
                if(0) fillHistWithValue( "NxPy_avPy2" ,      meanPtY*meanPtY      );
                if(0) fillHistWithValue( "NxPy_Nx_avPy" ,    w->_nX*meanPtY      );


                if(0) fillHistWithValue( "Nf_OVER_Ny" ,    w->_nF/(double)w->_nY      );
                fillHistWithValue( "Nb_OVER_Ny" ,    w->_nB/(double)w->_nY      );
                if(0) fillHistWithValue( "Nx_OVER_Ny" ,    w->_nX/(double)w->_nY      );

                if(0) fillHistWithValue( "Nf_avPy" ,   w->_nF*meanPtY       );
                if(0) fillHistWithValue( "Nb_avPy" ,   w->_nB*meanPtY       );
            }

            // x:
            if ( w->_nX > 0 )
            {
                fillHistWithValue( "x_Nevents" ,    1                  );

                if(0) fillHistWithValue( "PxNy_avPx" ,         meanPtX           );
                if(0) fillHistWithValue( "PxNy_Ny" ,        w->_nY            );
                if(0) fillHistWithValue( "PxNy_avPx2" ,     meanPtX*meanPtX      );
                if(0) fillHistWithValue( "PxNy_Ny2" ,      w->_nY*w->_nY      );
                if(0) fillHistWithValue( "PxNy_avPx_Ny" ,    meanPtX*w->_nY      );

                fillHistWithValue( "Nf_OVER_Nx" ,    w->_nF/(double)w->_nX      );
                if(0) fillHistWithValue( "Nb_OVER_Nx" ,    w->_nB/(double)w->_nX      );
                if(0) fillHistWithValue( "Ny_OVER_Nx" ,    w->_nY/(double)w->_nX      );

                if(0) fillHistWithValue( "Nf_avPx" ,   w->_nF*meanPtX       );
                if(0) fillHistWithValue( "Nb_avPx" ,   w->_nB*meanPtX       );
            }


            // xy:
            if ( w->_nX > 0 && w->_nY > 0 )
            {
                fillHistWithValue( "xy_Nevents" ,   1                  );
                if(1) fillHistWithValue( "PxPy_avPx" ,       meanPtX             );
                if(1) fillHistWithValue( "PxPy_avPy" ,       meanPtY             );
                if(0) fillHistWithValue( "PxPy_avPx2" ,     meanPtX*meanPtX      );
                if(0) fillHistWithValue( "PxPy_avPy2" ,     meanPtY*meanPtY      );
                if(1) fillHistWithValue( "PxPy_avPx_avPy" ,   meanPtX*meanPtY      );

                if(0) fillHistWithValue( "Nf_OVER_Nx_vs_avPy" ,   w->_nF/(double)w->_nX*meanPtY      );
                if(1) fillHistWithValue( "Nb_OVER_Ny_vs_avPx" ,   w->_nB/(double)w->_nY*meanPtX      );

                fillHistWithValue( "Nf_OVER_Nx_vs_Nb_OVER_Ny" ,   w->_nF/(double)w->_nX * w->_nB/(double)w->_nY      );
            }

            // fx:
            if ( w->_nF > 0 && w->_nX > 0 )
            {
                if(0) fillHistWithValue( "fx_Nevents" ,   1                  );

                if(0) fillHistWithValue( "avPf_avPx" ,   meanPtF*meanPtX      );

                if(0) fillHistWithValue( "Nb_OVER_Nf_vs_Ny_OVER_Nx" ,   w->_nB/(double)w->_nF * w->_nY/(double)w->_nX      );
            }
            // fy:
            if ( w->_nF > 0 && w->_nY > 0 )
            {
                fillHistWithValue( "fy_Nevents" ,   1                  );

                if(0) fillHistWithValue( "avPf_avPy" ,   meanPtF*meanPtY      );
                if(0) fillHistWithValue( "Nx_OVER_Nf_vs_avPy" ,   w->_nX/(double)w->_nF *meanPtY      );
                if(0) fillHistWithValue( "Nb_OVER_Ny_vs_avPf" ,   w->_nB/(double)w->_nY *meanPtF      );
            }
            // bx:
            if ( w->_nB > 0 && w->_nX > 0 )
            {
                if(0) fillHistWithValue( "bx_Nevents" ,   1                  );

                if(0) fillHistWithValue( "Pb_Px" ,   meanPtB*meanPtX      );
                if(0) fillHistWithValue( "Nf_OVER_Nx_vs_avPb" ,   w->_nF/(double)w->_nX *meanPtB      );
                if(0) fillHistWithValue( "Ny_OVER_Nb_vs_avPx" ,   w->_nY/(double)w->_nB *meanPtX      );
            }
            // by:
            if ( w->_nB > 0 && w->_nY > 0 )
            {
                if(0) fillHistWithValue( "by_Nevents" ,   1                  );
                if(0) fillHistWithValue( "avPb_avPy" ,   meanPtB*meanPtY      );

                if(0) fillHistWithValue( "Nf_OVER_Nb_vs_Nx_OVER_Ny" ,   w->_nF/(double)w->_nB * w->_nX/(double)w->_nY      );
            }


            w->resetEvent();
        } // end of loop over eta wins

    }   // end of finishEvent( int subId )



};


#endif












