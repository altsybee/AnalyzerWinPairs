#ifndef WinPairFilling_cxx
#define WinPairFilling_cxx


#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THn.h"

#include "TDirectory.h"
#include "TString.h"

#include "WinPairBase.h"

#include <iostream>
#include <map>
//#include <unordered_map>
#include <array>
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


    void addTrack( int pid, double eta, double phi, double pt, int charge, double weight = 1.0 )
    {
        int pidAbs = abs(pid);

        if ( pt < ptWin[0] || pt > ptWin[1] )
            return;

        // ##### backward window:
        // B
        if ( eta > eWin[0] && eta < eWin[1] && phi > phiWin[0] && phi < phiWin[1] )
            if ( partTypes[1]==0 || pidAbs == partTypes[1] )
                if ( partCharges[1]==0 || charge == partCharges[1] )
                {
                    //                    arr_pB[_nB] = pt;
                    //                    _nB++;
                    //                    _ptB += pt;
//  !!!! commented to speed up!   arr_pB[_int_counter_nB++] = pt *weight;
                    _nB += 1*weight;
                    _ptB += pt *weight;

                    _w_wMinus1_B += weight*(weight-1);
                }

        // Y
        if ( (!fullEtaForDenom && eta > eWin[0] && eta < eWin[1] && phi > phiWin[0] && phi < phiWin[1] )
             || (fullEtaForDenom && eta > etaForDenomMin && eta < etaForDenomMax ) )  // and full phi (2pi)
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
        if ( eta > eWin[2] && eta < eWin[3] && phi > phiWin[2] && phi < phiWin[3]  )
            if ( partTypes[0]==0 || pidAbs == partTypes[0] )
                if ( partCharges[0]==0 || charge == partCharges[0] )
                {
                    //                    arr_pF[_nF] = pt;
                    //                    _nF++;
                    //                    _ptF += pt;
     //  !!!! commented to speed up!                arr_pF[_int_counter_nF++] = pt *weight;
                    _nF += 1*weight ;
                    _ptF += pt *weight;

                    _w_wMinus1_F += weight*(weight-1);
                    _w_wMinus1_PF += weight*(weight-1)*pt;
//                    _w_wMinus1_PF2 += weight*(weight-1)*pt*pt;

                    _ptF2  += pt*pt *weight*weight;
                }
        // X
        if ( (!fullEtaForDenom && eta > eWin[2] && eta < eWin[3] && phi > phiWin[2] && phi < phiWin[3] )
             || (fullEtaForDenom && eta > etaForDenomMin && eta < etaForDenomMax) )  // and full phi (2pi)
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








const int MAX_N_WIN_PAIRS = 4096;//150;

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
    double phiSep[MAX_N_WIN_PAIRS];
    double phiFsize[MAX_N_WIN_PAIRS];
    double phiBsize[MAX_N_WIN_PAIRS];
    int nUnique_dEta;
    int nUnique_dPhi;

    int nWPs;
    TH3D *hAllWins; //!
    TH3D *hDetaDphi; //!
//    THnD *hDeltaEta; //!
    TString strAnLevel;

    // some vars to make flexible hist filling (i.e. don't fill if bin name is not in varNames array):
    int currentSubsampleId;
    int currentWinId;
//    double currentDeltaEta;
//    double currentDeltaPhi;
    int currentDetaDphiPairId;

    map<const char*,int> mapVar;   // helps checking which vars are not requested in varNames array
//    map<double*, int> map_DetaDphi_bin;  // correspondance between dEta-dPhi bins and binId in TH3D (y axis)
    map< array<int, 2>, int> map_DetaDphi_bin;  // correspondance between dEta-dPhi bins and binId in TH3D (y axis)
//    map< int, int> map_DetaDphi_bin;  // correspondance between dEta-dPhi bins and binId in TH3D (y axis)

    // ############
    WinPairWrapper()
    {
        nWPs = 0;
    }

    // ############
    void addWinPair( //const char* strPrefix,
                     int *_pTypes, int *_pCharges
                     , double *_eWins  //double _eMinB, double _eMaxB, double _eMinF, double _eMaxF
                     , double *_phiWins
                     , double _ptMin, double _ptMax
                     , bool _fullEtaForDenom = false, double _etaForDenomMin = -0.8, double _etaForDenomMax = 0.8 )
    {
        if ( nWPs >= MAX_N_WIN_PAIRS )
        {
            cout << "AHTUNG!!! nWPs >= MAX_N_WIN_PAIRS" << endl;
            int aa;
            cin >> aa;
        }
        int iWin = nWPs;
        WinPairFilling *_wp = new WinPairFilling;
        wp[iWin] = _wp;
        _wp->setParticleTypes( _pTypes, _pCharges );
        _wp->setWindows( /*_eMinB, _eMaxB, _eMinF, _eMaxF,*/_eWins, _phiWins, _ptMin, _ptMax, _fullEtaForDenom, _etaForDenomMin, _etaForDenomMax );
        _wp->init( 1500 );

//        double eBsize = _eMaxB-_eMinB;
//        double eFsize = _eMaxF-_eMinF;
        double eBpos = ( _eWins[1] + _eWins[0] )/2;
        double eFpos = ( _eWins[3] + _eWins[2] )/2;

        eSep[iWin] = round( (eFpos - eBpos)*100 ) / 100;
        eFsize[iWin] = round( (_eWins[3] - _eWins[2])*100 ) / 100;
        eBsize[iWin] = round( (_eWins[1] - _eWins[0])*100 ) / 100;

        double phiBpos = ( _phiWins[1] + _phiWins[0] )/2;
        double phiFpos = ( _phiWins[3] + _phiWins[2] )/2;
        phiSep[iWin] = round( (phiFpos - phiBpos)*100 ) / 100;
        phiFsize[iWin] = round( (_phiWins[3] - _phiWins[2])*100 ) / 100;
        phiBsize[iWin] = round( (_phiWins[3] - _phiWins[2])*100 ) / 100;

        nWPs++;
    }



    // ############
    void setHistAllWins( const char* strPrefix, int cBin, int nSub, const char* _varNames[], int _nVars )
    {
        if ( nWPs == 0 )
        {
            cout << "AHTUNG!!! No win pairs!" << endl;
            int tmpA;
            cin >> tmpA;
        }

        strAnLevel = Form("%s", strPrefix);

//        int nWPs = nEta;

        TString strHistName = Form("hAllWins_%s_cBin%d", strAnLevel.Data(), cBin);


        // all window pairs separately:
        hAllWins = new TH3D( strHistName, strHistName
                             //, nBinsX, -0.5, nBinsX-0.5
                             , _nVars, -0.5, _nVars-0.5
                             , nWPs, -0.5, nWPs-0.5
                             , nSub, -0.5, nSub-0.5
                             );

        //
        // count unique dEta-s:
        map <double, int> _mapUnique_dEta;
        for( int i=0; i < nWPs; i++)
            _mapUnique_dEta[ eSep[i] ]++;
        nUnique_dEta = _mapUnique_dEta.size();

        // count unique dEta & dPhi-s:
        map <double, int> _mapUnique_dPhi;
        for( int i=0; i < nWPs; i++)
            _mapUnique_dPhi[ phiSep[i] ]++;
        nUnique_dPhi = _mapUnique_dPhi.size();

        // QA:
        if(0)
        {
            for( int i=0; i<nWPs; i++)
                cout << "i = " << i << ", etaSep = " << eSep[i] << ": nWPs = " << _mapUnique_dEta[ eSep[i] ] << endl;

            for( int i=0; i<nWPs; i++)
                cout << "i = " << i << ", phiSep = " << phiSep[i] << ": nWPs = " << _mapUnique_dPhi[ phiSep[i] ] << endl;

            int aa;
            cin >> aa;
        }
        // find max dEta:
        double min_dEta = *min_element(eSep, eSep+nWPs);
        double max_dEta = *max_element(eSep, eSep+nWPs);
        double etaStep = ( max_dEta-min_dEta + eFsize[0] ) / nUnique_dEta;

        // find max dPhi:
        double min_dPhi = *min_element(phiSep, phiSep+nWPs);
        double max_dPhi = *max_element(phiSep, phiSep+nWPs);
        double phiStep = ( max_dPhi-min_dPhi + phiFsize[0] ) / nUnique_dPhi;



        TString str_dEta_dPhi_HistName = Form("hDetaDphi_%s_cBin%d", strAnLevel.Data(), cBin);
        hDetaDphi = new TH3D( str_dEta_dPhi_HistName, str_dEta_dPhi_HistName
                             , _nVars, -0.5, _nVars-0.5
//                             , nUnique_dEta, min_dEta - eFsize[0]/2, max_dEta + eFsize[0]/2
                             , nUnique_dEta*nUnique_dPhi, -0.5, nUnique_dEta*nUnique_dPhi-0.5
                             , nSub, -0.5, nSub-0.5
                             );

//        int arr_nBins[] = { _nVars, nUnique_dEta, nUnique_dPhi, nSub };
//        double arr_mins[] = { -0.5,        min_dEta - eFsize[0]/2,     min_dPhi - phiFsize[0]/2, -0.5 };
//        double arr_maxs[] = { _nVars-0.5,  max_dEta + eFsize[0]/2,     max_dPhi + phiFsize[0]/2, nSub-0.5 };

//        hDeltaEta = new THnD( str_dEta_HistName, str_dEta_HistName, 4
//                             , arr_nBins
//                             , arr_mins
//                             , arr_maxs
//                             );


        cout << "nUnique_dEta = " << nUnique_dEta << ", min_dEta = " << min_dEta << ", max_dEta = " << max_dEta << endl;
        cout << "nUnique_dPhi = " << nUnique_dPhi << endl;


        // set labels x and prepare a map with vars
        for( int i = 0; i < _nVars; i++ )
        {
            mapVar.insert( pair<const char*,int>( _varNames[i], i ) );
            hAllWins->GetXaxis()->SetBinLabel( i+1, _varNames[i] ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
            hDetaDphi->GetXaxis()->SetBinLabel( i+1, _varNames[i] ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
//            hDeltaEta->GetAxis(0)->SetBinLabel( i+1, _varNames[i] ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
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
        for( int winId = 0; winId < nWPs; winId++ )
            hAllWins->GetYaxis()->SetBinLabel( winId+1, Form("etaB_%.2f_%.2f_phiB_%.2f_%.2f_etaF_%.2f_%.2f_phiF_%.2f_%.2f"
                    , wp[winId]->eWin[0], wp[winId]->eWin[1],    wp[winId]->phiWin[0], wp[winId]->phiWin[1]
                    , wp[winId]->eWin[2], wp[winId]->eWin[3],    wp[winId]->phiWin[2], wp[winId]->phiWin[3]
                    ) );

        // set labels y for hDetaDphi
        int kWin = 0;
        for( int i = 0; i < nUnique_dEta; i++ )
            for( int j = 0; j < nUnique_dPhi; j++ )
            {
//                double pair_DetaDphi[] = { etaStep*i, phiStep*j };
//                int pair_DetaDphi[] = { (int)(etaStep*i *100), (int)(phiStep*j *100) }; // *100 to convert to int - to be stable wrt to double keys!
                array<int, 2> pair_DetaDphi = { (int)round(  (min_dEta + etaStep*i)*100 ), (int)round( (min_dPhi + phiStep*j) *100 ) };
//                int pair_DetaDphi = (int)(  (min_dEta + etaStep*i) * 100) *10000 + (int)(  (min_dPhi + phiStep*j) * 100);
                map_DetaDphi_bin[ pair_DetaDphi ] = kWin;
//                cout << "etaStep*i = " << min_dEta + etaStep*i << ", phiStep*j = " << min_dPhi + phiStep*j << ", pair_DetaDphi = " << pair_DetaDphi << ", map: " << map_DetaDphi_bin[ pair_DetaDphi ] << endl;
                hDetaDphi->GetYaxis()->SetBinLabel( kWin+1, Form("dEta_%.2f_dPhi_%.2f", min_dEta + etaStep*i, min_dPhi + phiStep*j ) );


//                cout << "min_dEta + etaStep*i, min_dPhi + phiStep*j: " << min_dEta + etaStep*i << " " << min_dPhi + phiStep*j << endl;

                kWin++;
            }
        cout << "N unique dEta dPhi = " << kWin << ", check map size: = " << map_DetaDphi_bin.size() << endl;
    }
    

    // ############
    void addTrack( int pid, double eta, double phi, double pt, int charge, double weight = 1.0 )
    {
        //        return;
        for( int eWin = 0; eWin < nWPs; eWin++ )
            wp[eWin]->addTrack( pid, eta, phi, pt, charge, weight );

        // if (Forward pid... charge..)
        // {
        //      histF_n->Fill(eta, phi, weight);
        //      histF_pt->Fill(eta, phi, weight);
        // }
        //
        // if (Backward pid... charge..)
        // {
        //      histB_n->Fill(eta, phi, weight);
        //      histB_pt->Fill(eta, phi, weight);
        // }
//        _w_wMinus1_F += weight*(weight-1);
//        _w_wMinus1_PF += weight*(weight-1)*pt;
////                    _w_wMinus1_PF2 += weight*(weight-1)*pt*pt;

//        _ptF2  += pt*pt *weight*weight;

    }



    void fillHistWithValue( const char *varName, double value, bool forceFilling = false )
    {
//        return;
        int varId = mapVar[varName];
        if( varId != 0 || forceFilling ) // because if varName is not in the varNames array, map will add a new key-value pair with value = 0.
        {
            hAllWins->Fill(  varId, currentWinId, currentSubsampleId, value  );
            hDetaDphi->Fill(  varId, currentDetaDphiPairId, currentSubsampleId, value  );
        }
        // here we are forcing filling "Nevents" (varId==0) by hand
    }


    // ############
    void finishEvent( int subId )
    {
        currentSubsampleId = subId;

        int sizeOfMap = map_DetaDphi_bin.size();
//        cout << "in finishEvent: " << nEta << endl;
        for( int wpId = 0; wpId < nWPs; wpId++ ) // loop over win pairs
        {
            currentWinId = wpId;

            array<int, 2> pair_DetaDphi = { (int)round(  eSep[wpId]*100 ), (int)round( phiSep[wpId] *100 ) };
            currentDetaDphiPairId = map_DetaDphi_bin[ pair_DetaDphi ];

            if ( currentDetaDphiPairId >= sizeOfMap )
                cout << "AHTUNG!!! currentDetaDphiPairId > map_DetaDphi_bin.size() !" << endl;

//            cout << "check map size: = " << map_DetaDphi_bin.size() << endl;
//            int aa;
//            cin >> aa;

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
        } // end of loop over win pairs

    }   // end of finishEvent



};


#endif












