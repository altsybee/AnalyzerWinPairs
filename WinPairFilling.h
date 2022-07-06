#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include "TDirectory.h"
#include "TString.h"

#include "WinPairBase.h"

#include <iostream>
#include <map>
#include <string>


using namespace std;


//bool MAKE_TEST_FIRST_FILL_FOR_NON_EXISTING_VALUES = true;

// available vars:
const char *availableVarNames[] = {    // now NOT really used!!! just to see the list as it was on July 6, 2022
    "Nevents",

    "f_Nevents",
    "b_Nevents",

    "x_Nevents",
    "y_Nevents",

    "fb_Nevents",
    "xy_Nevents",

    "fx_Nevents",
    "fy_Nevents",
    "bx_Nevents",
    "by_Nevents",

    // f, b
    "Nf",
    "Nb",
    "N2_f",
    "N2_b",
    "Nf_Nb",

    // x, y
    "Nx",
    "Ny",
    "Nx2",
    "Ny2",
    "Nx_Ny",

    //
    "NfPb_Nf",
    "NfPb_Pb",
    "NfPb_Nf2",
    "NfPb_Pb2",
    "NfPb_Nf_Pb",

    "PfNb_Pf",
    "PfNb_Nb",
    "PfNb_Pf2",
    "PfNb_Nb2",
    "PfNb_Pf_Nb",

    "PfPb_Pf",
    "PfPb_Pb",
    "PfPb_Pf2",
    "PfPb_Pb2",
    "PfPb_Pf_Pb",



    "NxPy_Nx",
    "NxPy_Py",
    "NxPy_Nx2",
    "NxPy_Py2",
    "NxPy_Nx_Py",

    "PxNy_Px",
    "PxNy_Ny",
    "PxNy_Px2",
    "PxNy_Ny2",
    "PxNy_Px_Ny",

    "PxPy_Px",
    "PxPy_Py",
    "PxPy_Px2",
    "PxPy_Py2",
    "PxPy_Px_Py",


    // ##### mixes:
    "Nf_Nx",
    "Nf_Ny",
    "Nb_Nx",
    "Nb_Ny",

    "Nf_Px",
    "Nf_Py",
    "Nb_Px",
    "Nb_Py",

    "Pf_Nx",
    "Pf_Ny",
    "Pb_Nx",
    "Pb_Ny",

    "Pf_Px",
    "Pf_Py",
    "Pb_Px",
    "Pb_Py",


    // RATIOS:
    "Nf_OVER_Nb",
    "Nf_OVER_Nx",
    "Nf_OVER_Ny",

    "Nb_OVER_Nf",
    "Nb_OVER_Nx",
    "Nb_OVER_Ny",

    "Nx_OVER_Nf",
    "Nx_OVER_Nb",
    "Nx_OVER_Ny",

    "Ny_OVER_Nf",
    "Ny_OVER_Nb",
    "Ny_OVER_Nx",

    // RATIO - RATIO:
    "Nf_OVER_Nb_vs_Nx_OVER_Ny",
    "Nf_OVER_Nx_vs_Nb_OVER_Ny",

    "Nb_OVER_Nf_vs_Ny_OVER_Nx",
    "Nx_OVER_Nf_vs_Ny_OVER_Nb",

    // RATIO - pT:
    "Nf_OVER_Nx_vs_Pb",
    "Nf_OVER_Nx_vs_Py",
    "Nb_OVER_Ny_vs_Pf",
    "Nb_OVER_Ny_vs_Px",

    "Nx_OVER_Nf_vs_Pb",
    "Nx_OVER_Nf_vs_Py",
    "Ny_OVER_Nb_vs_Pf",
    "Ny_OVER_Nb_vs_Px",


    // added on 04.04.2019: for dptdpt etc (from /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/AliForwardBackwardAnalysis.cxx:1581)
    "sumPtAllEvF"  ,
    "sumPtAllEvB"  ,
    "piFpjB"       ,
    "nF*PB"    ,
    "nB*PF"    ,
    "pipjF"        ,
    "(nF-1)*PF",
    "nF*(nF-1)"    ,
    "pipjB"        ,
    "(nB-1)*PB",
    "nB*(nB-1)"    ,


    // to check VV comparison ptpt vs dptdpt: special terms:
    "PF*PB",
    "PF",
    "PB",

    // for new ratio-pt - April 2021
    "nY*PF",
    "nB*PX",
    "nY*PX",

    "nX*PB",
    "nF*PY",
    "nX*PY",

    "Nb_OVER_Ny*PF",
    "Nb_OVER_Ny*PX",
    "Nf_OVER_Nx*PB",
    "Nf_OVER_Nx*PY",

    "sumPtAllEvX"  ,
    "sumPtAllEvY"  ,



};


//const int nVars = sizeof(varNames)/sizeof(*varNames);



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
    double _w_wMinus1_X; double _w_wMinus1_Y; // for corrections in one window!

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
        _nX(0), _nY(0),    _w_wMinus1_X(0), _w_wMinus1_Y(0),
        _ptX(0), _ptY(0)
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
                }

        //        cout << "test2" << endl;

    }

    void finishEvent()
    {
        // to calc C:
        if(1)  // 0 to speed up
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

        _w_wMinus1_X = 0;
        _w_wMinus1_Y = 0;


        // extra info:
        pipjF = 0; //sum_piF = 0;
        pipjB = 0; //sum_piB = 0;

        piFpjB = 0;

    }

    //    ClassDef(WinPairFilling, 1);


};








const int MAX_N_WIN_PAIRS = 150;

// ####################################################################################################################
// class which collects all window pairs and performs filling of TH3D histogram with e-by-e calculated values
// ####################################################################################################################

class WinPairWrapper
{
public:
    WinPairFilling *wp[MAX_N_WIN_PAIRS]; //! // win pairs
    int nEta;
    TH3D *hAllWins; //!
    TString strAnLevel;

    // some vars to make flexible hist filling (i.e. don't fill if bin name is not in varNames array):
    int currentSubsampleId;
    int currentWinId;

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
        nEta++;
    }



    // ############
    void setHistAllWins( const char* strPrefix, int cBin, int nSub, const char* _varNames[], int _nVars )
    {
        strAnLevel = Form("%s", strPrefix);

        int nWins = nEta;

        TString strHistName = Form("hAllWins_%s_cBin%d", strAnLevel.Data(), cBin);


        hAllWins = new TH3D( strHistName, strHistName
                             //, nBinsX, -0.5, nBinsX-0.5
                             , _nVars, -0.5, _nVars-0.5
                             , nWins, -0.5, nWins-0.5
                             , nSub, -0.5, nSub-0.5
                             );

        // set labels x and prepare a map with vars
        for( int i = 0; i < _nVars; i++ )
        {
            mapVar.insert( pair<const char*,int>( _varNames[i], i ) );
            hAllWins->GetXaxis()->SetBinLabel( i+1, _varNames[i] ); //h_wp->GetXaxis()->GetBinLabel( i+1 ) );
        }
//        cout << "check mapVar size: " << mapVar.size() << endl;

        // set labels y
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
        int varId = mapVar[varName];
        if( varId != 0 || forceFilling ) // because if varName is not in the varNames array, map will add a new key-value pair with value = 0.
            hAllWins->Fill(  varId, currentWinId, currentSubsampleId, value  );
        // here we force filling "Nevents" (varId==0) by hand
    }


    // ############
    void finishEvent( int subId )
    {
        currentSubsampleId = subId;

        for( int eWin = 0; eWin < nEta; eWin++ )
        {
            currentWinId = eWin;
//            varCounter = 0;

            //
            WinPairFilling *w = wp[eWin];
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
            fillHistWithValue( "N2_f"  ,   w->_nF * w->_nF      );
            fillHistWithValue( "N2_b"  ,   w->_nB * w->_nB      );
            fillHistWithValue( "Nf_Nb" ,   w->_nF * w->_nB      );

            // NxNy:
            fillHistWithValue( "Nx"      ,   w->_nX             );
            fillHistWithValue( "Ny"      ,   w->_nY             );
            fillHistWithValue( "Nx2"     ,   w->_nX * w->_nX - w->_w_wMinus1_X   ); // !fullEtaForDenom ? _nX*_nX  : _nX*_nX - _w_wMinus1_X   );
            fillHistWithValue( "Ny2"     ,   w->_nY * w->_nY - w->_w_wMinus1_Y   ); // !fullEtaForDenom ? _nY*_nY  : _nY*_nY - _w_wMinus1_Y   );
            fillHistWithValue( "Nx_Ny"   ,   w->_nX * w->_nY      );


            // Nf_Nx Nf_Ny Nb_Nx Nb_Ny:
            fillHistWithValue( "Nf_Nx"  ,   w->_nF * w->_nX       );
            fillHistWithValue( "Nf_Ny"  ,   w->_nF * w->_nY       );
            fillHistWithValue( "Nb_Nx"  ,   w->_nB * w->_nX       );
            fillHistWithValue( "Nb_Ny"  ,   w->_nB * w->_nY       );


            // FROM /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/AliForwardBackwardAnalysis.cxx:1583
            fillHistWithValue( "sumPtAllEvF" , w->_ptF );
            fillHistWithValue( "sumPtAllEvB" , w->_ptB );
            fillHistWithValue( "sumPtAllEvX" , w->_ptX );
            fillHistWithValue( "sumPtAllEvY" , w->_ptY );

            // for dptdpt:
            fillHistWithValue( "piFpjB" , w->piFpjB );
            fillHistWithValue( "nF*PB" ,  w->_nF * w->_ptB );
            fillHistWithValue( "nB*PF" ,  w->_nB * w->_ptF );

            // for C when av over pairs is OUTSIDE sum:
            fillHistWithValue( "pipjF" ,          w->pipjF );
            fillHistWithValue( "(nF-1)*PF" ,      (w->_nF-1)*w->_ptF ); // (n-1)*sumPt
            fillHistWithValue( "nF*(nF-1)" ,          (w->_nF-1)*w->_nF ); // (n-1)*n

            fillHistWithValue( "pipjB" ,         w->pipjB );
            fillHistWithValue( "(nB-1)*PB" , (w->_nB-1)*w->_ptB ); // (n-1)*sumPt
            fillHistWithValue( "nB*(nB-1)" ,     (w->_nB-1)*w->_nB ); // (n-1)*n



            // to check VV comparison ptpt vs dptdpt: special terms:
            fillHistWithValue( "PF*PB" ,      w->_ptF*w->_ptB );
            //        fillHistWithValue( "nB*PF" ,            _nB*_ptF );
            //        fillHistWithValue( "nF*PB" ,            _nF*_ptB );  //_nF*_nB*_ptB );
            fillHistWithValue( "PF" ,              w->_ptF );
            fillHistWithValue( "PB" ,              w->_ptB );


            // for new ratio-pt (April 2021)
            //        fillHistWithValue( "nB*PF" , _nB*_ptF ); // already have it above
            fillHistWithValue( "nY*PF" , w->_nY*w->_ptF );
            fillHistWithValue( "nB*PX" , w->_nB*w->_ptX );
            fillHistWithValue( "nY*PX" , w->_nY*w->_ptX );

            //        fillHistWithValue( "nF*PB" , _nF*_ptB );
            fillHistWithValue( "nX*PB" , w->_nX*w->_ptB );
            fillHistWithValue( "nF*PY" , w->_nF*w->_ptY );
            fillHistWithValue( "nX*PY" , w->_nX*w->_ptY );


            if ( w->_nY > 0 )
            {
                fillHistWithValue( "Nb_OVER_Ny*PF" , w->_nB/(double)w->_nY * w->_ptF );
                fillHistWithValue( "Nb_OVER_Ny*PX" , w->_nB/(double)w->_nY * w->_ptX );
            }
            if ( w->_nX > 0 )
            {
                fillHistWithValue( "Nf_OVER_Nx*PB" , w->_nF/(double)w->_nX * w->_ptB );
                fillHistWithValue( "Nf_OVER_Nx*PY" , w->_nF/(double)w->_nX * w->_ptY );
            }



            // NfPb:
            if ( w->_nB > 0 )
            {
                fillHistWithValue( "b_Nevents" ,    1                  );

                fillHistWithValue( "NfPb_Nf" ,       w->_nF             );
                fillHistWithValue( "NfPb_Pb" ,       meanPtB             );
                fillHistWithValue( "NfPb_Nf2" ,      w->_nF*w->_nF      );
                fillHistWithValue( "NfPb_Pb2" ,      meanPtB*meanPtB      );
                fillHistWithValue( "NfPb_Nf_Pb" ,    w->_nF*meanPtB      );

                fillHistWithValue( "Nf_OVER_Nb" ,    w->_nF/(double)w->_nB      );
                fillHistWithValue( "Nx_OVER_Nb" ,    w->_nX/(double)w->_nB      );
                fillHistWithValue( "Ny_OVER_Nb" ,    w->_nY/(double)w->_nB      );

                fillHistWithValue( "Pb_Nx" ,   meanPtB*w->_nX      );
                fillHistWithValue( "Pb_Ny" ,   meanPtB*w->_nY      );


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

                fillHistWithValue( "PfNb_Pf" ,         meanPtF           );
                fillHistWithValue( "PfNb_Nb" ,        w->_nB            );
                fillHistWithValue( "PfNb_Pf2" ,     meanPtF*meanPtF      );
                fillHistWithValue( "PfNb_Nb2" ,      w->_nB*w->_nB      );
                fillHistWithValue( "PfNb_Pf_Nb" ,    meanPtF*w->_nB      );

                fillHistWithValue( "Nb_OVER_Nf" ,    w->_nB/(double)w->_nF      );
                fillHistWithValue( "Nx_OVER_Nf" ,    w->_nX/(double)w->_nF      );
                fillHistWithValue( "Ny_OVER_Nf" ,    w->_nY/(double)w->_nF      );

                fillHistWithValue( "Pf_Nx" ,   meanPtF*w->_nX      );
                fillHistWithValue( "Pf_Ny" ,   meanPtF*w->_nY      );

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
                fillHistWithValue( "fb_Nevents" ,   1                  );

                fillHistWithValue( "PfPb_Pf" ,       meanPtF             );
                fillHistWithValue( "PfPb_Pb" ,       meanPtB             );
                fillHistWithValue( "PfPb_Pf2" ,     meanPtF*meanPtF      );
                fillHistWithValue( "PfPb_Pb2" ,     meanPtB*meanPtB      );
                fillHistWithValue( "PfPb_Pf_Pb" ,   meanPtF*meanPtB      );

                fillHistWithValue( "Nx_OVER_Nf_vs_Pb" ,   w->_nX/(double)w->_nF *meanPtB      );
                fillHistWithValue( "Ny_OVER_Nb_vs_Pf" ,   w->_nY/(double)w->_nB *meanPtF      );

                fillHistWithValue( "Nx_OVER_Nf_vs_Ny_OVER_Nb" ,   w->_nX/(double)w->_nF * w->_nY/(double)w->_nB      );

                //            // for dptdpt:
                //            fillHistWithValue( "piFpjB" ,     piFpjB );
                //            fillHistWithValue( "nF*sum_pB" ,  _nF*(_nB*_ptB) );
                //            fillHistWithValue( "nB*sum_pF" ,  _nB*(_nF*_ptF) );

            }


            // y:
            if ( w->_nY > 0 )
            {
                fillHistWithValue( "y_Nevents" ,    1                  );

                fillHistWithValue( "NxPy_Nx" ,        w->_nX             );
                fillHistWithValue( "NxPy_Py" ,       meanPtY             );
                fillHistWithValue( "NxPy_Nx2" ,      w->_nX*w->_nX      );
                fillHistWithValue( "NxPy_Py2" ,      meanPtY*meanPtY      );
                fillHistWithValue( "NxPy_Nx_Py" ,    w->_nX*meanPtY      );


                fillHistWithValue( "Nf_OVER_Ny" ,    w->_nF/(double)w->_nY      );
                fillHistWithValue( "Nb_OVER_Ny" ,    w->_nB/(double)w->_nY      );
                fillHistWithValue( "Nx_OVER_Ny" ,    w->_nX/(double)w->_nY      );

                fillHistWithValue( "Nf_Py" ,   w->_nF*meanPtY       );
                fillHistWithValue( "Nb_Py" ,   w->_nB*meanPtY       );
            }

            // x:
            if ( w->_nX > 0 )
            {
                fillHistWithValue( "x_Nevents" ,    1                  );

                fillHistWithValue( "PxNy_Px" ,         meanPtX           );
                fillHistWithValue( "PxNy_Ny" ,        w->_nY            );
                fillHistWithValue( "PxNy_Px2" ,     meanPtX*meanPtX      );
                fillHistWithValue( "PxNy_Ny2" ,      w->_nY*w->_nY      );
                fillHistWithValue( "PxNy_Px_Ny" ,    meanPtX*w->_nY      );

                fillHistWithValue( "Nf_OVER_Nx" ,    w->_nF/(double)w->_nX      );
                fillHistWithValue( "Nb_OVER_Nx" ,    w->_nB/(double)w->_nX      );
                fillHistWithValue( "Ny_OVER_Nx" ,    w->_nY/(double)w->_nX      );

                fillHistWithValue( "Nf_Px" ,   w->_nF*meanPtX       );
                fillHistWithValue( "Nb_Px" ,   w->_nB*meanPtX       );
            }


            // xy:
            if ( w->_nX > 0 && w->_nY > 0 )
            {
                fillHistWithValue( "xy_Nevents" ,   1                  );
                fillHistWithValue( "PxPy_Px" ,       meanPtX             );
                fillHistWithValue( "PxPy_Py" ,       meanPtY             );
                fillHistWithValue( "PxPy_Px2" ,     meanPtX*meanPtX      );
                fillHistWithValue( "PxPy_Py2" ,     meanPtY*meanPtY      );
                fillHistWithValue( "PxPy_Px_Py" ,   meanPtX*meanPtY      );

                fillHistWithValue( "Nf_OVER_Nx_vs_Py" ,   w->_nF/(double)w->_nX*meanPtY      );
                fillHistWithValue( "Nb_OVER_Ny_vs_Px" ,   w->_nB/(double)w->_nY*meanPtX      );

                fillHistWithValue( "Nf_OVER_Nx_vs_Nb_OVER_Ny" ,   w->_nF/(double)w->_nX * w->_nB/(double)w->_nY      );
            }

            // fx:
            if ( w->_nF > 0 && w->_nX > 0 )
            {
                fillHistWithValue( "fx_Nevents" ,   1                  );

                fillHistWithValue( "Pf_Px" ,   meanPtF*meanPtX      );

                fillHistWithValue( "Nb_OVER_Nf_vs_Ny_OVER_Nx" ,   w->_nB/(double)w->_nF * w->_nY/(double)w->_nX      );
            }
            // fy:
            if ( w->_nF > 0 && w->_nY > 0 )
            {
                fillHistWithValue( "fy_Nevents" ,   1                  );

                fillHistWithValue( "Pf_Py" ,   meanPtF*meanPtY      );
                fillHistWithValue( "Nx_OVER_Nf_vs_Py" ,   w->_nX/(double)w->_nF *meanPtY      );
                fillHistWithValue( "Nb_OVER_Ny_vs_Pf" ,   w->_nB/(double)w->_nY *meanPtF      );
            }
            // bx:
            if ( w->_nB > 0 && w->_nX > 0 )
            {
                fillHistWithValue( "bx_Nevents" ,   1                  );

                fillHistWithValue( "Pb_Px" ,   meanPtB*meanPtX      );
                fillHistWithValue( "Nf_OVER_Nx_vs_Pb" ,   w->_nF/(double)w->_nX *meanPtB      );
                fillHistWithValue( "Ny_OVER_Nb_vs_Px" ,   w->_nY/(double)w->_nB *meanPtX      );
            }
            // by:
            if ( w->_nB > 0 && w->_nY > 0 )
            {
                fillHistWithValue( "by_Nevents" ,   1                  );
                fillHistWithValue( "Pb_Py" ,   meanPtB*meanPtY      );

                fillHistWithValue( "Nf_OVER_Nb_vs_Nx_OVER_Ny" ,   w->_nF/(double)w->_nB * w->_nX/(double)w->_nY      );
            }


            w->resetEvent();
        } // end of loop over eta wins

    }   // end of finishEvent( int subId )



};














