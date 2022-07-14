#ifndef _WinPairBase
#define _WinPairBase


#include "TH1.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TString.h"



#include <iostream>
#include <map>
using namespace std;







// #############
struct WinPairBase
{
public:
    // acceptance of the windows:
    double eWin[4]; // eta for B (min, max), F (min, max),
    double ptWin[2];
    int partTypes[4];   // F,B, X,Y  // 0 - NO PID SELECTION
    int partCharges[4]; // F,B, X,Y  // 0 - NO CHARGE SELECTION

    // June 2021: new for ratio-ratio when denominator is in full acceptance
    bool fullEtaForDenom;
    double etaForDenomMin;
    double etaForDenomMax;

    WinPairBase() // int _subsId = -1)
    {
//        cout << "constructor WinPairBase" << endl;
        for( int i = 0; i < 4; i++ )
        {
            partTypes[i] = 0;
            partCharges[i] = 0;
        }

        //
        fullEtaForDenom = false;
        etaForDenomMin = -0.8;
        etaForDenomMax = 0.8;
    }

    ~WinPairBase()
    {
    }

    void setParticleTypes( int *_pTypes, int *_pCharges )
    {
        for( int i = 0; i < 4; i++ )
        {
            partTypes[i] = _pTypes[i];
            partCharges[i] = _pCharges[i];
        }
    }

    void setWindows( double _eMinB, double _eMaxB, double _eMinF, double _eMaxF, double _ptMin, double _ptMax
                     , bool _fullEtaForDenom = false, double _etaForDenomMin = -0.8, double _etaForDenomMax = 0.8 )
    {
        eWin[0] = _eMinB;
        eWin[1] = _eMaxB;
        eWin[2] = _eMinF;
        eWin[3] = _eMaxF;

        ptWin[0] = _ptMin;
        ptWin[1] = _ptMax;

//        TString strPID = Form( "pid_%d_%d_%d_%d_charge_%d_%d_%d_%d",
//                               partTypes[0], partTypes[1], partTypes[2], partTypes[3],
//                partCharges[0], partCharges[1], partCharges[2], partCharges[3]
//                );

        //        cout << "dir = " << dir << endl;
//        strAccumHistName = Form( "histAccumulatedValues_%s_%s_cW%d_cBin%d_eta_B_%.1f_%.1f_F_%.1f_%.1f_pt_%.1f_%.1f",
//                                strPrefix, strPID.Data(), iCW, cBin, eWin[0], eWin[1], eWin[2], eWin[3], ptWin[0], ptWin[1]
//                );
//        if ( subsampleId >=0 ) // i.e. work with subsamples
//            strAccumHistName.Append( Form( "_subs%d", subsampleId ) );

//        cout << strAccumHistName << endl;
//        if ( dir != 0x0 ) // get histogram from file
//            histAccumulatedValues = (TH1D*)dir->Get( strName );
//        else
//        {
// //            cout << "initializing window: " << strName.Data() << endl;
// //            histAccumulatedValues = new TH1D( strName //"histAccumulatedValues"
// //                                              , strName, nVars,-0.5,nVars-0.5);

// ////            TString gArrayMemberNames[nVars];
// //            for( int i=0; i < nVars; i++ )
// //                histAccumulatedValues->GetXaxis()->SetBinLabel( i+1, varNames[i] );
//        }

        fullEtaForDenom = _fullEtaForDenom;
        etaForDenomMin = _etaForDenomMin;
        etaForDenomMax = _etaForDenomMax;

    }


//    ClassDef(WinPairBase, 1);


};


#endif
