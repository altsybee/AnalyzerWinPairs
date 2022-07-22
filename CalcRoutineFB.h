#ifndef SimpleCalculations_cxx
#define SimpleCalculations_cxx


#include "TH1.h"
#include "TH1D.h"
#include "TH3D.h"
#include "THn.h"
#include "TDirectory.h"
#include "TString.h"
#include "TGraphErrors.h"

#include <iostream>
#include <string>
using namespace std;


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

    "avNF",
    "avNB",

    "avPtF",
    "avPtB",

    // Dec2019: bonus:

    "bcorr_XY",
    "nu_dyn_XY",
    "sigma_XY",

    "nu_dyn_FY",
    "sigma_FY",

    "nu_dyn_XB",
    "sigma_XB",


    "avX",
    "avY",



    "corr_rr_formula",
    "corr_rr_direct",

    "corr_R2_aa",
    "corr_R2_bb",
    "corr_R2_ab",
    "corr_R2_ba",


    // ##### meanPt-meanPt VS ETA (July 2022):
    "avPtF_avPtB_direct",
    "avPtF_avPtB_formula",


    // new ratio-ratio - June 2021: when denominator is in full acceptance
    "corr_rr_FULL_ETA_DENOM_formula",
    "corr_rr_FULL_ETA_DENOM_direct",

    "corr_rPt_formula",   // new ratio-meanPt - July 2022: when expansion for <pT> is done as <nB*PX>/.. + <nY*nX>/.. - <nB*nX>/.. - <nY*PX>/..
    "corr_rPt_direct",

    // new ratio-pt - April 2021: when pt is not averaged in each event
    "corr_rSumPt_formula",
    "corr_rSumPt_direct",


    //
    "coeff_dptdpt",
    "coeff_ptpt_CHECK_VV_formula",
    "coeff_ptpt_BOZEK_2017",
    "C_BOZEK_F",
    "C_BOZEK_B",
};
const int nObs = sizeof(obsNames)/sizeof(*obsNames); // number of observables




// ################################################
// check if particles are identical: F==B, X==Y
bool ifIdenticalParticlesFB( int *_pTypes, int *_pCharges )
{
    if( _pTypes[0] != _pTypes[1] ) return false;
    if( _pTypes[2] != _pTypes[3] ) return false;
    if( _pCharges[0] != _pCharges[1] ) return false;
    if( _pCharges[2] != _pCharges[3] ) return false;
    return true;
}





// ################################################
//const int MAX_N_WIN_PAIRS = 500;//150;
struct CalcWithSubsamples
{
    TH3D *h3D; // incoming hist with variables
    map<string, int> mapVarIdByName;

    // current status vars:
    TH1D *currentVarFullHist;
    int currenWinId;
    int currenSubId;
    double currenEtaSep;
    double currenPhiSep;

    TGraphErrors *grVsEta[nObs];
    double integrals[nObs];

    int nVars;
    int nWinPairs;
    int nSubs;

    TH3D *histCalcObs;   // hist with calculated observables
    map<string, int> mapObsIdByName;   // observable name-binId correspondance

    double eSizeNum;
    double eSizeDenom;
    bool if_Identical_FB_XY;

    CalcWithSubsamples()
    {
        h3D = 0x0;
    }

    TGraphErrors *getGraph( const char *strName )
    {
        int idObs = mapObsIdByName[strName];
        return grVsEta[idObs];
    }

    double getIntegral( const char *strName )
    {
        int idObs = mapObsIdByName[strName];
        return integrals[idObs];
    }

    void takeEtaSep( int _iWin )
    {
        double _etaSep = -1000;
        double _phiSep = -1000;

        TString strEtaBin;
        if ( h3D )
            strEtaBin = h3D->GetYaxis()->GetBinLabel( _iWin + 1 );

        if(0) cout << "strEtaBin: " << strEtaBin << endl;

        if ( strEtaBin.Contains("etaB") ) // i.e. we do pair-by-pair hist3D
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
            _etaSep = round( (eFpos - eBpos)*100 ) / 100;

            double phiBpos = ( phiBounds[0] + phiBounds[1] )/2;
            double phiFpos = ( phiBounds[2] + phiBounds[3] )/2;
            _phiSep = round( (phiFpos - phiBpos)*100 ) / 100;
        }
        else // if dEta-dPhi
        {
            float dEta, dPhi;
            sscanf( strEtaBin.Data(), "dEta_%f_dPhi_%f", &dEta, &dPhi );
//            dEta = h3D->GetYaxis()->GetBinCenter();
            _etaSep = round( dEta *100 ) / 100;
            _phiSep = round( dPhi *100 ) / 100;
            if(0)
                cout << "_etaSep = " << _etaSep << ", _phiSep = " << _phiSep << endl;
        }

        if ( _etaSep < -999 ) { cout << "AHTUNG!!!  _etaSep < -999!" << endl; int aa; cin>>aa; }
        if ( _phiSep < -999 ) { cout << "AHTUNG!!!  _phiSep < -999!" << endl; int aa; cin>>aa; }

        currenEtaSep = _etaSep;
        currenPhiSep = _phiSep;
    }


    // #####################
    void calc(  double _eSizeNum, double _eSizeDenom, bool _if_Identical_FB_XY )
    {
        eSizeNum = _eSizeNum;
        eSizeDenom = _eSizeDenom;
        if_Identical_FB_XY = _if_Identical_FB_XY;

        if(h3D)
        {
            nVars = h3D->GetNbinsX();
            nWinPairs = h3D->GetNbinsY();
            nSubs = h3D->GetNbinsZ();
        }

        createBinNameBinIdMap();


        // !!! check for empty bins:
        if(0)for( int i = 0; i < nVars; i++ )
            for( int j = 0; j < nWinPairs; j++ )
                for( int k = 0; k < nSubs; k++ )
                {
                    double binContent = h3D->GetBinContent(i+1, j+1, k+1);
                    if ( binContent < 0.00001 )
                    {
                        cout << "AHTUNG!!! binContent < 0.00001 !!" << " i=" << i << " j=" << j << " k=" << k <<
                             "bin label: " << h3D->GetXaxis()->GetBinLabel(i+1) << endl;
                    }

                }

        cout << "##### nObs = " << nObs << ", nVars = " << nVars << ", nSubs = " << nSubs << ", nWinPairs = " << nWinPairs << ", eSizeNum = " << eSizeNum << ", eSizeDenom = " << eSizeDenom << ", nSubs = " << nSubs << endl;

        // create graphs per each observable
        for( int iObs = 0; iObs < nObs; iObs++ )
        {
            grVsEta[iObs] = new TGraphErrors;
            integrals[iObs] = 0;
        }


        // calculated observables
        TString strHistObsTitle = Form( "%s_observables", h3D->GetName() );
        histCalcObs = new TH3D( strHistObsTitle, strHistObsTitle  //strAccumHistName , strAccumHistName //"histAccumulatedValues"
                                , nObs,-0.5, nObs-0.5
                                , nWinPairs, -0.5, nWinPairs-0.5
                                , nSubs+1, -0.5, nSubs-0.5 + 1  // +1 more for full calc (=sum of subsamples)
                                );


        // set x labels
        for( int i=0; i < nObs; i++ )
        {
            mapObsIdByName.insert( pair<string, int>( obsNames[i], i ) );
            histCalcObs->GetXaxis()->SetBinLabel( i+1, obsNames[i] );
        }


        // loop over win pairs:
        for ( int iWin = 0; iWin < nWinPairs; iWin++ )
        {
//            cout << "###### STARTING iWin = " << iWin << endl;

            // find out eSep, phiSep:
            takeEtaSep( iWin );
//            cout << "eSep = " << eSep << endl;


            TH1D *histVarValues = h3D->ProjectionX( "_px", iWin+1, iWin+1, 1, nSubs ); // project all subsamples on axis

            currentVarFullHist = histVarValues;
            currenWinId = iWin;
            currenSubId = nSubs; // i.e. calc for sum of all subsamples.
            finalCalc();

            // subsamples
            for ( int iSub = 0; iSub < nSubs; iSub++ )
            {
                // cout << "iSub=" << iSub << ", iType=" << iType << ", iCW=" << iCW << ", iWin=" << iWin << endl;
                currenSubId = iSub;
                finalCalc();
            }

            // calc mean and errors:
            for( int iObs = 0; iObs < nObs; iObs++ )
            {
                // mean
                double mean = histCalcObs->GetBinContent( iObs+1, iWin+1, nSubs+1 );


                // stdDev:
                double std_dev = 0;
                for ( int iSub = 0; iSub < nSubs; iSub++)
                {
                    double value = histCalcObs->GetBinContent( iObs+1, iWin+1, iSub+1 );
                    float diff = value - mean;
                    std_dev += diff*diff;
                }
                std_dev = sqrt(std_dev / nSubs / (nSubs-1) );

//                cout << "iObs = " << iObs << ", obs = " << obsNames[iObs] << ": mean +/- std_dev: " << mean << "   " << std_dev << endl;

                histCalcObs->SetBinError( iObs+1, iWin+1, nSubs+1, std_dev );

                //                cout << wpObs[iWin]->histCalcObs->GetXaxis()->GetBinLabel( bin+1 ) << " = " << mean << ", std_dev = " << std_dev << endl;

                // add points to graphs:
                grVsEta[iObs]->SetPoint( iWin, currenEtaSep, mean );
                grVsEta[iObs]->SetPointError( iWin, 0, std_dev );

                integrals[iObs] += mean; //*eSizeNum;
            }
        }  // end of loop over win pairs

    } // end of calc()





    // ###########
    void createBinNameBinIdMap() // const char* valueBinName, const char* nEvBinName )
    {
        for( int i=0; i < h3D->GetNbinsX(); i++ )
        {
            const char *strName = h3D->GetXaxis()->GetBinLabel(i+1);
            mapVarIdByName.insert( pair<string,int>( strName, i ) );
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
          std::cerr << " Obs " << obsName << ": out of range error: " << oor.what() << '\n';
        }
        if( obsId >= 0 )
            histCalcObs->SetBinContent( obsId+1, currenWinId+1, currenSubId+1, value );
        //        cout << "obsName = " << obsName << ", mapObsIdByName.at( obsName ) = " << mapObsIdByName.at( obsName ) << ", value = " << value << endl;
    }



    double getValueByName( const char* binName )
    {
//        cout << binName << endl;

        double value = -1000;
        int varId = -1;   //varIdByName( binName );

        try
        {
          varId = mapVarIdByName.at( binName );
        }
        catch (const std::out_of_range& oor) {
          std::cerr << " Var " << binName << ": out of range error: " << oor.what() << '\n';
        }


        if ( varId == -1 )
            return -1000;
//        cout << "getValueByName: binName = " << binName << ", mapVarIdByName[binName]+1 = " << mapVarIdByName[binName]+1 << endl;
        if ( currenSubId < nSubs )
            value = h3D->GetBinContent( varId+1, currenWinId+1, currenSubId+1 );
        else // i.e. we calc now the full hist (=sum of all subsamples)
            value = currentVarFullHist->GetBinContent(  varId+1 );

        return value;

    }
//    void setValueByName( const char* binName, double value )
//    {
//        int varId = varIdByName( binName );
//        if( varId >= 0 )
//            currentVarFullHist->SetBinContent(  varIdByName( binName )+1, value );
//    }



    // ###########
    void finalCalc() // /*int iWin,*/ double eSep, double phiSep ) //, bool if_Identical_FB_XY, double eSizeNum, double eSizeDenom )
    {
        double eSep = currenEtaSep;
        double phiSep = currenPhiSep;

        // we assume that eSep=0 means windows are completely overlapped!
        bool identical = ( eSep==0 && phiSep==0 && if_Identical_FB_XY );
//        cout << "finalCalc: " << eSep << endl;


        //  !!! ne nado average tut!!!      averageOverEvents( "sumPtAllEvF"  ,  "Nevents" );
        //  !!! ne nado average tut!!!      averageOverEvents( "sumPtAllEvB"  ,  "Nevents" );
        //  !!! ne nado average tut!!!      averageOverEvents( "piFpjB"       ,  "Nevents" );
        //  !!! ne nado average tut!!!      averageOverEvents( "nF*sum_pB"    ,  "Nevents" );
        //  !!! ne nado average tut!!!      averageOverEvents( "nB*sum_pF"    ,  "Nevents" );
        //  !!! ne nado average tut!!!      averageOverEvents( "pipjF"        ,  "Nevents" );
        //  !!! ne nado average tut!!!      averageOverEvents( "(nF-1)*sum_pF",  "Nevents" );
        //  !!! ne nado average tut!!!      averageOverEvents( "nF*(nF-1)"    ,  "Nevents" );
        //  !!! ne nado average tut!!!      averageOverEvents( "pipjB"        ,  "Nevents" );
        //  !!! ne nado average tut!!!      averageOverEvents( "(nB-1)*sum_pB",  "Nevents" );
        //  !!! ne nado average tut!!!      averageOverEvents( "nB*(nB-1)"    ,  "Nevents" );



        //        averageOverEvents( 0 /*corr type*/, "Nevents" );
        //        averageOverEvents( 1 /*corr type*/, "b_Nevents" );
        //        averageOverEvents( 2 /*corr type*/, "fb_Nevents" );

        //        cout << "mapVar[  Nf  ] = " << mapVar[ "Nf" ]  << endl;

        //        int aa;
        //        cin >> aa;
        //        calc_bCorr( 0, "Nx", "Ny", "Nx2", "Ny2", "Nx*Ny" );
//        calc_bCorr( 0, "Nf", "Nb", "Nf2", "Nb2", "Nf*Nb" );
//        calc_bCorr( 1, "NfPb_Nf", "NfPb_Pb", "NfPb_Nf2", "NfPb_Pb2", "NfPb_Nf_Pb" );
//        calc_bCorr( 2, "PfPb_Pf", "PfPb_Pb", "PfPb_Pf2", "PfPb_Pb2", "PfPb_Pf_Pb" );
        //        calc_bCorr( 1, mapVar[ "NfPb_Nf" ] );
        //        calc_bCorr( 2, mapVar[ "PfPb_Pf" ] );


//        fillHistWithValue( "bcorr_FB",      bcorr[0] );
//        fillHistWithValue( "bcorr_PtN",     bcorr[1] );
//        fillHistWithValue( "bcorr_PtPt",    bcorr[2] );
////        fillHistWithValue( "nu_dyn_FB",     nu_dyn[0] );
//        fillHistWithValue( "nu_dyn_PtN",    nu_dyn[1] );
//        fillHistWithValue( "nu_dyn_PtPt",   nu_dyn[2] );
////        fillHistWithValue( "sigma_FB",      sigma[0] );
//        fillHistWithValue( "sigma_PtN",     sigma[1] );
//        fillHistWithValue( "sigma_PtPt",    sigma[2] );


//        fillHistWithValue( "avNF", avF[0] );
//        fillHistWithValue( "avNB", avB[0] );

//        fillHistWithValue( "avPtF", avF[2] );
//        fillHistWithValue( "avPtB", avB[2] );




        // ##############################
        // now calc the observables

        double _nEvents     = getValueByName(  "Nevents" );

        double FB = getValueByName(  "Nf*Nb" ) / _nEvents;
        double XY = getValueByName(  "Nx*Ny" ) / _nEvents;
        double FY = getValueByName(  "Nf*Ny" ) / _nEvents;
        double XB = getValueByName(  "Nb*Nx" ) / _nEvents;

        double F = getValueByName(  "Nf" )  / _nEvents;
        double B = getValueByName(  "Nb" )  / _nEvents;
        double X = getValueByName(  "Nx" )  / _nEvents;
        double Y = getValueByName(  "Ny" )  / _nEvents;

        double F2 = getValueByName(  "Nf2" ) / _nEvents;
        double X2 = getValueByName(  "Nx2" ) / _nEvents;
        double B2 = getValueByName(  "Nb2" ) / _nEvents;
        double Y2 = getValueByName(  "Ny2" ) / _nEvents;



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
            fillHistWithValue( "corr_R2_ab", R2_aa );
            fillHistWithValue( "corr_R2_ba", R2_bb );
        }


        // ##############################
        // ##### coeff ratio-ratio
        // ##############################
        double rFB = getValueByName(  "Nf_OVER_Nx_vs_Nb_OVER_Ny" ) / getValueByName(  "xy_Nevents" );
        double ratioF = getValueByName(  "Nf_OVER_Nx" ) / getValueByName(  "xy_Nevents" );
        double ratioB = getValueByName(  "Nb_OVER_Ny" ) / getValueByName(  "xy_Nevents" );
        {


            double corr_rr_direct = F/eSizeNum * X/eSizeDenom / (F/eSizeNum + X/eSizeDenom) * ( rFB/ratioF/ratioB    - 1 );
            double corr_rr_formula = F/eSizeNum * X/eSizeDenom / (F/eSizeNum + X/eSizeDenom) * (FB/F/B + XY/X/Y - FY/F/Y - XB/X/B   );

            if (identical)
            {
                corr_rr_direct = F/eSizeNum * X/eSizeDenom / (F/eSizeNum + X/eSizeDenom) * ( rFB/ratioF/ratioB - 1/F - 1/X   - 1 );
                corr_rr_formula = F/eSizeNum * X/eSizeDenom / (F/eSizeNum + X/eSizeDenom) * (F2/F/F + X2/X/X - FY/F/Y - XB/X/B - 1/F - 1/X );
            }
            fillHistWithValue( "corr_rr_direct", corr_rr_direct );
            fillHistWithValue( "corr_rr_formula", corr_rr_formula );
            if(0)cout << "corr_rr_formula = " << corr_rr_formula
                 << ", F = " << F << ", X = " << X << ", Y = " << Y
                 << ", F2 = " << F2 << ", X2 = " << X2 << ", FY = " << FY << ", XB = " << XB
                 << ", _nEvents = " << _nEvents
                 << endl;

            // when denominator is in full acceptance:
            {
                double FX = FY;


                double corr_rr_FULL_ETA_DENOM_formula = F/eSizeNum * X/eSizeDenom / ( F/eSizeNum + X/eSizeDenom ) * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X  ) ;
                double corr_rr_FULL_ETA_DENOM_direct =  F/eSizeNum * X/eSizeDenom / ( F/eSizeNum + X/eSizeDenom ) * ( rFB/ratioF/ratioB - 1   -1/X  ) ;

                if (identical)
                {
                    corr_rr_FULL_ETA_DENOM_formula = F/eSizeNum * X/eSizeDenom / ( F/eSizeNum + X/eSizeDenom ) * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X  - 1/F ) ;
                    corr_rr_FULL_ETA_DENOM_direct =  F/eSizeNum * X/eSizeDenom / ( F/eSizeNum + X/eSizeDenom ) * ( rFB/ratioF/ratioB - 1   -1/X   - 1/F ) ;
                }

                //            fillHistWithValue( "corr_rr_FULL_ETA_DENOM_formula", (F+B)/2 * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X) );
                //            fillHistWithValue( "corr_rr_FULL_ETA_DENOM_formula", F*B/(F+B) * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X) );
                fillHistWithValue( "corr_rr_FULL_ETA_DENOM_formula", corr_rr_FULL_ETA_DENOM_formula );
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
            double PxPy_avPx      = getValueByName(  "PxPy_avPx" );
            double Nb_OVER_Ny_vs_Px = getValueByName(  "Nb_OVER_Ny_vs_avPx" );

            fillHistWithValue( "corr_rPt_direct", X/eSizeNum * ( Nb_OVER_Ny_vs_Px/PxPy_avPx/ratioB - 1 ) );

            // formula
            double nB_PX = getValueByName(  "nB*PX" ) / _nEvents;
            double nY_PX = getValueByName(  "nY*PX" ) / _nEvents;
            double PX      = getValueByName(  "PX" )  / _nEvents;

            // new ratio-meanPt formula - July 2022: when expansion for <pT> is done as <nB*PX>/.. + <nY*nX>/.. - <nB*nX>/.. - <nY*PX>/..
//            double _avSumPtInEvX      = getValueByName(  "sumPtAllEvX" )  / _nEvents;
            double rPt_full_formula = X/eSizeNum * (  nB_PX /B/PX  + XY/X/Y - XB/X/B - nY_PX /Y/PX  ); // here even if wins are identical, +1/X-1/X cancels!
            if (identical)
            {
//                double subtrX = identical ? 1/X : 0;
                double nX_PX = getValueByName(  "nX*PX" ) / _nEvents;
                rPt_full_formula = X/eSizeNum * (  nB_PX /B/PX  + X2/X/X - XB/X/B - nX_PX /X/PX );
            }
            fillHistWithValue( "corr_rPt_formula", rPt_full_formula );
        }




        // nu_dyn, sigma, avX, avY
        if(1)
        {
//            double numerator = XY - X * Y;
//            double denominator_bCorr = X2 - X*X;

            fillHistWithValue( "avX", F );
            fillHistWithValue( "avY", B );

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


        // ##### coeff avPt-avPt formula (new expansion, July 2022)
        {
            double PF_PB = getValueByName(  "PF*PB" ) / _nEvents;
            double nF_PB = getValueByName(  "nF*PB" ) / _nEvents;
            double nB_PF = getValueByName(  "nB*PF" ) / _nEvents;

            double PF = getValueByName(  "PF" ) / _nEvents;
            double PB = getValueByName(  "PB" ) / _nEvents;


            double avPtF_avPtB_formula = F/eSizeNum * ( PF_PB/PF/PB + FB/F/B - nF_PB/F/PB  - nB_PF/B/PF );

            if (identical)
            {
                double PF2 = getValueByName(  "PF2" ) / _nEvents;
                double nF_PF = getValueByName(  "nF*PF" ) / _nEvents;
                double piF2 = getValueByName(  "piF2" ) / _nEvents;

//                double subtrX = identical ? 1/X : 0;
                // DODELAT' PF_PB !!!!!
                avPtF_avPtB_formula = F/eSizeNum * ( (PF2 - piF2)/PF/PF + F2/F/F  - nF_PF/F/PF  - nF_PF/F/PF     - 1/F + 1/F + 1/F );
            }

            fillHistWithValue( "avPtF_avPtB_formula", avPtF_avPtB_formula );
        }

        // ##### coeff avPt-avPt direct
        {
            double Pf = getValueByName(  "PfPb_avPf" ) / getValueByName(  "fb_Nevents" );
            double Pb = getValueByName(  "PfPb_avPb" ) / getValueByName(  "fb_Nevents" );
            double Pf_Pb = getValueByName(  "PfPb_avPf_avPb" ) / getValueByName(  "fb_Nevents" );

            fillHistWithValue( "avPtF_avPtB_direct", F/eSizeNum * ( Pf_Pb/(Pf * Pb) - 1 ) );
//            cout << " >> Pf_Pb = " << Pf_Pb << ", Pf = " << Pf << ", Pb = " << Pb << endl;
        }







        // ##### dptdpt
        if(0)
        {
            double _avPtAllEvF     = getValueByName(  "sumPtAllEvF" )  / ( F*_nEvents );
            double _avPtAllEvB     = getValueByName(  "sumPtAllEvB" )  / ( B*_nEvents );
            double _piFpjB = getValueByName(  "piFpjB" );
            double _nF_sum_pB = getValueByName(  "nF*PB" ); // /2; // /2 is TMP!!!
            double _nB_sum_pF = getValueByName(  "nB*PF" ); // /2;

            // from /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/results/_common_files/routine.h:169
            fillHistWithValue( "coeff_dptdpt", (
                                   ( _piFpjB - _nF_sum_pB*_avPtAllEvF - _nB_sum_pF*_avPtAllEvB)/ (FB*_nEvents)
                                   + _avPtAllEvF*_avPtAllEvB
                                   ) / (_avPtAllEvF*_avPtAllEvB) );
        }

        //        }



        // ### coeff ptpt with means in denom
        if(0){

            fillHistWithValue( "coeff_ptpt_CHECK_VV_formula",
                               (getValueByName(  "Nf" )
                               *getValueByName(  "Nb" )
                    *getValueByName(  "PF*PB" )
                    -
                    getValueByName(  "Nb" )
                    *getValueByName(  "PF" )
                    *getValueByName(  "nF*PB" )// /2 // /2 is TMP!!!
                    -
                    getValueByName(  "Nf" )
                    *getValueByName(  "PB" )
                    *getValueByName(  "nB*PF" )// /2 // /2 is TMP!!!
                    +
                    getValueByName(  "PF" )
                    *getValueByName(  "PB" )
                    *getValueByName(  "Nf*Nb" )
                    ) / (
                        getValueByName(  "PF" )
                    *getValueByName(  "PB" )
                    *getValueByName(  "Nf*Nb" )
                    )
                    );
        }

        // ### coeff ptpt BOZEK_2017_USING_QM17_RESULTS_1704.02777.pdf
        if(0)
        {
            double _avPtAllEvF     = getValueByName(  "sumPtAllEvF" )  / ( F*_nEvents );
            double _avPtAllEvB     = getValueByName(  "sumPtAllEvB" )  / ( B*_nEvents );

            double Pf = getValueByName(  "PfPb_Pf" ) ;
            double Pb = getValueByName(  "PfPb_Pb" ) ;
            double Pf_Pb = getValueByName(  "PfPb_Pf_Pb" );

            double _pipjF = getValueByName(  "pipjF" );
            double _pipjB = getValueByName(  "pipjB" );

            double secondTermF     = getValueByName(  "(nF-1)*PF" );
            double secondTermB     = getValueByName(  "(nB-1)*PB" );

            double nn_1F     = getValueByName(  "nF*(nF-1)" );
            double nn_1B     = getValueByName(  "nB*(nB-1)" );


            double C_Bozek_F_nPairs_OUTSIDE_sum = ( _pipjF*2 - 2*_avPtAllEvF*secondTermF + _avPtAllEvF*_avPtAllEvF*nn_1F ) / nn_1F;
            double C_Bozek_B_nPairs_OUTSIDE_sum = ( _pipjB*2 - 2*_avPtAllEvB*secondTermB + _avPtAllEvB*_avPtAllEvB*nn_1B ) / nn_1B;

            //            double C_Bozek_F = pipjF_avPerEv/nEventsF*2 - 2*avPtAllEvF*secondTermFnew/nEventsF + avPtAllEvF*avPtAllEvF;
            //            double C_Bozek_B = pipjB_avPerEv/nEventsB*2 - 2*avPtAllEvB*secondTermBnew/nEventsB + avPtAllEvB*avPtAllEvB;

            fillHistWithValue( "coeff_ptpt_BOZEK_2017", (Pf_Pb - Pf*Pb)/sqrt( C_Bozek_F_nPairs_OUTSIDE_sum * C_Bozek_B_nPairs_OUTSIDE_sum ) );

            //            if(0)cout <<  ">>> C_Bozek_B_nPairs_OUTSIDE_sum: " << _pipjB << " "<< _avPtAllEvB << " " << _avPtAllEvB << " " << nn_1B << ", coeff_ptpt_BOZEK_2017=" << coeff_ptpt_BOZEK_2017 << endl;

            fillHistWithValue( "C_BOZEK_F", C_Bozek_F_nPairs_OUTSIDE_sum );
            fillHistWithValue( "C_BOZEK_B", C_Bozek_B_nPairs_OUTSIDE_sum );
        }



    } // end of final calc


}; // end of class



#endif

