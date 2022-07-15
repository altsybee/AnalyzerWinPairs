#ifndef SimpleCalculations_cxx
#define SimpleCalculations_cxx



#include "TH1.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TDirectory.h"
#include "TString.h"
#include "TGraphErrors.h"

//#include "SimpleWinPair_forPythiaAnalysis_v3.h"
//#include "WinPair2019.h"
//#include "WinPairBase.h"
#include "WinPairFilling.h"



#include <iostream>
using namespace std;

// from https://stackoverflow.com/questions/4157687/using-char-as-a-key-in-stdmap
struct cmp_str
{
    bool operator()(char const *a, char const *b) const
    {
        return std::strcmp(a, b) < 0;
    }
};











// names of final quantities
const char *obsNames[] =
{
//    "bcorr_FB",
//    "bcorr_PtN",
//    "bcorr_PtPt",

//    "nu_dyn_FB",
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

    "corr_rPt_formula",
    "corr_rPt_direct",
    "corr_rPt_full_formula",   // new ratio-meanPt - July 2022: when expansion for <pT> is done as <nB*PX>/.. + <nY*nX>/.. - <nB*nX>/.. - <nY*PX>/..

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




// #############
class SimpleCalculations //: public WinPairBase //WinPair2019
{
    static int global_obs_hist_counter;  // = 0;
public:
    map<const char*, double, cmp_str> mapData;
//    double eSep;  // we assume that eSep=0 means windows are completely overlapped!
    // final quantities
    double bcorr[3];
    double nu_dyn[3];
    double sigma[3];

    double avF[3];
    double avB[3];

    map<const char*,int> mapObs;   // helps checking var duplicates
//    bool firstFill; // to speed up (?) hist filling: skip if obs is not requested (but always try fill the first time because of the map usage)


    TH1D *histAccumulatedValues;   // //! accumulated values for observables
    TH1D *histCalcObs;   // //!   hist with calculated observables

    SimpleCalculations()
    {
//        firstFill = true;
//        eSep = -1000;
        histAccumulatedValues = 0x0;

        //        cout << "global_obs_hist_counter = " << global_obs_hist_counter << endl;
        histCalcObs = new TH1D( Form("observables_%d", global_obs_hist_counter++), "observables"  //strAccumHistName , strAccumHistName //"histAccumulatedValues"
                                , nObs,-0.5, nObs-0.5 );
        for( int i=0; i < nObs; i++ )
        {
            mapObs.insert( pair<const char*,int>( obsNames[i], i ) );
            histCalcObs->GetXaxis()->SetBinLabel( i+1, obsNames[i] );
        }


        if ( (int)mapObs.size() < nObs )
        {
            cout << "AHTUNG!!! there are duplicates in the varNames array!!!" << endl;
            int aa;
            cin >> aa;
        }

    }


    // ###########
    void createNameValueMap() // const char* valueBinName, const char* nEvBinName )
    {
        //        TString strTmp;
        for( int i=0; i < histAccumulatedValues->GetNbinsX(); i++ )
        {
            double value = histAccumulatedValues->GetBinContent( i+1 );
            const char *strName = histAccumulatedValues->GetXaxis()->GetBinLabel(i+1);
            if ( value == 0 )
            {
                //                cout << ">>> warning: strName=" << strName << " = " << value << endl;
                continue;
            }

            mapData[ strName ] = value;
            //            cout << "test: "
            //                 << strName //histAccumulatedValues->GetXaxis()->GetBinLabel(i+1)
            //                 << " len=" << strlen(strName) << ", last char=" << char(strName[ strlen(strName) ])
            //                 << " " << mapData[ strName ] << ", COMPARE: " << strcmp( strName, "Nevents" )
            //                 << ", value = " << value << ", mapValue = " << mapData[ strName ] << endl;
        }

        //        cout << "!!! TEST nEv: " << strlen( "Nevents" ) << " "
        //             << mapData[ "Nevents" ] << endl; //"Nevents" ] << endl;
    }



    // ###########
    void fillHistWithValue( const char *varName, double value )
    {
        // with mapObs, it was wrong! (corrupted logic..)
//        int varId = mapObs[varName];
//        if ( varId != 0 || firstFill ) // to avoid filling if no such obs name in array (== map returns zero)
//        {
//            cout << "varName = " << varName << ", varId = " << varId << endl;
//            histCalcObs->Fill( varId,   value );
        histCalcObs->Fill( varName,   value );
//            if( firstFill )
//                firstFill = false;
//        }
    }



    // ###########
    void averageOverEvents( const char* valueBinName, const char* nEvBinName )
    {
        long int nEv = mapData[ nEvBinName ];

        if ( nEv == 0 )
        {
            return; //
            cout << "ABNORMAL! nEv==0, valueBinName=" << valueBinName << ", nEvBinName=" << nEvBinName << endl;
            int tmpA;
            cin >> tmpA;
        }

        mapData[ valueBinName ] /= nEv;
        //        data[binId] /= nEv;

        if ( mapData[ valueBinName ] == 0 )
        {
            return; //
            cout << "ABNORMAL! mapData[ valueBinName ]==0, valueBinName = " << valueBinName << endl;
            int tmpA;
            cin >> tmpA;
        }
    }

    // ###########
    void calc_bCorr( int corrTypeId, const char* nameF, const char* nameB, const char* nameF2,
                     const char* nameB2, const char* nameFB )
    {
        //        double F = data[id] ;
        //        double B = data[id+1] ;
        //        double F2 = data[id+2];
        //        double B2 = data[id+3];
        //        double FB = data[id+4];

        double F = mapData[ nameF ] ;
        double B = mapData[ nameB ] ;
        double F2 = mapData[ nameF2 ];
        double B2 = mapData[ nameB2 ];
        double FB = mapData[ nameFB ];


        double numerator = FB - F * B;
        double denominator_bCorr = F2 - F*F;
        //        double denominator_bCorr = meanB2 - meanB*meanB;
        //        double denominator_bCorr = sqrt( (X2 - X*X)*(Y2 - Y*Y) );
        //    if (type == 0) // bCorr
        //        denominator = meanF2 - meanF*meanF;
        //    else if (type == 1) // C2
        //        denominator = meanF * meanB;

        avF[corrTypeId] = F;
        avB[corrTypeId] = B;


        //bcorr
        bcorr[corrTypeId] = -1000;
        if ( denominator_bCorr != 0 )
            bcorr[corrTypeId] = numerator / denominator_bCorr;
        else
            if(0) cout << "!!!!! WARNING denominator_bCorr=" << denominator_bCorr
                 << ", F2=" << F2 << ", F*F=" << F*F << endl;

        if ( F != 0 && B != 0 )
        {
            nu_dyn[corrTypeId] = ( F2 - F ) / F / F
                    + ( B2 - B ) / B / B
                    - 2*( FB/F/B);
            sigma[corrTypeId] = nu_dyn[corrTypeId] / (1./F + 1./B) + 1.;
        }
    }

    // ###########
    void finalCalc( double eSep, bool if_Identical_FB_XY, double eSizeNum, double eSizeDenom )
    {
        // we assume that eSep=0 means windows are completely overlapped!
        bool identical = ( eSep==0 && if_Identical_FB_XY );
//        cout << "finalCalc: " << eSep << endl;

        histCalcObs->Reset();

        createNameValueMap();
        //        mapData[ "Nevents" ] = histAccumulatedValues->GetBinContent( mapVar["Nevents"]+1 );

        // ### averaging over events:
        averageOverEvents( "Nf",    "Nevents" );
        averageOverEvents( "Nb",    "Nevents" );
        averageOverEvents( "Nf2",  "Nevents" );
        averageOverEvents( "Nb2",  "Nevents" );
        averageOverEvents( "Nf*Nb", "Nevents" );

        averageOverEvents( "NfPb_Nf",     "b_Nevents" );
        averageOverEvents( "NfPb_avPb",     "b_Nevents" );
        averageOverEvents( "NfPb_Nf2",    "b_Nevents" );
        averageOverEvents( "NfPb_avPb2",    "b_Nevents" );
        averageOverEvents( "NfPb_Nf_avPb",  "b_Nevents" );

        averageOverEvents( "PfNb_avPf",     "f_Nevents" );
        averageOverEvents( "PfNb_Nb",     "f_Nevents" );
        averageOverEvents( "PfNb_avPf2",    "f_Nevents" );
        averageOverEvents( "PfNb_Nb2",    "f_Nevents" );
        averageOverEvents( "PfNb_avPf_Nb",  "f_Nevents" );


        averageOverEvents( "PfPb_avPf",     "fb_Nevents" );
        averageOverEvents( "PfPb_avPb",     "fb_Nevents" );
        averageOverEvents( "PfPb_avPf2",    "fb_Nevents" );
        averageOverEvents( "PfPb_avPb2",    "fb_Nevents" );
        averageOverEvents( "PfPb_avPf_avPb",  "fb_Nevents" );



        // x, y
        averageOverEvents( "Nx",     "Nevents" );
        averageOverEvents( "Ny",     "Nevents" );
        averageOverEvents( "Nx2",    "Nevents" );
        averageOverEvents( "Ny2",    "Nevents" );
        averageOverEvents( "Nx*Ny",  "Nevents" );


        averageOverEvents( "NxPy_Nx",     "y_Nevents" );
        averageOverEvents( "NxPy_Py",     "y_Nevents" );
        averageOverEvents( "NxPy_Nx2",    "y_Nevents" );
        averageOverEvents( "NxPy_Py2",    "y_Nevents" );
        averageOverEvents( "NxPy_Nx_Py",  "y_Nevents" );


        averageOverEvents( "PxNy_Px",     "x_Nevents" );
        averageOverEvents( "PxNy_Ny",     "x_Nevents" );
        averageOverEvents( "PxNy_Px2",    "x_Nevents" );
        averageOverEvents( "PxNy_Ny2",    "x_Nevents" );
        averageOverEvents( "PxNy_Px_Ny",  "x_Nevents" );


        averageOverEvents( "PxPy_Px",     "xy_Nevents" );
        averageOverEvents( "PxPy_Py",     "xy_Nevents" );
        averageOverEvents( "PxPy_Px2",    "xy_Nevents" );
        averageOverEvents( "PxPy_Py2",    "xy_Nevents" );
        averageOverEvents( "PxPy_Px_Py",  "xy_Nevents" );





        // ##### mixes:

        averageOverEvents( "Nf*Nx",  "Nevents" );
        averageOverEvents( "Nf*Ny",  "Nevents" );
        averageOverEvents( "Nb*Nx",  "Nevents" );
        averageOverEvents( "Nb*Ny",  "Nevents" );


        averageOverEvents( "Nf_Px",  "x_Nevents" );
        averageOverEvents( "Nf_Py",  "y_Nevents" );
        averageOverEvents( "Nb_Px",  "x_Nevents" );
        averageOverEvents( "Nb_Py",  "y_Nevents" );


        averageOverEvents( "Pf_Nx",  "f_Nevents" );
        averageOverEvents( "Pf_Ny",  "f_Nevents" );
        averageOverEvents( "Pb_Nx",  "b_Nevents" );
        averageOverEvents( "Pb_Ny",  "b_Nevents" );



        averageOverEvents( "Pf_Px",  "fx_Nevents" );
        averageOverEvents( "Pf_Py",  "fy_Nevents" );
        averageOverEvents( "Pb_Px",  "bx_Nevents" );
        averageOverEvents( "Pb_Py",  "by_Nevents" );


        // RATIOS:
        averageOverEvents( "Nf_OVER_Nb",  "b_Nevents" );
        averageOverEvents( "Nf_OVER_Nx",  "x_Nevents" );
        averageOverEvents( "Nf_OVER_Ny",  "y_Nevents" );


        averageOverEvents( "Nb_OVER_Nf",  "f_Nevents" );
        averageOverEvents( "Nb_OVER_Nx",  "x_Nevents" );
        averageOverEvents( "Nb_OVER_Ny",  "xy_Nevents" );

        averageOverEvents( "Nx_OVER_Nf",  "f_Nevents" );
        averageOverEvents( "Nx_OVER_Nb",  "b_Nevents" );
        averageOverEvents( "Nx_OVER_Ny",  "y_Nevents" );

        averageOverEvents( "Ny_OVER_Nf",  "f_Nevents" );
        averageOverEvents( "Ny_OVER_Nb",  "b_Nevents" );
        averageOverEvents( "Ny_OVER_Nx",  "x_Nevents" );

        // RATIOS - pT:
        averageOverEvents( "Nf_OVER_Nx_vs_Pb",  "bx_Nevents" );
        averageOverEvents( "Nf_OVER_Nx_vs_Py",  "xy_Nevents" );
        averageOverEvents( "Nb_OVER_Ny_vs_Pf",  "fy_Nevents" );
        averageOverEvents( "Nb_OVER_Ny_vs_Px",  "xy_Nevents" );


        averageOverEvents( "Nx_OVER_Nf_vs_Pb",  "fb_Nevents" );
        averageOverEvents( "Nx_OVER_Nf_vs_Py",  "fy_Nevents" );
        averageOverEvents( "Ny_OVER_Nb_vs_Pf",  "fb_Nevents" );
        averageOverEvents( "Ny_OVER_Nb_vs_Px",  "bx_Nevents" );


        averageOverEvents( "Nf_OVER_Nb_vs_Nx_OVER_Ny",  "by_Nevents" );
        averageOverEvents( "Nf_OVER_Nx_vs_Nb_OVER_Ny",  "xy_Nevents" );
        averageOverEvents( "Nb_OVER_Nf_vs_Ny_OVER_Nx",  "fx_Nevents" );
        averageOverEvents( "Nx_OVER_Nf_vs_Ny_OVER_Nb",  "fb_Nevents" );

        // new ratio-pT, April 2021
        averageOverEvents( "Nb_OVER_Ny*PF",  "y_Nevents" );
        averageOverEvents( "Nb_OVER_Ny*PX",  "y_Nevents" );
        averageOverEvents( "Nf_OVER_Nx*PB",  "x_Nevents" );
        averageOverEvents( "Nf_OVER_Nx*PY",  "x_Nevents" );

//        averageOverEvents( "nB*PF",  "Nevents" );
//        averageOverEvents( "nY*PF",  "Nevents" );
//        averageOverEvents( "nB*PX",  "Nevents" );
//        averageOverEvents( "nY*PX",  "Nevents" );


        averageOverEvents( "PxPy_avPx",  "Nevents" );
        averageOverEvents( "Nb_OVER_Ny_vs_avPx",  "xy_Nevents" );






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
        calc_bCorr( 0, "Nf", "Nb", "Nf2", "Nb2", "Nf*Nb" );
        calc_bCorr( 1, "NfPb_Nf", "NfPb_Pb", "NfPb_Nf2", "NfPb_Pb2", "NfPb_Nf_Pb" );
        calc_bCorr( 2, "PfPb_Pf", "PfPb_Pb", "PfPb_Pf2", "PfPb_Pb2", "PfPb_Pf_Pb" );
        //        calc_bCorr( 1, mapVar[ "NfPb_Nf" ] );
        //        calc_bCorr( 2, mapVar[ "PfPb_Pf" ] );


        fillHistWithValue( "bcorr_FB",      bcorr[0] );
        fillHistWithValue( "bcorr_PtN",     bcorr[1] );
        fillHistWithValue( "bcorr_PtPt",    bcorr[2] );
//        fillHistWithValue( "nu_dyn_FB",     nu_dyn[0] );
        fillHistWithValue( "nu_dyn_PtN",    nu_dyn[1] );
        fillHistWithValue( "nu_dyn_PtPt",   nu_dyn[2] );
//        fillHistWithValue( "sigma_FB",      sigma[0] );
        fillHistWithValue( "sigma_PtN",     sigma[1] );
        fillHistWithValue( "sigma_PtPt",    sigma[2] );


        fillHistWithValue( "avNF", avF[0] );
        fillHistWithValue( "avNB", avB[0] );

        fillHistWithValue( "avPtF", avF[2] );
        fillHistWithValue( "avPtB", avB[2] );




        // ##############################
        // now calc the observables

        double _nEvents     = mapData[ "Nevents" ];

        double FB = mapData[ "Nf*Nb" ];
        double XY = mapData[ "Nx*Ny" ];
        double FY = mapData[ "Nf*Ny" ];
        double XB = mapData[ "Nb*Nx" ];

        double F = mapData[ "Nf" ];
        double B = mapData[ "Nb" ];
        double X = mapData[ "Nx" ];
        double Y = mapData[ "Ny" ];

        double F2 = mapData[ "Nf2" ];
        double X2 = mapData[ "Nx2" ];
        double B2 = mapData[ "Nb2" ];
        double Y2 = mapData[ "Ny2" ];


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
        double rFB = mapData[ "Nf_OVER_Nx_vs_Nb_OVER_Ny" ];
        double ratioF = mapData[ "Nf_OVER_Nx" ];
        double ratioB = mapData[ "Nb_OVER_Ny" ];
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
//            double Nb_OVER_Ny_PX = mapData[ "Nb_OVER_Ny*PX" ];
            double PxPy_avPx      = mapData[ "PxPy_avPx" ];
            double Nb_OVER_Ny_vs_Px = mapData[ "Nb_OVER_Ny_vs_avPx" ];

            fillHistWithValue( "corr_rPt_direct", X/eSizeNum * ( Nb_OVER_Ny_vs_Px/PxPy_avPx/ratioB - 1 ) );

            // formula
            double nB_PX = mapData[ "nB*PX" ] / _nEvents;
            double nY_PX = mapData[ "nY*PX" ] / _nEvents;
            double PX      = mapData[ "PX" ]  / _nEvents;

            // new ratio-meanPt formula - July 2022: when expansion for <pT> is done as <nB*PX>/.. + <nY*nX>/.. - <nB*nX>/.. - <nY*PX>/..
//            double _avSumPtInEvX      = mapData[ "sumPtAllEvX" ]  / _nEvents;
            double rPt_full_formula = X/eSizeNum * (  nB_PX /B/PX  + XY/X/Y - XB/X/B - nY_PX /Y/PX  ); // here even if wins are identical, +1/X-1/X cancels!
            if (identical)
            {
//                double subtrX = identical ? 1/X : 0;
                double nX_PX = mapData[ "nX*PX" ] / _nEvents;
                rPt_full_formula = X/eSizeNum * (  nB_PX /B/PX  + X2/X/X - XB/X/B - nX_PX /X/PX );
            }
            fillHistWithValue( "corr_rPt_full_formula", rPt_full_formula );
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
            double PF_PB = mapData[ "PF*PB" ] / mapData[ "Nevents" ];
            double nF_PB = mapData[ "nF*PB" ] / mapData[ "Nevents" ];
            double nB_PF = mapData[ "nB*PF" ] / mapData[ "Nevents" ];

            double PF = mapData[ "PF" ] / mapData[ "Nevents" ];
            double PB = mapData[ "PB" ] / mapData[ "Nevents" ];


            double avPtF_avPtB_formula = F/eSizeNum * ( PF_PB/PF/PB + FB/F/B - nF_PB/F/PB  - nB_PF/B/PF );

            if (identical)
            {
                double PF2 = mapData[ "PF2" ] / mapData[ "Nevents" ];
                double nF_PF = mapData[ "nF*PF" ] / mapData[ "Nevents" ];
                double piF2 = mapData[ "piF2" ] / mapData[ "Nevents" ];

//                double subtrX = identical ? 1/X : 0;
                // DODELAT' PF_PB !!!!!
                avPtF_avPtB_formula = F/eSizeNum * ( (PF2 - piF2)/PF/PF + F2/F/F  - nF_PF/F/PF  - nF_PF/F/PF     - 1/F + 1/F + 1/F );
            }

            fillHistWithValue( "avPtF_avPtB_formula", avPtF_avPtB_formula );
        }

        // ##### coeff avPt-avPt direct
        {
            double Pf = mapData[ "PfPb_avPf" ] ;
            double Pb = mapData[ "PfPb_avPb" ] ;
            double Pf_Pb = mapData[ "PfPb_avPf_avPb" ];

            fillHistWithValue( "avPtF_avPtB_direct", F/eSizeNum * ( Pf_Pb/(Pf * Pb) - 1 ) );
//            cout << " >> Pf_Pb = " << Pf_Pb << ", Pf = " << Pf << ", Pb = " << Pb << endl;
        }







        // ##### dptdpt
        if(0)
        {
            double _avPtAllEvF     = mapData[ "sumPtAllEvF" ]  / ( F*_nEvents );
            double _avPtAllEvB     = mapData[ "sumPtAllEvB" ]  / ( B*_nEvents );
            double _piFpjB = mapData[ "piFpjB" ];
            //        double _nF_sum_pB = mapData[ "nF*sum_pB" ];
            //        double _nB_sum_pF = mapData[ "nB*sum_pF" ];
            double _nF_sum_pB = mapData[ "nF*PB" ]; // /2; // /2 is TMP!!!
            double _nB_sum_pF = mapData[ "nB*PF" ]; // /2;
            //        cout << " >> " << _piFpjB << " "<< _nF_sum_pB << " " << _nB_sum_pF << endl;
            //        cout << " >>>> " << ( _piFpjB - _nF_sum_pB*_avPtAllEvF - _nB_sum_pF*_avPtAllEvB)/ (FB*_nEvents) << endl;

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
                               (mapData[ "Nf" ]
                               *mapData[ "Nb" ]
                    *mapData[ "PF*PB" ]
                    -
                    mapData[ "Nb" ]
                    *mapData[ "PF" ]
                    *mapData[ "nF*PB" ]// /2 // /2 is TMP!!!
                    -
                    mapData[ "Nf" ]
                    *mapData[ "PB" ]
                    *mapData[ "nB*PF" ]// /2 // /2 is TMP!!!
                    +
                    mapData[ "PF" ]
                    *mapData[ "PB" ]
                    *mapData[ "Nf*Nb" ]
                    ) / (
                        mapData[ "PF" ]
                    *mapData[ "PB" ]
                    *mapData[ "Nf*Nb" ]
                    )
                    );
        }

        // ### coeff ptpt BOZEK_2017_USING_QM17_RESULTS_1704.02777.pdf
        if(0)
        {
            double _avPtAllEvF     = mapData[ "sumPtAllEvF" ]  / ( F*_nEvents );
            double _avPtAllEvB     = mapData[ "sumPtAllEvB" ]  / ( B*_nEvents );

            double Pf = mapData[ "PfPb_Pf" ] ;
            double Pb = mapData[ "PfPb_Pb" ] ;
            double Pf_Pb = mapData[ "PfPb_Pf_Pb" ];

            double _pipjF = mapData[ "pipjF" ];
            double _pipjB = mapData[ "pipjB" ];

            double secondTermF     = mapData[ "(nF-1)*PF" ];
            double secondTermB     = mapData[ "(nB-1)*PB" ];

            double nn_1F     = mapData[ "nF*(nF-1)" ];
            double nn_1B     = mapData[ "nB*(nB-1)" ];


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


    double getValue( TString strBinName ) //char *strBinName )
    {
        //        cout << "histCalcObs: " << histCalcObs << endl;
        //        cout << histCalcObs->GetName() << endl;
        //        cout << histCalcObs->GetXaxis()->FindBin( strBinName ) << endl;
        return histCalcObs->GetBinContent( histCalcObs->GetXaxis()->FindBin( strBinName ) );
        //        return histCalcObs->GetBinContent( 1 );
    }

    void getValueWithError( TGraphErrors *gr, int pointId, double x, TString strBinName ) //char *strBinName )
    {
        //        cout << "histCalcObs: " << histCalcObs << endl;
        //        cout << histCalcObs->GetName() << endl;
        //        cout << histCalcObs->GetXaxis()->FindBin( strBinName ) << endl;
        double y = histCalcObs->GetBinContent( histCalcObs->GetXaxis()->FindBin( strBinName ) );
//        if( strBinName.CompareTo( "sigma_FB" ) )
//            cout << "y = " << y << endl;
        double yerr = histCalcObs->GetBinError( histCalcObs->GetXaxis()->FindBin( strBinName ) );
        //        return histCalcObs->GetBinContent( 1 );
        gr->SetPoint( pointId, x, y );
        gr->SetPointError( pointId, 0, yerr );
    }
};

int SimpleCalculations::global_obs_hist_counter = 0;


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
    TH3D *h3D;
    SimpleCalculations **wpObs;
//    double eSep[MAX_N_WIN_PAIRS];


    double takeEtaSep( TH3D *_h3D, int _iEta )
    {
        TString strEtaBin = _h3D->GetYaxis()->GetBinLabel( _iEta + 1 );
//        cout << "strEtaBin: " << strEtaBin << endl;

        if ( strEtaBin.Contains("eB") ) // i.e. we do pair-by-pair hist3D
        {
            float eBounds[4];
            sscanf( strEtaBin.Data(), "eB_%f_%f_eF_%f_%f", &eBounds[0], &eBounds[1], &eBounds[2], &eBounds[3] );
            if(0)for ( int k = 0; k < 4; k++ )
            {
//                    eBounds[k] = round( eBounds[k]*100 ) / 100;
                cout << "eBounds[" << k << "] = " << eBounds[k]  << endl;
            }
            double eBpos = ( eBounds[0] + eBounds[1] )/2;
            double eFpos = ( eBounds[2] + eBounds[3] )/2;
            return round( (eFpos - eBpos)*100 ) / 100;
        }
        else
            return round( _h3D->GetYaxis()->GetBinCenter( _iEta+1 ) *100 ) / 100;
    }

    void calc(  double eSizeNum, double eSizeDenom, bool if_Identical_FB_XY )
    {
        int nEta = h3D->GetNbinsY();
        int nSubs = h3D->GetNbinsZ();
        //        int nObs = h3D->GetNbinsX();

        wpObs = new SimpleCalculations*[nEta];

        cout << "##### nObs = " << nObs << endl;
        SimpleCalculations *wpSubsamples = new SimpleCalculations[nSubs];

        for ( int iEta = 0; iEta < h3D->GetNbinsY(); iEta++ )
        {
//            cout << "###### STARTING iEta = " << iEta << endl;
            wpObs[iEta] = new SimpleCalculations;
            // continue;


            // mean values
            wpObs[iEta]->histAccumulatedValues = h3D->ProjectionX( "_px", iEta+1, iEta+1, 1, nSubs ); // project all subsamples on axis

            // find out eSep:
            double eSep = takeEtaSep( h3D, iEta );
//            cout << "eSep = " << eSep << endl;

            wpObs[iEta]->finalCalc( eSep, if_Identical_FB_XY, eSizeNum, eSizeDenom );


            //            cout << "wpObs[iEta].getValue( corr_rr_formula): " << wpObs[iEta]->getValue( "corr_rr_formula")  << endl;

            // subsamples
            for ( int iSub = 0; iSub < nSubs; iSub++ )
            {
                //                cout << "iSub=" << iSub << ", iType=" << iType << ", iCW=" << iCW << ", iEta=" << iEta << endl;
                //                cout << "iSub=" << iSub << ", iEta=" << iEta << endl;

                SimpleCalculations *wp = &wpSubsamples[iSub];

                wp->histAccumulatedValues = h3D->ProjectionX( "_px", iEta+1, iEta+1, iSub+1, iSub+1 );
                wp->finalCalc( eSep, if_Identical_FB_XY, eSizeNum, eSizeDenom );
            } // end of subsamples


            // subsampling for each observable:
            for( int bin = 0; bin < nObs; bin++ )
            {
                // mean
                double mean = wpObs[iEta]->histCalcObs->GetBinContent( bin+1);


                // stdDev:
                double std_dev = 0;
                for ( int iSub = 0; iSub < nSubs; iSub++)
                {
                    double value = wpSubsamples[iSub].histCalcObs->GetBinContent( bin+1);
                    float diff = value - mean;
                    std_dev += diff*diff;
                }
                std_dev = sqrt(std_dev / nSubs / (nSubs-1) );

                wpObs[iEta]->histCalcObs->SetBinError( bin+1, std_dev );

                //                cout << wpObs[iEta]->histCalcObs->GetXaxis()->GetBinLabel( bin+1 ) << " = " << mean << ", std_dev = " << std_dev << endl;

            }
        }
        delete [] wpSubsamples;
    }
};







#endif











