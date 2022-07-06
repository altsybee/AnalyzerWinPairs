#include "TH1.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TDirectory.h"
#include "TString.h"
#include "TGraphErrors.h"

//#include "SimpleWinPair_forPythiaAnalysis_v3.h"
//#include "WinPair2019.h"
#include "WinPairBase.h"



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
    "bcorr_FB",
    "bcorr_PtN",
    "bcorr_PtPt",

    "nu_dyn_NN",
    "nu_dyn_PtN",
    "nu_dyn_PtPt",

    "sigma_FB",
    "sigma_PtN",
    "sigma_PtPt",

    "avNF",
    "avNB",

    "avPtF",
    "avPtB",

    // Dec2019: bonus:
    "sigma_approx_FB",

    "bcorr_XY",
    "nu_dyn_XY",
    "sigma_XY",
    "sigma_approx_XY",

    "nu_dyn_FY",
    "sigma_FY",
    "sigma_approx_FY",

    "nu_dyn_XB",
    "sigma_XB",
    "sigma_approx_XB",


    "avX",
    "avY",



    "corr_rr_formula",
    "corr_rr_direct",

    "corr_R2_aa",
    "corr_R2_bb",
    "corr_R2_ab",
    "corr_R2_ba",


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
    "coeff_ptpt_with_means_in_denom",
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

    // final quantities
    double bcorr[3];
    double nu_dyn[3];
    double sigma[3];

    double avF[3];
    double avB[3];


    TH1D *histAccumulatedValues;   // //! accumulated values for observables
    TH1D *histCalcObs;   // //!   hist with calculated observables

    SimpleCalculations()
    {
        histAccumulatedValues = 0x0;

        //        cout << "global_obs_hist_counter = " << global_obs_hist_counter << endl;
        histCalcObs = new TH1D( Form("observables_%d", global_obs_hist_counter++), "observables"  //strAccumHistName , strAccumHistName //"histAccumulatedValues"
                                , nObs,-0.5, nObs-0.5 );
        for( int i=0; i < nObs; i++ )
            histCalcObs->GetXaxis()->SetBinLabel( i+1, obsNames[i] );

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
            cout << "!!!!! WARNING denominator_bCorr=" << denominator_bCorr
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
    void finalCalc( double eSizeNum, double eSizeDenom )
    {
        histCalcObs->Reset();

        createNameValueMap();
        //        mapData[ "Nevents" ] = histAccumulatedValues->GetBinContent( mapVar["Nevents"]+1 );

        // ### averaging over events:
        averageOverEvents( "Nf",    "Nevents" );
        averageOverEvents( "Nb",    "Nevents" );
        averageOverEvents( "N2_f",  "Nevents" );
        averageOverEvents( "N2_b",  "Nevents" );
        averageOverEvents( "Nf_Nb", "Nevents" );

        averageOverEvents( "NfPb_Nf",     "b_Nevents" );
        averageOverEvents( "NfPb_Pb",     "b_Nevents" );
        averageOverEvents( "NfPb_Nf2",    "b_Nevents" );
        averageOverEvents( "NfPb_Pb2",    "b_Nevents" );
        averageOverEvents( "NfPb_Nf_Pb",  "b_Nevents" );

        averageOverEvents( "PfNb_Pf",     "f_Nevents" );
        averageOverEvents( "PfNb_Nb",     "f_Nevents" );
        averageOverEvents( "PfNb_Pf2",    "f_Nevents" );
        averageOverEvents( "PfNb_Nb2",    "f_Nevents" );
        averageOverEvents( "PfNb_Pf_Nb",  "f_Nevents" );


        averageOverEvents( "PfPb_Pf",     "fb_Nevents" );
        averageOverEvents( "PfPb_Pb",     "fb_Nevents" );
        averageOverEvents( "PfPb_Pf2",    "fb_Nevents" );
        averageOverEvents( "PfPb_Pb2",    "fb_Nevents" );
        averageOverEvents( "PfPb_Pf_Pb",  "fb_Nevents" );
        // KOSTYL'!!!!
        //        averageOverEvents( "PfPb_Pf",     "Nevents" );
        //        averageOverEvents( "PfPb_Pb",     "Nevents" );
        //        averageOverEvents( "PfPb_Pf2",    "Nevents" );
        //        averageOverEvents( "PfPb_Pb2",    "Nevents" );
        //        averageOverEvents( "PfPb_Pf_Pb",  "Nevents" );



        // x, y
        averageOverEvents( "Nx",     "Nevents" );
        averageOverEvents( "Ny",     "Nevents" );
        averageOverEvents( "Nx2",    "Nevents" );
        averageOverEvents( "Ny2",    "Nevents" );
        averageOverEvents( "Nx_Ny",  "Nevents" );


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

        averageOverEvents( "Nf_Nx",  "Nevents" );
        averageOverEvents( "Nf_Ny",  "Nevents" );
        averageOverEvents( "Nb_Nx",  "Nevents" );
        averageOverEvents( "Nb_Ny",  "Nevents" );


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

        averageOverEvents( "nB*PF",  "Nevents" );
        averageOverEvents( "nY*PF",  "Nevents" );
        averageOverEvents( "nB*PX",  "Nevents" );
        averageOverEvents( "nY*PX",  "Nevents" );



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





        // ### !!! check that all data have been writen successfully:
        //        for( int i = 11; i < nVars; i++ )
        //            if ( data[i] == -1000 )
        //            {
        ////                cout << "ABNORMAL! data[i] == -1000 for i = " << i << ", varName = " << varNames[i] << endl;
        ////                int tmpA;
        ////                cin >> tmpA;
        //            }



        //        averageOverEvents( 0 /*corr type*/, "Nevents" );
        //        averageOverEvents( 1 /*corr type*/, "b_Nevents" );
        //        averageOverEvents( 2 /*corr type*/, "fb_Nevents" );

        //        cout << "mapVar[  Nf  ] = " << mapVar[ "Nf" ]  << endl;

        //        int aa;
        //        cin >> aa;
        //        calc_bCorr( 0, "Nx", "Ny", "Nx2", "Ny2", "Nx_Ny" );
        calc_bCorr( 0, "Nf", "Nb", "N2_f", "N2_b", "Nf_Nb" );
        calc_bCorr( 1, "NfPb_Nf", "NfPb_Pb", "NfPb_Nf2", "NfPb_Pb2", "NfPb_Nf_Pb" );
        calc_bCorr( 2, "PfPb_Pf", "PfPb_Pb", "PfPb_Pf2", "PfPb_Pb2", "PfPb_Pf_Pb" );
        //        calc_bCorr( 1, mapVar[ "NfPb_Nf" ] );
        //        calc_bCorr( 2, mapVar[ "PfPb_Pf" ] );


        histCalcObs->Fill( "bcorr_FB",      bcorr[0] );
        histCalcObs->Fill( "bcorr_PtN",     bcorr[1] );
        histCalcObs->Fill( "bcorr_PtPt",    bcorr[2] );
        histCalcObs->Fill( "nu_dyn_NN",     nu_dyn[0] );
        histCalcObs->Fill( "nu_dyn_PtN",    nu_dyn[1] );
        histCalcObs->Fill( "nu_dyn_PtPt",   nu_dyn[2] );
        histCalcObs->Fill( "sigma_FB",      sigma[0] );
        histCalcObs->Fill( "sigma_PtN",     sigma[1] );
        histCalcObs->Fill( "sigma_PtPt",    sigma[2] );


        histCalcObs->Fill( "avNF", avF[0] );
        histCalcObs->Fill( "avNB", avB[0] );

        histCalcObs->Fill( "avPtF", avF[2] );
        histCalcObs->Fill( "avPtB", avB[2] );





        // ##############################
        // ##### coeff ratio-ratio
        // ##############################

        // formula:
        double FB = mapData[ "Nf_Nb" ];
        double XY = mapData[ "Nx_Ny" ];
        double FY = mapData[ "Nf_Ny" ];
        double XB = mapData[ "Nb_Nx" ];

        double F = mapData[ "Nf" ];
        double B = mapData[ "Nb" ];
        double X = mapData[ "Nx" ];
        double Y = mapData[ "Ny" ];

        histCalcObs->Fill( "corr_R2_aa", FB/F/B - 1 );
        histCalcObs->Fill( "corr_R2_bb", XY/X/Y - 1 );
        histCalcObs->Fill( "corr_R2_ab", FY/F/Y - 1 );
        histCalcObs->Fill( "corr_R2_ba", XB/X/B - 1 );

        //        histCalcObs->Fill( "corr_rr_formula", (F+B)/2 * (FB/F/B + XY/X/Y - FY/F/Y - XB/X/B) );
        //        histCalcObs->Fill( "corr_rr_formula", F*B/(F+B) * (FB/F/B + XY/X/Y - FY/F/Y - XB/X/B) );
        histCalcObs->Fill( "corr_rr_formula", F/eSizeNum * X/eSizeDenom / (F/eSizeNum + X/eSizeDenom) * (FB/F/B + XY/X/Y - FY/F/Y - XB/X/B) );


        if(0)
        {
            cout << "F=" << F << ", B=" << B;
            cout << ", X=" << X << ", Y=" << Y;
            cout << ", FB=" << FB << "    XY=" << XY << ", FY=" << FY <<  ", XB=" << XB;
            cout << ", bcorr_ratio=" << FB/F/B + XY/X/Y - FY/F/Y - XB/X/B << endl;
        }

        // direct:
        double rFB = mapData[ "Nf_OVER_Nx_vs_Nb_OVER_Ny" ];
        double ratioF = mapData[ "Nf_OVER_Nx" ];
        double ratioB = mapData[ "Nb_OVER_Ny" ];

        //        histCalcObs->Fill( "corr_rr_direct", (F+B)/2 * ( rFB/ratioF/ratioB - 1 ) );
        //        histCalcObs->Fill( "corr_rr_direct", F*B/(F+B) * ( rFB/ratioF/ratioB - 1 ) );
        histCalcObs->Fill( "corr_rr_direct", F/eSizeNum * X/eSizeDenom / (F/eSizeNum + X/eSizeDenom) * ( rFB/ratioF/ratioB - 1 ) );

        //        cout <<

        // new ratio-ratio - June 2021: when denominator is in full acceptance
        {
            double X2 = mapData[ "Nx2" ];
            double FX = FY;
            //            histCalcObs->Fill( "corr_rr_FULL_ETA_DENOM_formula", (F+B)/2 * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X) );
            //            histCalcObs->Fill( "corr_rr_FULL_ETA_DENOM_formula", F*B/(F+B) * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X) );
            histCalcObs->Fill( "corr_rr_FULL_ETA_DENOM_formula", F/eSizeNum * X/eSizeDenom / ( F/eSizeNum + X/eSizeDenom ) * (FB/F/B + X2/X/X - FX/F/X - XB/X/B    -1/X) );


            histCalcObs->Fill( "corr_rr_FULL_ETA_DENOM_direct", F/eSizeNum * X/eSizeDenom / ( F/eSizeNum + X/eSizeDenom ) * ( rFB/ratioF/ratioB - 1   -1/X) );

        }

        //        cout << ">>> F/eSizeNum = " << F/eSizeNum << ", X/eSizeNum = " << X/eSizeNum << ", X/eSizeDenom = " << X/eSizeDenom << endl;
        //        cout << ">>> prefactor in corr_rr: " << F/eSizeNum * X/eSizeNum / (F/eSizeNum + X/eSizeNum) << "  ";
        //        cout << ">>> prefactor in corr_rr_FULL_ETA_DENOM: " << F/eSizeNum * X/eSizeDenom / ( F/eSizeNum + X/eSizeDenom ) << endl;





        // ##############################
        // ##### coeff ratio-pT
        // ##############################

        // formula:
        //        double pXB = mapData[ "Nb_Px" ];
        //        double pXY = mapData[ "PxNy_Px_Ny" ];
        //        double ptX = mapData[ "PxNy_Px"  ];         // TAKE pT for x!!!

        //        corr_rPt_formula = B*( pXB/ptX/B - pXY/ptX/Y);

        double pYF = mapData[ "Nf_Py"       ];
        double pYX = mapData[ "NxPy_Nx_Py"  ];
        double ptY = mapData[ "NxPy_Py"     ];         // TAKE pT for x!!!

        // April 2021:
        double pFB = mapData[ "PfNb_Pf_Nb" ];
        double pFY = mapData[ "Pf_Ny" ];
        //        double ptB = mapData[ "NfPb_Pb"  ];
        double ptF = mapData[ "PfNb_Pf"  ];



        //        corr_rPt_formula = F*( pYF/ptY/F - pYX/ptY/X);
        histCalcObs->Fill( "corr_rPt_formula", X/eSizeNum * ( pFB/ptF/B - pFY/ptF/Y) );


        if(0)cout << " >> pYF=" << pYF << ", pFB=" << pFB
                  << " >> ptY=" << ptY << ", ptF=" << ptF
                  << " >> pYX=" << pYX << ", pFY=" << pFY
                  << endl;


        double eff = 1;
        //        double eff = 0.8;

        // direct:
        //        double BoverY_vs_Px = mapData[ "Nb_OVER_Ny_vs_Px" ]; // TAKE pT for x!!!
        double FoverX_vs_Pb = mapData[ "Nf_OVER_Nx_vs_Pb" ]; //"Nb_OVER_Ny_vs_Px" ]; // TAKE pT for x!!!
        double BoverY_vs_Pf = mapData[ "Nb_OVER_Ny_vs_Pf" ];

        //        corr_rPt_direct = B*( BoverY_vs_Px/ptX/ratioB - 1 );
        //        corr_rPt_direct = F*( FoverX_vs_Pb/ptY/ratioF - 1 );
        histCalcObs->Fill( "corr_rPt_direct", X/eSizeNum * ( BoverY_vs_Pf/ptF/ratioB - 1 ) );
//        cout << " BoverY_vs_Pf = " << BoverY_vs_Pf << ", ratioB = " << ratioB << ", ptF = " << ptF << ", X/eSizeNum = " << X/eSizeNum << endl;

        // new ratio-pt - April 2021: when pt is not averaged in each event
        double _nEvents     = mapData[ "Nevents" ];
        double _avPtAllEvF      = mapData[ "sumPtAllEvF" ]  / ( F*_nEvents );
        //        double _avPtAllEvB     = mapData[ "sumPtAllEvB" ]  / ( B*_nEvents );
        double _avPtAllEvX      = mapData[ "sumPtAllEvX" ]  / ( X*_nEvents );
        {
            // B/Y vs pT_F
            double av_nB_PF = mapData[ "nB*PF" ]; // / _nEvents;
            double av_nY_PF = mapData[ "nY*PF" ]; // / _nEvents;
            double Nb_OVER_Ny_PF = mapData[ "Nb_OVER_Ny*PF" ];

            //            corr_rSumPt_formula = F/eff*( av_nB_PF /B/F/_avPtAllEvF - av_nY_PF /Y/F/_avPtAllEvF );//- 1 );
            //            corr_rSumPt_direct  = F/eff*( Nb_OVER_Ny_PF/F/ratioB/_avPtAllEvF - 1 );

            // B/Y vs pT_X
            double av_nB_PX = mapData[ "nB*PX" ];   // / _nEvents;
            double av_nY_PX = mapData[ "nY*PX" ];   // / _nEvents;
            double Nb_OVER_Ny_PX = mapData[ "Nb_OVER_Ny*PX" ];

            histCalcObs->Fill( "corr_rSumPt_formula", X/eSizeNum */*X/eff**/( av_nB_PX /B/X/_avPtAllEvX - av_nY_PX /Y/X/_avPtAllEvX ) );//- 1 );
            histCalcObs->Fill( "corr_rSumPt_direct", X/eSizeNum */*X/eff**/( Nb_OVER_Ny_PX/X/ratioB/_avPtAllEvX - 1 ) );


            // new ratio-meanPt - July 2022: when expansion for <pT> is done as <nB*PX>/.. + <nY*nX>/.. - <nB*nX>/.. - <nY*PX>/..
            double _avSumPtInEvX      = mapData[ "sumPtAllEvX" ]  / _nEvents;
            histCalcObs->Fill( "corr_rPt_full_formula", X/eSizeNum * (  av_nB_PX /B/_avSumPtInEvX  + XY/X/Y - XB/X/B - av_nY_PX /Y/_avSumPtInEvX ) );


            if(0)
            {
                cout << " >> _avPtAllEvF=" << _avPtAllEvF << ", _avPtAllEvX=" << _avPtAllEvX << endl;
                cout << " >> Nb_OVER_Ny_PF=" << Nb_OVER_Ny_PF << ", Nb_OVER_Ny_PX=" << Nb_OVER_Ny_PX
                     << " >> F=" << F << ", X=" << X
                     << " >> ratioB=" << ratioB << ", ratioB=" << ratioB
                     << endl;
            }
        }

        if(0)
        {
            //            cout << ">>> pXB=" << pXB << ", pXY=" << pXY;
            //            cout << ", ptX=" << ptX;
            cout << ">>> pYF=" << pYF << ", pYX=" << pYX;
            cout << ", ptY=" << ptY;
            //            cout << " >> corr_rPt_formula=" << corr_rPt_formula << ", corr_rPt_direct=" << corr_rPt_direct << endl;
        }
        //        cout << ">>> pYF=" << pYF << ", pYX=" << pYX;
        //        cout << ", ptY=" << ptY;
        //        cout << " >> FoverX_vs_Pb=" << FoverX_vs_Pb << ", BoverY_vs_Pf=" << BoverY_vs_Pf
        //             << " >> ptY=" << ptY << ", ptF=" << ptF
        //             << " >> ratioF=" << ratioF << ", ratioB=" << ratioB
        //             << endl;


        //        cout << " >> corr_rPt_formula=" << corr_rPt_formula << ", corr_rPt_direct=" << corr_rPt_direct << endl;




        // nn for X, Y:
        if(1)
        {
            double X2 = mapData[ "Nx2" ];
            double Y2 = mapData[ "Ny2" ];

            double F2 = mapData[ "N2_f" ];
            double B2 = mapData[ "N2_b" ];


            double numerator = XY - X * Y;
            double denominator_bCorr = X2 - X*X;

            histCalcObs->Fill( "avX", F );
            histCalcObs->Fill( "avY", B );


            //bcorr
            histCalcObs->Fill( "bcorr_XY", -1000 );
            if ( denominator_bCorr != 0 )
                histCalcObs->Fill( "bcorr_XY", numerator / denominator_bCorr );
            else
                cout << "!!!!! WARNING denominator_bCorr=" << denominator_bCorr
                     << ", X2=" << X2 << ", X*X=" << X*X << endl;

            if ( X != 0 && Y != 0 )
            {
                histCalcObs->Fill( "sigma_approx_FB", ( 2 - 2*( FB/F/B) ) / (1./F + 1./B) + 1. );

                // XY
                double nu_dyn_XY = ( X2 - X ) / X / X
                        + ( Y2 - Y ) / Y / Y
                        - 2*( XY/X/Y);
                histCalcObs->Fill( "nu_dyn_XY", nu_dyn_XY );
                histCalcObs->Fill( "sigma_XY", nu_dyn_XY / (1./X + 1./Y) + 1. );
                histCalcObs->Fill( "sigma_approx_XY", ( 2 - 2*( XY/X/Y) ) / (1./X + 1./Y) + 1. );

                // FY
                double nu_dyn_FY = ( F2 - F ) / F / F
                        + ( Y2 - Y ) / Y / Y
                        - 2*( FY/F/Y);
                histCalcObs->Fill( "nu_dyn_FY", nu_dyn_FY );
                histCalcObs->Fill( "sigma_FY", nu_dyn_FY / (1./F + 1./Y) + 1.);
                histCalcObs->Fill( "sigma_approx_FY", ( 2 - 2*( FY/F/Y) ) / (1./F + 1./Y) + 1. );

                // XB
                double nu_dyn_XB = ( X2 - X ) / X / X
                        + ( B2 - B ) / B / B
                        - 2*( XB/X/B);
                histCalcObs->Fill( "nu_dyn_XB", nu_dyn_XB );
                histCalcObs->Fill( "sigma_XB", nu_dyn_XB / (1./X + 1./B) + 1. );
                histCalcObs->Fill( "sigma_approx_XB", ( 2 - 2*( XB/X/B) ) / (1./X + 1./B) + 1. );
            }
        }




        // ##### coeff ratio-sum pT (April 2021)
        {

        }







        // ##### dptdpt
        //        {
        //            double _nEvents     = mapData[ "Nevents" ];
        //            double _avPtAllEvF     = mapData[ "sumPtAllEvF" ]  / ( F*_nEvents );
        double _avPtAllEvB     = mapData[ "sumPtAllEvB" ]  / ( B*_nEvents );
        if(0)cout << F << " "<< B << " " << _nEvents << " " << _avPtAllEvF << " " << _avPtAllEvB << endl;

        double _piFpjB = mapData[ "piFpjB" ];
        //        double _nF_sum_pB = mapData[ "nF*sum_pB" ];
        //        double _nB_sum_pF = mapData[ "nB*sum_pF" ];
        double _nF_sum_pB = mapData[ "nF*PB" ]; // /2; // /2 is TMP!!!
        double _nB_sum_pF = mapData[ "nB*PF" ]; // /2;
        //        cout << " >> " << _piFpjB << " "<< _nF_sum_pB << " " << _nB_sum_pF << endl;
        //        cout << " >>>> " << ( _piFpjB - _nF_sum_pB*_avPtAllEvF - _nB_sum_pF*_avPtAllEvB)/ (FB*_nEvents) << endl;

        // from /Volumes/OptibaySSD/ALICE_analysis/AliceTaskGetEventTreeIA/task_FB_and_DptDpt_analysis/results/_common_files/routine.h:169
        histCalcObs->Fill( "coeff_dptdpt", (
                               ( _piFpjB - _nF_sum_pB*_avPtAllEvF - _nB_sum_pF*_avPtAllEvB)/ (FB*_nEvents)
                               + _avPtAllEvF*_avPtAllEvB
                               ) / (_avPtAllEvF*_avPtAllEvB) );

        //        }
        // ### coeff ptpt with means in denom
        {
            double Pf = mapData[ "PfPb_Pf" ] ;
            double Pb = mapData[ "PfPb_Pb" ] ;
            double Pf_Pb = mapData[ "PfPb_Pf_Pb" ];

            histCalcObs->Fill( "coeff_ptpt_with_means_in_denom",  Pf_Pb/(Pf * Pb) - 1 );
        }


        // ### coeff ptpt with means in denom
        {

            histCalcObs->Fill( "coeff_ptpt_CHECK_VV_formula",
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
                    *mapData[ "Nf_Nb" ]
                    ) / (
                        mapData[ "PF" ]
                    *mapData[ "PB" ]
                    *mapData[ "Nf_Nb" ]
                    )
                    );
        }

        // ### coeff ptpt BOZEK_2017_USING_QM17_RESULTS_1704.02777.pdf
        {
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

            histCalcObs->Fill( "coeff_ptpt_BOZEK_2017", (Pf_Pb - Pf*Pb)/sqrt( C_Bozek_F_nPairs_OUTSIDE_sum * C_Bozek_B_nPairs_OUTSIDE_sum ) );

            //            if(0)cout <<  ">>> C_Bozek_B_nPairs_OUTSIDE_sum: " << _pipjB << " "<< _avPtAllEvB << " " << _avPtAllEvB << " " << nn_1B << ", coeff_ptpt_BOZEK_2017=" << coeff_ptpt_BOZEK_2017 << endl;

            histCalcObs->Fill( "C_BOZEK_F", C_Bozek_F_nPairs_OUTSIDE_sum );
            histCalcObs->Fill( "C_BOZEK_B", C_Bozek_B_nPairs_OUTSIDE_sum );
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
        double yerr = histCalcObs->GetBinError( histCalcObs->GetXaxis()->FindBin( strBinName ) );
        //        return histCalcObs->GetBinContent( 1 );
        gr->SetPoint( pointId, x, y );
        gr->SetPointError( pointId, 0, yerr );
    }
};


int SimpleCalculations::global_obs_hist_counter = 0;




struct CalcWithSubsamples
{
    TH3D *h3D;
    SimpleCalculations **wpObs;

    void calc(  double eSizeNum, double eSizeDenom )
    {
        int nEta = h3D->GetNbinsY();
        int nSubs = h3D->GetNbinsZ();
        //        int nObs = h3D->GetNbinsX();

        wpObs = new SimpleCalculations*[nEta];

        cout << "##### nObs = " << nObs << endl;
        SimpleCalculations *wpSubsamples = new SimpleCalculations[nSubs];

        for ( int iEta = 0; iEta < h3D->GetNbinsY(); iEta++ )
        {
            cout << "###### STARTING iEta=" << iEta << endl;
            wpObs[iEta] = new SimpleCalculations;
            // continue;


            // mean values
            wpObs[iEta]->histAccumulatedValues = h3D->ProjectionX("_px", iEta+1, iEta+1, 1, nSubs ); // project all subsamples on axis
            wpObs[iEta]->finalCalc( eSizeNum, eSizeDenom );

            //            cout << "wpObs[iEta].getValue( corr_rr_formula): " << wpObs[iEta]->getValue( "corr_rr_formula")  << endl;

            // subsamples
            for ( int iSub = 0; iSub < nSubs; iSub++ )
            {
                //                cout << "iSub=" << iSub << ", iType=" << iType << ", iCW=" << iCW << ", iEta=" << iEta << endl;
                //                cout << "iSub=" << iSub << ", iEta=" << iEta << endl;

                SimpleCalculations *wp = &wpSubsamples[iSub];

                wp->histAccumulatedValues = h3D->ProjectionX( "_px", iEta+1, iEta+1, iSub+1, iSub+1 );
                wp->finalCalc( eSizeNum, eSizeDenom );
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


















