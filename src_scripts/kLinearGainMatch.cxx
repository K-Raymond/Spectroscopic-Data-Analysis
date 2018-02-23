// To compile:
// Note: GRSISort looks for .cxx extensions when compiling (for example it looks
// in the myAnalysis directory)
// Alternatively you may use the following to compile:
// g++ GainMatch.C -o MyGainMatch -std=c++0x -I$GRSISYS/GRSISort/include/
// `grsi-config --cflags --all-libs --root`

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <math.h> // round, floor, ceil, trunc
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "THStack.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
//#include "TFileIter.h"
#include "TApplication.h"
#include "TF1.h"
#include "TGainMatch.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TPPG.h"
#include "TPeak.h"
#include "TROOT.h"
#include "TScaler.h"
#include "TSpectrum.h"
#include "TStyle.h"
#include "TTree.h"

#ifndef __CINT__
#include "TGriffin.h"
#include "TSceptar.h"
#endif

///// GLOBAL VARIABLES /////

// When using fragment trees, set gIsCalibration = 1
// Otherwise gIsCalibration = 0
const bool gIsCalibration = 0;

// Input the two peaks for fitting
//  {peak, fit window}
const double_t gCalPeaks[2][2] = { {315.42, 20}, {1864.89, 20} };

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage is: %s <fragment or analysis tree file> ).\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    TFile *pFile = new TFile(argv[1]);

    if( pFile == nullptr ) {
        printf("Failed to open file '%s'!\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    TTree *pTree = nullptr;
    if( gIsCalibration )
        pTree = static_cast<TTree*>( pFile->Get("FragmentTree") );
    else
        pTree = static_cast<TTree*>( pFile->Get("AnalysisTree") );

    if (!pTree) {
        printf("Failed to find TTree\n");
        exit(EXIT_FAILURE);
    }

    double_t MeasuredPeaks[2];
    
    TChannel* pChannel = nullptr;
    TChannel::ReadCalFromTree(pTree);

    // Keep statistics at 1 keV/bin
    TH2D *mat_en = new TH2D("mat_en", "", 64, 0, 64, 5000, 0, 5000);

    if (gIsCalibration)
        pTree->Project("mat_en",
                        "TFragment.GetEnergy():TFragment.GetChannelNumber()");
	else
        pTree->Project("mat_en", "TGriffin.fGriffinLowGainHits.GetEnergy():"
                                  "TGriffin.fGriffinLowGainHits.GetChannel()."
                                  "fNumber"); // for cal files
    for (int i = 0; i < 64 ; i++ ) {
        pChannel = TChannel::GetChannelByNumber(i);
        TH1D *h_en = mat_en->ProjectionY(Form("h_%.2i", i), i + 1, i + 1);

        for ( int k = 0; k < 2 ; k++ ) {
            TSpectrum s;
            h_en->GetXaxis()->SetRangeUser( 
                    gCalPeaks[k][0] - gCalPeaks[k][1],
                    gCalPeaks[k][0] + gCalPeaks[k][1] );
            s.Search(h_en, 2, "", 0.25);
            double_t SpecPeak = s.GetPositionX()[0];
            h_en->GetXaxis()->UnZoom();

            TPeak* TempP = new TPeak(SpecPeak, SpecPeak - gCalPeaks[k][1],
                    SpecPeak + gCalPeaks[k][1] );

            TempP->Fit(h_en,"MQ");
            MeasuredPeaks[k] = TempP->GetCentroid();
            delete TempP;
        }

        double_t oldSlope = pChannel->GetENGCoeff()[1];
        double_t oldOffset = pChannel->GetENGCoeff()[0];
        double_t ChargePeaks[2];

        for ( int k = 0; k < 2 ; k++ ) {
            ChargePeaks[k] = ( MeasuredPeaks[k] - oldOffset ) / oldSlope;
        }

        double_t newSlope = (gCalPeaks[1][0] - gCalPeaks[0][0])
            / (ChargePeaks[1] - ChargePeaks[0]);
        double_t newOffset = gCalPeaks[0][0] - newSlope*ChargePeaks[0];

        pChannel->DestroyENGCal();
        pChannel->AddENGCoefficient( static_cast<Float_t>( newOffset ) );
        pChannel->AddENGCoefficient( static_cast<Float_t>( newSlope ) );
    }

    TChannel::WriteToRoot();
    //TChannel::WriteCalFile("./newtest.cal");

    exit(EXIT_SUCCESS);
}
