/*
 * Magical Nonlinear finder and TSpline Generator
 *
 * Created by Kurtis Raymond (kraymond@sfu.ca)
 *
 * The analysis tree loaded in should have its channels already
 * gain matched.
 *
 * This script takes in an analsys tree from the provided root file,
 * measures the amount the peaks specified (in gPeaks) from the measured data.
 *
 * The resulting TGraphs are stored in a directory "Energy_Residuals" in the
 * provided root file.
 *
 */

/* 
 * Example loading script. Put at the begining of your other analysis scripts
 * such as LeanMatricies.cxx, after the TGriffin address has been loaded.
 *
 * Change pGriff to the TGriffin address pointer variable.
 *
 * Note that the TGraphs are stored from 1->64, but they are loaded in through
 * LoadEnergyResidual from 0->63.

    if( gFile->cd("Energy_Residuals") ) {
        printf("Energy residual calibration data found. Loading...\n");

        TGraph* TempResidual;
        for( int i = 1 ; i <= 64 ; i++ ) {
            gDirectory->GetObject(Form("Graph;%d",i), TempResidual);
            pGriff->LoadEnergyResidual( i-1, TempResidual);
        }
        gFile->cd(); // Return to the top directory
        printf("Done.\n");
    }
 */

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

#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF1.h"
#include "TFile.h"
#include "TGainMatch.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TPPG.h"
#include "TPad.h"
#include "TPeak.h"
#include "TROOT.h"
#include "TScaler.h"
#include "TSpectrum.h"
#include "TSpline.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

#ifndef __CINT__
#include "TGriffin.h"
#include "TSceptar.h"
#endif

//=== Global Variables ===//
// 1 for fragment tree
// 0 for analysis tree
const bool gIsFragmentFile = 0;

// Input calibration peaks for finding nonlinearities
/*
const std::vector<double_t> gPeaks = {121.7817, 244.6974, 964.057,
                                      1085.837, 1112.076};

const std::vector<double_t> gWidths = {16, 20, 20, 13, 13};
*/

// Key peaks for 129Sn
const std::vector<double_t> gPeaks = {
    315.42,
    511.0,
    570.41,
    645.2,
    728.53,
    769.31,
//    907.34, // Double Peak
    1008.53,
    1054.3,
//    1222.51,
//    1781.54,
    1864.89,
    2118.26,
    2546.61};

const std::vector<double_t> gWidths = {
    20,
    20,
    20,
    20,
    20,
    20,
//    20,
    20,
    15,
//    20,
//    20,
    20,
    20,
    20};

// Used for displaying individual peak fitting information
const bool gPrintFlag = true;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s <fragment or analysis tree file>.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    TFile *pFile = new TFile(argv[1], "update");

    if (!pFile->IsOpen()) {
        printf("Failed to open file '%s'\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    if( gPeaks.size() != gWidths.size() ) {
        printf("Number of peaks and widths are different!\n");
        exit(EXIT_FAILURE);
    }

    TTree *pTree = nullptr;

    if (gIsFragmentFile)
        pTree = (TTree *)pFile->Get("FragmentTree");
    else
        pTree = (TTree *)pFile->Get("AnalysisTree");

    TChannel::ReadCalFromTree(pTree);

    if (pTree == nullptr) {
        printf("Failed to find fragment or analysis tree in file '%s'.\n",
               argv[1]);
        exit(EXIT_FAILURE);
    }

    // Setup TGriffin
    TGriffin *pGriff = nullptr;
    pTree->SetBranchAddress("TGriffin", &pGriff);

    //TChannel *pChannel = nullptr;
    //TChannel::ReadCalFromTree(pTree);

    printf("Generating empty matrix");
    // Load energy matrix
    TH2D *mat_en = new TH2D("mat_en", "", 64, 0, 64, 5000, 0, 5000);

    printf("Filling energy matrix");
    // Load in Energy data
    if (gIsFragmentFile)
        pTree->Project("mat_en", "TFragment.GetEnergy():"
                                 "TFragment.GetChannelNumber()");
    else
        pTree->Project("mat_en",
                       "TGriffin.fGriffinLowGainHits.GetEnergy():"
                       "TGriffin.fGriffinLowGainHits.GetChannel().fNumber");

    // Make a list that will store all the energy differences
    TList* LNonlinearitiesGraphs = new TList();
    TList* LNonlinearitiesGraphsErr = new TList();

    int nPeaks = gPeaks.size();

    for (int i = 0; i < 64; i++) {
        printf("Starting new channel %d:\n", i);

        // Project the project the histogram for the current channel
        TH1D *h_en = mat_en->ProjectionY(Form("h_%.2i", i), i + 1, i + 1);

        std::vector<double_t> EngDiff = {};
        std::vector<double_t> EngDiffErr = {};
        std::vector<double_t> EngX = {};

        // Fit all the peaks in our calibration and collect their centroids
        for (int k = 0; k < nPeaks; k++) {
            double_t CalPeak, DataPeak, CalWidth;
            double_t DataPeakErr;
            CalPeak = gPeaks[k];
            CalWidth = gWidths[k];

            // We use TSpectrum::Search() to grab all the peaks. Output is
            // ordered from the most intense peak to the least.
            if (gPrintFlag)
                printf("Fitting peak %g .", gPeaks[k]);
            TSpectrum s;
            h_en->GetXaxis()->SetRangeUser(CalPeak - CalWidth,
                                           CalPeak + CalWidth);
            s.Search(h_en, 2, "", 0.15);    // Hist, Sigma, Opt, Threshold
            DataPeak = s.GetPositionX()[0]; // Grab most intense peak

            if ( DataPeak < 1 ) { // Errors commonly fall under this condition
                printf("Could not find Peak, Skipping.\n");
                continue;
            }

            if (gPrintFlag)
                printf("Roughly at %g ", DataPeak);
            h_en->GetXaxis()->UnZoom();

            // Fit the peak
            TPeak *CurPeak =
                new TPeak(DataPeak, DataPeak - CalWidth, DataPeak + CalWidth);
            CurPeak->Fit(h_en, "MQ+"); // Quiet Flag
            DataPeak = CurPeak->GetCentroid();
            DataPeakErr = CurPeak->GetCentroidErr();

            // Report the peak
            if (gPrintFlag)
                printf("... found at %g, ", DataPeak);
            // We compute the quantanty that will be subtracted by the data
            EngDiff.push_back( DataPeak - CalPeak );
            EngDiffErr.push_back(DataPeakErr);
            if (gPrintFlag)
                printf(" difference of %g\n", EngDiff.back());
            EngX.push_back(DataPeak);

            // Loop Cleanup
            delete CurPeak;
        }

        // Boundary conditions to prevent too much exterpolation
        EngX.push_back(EngX.back() + 10);
        EngDiffErr.push_back(0.0);
        EngDiff.push_back(0.0);
        EngX.push_back(EngX.back() + 20);
        EngDiffErr.push_back(0.0);
        EngDiff.push_back(0.0);

        // Make a TGraph that can be used for interpolating the values
        TGraph* TempGraph = new TGraph(EngX.size(), EngX.data(), EngDiff.data());
        TempGraph->SetTitle("");
        LNonlinearitiesGraphs->Add(TempGraph);

        // Make a graph that represents this channel's offsets and errors
        TGraphErrors* TempGraphErr = new TGraphErrors(EngX.size(), EngX.data(),
                EngDiff.data(), EngDiffErr.data(), EngDiffErr.data() );
        LNonlinearitiesGraphsErr->Add(TempGraphErr);
    }

    // We want to make an extra directory to store all of our energy
    // residuals in. The assumption is made that these TGraphs are
    // written in order.
    //
    // See sample_load_code.cxx for details on how to load these.

    printf("Writing Non-Linarities\n");
    TDirectory* NonLinearDirectory;
    TIter next = TIter(LNonlinearitiesGraphs->MakeIterator());
    if ( pFile->cd("Energy_Residuals") ) {
        NonLinearDirectory = gDirectory; // Set to current directory
        NonLinearDirectory->Delete("Graph;*");
    }
    else
        NonLinearDirectory = pFile->mkdir("Energy_Residuals");

    NonLinearDirectory->cd();
    LNonlinearitiesGraphs->Write();
    pFile->cd();

    printf("Writing Graphs\n");
    // Compose a graph that gives an overview of all the residuals.
    TCanvas * c1 = new TCanvas("Residuals", "Residuals", 800, 800);
    c1->SetFrameBorderMode(0);
    c1->Divide(4,4);
    next = TIter(LNonlinearitiesGraphsErr->MakeIterator());
    for( int i = 1 ; i <= 16 ; i++ ) {
        // cd(0) is the canvas itself
        // Iterate through all the TPads
        c1->cd(i);
        gPad->SetTitle(Form("Digitizer %d", i));
        gPad->SetLeftMargin(0.05);
        gPad->SetBottomMargin(0.05);
        gPad->SetRightMargin(0.00);
        gPad->SetTopMargin(0.00);
        TMultiGraph* mg = new TMultiGraph();
        for(int k = 0; k < 4; k++ ) {
            // Iterate through the list to produce four data sets in
            // the multigraph 
            TGraphErrors* obj = (TGraphErrors*) next();
            obj->SetMarkerStyle(31); // kStar
            obj->SetLineColor(k+1);
            obj->SetMarkerColor(k+1); // Nice marker colors
            mg->Add(obj);
        }
        mg->Draw("APL");
    }
    c1->Update();
    c1->Draw();
    c1->Write("Non-Linearity Overview");

    // Project Matrix
    // Cleanup
    delete mat_en;
    pFile->Close();
}
