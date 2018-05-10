// g++ LeanMatrices.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/libraries
// -lAnalysisTreeBuilder -lGriffin -lSceptar -lDescant -lPaces -lGRSIDetector
// -lTGRSIFit -lTigress -lSharc -lCSM -lTriFoil -lTGRSIint -lGRSILoop
// -lMidasFormat -lGRSIRootIO -lDataParser -lGRSIFormat -lMidasFormat
// -lXMLParser -lXMLIO -lProof -lGuiHtml `grsi-config --cflags --libs`
// `root-config --cflags --libs`  -lTreePlayer -lGROOT -lX11 -lXpm -lSpectrum

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sys/stat.h>
#include <utility>
#include <vector>

#include "Globals.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGRSIRunInfo.h"
#include "TGRSISortInfo.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "TPPG.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TVectorD.h"
#include "TVirtualIndex.h"
#include "TGRSIOptions.h"
#include "THnSparse.h"
#include "TSpline.h"
#include "TMVA/TSpline1.h"

#ifndef __CINT__
#include "TGriffin.h"
#include "TSceptar.h"
#include "TGRSISelector.h"
#endif

std::vector<TSpline*> ResidualVec;

// This function gets run if running interpretively
// Not recommended for the analysis scripts
#ifdef __CINT__
void LeanMatrices() {
    if (!AnalysisTree) {
        printf("No analysis tree found!\n");
        return;
    }
    // coinc window = 0-20, bg window 40-60, 6000 bins from 0. to 6000. (default
    // is 4000)
    TList *list = LeanMatrices(AnalysisTree, TPPG, TGRSIRunInfo, 0.);

    TFile *outfile = new TFile("output.root", "recreate");
    list->Write();
}
#endif

TList *LeanMatrices(TTree *tree, TPPG *ppg, TGRSIRunInfo *runInfo,
                    long maxEntries = 0, TStopwatch *w = NULL) {
    if (runInfo == NULL) {
        return NULL;
    }
    
    
    ///////////////////////////////////// SETUP
    //////////////////////////////////////////
    // gamma histogram limits
    Double_t low = 0;
    Double_t high = 6000;    // 10000
    Double_t nofBins = 6000; // 10000

    // Coincidence Parameters
    // times in ns
    Double_t ggTlow = 0.;
    Double_t ggThigh = 350.;

    Double_t ggBGlow = 1000.;
    Double_t ggBGhigh = 2000.;

    if (w == NULL) {
        w = new TStopwatch;
        w->Start();
    }

    TList *list = new TList;

    // TH1D* ggTimeDiff = new TH1D("ggTimeDiff", "#gamma-#gamma time
    // difference", 2000,0,4000); list->Add(ggTimeDiff);
    // TH1D* aaTimeDiff = new TH1D("aaTimeDiff", "#gamma-#gamma time
    // difference", 300,0,300); list->Add(aaTimeDiff);

    // Energy matrices for individual crystals and

    TH1D *Singles_total = new TH1D(
        "Singles_total", "#gamma singles for all crystals", nofBins, low, high);
    list->Add(Singles_total);
    TH2D *Singles_vs_Crystal =
        new TH2D("Singles_vs_Crystal", "#gamma singles for each crystal", 64, 0,
                 64, nofBins, low, high);
    list->Add(Singles_vs_Crystal);
    TH1D *Addback_total =
        new TH1D("Addback_total", "#gamma addback for all detectors", nofBins,
                 low, high);
    list->Add(Addback_total);
    TH2D *Addback_vs_Detector =
        new TH2D("Addback_vs_Detector", "#gamma addback for each detector", 16,
                 0, 16, nofBins, low, high);
    list->Add(Addback_vs_Detector);
    TH1D *Singles[64];
    TH1D *Addback[16];

    TH2D *ggsummat = new TH2D("ggsummat","#gamma-#gamma matrix 180 degrees",nofBins, low, high,nofBins, low, high); list->Add(ggsummat);


    for (int i = 0; i < 64; i++) {
        char name[128];

        sprintf(name, "Singles_%.2d", i);

        Singles[i] = new TH1D(name, Form("Singles for crystal %d", i), nofBins,
                              low, high);
        list->Add(Singles[i]);
    }
    for (int i = 0; i < 16; i++) {
        char name[128];

        sprintf(name, "Addback_%.2d", i);

        Addback[i] = new TH1D(name, Form("Addback for detector %d", i), nofBins,
                              low, high);
        list->Add(Addback[i]);
    }

    // Check what the times are
    TH1F *timeinrun = new TH1F(
        "timeinrun", "Time within run for #gamma singles", 1000, 0., 6.0e12);
    list->Add(timeinrun);

    TGriffin *grif = 0;
    tree->SetBranchAddress("TGriffin",
                           &grif); // We assume we always have a Griffin branch
            TGRSIOptions::AnalysisOptions()->SetCorrectCrossTalk(true);
        //TGRSIOptions::AnalysisOptions()->IsCorrectingCrossTalk();
    //grif->ResetFlags();
    if( ResidualVec.size() == 64 ) {
        printf("Loading in energy residuals\n");
        for (int k = 0; k < 64; k++) {
            grif->LoadEnergyResidual(k+1, ResidualVec[k]);
        }
    }

    // Indices of the two hits being compared
    int one;
    int two;

    std::cout << std::fixed
              << std::setprecision(
                     1); // This just make outputs not look terrible

    if (maxEntries == 0 || maxEntries > tree->GetEntries()) {
        maxEntries = tree->GetEntries();
    }

    int entry;
    // Only loop over the set number of entries
    // I'm starting at entry 1 because of the weird high stamp of 4
    for (entry = 1; entry < maxEntries; ++entry) {

        tree->GetEntry(entry);

        grif->ResetAddback();

        // loop over the gammas in the event packet
        for (one = 0; one < (int)grif->GetMultiplicity(); ++one) {
                
            timeinrun->Fill(grif->GetHit(one)->GetTime());
                if(grif->GetGriffinHit(one)->GetKValue() != 700) continue;
            int crystal = grif->GetGriffinHit(one)->GetCrystal() +
                          (grif->GetGriffinHit(one)->GetDetector() - 1) * 4;

            // We want to put every gamma ray in this event into the singles
            Singles_total->Fill(grif->GetGriffinHit(one)->GetEnergy());
            Singles_vs_Crystal->Fill(crystal,
                                     grif->GetGriffinHit(one)->GetEnergy());

            for (int i = 0; i < 64; i++) {

                if (crystal == i)
                    Singles[i]->Fill(grif->GetGriffinHit(one)->GetEnergy());

            } // crystal loop

                     // We now want to loop over any other gammas in this
               //packet
                     for(two = 0; two < (int) grif->GetMultiplicity(); ++two) {
                        if(two == one) continue; //If we are looking at the same
               //gamma we don't want to call it a coincidence

                        //ggTimeDiff->Fill(TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()));

                        if(ggTlow <=
               TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime())
               &&
               TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime())
               < ggThigh) {
                           //If they are close enough in time, fill the
               //gamma-gamma matrix. This will be symmetric because we are doing a
               //double loop over gammas
               
               if(grif->GetGriffinHit(one)->GetPosition()               //180 sum coincidence
                             .Angle(grif->GetGriffinHit(two)->GetPosition()) 
                           > 3.13)
               
                           ggsummat->Fill(grif->GetGriffinHit(one)->GetEnergy(),
               grif->GetGriffinHit(two)->GetEnergy());

                        } // gg prompt coincidences

                //        if(ggBGlow <=
              // TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime())
             //  &&
             //  TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime())
             //  < ggBGhigh) {
                           //If they are not close enough in time, fill the
               //time-random gamma-gamma matrix. This will be symmetric because we
               //are doing a double loop over gammas
                        } // gg time random coincidences

                     } // gg second gamma loop
            

         // g first gamma loop

        // loop over the addbacks in the event packet
        for (one = 0; one < (int)grif->GetAddbackMultiplicity(); ++one) {

            int detector = grif->GetAddbackHit(one)->GetDetector() - 1;
                if(grif->GetAddbackHit(one)->GetKValue() != 700) continue;
            // We want to put every gamma ray in this event into the singles
            Addback_total->Fill(grif->GetAddbackHit(one)->GetEnergy());
            Addback_vs_Detector->Fill(detector,
                                      grif->GetAddbackHit(one)->GetEnergy());

            for (int i = 0; i < 16; i++) {

                if (detector == i)
                    Addback[i]->Fill(grif->GetAddbackHit(one)->GetEnergy());

            } // detectorloop

            // We now want to loop over any other gammas in this packet
            /*         for(two = 0; two < (int) grif->GetAddbackMultiplicity();
               ++two) {
                                                    if(two == one) continue;
               //If we are looking at the same gamma we don't want to call it a
               coincidence


                        //VINZENZ I THINK THIS IS BREAKING THIS FOR SOME REASON.
                        aaTimeDiff->Fill(TMath::Abs(grif->GetAddbackHit(two)->GetTime()-grif->GetAddbackHit(one)->GetTime()));

                     } // aa multiplicity loop
            */
        } // a loop

        if ((entry % 10000) == 0) {
            printf("Completed %d of %ld \r", entry, maxEntries);
        }

    } // entry loop

    std::cout << "creating histograms done after " << w->RealTime()
              << " seconds" << std::endl;
    w->Continue();
    return list;
}

// This function gets run if running in compiled mode
#ifndef __CINT__
int main(int argc, char **argv) {
    if (argc != 4 && argc != 3 && argc != 2) {
        printf("try again (usage: %s <analysis tree file> <optional: output "
               "file> <max entries>).\n",
               argv[0]);
        return 0;
    }

    // We use a stopwatch so that we can watch progress
    TStopwatch w;
    w.Start();

    TFile *file = new TFile(argv[1]);
    //TGRSIOptions::Get()->ReadFromFile(file);
    //TGRSIOptions::AnalysisOptions()->SetCorrectCrossTalk(true);
    if (file == NULL) {
        printf("Failed to open file '%s'!\n", argv[1]);
        return 1;
    }
    if (!file->IsOpen()) {
        printf("Failed to open file '%s'!\n", argv[1]);
        return 1;
    }
    printf("Sorting file:" DBLUE " %s" RESET_COLOR "\n", file->GetName());

    // Get PPG from File
    TPPG *myPPG = (TPPG *)file->Get("TPPG");
    /*      if(myPPG == NULL) {
            printf("Failed to find PPG information in file '%s'!\n",argv[1]);
            return 1;
           }*/

    TFile *pResFile = nullptr;
    if ( argc > 2 ) // Check if the extra file is tacked on
        pResFile = new TFile( argv[2], "READ");
    if ( pResFile != nullptr )
    {
        pResFile->cd();
        if (pResFile->cd("Energy_Residuals")) {
            printf("Energy residuals found, loading...\n");
            TGraph* TempGraph;
            for (int k = 0 ; k < 64; k++) {
                gDirectory->GetObject(Form("Graph;%d", k + 1), TempGraph);
                ResidualVec.push_back( new TMVA::TSpline1("", TempGraph );
            }
        } else {
            printf("No energy residuals found\n");
        }
        pResFile->Close();
    }

    // Get run info from File
    TGRSIRunInfo *runInfo = (TGRSIRunInfo *)file->Get("TGRSIRunInfo");
    TGRSIRunInfo::Get()->SetRunInfo(runInfo);
    if (runInfo == NULL) {
        printf("Failed to find run information in file '%s'!\n", argv[1]);
        return 1;
    }

    TTree *tree = (TTree *)file->Get("AnalysisTree");
    TChannel::ReadCalFromTree(tree);
    if (tree == NULL) {
        printf("Failed to find analysis tree in file '%s'!\n", argv[1]);
        return 1;
    }

    // Get the TGRSIRunInfo from the analysis Tree.


    TList *list; // We return a list because we fill a bunch of TH1's and shove
                 // them into this list.
                 
    TFile *outfile;
    if (runInfo == nullptr) {
        printf(
            "Could not find run info, please provide output file name\n");
        return 0;
    }
    int runnumber = runInfo->RunNumber();
    int subrunnumber = runInfo->SubRunNumber();
    outfile = new TFile(
        Form("matrix%05d_%03d.root", runnumber, subrunnumber), "recreate");             
                 
/*    TFile *outfile;
    if (argc < 3) {
        if (!runInfo) {
            printf(
                "Could not find run info, please provide output file name\n");
            return 0;
        }
        int runnumber = runInfo->RunNumber();
        int subrunnumber = runInfo->SubRunNumber();
        outfile = new TFile(
            Form("calmat%05d_%03d.root", runnumber, subrunnumber), "recreate");
    } else {
        outfile = new TFile(argv[2], "recreate");
    }
*/
    std::cout << argv[0] << ": starting Analysis after " << w.RealTime()
              << " seconds" << std::endl;
    w.Continue();
    if (argc < 4) {
        list = LeanMatrices(tree, myPPG, runInfo, 0, &w);
    } else {
        int entries = atoi(argv[3]);
        std::cout << "Limiting processing of analysis tree to " << entries
                  << " entries!" << std::endl;
        list = LeanMatrices(tree, myPPG, runInfo, entries, &w);
    }
    if (list == NULL) {
        std::cout << "LeanMatrices returned TList* NULL!\n" << std::endl;
        return 1;
    }

    printf("Writing to File: " DYELLOW "%s" RESET_COLOR "\n",
           outfile->GetName());
    list->Write();

    outfile->Close();

    std::cout << argv[0] << " done after " << w.RealTime() << " seconds"
              << std::endl
              << std::endl;

    return 0;
}

#endif
