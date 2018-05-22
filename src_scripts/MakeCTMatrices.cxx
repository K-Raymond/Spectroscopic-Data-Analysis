// g++ LeanMatrices.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/libraries
// -lAnalysisTreeBuilder -lGriffin -lSceptar -lDescant -lPaces -lGRSIDetector
// -lTGRSIFit -lTigress -lSharc -lCSM -lTriFoil -lTGRSIint -lGRSILoop
// -lMidasFormat -lGRSIRootIO -lDataParser -lGRSIFormat -lMidasFormat
// -lXMLParser -lXMLIO -lProof -lGuiHtml `grsi-config --cflags --libs`
// `root-config --cflags --libs`  -lTreePlayer -lGROOT -lX11 -lXpm -lSpectrum

//

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

std::vector<TMVA::TSpline1*> ResidualVec;

// This function gets run if running interpretively
// Not recommended for the analysis scripts
#ifdef __CINT__
void LeanMatrices() {
    if (!AnalysisTree) {
        printf("No analysis tree found!\n");
        return;
    }
    // coinc window = 0-20, bg window 40-60, 6000 bins from 0. to 6000. (default is 4000)
    TList *list = LeanMatrices(AnalysisTree, 0.);

    TFile *outfile = new TFile("output.root", "recreate");
    list->Write();
}
#endif

TList *LeanMatrices(TTree *tree, long maxEntries = 0, TStopwatch *w = NULL)
{

  ///////////////////////////////////// SETUP
  //////////////////////////////////////////
  // gamma histogram limits
  Double_t low = 0;
  Double_t high = 2000;    // 10000
  Double_t nofBins = 2000; // 10000

  if (w == NULL) {
      w = new TStopwatch;
      w->Start();
  }

  TList *list = new TList;

  TH2F *histos[96];

	for(int det_num=0; det_num<16; det_num++){
		histos[det_num*6]   = new TH2F(Form("det_%d_0_1",det_num+1),"",nofBins,low,high,nofBins,low,high);
		histos[det_num*6+1] = new TH2F(Form("det_%d_0_2",det_num+1),"",nofBins,low,high,nofBins,low,high);
		histos[det_num*6+2] = new TH2F(Form("det_%d_0_3",det_num+1),"",nofBins,low,high,nofBins,low,high);
		histos[det_num*6+3] = new TH2F(Form("det_%d_1_2",det_num+1),"",nofBins,low,high,nofBins,low,high);
		histos[det_num*6+4] = new TH2F(Form("det_%d_1_3",det_num+1),"",nofBins,low,high,nofBins,low,high);
		histos[det_num*6+5] = new TH2F(Form("det_%d_2_3",det_num+1),"",nofBins,low,high,nofBins,low,high);
	}

	for( int i=0; i<96; i++) { list->Add(histos[i]);}


  TGriffin *grif = 0;
  tree->SetBranchAddress("TGriffin", &grif); // We assume we always have a Griffin branch
	// try to make sure "cross talk" isn't being used (even though it shouldn't exist now anyway)
  TGRSIOptions::AnalysisOptions()->SetCorrectCrossTalk(false);

  grif->ResetFlags();

  if( ResidualVec.size() == 64 ) {
      printf("Loading in energy residuals\n");
      for (int k = 0; k < 64; k++) {
          grif->LoadEnergyResidual(k+1, ResidualVec[k]);
      }
  }

  // Indices of the two hits being compared
  int one;
  int two;

  std::cout << std::fixed << std::setprecision(1); // This just make outputs not look terrible

  if (maxEntries == 0 || maxEntries > tree->GetEntries()) { maxEntries = tree->GetEntries(); }



/*	// ****** Example from CrossTalk.C selector script
// steps: 1. distribute total multiplicity among each detector (1-16)
// 2. first loop through hits in event (through entire multiplicity, across different detectors)
// 3. reject pileup of current first hit index
// 4. second loop through subsequent hits in event (only higher values of multiplicity [not symmetrized], across different detectors)
// 5. reject pileup of current second hit index
// 6. check if (1) current first indexed detector has mult = 2, (2) if both indices have same detector, and (3) prompt time difference (300)
// 7. order the indices by crystal number (0-2 for "low" and 1-3 for "high")
// 8. exclude same-crystal hits
// 9. fill histogram with crystal ordering
void CrossTalk::FillHistograms() {
	//find the multiplicity in each clover over the entire event
   //we do this because we want to force a multiplicity of 2
	Int_t det_multiplicity[17] = {0};
	for(auto gr1 = 0; gr1 < fGrif->GetMultiplicity(); ++gr1){
		++(det_multiplicity[fGrif->GetGriffinHit(gr1)->GetDetector()]);
	}
   for(auto gr1 = 0; gr1 < fGrif->GetMultiplicity(); ++gr1){
		if(pileup_reject && (fGrif->GetGriffinHit(gr1)->GetKValue() != 700)) continue; //This pileup number might have to change for other expmnts

		//fH1[Form("gEdet%d",fGrif->GetGriffinHit(gr1)->GetDetector())]->Fill(fGrif->GetGriffinHit(gr1)->GetEnergy());
		//fH2["gE_chan"]->Fill(fGrif->GetGriffinHit(gr1)->GetArrayNumber(),fGrif->GetGriffinHit(gr1)->GetEnergy());
		//fH1["gE"]->Fill(fGrif->GetGriffinHit(gr1)->GetEnergy());
		//fH1["gEnoCT"]->Fill(fGrif->GetGriffinHit(gr1)->GetNoCTEnergy());

		for(auto gr2 = gr1 + 1; gr2 < fGrif->GetMultiplicity(); ++gr2){
			if(pileup_reject && fGrif->GetGriffinHit(gr2)->GetKValue() != 700) continue; //This pileup number might have to change for other expmnts

			if((det_multiplicity[fGrif->GetGriffinHit(gr1)->GetDetector()] == 2) && Addback(*(fGrif->GetGriffinHit(gr1)), *(fGrif->GetGriffinHit(gr2)))){

				TGriffinHit *low_crys_hit, *high_crys_hit;

				if(fGrif->GetGriffinHit(gr1)->GetCrystal() < fGrif->GetGriffinHit(gr2)->GetCrystal()){
					low_crys_hit = fGrif->GetGriffinHit(gr1);
					high_crys_hit = fGrif->GetGriffinHit(gr2);
				}
				else{
					low_crys_hit = fGrif->GetGriffinHit(gr2);
					high_crys_hit = fGrif->GetGriffinHit(gr1);
				}
				if(low_crys_hit->GetCrystal() != high_crys_hit->GetCrystal()){
					//fH2[Form("det_%d_%d_%d",low_crys_hit->GetDetector(),low_crys_hit->GetCrystal(),high_crys_hit->GetCrystal())]->Fill(low_crys_hit->GetNoCTEnergy(),high_crys_hit->GetNoCTEnergy());
					fH2[Form("det_%d_%d_%d",low_crys_hit->GetDetector(),low_crys_hit->GetCrystal(),high_crys_hit->GetCrystal())]->Fill(low_crys_hit->GetEnergy(),high_crys_hit->GetEnergy());
				}
			}
		}
	}

//   for(auto gr1 = 0; gr1 < fGrif->GetAddbackMultiplicity(); ++gr1) {
//      if(pileup_reject && (fGrif->GetAddbackHit(gr1)->GetKValue() != 700))
//         continue; // This pileup number might have to change for other expmnts
//      fH1["aE"]->Fill(fGrif->GetAddbackHit(gr1)->GetEnergy());
//      fH1[Form("aEdet%d", fGrif->GetAddbackHit(gr1)->GetDetector())]->Fill(fGrif->GetAddbackHit(gr1)->GetEnergy());
//      fH1["aMult"]->Fill(fGrif->GetNAddbackFrags(gr1));
//      if(fGrif->GetNAddbackFrags(gr1) == 2)
//         fH1[Form("aE2det%d", fGrif->GetAddbackHit(gr1)->GetDetector())]->Fill(fGrif->GetAddbackHit(gr1)->GetEnergy());
//   }
}

*/	// ********


  int entry;
  // Only loop over the set number of entries
  // I'm starting at entry 1 because of the weird high stamp of 4
  for (entry = 1; entry < maxEntries; ++entry)
	{

		tree->GetEntry(entry);

		int det_mult[17] = {0};

		// 1. distribute total multiplicity among each detector (1-16)
		for( one = 0; one < grif->GetMultiplicity(); one++)
		{
			++(det_mult[grif->GetGriffinHit(one)->GetDetector()]);
		}

		// 2. first loop through hits in event (through entire multiplicity, across different detectors)
		for (one = 0; one < grif->GetMultiplicity(); one++)
		{

			// 3. reject pileup of current first hit index
			if(grif->GetGriffinHit(one)->GetKValue() != 700) continue;

			// 4. second loop through subsequent hits in event (only higher values of multiplicity [not symmetrized], across different detectors)
			for( two = one+1; two < grif->GetMultiplicity(); two++)
			{
				// 5. reject pileup of current second hit index
				if( grif->GetGriffinHit(two)->GetKValue() != 700 ) continue;

				// 6. check if (1) current first indexed detector has mult = 2, (2) if both indices have same detector, and (3) prompt time difference (300)
				// 8. exclude same-crystal hits
				if( det_mult[grif->GetGriffinHit(one)->GetDetector()] !=2 ) continue;
				if( grif->GetGriffinHit(one)->GetDetector() != grif->GetGriffinHit(two)->GetDetector() ) continue;
				if( std::fabs( grif->GetGriffinHit(one)->GetTime() - grif->GetGriffinHit(two)->GetTime() ) > 300. ) continue;
				if( grif->GetGriffinHit(one)->GetCrystal() == grif->GetGriffinHit(two)->GetCrystal() ) continue;

				// 7. order the indices by crystal number (0-2 for "low" and 1-3 for "high")
				TGriffinHit* low_crys_hit;
				TGriffinHit* high_crys_hit;

				if( grif->GetGriffinHit(one)->GetCrystal() < grif->GetGriffinHit(two)->GetCrystal() )
				{
					low_crys_hit = grif->GetGriffinHit(one);
					high_crys_hit = grif->GetGriffinHit(two);
				}
				else
				{
					low_crys_hit = grif->GetGriffinHit(two);
					high_crys_hit = grif->GetGriffinHit(one);
				}

				// 9. fill histogram with crystal ordering
				int hist_index = -1;

				if( low_crys_hit->GetCrystal() == 0 && high_crys_hit->GetCrystal() == 1) { hist_index = 0; }
				if( low_crys_hit->GetCrystal() == 0 && high_crys_hit->GetCrystal() == 2) { hist_index = 1; }
				if( low_crys_hit->GetCrystal() == 0 && high_crys_hit->GetCrystal() == 3) { hist_index = 2; }
				if( low_crys_hit->GetCrystal() == 1 && high_crys_hit->GetCrystal() == 2) { hist_index = 3; }
				if( low_crys_hit->GetCrystal() == 1 && high_crys_hit->GetCrystal() == 3) { hist_index = 4; }
				if( low_crys_hit->GetCrystal() == 2 && high_crys_hit->GetCrystal() == 3) { hist_index = 5; }

				hist_index += ( 6*( grif->GetGriffinHit(one)->GetDetector()-1 ) );

				// Fill( low crystal, high crystal );
				histos[hist_index]->Fill(low_crys_hit->GetEnergy(),high_crys_hit->GetEnergy());

			}	// second gamma loop

		}	// first gamma loop

    if ((entry % 10000) == 0) { printf("Completed %d of %ld \r", entry, maxEntries); }

	} // entry loop

  std::cout << "creating histograms done after " << w->RealTime() << " seconds" << std::endl;
  w->Continue();
  return list;
}



// This function gets run if running in compiled mode
#ifndef __CINT__
int main(int argc, char **argv) {
    if ( argc != 3 && argc != 2) {
        printf("try again (usage: %s <analysis tree file> <optional:residuals file>).\n", argv[0]);
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
                ResidualVec.push_back( new TMVA::TSpline1("", TempGraph ) );
            }
        } else {
            printf("No energy residuals found\n");
        }
        pResFile->Close();
    }


    TTree *tree = (TTree *)file->Get("AnalysisTree");
    TChannel::ReadCalFromTree(tree);
    if (tree == NULL) {
        printf("Failed to find analysis tree in file '%s'!\n", argv[1]);
        return 1;
    }

    TList *list; // We return a list because we fill a bunch of TH1's and shove them into this list.

    TFile *outfile = new TFile("CrossTalk_histos.root","recreate");

    std::cout << argv[0] << ": starting Analysis after " << w.RealTime() << " seconds" << std::endl;
    w.Continue();

    if (argc < 4) {
        list = LeanMatrices(tree, 0, &w);
    } else {
        int entries = atoi(argv[3]);
        std::cout << "Limiting processing of analysis tree to " << entries << " entries!" << std::endl;
        list = LeanMatrices(tree, entries, &w);
    }
    if (list == NULL) {
        std::cout << "LeanMatrices returned TList* NULL!\n" << std::endl;
        return 1;
    }

    printf("Writing to File: " DYELLOW "%s" RESET_COLOR "\n",outfile->GetName());
    list->Write();

    outfile->Close();

    std::cout << argv[0] << " done after " << w.RealTime() << " seconds" << std::endl << std::endl;

    return 0;
}

#endif
