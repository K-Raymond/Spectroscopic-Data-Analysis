//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 25 13:18:27 2016 by ROOT version 5.34/24
// from TTree FragmentTree/FragmentTree
// found on file: fragment07844_000.root
//////////////////////////////////////////////////////////

#ifndef CrossTalk_h
#define CrossTalk_h

#include "TChain.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"

// Header file for the classes stored in the TTree if any.
#include "TGriffin.h"
#include "TSceptar.h"
#include "TGRSISelector.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class CrossTalk : public TGRSISelector {

public:
   TGriffin* fGrif;
   TSceptar* fScep;
   std::vector<TGraph*> fResidualVec;

   CrossTalk(TTree* /*tree*/ = 0) : TGRSISelector(), fGrif(0), fScep(0) { SetOutputPrefix("Crosstalk"); }
   virtual ~CrossTalk() {}
   virtual Int_t Version() const { return 2; }
   void          CreateHistograms();
   void          FillHistograms();
   void InitializeBranches(TTree* tree);

   ClassDef(CrossTalk, 2);
};

#endif

#ifdef CrossTalk_cxx
void CrossTalk::InitializeBranches(TTree* tree)
{
   if(!tree) return;
   tree->SetBranchAddress("TGriffin", &fGrif);
   //tree->SetBranchAddress("TSceptar", &fScep);
    TFile *pResFile = nullptr;
    if ( argc > 2 )
        pResFile = new TFile( "residuals.root", "READ");
    if ( pResFile != nullptr ) {
        pResFile->cd();
        if ( pResFile->cd("Energy_Residuals") ) {
            printf("Energy residuals found, loading...\n");
            TGraph* TempGraph;
            for (int k = 0 ; k < 64; k++) {
                gDirectory->GetObject(Form("Graph;%d", k + 1), TempGraph);
                ResidualVec.push_back( TempGraph );
            }
        } else {
            printf("No energy residuals found\n");
        }
        pResFile->Close();
    }


}

#endif // #ifdef CrossTalk_cxx
