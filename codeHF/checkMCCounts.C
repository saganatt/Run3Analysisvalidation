#include <TH1.h>
#include <TFile.h>
#include <iostream>

int getBinNumber(Int_t pdg)
{
  switch (pdg) {
    case 0:
      return 1;
    case 1:
      return 2;
    case 2:
      return 3;
    case 3:
      return 4;
    default:
      return -1;
  }
}

int checkMCCounts(TString filename = "O2_AOD.root", bool isAOD)
{
  TFile* f = new TFile(filename.Data());
  if (f->IsZombie()) {
    printf("Failed to open file %s\n", filename.Data());
    return 1;
  }

  TH1F* histCounts = new TH1F("PDG counts", "PDG counts", 4, 0, 4);
  histCounts->GetXaxis()->SetBinLabel(1, "el");
  histCounts->GetXaxis()->SetBinLabel(2, "ka");
  histCounts->GetXaxis()->SetBinLabel(3, "pi");
  histCounts->GetXaxis()->SetBinLabel(4, "pr");

  if (isAOD) {
  } else {
    TTree* tracks = (TTree*)f->get("esdTree/Tracks");
    Int_t pdg;
    Bool_t isPhysicalPrimary;
    tracks->SetBranchAddress("PdgCode", &pdg);
    tracks->SetBranchAddress("IsPhysicalPrimary", &isPhysicalPrimary);
    for (int i = 0; i < tracks->GetEntries(); i++) {
      tracks->GetEntry(i);
      std::cout << "PDG: " << pdg << " physical primary: " << isPhysicalPrimary << std::endl;
      if (isPhysicalPrimary) {
        int bin = getBinNumber(pdg);
        if (bin != -1) {
          histCounts->Fill(bin);
        }
      }
    }
  }
}
