#include "AliESDEvent.h"
#include "AliVParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliVEventHandler.h"

#include "AliAnalysisTaskCheckMCCountsESD.h"

ClassImp(AliAnalysisTaskTrackingEffPID);

AliAnalysisTaskCheckMCCountsESD::AliAnalysisTaskCheckMCCountsESD() : AliAnalysisTaskSE("CheckMCCountsESD"),
                                                                     fOutputList{0x0},
                                                                     fHistCounts{0x0},
                                                                     fInd(0)
{
  fDebug.open("debug_esd.txt");
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskCheckMCCountsESD::~AliAnalysisTaskCheckMCCountsESD()
{
  if (fOutputList)
    delete fOutputList;
}

void AliAnalysisTaskCheckMCCountsESD::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  fHistCounts = new TH1F("PDG counts ESD", "PDG counts ESD", 4, 0, 4);
  fHistCounts->GetXaxis()->SetBinLabel(1, "el");
  fHistCounts->GetXaxis()->SetBinLabel(2, "ka");
  fHistCounts->GetXaxis()->SetBinLabel(3, "pi");
  fHistCounts->GetXaxis()->SetBinLabel(4, "pr");
  fOutputList->Add(fHistCounts);

  PostData(1, fOutputList);
}

int AliAnalysisTaskCheckMCCountsESD::GetBinNumber(Int_t pdg)
{
  switch (pdg) {
    case kElectron:
      return 0;
    case kKPlus:
      return 1;
    case kPiPlus:
      return 2;
    case kProton:
      return 3;
    default:
      return -1;
  }
}

void AliAnalysisTaskCheckMCCountsESD::UserExec(Option_t*)
{
  AliVEvent* ev = fInputEvent;
  if (!ev) {
    AliFatal("NO EVENT FOUND!");
    return;
  }
  if (!fMCEvent) {
    AliFatal("NO MC INFO FOUND");
    return;
  }

  for (int iMC = 0; iMC < fMCEvent->GetNumberOfTracks(); ++iMC) {
    AliVParticle* part = (AliVParticle*)fMCEvent->GetTrack(iMC);
    const int pdg = part->PdgCode();
    bool isPhysicalPrimary = fMCEvent->IsPhysicalPrimary(iMC);
    fDebug << "Particle " << fInd << " PDG: " << pdg << " physical primary: " << isPhysicalPrimary << std::endl;
    if (isPhysicalPrimary) {
      int bin = GetBinNumber(pdg);
      if (bin != -1) {
        fDebug << "Particle " << fInd << " bin: " << bin << " PDG code: " << pdg << " kPion: " << kPiPlus << " kProton: " << kProton << " kElectron: " << kElectron << " kKPlus: " << kKPlus << std::endl;
        fHistCounts->Fill(bin);
      }
    }
    fInd++;
  }

  PostData(1, fOutputList);
}

void AliAnalysisTaskCheckMCCountsESD::Terminate(Option_t*)
{
  fDebug.close();

  TCanvas* canHis = new TCanvas("Counts", "Counts", 800, 600);
  canHis->cd();
  fHistCounts->Draw();

  TLatex* latexa = new TLatex();
  latexa->SetTextSize(0.04);
  latexa->SetTextFont(42);
  latexa->SetTextAlign(22);
  for (int i = 1; i <= fHistCounts->GetNbinsX(); i++) {
    std::string binstr = std::to_string((int)(fHistCounts->GetBinContent(i)));
    float bincen = 0.1f + fHistCounts->GetXaxis()->GetBinCenter(i) / (fHistCounts->GetNbinsX() + 1);
    latexa->DrawLatexNDC(bincen, 0.2f, binstr.c_str());
  }

  canHis->SaveAs("mccounts_esd.png");
  delete canHis;
  return;
}
