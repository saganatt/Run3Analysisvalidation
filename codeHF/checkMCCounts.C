#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TPDGCode.h>
#include <iostream>

#include "AliRunLoader.h"
#include "AliStack.h"

int getBinNumber(Int_t pdg)
{
  switch (pdg) {
    case kElectron:
      return 1;
    case kKPlus:
      return 2;
    case kPiPlus:
      return 3;
    case kProton:
      return 4;
    default:
      return -1;
  }
}

int checkCountsAOD(TString filename = "AO2D.root", TH1F* histCounts = nullptr)
{
  TFile* f = new TFile(filename.Data());
  if (f->IsZombie()) {
    printf("Failed to open file %s\n", filename.Data());
    return 1;
  }
  TTree* tracks = (TTree*)f->Get("O2tracks");
  std::cout << "Not implemented yet!" << std::endl;
  return 1;
}

int checkCountsESD(TString filepath = "/data2/data/Run2/pp_13TeV/sim/LHC20f4a/294774/001/", TH1F* histCounts = nullptr)
{
  AliRunLoader* RunLoader = AliRunLoader::Open(TString::Format("%s/%s", filepath.Data(), "galice.root"), "read");
  if (RunLoader == 0x0) {
    printf("Error opening %s file \n", filepath.Data());
    return 1;
  }
  RunLoader->LoadHeader();
  RunLoader->LoadKinematics();
  int nevents = RunLoader->GetNumberOfEvents();

  for (int ievent = 0; ievent < nevents; ievent++) {
    RunLoader->GetEvent(ievent);
    AliStack* stack = RunLoader->Stack();
    Int_t nparticles = stack->GetNtrack();
    for (int i = 0; i < nparticles; i++) {
      TParticle* p = stack->Particle(i);
      //std::cout << "Particle PDG: " << p->GetPdgCode() << " physical primary: " << p->IsPrimary() << std::endl;
      if (p->IsPrimary()) {
        int bin = getBinNumber(p->GetPdgCode());
        if (bin != -1) {
          histCounts->Fill(bin);
        } else {
          std::cout << "Unidentified PDG: " << p->GetPdgCode() << std::endl;
        }
      }
    }
  }
  return 0;
}

//int checkCountsESDManager(TString txtfile = "./list_ali.txt")
//{
//  AliAnalysisManager* mgr = new AliAnalysisManager("testAnalysis");
//  TChain* chainESD = CreateLocalChain(txtfile.Data());
//  if (!chainESD) {
//    Error("CreateLocalChain", "Failed to create chain from file %s", txtfile.Data());
//    return -1;
//  }
//
//  // Create and configure the alien handler plugin
//  AliESDInputHandler* esdH = new AliESDInputHandler();
//  //  esdH->SetNeedField(kTRUE);
//  mgr->SetInputEventHandler(esdH);
//
//  AliMCEventHandler* handler = new AliMCEventHandler;
//  handler->SetReadTR(kFALSE);
//  mgr->SetMCtruthEventHandler(handler);
//
//  AnliAnalysisTaskCheckMCCountsESD* taskeff = reinterpret_cast<AliAnalysisTaskCheckMCCountsESD*>(gInterpreter->ProcessLine(".x AddTaskCheckMCCountsESD()", "CheckMCCountsESD.C"));
//
//  mgr->InitAnalysis();
//  mgr->PrintStatus();
//  return mgr->StartAnalysis("local", chainESD);
//}

int checkMCCounts(TString filename = "O2_AOD.root", bool isAOD = true)
{
  gROOT->SetBatch(kTRUE);

  //TFile* output = TFile::Open("counts.root", "RECREATE");
  TCanvas* canHis = new TCanvas("Counts", "Counts", 800, 600);

  TH1F* histCounts = new TH1F("PDG counts", "PDG counts", 4, 0, 4);
  histCounts->GetXaxis()->SetBinLabel(1, "el");
  histCounts->GetXaxis()->SetBinLabel(2, "ka");
  histCounts->GetXaxis()->SetBinLabel(3, "pi");
  histCounts->GetXaxis()->SetBinLabel(4, "pr");

  int ret = 0;
  if (isAOD) {
    ret = checkCountsAOD(filename, histCounts);
  } else {
    ret = checkCountsESD(filename, histCounts);
  }

  if (!ret) {
    return ret;
  }

  canHis->SaveAs("mccounts.png");
  return 0;
}
