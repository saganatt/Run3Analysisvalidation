#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TParticle.h>
#include <TPDGCode.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>

#include <fstream>
#include <iostream>

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliRunLoader.h"
#include "AliStack.h"

#include "Framework/DataTypes.h"

R__ADD_INCLUDE_PATH(.)
#include "AliAnalysisTaskCheckMCCountsESD.h"
#endif

int getBinNumber(Int_t pdg)
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

TChain* CreateLocalChain(const char* txtfile)
{
  // Open the file
  ifstream in;
  in.open(txtfile);
  Int_t count = 0;
  // Read the input list of files and add them to the chain
  TString line;
  TChain* chain = new TChain("esdTree");
  while (in.good()) {
    in >> line;
    if (line.IsNull() || line.BeginsWith("#"))
      continue;
    TString esdFile(line);
    TFile* file = TFile::Open(esdFile);
    if (file && !file->IsZombie()) {
      chain->Add(esdFile);
      file->Close();
    } else {
      Error("CreateLocalChain", "Skipping un-openable file: %s", esdFile.Data());
    }
  }
  in.close();
  if (!chain->GetListOfFiles()->GetEntries()) {
    Error("CreateLocalChain", "No file from %s could be opened", txtfile);
    delete chain;
    return nullptr;
  }
  return chain;
}

int checkCountsAOD(TString txtfile = "list_o2.txt")
{
  TChain* chainAOD = CreateLocalChain(txtfile.Data());
  if (!chainAOD) {
    std::cout << "Failed to create chain from file " << txtfile.Data() << std::endl;
    return 1;
  }

  TFile* output = TFile::Open("CountsAOD.root", "RECREATE");
  std::ofstream debug;
  debug.open("debug_aod.txt");

  TH1F* hCounts = new TH1F("PDG counts AOD", "PDG counts AOD", 4, 0, 4);
  hCounts->GetXaxis()->SetBinLabel(1, "el");
  hCounts->GetXaxis()->SetBinLabel(2, "ka");
  hCounts->GetXaxis()->SetBinLabel(3, "pi");
  hCounts->GetXaxis()->SetBinLabel(4, "pr");

  int ind = 0;

  TObjArray* files = chainAOD->GetListOfFiles();
  TIter ifile(files);
  TChainElement* chEl = nullptr;
  while ((chEl = (TChainElement*)ifile())) {
    TFile f(chEl->GetTitle());
    TList* l = f.GetListOfKeys();
    TIter it(l);
    TObject* obj = nullptr;
    while ((obj = it())) {
      if (strcmp(obj->GetName(), "metaData") == 0) {
        continue;
      }
      TString treeName = TString::Format("%s/O2mcparticle_001", obj->GetName());
      TTree* mc = (TTree*)f.Get(treeName);
      if (!mc) {
        std::cout << "Tree is null: " << treeName.Data() << std::endl;
        return 1;
      }

      int pdg;
      mc->SetBranchAddress("fPdgCode", &pdg);
      uint8_t flags;
      mc->SetBranchAddress("fFlags", &flags);

      for (int i = 0; i < mc->GetEntries(); i++) {
        mc->GetEntry(i);
        bool isPhysicalPrimary = (flags & o2::aod::mcparticle::enums::PhysicalPrimary) == o2::aod::mcparticle::enums::PhysicalPrimary;
        debug << "Particle " << ind << " PDG: " << pdg << " physical primary: " << isPhysicalPrimary << std::endl;
        if (isPhysicalPrimary) {
          int bin = getBinNumber(pdg);
          if (bin != -1) {
            debug << "Particle " << ind << " bin: " << bin << " PDG code: " << pdg << " kPion: " << kPiPlus << " kProton: " << kProton << " kElectron: " << kElectron << " kKPlus: " << kKPlus << std::endl;
            hCounts->Fill(bin);
          }
        }
        ind++;
      }
    }
  }

  debug.close();

  TCanvas* canHis = new TCanvas("Counts", "Counts", 800, 600);
  canHis->cd();
  hCounts->Draw();

  TLatex* latexa = new TLatex();
  latexa->SetTextSize(0.04);
  latexa->SetTextFont(42);
  latexa->SetTextAlign(22);
  for (int i = 1; i <= hCounts->GetNbinsX(); i++) {
    std::string binstr = std::to_string((int)(hCounts->GetBinContent(i)));
    float bincen = 0.1f + hCounts->GetXaxis()->GetBinCenter(i) / (hCounts->GetNbinsX() + 1);
    latexa->DrawLatexNDC(bincen, 0.2f, binstr.c_str());
  }

  canHis->SaveAs("mccounts_aod.png");
  output->Write();
  delete canHis;

  return 0;
}

int checkCountsESDGAlice(TString filepath = "/data2/data/Run2/pp_13TeV/sim/LHC20f4a/294774/001/", TH1F* hCounts = nullptr)
{
  AliRunLoader* runLoader = AliRunLoader::Open(TString::Format("%s/%s", filepath.Data(), "galice.root"), "read");
  if (runLoader == 0x0) {
    printf("Error opening %s file \n", filepath.Data());
    return 1;
  }
  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  int nevents = runLoader->GetNumberOfEvents();

  for (int ievent = 0; ievent < nevents; ievent++) {
    runLoader->GetEvent(ievent);
    AliStack* stack = runLoader->Stack();
    Int_t nparticles = stack->GetNtrack();
    for (int i = 0; i < nparticles; i++) {
      TParticle* p = stack->Particle(i);
      if (p->IsPrimary()) {
        int bin = getBinNumber(p->GetPdgCode());
        if (bin != -1) {
          std::cout << "Bin: " << bin << " PDG code: " << p->GetPdgCode() << " kPion: " << kPiPlus << " kProton: " << kProton << " kElectron: " << kElectron << " kKPlus: " << kKPlus << std::endl;
          hCounts->Fill(bin);
        }
      }
    }
  }
  return 0;
}

int checkCountsESD(TString txtfile = "list_ali.txt")
{
  AliAnalysisManager* mgr = new AliAnalysisManager("testAnalysis");
  TChain* chainESD = CreateLocalChain(txtfile.Data());
  if (!chainESD) {
    std::cout << "Failed to create chain from file " << txtfile.Data() << std::endl;
    return 1;
  }

  // Create and configure the alien handler plugin
  AliESDInputHandler* esdH = new AliESDInputHandler();
  //  esdH->SetNeedField(kTRUE);
  mgr->SetInputEventHandler(esdH);

  AliMCEventHandler* handler = new AliMCEventHandler;
  handler->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(handler);

  AliAnalysisTaskCheckMCCountsESD* task = new AliAnalysisTaskCheckMCCountsESD();
  mgr->AddTask(task);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":CountsMC";
  TString listname = "listCountsMC";

  AliAnalysisDataContainer* coutput = mgr->CreateContainer(listname,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputFileName.Data());

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput);

  mgr->InitAnalysis();
  mgr->PrintStatus();
  return mgr->StartAnalysis("local", chainESD);
}

int checkMCCounts(TString filename = "O2_AOD.root", bool isAOD = true)
{
  gROOT->SetBatch(kTRUE);

  if (isAOD) {
    return checkCountsAOD(filename);
  } else {
    return checkCountsESD(filename);
  }

  return 0;
}
