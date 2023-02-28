#ifndef ALIANALYSISTASKCHECKMCCOUNTSESD
#define ALIANALYSISTASKCHECKMCCOUNTSESD

#include <TH1.h>

#include <fstream>

#include "AliAnalysisTaskSE.h"

class TList;

class AliAnalysisTaskCheckMCCountsESD : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskCheckMCCountsESD();
  virtual ~AliAnalysisTaskCheckMCCountsESD();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t*);
  virtual void Terminate(Option_t*);

  int GetBinNumber(Int_t);

 private:
  AliAnalysisTaskCheckMCCountsESD(const AliAnalysisTaskCheckMCCountsESD& source);
  AliAnalysisTaskCheckMCCountsESD& operator=(const AliAnalysisTaskCheckMCCountsESD& source);
  TList* fOutputList;
  TH1F* fHistCounts;
  std::ofstream fDebug;
  int fInd;

  ClassDef(AliAnalysisTaskCheckMCCountsESD, 1);
};
#endif

