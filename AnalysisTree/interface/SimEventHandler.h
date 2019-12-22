#ifndef SIMEVENTHANDLER_H
#define SIMEVENTHANDLER_H

#include <vector>
#include <unordered_map>
#include "IvyBase.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "SystematicVariations.h"
#include "TH1F.h"


class SimEventHandler : public IvyBase{
public:
  enum EventRandomNumberType{
    kDataPeriod,
    kGenMETSmear
  };

protected:
  std::unordered_map< TString, std::vector<TH1F*> > map_DataPeriod_PUHistList;

  std::unordered_map<EventRandomNumberType, unsigned long long> product_rnds;
  TString theChosenDataPeriod;
  float pileupWeight;

  void setupPUHistograms();
  void clearPUHistograms();

  bool constructRandomNumbers();
  bool constructPUWeight(SystematicsHelpers::SystematicVariationTypes const& syst);

  void clear();

public:
  SimEventHandler();
  ~SimEventHandler();

  bool constructSimEvent(SystematicsHelpers::SystematicVariationTypes const& syst);

  void bookBranches(BaseTree* intree);

  TString const& getChosenDataPeriod() const;
  float const& getPileUpWeight() const{ return pileupWeight; }

  unsigned long long const& getRandomNumberSeed(SimEventHandler::EventRandomNumberType type) const;

};


#endif
