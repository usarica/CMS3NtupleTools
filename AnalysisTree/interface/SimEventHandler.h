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
  // There are two versions of the data period random numbers in order to save headache:
  // - kDataPeriod_global, defined with respect to the entire year.
  // - kDataPeriod_local, defined with respect to the chosen data period.
  enum EventRandomNumberType{
    kDataPeriod_global,
    kDataPeriod_local,
    kGenMETSmear
  };

protected:
  std::unordered_map< TString, std::vector<TH1F*> > map_DataPeriod_PUHistList;
  std::unordered_map<TString, TH1F*> map_exceptionalPUHistList;

  bool hasPUException;
  TString theChosenDataPeriod;
  // These are the seeds, not the random numbers themselves:
  std::unordered_map<EventRandomNumberType, unsigned long long> product_rnds;
  // These are in fact the random numbers:
  std::unordered_map<EventRandomNumberType, double> product_rndnums;
  bool hasHEM2018Issue;
  std::vector<float> pileupWeights; // sNominal, ePUDn, ePUUp
  std::vector<float const*> l1prefiringWeights; // sNominal, eL1PrefiringDn, eL1PrefiringUp

  void setupPUHistograms();
  void clearPUHistograms();

  bool constructRandomNumbers();
  bool constructPUWeight();
  bool constructL1PrefiringWeight();

  void clear();

public:
  SimEventHandler();
  ~SimEventHandler();

  bool constructSimEvent();

  bool wrapTree(BaseTree* tree);

  void bookBranches(BaseTree* intree);

  TString const& getChosenDataPeriod() const;
  int getChosenRunNumber() const;
  bool const& getHasHEM2018Issue() const{ return hasHEM2018Issue; }

  float const& getPileUpWeight(SystematicsHelpers::SystematicVariationTypes const& syst) const;
  float getL1PrefiringWeight(SystematicsHelpers::SystematicVariationTypes const& syst) const;

  unsigned long long const& getRandomNumberSeed(SimEventHandler::EventRandomNumberType type) const;
  double const& getRandomNumber(SimEventHandler::EventRandomNumberType type) const;

};


#endif
