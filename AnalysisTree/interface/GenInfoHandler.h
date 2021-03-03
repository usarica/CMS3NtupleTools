#ifndef GENINFOHANDLER_H
#define GENINFOHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "SampleExceptions.h"
#include "GenInfoObject.h"
#include "LHEParticleObject.h"
#include "GenParticleObject.h"
#include "GenJetObject.h"
#include "SystematicVariations.h"


class GenInfoHandler : public IvyBase{
protected:
  std::unordered_map< BaseTree*, std::vector<TString> > tree_kfactorlist_map;
  std::unordered_map< BaseTree*, std::vector<TString> > tree_MElist_map;
  std::unordered_map< BaseTree*, bool > tree_lheparticles_present_map;
  std::unordered_map< BaseTree*, float > abs_genWeight_default_thr_map;

  bool acquireCoreGenInfo;
  bool acquireLHEMEWeights;
  bool acquireLHEParticles;
  bool acquireGenParticles;
  bool acquireGenAK4Jets;
  bool acquireGenAK8Jets;
  bool doGenJetsCleaning; // Needs gen. particles and LHE taus
  bool doGenJetsVDecayCleaning; // Needs LHE particles
  bool allowLargeGenWeightRemoval;

  SampleHelpers::GenWeightExceptionType genWeightException;
  float const* abs_genWeight_default_thr;

  GenInfoObject* genInfo;
  std::vector<LHEParticleObject*> lheparticles;
  std::vector<GenParticleObject*> genparticles;
  std::vector<GenJetObject*> genak4jets;
  std::vector<GenJetObject*> genak8jets;

  bool constructCoreGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool constructLHEParticles();
  bool constructGenParticles();
  bool constructGenAK4Jets();
  bool constructGenAK8Jets();

  bool determineWeightThresholds();

  bool applyGenJetCleaning();

  void clear();

public:
  static const std::string colName_lheparticles;
  static const std::string colName_genparticles;
  static const std::string colName_genak4jets;
  static const std::string colName_genak8jets;

  GenInfoHandler();
  ~GenInfoHandler(){ clear(); }

  bool constructGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst);

  GenInfoObject* const& getGenInfo() const{ return genInfo; }
  std::vector<LHEParticleObject*> const& getLHEParticles() const{ return lheparticles; }
  std::vector<GenParticleObject*> const& getGenParticles() const{ return genparticles; }
  std::vector<GenJetObject*> const& getGenAK4Jets() const{ return genak4jets; }
  std::vector<GenJetObject*> const& getGenAK8Jets() const{ return genak8jets; }

  void setAcquireCoreGenInfo(bool flag){ acquireCoreGenInfo=flag; }
  void setAcquireLHEMEWeights(bool flag){ acquireLHEMEWeights=flag; }
  void setAcquireLHEParticles(bool flag){ acquireLHEParticles=flag; }
  void setAcquireGenParticles(bool flag){ acquireGenParticles=flag; }
  void setAcquireGenAK4Jets(bool flag){ acquireGenAK4Jets=flag; }
  void setAcquireGenAK8Jets(bool flag){ acquireGenAK8Jets=flag; }
  void setDoGenJetsCleaning(bool flag){ doGenJetsCleaning=flag; }
  void setDoGenJetsVDecayCleaning(bool flag){ doGenJetsVDecayCleaning=flag; }
  void setAllowLargeGenWeightRemoval(bool flag){ allowLargeGenWeightRemoval=flag; }

  bool wrapTree(BaseTree* tree);

  void bookBranches(BaseTree* tree);

};


#endif
