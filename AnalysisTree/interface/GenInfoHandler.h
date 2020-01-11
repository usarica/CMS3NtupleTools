#ifndef GENINFOHANDLER_H
#define GENINFOHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "GenInfoObject.h"
#include "LHEParticleObject.h"
#include "SystematicVariations.h"


class GenInfoHandler : public IvyBase{
protected:
  std::unordered_map< BaseTree*, std::vector<TString> > tree_MElist_map;
  std::unordered_map< BaseTree*, bool > tree_lheparticles_present_map;

  bool acquireCoreGenInfo;
  bool acquireLHEMEWeights;
  bool acquireLHEParticles;
  bool acquireGenParticles;

  GenInfoObject* genInfo;
  std::vector<LHEParticleObject*> lheparticles;

  bool constructCoreGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool constructLHEParticles();

  void clear(){ delete genInfo; genInfo=nullptr; for (auto*& part:lheparticles) delete part; lheparticles.clear(); }

public:
  static const std::string colName_lheparticles;
  static const std::string colName_genparticles;

  GenInfoHandler();
  ~GenInfoHandler(){ clear(); }

  bool constructGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst);

  GenInfoObject* const& getGenInfo() const{ return genInfo; }
  std::vector<LHEParticleObject*> const& getLHEParticles() const{ return lheparticles; }

  void setAcquireCoreGenInfo(bool flag){ acquireCoreGenInfo=flag; }
  void setAcquireLHEMEWeights(bool flag){ acquireLHEMEWeights=flag; }
  void setAcquireLHEParticles(bool flag){ acquireLHEParticles=flag; }
  void setAcquireGenParticles(bool flag){ acquireGenParticles=flag; }

  void bookBranches(BaseTree* tree);

};


#endif
