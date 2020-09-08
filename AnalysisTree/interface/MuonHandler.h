#ifndef MUONHANDLER_H
#define MUONHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "MuonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "PFCandidateObject.h"
#include "OverlapMapHandler.h"
#include "ParticleDisambiguator.h"
#include "SystematicVariations.h"


class MuonHandler : public IvyBase{
public:
  typedef MuonObject ProductType_t;
  static const std::string colName;

protected:
  friend class ParticleDisambiguator;

  bool has_precomputed_timing;
  bool has_genmatching;
  bool hasOverlapMaps;
  OverlapMapHandler<MuonObject, AK4JetObject>* overlapMap_muons_ak4jets;
  OverlapMapHandler<MuonObject, AK8JetObject>* overlapMap_muons_ak8jets;

  std::vector<ProductType_t*> productList;

  void clear(){ this->resetCache(); for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

  bool constructMuonObjects(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool associatePFCandidates(std::vector<PFCandidateObject*> const* pfcandidates) const;
  bool linkOverlapElements() const;

  static void checkOptionalInfo(BaseTree* tree, bool& flag_precomputed_timing, bool& flag_genmatching);

public:
  // Constructors
  MuonHandler();

  // Destructors
  ~MuonHandler(){ clear(); }

  bool constructMuons(SystematicsHelpers::SystematicVariationTypes const& syst, std::vector<PFCandidateObject*> const* pfcandidates);

  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  bool wrapTree(BaseTree* tree);

  void bookBranches(BaseTree* tree);

  void registerOverlapMaps(
    OverlapMapHandler<MuonObject, AK4JetObject>& overlapMap_muons_ak4jets_,
    OverlapMapHandler<MuonObject, AK8JetObject>& overlapMap_muons_ak8jets_
  );

};


#endif
