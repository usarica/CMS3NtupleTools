#ifndef ELECTRONHANDLER_H
#define ELECTRONHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "ElectronObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "PFCandidateObject.h"
#include "OverlapMapHandler.h"
#include "ParticleDisambiguator.h"
#include "SystematicVariations.h"


class ElectronHandler : public IvyBase{
public:
  typedef ElectronObject ProductType_t;
  static const std::string colName;

protected:
  friend class ParticleDisambiguator;

  bool has_mvaid_extras;
  bool has_genmatching;
  bool hasOverlapMaps;
  OverlapMapHandler<ElectronObject, AK4JetObject>* overlapMap_electrons_ak4jets;
  OverlapMapHandler<ElectronObject, AK8JetObject>* overlapMap_electrons_ak8jets;

  std::vector<ProductType_t*> productList;

  void clear(){ this->resetCache(); for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

  bool constructElectronObjects(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool associatePFCandidates(std::vector<PFCandidateObject*> const* pfcandidates) const;
  bool linkOverlapElements() const;

  static void checkOptionalInfo(BaseTree* tree, bool& flag_mvaid_extras, bool& flag_genmatching);

public:
  // Constructors
  ElectronHandler();

  // Destructors
  ~ElectronHandler(){ clear(); }

  bool constructElectrons(SystematicsHelpers::SystematicVariationTypes const& syst, std::vector<PFCandidateObject*> const* pfcandidates);

  bool wrapTree(BaseTree* tree);

  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  void bookBranches(BaseTree* tree);

  void registerOverlapMaps(
    OverlapMapHandler<ElectronObject, AK4JetObject>& overlapMap_electrons_ak4jets_,
    OverlapMapHandler<ElectronObject, AK8JetObject>& overlapMap_electrons_ak8jets_
  );

};


#endif
