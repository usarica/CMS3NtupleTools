#ifndef PHOTONHANDLER_H
#define PHOTONHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "PFCandidateObject.h"
#include "OverlapMapHandler.h"
#include "ParticleDisambiguator.h"
#include "SystematicVariations.h"


class PhotonHandler : public IvyBase{
public:
  typedef PhotonObject ProductType_t;
  static const std::string colName;

protected:
  friend class ParticleDisambiguator;

  bool has_genmatching;
  bool hasOverlapMaps;
  OverlapMapHandler<PhotonObject, AK4JetObject>* overlapMap_photons_ak4jets;
  OverlapMapHandler<PhotonObject, AK8JetObject>* overlapMap_photons_ak8jets;

  std::vector<ProductType_t*> productList;

  void clear(){ this->resetCache(); for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

  bool constructPhotonObjects(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool associatePFCandidates(std::vector<PFCandidateObject*> const* pfcandidates) const;
  bool linkOverlapElements() const;

  static void checkOptionalInfo(BaseTree* tree, bool& flag_genmatching);

public:
  // Constructors
  PhotonHandler();

  // Destructors
  ~PhotonHandler(){ clear(); }

  bool constructPhotons(SystematicsHelpers::SystematicVariationTypes const& syst, std::vector<PFCandidateObject*> const* pfcandidates);

  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  bool wrapTree(BaseTree* tree);

  void bookBranches(BaseTree* tree);

  void registerOverlapMaps(
    OverlapMapHandler<PhotonObject, AK4JetObject>& overlapMap_photons_ak4jets_,
    OverlapMapHandler<PhotonObject, AK8JetObject>& overlapMap_photons_ak8jets_
  );

};


#endif
