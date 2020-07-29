#ifndef PHOTONHANDLER_H
#define PHOTONHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "ParticleDisambiguator.h"
#include "SystematicVariations.h"


class PhotonHandler : public IvyBase{
public:
  typedef PhotonObject ProductType_t;
  static const std::string colName;

protected:
  friend class ParticleDisambiguator;

  bool has_genmatching;

  std::vector<ProductType_t*> productList;

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

  static void checkOptionalInfo(BaseTree* tree, bool& flag_genmatching);

public:
  // Constructors
  PhotonHandler();

  // Destructors
  ~PhotonHandler(){ clear(); }

  bool constructPhotons(SystematicsHelpers::SystematicVariationTypes const& syst);
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  bool wrapTree(BaseTree* tree);

  void bookBranches(BaseTree* tree);

};


#endif
