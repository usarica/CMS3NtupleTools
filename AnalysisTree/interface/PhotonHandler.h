#ifndef PHOTONHANDLER_H
#define PHOTONHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "SystematicVariations.h"


class PhotonHandler : public IvyBase{
public:
  typedef PhotonObject ProductType_t;
  static const std::string colName;

protected:
  std::vector<ProductType_t*> productList;

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

  bool applyCleaning(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons);

public:
  // Constructors
  PhotonHandler();

  // Destructors
  ~PhotonHandler(){ clear(); }

  bool constructPhotons(SystematicsHelpers::SystematicVariationTypes const& syst, std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons);
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  static void bookBranches(BaseTree* tree);

};


#endif
