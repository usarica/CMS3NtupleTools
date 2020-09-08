#ifndef ISOTRACKHANDLER_H
#define ISOTRACKHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "IsotrackObject.h"
#include "ParticleDisambiguator.h"
#include "SystematicVariations.h"


class IsotrackHandler : public IvyBase{
public:
  typedef IsotrackObject ProductType_t;
  static const std::string colName;

protected:
  friend class ParticleDisambiguator;

  std::vector<ProductType_t*> productList;

  void clear(){ this->resetCache(); for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

  bool applyCleaning(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons);

public:
  // Constructors
  IsotrackHandler();

  // Destructors
  ~IsotrackHandler(){ clear(); }

  bool constructIsotracks(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons);

  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  static void bookBranches(BaseTree* tree);

};


#endif
