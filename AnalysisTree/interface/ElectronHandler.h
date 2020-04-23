#ifndef ELECTRONHANDLER_H
#define ELECTRONHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "ElectronObject.h"
#include "ParticleDisambiguator.h"
#include "SystematicVariations.h"


class ElectronHandler : public IvyBase{
public:
  typedef ElectronObject ProductType_t;
  static const std::string colName;

protected:
  friend class ParticleDisambiguator;

  bool has_mvaid_extras;
  std::vector<ProductType_t*> productList;

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

public:
  // Constructors
  ElectronHandler();

  // Destructors
  ~ElectronHandler(){ clear(); }

  bool constructElectrons(SystematicsHelpers::SystematicVariationTypes const& syst);

  bool wrapTree(BaseTree* tree);

  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  void bookBranches(BaseTree* tree);

};


#endif
