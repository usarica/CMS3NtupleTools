#ifndef SUPERCLUSTERHANDLER_H
#define SUPERCLUSTERHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "SuperclusterObject.h"
#include "ParticleDisambiguator.h"


class SuperclusterHandler : public IvyBase{
public:
  typedef SuperclusterObject ProductType_t;
  static const std::string colName;

protected:
  friend class ParticleDisambiguator;

  std::vector<ProductType_t*> productList;

  void clear(){ this->resetCache(); for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

public:
  // Constructors
  SuperclusterHandler();

  // Destructors
  ~SuperclusterHandler(){ clear(); }

  bool constructSuperclusters(SystematicsHelpers::SystematicVariationTypes const& syst);
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  static void bookBranches(BaseTree* tree);

};


#endif
