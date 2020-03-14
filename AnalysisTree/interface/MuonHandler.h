#ifndef MUONHANDLER_H
#define MUONHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "MuonObject.h"
#include "ParticleDisambiguator.h"
#include "SystematicVariations.h"


class MuonHandler : public IvyBase{
public:
  typedef MuonObject ProductType_t;
  static const std::string colName;

protected:
  friend class ParticleDisambiguator;

  bool has_precomputed_timing_flag;
  std::vector<ProductType_t*> productList;

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

public:
  // Constructors
  MuonHandler();

  // Destructors
  ~MuonHandler(){ clear(); }

  bool constructMuons(SystematicsHelpers::SystematicVariationTypes const& syst);
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  bool wrapTree(BaseTree* tree);

  void bookBranches(BaseTree* tree);

};


#endif
