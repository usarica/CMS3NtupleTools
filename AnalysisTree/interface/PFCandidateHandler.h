#ifndef PFCANDIDATEHANDLER_H
#define PFCANDIDATEHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "PFCandidateObject.h"
#include "SystematicVariations.h"


class PFCandidateHandler : public IvyBase{
public:
  typedef PFCandidateObject ProductType_t;
  static const std::string colName;

protected:
  std::vector<ProductType_t*> productList;

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

public:
  // Constructors
  PFCandidateHandler();

  // Destructors
  ~PFCandidateHandler(){ clear(); }

  bool constructPFCandidates(SystematicsHelpers::SystematicVariationTypes const& syst);
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  static void bookBranches(BaseTree* tree);

};


#endif
