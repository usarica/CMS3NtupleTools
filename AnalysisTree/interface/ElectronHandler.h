#ifndef ELECTRONHANDLER_H
#define ELECTRONHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "ElectronObject.h"
#include "SystematicVariations.h"


class ElectronHandler : public IvyBase{
public:
  typedef ElectronObject ProductType_t;
  static const std::string colName;

protected:
  std::vector<ProductType_t*> productList;

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

public:
  // Constructors
  ElectronHandler();

  // Destructors
  ~ElectronHandler(){ clear(); }

  bool constructElectrons(SystematicsHelpers::SystematicVariationTypes const& syst);
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  static void bookBranches(BaseTree* tree);

};


#endif
