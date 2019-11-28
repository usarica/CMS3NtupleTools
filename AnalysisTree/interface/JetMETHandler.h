#ifndef JETMETHANDLER_H
#define JETMETHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "SystematicVariations.h"


class JetMETHandler : public IvyBase{
public:
  static const std::string colName_ak4jets;
  static const std::string colName_ak8jets;

protected:
  std::vector<AK4JetObject*> ak4jets;
  std::vector<AK8JetObject*> ak8jets;

  void clear();
  //void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

public:
  // Constructors
  JetMETHandler();

  // Destructors
  ~JetMETHandler(){ clear(); }

  bool constructProducts(SystematicsHelpers::SystematicVariationTypes const& syst);
  std::vector<AK4JetObject*> const& getAK4Jets() const{ return ak4jets; }
  std::vector<AK8JetObject*> const& getAK8Jets() const{ return ak8jets; }

  static void bookBranches(BaseTree* tree);

};


#endif
