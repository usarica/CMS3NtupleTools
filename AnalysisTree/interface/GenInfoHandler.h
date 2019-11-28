#ifndef GENINFOHANDLER_H
#define GENINFOHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "GenInfoObject.h"
#include "SystematicVariations.h"


class GenInfoHandler : public IvyBase{
protected:
  GenInfoObject* genInfo;

  void clear(){ delete genInfo; genInfo=nullptr; }

public:
  GenInfoHandler();
  ~GenInfoHandler(){ clear(); }

  bool constructGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst);
  GenInfoObject* const& getGenInfo() const{ return genInfo; }

  static void bookBranches(BaseTree* tree);

};


#endif
