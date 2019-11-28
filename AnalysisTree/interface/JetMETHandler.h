#ifndef JETMETHANDLER_H
#define JETMETHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "SystematicVariations.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "METObject.h"


class JetMETHandler : public IvyBase{
public:
  static const std::string colName_ak4jets;
  static const std::string colName_ak8jets;
  static const std::string colName_pfchsmet;
  static const std::string colName_pfpuppimet;

protected:
  std::vector<AK4JetObject*> ak4jets;
  std::vector<AK8JetObject*> ak8jets;
  METObject* pfchsmet;
  METObject* pfpuppimet;

  void clear();

public:
  // Constructors
  JetMETHandler();

  // Destructors
  ~JetMETHandler(){ clear(); }

  bool constructJetMET(SystematicsHelpers::SystematicVariationTypes const& syst);

  std::vector<AK4JetObject*> const& getAK4Jets() const{ return ak4jets; }
  std::vector<AK8JetObject*> const& getAK8Jets() const{ return ak8jets; }
  METObject* const& getPFCHSMET() const{ return pfchsmet; }
  METObject* const& getPFPUPPIMET() const{ return pfpuppimet; }

  static void bookBranches(BaseTree* tree);

};


#endif
