#ifndef JETMETHANDLER_H
#define JETMETHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "SystematicVariations.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
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

  bool constructAK4Jets(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool constructAK8Jets(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool constructMET(SystematicsHelpers::SystematicVariationTypes const& syst);

  bool applyJetCleaning(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons);

public:
  // Constructors
  JetMETHandler();

  // Destructors
  ~JetMETHandler(){ clear(); }

  bool constructJetMET(SystematicsHelpers::SystematicVariationTypes const& syst, std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons);

  std::vector<AK4JetObject*> const& getAK4Jets() const{ return ak4jets; }
  std::vector<AK8JetObject*> const& getAK8Jets() const{ return ak8jets; }
  METObject* const& getPFCHSMET() const{ return pfchsmet; }
  METObject* const& getPFPUPPIMET() const{ return pfpuppimet; }

  static void bookBranches(BaseTree* tree);

};


#endif
