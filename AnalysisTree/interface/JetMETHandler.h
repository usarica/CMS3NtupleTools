#ifndef JETMETHANDLER_H
#define JETMETHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "SystematicVariations.h"
#include "PFCandidateObject.h"
#include "OverlapMapHandler.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "METObject.h"


class JetMETHandler : public IvyBase{
public:
  static const std::string& colName_ak4jets;
  static const std::string& colName_ak8jets;
  static const std::string colName_pfmet;
  static const std::string colName_pfpuppimet;
  static const std::string colName_vertices;

protected:
  bool hasOverlapMaps;
  OverlapMapHandler<MuonObject, AK4JetObject>* overlapMap_muons_ak4jets;
  OverlapMapHandler<MuonObject, AK8JetObject>* overlapMap_muons_ak8jets;
  OverlapMapHandler<ElectronObject, AK4JetObject>* overlapMap_electrons_ak4jets;
  OverlapMapHandler<ElectronObject, AK8JetObject>* overlapMap_electrons_ak8jets;
  OverlapMapHandler<PhotonObject, AK4JetObject>* overlapMap_photons_ak4jets;
  OverlapMapHandler<PhotonObject, AK8JetObject>* overlapMap_photons_ak8jets;

  std::vector<AK4JetObject*> ak4jets;
  std::vector<AK8JetObject*> ak8jets;
  METObject* pfmet;
  METObject* pfpuppimet;

  float pfmet_XYcorr_xCoeffA; float pfmet_XYcorr_xCoeffB;
  float pfmet_XYcorr_yCoeffA; float pfmet_XYcorr_yCoeffB;

  void clear();

  bool constructAK4Jets(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool constructAK8Jets(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool associatePFCandidates(std::vector<PFCandidateObject*> const* pfcandidates) const;
  bool linkOverlapElements() const;
  bool applyJetCleaning(bool usePFCandidates, std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons);

  bool constructMET(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool assignMETXYShifts(SystematicsHelpers::SystematicVariationTypes const& syst);
  bool applyMETParticleShifts(bool usePFCandidates, std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons);

public:
  // Constructors
  JetMETHandler();

  // Destructors
  ~JetMETHandler(){ clear(); }

  bool constructJetMET(
    SystematicsHelpers::SystematicVariationTypes const& syst,
    std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons,
    std::vector<PFCandidateObject*> const* pfcandidates
  );

  std::vector<AK4JetObject*> const& getAK4Jets() const{ return ak4jets; }
  std::vector<AK8JetObject*> const& getAK8Jets() const{ return ak8jets; }
  METObject* const& getPFMET() const{ return pfmet; }
  METObject* const& getPFPUPPIMET() const{ return pfpuppimet; }

  bool wrapTree(BaseTree* tree);

  void bookBranches(BaseTree* tree);

  void registerOverlapMaps(
    OverlapMapHandler<MuonObject, AK4JetObject>& overlapMap_muons_ak4jets_,
    OverlapMapHandler<MuonObject, AK8JetObject>& overlapMap_muons_ak8jets_,
    OverlapMapHandler<ElectronObject, AK4JetObject>& overlapMap_electrons_ak4jets_,
    OverlapMapHandler<ElectronObject, AK8JetObject>& overlapMap_electrons_ak8jets_,
    OverlapMapHandler<PhotonObject, AK4JetObject>& overlapMap_photons_ak4jets_,
    OverlapMapHandler<PhotonObject, AK8JetObject>& overlapMap_photons_ak8jets_
  );

};


#endif
