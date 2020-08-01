#ifndef MUONSCALEFACTORHANDLER_H
#define MUONSCALEFACTORHANDLER_H

#include <unordered_map>
#include "ScaleFactorHandlerBase.h"
#include "SystematicVariations.h"
#include "MuonObject.h"


class MuonScaleFactorHandler : public ScaleFactorHandlerBase{
public:
  enum EfficiencyType{
    kIdEff = 0,
    kLooseIsoEff,
    kTightIsoEff,
    kAllEffs,

    nEfficiencyTypes
  };

protected:
  ExtendedHistogram_2D eff_mc_id_hists;
  ExtendedHistogram_2D eff_mc_iso_loose_hists;
  ExtendedHistogram_2D eff_mc_iso_tight_hists;

  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_SF_id_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_SF_iso_loose_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_SF_iso_tight_map;

  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& pt, float const& eta, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const;

public:
  MuonScaleFactorHandler();
  ~MuonScaleFactorHandler();

  bool setup();
  void reset();

  void getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& eta, MuonScaleFactorHandler::EfficiencyType type, float& val, float* effval) const;
  void getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, MuonObject const* obj, MuonScaleFactorHandler::EfficiencyType type, float& val, float* effval) const;

};



#endif
