#ifndef MUONSCALEFACTORHANDLER_H
#define MUONSCALEFACTORHANDLER_H

#include <unordered_map>
#include "ScaleFactorHandlerBase.h"
#include "SystematicVariations.h"
#include "MuonObject.h"


class MuonScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D_f> syst_eff_mc_id_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D_f> syst_eff_mc_iso_loose_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D_f> syst_eff_mc_iso_tight_map;

  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D_f> syst_SF_id_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D_f> syst_SF_iso_loose_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D_f> syst_SF_iso_tight_map;

  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& pt, float const& eta, ExtendedHistogram_2D_f const& hist, bool etaOnY, bool useAbsEta) const;

public:
  MuonScaleFactorHandler();
  ~MuonScaleFactorHandler();

  bool setup();
  void reset();

  void getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& eta, bool const& passId, bool const& passLooseIso, bool const& passTightIso, float& val, float* effval) const;
  void getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, MuonObject const* obj, float& val, float* effval) const;

};



#endif
