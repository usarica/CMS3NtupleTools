#ifndef TRIGGERSCALEFACTORHANDLER_H
#define TRIGGERSCALEFACTORHANDLER_H

#include <unordered_map>
#include "ScaleFactorHandlerBase.h"
#include "SystematicVariations.h"
#include "ParticleObject.h"


class TriggerScaleFactorHandler : public ScaleFactorHandlerBase{
public:
  enum CombinedTriggerType{
    kCombined_Dilepton,
    kCombined_SingleLepton,

    nCombinedTriggerTypes
  };

protected:
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_eff_mc_Dilepton_SingleLepton_mumu_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_eff_mc_Dilepton_SingleLepton_ee_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_eff_mc_Dilepton_SingleLepton_mue_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D > syst_eff_mc_SingleMuon_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D > syst_eff_mc_SingleElectron_map;

  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_SF_Dilepton_SingleLepton_mumu_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_SF_Dilepton_SingleLepton_ee_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_SF_Dilepton_SingleLepton_mue_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D > syst_SF_SingleMuon_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D > syst_SF_SingleElectron_map;

  void evalScaleFactorFromHistogram_PtEta(float& theSF, float const& pt, float const& eta, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const;
  void evalScaleFactorFromHistogram_PtPt(float& theSF, float const& pt1, float const& pt2, ExtendedHistogram_2D const& hist) const;

public:
  TriggerScaleFactorHandler();
  ~TriggerScaleFactorHandler();

  bool setup();
  void reset();

  void getCombinedDileptonSFAndEff(
    SystematicsHelpers::SystematicVariationTypes const& syst,
    float pt1, float eta1, cms3_id_t id1,
    float pt2, float eta2, cms3_id_t id2,
    bool passTrigger,
    float& val, float* effval
  ) const;
  void getCombinedDileptonSFAndEff(
    SystematicsHelpers::SystematicVariationTypes const& syst,
    ParticleObject const* obj1, ParticleObject const* obj2,
    bool passTrigger,
    float& val, float* effval
  ) const;

  void getCombinedSingleLeptonSFAndEff(
    SystematicsHelpers::SystematicVariationTypes const& syst,
    float const& pt, float const& eta, cms3_id_t const& partId,
    bool passTrigger,
    float& val, float* effval
  ) const;
  void getCombinedSingleLeptonSFAndEff(
    SystematicsHelpers::SystematicVariationTypes const& syst,
    ParticleObject const* obj,
    bool passTrigger,
    float& val, float* effval
  ) const;

};



#endif
