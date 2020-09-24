#ifndef ELECTRONSCALEFACTORHANDLER_H
#define ELECTRONSCALEFACTORHANDLER_H

#include <unordered_map>
#include "ScaleFactorHandlerBase.h"
#include "SystematicVariations.h"
#include "ElectronObject.h"


class ElectronScaleFactorHandler : public ScaleFactorHandlerBase{
public:
  enum EfficiencyType{
    kTrackingEff = 0,
    kIdEff,
    kLooseIsoEff,
    kTightIsoEff,
    kAllEffs,

    nEfficiencyTypes
  };

protected:
  std::vector<ExtendedHistogram_2D> eff_mc_reco_hists;
  std::vector<ExtendedHistogram_2D> eff_mc_id_hists;
  std::vector<ExtendedHistogram_2D> eff_mc_iso_loose_hists;
  std::vector<ExtendedHistogram_2D> eff_mc_iso_tight_hists;

  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_SF_reco_map;
  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_SF_id_map;
  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_SF_iso_loose_map;
  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_SF_iso_tight_map;

  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& pt, float const& etaSC, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const;

public:
  ElectronScaleFactorHandler();
  ~ElectronScaleFactorHandler();

  bool setup();
  void reset();

  // idx_gap==0: Non-gap, ==1: gap, ==2: combined
  void getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& etaSC, unsigned short const& idx_gap, bool const& passId, bool const& passLooseIso, bool const& passTightIso, float& val, float* effval) const;
  void getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, ElectronObject const* obj, float& val, float* effval) const;

};



#endif
