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
  // The map values are vectors for nongap, gap, nongap_gap (combined) histograms separated and in this order.
  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_eff_mc_reco_map;
  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_eff_mc_id_map;
  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_eff_mc_iso_loose_map;
  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_eff_mc_iso_tight_map;

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
