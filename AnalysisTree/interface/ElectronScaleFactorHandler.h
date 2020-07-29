#ifndef ELECTRONSCALEFACTORHANDLER_H
#define ELECTRONSCALEFACTORHANDLER_H

#include "ScaleFactorHandlerBase.h"
#include "SystematicVariations.h"
#include "ElectronObject.h"


class ElectronScaleFactorHandler : public ScaleFactorHandlerBase{
public:
  enum EfficiencyType{
    kTrackingEff,
    kIdEff,
    kLooseIsoEff,
    kTightIsoEff,
    kAll
  };

protected:
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_eff_mc_id_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_eff_mc_iso_loose_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_eff_mc_iso_tight_map;

  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_eff_data_id_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_eff_data_iso_loose_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_eff_data_iso_tight_map;

  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_SF_id_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_SF_iso_loose_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, ExtendedHistogram_2D> syst_SF_iso_tight_map;

  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& pt, float const& etaSC, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const;
  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, ElectronObject const* obj, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const;

public:
  ElectronScaleFactorHandler();
  ~ElectronScaleFactorHandler();

  bool setup();
  void reset();

  void getEffAndError(float& theEff, float& theEffRelErr, float const& pt, float const& etaSC, bool isData, bool useFastSim, ElectronScaleFactorHandler::EfficiencyType type) const;
  void getSFAndError(float& theSF, float& theSFRelErr, float const& pt, float const& etaSC, bool useFastSim, ElectronScaleFactorHandler::EfficiencyType type) const;

  void getEffAndError(float& theEff, float& theEffRelErr, ElectronObject const* obj, bool isData, bool useFastSim, ElectronScaleFactorHandler::EfficiencyType type) const;
  void getSFAndError(float& theSF, float& theSFRelErr, ElectronObject const* obj, bool useFastSim, ElectronScaleFactorHandler::EfficiencyType type) const;

};



#endif
