#ifndef PUJETIDSCALEFACTORHANDLER_H
#define PUJETIDSCALEFACTORHANDLER_H

#include <unordered_map>
#include <vector>
#include <CMSDataTools/AnalysisTree/interface/ExtendedHistogram_2D.h>
#include "AK4JetObject.h"
#include "ScaleFactorHandlerBase.h"
#include "SystematicVariations.h"


class PUJetIdScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_pujetidwp_effs_map_mistagged;
  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::vector<ExtendedHistogram_2D> > syst_pujetidwp_effs_map_matched;
  std::vector<ExtendedHistogram_2D> pujetidwp_SFs_map_mistagged;
  std::vector<ExtendedHistogram_2D> pujetidwp_SFs_map_matched;

  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& pt, float const& eta, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const;

public:
  PUJetIdScaleFactorHandler();
  ~PUJetIdScaleFactorHandler();

  bool setup();
  void reset();

  void getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& eta, bool const& isMatched, bool const& isLoose, bool const& isMedium, bool const& isTight, float& val, float* effval) const;
  void getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, AK4JetObject const* obj, float& val, float* effval) const;

};



#endif
