#ifndef BTAGSCALEFACTORHANDLER_H
#define BTAGSCALEFACTORHANDLER_H

#include <unordered_map>
#include <CMSDataTools/AnalysisTree/interface/ExtendedHistogram_2D.h>
#include "BTagCalibrationStandalone.h"
#include "ScaleFactorHandlerBase.h"
#include "AK4JetObject.h"
#include "SystematicVariations.h"
#include "BtagHelpers.h"


class BtagScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  std::unordered_map< SystematicsHelpers::SystematicVariationTypes, std::unordered_map< BTagEntry::JetFlavor, std::vector<std::vector<ExtendedHistogram_2D_f>> > > syst_flav_pujetid_WP_mceffhist_map;

  std::unordered_map<BtagHelpers::BtagWPType, BTagCalibration*> WP_calib_map;
  std::unordered_map<BtagHelpers::BtagWPType, BTagCalibrationReader*> WP_calibreader_map_nominal;
  std::unordered_map<BtagHelpers::BtagWPType, BTagCalibrationReader*> WP_calibreader_map_dn;
  std::unordered_map<BtagHelpers::BtagWPType, BTagCalibrationReader*> WP_calibreader_map_up;

  void evalEfficiencyFromHistogram(float& val, float const& pt, float const& eta, ExtendedHistogram_2D_f const& hist, bool etaOnY, bool useAbsEta) const;

  static void loadBTagCalibrations(BTagCalibrationReader* const& reader, BTagCalibration* const& calibration);

  float getSFFromBTagCalibrationReader(
    BTagCalibrationReader const* calibReader, BTagCalibrationReader const* calibReader_Nominal,
    SystematicsHelpers::SystematicVariationTypes const& syst, BTagEntry::JetFlavor const& flav, float const& pt, float const& eta
  ) const;

public:
  BtagScaleFactorHandler();
  ~BtagScaleFactorHandler();

  bool setup();
  void reset();

  void getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& eta, unsigned short const& pujetidcat, BTagEntry::JetFlavor const& flav, float const& btagval, float& val, float* effval) const;
  void getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, AK4JetObject const* obj, float& val, float* effval) const;

};



#endif
