#ifndef PHOTONSCALEFACTORHANDLER_H
#define PHOTONSCALEFACTORHANDLER_H

#include "ExtendedHistogram_2D.h"
#include "ScaleFactorHandlerBase.h"
#include "PhotonObject.h"
#include "PhotonSelectionHelpers.h"


class PhotonScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  ExtendedHistogram_2D_f h_eff_mc_tampon;
  ExtendedHistogram_2D_f h_eff_mc_tight;

  ExtendedHistogram_2D_f h_SF_tampon;
  ExtendedHistogram_2D_f h_SF_tight;

  void evalScaleFactorFromHistogram(float& val, float& relerr, float const& pt, float const& etaSC, ExtendedHistogram_2D_f const& hist, bool etaOnY, bool useAbsEta) const;

  static TString getScaleFactorFileName(PhotonSelectionHelpers::SelectionBits const& preselectionBit, int const& year);

public:
  PhotonScaleFactorHandler();
  ~PhotonScaleFactorHandler();

  bool setup();
  void reset();

  void getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& etaSC, bool const& isTight, bool const& isTampon, float& val, float* effval) const;
  void getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, PhotonObject const* obj, float& val, float* effval) const;

};



#endif
