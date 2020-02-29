#ifndef ELECTRONSCALEFACTORHANDLER_H
#define ELECTRONSCALEFACTORHANDLER_H

#include "TH2F.h"
#include "TFile.h"
#include "ScaleFactorHandlerBase.h"
#include "ElectronObject.h"


class ElectronScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  ExtendedHistogram_2D h_eff_mc_id;
  ExtendedHistogram_2D h_eff_mc_iso;

  ExtendedHistogram_2D h_eff_data_id;
  ExtendedHistogram_2D h_eff_data_iso;

  ExtendedHistogram_2D h_SF_id;
  ExtendedHistogram_2D h_SF_iso;

  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& pt, float const& etaSC, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const;
  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, ElectronObject const* obj, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const;

public:
  ElectronScaleFactorHandler();
  ~ElectronScaleFactorHandler();

  bool setup();
  void reset();

  void getIdIsoEffAndError(float& theEff, float& theEffRelErr, float const& pt, float const& etaSC, bool isData, bool useFastSim) const;
  void getIdIsoSFAndError(float& theSF, float& theSFRelErr, float const& pt, float const& etaSC, bool useFastSim) const;

  void getIdIsoEffAndError(float& theEff, float& theEffRelErr, ElectronObject const* obj, bool isData, bool useFastSim) const;
  void getIdIsoSFAndError(float& theSF, float& theSFRelErr, ElectronObject const* obj, bool useFastSim) const;

};



#endif
