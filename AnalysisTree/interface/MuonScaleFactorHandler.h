#ifndef MUONSCALEFACTORHANDLER_H
#define MUONSCALEFACTORHANDLER_H

#include "TH2F.h"
#include "TFile.h"
#include "ScaleFactorHandlerBase.h"
#include "MuonObject.h"


class MuonScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  TFile* finput_SF_id;
  TFile* finput_SF_iso;

  TH2F* h_SF_id;
  TH2F* h_SF_iso;

  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, MuonObject const* obj, TH2F const* hist, bool etaOnY, bool useAbsEta) const;

public:
  MuonScaleFactorHandler();
  ~MuonScaleFactorHandler();

  bool setup();
  void reset();

  void getIdIsoSFAndError(float& theSF, float& theSFRelErr, MuonObject const* obj, bool useFastSim) const;

};



#endif
