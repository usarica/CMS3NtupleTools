#ifndef METCORRECTIONHANDLER_H
#define METCORRECTIONHANDLER_H

#include <vector>
#include <unordered_map>
#include <utility>
#include "METObject.h"
#include "SystematicVariations.h"
#include "SimEventHandler.h"
#include "ScaleFactorHandlerBase.h"


class METCorrectionParameters{
public:
  std::vector<float> sigmas_nominal;
  std::vector<float> sigmas_dn;
  std::vector<float> sigmas_up;
  std::vector<float> fracs;

  METCorrectionParameters();
  METCorrectionParameters(METCorrectionParameters const& other);
  ~METCorrectionParameters(){}

  METCorrectionParameters& operator=(METCorrectionParameters const& other);
  void swap(METCorrectionParameters& other);

  void translateFractions();

};

class METCorrectionHandler : public ScaleFactorHandlerBase{
protected:
  std::unordered_map<
    TString,
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
  > pfmet_XY_JER_map;
  std::unordered_map<
    TString,
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
  > pfmet_JER_map;
  std::unordered_map<
    TString,
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
  > pfmet_XY_map;
  std::unordered_map<
    TString,
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
  > pfmet_map;
  std::unordered_map<
    TString,
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
  > puppimet_map;

  void readFile(TString const& strinput, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>& pars);

public:
  METCorrectionHandler();
  ~METCorrectionHandler();

  bool setup();
  void reset();

  void applyCorrections(
    TString const& effDataPeriod,
    float const& genMET, float const& genMETPhi,
    METObject* obj, bool isPFMET
  ) const;

  void printParameters() const;

};


#endif