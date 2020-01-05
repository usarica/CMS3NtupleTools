#ifndef SYSTEMATICVARIATIONS_H
#define SYSTEMATICVARIATIONS_H

#include <string>


namespace SystematicsHelpers{

  enum SystematicVariationTypes{
    sNominal,

    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp,
    tPythiaScaleDn, tPythiaScaleUp,
    tPythiaTuneDn, tPythiaTuneUp,

    eEleEffDn, eEleEffUp,
    eEleScaleDn, eEleScaleUp,
    eEleResDn, eEleResUp,

    eMuEffDn, eMuEffUp,
    eMuScaleDn, eMuScaleUp,
    eMuResDn, eMuResUp,

    ePhoEffDn, ePhoEffUp,
    ePhoScaleDn, ePhoScaleUp,
    ePhoResDn, ePhoResUp,

    eMETDn, eMETUp,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp,
    eBTagSFDn, eBTagSFUp,

    nSystematicVariations,

    sUncorrected // For checks
  };

  std::string getSystCoreName(SystematicsHelpers::SystematicVariationTypes const& type);
  std::string getSystName(SystematicsHelpers::SystematicVariationTypes const& type);

}

#endif
