#ifndef SYSTEMATICVARIATIONS_H
#define SYSTEMATICVARIATIONS_H


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

    eMETDn, eMETUp,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp,
    eBTagSFDn, eBTagSFUp,

    nSystematicVariations
  };

}

#endif
