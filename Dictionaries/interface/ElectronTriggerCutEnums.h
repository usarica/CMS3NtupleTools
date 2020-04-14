#ifndef ELECTRONCUTENUMS_H
#define ELECTRONCUTENUMS_H

namespace ElectronTriggerCutEnums{

  enum ElectronTriggerCutTypes{
    kCuts_CaloIdL_TrackIdL,

    kCuts_CaloIdL_TrackIdL_IsoVL_v1,
    kCuts_CaloIdL_TrackIdL_IsoVL_v2,
    kCuts_CaloIdL_TrackIdL_IsoVL_v3,

    kCuts_CaloIdL_GsfTrkIdVL,
    kCuts_CaloIdL_MW,

    kCuts_WPLoose,
    kCuts_WPTight_v1,
    kCuts_WPTight_v2,

    kCuts_DoublePhoton,

    kCuts_CustomTriggerSafe_Id_v1,
    kCuts_CustomTriggerSafe_Iso_v1,

    kCuts_CustomTriggerSafe_Id_v2,
    kCuts_CustomTriggerSafe_Iso_v2,

    nElectronTriggerCutTypes
  };

}

#endif
