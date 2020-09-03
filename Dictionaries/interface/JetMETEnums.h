#ifndef JETMETENUMS_H
#define JETMETENUMS_H

namespace JetMETEnums{
  enum METShiftType{
    kMETShift_JECNominal,
    kMETShift_JECDn,
    kMETShift_JECUp,
    kMETShift_JECNominal_JERNominal,
    kMETShift_JECNominal_JERDn,
    kMETShift_JECNominal_JERUp,
    kMETShift_JECDn_JERNominal,
    kMETShift_JECUp_JERNominal,

    nMETShiftTypes
  };

}

#endif
