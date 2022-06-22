#ifndef CMS3_COMMONTYPEDEFS_H
#define CMS3_COMMONTYPEDEFS_H

#include <cstdint>
#include "DataFormats/Provenance/interface/RunLumiEventNumber.h"
#include "IvyFramework/IvyDataTools/interface/IvyDataTypes.h"


typedef ivy_listIndex_short_t cms3_listIndex_short_t;
typedef ivy_listIndex_long_t cms3_listIndex_long_t;
typedef ivy_listSize_t cms3_listSize_t;

typedef ivy_listIndex_signed_short_t cms3_listIndex_signed_short_t;
typedef ivy_listIndex_signed_long_t cms3_listIndex_signed_long_t;

typedef cms3_listSize_t cms3_vtxs_nvtxs_t;

typedef short cms3_triggertype_t;
typedef cms3_listIndex_short_t cms3_triggerIndex_t;

typedef char cms3_charge_t;

typedef ivy_id_t cms3_id_t;
typedef ivy_absid_t cms3_absid_t;

typedef short cms3_genstatus_t;
typedef unsigned short cms3_nShowerGluons_t;

typedef short cms3_refkey_t;

typedef unsigned int cms3_muon_pogselectorbits_t;
typedef unsigned char cms3_muon_cutbasedbits_h4l_t;

typedef uint16_t cms3_egamma_fid_type_mask_t;

typedef unsigned char cms3_electron_charge_consistency_bits_t;

typedef unsigned short cms3_electron_cutbasedbits_t;
typedef unsigned char cms3_electron_mvacat_t;
typedef unsigned short cms3_electron_missinghits_t; // Could be unsigned char really...
typedef unsigned char cms3_electron_cutbasedbits_egPFElectron_t;
typedef unsigned int cms3_electron_cutbasedbits_triggeremulation_t;

typedef unsigned short cms3_photon_cutbasedbits_t; // Could be unsigned char, but better to keep as short-type in case new bits are added.
typedef unsigned char cms3_photon_mvacat_t;
typedef bool cms3_photon_cutbasedbits_hgg_t;
typedef unsigned char cms3_photon_cutbasedbits_egPFPhoton_t;

typedef char cms3_jet_genflavor_t;
typedef unsigned char cms3_jet_pujetid_t;
typedef uint8_t cms3_metsafety_t;

typedef uint16_t cms3_pfcand_qualityflag_t;

typedef edm::RunNumber_t cms3_runnumber_t;
typedef edm::LuminosityBlockNumber_t cms3_lumiblock_t;
typedef edm::EventNumber_t cms3_eventnumber_t;


#endif
