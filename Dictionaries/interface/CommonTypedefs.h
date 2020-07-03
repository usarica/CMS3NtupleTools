#ifndef CMS3_COMMONTYPEDEFS_H
#define CMS3_COMMONTYPEDEFS_H

#include <cstdint>


typedef char cms3_charge_t;

typedef int cms3_id_t;
typedef unsigned int cms3_absid_t;

typedef short cms3_genstatus_t;
typedef unsigned short cms3_nShowerGluons_t;

typedef unsigned int cms3_muon_pogselectorbits_t;

typedef uint16_t cms3_egamma_fid_type_mask_t;

typedef unsigned char cms3_electron_charge_consistency_bits_t;

typedef unsigned short cms3_electron_cutbasedbits_t;
typedef unsigned char cms3_electron_mvacat_t;
typedef unsigned int cms3_electron_missinghits_t; // Could be unsigned char really...
typedef unsigned int cms3_electron_cutbasedbits_triggeremulation_t;

typedef unsigned int cms3_photon_cutbasedbits_t;
typedef unsigned int cms3_photon_mvacat_t;
typedef bool cms3_photon_cutbasedbits_hgg_t;

typedef char cms3_jet_genflavor_t;
typedef unsigned char cms3_jet_pujetid_t;

#endif
