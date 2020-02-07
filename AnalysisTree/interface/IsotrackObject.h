#ifndef ISOTRACKOBJECT_H
#define ISOTRACKOBJECT_H

#include "ParticleObject.h"


#define ISOTRACK_VARIABLES \
ISOTRACK_VARIABLE(bool, fromPV, false) \
ISOTRACK_VARIABLE(bool, is_pfCand, false) \
ISOTRACK_VARIABLE(bool, is_lostTrack, false) \
ISOTRACK_VARIABLE(bool, is_highPurityTrack, false) \
ISOTRACK_VARIABLE(bool, is_tightTrack, false) \
ISOTRACK_VARIABLE(int, nearestPFcand_id, 0) \
ISOTRACK_VARIABLE(float, nearestPFcand_deltaR, 0) \
ISOTRACK_VARIABLE(float, pfIso03_ch, 0) \
ISOTRACK_VARIABLE(float, pfIso03_comb_nofsr, 0) \
ISOTRACK_VARIABLE(float, miniIso_ch, 0) \
ISOTRACK_VARIABLE(float, miniIso_comb_nofsr, 0) \
ISOTRACK_VARIABLE(float, dxy, 0) \
ISOTRACK_VARIABLE(float, dxyerr, 0) \
ISOTRACK_VARIABLE(float, dz, 0) \
ISOTRACK_VARIABLE(float, dzerr, 0)


class IsotrackVariables{
public:
#define ISOTRACK_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  ISOTRACK_VARIABLES;
#undef ISOTRACK_VARIABLE

  IsotrackVariables();
  IsotrackVariables(IsotrackVariables const& other);
  IsotrackVariables& operator=(const IsotrackVariables& other);

  void swap(IsotrackVariables& other);

};

class IsotrackObject : public ParticleObject{
public:
  IsotrackVariables extras;

  IsotrackObject();
  IsotrackObject(int id);
  IsotrackObject(int id, LorentzVector_t const& mom_);
  IsotrackObject(const IsotrackObject& other);
  IsotrackObject& operator=(const IsotrackObject& other);
  ~IsotrackObject();

  void swap(IsotrackObject& other);

};

#endif
