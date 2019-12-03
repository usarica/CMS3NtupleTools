#ifndef DILEPTONOBJECT_H
#define DILEPTONOBJECT_H

#include "ParticleObject.h"


#define DILEPTON_VARIABLES \
DILEPTON_VARIABLE(bool, isValid, false) \
DILEPTON_VARIABLE(bool, isOS, false) \
DILEPTON_VARIABLE(bool, isSF, false) \
DILEPTON_VARIABLE(unsigned char, nTightDaughters, 0)


class DileptonVariables{
public:
#define DILEPTON_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  DILEPTON_VARIABLES;
#undef DILEPTON_VARIABLE

  DileptonVariables();
  DileptonVariables(DileptonVariables const& other);
  DileptonVariables& operator=(const DileptonVariables& other);

  void swap(DileptonVariables& other);

};

class DileptonObject : public ParticleObject{
public:
  DileptonVariables extras;

  DileptonObject();
  DileptonObject(int id_);
  DileptonObject(int id_, LorentzVector_t const& mom_);
  DileptonObject(const DileptonObject& other);
  DileptonObject& operator=(const DileptonObject& other);
  ~DileptonObject();

  void swap(DileptonObject& other);

  void configure();

  bool isValid() const{ return extras.isValid; }
  bool isOS() const{ return extras.isOS; }
  bool isSF() const{ return extras.isSF; }
  unsigned char nTightDaughters() const{ return extras.nTightDaughters; }

};

#endif
