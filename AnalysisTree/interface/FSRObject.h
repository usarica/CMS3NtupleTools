#ifndef FSROBJECT_H
#define FSROBJECT_H

#include "ParticleObject.h"
#include "PhotonObject.h"


#define FSR_VECTOR_VARIABLES \
FSR_VECTOR_VARIABLE(std::vector<cms3_listIndex_short_t>, fsrMatch_muon_index_list) \
FSR_VECTOR_VARIABLE(std::vector<cms3_listIndex_short_t>, fsrMatch_electron_index_list) \
FSR_VECTOR_VARIABLE(std::vector<cms3_listIndex_short_t>, fsrMatch_ak4jet_index_list) \
FSR_VECTOR_VARIABLE(std::vector<cms3_listIndex_short_t>, fsrMatch_ak8jet_index_list) \
FSR_VECTOR_VARIABLE(std::vector<cms3_listIndex_short_t>, photonVeto_index_list)

#define FSR_VARIABLES \
FSR_VECTOR_VARIABLES


class FSRVariables{
public:
#define FSR_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
#define FSR_VECTOR_VARIABLE(TYPE, NAME) TYPE NAME;
  FSR_VARIABLES;
#undef FSR_VECTOR_VARIABLE
#undef FSR_VARIABLE

  FSRVariables();
  FSRVariables(FSRVariables const& other);
  FSRVariables& operator=(const FSRVariables& other);

  void swap(FSRVariables& other);

};

class FSRObject : public ParticleObject{
public:
  FSRVariables extras;

  FSRObject();
  FSRObject(LorentzVector_t const& mom_);
  FSRObject(const FSRObject& other);
  FSRObject& operator=(const FSRObject& other);
  ~FSRObject();

  void swap(FSRObject& other);

  void setAssociatedPhoton(PhotonObject* part){ this->addMother(part); }

};

#endif
