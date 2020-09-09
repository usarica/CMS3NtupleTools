#ifndef VERTEXOBJECT_H
#define VERTEXOBJECT_H

#include "SystematicVariations.h"
#include "ParticleObject.h"


// This is all of what is needed
#define VERTEX_VARIABLES \
VERTEX_VARIABLE(bool, is_good, false)

// These are defined for the entire event
#define VERTEX_EVENT_VARIABLES \
VERTEX_EVENT_VARIABLE(cms3_vtxs_nvtxs_t, nvtxs, 0) \
VERTEX_EVENT_VARIABLE(cms3_vtxs_nvtxs_t, nvtxs_good, 0) \
VERTEX_EVENT_VARIABLE(cms3_vtxs_nvtxs_t, nvtxs_good_JEC, 0)


class VertexVariables{
public:
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  VERTEX_VARIABLES;
#undef VERTEX_VARIABLE

  VertexVariables();
  VertexVariables(VertexVariables const& other);
  VertexVariables& operator=(const VertexVariables& other);

  void swap(VertexVariables& other);

};

class VertexObject{
public:
  VertexVariables extras;

public:
  VertexObject();
  VertexObject(const VertexObject& other);
  VertexObject& operator=(const VertexObject& other);
  ~VertexObject();

  void swap(VertexObject& other);

};

#endif
