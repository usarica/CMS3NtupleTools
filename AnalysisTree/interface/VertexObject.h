#ifndef VERTEXOBJECT_H
#define VERTEXOBJECT_H

#include "SystematicVariations.h"
#include "ParticleObject.h"


// This is all of what is needed
#define VERTEX_VARIABLES \
VERTEX_VARIABLE(bool, is_good, false)


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
