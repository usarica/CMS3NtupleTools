#include <algorithm>
#include <utility>
#include <cmath>
#include "VertexObject.h"


VertexVariables::VertexVariables(){
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  VERTEX_VARIABLES;
#undef VERTEX_VARIABLE
}
VertexVariables::VertexVariables(VertexVariables const& other){
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  VERTEX_VARIABLES;
#undef VERTEX_VARIABLE
}
void VertexVariables::swap(VertexVariables& other){
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  VERTEX_VARIABLES;
#undef VERTEX_VARIABLE
}
VertexVariables& VertexVariables::operator=(const VertexVariables& other){
  VertexVariables tmp(other);
  swap(tmp);
  return *this;
}

VertexObject::VertexObject() :
  extras()
{}
VertexObject::VertexObject(const VertexObject& other) :
  extras(other.extras)
{}
void VertexObject::swap(VertexObject& other){
  extras.swap(other.extras);
}
VertexObject& VertexObject::operator=(const VertexObject& other){
  VertexObject tmp(other);
  swap(tmp);
  return *this;
}
VertexObject::~VertexObject(){}
