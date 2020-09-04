#ifndef VERTEXHANDLER_H
#define VERTEXHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "VertexObject.h"


class VertexHandler : public IvyBase{
public:
  typedef VertexObject ProductType_t;
  static const std::string colName;

protected:
  // These two variables are not the same as the size of the vertex collection stored
#define VERTEX_EVENT_VARIABLE(TYPE, NAME, DEFVAL) TYPE product_##NAME;
  VERTEX_EVENT_VARIABLES;
#undef VERTEX_EVENT_VARIABLE
  bool product_hasGoodPrimaryVertex;

  std::vector<ProductType_t*> productList;

  void clear();

public:
  // Constructors
  VertexHandler();

  // Destructors
  ~VertexHandler(){ clear(); }

  bool constructVertices();
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }
  cms3_listSize_t const& getNVertices() const{ return product_n_vtxs; }
  cms3_listSize_t const& getNGoodVertices() const{ return product_n_vtxs_good; }
  cms3_listSize_t const& getNGoodJECVertices() const{ return product_n_vtxs_good_JEC; }
  bool hasVertex() const{ return product_n_vtxs>0; }
  bool hasGoodVertex() const{ return product_n_vtxs_good>0; }
  bool hasGoodJECVertex() const{ return product_n_vtxs_good_JEC>0; }
  bool const& hasGoodPrimaryVertex() const{ return product_hasGoodPrimaryVertex; }

  static void bookBranches(BaseTree* tree);

};


#endif
