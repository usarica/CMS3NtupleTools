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
  unsigned int product_nvtxs;
  unsigned int product_nvtxs_good;
  bool product_hasGoodPrimaryVertex;
  std::vector<ProductType_t*> productList;

  void clear(){ product_nvtxs=product_nvtxs_good=0; product_hasGoodPrimaryVertex=false; for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

public:
  // Constructors
  VertexHandler();

  // Destructors
  ~VertexHandler(){ clear(); }

  bool constructVertices();
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }
  unsigned int const& getNVertices() const{ return product_nvtxs; }
  unsigned int const& getNGoodVertices() const{ return product_nvtxs_good; }
  bool hasGoodVertex() const{ return product_nvtxs_good>0; }
  bool const& hasGoodPrimaryVertex() const{ return product_hasGoodPrimaryVertex; }

  static void bookBranches(BaseTree* tree);

};


#endif
