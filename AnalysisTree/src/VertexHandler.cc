#include <cassert>
#include "VertexHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define VECTOR_ITERATOR_HANDLER_DIRECTIVES \
VERTEX_VARIABLES


const std::string VertexHandler::colName = "vtxs";

VertexHandler::VertexHandler() :
  IvyBase(),
  product_nvtxs(0),
  product_nvtxs_good(0),
  product_hasGoodPrimaryVertex(false)
{
  this->addConsumed<unsigned int>(VertexHandler::colName + "_nvtxs");
  this->addConsumed<unsigned int>(VertexHandler::colName + "_nvtxs_good");
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(VertexHandler::colName + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef VERTEX_VARIABLE
}

void VertexHandler::clear(){
  product_nvtxs = product_nvtxs_good = 0;
  product_hasGoodPrimaryVertex=false;
  
  for (ProductType_t*& prod:productList) delete prod;
  productList.clear();
}

bool VertexHandler::constructVertices(){
  clear();
  if (!currentTree) return false;

#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef VERTEX_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(VertexHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef VERTEX_VARIABLE
  allVariablesPresent &= this->getConsumedValue(VertexHandler::colName + "_nvtxs", product_nvtxs);
  allVariablesPresent &= this->getConsumedValue(VertexHandler::colName + "_nvtxs_good", product_nvtxs_good);
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "VertexHandler::constructVertices: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "VertexHandler::constructVertices: All variables are set up!" << endl;

  if (itBegin_is_good == itEnd_is_good) return true; // Construction is successful, it is just that no vertices exist.

  size_t nProducts = (itEnd_is_good - itBegin_is_good);
  productList.reserve(nProducts);
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef VERTEX_VARIABLE
  {
    size_t ip=0;
    while (it_is_good != itEnd_is_good){
      if (this->verbosity>=TVar::DEBUG) MELAout << "VertexHandler::constructVertices: Attempting vertex " << ip << "..." << endl;

      productList.push_back(new VertexObject());
      VertexObject*& obj = productList.back();

      // Set extras
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      VERTEX_VARIABLES;
#undef VERTEX_VARIABLE

      // Set the product_hasGoodPrimaryVertex flag for the first vertex
      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- The vertex is " << (*it_is_good ? "good" : "bad") << "." << endl;
      if (ip==0) product_hasGoodPrimaryVertex = (*it_is_good);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef VERTEX_VARIABLE
    }
  }

  return true;
}

void VertexHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  tree->bookBranch<unsigned int>(VertexHandler::colName + "_nvtxs", 0);
  tree->bookBranch<unsigned int>(VertexHandler::colName + "_nvtxs_good", 0);
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(VertexHandler::colName + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef VERTEX_VARIABLE
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
