#ifndef OVERLAPMAPHANDLER_H
#define OVERLAPMAPHANDLER_H

#include <vector>

#include "IvyBase.h"
#include "OverlapMapElement.h"


template<typename T, typename U> class OverlapMapHandler : public IvyBase{
public:
  typedef OverlapMapElement<T, U> ProductType_t;
  static const std::string colName;

protected:
  std::vector<ProductType_t*> productList;

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

public:
  // Constructors
  OverlapMapHandler();

  // Destructors
  ~OverlapMapHandler(){ clear(); }

  bool constructOverlapMaps();

  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  void bookBranches(BaseTree* tree);

};


#define OVERLAPMAP_SPECIALIZATIONS \
OVERLAPMAP_SPECIALIZATION(MuonObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(MuonObject, AK8JetObject) \
OVERLAPMAP_SPECIALIZATION(ElectronObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(ElectronObject, AK8JetObject) \
OVERLAPMAP_SPECIALIZATION(PhotonObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(PhotonObject, AK8JetObject)

#define OVERLAPMAP_SPECIALIZATION(T1, T2) \
template<> OverlapMapHandler<T1, T2>::OverlapMapHandler(); \
template<> bool OverlapMapHandler<T1, T2>::constructOverlapMaps(); \
template<> void OverlapMapHandler<T1, T2>::bookBranches(BaseTree* tree);

OVERLAPMAP_SPECIALIZATIONS;

#undef OVERLAPMAP_SPECIALIZATION
#undef OVERLAPMAP_SPECIALIZATIONS


#endif
