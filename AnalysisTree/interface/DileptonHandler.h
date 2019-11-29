#ifndef DILEPTONHANDLER_H
#define DILEPTONHANDLER_H

#include <vector>
#include "MuonObject.h"
#include "ElectronObject.h"


class DileptonHandler{
protected:
  std::vector<ParticleObject*> productList;

  void clear(){ for (auto& prod:productList) delete prod; productList.clear(); }

  bool constructSSDileptons(
    std::vector<MuonObject*> const* muons,
    std::vector<ElectronObject*> const* electrons
  );
  bool constructOSDileptons(
    std::vector<MuonObject*> const* muons,
    std::vector<ElectronObject*> const* electrons
  );

public:
  DileptonHandler();
  ~DileptonHandler(){ clear(); }

  bool constructDileptons(
    std::vector<MuonObject*> const* muons,
    std::vector<ElectronObject*> const* electrons
  );
  std::vector<ParticleObject*> const& getProducts() const{ return productList; }

};


#endif
