#ifndef DILEPTONHANDLER_H
#define DILEPTONHANDLER_H

#include <vector>
#include "MuonObject.h"
#include "ElectronObject.h"
#include "DileptonObject.h"


class DileptonHandler{
protected:
  std::vector<DileptonObject*> productList;

  void clear(){ for (auto& prod:productList) delete prod; productList.clear(); }

  bool constructSSDileptons(
    std::vector<MuonObject*> const* muons,
    std::vector<ElectronObject*> const* electrons
  );
  bool constructOSDileptons(
    std::vector<MuonObject*> const* muons,
    std::vector<ElectronObject*> const* electrons
  );
  void setDileptonFlags() const{ for (auto* prod:productList) prod->configure(); }

public:
  DileptonHandler();
  ~DileptonHandler(){ clear(); }

  bool constructDileptons(
    std::vector<MuonObject*> const* muons,
    std::vector<ElectronObject*> const* electrons
  );
  std::vector<DileptonObject*> const& getProducts() const{ return productList; }

};


#endif
