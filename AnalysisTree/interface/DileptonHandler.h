#ifndef DILEPTONHANDLER_H
#define DILEPTONHANDLER_H

#include <vector>
#include "VerbosityLevel.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "DileptonObject.h"


class DileptonHandler{
protected:
  MiscUtils::VerbosityLevel verbosity;

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

  void setVerbosity(MiscUtils::VerbosityLevel verbosity_){ verbosity=verbosity_; }

  bool constructDileptons(
    std::vector<MuonObject*> const* muons,
    std::vector<ElectronObject*> const* electrons
  );
  std::vector<DileptonObject*> const& getProducts() const{ return productList; }

};


#endif
