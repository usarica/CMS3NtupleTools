#include <cassert>
#include "SystematicVariations.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


std::string SystematicsHelpers::getSystCoreName(SystematicsHelpers::SystematicVariationTypes const& type){
  switch (type){
  case sNominal:
    return "Nominal";

  case tPDFScaleDn:
  case tPDFScaleUp:
    return "FacScale";
  case tQCDScaleDn:
  case tQCDScaleUp:
    return "RenScale";
  case tAsMZDn:
  case tAsMZUp:
    return "AsMZ";
  case tPDFReplicaDn:
  case tPDFReplicaUp:
    return "PDFVariation";
  case tPythiaScaleDn:
  case tPythiaScaleUp:
    return "PythiaScale";
  case tPythiaTuneDn:
  case tPythiaTuneUp:
    return "PythiaTune";

  case eEleEffDn:
  case eEleEffUp:
    return "ElectronEff";
  case eEleEffStatDn:
  case eEleEffStatUp:
    return "ElectronEff_Stat";
  case eEleEffSystDn:
  case eEleEffSystUp:
    return "ElectronEff_Syst";
  case eEleEffAltMCDn:
  case eEleEffAltMCUp:
    return "ElectronEff_AltMC";
  case eEleScaleDn:
  case eEleScaleUp:
    return "ElectronScale";
  case eEleResDn:
  case eEleResUp:
    return "ElectronRes";

  case eMuEffDn:
  case eMuEffUp:
    return "MuonEff";
  case eMuEffStatDn:
  case eMuEffStatUp:
    return "MuonEff_Stat";
  case eMuEffSystDn:
  case eMuEffSystUp:
    return "MuonEff_Syst";
  case eMuEffAltMCDn:
  case eMuEffAltMCUp:
    return "MuonEff_AltMC";
  case eMuScaleDn:
  case eMuScaleUp:
    return "MuonScale";
  case eMuResDn:
  case eMuResUp:
    return "MuonRes";

  case ePhoEffDn:
  case ePhoEffUp:
    return "PhotonEff";
  case ePhoScaleDn:
  case ePhoScaleUp:
    return "PhotonScale";
  case ePhoResDn:
  case ePhoResUp:
    return "PhotonRes";

  case eMETDn:
  case eMETUp:
    return "MET";
  case eJECDn:
  case eJECUp:
    return "JEC";
  case eJERDn:
  case eJERUp:
    return "JER";
  case ePUDn:
  case ePUUp:
    return "PU";
  case eBTagSFDn:
  case eBTagSFUp:
    return "BTagSF";

  case eL1PrefiringDn:
  case eL1PrefiringUp:
    return "L1Prefiring";

  case sUncorrected:
    return "Uncorrected";

  default:
    MELAerr << "SystematicsHelpers::getSystCoreName: Systematic " << type << " is not defined!" << endl;
    assert(0);
    return "";
  }
}
bool SystematicsHelpers::isDownSystematic(SystematicsHelpers::SystematicVariationTypes const& type){
  if (type<nSystematicVariations && type!=sNominal) return (((int) type)%2 == 1);
  else return false;
}
bool SystematicsHelpers::isUpSystematic(SystematicsHelpers::SystematicVariationTypes const& type){
  if (type<nSystematicVariations && type!=sNominal) return (((int) type)%2 == 0);
  else return false;
}
SystematicsHelpers::SystematicVariationTypes SystematicsHelpers::getSystComplement(SystematicsHelpers::SystematicVariationTypes const& type){
  bool isDn = SystematicsHelpers::isDownSystematic(type);
  bool isUp = SystematicsHelpers::isUpSystematic(type);
  if (isDn == isUp) return type;
  else if (isDn) return static_cast<SystematicsHelpers::SystematicVariationTypes>(static_cast<int>(type)+1);
  else /*if (isUp)*/ return static_cast<SystematicsHelpers::SystematicVariationTypes>(static_cast<int>(type)-1);
}
std::string SystematicsHelpers::getSystName(SystematicsHelpers::SystematicVariationTypes const& type){
  std::string res = getSystCoreName(type);
  if (SystematicsHelpers::isDownSystematic(type)) res = res + "Dn";
  else if (SystematicsHelpers::isUpSystematic(type)) res = res + "Up";
  return res;
}
