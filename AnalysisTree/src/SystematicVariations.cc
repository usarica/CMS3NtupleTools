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
  case eEleScaleDn:
  case eEleScaleUp:
    return "ElectronScale";
  case eEleResDn:
  case eEleResUp:
    return "ElectronRes";

  case eMuEffDn:
  case eMuEffUp:
    return "MuonEff";
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

  case sUncorrected:
    return "Uncorrected";

  default:
    MELAerr << "SystematicsHelpers::getSystCoreName: Systematic " << type << " is not defined!" << endl;
    assert(0);
    return "";
  }
}
std::string SystematicsHelpers::getSystName(SystematicsHelpers::SystematicVariationTypes const& type){
  std::string res = getSystCoreName(type);
  switch (type){
  case tPDFScaleDn:
  case tQCDScaleDn:
  case tAsMZDn:
  case tPDFReplicaDn:
  case tPythiaScaleDn:
  case tPythiaTuneDn:
  case eEleEffDn:
  case eEleScaleDn:
  case eEleResDn:
  case eMuEffDn:
  case eMuScaleDn:
  case eMuResDn:
  case ePhoEffDn:
  case ePhoScaleDn:
  case ePhoResDn:
  case eMETDn:
  case eJECDn:
  case eJERDn:
  case ePUDn:
  case eBTagSFDn:
    res = res + "Dn";
    break;

  case tPDFScaleUp:
  case tQCDScaleUp:
  case tAsMZUp:
  case tPDFReplicaUp:
  case tPythiaScaleUp:
  case tPythiaTuneUp:
  case eEleEffUp:
  case eEleScaleUp:
  case eEleResUp:
  case eMuEffUp:
  case eMuScaleUp:
  case eMuResUp:
  case ePhoEffUp:
  case ePhoScaleUp:
  case ePhoResUp:
  case eMETUp:
  case eJECUp:
  case eJERUp:
  case ePUUp:
  case eBTagSFUp:
    res = res + "Up";

  case sNominal:
  case sUncorrected:
    break;

  default:
    MELAerr << "SystematicsHelpers::getSystName: Systematic " << type << " is not defined!" << endl;
    assert(0);
  }
  return res;
}
