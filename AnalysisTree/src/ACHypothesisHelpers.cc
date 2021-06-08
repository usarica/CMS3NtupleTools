#include <cassert>
#include "ACHypothesisHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


bool ACHypothesisHelpers::isOnshellDecay(DecayType dktype){
  return (dktype == kZZ4l_onshell);
}

TString ACHypothesisHelpers::getACHypothesisName(ACHypothesisHelpers::ACHypothesis hypo){
  switch (hypo){
  case kSM:
    return "SM";
  case kL1:
    return "L1";
  case kA2:
    return "a2";
  case kA3:
    return "a3";
  case kL1ZGs:
    return "L1ZGs";
  default:
    return "";
  };
}
TString ACHypothesisHelpers::getACHypothesisLabel(ACHypothesisHelpers::ACHypothesis hypo){
  TString res;
  switch (hypo){
  case ACHypothesisHelpers::kL1:
    res="#Lambda_{1}";
    break;
  case ACHypothesisHelpers::kA2:
    res="a_{2}";
    break;
  case ACHypothesisHelpers::kA3:
    res="a_{3}";
    break;
  case ACHypothesisHelpers::kL1ZGs:
    res="#Lambda_{1}^{Z#gamma}";
    break;
  default:
    break;
  };
  return res;
}
TString ACHypothesisHelpers::getACHypothesisFLabel(ACHypothesisHelpers::ACHypothesis hypo){
  TString res;
  switch (hypo){
  case ACHypothesisHelpers::kL1:
    res="f_{#Lambda1}";
    break;
  case ACHypothesisHelpers::kA2:
    res="f_{a2}";
    break;
  case ACHypothesisHelpers::kA3:
    res="f_{a3}";
    break;
  case ACHypothesisHelpers::kL1ZGs:
    res="f_{#Lambda1}^{Z#gamma}";
    break;
  default:
    break;
  };
  return res;
}

float ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::ACHypothesis hypo){
  if (hypo==kL1 || hypo==kL1ZGs) return 1e4;
  else return 1;
}
float ACHypothesisHelpers::getACHypothesisHZZGVal(ACHypothesisHelpers::ACHypothesis hypo){
  switch (hypo){
  case kL1:
    return -1.211020e4;
  case kA2:
    return 1.663195;
  case kA3:
    return 2.55502;
  case kL1ZGs:
    return -7613.351;
  default:
    return 1;
  };
}

std::vector<TString> ACHypothesisHelpers::getACHypothesisKDNameSet(ACHypothesisHelpers::ACHypothesis hypo, ACHypothesisHelpers::ProductionType prod_type, ACHypothesisHelpers::DecayType decay_type){
  std::vector<DiscriminantClasses::Type> KDset = ACHypothesisHelpers::getACHypothesisKDSet(hypo, prod_type, decay_type);
  vector<TString> res;
  for (auto& type:KDset) res.push_back(DiscriminantClasses::getKDName(type));
  return res;
}
std::vector<DiscriminantClasses::Type> ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::ACHypothesis hypo, ACHypothesisHelpers::ProductionType prod_type, ACHypothesisHelpers::DecayType decay_type){
  std::vector<DiscriminantClasses::Type> res;
  if (decay_type==kZZ4l_offshell){
    // First dimension, "mass", should be added separately.
    if (prod_type==kGG){
      res.push_back(DiscriminantClasses::kDbkgkin);
      switch (hypo){
      case kSM:
        res.push_back(DiscriminantClasses::kCggint);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDL1dec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDa2dec);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDa3dec);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDL1ZGsdec);
        break;
      default:
        break;
      };
    }
    else if (prod_type==kVBF){
      res.push_back(DiscriminantClasses::kDbkgjjEWQCD);
      switch (hypo){
      case kSM:
        res.push_back(DiscriminantClasses::kCjjVBFint);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDL1jjVBFdec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDa2jjVBFdec);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDa3jjVBFdec);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDL1ZGsjjVBFdec);
        break;
      default:
        break;
      };
    }
    else if (prod_type==kHadVH){
      res.push_back(DiscriminantClasses::kDbkgjjEWQCD);
      switch (hypo){
      case kSM:
        res.push_back(DiscriminantClasses::kCjjVHint);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDL1jjVHdec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDa2jjVHdec);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDa3jjVHdec);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDL1ZGsjjVHdec);
        break;
      default:
        break;
      };
    }
  }
  else if (decay_type==kZZ4l_onshell){
    if (prod_type==kGG){
      switch (hypo){
      case kSM:
        // First dimension, "mass", should be added separately.
        res.push_back(DiscriminantClasses::kDbkgkin);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDbkgm4l);
        res.push_back(DiscriminantClasses::kDL1dec);
        res.push_back(DiscriminantClasses::kDa2dec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDbkgm4l);
        res.push_back(DiscriminantClasses::kDa2dec);
        res.push_back(DiscriminantClasses::kCa2decint);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDbkgm4l);
        res.push_back(DiscriminantClasses::kDa3dec);
        res.push_back(DiscriminantClasses::kCa3decint);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDbkgm4l);
        res.push_back(DiscriminantClasses::kDL1ZGsdec);
        res.push_back(DiscriminantClasses::kDa2dec);
        //res.push_back(DiscriminantClasses::kCL1ZGsdecint);
        break;
      default:
        break;
      };
    }
    else if (prod_type==kVBF){
      switch (hypo){
      case kSM:
        // First dimension, "mass", should be added separately.
        res.push_back(DiscriminantClasses::kDbkgjjEWQCD);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDL1jjVBFdec);
        res.push_back(DiscriminantClasses::kDa2jjVBFdec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDa2jjVBFdec);
        res.push_back(DiscriminantClasses::kCa2jjVBFint);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDa3jjVBFdec);
        res.push_back(DiscriminantClasses::kCa3jjVBFint);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDL1ZGsjjVBFdec);
        res.push_back(DiscriminantClasses::kDa2jjVBFdec);
        //res.push_back(DiscriminantClasses::kCL1ZGsjjVBFint);
        break;
      default:
        break;
      };
    }
    else if (prod_type==kHadVH){
      switch (hypo){
      case kSM:
        // First dimension, "mass", should be added separately.
        res.push_back(DiscriminantClasses::kDbkgjjEWQCD);
        break;
      case kL1:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDL1jjVHdec);
        res.push_back(DiscriminantClasses::kDa2jjVHdec);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDa2jjVHdec);
        res.push_back(DiscriminantClasses::kCa2jjVHint);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDa3jjVHdec);
        res.push_back(DiscriminantClasses::kCa3jjVHint);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDbkgm4ljjEWQCD);
        res.push_back(DiscriminantClasses::kDL1ZGsjjVHdec);
        res.push_back(DiscriminantClasses::kDa2jjVHdec);
        //res.push_back(DiscriminantClasses::kCL1ZGsjjVHint);
        break;
      default:
        break;
      };
    }
  }
  else if (decay_type==kZZ2l2nu_offshell){
    // First dimension, "mass", should be added separately.
    if (prod_type==kVBF){
      res.push_back(DiscriminantClasses::kDjjVBF); // Always used
      switch (hypo){
      case kL1:
        res.push_back(DiscriminantClasses::kDjjVBFL1);
        break;
      case kSM:
      case kA2:
        res.push_back(DiscriminantClasses::kDjjVBFa2);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDjjVBFa3);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDjjVBFL1ZGs);
        break;
      default:
        break;
      };
    }
    /*
    else if (prod_type==kHadVH){
      // Always used
      res.push_back(DiscriminantClasses::kDjjZH);
      res.push_back(DiscriminantClasses::kDjjWH);
      switch (hypo){
      case kL1:
        res.push_back(DiscriminantClasses::kDjjZHL1);
        res.push_back(DiscriminantClasses::kDjjWHL1);
        break;
      case kA2:
        res.push_back(DiscriminantClasses::kDjjZHa2);
        res.push_back(DiscriminantClasses::kDjjWHa2);
        break;
      case kA3:
        res.push_back(DiscriminantClasses::kDjjZHa3);
        res.push_back(DiscriminantClasses::kDjjWHa3);
        break;
      case kL1ZGs:
        res.push_back(DiscriminantClasses::kDjjZHL1ZGs);
        break;
      default:
        break;
      };
    }
    */
  }
  else if (decay_type==kZW3l1nu){
    // First dimension, "mass", should be added separately.
    // There are no discriminants involved in the ZW analysis.
  }

  return res;
}

TString ACHypothesisHelpers::getDecayFinalStateLabel(ACHypothesisHelpers::DecayType dktype){
  switch (dktype){
  case kZZ4l_onshell:
  case kZZ4l_offshell:
    return "4l";
  case kZZ2l2nu_offshell:
    return "2l2#nu";
  case kZW3l1nu:
    return "3l1#nu";
  default:
    MELAerr << "ACHypothesisHelpers::getDecayFinalStateLabel: Decay type " << dktype << " is not recognized. Please fix the implementation." << endl;
    assert(0);
    return "";
  }
}
