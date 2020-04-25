#include <cassert>
#include <sstream>
#include <FWCore/Utilities/interface/Exception.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include <CMS3/NtupleMaker/interface/CMS3ObjectHelpers.h>
#include <CMS3/NtupleMaker/interface/KFactorHelpers.h>
#include <CMSDataTools/AnalysisTree/interface/HostHelpersCore.h>
#include "MELAStreamHelpers.hh"
#include "TDirectory.h"


using namespace std;
using namespace MELAStreamHelpers;


namespace KFactorHelpers{

  float getBeamEnergy(int const& year){
    float beamEnergy=0;
    switch (year){
    case 2011:
      beamEnergy=3500;
      break;
    case 2012:
      beamEnergy=4000;
      break;
    case 2015:
    case 2016:
    case 2017:
    case 2018:
      beamEnergy=6500;
      break;

    default:
      throw cms::Exception("UnknownYear") << "KFactorHelpers::getBeamEnergy: Beam energy for year " << year << " is not implemented.";
      break;
    }
    return beamEnergy;
  }

  void getVVTopology(
    VVFinalStateType const& type,
    std::vector<reco::GenParticle> const& genparticles,
    std::vector<reco::GenParticle const*>& incomingQuarks,
    std::vector<reco::GenParticle const*>& incomingGluons,
    std::vector<reco::GenParticle const*>& outgoingQuarks,
    std::vector<reco::GenParticle const*>& outgoingGluons,
    std::pair<reco::GenParticle const*, reco::GenParticle const*>& V1pair,
    std::pair<reco::GenParticle const*, reco::GenParticle const*>& V2pair
  ){
    static const float mZsq = PDGHelpers::Zmass * PDGHelpers::Zmass;
    static const float mWsq = PDGHelpers::Wmass * PDGHelpers::Wmass;
    static const float mZGaZ = PDGHelpers::Zmass * PDGHelpers::Zwidth;
    static const float mWGaW = PDGHelpers::Wmass * PDGHelpers::Wwidth;
    static const float mZGaZsq = mZGaZ*mZGaZ;
    static const float mWGaWsq = mWGaW*mWGaW;

    bool const doZZ = (type == kZZ);
    bool const doWZ = (type == kWZ);
    bool const doWW = (type == kWW);
    if (!doZZ && !doWZ && !doWW){
      throw cms::Exception("InvalidTopology") << "KFactorHelpers::getVVTopology: Invalid VV topology flag " << type;
      return;
    }

    incomingQuarks.clear(); incomingQuarks.reserve(2);
    incomingGluons.clear(); incomingGluons.reserve(2);
    outgoingQuarks.clear(); outgoingQuarks.reserve(2);
    outgoingGluons.clear(); outgoingGluons.reserve(2);
    std::vector<reco::GenParticle const*> lepMinusPlus[3][2]; // l-, l+
    std::vector<reco::GenParticle const*> lepNuNubar[3][2]; // nu, nub
    std::vector<reco::GenParticle const*> quarkAntiquark[6][2]; // q, qb
    for (auto it_part = genparticles.cbegin(); it_part != genparticles.cend(); it_part++){
      reco::GenParticle const* part = &(*it_part);
      if (!part->isHardProcess()) continue;
      int id = part->pdgId();
      int status = part->status();
      //MELAout << "Gen particle (id, status) = (" << id << ", " << status << ")" << endl;
      if (PDGHelpers::isAQuark(id) && status == 21) incomingQuarks.push_back(part);
      else if (PDGHelpers::isAGluon(id) && status == 21) incomingGluons.push_back(part);
      else if (PDGHelpers::isAQuark(id)) quarkAntiquark[std::abs(id)-1][(id>0 ? 0 : 1)].push_back(part);
      else if (PDGHelpers::isALepton(id)) lepMinusPlus[(std::abs(id)-11)/2][(id>0 ? 0 : 1)].push_back(part);
      else if (PDGHelpers::isANeutrino(id)) lepNuNubar[(std::abs(id)-12)/2][(id>0 ? 0 : 1)].push_back(part);
      else if (PDGHelpers::isAGluon(id)) outgoingGluons.push_back(part);
    }
    CMS3ObjectHelpers::sortByGreaterPz(incomingQuarks);
    CMS3ObjectHelpers::sortByGreaterPz(incomingGluons);

    /*
    for (auto const& part:incomingQuarks) MELAout << "Incoming quark (id, status) = ( " << part->pdgId() << " , " << part->status() << " )" << endl;
    for (auto const& part:incomingGluons) MELAout << "Incoming gluon (id, status) = ( " << part->pdgId() << " , " << part->status() << " )" << endl;
    for (unsigned short c=0; c<3; c++){
      for (unsigned short d=0; d<2; d++){
        for (auto const& part:lepMinusPlus[c][d]) MELAout << "Lepton[" << c << "][" << d << "] (id, status) = ( " << part->pdgId() << " , " << part->status() << " )" << endl;
      }
    }
    for (unsigned short c=0; c<3; c++){
      for (unsigned short d=0; d<2; d++){
        for (auto const& part:lepNuNubar[c][d]) MELAout << "Neutrino[" << c << "][" << d << "] (id, status) = ( " << part->pdgId() << " , " << part->status() << " )" << endl;
      }
    }
    for (unsigned short c=0; c<6; c++){
      for (unsigned short d=0; d<2; d++){
        for (auto const& part:quarkAntiquark[c][d]) MELAout << "Quark[" << c << "][" << d << "] (id, status) = ( " << part->pdgId() << " , " << part->status() << " )" << endl;
      }
    }
    */

    // First, construct all possible Vs
    std::vector< std::pair<reco::GenParticle const*, reco::GenParticle const*> > tmpVhandle;
    if (doZZ || doWZ){ // Find Zs
      for (unsigned char c=0; c<3; c++){ for (auto const& F1:lepMinusPlus[c][0]){ for (auto const& F2:lepMinusPlus[c][1]) tmpVhandle.emplace_back(F1, F2); } }
      for (unsigned char c=0; c<3; c++){ for (auto const& F1:lepNuNubar[c][0]){ for (auto const& F2:lepNuNubar[c][1]) tmpVhandle.emplace_back(F1, F2); } }
      for (unsigned char c=0; c<6; c++){ for (auto const& F1:quarkAntiquark[c][0]){ for (auto const& F2:quarkAntiquark[c][1]) tmpVhandle.emplace_back(F1, F2); } }
    }
    if (doWW || doWZ){ // Find Ws
      for (unsigned char c=0; c<3; c++){
        for (unsigned char signW=0; signW<2; signW++){ // ==0: W+, ==1: W-
          for (auto const& F1:lepMinusPlus[c][1-signW]){
            for (auto const& F2:lepNuNubar[c][signW]) tmpVhandle.emplace_back(F1, F2);
          }
        }
      }
      for (unsigned char c=0; c<6; c++){
        for (unsigned char d=0; d<6; d++){
          if (d%2 == c%2) continue;
          for (auto const& F1:quarkAntiquark[c][0]){
            for (auto const& F2:quarkAntiquark[d][1]) tmpVhandle.emplace_back(F1, F2);
          }
        }
      }
    }
    for (auto& tmppair:tmpVhandle){ if (tmppair.first->pdgId()<0) std::swap(tmppair.first, tmppair.second); }
    //for (auto const& tmppair:tmpVhandle) MELAout << "Intermediate V (m=" << (tmppair.first->p4() + tmppair.second->p4()).M() << ") found. (id1, st1) = " << tmppair.first->pdgId() << ", " << tmppair.first->status() << ", (id2, st2) = " << tmppair.second->pdgId() << ", " << tmppair.second->status() << endl;

    // Determine best V1 and V2
    std::pair<reco::GenParticle const*, reco::GenParticle const*> const* bestV1pair=nullptr;
    std::pair<reco::GenParticle const*, reco::GenParticle const*> const* bestV2pair=nullptr;
    float bestDiffMass=-1;
    for (auto it1 = tmpVhandle.cbegin(); it1!=tmpVhandle.cend(); it1++){
      auto const& Vi = *it1;
      for (auto it2 = it1; it2!=tmpVhandle.cend(); it2++){
        if (it2==it1) continue;
        auto const& Vj = *it2;
        if (Vi.first==Vj.first || Vi.first==Vj.second || Vi.second==Vj.first || Vi.second==Vj.second) continue;
        bool isViZ = Vi.first->pdgId() == -Vi.second->pdgId();
        bool isVjZ = Vj.first->pdgId() == -Vj.second->pdgId();
        if (doZZ && !(isViZ && isVjZ)) continue;
        if (doWZ && (isViZ == isVjZ)) continue;
        if (doWW){
          if (isViZ || isVjZ) continue;
          // Additional check to ensure that the incoming single quark flavor is conserved
          if (incomingQuarks.size()==1){
            int qflav = incomingQuarks.front()->pdgId();
            unsigned int n_outgoing_sameflavors=0;
            for (auto const& part:quarkAntiquark[std::abs(qflav)-1][qflav>0 ? 0 : 1]){
              if (part==Vi.first || part==Vi.second || part==Vj.first || part==Vj.second) continue;
              n_outgoing_sameflavors++;
            }
            if (n_outgoing_sameflavors%2 == 0) continue; // Make sure there is an odd number of them
          }
        }
        auto pVi = Vi.first->p4() + Vi.second->p4();
        auto pVj = Vj.first->p4() + Vj.second->p4();
        float diffMass = std::pow((std::pow(pVi.M(), 2) - (isViZ ? mZsq : mWsq)), 2) + (isViZ ? mZGaZsq : mWGaWsq);
        diffMass *= std::pow((std::pow(pVj.M(), 2) - (isVjZ ? mZsq : mWsq)), 2) + (isVjZ ? mZGaZsq : mWGaWsq);
        if (bestDiffMass<0.f || diffMass<bestDiffMass){
          bestDiffMass = diffMass;
          bestV1pair = &Vi;
          bestV2pair = &Vj;
        }
      }
    }

    if (!bestV1pair) throw cms::Exception("InvalidTopology") << "KFactorHelpers::getVVTopology: Best V1 cannot be determined!";
    if (!bestV2pair) throw cms::Exception("InvalidTopology") << "KFactorHelpers::getVVTopology: Best V2 cannot be determined!";

    // Return best-Z, second-best-Z in ZZ, ZW in WZ, and W+W- in WW
    bool doSwap = false;
    auto pBestV1 = bestV1pair->first->p4() + bestV1pair->second->p4();
    auto pBestV2 = bestV2pair->first->p4() + bestV2pair->second->p4();
    if (doWW) doSwap = (PDGHelpers::isANeutrino(bestV2pair->first->pdgId()) || PDGHelpers::isUpTypeQuark(bestV2pair->first->pdgId()));
    else if (doWZ) doSwap = (bestV2pair->first->pdgId() == -bestV2pair->second->pdgId());
    else doSwap = (std::abs(pBestV1.M() - PDGHelpers::Zmass)>std::abs(pBestV2.M() - PDGHelpers::Zmass));
    if (doSwap) std::swap(bestV1pair, bestV2pair);

    //MELAout << "m12 = " << (bestV1pair->first->p4() + bestV1pair->second->p4()).M() << endl;
    //MELAout << "m34 = " << (bestV2pair->first->p4() + bestV2pair->second->p4()).M() << endl;

    V1pair = *bestV1pair;
    V2pair = *bestV2pair;
    
    // Finally, determine the associated outgoing quarks
    for (unsigned char c=0; c<6; c++){
      for (unsigned char d=0; d<1; d++){
        for (auto const& part:quarkAntiquark[c][d]){
          if (part==bestV1pair->first || part==bestV1pair->second || part==bestV2pair->first || part==bestV2pair->second) continue;
          outgoingQuarks.push_back(part);
        }
      }
    }
  }

  double xsecRatio_qqZZ4f_QCD_NLO_NNLO_13TeV_byMass(double mass, bool is_4f_2f2fp, unsigned char order){ // is_4f_2f2fp: 0=4l, 1=2l2l'; order: 0=NLO/LO, 1=NNLO/NLO
    // Xsecs come from CJLST/ZZAnalysis:AnalysisStep/src/kFactors.C::xsec_qqZZ_qcd_M
    static const std::vector< std::vector< std::vector<double> > > xseclist={
      {
        { 0, 1.095828532, 1.623089848, 2.006355102 },
        { 25, 1.70498411, 2.55680607, 3.00553358 },
        { 50, 0.7214786, 1.07670322, 1.26022261 },
        { 75, 9.93585169, 12.49651341, 12.8890551 },
        { 100, 1.88979055, 2.32301901, 2.44580392 },
        { 125, 0.80689284, 1.06246962, 1.18239288 },
        { 150, 0.420080936, 0.576994967, 0.65408981 },
        { 175, 0.573751417, 0.755001119, 0.83318604 },
        { 200, 0.6426366, 0.84359613, 0.92841113 },
        { 225, 0.460870086, 0.6155138, 0.68305214 },
        { 250, 0.325131608, 0.440635023, 0.493815795 },
        { 275, 0.234474317, 0.321194629, 0.358418201 },
        { 300, 0.173099357, 0.239157872, 0.272415573 },
        { 325, 0.130327192, 0.18182296, 0.208832023 },
        { 350, 0.100194403, 0.140405709, 0.160927729 },
        { 375, 0.078065965, 0.109969337, 0.125995342 },
        { 400, 0.061860538, 0.087368406, 0.099450959 },
        { 425, 0.049470606, 0.070200812, 0.081096816 },
        { 450, 0.040036726, 0.057035286, 0.064837612 },
        { 475, 0.03272652, 0.046808668, 0.053386015 }
      },
      {
        { 0, 5.80657303, 8.38413944, 10.48809451 },
        { 25, 4.62863121, 6.84873304, 8.38692118 },
        { 50, 1.61514729, 2.39428421, 2.85607864 },
        { 75, 19.10610496, 24.00530313, 25.10894849 },
        { 100, 3.86816714, 4.7520914, 5.14762763 },
        { 125, 1.64560534, 2.16841494, 2.38514732 },
        { 150, 0.85582638, 1.17224019, 1.3679862 },
        { 175, 1.1629688, 1.52672919, 1.68549457 },
        { 200, 1.29516749, 1.69836336, 1.87826529 },
        { 225, 0.92621932, 1.23581661, 1.36793012 },
        { 250, 0.65489237, 0.88543635, 0.98456032 },
        { 275, 0.470967807, 0.64402568, 0.73111461 },
        { 300, 0.347254403, 0.480862926, 0.538062005 },
        { 325, 0.261981026, 0.363080611, 0.413543401 },
        { 350, 0.201165289, 0.280556435, 0.323941619 },
        { 375, 0.157331951, 0.220998012, 0.259322746 },
        { 400, 0.124083729, 0.174598837, 0.209681592 },
        { 425, 0.099317367, 0.139860047, 0.166315351 },
        { 450, 0.080436156, 0.114091349, 0.135250739 },
        { 475, 0.065390164, 0.093794404, 0.104784251 }
      }
    };

    if (mass<0.) mass=0;
    std::vector< std::vector<double> > const& xsec = xseclist.at(is_4f_2f2fp);
    auto it_xsec = xsec.cbegin();
    auto itEnd_xsec = xsec.cend();
    for (; it_xsec != itEnd_xsec; it_xsec++){
      auto itNext_xsec = it_xsec; itNext_xsec++;
      if (mass>=it_xsec->front() && (itNext_xsec==itEnd_xsec || mass<itNext_xsec->front())) break;
    }
    if (static_cast<size_t>(order+2)>=it_xsec->size()) throw cms::Exception("UnknownOrder") << "KFactorHelpers::xsecRatio_qqZZ4f_QCD_NLO_NNLO_13TeV_byMass: Order " << order << " is not implemented.";
    return it_xsec->at(order+2) / it_xsec->at(order+1);
  }

  double xsecRatio_qqWW4f_QCD_NLO_NNLO_13TeV_byPt(double pt, unsigned char order){ // order: 0=NLO/LO, 1=NNLO/NLO
    // Xsecs come from latinos/LatinoAnalysis commit 9864c75c7d532b57fa399d0631634e98e87108d4. Later commit is totally useless because it tells only about NNLL corrections on NNLO, and nothing about the actual NNLO/NLO or NLO/LO ratios.
    // NLO = Gardener/python/data/wwresum/powheg_2l2nu_nlo.dat (multiplied by 1000 to match nnlo_central.dat)
    // NNLO = Gardener/python/data/wwresum/nnlo_central.dat
    static const std::vector< std::pair<double, double> > pt_dxsecNLOdpt_pairs={
      { 0, 528.447 },
      { 1, 1813.59 },
      { 2, 2783.3 },
      { 3, 3336.31 },
      { 4, 3537.45 },
      { 5, 3544.56 },
      { 6, 3439.11 },
      { 7, 3271.47 },
      { 8, 3108.5 },
      { 9, 2943.3 },
      { 10, 2792.09 },
      { 11, 2638.71 },
      { 12, 2505.28 },
      { 13, 2360.85 },
      { 14, 2263.21 },
      { 15, 2125.06 },
      { 16, 2028.94 },
      { 17, 1912.81 },
      { 18, 1828.56 },
      { 19, 1742.14 },
      { 20, 1667.8 },
      { 21, 1585.01 },
      { 22, 1503.15 },
      { 23, 1454.35 },
      { 24, 1367.93 },
      { 25, 1328.79 },
      { 26, 1269.86 },
      { 27, 1214.83 },
      { 28, 1164.57 },
      { 29, 1119.84 },
      { 30, 1065.84 },
      { 31, 1038.89 },
      { 32, 1001.49 },
      { 33, 973.184 },
      { 34, 920.865 },
      { 35, 902.811 },
      { 36, 867.57 },
      { 37, 828.589 },
      { 38, 809.721 },
      { 39, 775.023 },
      { 40, 765.698 },
      { 41, 725.523 },
      { 42, 700.855 },
      { 43, 680.903 },
      { 44, 660.192 },
      { 45, 648.265 },
      { 46, 624.627 },
      { 47, 601.097 },
      { 48, 586.892 },
      { 49, 562.169 },
      { 50, 541.242 },
      { 51, 533.055 },
      { 52, 520.531 },
      { 53, 496.676 },
      { 54, 496.514 },
      { 55, 475.207 },
      { 56, 461.49 },
      { 57, 451.405 },
      { 58, 438.448 },
      { 59, 430.912 },
      { 60, 417.466 },
      { 61, 408.032 },
      { 62, 399.683 },
      { 63, 392.581 },
      { 64, 378.864 },
      { 65, 367.858 },
      { 66, 365.473 },
      { 67, 347.581 },
      { 68, 334.786 },
      { 69, 335.707 },
      { 70, 320.202 },
      { 71, 320.256 },
      { 72, 323.835 },
      { 73, 300.521 },
      { 74, 306.539 },
      { 75, 294.72 },
      { 76, 284.961 },
      { 77, 280.082 },
      { 78, 271.244 },
      { 79, 263.275 },
      { 80, 261.16 },
      { 81, 254.221 },
      { 82, 250.154 },
      { 83, 238.931 },
      { 84, 240.016 },
      { 85, 233.456 },
      { 86, 222.721 },
      { 87, 221.907 },
      { 88, 220.498 },
      { 89, 217.895 },
      { 90, 209.112 },
      { 91, 202.444 },
      { 92, 199.245 },
      { 93, 193.932 },
      { 94, 190.191 },
      { 95, 194.366 },
      { 96, 190.841 },
      { 97, 181.787 },
      { 98, 175.227 },
      { 99, 169.426 },
      { 100, 169.697 },
      { 101, 164.438 },
      { 102, 166.227 },
      { 103, 159.451 },
      { 104, 162.703 },
      { 105, 157.065 },
      { 106, 155.981 },
      { 107, 150.613 },
      { 108, 148.878 },
      { 109, 141.396 },
      { 110, 142.047 },
      { 111, 136.788 },
      { 112, 132.017 },
      { 113, 127.137 },
      { 114, 128.059 },
      { 115, 134.727 },
      { 116, 128.926 },
      { 117, 130.065 },
      { 118, 119.113 },
      { 119, 121.445 },
      { 120, 117.162 },
      { 121, 118.951 },
      { 122, 112.77 },
      { 123, 112.119 },
      { 124, 108.649 },
      { 125, 106.047 },
      { 126, 104.204 },
      { 127, 101.927 },
      { 128, 103.282 },
      { 129, 99.3244 },
      { 130, 92.2762 },
      { 131, 95.4207 },
      { 132, 96.234 },
      { 133, 99.4871 },
      { 134, 91.8425 },
      { 135, 90.3787 },
      { 136, 85.6618 },
      { 137, 84.1438 },
      { 138, 86.1498 },
      { 139, 84.6317 },
      { 140, 82.7884 },
      { 141, 80.6739 },
      { 142, 81.7583 },
      { 143, 75.7403 },
      { 144, 76.5535 },
      { 145, 80.6197 },
      { 146, 73.1921 },
      { 147, 72.0536 },
      { 148, 72.0536 },
      { 149, 71.6198 },
      { 150, 69.4512 },
      { 151, 71.1861 },
      { 152, 69.2343 },
      { 153, 63.9753 },
      { 154, 64.3548 },
      { 155, 64.2464 },
      { 156, 61.644 },
      { 157, 62.6199 },
      { 158, 60.1802 },
      { 159, 58.3368 }
    };
    static const std::vector< std::pair<double, double> > pt_dxsecNNLOdpt_pairs={
      { 0, 1352.7 },
      { 2, 3026.91 },
      { 4, 3474.95 },
      { 6, 3372.37 },
      { 8, 3124.28 },
      { 10, 2817.84 },
      { 12, 2527.44 },
      { 14, 2275.05 },
      { 16, 2054.14 },
      { 18, 1866.09 },
      { 20, 1716.01 },
      { 22, 1533.89 },
      { 24, 1412.98 },
      { 26, 1296.62 },
      { 28, 1197.47 },
      { 30, 1108.98 },
      { 32, 1028.53 },
      { 34, 950.774 },
      { 36, 894.662 },
      { 38, 828.648 },
      { 40, 784.771 },
      { 42, 728.839 },
      { 44, 686.004 },
      { 46, 652.744 },
      { 48, 611.308 },
      { 50, 578.607 },
      { 52, 548.151 },
      { 54, 513.692 },
      { 56, 487.203 },
      { 58, 464.965 },
      { 60, 440.997 },
      { 62, 424.509 },
      { 64, 395.259 },
      { 66, 384.622 },
      { 68, 368.326 },
      { 70, 350.313 },
      { 72, 331.228 },
      { 74, 319.046 },
      { 76, 303.94 },
      { 78, 288.201 },
      { 80, 279.139 },
      { 82, 266.371 },
      { 84, 254.341 },
      { 86, 247.346 },
      { 88, 231.843 },
      { 90, 226.656 },
      { 92, 216.102 },
      { 94, 207.329 },
      { 96, 203.623 },
      { 98, 191.013 },
      { 100, 185.838 },
      { 102, 178.242 },
      { 104, 169.752 },
      { 106, 163.426 },
      { 108, 157.387 },
      { 110, 153.889 },
      { 112, 147.664 },
      { 114, 142.535 },
      { 116, 136.604 },
      { 118, 132.476 },
      { 120, 127.013 },
      { 122, 124.151 },
      { 124, 116.889 },
      { 126, 115.057 },
      { 128, 110.377 },
      { 130, 106.281 },
      { 132, 101.908 },
      { 134, 99.4697 },
      { 136, 95.9682 },
      { 138, 91.841 },
      { 140, 89.0867 },
      { 142, 86.8707 },
      { 144, 84.005 },
      { 146, 80.9977 },
      { 148, 79.1574 },
      { 150, 77.4441 },
      { 152, 73.1775 },
      { 154, 71.7693 },
      { 156, 69.2995 },
      { 158, 67.2147 },
      { 160, 64.5131 },
      { 162, 63.2552 },
      { 164, 60.6735 },
      { 166, 58.782 },
      { 168, 57.1164 },
      { 170, 55.0556 },
      { 172, 53.1484 },
      { 174, 52.613 },
      { 176, 50.4792 },
      { 178, 49.9026 },
      { 180, 46.8545 },
      { 182, 46.553 },
      { 184, 44.1068 },
      { 186, 44.5581 },
      { 188, 42.2943 },
      { 190, 41.1089 },
      { 192, 39.4538 },
      { 194, 37.9085 },
      { 196, 38.823 },
      { 198, 36.3036 },
      { 200, 35.5537 },
      { 202, 34.0922 },
      { 204, 34.1132 },
      { 206, 31.9836 },
      { 208, 32.4309 },
      { 210, 29.7868 },
      { 212, 30.4933 },
      { 214, 28.9013 },
      { 216, 27.6841 },
      { 218, 27.935 },
      { 220, 26.4002 },
      { 222, 25.3582 },
      { 224, 25.6734 },
      { 226, 24.3061 },
      { 228, 24.1658 },
      { 230, 23.3597 },
      { 232, 22.3155 },
      { 234, 22.3404 },
      { 236, 21.4867 },
      { 238, 20.916 },
      { 240, 20.9879 },
      { 242, 19.1895 },
      { 244, 19.002 },
      { 246, 19.1301 },
      { 248, 18.334 },
      { 250, 17.4082 },
      { 252, 17.5854 },
      { 254, 16.4798 },
      { 256, 16.4874 },
      { 258, 16.57 },
      { 260, 15.7267 },
      { 262, 14.7291 },
      { 264, 15.0447 },
      { 266, 14.5847 },
      { 268, 14.1532 },
      { 270, 14.0834 },
      { 272, 13.2301 },
      { 274, 13.3985 },
      { 276, 12.6008 },
      { 278, 12.6924 },
      { 280, 12.0855 },
      { 282, 12.0805 },
      { 284, 11.8038 },
      { 286, 11.1154 },
      { 288, 11.4354 },
      { 290, 10.6224 },
      { 292, 10.612 },
      { 294, 10.3012 },
      { 296, 10.2084 },
      { 298, 10.1103 },
      { 300, 9.46106 },
      { 302, 9.42522 },
      { 304, 9.10904 },
      { 306, 8.77054 },
      { 308, 8.61 },
      { 310, 8.42472 },
      { 312, 8.35905 },
      { 314, 8.08369 },
      { 316, 7.96076 },
      { 318, 8.14242 },
      { 320, 7.55919 },
      { 322, 7.35281 },
      { 324, 7.59602 },
      { 326, 7.05912 },
      { 328, 7.23618 },
      { 330, 6.60041 },
      { 332, 6.80337 },
      { 334, 6.59755 },
      { 336, 6.26821 },
      { 338, 6.38093 },
      { 340, 6.08783 },
      { 342, 6.08606 },
      { 344, 5.72672 },
      { 346, 5.91723 },
      { 348, 5.59435 },
      { 350, 5.27735 },
      { 352, 5.33953 },
      { 354, 5.25781 },
      { 356, 4.94777 },
      { 358, 5.09145 },
      { 360, 4.90733 },
      { 362, 4.9163 },
      { 364, 4.73315 },
      { 366, 4.5936 },
      { 368, 4.57343 },
      { 370, 4.58277 },
      { 372, 4.39141 },
      { 374, 4.07482 },
      { 376, 4.30883 },
      { 378, 3.96164 },
      { 380, 3.98015 },
      { 382, 3.88802 },
      { 384, 3.92449 },
      { 386, 3.8463 },
      { 388, 3.57269 },
      { 390, 3.59467 },
      { 392, 3.54517 },
      { 394, 3.54852 },
      { 396, 3.44569 },
      { 398, 3.40717 },
      { 400, 3.07584 },
      { 402, 3.33014 },
      { 404, 2.99553 },
      { 406, 2.95431 },
      { 408, 3.11544 },
      { 410, 2.98542 },
      { 412, 2.86412 },
      { 414, 2.81536 },
      { 416, 2.83647 },
      { 418, 2.76095 },
      { 420, 2.74509 },
      { 422, 2.69556 },
      { 424, 2.63706 },
      { 426, 2.4831 },
      { 428, 2.60245 },
      { 430, 2.43423 },
      { 432, 2.22897 },
      { 434, 2.55714 },
      { 436, 2.15999 },
      { 438, 2.26069 },
      { 440, 2.35289 },
      { 442, 2.02684 },
      { 444, 2.29214 },
      { 446, 2.05598 },
      { 448, 2.06804 },
      { 450, 2.09245 },
      { 452, 2.04225 },
      { 454, 2.07283 },
      { 456, 1.88161 },
      { 458, 1.91933 },
      { 460, 1.85519 },
      { 462, 1.79867 },
      { 464, 1.70363 },
      { 466, 1.78892 },
      { 468, 1.84855 },
      { 470, 1.70237 },
      { 472, 1.6893 },
      { 474, 1.64749 },
      { 476, 1.56984 },
      { 478, 1.55488 },
      { 480, 1.62439 },
      { 482, 1.44774 },
      { 484, 1.45588 },
      { 486, 1.51294 },
      { 488, 1.4278 },
      { 490, 1.37181 },
      { 492, 1.46301 },
      { 494, 1.32757 },
      { 496, 1.42946 },
      { 498, 1.24183 }
    };
    if (order == 0){
      return 1.549; // From MATRIX manual, Table 6, process ppwxw02)
    }
    else if (order == 1){
      // Protection for different pT ranges
      pt = std::min(pt, std::min(pt_dxsecNLOdpt_pairs.back().first, pt_dxsecNNLOdpt_pairs.back().first));

      auto it_pt_xsecratio_pair = pt_dxsecNLOdpt_pairs.cbegin();
      auto itEnd_pt_xsecratio_pair = pt_dxsecNLOdpt_pairs.cend();
      for (; it_pt_xsecratio_pair != itEnd_pt_xsecratio_pair; it_pt_xsecratio_pair++){
        auto itNext_pt_xsecratio_pair = it_pt_xsecratio_pair; itNext_pt_xsecratio_pair++;
        if (pt>=it_pt_xsecratio_pair->first && (itNext_pt_xsecratio_pair==itEnd_pt_xsecratio_pair || pt<itNext_pt_xsecratio_pair->first)) break;
      }
      double const& xsec_NLO = it_pt_xsecratio_pair->second;

      it_pt_xsecratio_pair = pt_dxsecNNLOdpt_pairs.cbegin();
      itEnd_pt_xsecratio_pair = pt_dxsecNNLOdpt_pairs.cend();
      for (; it_pt_xsecratio_pair != itEnd_pt_xsecratio_pair; it_pt_xsecratio_pair++){
        auto itNext_pt_xsecratio_pair = it_pt_xsecratio_pair; itNext_pt_xsecratio_pair++;
        if (pt>=it_pt_xsecratio_pair->first && (itNext_pt_xsecratio_pair==itEnd_pt_xsecratio_pair || pt<itNext_pt_xsecratio_pair->first)) break;
      }
      double const& xsec_NNLO = it_pt_xsecratio_pair->second;

      return xsec_NNLO/xsec_NLO;
    }
    else{
      throw cms::Exception("UnknownOrder") << "KFactorHelpers::xsecRatio_qqWW4f_QCD_NLO_NNLO_13TeV_byPt: Order " << order << " is not implemented.";
      return 1;
    }
  }

  double xsecRatio_qqWZ4f_QCD_NLO_NNLO_13TeV_flat(bool isWpZ, unsigned char order){ // order: 0=NLO/LO, 1=NNLO/NLO
    // NLO/LO ratios are from the 2l2nu framework, NNLO/NLO are from the MATRIX manual, Table 6
    return (order == 0 ? (isWpZ ? 28.55/15.51 : 18.19/9.53) : (isWpZ ? 1.106 : 1.111));
  }


  const std::string KFactorHandler_QCD_ggVV_Sig::KFactorArgName = "KFactor_QCD_ggVV_Sig_arg";
  KFactorHandler_QCD_ggVV_Sig::KFactorHandler_QCD_ggVV_Sig(int const& year) :
    KFactorHandlerBase(),
    kfFile_NNLO(nullptr),
    kfFile_NLO(nullptr)
  {
    TString strKFactorDir = "${CMSSW_BASE}/src/CMS3/NtupleMaker/data/Kfactors/";
    HostHelpers::ExpandEnvironmentVariables(strKFactorDir);

    TString strSqrts="";
    switch (year){
    case 2015:
    case 2016:
    case 2017:
    case 2018:
      strSqrts="13TeV";
      break;
    default:
      throw cms::Exception("UnknownYear") << "KFactorHandler_QCD_ggVV_Sig::KFactorHandler_QCD_ggVV_Sig: Year " << year << " is not implemented.";
      break;
    }
    strKFactorDir += strSqrts;

    TString strKFactorNNLOFile = strKFactorDir + "/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root";
    TString strKFactorNLOFile = strKFactorDir + "/Kfactor_Collected_ggHZZ_2l2l_NLO_NNPDF_NarrowWidth_13TeV.root";

    TDirectory* curdir = gDirectory;
    kfFile_NNLO = std::make_shared<TFile>(strKFactorNNLOFile, "read"); curdir->cd();
    kfFile_NLO = std::make_shared<TFile>(strKFactorNLOFile, "read"); curdir->cd();

    this->setup();
  }
  KFactorHandler_QCD_ggVV_Sig::KFactorHandler_QCD_ggVV_Sig(KFactorHandler_QCD_ggVV_Sig const& other) :
    KFactorHandlerBase(other),
    kfFile_NNLO(other.kfFile_NNLO),
    kfFile_NLO(other.kfFile_NLO)
  {
    // Do not copy vectors of splines, get them from scratch again
    this->setup();
  }
  void KFactorHandler_QCD_ggVV_Sig::setup(){
    TDirectory* curdir = gDirectory;

    kfFile_NNLO->cd();
    sp_NNLO.reserve(9);
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_Nominal"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_QCDScaleDn"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_QCDScaleUp"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_PDFScaleDn"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_PDFScaleUp"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_PDFReplicaDn"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_PDFReplicaUp"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_AsDn"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_AsUp"));
    for (size_t ikf=0; ikf<sp_NNLO.size(); ikf++){
      if (!sp_NNLO.at(ikf)) throw cms::Exception(Form("KFactorHandler_QCD_ggVV_Sig::setup: NNLO K factor at location %lu cannot be found.", ikf));
    }

    kfFile_NLO->cd();
    sp_NLO.reserve(9);
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_Nominal"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_QCDScaleDn"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_QCDScaleUp"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_PDFScaleDn"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_PDFScaleUp"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_PDFReplicaDn"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_PDFReplicaUp"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_AsDn"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_AsUp"));
    for (size_t ikf=0; ikf<sp_NLO.size(); ikf++){
      if (!sp_NLO.at(ikf)) throw cms::Exception(Form("KFactorHandler_QCD_ggVV_Sig::setup: NLO K factor at location %lu cannot be found.", ikf));
    }

    curdir->cd();
  }
  void KFactorHandler_QCD_ggVV_Sig::eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType denominator, std::unordered_map<std::string, float>& kfactors_map) const{
    static const std::vector<std::string> kfactornames{
      "KFactor_QCD_NNLO_ggVV_Sig_Nominal",
      "KFactor_QCD_NNLO_ggVV_Sig_QCDScaleDn",
      "KFactor_QCD_NNLO_ggVV_Sig_QCDScaleUp",
      "KFactor_QCD_NNLO_ggVV_Sig_PDFScaleDn",
      "KFactor_QCD_NNLO_ggVV_Sig_PDFScaleUp",
      "KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaDn",
      "KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaUp",
      "KFactor_QCD_NNLO_ggVV_Sig_AsDn",
      "KFactor_QCD_NNLO_ggVV_Sig_AsUp",
      "KFactor_QCD_NLO_ggVV_Sig_Nominal",
      "KFactor_QCD_NLO_ggVV_Sig_QCDScaleDn",
      "KFactor_QCD_NLO_ggVV_Sig_QCDScaleUp",
      "KFactor_QCD_NLO_ggVV_Sig_PDFScaleDn",
      "KFactor_QCD_NLO_ggVV_Sig_PDFScaleUp",
      "KFactor_QCD_NLO_ggVV_Sig_PDFReplicaDn",
      "KFactor_QCD_NLO_ggVV_Sig_PDFReplicaUp",
      "KFactor_QCD_NLO_ggVV_Sig_AsDn",
      "KFactor_QCD_NLO_ggVV_Sig_AsUp"
    };
    assert(kfactornames.size()==(this->sp_NNLO.size() + this->sp_NLO.size()));

    auto it_arg = kfactors_map.find(KFactorHandler_QCD_ggVV_Sig::KFactorArgName);
    if (it_arg == kfactors_map.end()) throw cms::Exception(Form("KFactorHandler_QCD_ggVV_Sig::eval: K factor evaluation argument, candidate mass with name %s, cannot be found.", KFactorHandler_QCD_ggVV_Sig::KFactorArgName.data()));
    float const& kfactor_arg = it_arg->second;

    float kfactor_denominator = 1;
    switch (denominator){
    case kf_QCD_NNLO_GGVV_SIG:
      kfactor_denominator = sp_NNLO.front()->Eval(kfactor_arg);
      break;
    case kf_QCD_NLO_GGVV_SIG:
      kfactor_denominator = sp_NLO.front()->Eval(kfactor_arg);
      break;
    default:
      break;
    }

    size_t ikf=0;
    if (type == kf_QCD_NNLO_GGVV_SIG){
      for (auto const& spkf:sp_NNLO){
        kfactors_map[kfactornames.at(ikf)] = spkf->Eval(kfactor_arg) / kfactor_denominator;
        ikf++;
      }
    }
    ikf = sp_NNLO.size();
    if (type == kf_QCD_NLO_GGVV_SIG){
      for (auto const& spkf:sp_NLO){
        kfactors_map[kfactornames.at(ikf)] = spkf->Eval(kfactor_arg) / kfactor_denominator;
        ikf++;
      }
    }
  }


  KFactorHandler_EW_qqVV_Bkg::KFactorHandler_EW_qqVV_Bkg() :
    KFactorHandlerBase(),
    type(nKFactorTypes),
    beamEnergy(0),

    fcn_Kfactor_QCD_qqZZ(nullptr),
    fcn_Kfactor_QCD_qqWZ(nullptr),
    fcn_Kfactor_QCD_qqWW(nullptr)
  {}
  KFactorHandler_EW_qqVV_Bkg::KFactorHandler_EW_qqVV_Bkg(int const& year, KFactorHelpers::KFactorType const& type_) :
    KFactorHandlerBase(),
    type(type_),
    beamEnergy(0),

    fcn_Kfactor_QCD_qqZZ(nullptr),
    fcn_Kfactor_QCD_qqWZ(nullptr),
    fcn_Kfactor_QCD_qqWW(nullptr)
  {
    TString strKFactorDir = "${CMSSW_BASE}/src/CMS3/NtupleMaker/data/Kfactors/";
    HostHelpers::ExpandEnvironmentVariables(strKFactorDir);
    unsigned int sqrts=0;
    switch(year){
    case 2011:
      sqrts=7;
      throw cms::Exception("UnknownKFactorQCD") << "KFactorHandler_EW_qqVV_Bkg::KFactorHandler_EW_qqVV_Bkg: sqrts for year " << year << " does not have QCD K factor functions implemented.";
      break;
    case 2012:
      sqrts=8;
      throw cms::Exception("UnknownKFactorQCD") << "KFactorHandler_EW_qqVV_Bkg::KFactorHandler_EW_qqVV_Bkg: sqrts for year " << year << " does not have QCD K factor functions implemented.";
      break;
    case 2015:
    case 2016:
    case 2017:
    case 2018:
      sqrts=13;
      fcn_Kfactor_QCD_qqZZ = &(KFactorHelpers::xsecRatio_qqZZ4f_QCD_NLO_NNLO_13TeV_byMass);
      fcn_Kfactor_QCD_qqWZ = &(KFactorHelpers::xsecRatio_qqWZ4f_QCD_NLO_NNLO_13TeV_flat);
      fcn_Kfactor_QCD_qqWW = &(KFactorHelpers::xsecRatio_qqWW4f_QCD_NLO_NNLO_13TeV_byPt);
      break;
    default:
      throw cms::Exception("UnknownYear") << "KFactorHandler_EW_qqVV_Bkg::KFactorHandler_EW_qqVV_Bkg: sqrts for year " << year << " is not implemented.";
      break;
    }
    strKFactorDir += Form("%uTeV", sqrts);

    TString strKFactorFile;
    if (type == kf_EW_NLO_QQZZ_BKG) strKFactorFile = strKFactorDir + "/ZZ_EWCorrections.dat";
    else if (type == kf_EW_NLO_QQWZ_BKG) strKFactorFile = strKFactorDir + "/WZ_EWCorrections.dat";
    else if (type == kf_EW_NLO_QQWW_BKG) strKFactorFile = strKFactorDir + "/WW_EWCorrections.dat";
    else throw cms::Exception("UnknownType") << "KFactorHandler_EW_qqVV_Bkg::KFactorHandler_EW_qqVV_Bkg: Correction type " << type << " is not implemented.";
    
    readTableFromFile(strKFactorFile, table_VV);
    beamEnergy = getBeamEnergy(year);

    this->setup();
  }
  KFactorHandler_EW_qqVV_Bkg::KFactorHandler_EW_qqVV_Bkg(KFactorHandler_EW_qqVV_Bkg const& other) :
    KFactorHandlerBase(other),
    type(other.type),
    beamEnergy(other.beamEnergy),
    fcn_Kfactor_QCD_qqZZ(other.fcn_Kfactor_QCD_qqZZ),
    fcn_Kfactor_QCD_qqWZ(other.fcn_Kfactor_QCD_qqWZ),
    fcn_Kfactor_QCD_qqWW(other.fcn_Kfactor_QCD_qqWW),
    table_VV(other.table_VV)
  {
    this->setup();
  }
  
  void KFactorHandler_EW_qqVV_Bkg::readTableFromFile(TString const& fname, std::vector< std::vector<double> >& table){
    //MELAout << "KFactorHandler_EW_qqVV_Bkg::readTableFromFile: Reading table " << fname << endl;
    if (!HostHelpers::FileReadable(fname.Data())) throw cms::Exception("IOError") << "KFactorHandler_EW_qqVV_Bkg::readTableFromFile " << fname << " does not exist!";
    ifstream fin(fname.Data(), ios_base::in);
    while (!fin.eof()){
      std::string strline;
      std::getline(fin, strline);
      if (strline.empty()) continue;

      //MELAout << "=> " << strline << endl;

      table.push_back(std::vector<double>());
      std::vector<double>& tbl_line = table.back(); tbl_line.reserve(5);

      std::stringstream ss(strline);
      std::string strtmp;
      while (ss >> strtmp){
        if (strtmp.empty()) continue;
        tbl_line.push_back(std::stod(strtmp));
      }
      //MELAout << "==>> " << tbl_line << endl;
    }
    //MELAout << "KFactorHandler_EW_qqVV_Bkg::readTableFromFile: Done reading table" << endl;
  }
  std::vector<double> KFactorHandler_EW_qqVV_Bkg::findTableEntry(double const& mhat, double const& that) const{
    std::vector<std::vector<double>>::const_iterator const table_end = table_VV.cend();

    double bestShatDiff = -1;
    std::vector<std::vector<double>>::const_iterator it_bestShat = table_end;
    std::vector<std::vector<double>>::const_iterator itNext_bestShat = table_end;
    for (size_t is=0; is<table_sqrtShatBegin_VV.size()-1;is++){
      auto const& itFirst = table_sqrtShatBegin_VV.at(is);
      auto const& itSecond = table_sqrtShatBegin_VV.at(is+1);
      double tmpShatDiff = std::abs(itFirst->front() - mhat);
      if (bestShatDiff<0. || tmpShatDiff < bestShatDiff){
        bestShatDiff = tmpShatDiff;
        it_bestShat = itFirst;
        itNext_bestShat = itSecond;
      }
    }
    assert(it_bestShat != table_end);

    std::vector<std::vector<double>>::const_iterator it_That = it_bestShat;
    std::vector<std::vector<double>>::const_iterator it_bestThat = table_end;
    double bestThatDiff = -1;
    while (it_That != itNext_bestShat){
      double tmpThatDiff = std::abs(it_That->at(1) - that);
      if (bestThatDiff<0. || tmpThatDiff<bestThatDiff){
        bestThatDiff = tmpThatDiff;
        it_bestThat = it_That;
      }
      it_That++;
    }
    assert(it_bestThat != table_end);

    double const& val_uc = it_bestThat->at(2);
    double const& val_ds = it_bestThat->at(3);
    double const& val_b = it_bestThat->at(4);
    return std::vector<double>{val_uc, val_ds, val_b};
  }

  void KFactorHandler_EW_qqVV_Bkg::setup(){
    double sqrtShat = -1;
    for (std::vector<std::vector<double>>::const_iterator it = table_VV.cbegin(); it!=table_VV.cend(); it++){
      double const& tmpSqrtShat = it->front();
      if (tmpSqrtShat!=sqrtShat){
        sqrtShat = tmpSqrtShat;
        table_sqrtShatBegin_VV.push_back(it);
      }
    }
    table_sqrtShatBegin_VV.push_back(table_VV.cend());
  }

  void KFactorHandler_EW_qqVV_Bkg::eval(
    GenEventInfoProduct const& eventInfo,
    std::vector<reco::GenParticle const*> const& incomingQuarks,
    std::pair<reco::GenParticle const*, reco::GenParticle const*> const& bestV1pair,
    std::pair<reco::GenParticle const*, reco::GenParticle const*> const& bestV2pair,
    std::unordered_map<std::string, float>& kfactors_map
  ) const{
    static const float mZsq = PDGHelpers::Zmass * PDGHelpers::Zmass;
    static const float mWsq = PDGHelpers::Wmass * PDGHelpers::Wmass;
    static const std::vector<std::string> kfactornames{
      "KFactor_EW_NLO_qqVV_Bkg_Nominal",
      "KFactor_EW_NLO_qqVV_Bkg_EWDn",
      "KFactor_EW_NLO_qqVV_Bkg_EWUp"
    };
    static const std::vector<std::string> kfactorargnames{
      "KFactor_EW_NLO_qqVV_Bkg_arg_mass",
      "KFactor_EW_NLO_qqVV_Bkg_arg_that",
      "KFactor_EW_NLO_qqVV_Bkg_arg_pthat",
      "KFactor_EW_NLO_qqVV_Bkg_arg_rho"
    };

    if (!eventInfo.pdf()) throw cms::Exception("InvalidGenEventInfoProduct") << "KFactorHandler_EW_qqVV_Bkg::eval: eventInfo is invalid.";
    if (!(bestV1pair.first && bestV1pair.second)) throw cms::Exception("InvalidKinematics") << "KFactorHandler_EW_qqVV_Bkg::eval: V1 is invalid.";
    if (!(bestV2pair.first && bestV2pair.second)) throw cms::Exception("InvalidKinematics") << "KFactorHandler_EW_qqVV_Bkg::eval: V2 is invalid.";

    if (
      !eventInfo.pdf()
      ||
      !(bestV1pair.first && bestV1pair.second && bestV2pair.first && bestV2pair.second)
      ){
      for (auto const& kfactorname:kfactornames) kfactors_map[kfactorname]=1;
      for (auto const& kfactorname:kfactorargnames) kfactors_map[kfactorname]=0;
      return;
    }
    auto const& x1 = eventInfo.pdf()->x.first;
    auto const& x2 = eventInfo.pdf()->x.second;

    float sum_mV1_mV2 = 0;
    switch (type){
    case kf_EW_NLO_QQZZ_BKG:
      sum_mV1_mV2 = 2.*PDGHelpers::Zmass;
      break;
    case kf_EW_NLO_QQWZ_BKG:
      sum_mV1_mV2 = PDGHelpers::Wmass + PDGHelpers::Zmass;
      break;
    case kf_EW_NLO_QQWW_BKG:
      sum_mV1_mV2 = 2.*PDGHelpers::Wmass;
      break;
    default:
      // Do nothing
      break;
    }

    auto pBestV1 = bestV1pair.first->p4() + bestV1pair.second->p4();
    auto pBestV2 = bestV2pair.first->p4() + bestV2pair.second->p4();
    auto pVV = pBestV1 + pBestV2;
    TLorentzVector pVV_tlv(pVV.px(), pVV.py(), pVV.pz(), pVV.energy());
    double m_hat = pVV.M();
    double pt_hat = pVV.Pt();
    if (m_hat<sum_mV1_mV2 || incomingQuarks.empty()){
      for (auto const& kfactorname:kfactornames) kfactors_map[kfactorname]=1;
      kfactors_map[kfactorargnames.at(0)]=m_hat;
      kfactors_map[kfactorargnames.at(1)]=0;
      kfactors_map[kfactorargnames.at(2)]=pt_hat;
      kfactors_map[kfactorargnames.at(3)]=0;
      return;
    }
    double s_hat = std::pow(m_hat, 2);
    double rho = pVV.Pt() / (bestV1pair.first->pt() + bestV1pair.second->pt() + bestV2pair.first->pt() + bestV2pair.second->pt());

    TLorentzVector p1_b, p2_b, pV1_b;
    TVector3 boost_VV = -pVV_tlv.BoostVector();
    pV1_b.SetXYZT(pBestV1.px(), pBestV1.py(), pBestV1.pz(), pBestV1.energy()); pV1_b.Boost(boost_VV);
    p1_b.SetXYZT(0., 0., x1*beamEnergy, x1*beamEnergy); p1_b.Boost(boost_VV);
    p2_b.SetXYZT(0., 0., -x2*beamEnergy, x2*beamEnergy); p2_b.Boost(boost_VV);
    TVector3 nV1_b = pV1_b.Vect().Unit();
    TVector3 n1_b = p1_b.Vect().Unit();
    TVector3 n2_b = p2_b.Vect().Unit();
    TVector3 n_beam = n1_b - n2_b; n_beam = n_beam.Unit();
    double cos_theta = n_beam.Dot(nV1_b);
    double t_hat=0;
    if (type == kf_EW_NLO_QQWW_BKG) t_hat = mWsq - 0.5*s_hat + cos_theta * std::sqrt(0.25*s_hat*s_hat - mWsq*s_hat);
    else if (type == kf_EW_NLO_QQZZ_BKG) t_hat = mZsq - 0.5*s_hat + cos_theta * std::sqrt(0.25*s_hat*s_hat - mZsq*s_hat);
    else{
      double b = 0.5/m_hat * std::sqrt(std::pow(s_hat - mZsq - mWsq, 2) - 4.*mZsq*mWsq);
      double a = std::sqrt(b*b + mZsq);
      t_hat = mZsq - m_hat * (a - b * cos_theta); // Awful calculation, needed to put ourselves to the center-of-mass frame with the 2 particles having a different mass!
    }

    unsigned char qtype = std::abs(incomingQuarks.front()->pdgId());
    std::vector<double> vCorr = findTableEntry(m_hat, t_hat);

    double kfactor = 1. + vCorr.at(1*(qtype==1 || qtype==3) + 2*(qtype==5));
    double rel_error = 0;
    double kfactor_QCD = 1;
    double kfactor_virtphoton = 1;
    if (type == kf_EW_NLO_QQZZ_BKG){ // ZZ
      kfactor_QCD = fcn_Kfactor_QCD_qqZZ(m_hat, !(std::abs(bestV1pair.first->pdgId())==std::abs(bestV2pair.first->pdgId())), 0);
    }
    else if (type == kf_EW_NLO_QQWW_BKG){ // WW
      kfactor_QCD = fcn_Kfactor_QCD_qqWW(pt_hat, 0);
    }
    else if (PDGHelpers::isANeutrino(bestV2pair.first->pdgId()) || PDGHelpers::isUpTypeQuark(bestV2pair.first->pdgId())){ // W+Z
      kfactor_QCD = fcn_Kfactor_QCD_qqWZ(true, 0);
      kfactor_virtphoton = (1. + 0.00559445 - 5.17082e-6 * m_hat + 3.63331e-8 * s_hat);
    }
    else{ // W-Z
      kfactor_QCD = fcn_Kfactor_QCD_qqWZ(false, 0);
      kfactor_virtphoton = (1. + 0.00174737 + 1.70668e-5 * m_hat + 2.26398e-8 * s_hat);
    }

    if (rho<0.3) rel_error = std::abs((kfactor-1.)*(kfactor_QCD - 1.)); //If rho is small, only corrections in QCD X EWK
    else rel_error = std::abs(kfactor-1.); // If rho is large, 100% because of large colinear gluon radiations

    kfactor *= kfactor_virtphoton; // No associated unc.

    kfactors_map[kfactornames.at(0)]=kfactor;
    kfactors_map[kfactornames.at(1)]=kfactor*(1.f - rel_error);
    kfactors_map[kfactornames.at(2)]=kfactor*(1.f + rel_error);

    kfactors_map[kfactorargnames.at(0)]=m_hat;
    kfactors_map[kfactorargnames.at(1)]=t_hat;
    kfactors_map[kfactorargnames.at(2)]=pt_hat;
    kfactors_map[kfactorargnames.at(3)]=rho;
  }


  KFactorHandler_QCD_qqVV_Bkg::KFactorHandler_QCD_qqVV_Bkg() :
    KFactorHandlerBase(),

    fcn_Kfactor_QCD_qqZZ(nullptr),
    fcn_Kfactor_QCD_qqWZ(nullptr),
    fcn_Kfactor_QCD_qqWW(nullptr)
  {}
  KFactorHandler_QCD_qqVV_Bkg::KFactorHandler_QCD_qqVV_Bkg(int const& year) :
    KFactorHandlerBase(),

    fcn_Kfactor_QCD_qqZZ(nullptr),
    fcn_Kfactor_QCD_qqWZ(nullptr),
    fcn_Kfactor_QCD_qqWW(nullptr)
  {
    switch (year){
    case 2015:
    case 2016:
    case 2017:
    case 2018:
      fcn_Kfactor_QCD_qqZZ = &(KFactorHelpers::xsecRatio_qqZZ4f_QCD_NLO_NNLO_13TeV_byMass);
      fcn_Kfactor_QCD_qqWZ = &(KFactorHelpers::xsecRatio_qqWZ4f_QCD_NLO_NNLO_13TeV_flat);
      fcn_Kfactor_QCD_qqWW = &(KFactorHelpers::xsecRatio_qqWW4f_QCD_NLO_NNLO_13TeV_byPt);
      break;
    default:
      throw cms::Exception("UnknownYear") << "KFactorHandler_QCD_qqVV_Bkg::KFactorHandler_QCD_qqVV_Bkg: Year " << year << " is not implemented.";
      break;
    }
  }
  KFactorHandler_QCD_qqVV_Bkg::KFactorHandler_QCD_qqVV_Bkg(KFactorHandler_QCD_qqVV_Bkg const& other) :
    KFactorHandlerBase(other),
    fcn_Kfactor_QCD_qqZZ(other.fcn_Kfactor_QCD_qqZZ),
    fcn_Kfactor_QCD_qqWZ(other.fcn_Kfactor_QCD_qqWZ),
    fcn_Kfactor_QCD_qqWW(other.fcn_Kfactor_QCD_qqWW)
  {}
  void KFactorHandler_QCD_qqVV_Bkg::eval(
    KFactorHelpers::KFactorType type,
    std::pair<reco::GenParticle const*, reco::GenParticle const*> const& bestV1pair,
    std::pair<reco::GenParticle const*, reco::GenParticle const*> const& bestV2pair,
    std::unordered_map<std::string, float>& kfactors_map
  ) const{
    static const std::vector<std::string> kfactornames{
      "KFactor_QCD_NNLO_qqVV_Bkg_Nominal"
    };
    static const std::vector<std::string> kfactorargnames{
      "KFactor_QCD_NNLO_qqVV_Bkg_arg_mass",
      "KFactor_QCD_NNLO_qqVV_Bkg_arg_pthat"
    };

    if (!(bestV1pair.first && bestV1pair.second)) throw cms::Exception("InvalidKinematics") << "KFactorHandler_QCD_qqVV_Bkg::eval: V1 is invalid.";
    if (!(bestV2pair.first && bestV2pair.second)) throw cms::Exception("InvalidKinematics") << "KFactorHandler_QCD_qqVV_Bkg::eval: V2 is invalid.";

    auto pBestV1 = bestV1pair.first->p4() + bestV1pair.second->p4();
    auto pBestV2 = bestV2pair.first->p4() + bestV2pair.second->p4();
    auto pVV = pBestV1 + pBestV2;
    double m_hat = pVV.M();
    double pt_hat = pVV.Pt();

    double kfactor=1;
    if (type == kf_QCD_NNLO_QQZZ_BKG){ // ZZ
      kfactor = fcn_Kfactor_QCD_qqZZ(m_hat, !(std::abs(bestV1pair.first->pdgId())==std::abs(bestV2pair.first->pdgId())), 1);
    }
    else if (type == kf_QCD_NNLO_QQWW_BKG){ // WW
      kfactor = fcn_Kfactor_QCD_qqWW(pt_hat, 1);
    }
    else if (type == kf_QCD_NNLO_QQWZ_BKG){ // WZ
      if (PDGHelpers::isANeutrino(bestV2pair.first->pdgId()) || PDGHelpers::isUpTypeQuark(bestV2pair.first->pdgId())){ // W+Z
        kfactor = fcn_Kfactor_QCD_qqWZ(true, 1);
      }
      else{ // W-Z
        kfactor = fcn_Kfactor_QCD_qqWZ(false, 1);
      }
    }

    kfactors_map[kfactornames.at(0)]=kfactor;

    kfactors_map[kfactorargnames.at(0)]=m_hat;
    kfactors_map[kfactorargnames.at(1)]=pt_hat;
  }

}
