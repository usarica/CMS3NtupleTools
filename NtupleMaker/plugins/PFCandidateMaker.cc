#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// #include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CMS3/NtupleMaker/interface/plugins/PFCandidateMaker.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

PFCandidateMaker::PFCandidateMaker(const edm::ParameterSet& iConfig){

  pfCandidatesToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesTag"));
  minPt_            = iConfig.getParameter<double>          ("minPt"              );

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<vector<LorentzVector>> (branchprefix+"p4"                         ).setBranchAlias(aliasprefix_+"_p4"                         );
  produces<vector<float>>         (branchprefix+"mass"                       ).setBranchAlias(aliasprefix_+"_mass"                       );
  produces<vector<float>>         (branchprefix+"dz"                         ).setBranchAlias(aliasprefix_+"_dz"                         );
  produces<vector<float>>         (branchprefix+"dxy"                        ).setBranchAlias(aliasprefix_+"_dxy"                        );
  produces<vector<float>>         (branchprefix+"dzError"                    ).setBranchAlias(aliasprefix_+"_dzError"                    );
  produces<vector<float>>         (branchprefix+"dxyError"                   ).setBranchAlias(aliasprefix_+"_dxyError"                   );
  produces<vector<int>>           (branchprefix+"charge"                     ).setBranchAlias(aliasprefix_+"_charge"                     );
  produces<vector<int>>           (branchprefix+"particleId"                 ).setBranchAlias(aliasprefix_+"_particleId"                 );
  produces<vector<uint8_t>>       (branchprefix+"fromPV"                     ).setBranchAlias(aliasprefix_+"_fromPV"                     );
  produces<vector<bool>>          (branchprefix+"isStandAloneMuon"           ).setBranchAlias(aliasprefix_+"_isStandAloneMuon"           );
  produces<vector<bool>>          (branchprefix+"isGlobalMuon"               ).setBranchAlias(aliasprefix_+"_isGlobalMuon"               );
  produces<vector<uint8_t>>       (branchprefix+"pvAssociationQuality"       ).setBranchAlias(aliasprefix_+"_pvAssociationQuality"       );
  produces<vector<int>>           (branchprefix+"IdAssociatedPV"             ).setBranchAlias(aliasprefix_+"_IdAssociatedPV"             );
  produces<vector<float>>         (branchprefix+"dzAssociatedPV"             ).setBranchAlias(aliasprefix_+"_dzAssociatedPV"             );
  produces<vector<float>>         (branchprefix+"puppiWeight"                ).setBranchAlias(aliasprefix_+"_puppiWeight"                );
  produces<vector<float>>         (branchprefix+"puppiWeightNoLep"           ).setBranchAlias(aliasprefix_+"_puppiWeightNoLep"           );
  produces<vector<float>>         (branchprefix+"trackIso"                   ).setBranchAlias(aliasprefix_+"_trackIso"                   );
  produces<vector<float>>         (branchprefix+"miniTrackIso"               ).setBranchAlias(aliasprefix_+"_miniTrackIso"               );
  // produces<vector<uint8_t>>       (branchprefix+"packedHits"                 ).setBranchAlias(aliasprefix_+"_packedHits"                 );
  // produces<vector<uint8_t>>       (branchprefix+"packedLayers"               ).setBranchAlias(aliasprefix_+"_packedLayers"               );
  // produces<vector<uint16_t>>      (branchprefix+"qualityFlags"               ).setBranchAlias(aliasprefix_+"_qualityFlags"               );
  // produces<vector<uint8_t>>       (branchprefix+"normalizedChi2"             ).setBranchAlias(aliasprefix_+"_normalizedChi2"             );
  produces<vector<int>>           (branchprefix+"numberOfPixelHits"          ).setBranchAlias(aliasprefix_+"_numberOfPixelHits"          );
  produces<vector<int>>           (branchprefix+"numberOfHits"               ).setBranchAlias(aliasprefix_+"_numberOfHits"               );
  produces<vector<int>>           (branchprefix+"pixelLayersWithMeasurement" ).setBranchAlias(aliasprefix_+"_pixelLayersWithMeasurement" );
  produces<vector<int>>           (branchprefix+"stripLayersWithMeasurement" ).setBranchAlias(aliasprefix_+"_stripLayersWithMeasurement" );
  produces<vector<bool>>          (branchprefix+"trackHighPurity"            ).setBranchAlias(aliasprefix_+"_trackHighPurity"            );

}

PFCandidateMaker::~PFCandidateMaker(){}
void  PFCandidateMaker::beginRun(const edm::Run&, const edm::EventSetup& es){}
void PFCandidateMaker::beginJob() {}
void PFCandidateMaker::endJob()   {}

void PFCandidateMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  unique_ptr<vector<LorentzVector>> pfcands_p4                         (new vector<LorentzVector> );
  unique_ptr<vector<float>>         pfcands_mass                       (new vector<float>         );
  unique_ptr<vector<float>>         pfcands_dz                         (new vector<float>         );
  unique_ptr<vector<float>>         pfcands_dxy                        (new vector<float>         );
  unique_ptr<vector<float>>         pfcands_dzError                    (new vector<float>         );
  unique_ptr<vector<float>>         pfcands_dxyError                   (new vector<float>         );
  unique_ptr<vector<int>>           pfcands_charge                     (new vector<int>           );
  unique_ptr<vector<int>>           pfcands_particleId                 (new vector<int>           );
  unique_ptr<vector<bool>>          pfcands_isStandAloneMuon           (new vector<bool>          );
  unique_ptr<vector<bool>>          pfcands_isGlobalMuon               (new vector<bool>          );
  unique_ptr<vector<uint8_t>>       pfcands_fromPV                     (new vector<uint8_t>       );
  unique_ptr<vector<uint8_t>>       pfcands_pvAssociationQuality       (new vector<uint8_t>       );
  unique_ptr<vector<int>>           pfcands_IdAssociatedPV             (new vector<int>           );
  unique_ptr<vector<float>>         pfcands_dzAssociatedPV             (new vector<float>         );
  unique_ptr<vector<float>>         pfcands_puppiWeight                (new vector<float>         );
  unique_ptr<vector<float>>         pfcands_puppiWeightNoLep           (new vector<float>         );
  unique_ptr<vector<float>>         pfcands_trackIso                   (new vector<float>         );
  unique_ptr<vector<float>>         pfcands_miniTrackIso               (new vector<float>         );
  unique_ptr<vector<int>>           pfcands_numberOfPixelHits          (new vector<int>           );
  unique_ptr<vector<int>>           pfcands_numberOfHits               (new vector<int>           );
  unique_ptr<vector<int>>           pfcands_pixelLayersWithMeasurement (new vector<int>           );
  unique_ptr<vector<int>>           pfcands_stripLayersWithMeasurement (new vector<int>           );
  unique_ptr<vector<bool>>          pfcands_trackHighPurity            (new vector<bool>          );

  // unique_ptr<vector<uint8_t>>       pfcands_packedHits          (new vector<uint8_t>       );
  // unique_ptr<vector<uint8_t>>       pfcands_packedLayers        (new vector<uint8_t>       );
  // unique_ptr<vector<uint16_t>>      pfcands_qualityFlags        (new vector<uint16_t>      );
  // unique_ptr<vector<uint8_t>>       pfcands_normalizedChi2      (new vector<uint8_t>       );

  unique_ptr<vector<float>>         pfcands_helperPhi           (new vector<float>         );
  unique_ptr<vector<float>>         pfcands_helperEta           (new vector<float>         );
  unique_ptr<vector<bool>>          pfcands_helperPVinfo        (new vector<bool>          );

  //get pfcandidates
  Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(pfCandidatesToken, pfCandidatesHandle);
  pfCandidates  = pfCandidatesHandle.product();

  for (pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++) {
    pfcands_helperPhi    -> push_back( pf_it->p4().Phi() );
    pfcands_helperEta    -> push_back( pf_it->p4().Eta() );
    pfcands_helperPVinfo -> push_back( (!pf_it->vertexRef().isNull() && fabs(pf_it->dz()) < 0.1) || pf_it->fromPV() > 1 );
    if (pf_it->p4().pt() < minPt_ && abs(pf_it->pdgId()) != 11 && abs(pf_it->pdgId()) != 13) continue;
    pfcands_p4   -> push_back( LorentzVector(pf_it->p4()) );
    pfcands_mass -> push_back( pf_it->mass()              );

    if (!pf_it->vertexRef().isNull()){
      pfcands_dz                   -> push_back( pf_it->dz()                   );
      pfcands_dxy                  -> push_back( pf_it->dxy()                  );
      if (pf_it->hasTrackDetails()) {
          pfcands_dzError              -> push_back( pf_it->dzError()              );
          pfcands_dxyError             -> push_back( pf_it->dxyError()             );
      } else {
          pfcands_dzError              -> push_back(0.);
          pfcands_dxyError             -> push_back(0.);
      }
      pfcands_pvAssociationQuality -> push_back( pf_it->pvAssociationQuality() );
      pfcands_dzAssociatedPV       -> push_back( pf_it->dzAssociatedPV()       );
      pfcands_IdAssociatedPV       -> push_back( pf_it->vertexRef().key()      );
    } else {
      pfcands_dz                   -> push_back( -9999.                        );
      pfcands_pvAssociationQuality -> push_back( 0                             );
      pfcands_dzAssociatedPV       -> push_back( -9999.                        );
      pfcands_IdAssociatedPV       -> push_back( -9999                         );
    }

    pfcands_charge                     -> push_back( pf_it->charge()                     );
    pfcands_particleId                 -> push_back( pf_it->pdgId()                      );
    pfcands_fromPV                     -> push_back( pf_it->fromPV()                     );
    pfcands_puppiWeight                -> push_back( pf_it->puppiWeight()                );
    pfcands_puppiWeightNoLep           -> push_back( pf_it->puppiWeightNoLep()           );
    pfcands_isStandAloneMuon           -> push_back( pf_it->isStandAloneMuon()           );
    pfcands_isGlobalMuon               -> push_back( pf_it->isGlobalMuon()               );
    pfcands_numberOfPixelHits          -> push_back( pf_it->numberOfPixelHits()          );
    pfcands_numberOfHits               -> push_back( pf_it->numberOfHits()               );
    pfcands_pixelLayersWithMeasurement -> push_back( pf_it->pixelLayersWithMeasurement() );
    pfcands_stripLayersWithMeasurement -> push_back( pf_it->stripLayersWithMeasurement() );
    pfcands_trackHighPurity            -> push_back( pf_it->trackHighPurity()            );

  }//loop over candidate collection

  for (auto ibegin = pfCandidates->begin(), pf_it1 = ibegin; pf_it1 != pfCandidates->end(); pf_it1++) {
    if (pf_it1->p4().pt() < minPt_ && abs(pf_it1->pdgId()) != 11 && abs(pf_it1->pdgId()) != 13) continue;
    // before Isolation are included in miniAOD, do all the iso calculations in 1 loop here
    float absIso = 0.0;
    float miniTkIso = 0.0;
    float miniDR = (pf_it1->pt() > 50)? (pf_it1->pt() > 200)? 0.05 : 10./pf_it1->pt() : 0.2; // compact version of miniDR calculation
    float baseDR = max(0.3F, miniDR);
    float thisPhi = (*pfcands_helperPhi)[pf_it1-ibegin];
    float thisEta = (*pfcands_helperEta)[pf_it1-ibegin];
    for (auto pf_it2 = ibegin; pf_it2 != pfCandidates->end(); pf_it2++) {
      if (pf_it2 == pf_it1) continue;
      // if (pf_it2->charge() == 0) continue; // skip neutrals
      if (abs(pf_it2->pdgId()) != 211) continue; // not considering leptons
      float deltaPhi = fabs((*pfcands_helperPhi)[pf_it2-ibegin] - thisPhi);
      if (deltaPhi > M_PI) deltaPhi = fabs(deltaPhi - 2.0*M_PI);
      if (deltaPhi > baseDR) continue;
      float deltaEta = fabs((*pfcands_helperEta)[pf_it2-ibegin] - thisEta);
      if (deltaEta > baseDR) continue;
      float dr = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
      if (dr > baseDR) continue;

      // use dz acceptance at 0.1, and accept if from PV
      if (pf_it2->pt() >= 0.0 && (*pfcands_helperPVinfo)[pf_it2-ibegin]) {
        // float dr = ROOT::Math::VectorUtil::DeltaR( pf_it2->p4(), pf_it1->p4() );
        if (dr < 0.3) absIso += pf_it2->pt();            // normal trackIso with cone size 0.3
        if (dr < miniDR) miniTkIso += pf_it2->pt();      // miniTrackIso
      }
    }
    pfcands_trackIso     -> push_back( absIso    );
    pfcands_miniTrackIso -> push_back( miniTkIso );
  }

  //define the phi bins
  vector<float> phibins;
  for (int i=0;i<10;i++) phibins.push_back(-TMath::Pi()+(2*i+1)*TMath::TwoPi()/20.);

  //define the eta bins
  vector<float> etabins_ctr;
  for (int i=0;i<8;++i) etabins_ctr.push_back(-2.1+0.6*i);
  vector<float> etabins_fwd;
  for (int i=0;i<10;++i) {
    if (i<5) etabins_fwd.push_back(-5.1+0.6*i);
    else etabins_fwd.push_back(2.7+0.6*(i-5));
  }
  vector<float> etabins_all;
  for (int i=0;i<18;++i) etabins_all.push_back(-5.1+0.6*i);

  //compute it
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  //Keep it
  iEvent.put(std::move(pfcands_p4                        ), branchprefix+"p4"                         );
  iEvent.put(std::move(pfcands_mass                      ), branchprefix+"mass"                       );
  iEvent.put(std::move(pfcands_dz                        ), branchprefix+"dz"                         );
  iEvent.put(std::move(pfcands_dxy                       ), branchprefix+"dxy"                        );
  iEvent.put(std::move(pfcands_dzError                   ), branchprefix+"dzError"                    );
  iEvent.put(std::move(pfcands_dxyError                  ), branchprefix+"dxyError"                   );
  iEvent.put(std::move(pfcands_charge                    ), branchprefix+"charge"                     );
  iEvent.put(std::move(pfcands_particleId                ), branchprefix+"particleId"                 );
  iEvent.put(std::move(pfcands_isGlobalMuon              ), branchprefix+"isGlobalMuon"               );
  iEvent.put(std::move(pfcands_isStandAloneMuon          ), branchprefix+"isStandAloneMuon"           );
  iEvent.put(std::move(pfcands_fromPV                    ), branchprefix+"fromPV"                     );
  iEvent.put(std::move(pfcands_pvAssociationQuality      ), branchprefix+"pvAssociationQuality"       );
  iEvent.put(std::move(pfcands_IdAssociatedPV            ), branchprefix+"IdAssociatedPV"             );
  iEvent.put(std::move(pfcands_dzAssociatedPV            ), branchprefix+"dzAssociatedPV"             );
  iEvent.put(std::move(pfcands_puppiWeight               ), branchprefix+"puppiWeight"                );
  iEvent.put(std::move(pfcands_puppiWeightNoLep          ), branchprefix+"puppiWeightNoLep"           );
  iEvent.put(std::move(pfcands_trackIso                  ), branchprefix+"trackIso"                   );
  iEvent.put(std::move(pfcands_miniTrackIso              ), branchprefix+"miniTrackIso"               );
  iEvent.put(std::move(pfcands_numberOfPixelHits         ), branchprefix+"numberOfPixelHits"          );
  iEvent.put(std::move(pfcands_numberOfHits              ), branchprefix+"numberOfHits"               );
  iEvent.put(std::move(pfcands_pixelLayersWithMeasurement), branchprefix+"pixelLayersWithMeasurement" );
  iEvent.put(std::move(pfcands_stripLayersWithMeasurement), branchprefix+"stripLayersWithMeasurement" );
  iEvent.put(std::move(pfcands_trackHighPurity           ), branchprefix+"trackHighPurity"            );

}

float PFCandidateMaker::getFixGridRho(std::vector<float>& etabins,std::vector<float>& phibins) {
  float etadist = etabins[1]-etabins[0];
  float phidist = phibins[1]-phibins[0];
  float etahalfdist = (etabins[1]-etabins[0])/2.;
  float phihalfdist = (phibins[1]-phibins[0])/2.;
  vector<float> sumPFNallSMDQ;
  sumPFNallSMDQ.reserve(etabins.size()*phibins.size());
  for (unsigned int ieta=0;ieta<etabins.size();++ieta) {
    for (unsigned int iphi=0;iphi<phibins.size();++iphi) {
      float pfniso_ieta_iphi = 0;
      for (pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++) {
        if (fabs(etabins[ieta]-pf_it->eta())>etahalfdist) continue;
        if (fabs(reco::deltaPhi(phibins[iphi],pf_it->phi()))>phihalfdist) continue;
        pfniso_ieta_iphi+=pf_it->pt();
      }
      sumPFNallSMDQ.push_back(pfniso_ieta_iphi);
    }
  }
  float evt_smdq = 0;
  sort(sumPFNallSMDQ.begin(),sumPFNallSMDQ.end());
  if (sumPFNallSMDQ.size()%2) evt_smdq = sumPFNallSMDQ[(sumPFNallSMDQ.size()-1)/2];
  else evt_smdq = (sumPFNallSMDQ[sumPFNallSMDQ.size()/2]+sumPFNallSMDQ[(sumPFNallSMDQ.size()-2)/2])/2.;
  return evt_smdq/(etadist*phidist);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFCandidateMaker);
