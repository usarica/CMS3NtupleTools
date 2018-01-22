#include "CMS3/NtupleMaker/interface/MatchUtilities.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Utilities/interface/Exception.h"

typedef math::XYZTLorentzVectorF LorentzVector;

MatchUtilities::MatchUtilities() {
}

MatchUtilities::~MatchUtilities() {
}



//----------------------------------------------------------------------------------------------
const reco::GenParticle* MatchUtilities::matchCandToGen(const reco::Candidate& cand, 
							const std::vector<reco::GenParticle>* genParticles, 
							int& genidx, int status, const std::vector<int> v_PIDsToExclude) {

  const reco::GenParticle* output = 0;
  double dRmin = 0.2;
  unsigned int i = 0;
  genidx = -9999;
  
  std::vector<reco::GenParticle>::const_iterator itPartEnd = genParticles->end();
  for(std::vector<reco::GenParticle>::const_iterator itPart=genParticles->begin(); itPart!=itPartEnd; ++itPart, ++i) {

    if ( status != 999 && itPart->status() != status ) continue;
    if ( find(v_PIDsToExclude.begin(), v_PIDsToExclude.end(), abs(itPart->pdgId()) ) != v_PIDsToExclude.end() ) 
      continue;

    const math::XYZVector v1(itPart->momentum().x(), itPart->momentum().y(), itPart->momentum().z());

    double dR = ROOT::Math::VectorUtil::DeltaR(v1,cand.p4());
    if (dR < dRmin) {
      dRmin = dR;
      output = &(*itPart);
      genidx = i;
    }
    
  }

  return output;
}

//----------------------------------------------------------------------------------------------
const reco::GenParticle* MatchUtilities::matchCandToGen(const reco::Track& track, 
							const std::vector<reco::GenParticle>* genParticles, 
							int& genidx, int status, const std::vector<int> v_PIDsToExclude) {

  const reco::GenParticle* output = 0;
  double dRmin = 0.2;
  unsigned int i = 0;
  genidx = -9999;
  
  std::vector<reco::GenParticle>::const_iterator itPartEnd = genParticles->end();
  for(std::vector<reco::GenParticle>::const_iterator itPart=genParticles->begin(); 
      itPart!=itPartEnd; ++itPart, ++i) {

    if ( status != 999 && itPart->status() != status ) continue;
    if ( find(v_PIDsToExclude.begin(), v_PIDsToExclude.end(), abs(itPart->pdgId()) ) != v_PIDsToExclude.end() ) 
      continue;

    const math::XYZVector v1(itPart->momentum().x(), itPart->momentum().y(), itPart->momentum().z());

    LorentzVector cand( track.px(), track.py(), track.pz(), track.p() );

    double dR = ROOT::Math::VectorUtil::DeltaR(v1,cand);

    if (dR < dRmin) {
      dRmin = dR;
      output = &(*itPart);
      genidx = i;
    }//END find minimum delta R loop
  
  }//END loop over genParticles

  return output;
}


//----------------------------------------------------------------------------------------------
const reco::GenParticle* MatchUtilities::matchCandToGen(const LorentzVector& candp4, 
							const std::vector<reco::GenParticle>* genParticles, 
							int& genidx, int status, const std::vector<int> v_PIDsToExclude) {

  const reco::GenParticle* output = 0;
  double dRmin = 0.2;
  unsigned int i = 0;
  genidx = -9999;
  
  std::vector<reco::GenParticle>::const_iterator itPartEnd = genParticles->end();
  for(std::vector<reco::GenParticle>::const_iterator itPart=genParticles->begin(); itPart!=itPartEnd; ++itPart, ++i) {

    if ( status != 999 && itPart->status() != status ) continue;
    if ( find(v_PIDsToExclude.begin(), v_PIDsToExclude.end(), abs(itPart->pdgId()) ) != v_PIDsToExclude.end() ) 
      continue;

    const math::XYZVector v1(itPart->momentum().x(), itPart->momentum().y(), itPart->momentum().z());

    double dR = ROOT::Math::VectorUtil::DeltaR(v1,candp4);

    if (dR < dRmin) {
      dRmin = dR;
      output = &(*itPart);
      genidx = i;
    }//END find minimum delta R loop
  
  }//END loop over genParticles

  return output;
}
 
//----------------------------------------------------------------------------------------------
// THIS IS THE BAD GUY
const pat::PackedGenParticle* MatchUtilities::matchCandToGen(const LorentzVector& candp4, 
							const std::vector<pat::PackedGenParticle>* genParticles, 
							int& genidx, int status, const std::vector<int> v_PIDsToExclude) {

  const pat::PackedGenParticle* output = 0;
  double dRmin = 0.2;
  unsigned int i = 0;
  genidx = -9999;

  double phi = candp4.Phi();
  double eta = candp4.Eta();
  double pi = M_PI;

  
  std::vector<pat::PackedGenParticle>::const_iterator itPartEnd = genParticles->end();
  for(std::vector<pat::PackedGenParticle>::const_iterator itPart=genParticles->begin(); itPart!=itPartEnd; ++itPart, ++i) {

    if ( status != 999 && itPart->status() != status ) continue;

    int id = abs(itPart->pdgId());
    if (id == 12) continue;
    if (id == 14) continue;
    if (id == 16) continue;
    if (id == 18) continue;
    if (id == 1000022) continue;

    double deltaPhi = phi-itPart->phi();
    if ( deltaPhi > pi ) deltaPhi -= 2.0*pi;
    else if ( deltaPhi <= -pi ) deltaPhi += 2.0*pi;
    deltaPhi = fabs(deltaPhi);
    if (deltaPhi > dRmin) continue;
    double deltaEta = fabs(itPart->eta()-eta);
    if (deltaEta > dRmin) continue;
    double dR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

    if (dR < dRmin) {
      dRmin = dR;
      output = &(*itPart);
      genidx = i;
    }//END find minimum delta R loop
  
  }//END loop over genParticles

  return output;
}

//----------------------------------------------------------------------------------------------
const reco::Candidate* MatchUtilities::matchGenToCand(const reco::GenJet& genJet,
						      std::vector<const reco::Candidate*> cand) {

  const reco::Candidate* output = 0;
  double dRmin = 0.2;

  std::vector<const reco::Candidate*>::const_iterator itCand;
  
  for(itCand=cand.begin(); itCand!=cand.end(); ++itCand) {

    const math::XYZVector v1(genJet.momentum().x(), genJet.momentum().y(), genJet.momentum().z());

    double dR = ROOT::Math::VectorUtil::DeltaR(v1, (*itCand)->p4());
    if (dR < dRmin) {
      dRmin = dR;
      output = *itCand;
    }
  }

  return output;
}


//----------------------------------------------------------------------------------------------
const reco::GenJet* MatchUtilities::matchCandToGenJet(const reco::Candidate& jet, 
						      const std::vector<reco::GenJet>* genJets) { 
  
  const reco::GenJet* output = 0;
  double dRmin = 0.3;
  
  std::vector<reco::GenJet>::const_iterator itJetEnd = genJets->end();
  for(std::vector<reco::GenJet>::const_iterator itJet=genJets->begin(); itJet!=itJetEnd; ++itJet) {

    const math::XYZVector v1(itJet->momentum().x(), itJet->momentum().y(), itJet->momentum().z());

    double dR = ROOT::Math::VectorUtil::DeltaR(v1, jet.p4());
    if (dR < dRmin) {
      dRmin = dR;
      output = &(*itJet);
    }
  }

  return output;
}
//----------------------------------------------------------------------------------------------
const reco::GenJet* MatchUtilities::matchCandToGenJet(const LorentzVector& jetp4, 
						      const std::vector<reco::GenJet>* genJets,
						      int &genidx) { 
  
  const reco::GenJet* output = 0;
  double dRmin = 0.3;
  int i = 0;
  genidx = -9999;
  
  std::vector<reco::GenJet>::const_iterator itJetEnd = genJets->end();
  for(std::vector<reco::GenJet>::const_iterator itJet=genJets->begin(); itJet!=itJetEnd; ++itJet, ++i) {

    const math::XYZVector v1(itJet->momentum().x(), itJet->momentum().y(), itJet->momentum().z());

    double dR = ROOT::Math::VectorUtil::DeltaR(v1, jetp4);
    if (dR < dRmin) {
      dRmin = dR;
      output = &(*itJet);
      genidx = i;
    }
  }

  return output;
}



//----------------------------------------------------------------------------------------------
const reco::Candidate* MatchUtilities::matchGenToCand(const reco::GenParticle& p, 
						      std::vector<const reco::Candidate*> cand) {

  const reco::Candidate* output = 0;
  double dRmin = 0.2;

  std::vector<const reco::Candidate*>::const_iterator itCand;
  
  for(itCand=cand.begin(); itCand!=cand.end(); ++itCand) {

    const math::XYZVector v1(p.momentum().x(), p.momentum().y(), p.momentum().z());

    double dR = ROOT::Math::VectorUtil::DeltaR(v1, (*itCand)->p4());
    if (dR < dRmin) {
      dRmin = dR;
      output = *itCand;
    }
  }

  return output;
}
//----------------------------------------------------------------------------------------------

const int MatchUtilities::getMatchedGenIndex(const reco::GenParticle& p, const std::vector<reco::GenParticle>* genParticles, 
					     int status, const std::vector<int> v_PIDsToExclude) {

  double dRmin = 0.2; 
  std::vector<reco::GenParticle>::const_iterator itCand;
  int idx = -9999;
  int temp = 0;
  math::XYZVector v1(p.momentum().x(), p.momentum().y(), p.momentum().z());
  for(itCand = genParticles->begin(); itCand != genParticles->end(); itCand++, temp++) {
    
    if(status != 999 && itCand->status() != status)
      continue;
    if ( find(v_PIDsToExclude.begin(), v_PIDsToExclude.end(), abs(itCand->pdgId()) ) != v_PIDsToExclude.end() ) 
      continue;

    double dR = ROOT::Math::VectorUtil::DeltaR(v1, itCand->p4());
    
    if(dR < dRmin) {
      idx = temp;
      dRmin = dR;
    }
  }

  return idx;
}


//----------------------------------------------------------------------------------------------

const void MatchUtilities::alignRecoPatJetCollections(const std::vector<reco::CaloJet>& v_ref,
						      std::vector<pat::Jet>& v_toAllign) {
  alignCollections(v_ref, v_toAllign);
}
//----------------------------------------------------------------------------------------------
const void MatchUtilities::alignRecoPatElectronCollections(const std::vector<reco::GsfElectron>& v_ref,
							   std::vector<pat::Electron>& v_toAllign) {
  alignCollections(v_ref, v_toAllign);
}
//----------------------------------------------------------------------------------------------
const void MatchUtilities::alignRecoPatMuonCollections(const std::vector<reco::Muon>& v_ref,
						       std::vector<pat::Muon>& v_toAllign) {
  alignCollections(v_ref, v_toAllign);
}
//----------------------------------------------------------------------------------------------
const void MatchUtilities::alignJPTcaloJetCollections(const std::vector<reco::CaloJet>& v_ref,
						      std::vector<reco::CaloJet>& v_toAllign) {
  alignCollections(v_ref, v_toAllign);
}

  
const std::tuple<float,float,int,float> MatchUtilities::getLepMVAInfo(edm::Ptr<reco::Candidate> lep, edm::Handle<edm::View<pat::Jet> > pfJetsHandle, const reco::Vertex &vtx) {
    // Taken from https://github.com/cms-sw/cmssw/blob/70922db3413f0934414865c37e735a97add7daef/PhysicsTools/NanoAOD/plugins/LeptonJetVarProducer.cc#L175
    int jetNDauChargedMVASel = 0;
    float ptratio = 1.0;
    float ptrel = 0.0;
    float disc = -99.; // -99 if not found: https://github.com/cms-sw/cmssw/blob/1fed0702c9db78ad3f8c2b1e53967c2b6c8ae939/PhysicsTools/NanoAOD/python/electrons_cff.py#L123
    for (unsigned int ij = 0; ij<pfJetsHandle->size(); ij++){
        auto jet = pfJetsHandle->ptrAt(ij);

        auto rawp4 = jet->correctedP4("Uncorrected");
        auto lepp4 = lep->p4();
        if ((rawp4-lepp4).R()<1e-4) break;

        if(matchByCommonSourceCandidatePtr(*lep,*jet)) {
            unsigned int jndau = 0;
            for(const auto _d : jet->daughterPtrVector()) {
                const auto d = dynamic_cast<const pat::PackedCandidate*>(_d.get());
                if (d->charge()==0) continue;
                if (d->fromPV()<=1) continue;
                if (reco::deltaR(*d,*lep)>0.4) continue;
                if (!(d->hasTrackDetails())) continue;
                auto tk = d->pseudoTrack();
                if(tk.pt()>1 &&
                        tk.hitPattern().numberOfValidHits()>=8 &&
                        tk.hitPattern().numberOfValidPixelHits()>=2 &&
                        tk.normalizedChi2()<5 &&
                        fabs(tk.dxy(vtx.position()))<0.2 &&
                        fabs(tk.dz(vtx.position()))<17
                  ) jndau++;
            }
            jetNDauChargedMVASel = jndau;
            auto jetp4 = (rawp4 - lepp4*(1.0/jet->jecFactor("L1FastJet")))*(jet->pt()/rawp4.pt())+lepp4;
            ptratio = lepp4.pt()/jetp4.pt();
            ptrel = lepp4.Vect().Cross((jetp4-lepp4).Vect().Unit()).R();
            disc = jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
            break; // take leading jet with shared source candidates
        }
    }
    return std::tuple<float,float,int,float>(ptratio,ptrel,jetNDauChargedMVASel,disc);
}
