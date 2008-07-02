#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>

MatchUtilities::MatchUtilities() {
}

MatchUtilities::~MatchUtilities() {
}


const reco::Candidate* MatchUtilities::matchGenToCand(const reco::GenJet& genJet, std::vector<const reco::Candidate*> cand) {

  const reco::Candidate* output = 0;
  double dRmin = 0.25;

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

const reco::GenJet* MatchUtilities::matchCandToGenJet(const reco::Candidate& jet, const std::vector<reco::GenJet>* genJets) { 
  
  const reco::GenJet* output = 0;
  double dRmin = 0.25;
  
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


const reco::GenParticle* MatchUtilities::matchCandToGen(const reco::Candidate& cand, const std::vector<reco::GenParticle>* genParticles) {

  const reco::GenParticle* output = 0;
  double dRmin = 0.1;
  
  std::vector<reco::GenParticle>::const_iterator itPartEnd = genParticles->end();
  for(std::vector<reco::GenParticle>::const_iterator itPart=genParticles->begin(); itPart!=itPartEnd; ++itPart) {

    if ( itPart->status() != 3 ) {

      const math::XYZVector v1(itPart->momentum().x(), itPart->momentum().y(), itPart->momentum().z());

      double dR = ROOT::Math::VectorUtil::DeltaR(v1,cand.p4());
      if (dR < dRmin) {
	dRmin = dR;
	output = &(*itPart);
      }
    }
  }

  return output;
}

const reco::GenParticle* MatchUtilities::matchCandToGen(const reco::Candidate& cand, const std::vector<reco::GenParticle>* genParticles, int& genidx) {

  const reco::GenParticle* output = 0;
  double dRmin = 0.1;
  unsigned int i = 0;

  std::vector<reco::GenParticle>::const_iterator itPartEnd = genParticles->end();
  for(std::vector<reco::GenParticle>::const_iterator itPart=genParticles->begin(); itPart!=itPartEnd; ++itPart, ++i) {

    if ( itPart->status() != 3 ) {

      const math::XYZVector v1(itPart->momentum().x(), itPart->momentum().y(), itPart->momentum().z());

      double dR = ROOT::Math::VectorUtil::DeltaR(v1,cand.p4());
      if (dR < dRmin) {
	dRmin = dR;
	output = &(*itPart);
	genidx = i;
      }
    }
  }

  return output;
}

const reco::GenParticle* MatchUtilities::matchCandToGen(const reco::Track& track, const std::vector<reco::GenParticle>* genParticles, int& genidx) {

  const reco::GenParticle* output = 0;
  double dRmin = 0.1;
  unsigned int i = 0;
  
  std::vector<reco::GenParticle>::const_iterator itPartEnd = genParticles->end();
  for(std::vector<reco::GenParticle>::const_iterator itPart=genParticles->begin(); itPart!=itPartEnd; ++itPart, ++i) {

    if ( itPart->status() != 3 ) {

      const math::XYZVector v1(itPart->momentum().x(), itPart->momentum().y(), itPart->momentum().z());

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > cand( track.px(), track.py(), track.pz(), track.p() );

      double dR = ROOT::Math::VectorUtil::DeltaR(v1,cand);

      if (dR < dRmin) {
	dRmin = dR;
	output = &(*itPart);
	genidx = i;
      }//END find minimum delta R loop
    }//END loop over status != 3 particles
  }//END loop over genParticles

  return output;
}

const reco::Candidate* MatchUtilities::matchGenToCand(const reco::GenParticle& p, std::vector<const reco::Candidate*> cand) {

  const reco::Candidate* output = 0;
  double dRmin = 0.1;

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
