#include "CMS2/NtupleMaker/interface/MCUtilities.h"
#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace std;


MCUtilities::MCUtilities() {
}

MCUtilities::~MCUtilities() {
}

const reco::GenParticle* MCUtilities::motherID(const reco::GenParticle& gp) {
  //extract motherid from a GenParticle by walking on the left side of the mother graph
  //inlining it here to avoid cyclic depce
    
  //yuck; this is all because mother access in the candidate is not virtual
  const reco::GenParticle* mom = &gp;
  while( mom->numberOfMothers() > 0 ) {
    for(uint j = 0; j < mom->numberOfMothers(); ++j) {

      mom = dynamic_cast<const reco::GenParticle*>( mom->mother(j) );

      if( mom->pdgId()!=gp.pdgId() )
	return mom;
    }
  }

  return mom;
}

void MCUtilities::writeDaughter( const reco::GenParticle& gp, int idx, vector<int>& genps_ld_id,
				 vector<int>& genps_ld_idx, vector<LorentzVector>& genps_ld_p4) {
  //call this for the status 3 particles to add all of their status 1 (not 2) daughters ( and grand daughters and great grand daughters ... )

  for( unsigned int i = 0; i < gp.numberOfDaughters(); i++ ) {

    if( gp.daughter(i)->status() == 1 ) {

      genps_ld_id.push_back( gp.daughter(i)->pdgId() );
      genps_ld_idx.push_back( idx                     );
      genps_ld_p4.push_back( LorentzVector(gp.daughter(i)->p4().px(),
					     gp.daughter(i)->p4().py(),
					     gp.daughter(i)->p4().pz(),
					     gp.daughter(i)->p4().e() ) );
    }
    else { 
      const reco::GenParticle* dau = dynamic_cast<const reco::GenParticle*>( gp.daughter(i) );
      MCUtilities::writeDaughter( *dau, idx, genps_ld_id, genps_ld_idx, genps_ld_p4 );
    }
  }
}

