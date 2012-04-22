#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CMS2/NtupleMaker/interface/MVAJetIdMaker.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "RecoMET/METAlgorithms/interface/PFSpecificAlgo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Utilities/interface/Exception.h"

typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
MVAJetIdMaker::MVAJetIdMaker(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;

  // product of this EDProducer
  produces<vector<float> >         ( "pfjetsmvavalue" ).setBranchAlias( "pfjets_mvavalue"	);
  
  //
  fVertexNameTag_   = iConfig.getParameter<InputTag>	( "VertexName" 		);
  fCorrJetName    	= iConfig.getParameter<InputTag>	( "CorrJetName"		);
  fUnCorrJetName  	= iConfig.getParameter<InputTag>	( "JetName"			);

  // 
  fPUJetIdAlgo    	= new PileupJetIdAlgo(iConfig); 

}

// Destructor
MVAJetIdMaker::~MVAJetIdMaker(){}

// ------------ method called once each job just before starting event loop  ------------
void MVAJetIdMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void MVAJetIdMaker::endJob() {}

// ------------ method called to produce the data  ------------
void MVAJetIdMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;
 
  // create containers
  auto_ptr<vector<float> >         pfjets_mvavalue                (new vector<float>          );  


  //Uncorrected Jets
  Handle<PFJetCollection>       lHUCJets;
  iEvent.getByLabel(fUnCorrJetName, lHUCJets);
  PFJetCollection               lUCJets = *lHUCJets;

  //Corrected Jets
  Handle<PFJetCollection>       lHCJets;
  iEvent.getByLabel(fCorrJetName  , lHCJets);
  PFJetCollection               lCJets = *lHCJets;

  // vertices    
  Handle<reco::VertexCollection> lHVertices;
  iEvent.getByLabel(fVertexNameTag_      , lHVertices); 
  VertexCollection lVertices = *lHVertices;
  Vertex *lPV = 0; if(lVertices.size() > 0) lPV = &lVertices[0]; 


  for(int i0   = 0; i0 < (int) lUCJets.size(); i0++) {   // uncorrecte jets collection                                           
	  const PFJet       *pUCJet = &(lUCJets.at(i0));
	  for(int i1 = 0; i1 < (int) lCJets .size(); i1++) {   // corrected jets collection                                         

		  const PFJet     *pCJet  = &(lCJets.at(i1));
		  if(       pUCJet->jetArea() != pCJet->jetArea()                  ) continue;
		  if( fabs(pUCJet->eta() - pCJet->eta())         > 0.01            ) continue;
		  double lJec = pCJet ->pt()/pUCJet->pt();
		  PileupJetIdentifier lPUJetId =  fPUJetIdAlgo->computeIdVariables(pCJet,lJec,&lVertices[0],lVertices,true);
	
		  //
		  // fill branch
		  //
		  pfjets_mvavalue              	->push_back( lPUJetId.mva()              );
	
		  // debug
		  if(false)
		  {
			  std::cout << "Debug Jet MVA: "
				  << lPUJetId.nvtx()      << " "
				  << lPUJetId.jetPt()     << " "
				  << lPUJetId.jetEta()     << " "
				  << lPUJetId.jetPhi()     << " "
				  << lPUJetId.d0()       << " "
				  << lPUJetId.dZ()      << " "
				  << lPUJetId.beta()      << " "
				  << lPUJetId.betaStar()  << " "
				  << lPUJetId.nCharged()  << " "
				  << lPUJetId.nNeutrals() << " "
				  << lPUJetId.dRMean()    << " "
				  << lPUJetId.frac01()    << " "
				  << lPUJetId.frac02()    << " "
				  << lPUJetId.frac03()    << " "
				  << lPUJetId.frac04()    << " "
				  << lPUJetId.frac05()
				  << " === : === "
				  << lPUJetId.mva() << " " << endl;
		  }

		  break;
	  }
  }



  // 
  // 
  // 
  iEvent.put(pfjets_mvavalue,             	 "pfjetsmvavalue"                 );
}



//define this as a plug-in
DEFINE_FWK_MODULE(MVAJetIdMaker);
