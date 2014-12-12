//  //-*- C++ -*-
//  //
//  // Package:    NtupleMaker
//  // Class:      MITConversionMaker
//  // 
//  /**\class MITConversionMaker MITConversionMaker.cc CMS2/NtupleMaker/src/MITConversionMaker.cc
//  
//     Description: <one line class summary>
//  
//     Implementation:
//     <Notes on implementation>
//  */
//  //
//  // Original Author:  Puneeth Kalavase
//  //         Created:  Tues Sep  1 11:07:38 CDT 2009
//  // $Id: MITConversionMaker.cc,v 1.4 2011/02/17 23:12:23 dbarge Exp $
//  //
//  
//  // system include files
//  #include <memory>
//  
//  // user include files
//  #include "FWCore/Framework/interface/Frameworkfwd.h"
//  #include "FWCore/Framework/interface/EDProducer.h"
//  #include "FWCore/Framework/interface/Event.h"
//  #include "FWCore/Framework/interface/MakerMacros.h"
//  #include "FWCore/ParameterSet/interface/ParameterSet.h"
//  
//  #include "CMS2/NtupleMaker/interface/MITConversionMaker.h" 
//  #include "DataFormats/Common/interface/Handle.h"
//  #include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//  #include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//  #include "DataFormats/TrackReco/interface/Track.h"
//  #include "DataFormats/TrackReco/interface/TrackFwd.h"
//  #include "MitEdm/ConversionRejection/interface/ConversionMatcher.h"
//  #include "MitEdm/DataFormats/interface/StablePart.h"
//  #include "MitEdm/DataFormats/interface/StableData.h"
//  
//  #include "DataFormats/Math/interface/LorentzVector.h"
//  #include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
//  #include "DataFormats/Common/interface/RefToPtr.h"
//  
//  #include "CMS2/NtupleMaker/interface/MITConversionMaker.h" 
//  #include "DataFormats/Common/interface/Handle.h"
//  #include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//  #include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//  #include "DataFormats/TrackReco/interface/Track.h"
//  #include "DataFormats/TrackReco/interface/TrackFwd.h"
//  #include "MitEdm/ConversionRejection/interface/ConversionMatcher.h"
//  #include "MitEdm/DataFormats/interface/StablePart.h"
//  #include "DataFormats/Math/interface/Point3D.h"
//  #include "DataFormats/BeamSpot/interface/BeamSpot.h"
//  
//  
//  #include <TMath.h>
//  #include <Math/VectorUtil.h>
//  typedef math::XYZTLorentzVectorF LorentzVector;
//  
//  //
//  // constructors and destructor
//  //
//  
//  
//  bool matchesConversion(const reco::TrackRef& tk, const mitedm::StablePart *sp) {
//    
//    
//    if ( sp->trackPtr() == edm::Ptr<reco::Track>(refToPtr(tk)))
//      return true;
//    
//    return false;
//  }
//  
//  
//  bool matchesConversion(const reco::GsfElectron &ele, const mitedm::StablePart *sp) {
//    
//    
//    if ( sp->trackPtr() == edm::Ptr<reco::Track>(refToPtr(ele.gsfTrack())) ) 
//      return true;
//    
//    return false;
//  }
//  
//  
//  MITConversionMaker::MITConversionMaker(const edm::ParameterSet& iConfig) {
//  
//    using namespace edm;
//    using namespace std;
//  
//    elsInputTag_		= iConfig.getParameter<InputTag>("elsInputTag");
//    mitConversionsTag_	= iConfig.getParameter<InputTag>("mitConversionsInputTag");
//    ctfTrksInputTag_	= iConfig.getParameter<InputTag>("ctfTrksInputTag");
//    beamSpotTag_          = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
//  
//    produces<vector<LorentzVector> >		("mitconvposp4"			).setBranchAlias("mitconv_pos_p4"		);
//    produces<vector<float> >			("mitconvlxy"			).setBranchAlias("mitconv_lxy"			);
//    produces<vector<float> >			("mitconvlz"			).setBranchAlias("mitconv_lz"			);
//    produces<vector<double> >			("mitconvprob"			).setBranchAlias("mitconv_prob"			);
//    produces<vector<float> >			("mitconvndof"			).setBranchAlias("mitconv_ndof"			);
//    produces<vector<float> >			("mitconvchi2"			).setBranchAlias("mitconv_chi2"			);
//    produces<vector<int> >			("mitconvisGoodConversion"	).setBranchAlias("mitconv_isGoodConversion"	);
//    produces<vector<vector<int> > >		("mitconvtkalgo"                ).setBranchAlias("mitconv_tkalgo"               );
//    produces<vector<vector<int> > >		("mitconvtkidx"                 ).setBranchAlias("mitconv_tkidx"                );  
//    produces<vector<vector<LorentzVector> > >	("mitconvtkp4"                  ).setBranchAlias("mitconv_tkp4"                 );
//    produces<vector<vector<float> > >		("mitconvtkd0corr"              ).setBranchAlias("mitconv_tkd0corr"             );
//    produces<vector<vector<float> > >		("mitconvtkd0"                  ).setBranchAlias("mitconv_tkd0"                 );
//    produces<vector<vector<int> > >               ("mitconvtkq"                   ).setBranchAlias("mitconv_tkq"                  ); 
//      
//  }
//  
//  MITConversionMaker::~MITConversionMaker()
//  {
//  }
//  
//  void  MITConversionMaker::beginJob()
//  {
//  }
//  
//  void MITConversionMaker::endJob()
//  {
//  }
//  
//  // ------------ method called to produce the data  ------------
//  void MITConversionMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
//  
//    using namespace std;
//    using namespace mitedm;
//    using namespace edm;
//    using namespace reco;
//  
//    auto_ptr<vector<LorentzVector> >		mitconv_pos_p4			(new vector<LorentzVector>		);
//    auto_ptr<vector<float> >			mitconv_lxy			(new vector<float>			);
//    auto_ptr<vector<float> >			mitconv_lz			(new vector<float>			);
//    auto_ptr<vector<double> >			mitconv_prob			(new vector<double>			);
//    auto_ptr<vector<float> >			mitconv_ndof			(new vector<float >			);
//    auto_ptr<vector<float> >			mitconv_chi2			(new vector<float>			);
//    auto_ptr<vector<int> >			mitconv_isGoodConversion	(new vector<int>			);
//    auto_ptr<vector<vector<int> > >		mitconv_tkalgo                  (new vector<vector<int> >		);
//    auto_ptr<vector<vector<int> > >		mitconv_tkidx                   (new vector<vector<int> >		);
//    auto_ptr<vector<vector<LorentzVector> > >	mitconv_tkp4			(new vector<vector<LorentzVector> >	);
//    auto_ptr<vector<vector<float> > >		mitconv_tkd0corr		(new vector<vector<float> >		);
//    auto_ptr<vector<vector<float> > >		mitconv_tkd0			(new vector<vector<float> >		);
//    auto_ptr<vector<vector<int> > >		mitconv_tkq			(new vector<vector<int> > 		);
//  
//  
//    
//    //get electron collection
//    Handle<reco::GsfElectronCollection> hElectrons;
//    iEvent.getByLabel(elsInputTag_, hElectrons);
//    const reco::GsfElectronCollection *v_electronCol = hElectrons.product();
//    
//    //get collection of reconstructed conversions
//    edm::Handle<std::vector<mitedm::DecayPart> > hConversions;
//    iEvent.getByLabel(mitConversionsTag_, hConversions);
//    const vector<mitedm::DecayPart> *v_mitconversions = hConversions.product();
//  
//    edm::Handle<TrackCollection> hTracks;
//    iEvent.getByLabel(ctfTrksInputTag_, hTracks);
//    const vector<Track> *v_trks = hTracks.product();
//  
//    Handle<BeamSpot> beamSpotH;
//    iEvent.getByLabel(beamSpotTag_, beamSpotH);
//  
//    for(vector<mitedm::DecayPart>::const_iterator it_conv = v_mitconversions->begin();
//        it_conv != v_mitconversions->end(); it_conv++) {
//  
//      //initialize with all the default cuts, except the r of conv cut which has been relaxed to 10 cm
//      mitedm::ConversionMatcher convMatcher(-10, 0.0, 0.0, 1e-6, 1e-6,0); 
//  
//      mitconv_pos_p4->push_back(LorentzVector(it_conv->position().x(), it_conv->position().y(), it_conv->position().z(), 0.0));
//      mitconv_lxy->push_back(it_conv->lxy());
//      mitconv_lz->push_back(it_conv->lz());
//      mitconv_prob->push_back(TMath::Prob(it_conv->chi2(),it_conv->ndof()));    
//      mitconv_ndof->push_back(it_conv->ndof());
//      mitconv_chi2->push_back(it_conv->chi2());
//      mitconv_isGoodConversion->push_back(convMatcher.isGoodConversion(*it_conv));
//  
//      vector<int> vtemp_tkalgo;
//      vector<int> vtemp_tkidx;    
//      vector<LorentzVector> vtemp_tkp4;    
//      vector<float> vtemp_tkd0corr;
//      vector<float> vtemp_tkd0;
//      vector<int> vtemp_tkq;
//  
//      //now look for the tracks belonging to this vertex
//      for(int i = 0; i < it_conv->nStableChild(); i++) {
//        const StableData &sd = it_conv->getStableData(i); 
//        const StablePart *sp = dynamic_cast<const StablePart*>(sd.originalPtr().get());
//        const reco::Track *trk = sp->track();
//  
//        double corrd0 = beamSpotH.isValid() ? -1 * ( trk->dxy(beamSpotH->position()) ) : trk->d0();		           
//        vtemp_tkalgo.push_back(trk->algo());
//        vtemp_tkp4.push_back(LorentzVector(trk->px(), trk->py(), trk->pz(), trk->p()));
//        vtemp_tkd0corr.push_back( corrd0                                                   );
//        vtemp_tkd0.push_back(trk->d0());
//        vtemp_tkq.push_back(trk->charge());
//  
//  
//        bool foundtrack = false;
//        for(vector<GsfElectron>::const_iterator el_it = v_electronCol->begin();
//  	  el_it != v_electronCol->end(); el_it++) {
//  	if(matchesConversion(*el_it, sp)) {	    
//  	  if(trk->algo() != 29)
//  	    throw cms::Exception("MITConversionMaker::produce(): The candidate track matches an electron track, but the track algo does not match! Needs Investigation!");
//  	  vtemp_tkidx.push_back(el_it - v_electronCol->begin());
//  	  foundtrack = true;
//  	  break;
//  	}
//        }
//        if(foundtrack)
//  	continue;
//        
//        for(vector<Track>::const_iterator tk_it = v_trks->begin();
//  	  tk_it != v_trks->end(); tk_it++) {
//  	TrackRef tempRef = TrackRef(hTracks, tk_it - v_trks->begin());
//  	if(matchesConversion(tempRef, sp) ) {	    
//  	  if(!(trk->algo() > 0 && trk->algo() < 15)) 
//  	    throw cms::Exception("MITConversionMaker::produce(): The candidate track matches a ctf track, but the track algo does not match! Needs Investigation!");
//  	  vtemp_tkidx.push_back(tk_it - v_trks->begin());
//  	  foundtrack = true;
//  	  break;
//  	}
//        }
//        
//        if(foundtrack)
//  	continue;
//        
//        vtemp_tkidx.push_back(-9999);
//      }
//      
//      mitconv_tkidx	->push_back(vtemp_tkidx);
//      mitconv_tkalgo	->push_back(vtemp_tkalgo);
//      mitconv_tkp4	->push_back(vtemp_tkp4);
//      mitconv_tkd0corr	->push_back(vtemp_tkd0corr);
//      mitconv_tkd0	->push_back(vtemp_tkd0);
//      mitconv_tkq		->push_back(vtemp_tkq);
//      
//    }//loop over conversions
//    
//  
//    iEvent.put(mitconv_pos_p4,		"mitconvposp4"			);
//    iEvent.put(mitconv_lxy,		"mitconvlxy"			);
//    iEvent.put(mitconv_lz,		"mitconvlz"			);
//    iEvent.put(mitconv_prob,		"mitconvprob"			);
//    iEvent.put(mitconv_ndof,		"mitconvndof"			);
//    iEvent.put(mitconv_chi2,		"mitconvchi2"			);
//    iEvent.put(mitconv_isGoodConversion,	"mitconvisGoodConversion"	);
//    iEvent.put(mitconv_tkalgo,		"mitconvtkalgo"			);
//    iEvent.put(mitconv_tkidx,		"mitconvtkidx"			);
//    iEvent.put(mitconv_tkp4,		"mitconvtkp4"			);
//    iEvent.put(mitconv_tkd0corr,		"mitconvtkd0corr"		);
//    iEvent.put(mitconv_tkd0,		"mitconvtkd0"			);
//    iEvent.put(mitconv_tkq,		"mitconvtkq"			); 
//    
//    
//  
//    
//  }
//  
//  
//  
//  
//  
//  //define this as a plug-in
//  DEFINE_FWK_MODULE(MITConversionMaker);
