#include "CMS2/NtupleMaker/interface/ElUtilities.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"

using namespace std;
using namespace reco;

ElUtilities::ElUtilities() {
}

ElUtilities::~ElUtilities() {
}


//--------------------------------------------------------------------------------
//get the electron collection
//--------------------------------------------------------------------------------
vector<const GsfElectron*> ElUtilities::getElectrons(const edm::Event& iEvent, 
						     const edm::InputTag electronsInputTag) {

  edm::Handle<edm::View<GsfElectron> > electron_h;
  iEvent.getByLabel(electronsInputTag, electron_h);
  vector<const GsfElectron*> collection;

  for(edm::View<reco::GsfElectron>::const_iterator electron = electron_h->begin();
      electron != electron_h->end(); ++electron){
    collection.push_back(&*electron);
  }

  return collection;

}


//--------------------------------------------------------------------------------------------
//remove electrons that have the same SC
//--------------------------------------------------------------------------------------------
void ElUtilities::removeElectrons(const vector<const GsfElectron*>* collection) {

  vector<const reco::GsfElectron*>* newEl = const_cast<vector<const GsfElectron*>* >(collection);
  vector<const reco::GsfElectron*> copy = *newEl;

  newEl->clear();

  vector<const GsfElectron*>::iterator it1, it2;

  for(it1=copy.begin(); it1!=copy.end(); ++it1) {

    bool isRemoved = false;
    for(it2=copy.begin(); it2!=copy.end(); ++it2) {
      if (it1 == it2)
        continue;
      if (((**it1).superCluster().id() == (**it2).superCluster().id()) &&
          ((**it1).superCluster().index() == (**it2).superCluster().index())) {
	cout << __FILE__ << endl;
        float deltaEp1 = fabs((**it1).eSuperClusterOverP() - 1.);
        float deltaEp2 = fabs((**it2).eSuperClusterOverP() - 1.);
        if (deltaEp1 > deltaEp2) {
          isRemoved = true;
          break;
        }
      }
    }

    if (!isRemoved)
      newEl->push_back(*it1);

  }
}
//--------------------------------------------------------------------------------
