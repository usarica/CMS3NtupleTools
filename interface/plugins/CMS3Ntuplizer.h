#include <cassert>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include <memory>

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/Run.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <FWCore/ServiceRegistry/interface/Service.h>

#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include <SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

#include <CMSDataTools/AnalysisTree/interface/SimpleEntry.h>
#include <CMSDataTools/AnalysisTree/interface/BaseTree.h>


class CMS3Ntuplizer : public edm::EDAnalyzer{
public:
  explicit CMS3Ntuplizer(const edm::ParameterSet&);
  ~CMS3Ntuplizer();

protected:
  const edm::ParameterSet pset;
  BaseTree* outtree;

  TString treename;
  //TString outfilename;


private:
  virtual void beginJob();
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);

  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

};
