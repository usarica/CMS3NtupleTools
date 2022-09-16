<use   name="FWCore/Framework"/>
<use   name="FWCore/PluginManager"/>
<use   name="FWCore/ParameterSet"/>
<use   name="FWCore/Utilities"/>
<use   name="CondFormats/JetMETObjects"/>
<use   name="DataFormats/TauReco"/>
<use   name="DataFormats/MuonReco"/>
<use   name="DataFormats/JetReco"/>
<use   name="DataFormats/GsfTrackReco"/>
<use   name="DataFormats/EgammaReco"/>
<use   name="DataFormats/EgammaCandidates"/>
<use   name="DataFormats/EcalDetId"/>
<use   name="DataFormats/TrackReco"/>
<use   name="DataFormats/Common"/>
<use   name="DataFormats/Math"/>
<use   name="DataFormats/PatCandidates"/>
<use   name="DataFormats/METReco"/>
<use   name="Geometry/Records"/>
<use   name="Geometry/CommonDetUnit"/>
<use   name="Geometry/TrackerGeometryBuilder"/>
<use   name="RecoParticleFlow/PFProducer"/>
<use   name="CommonTools/ParticleFlow"/>
<use   name="RecoEcal/EgammaCoreTools"/>
<use   name="RecoEgamma/EgammaTools"/>
<use   name="JetMETCorrections/Objects"/>
<use   name="JetMETCorrections/Modules"/>
<use   name="HLTrigger/HLTcore"/>
<use   name="EgammaAnalysis/ElectronTools"/>
<use   name="PhysicsTools/NanoAOD"/>
<use   name="NNKit/FatJetNN"/>
<use   name="CommonLHETools/LHEHandler"/>
<use   name="IvyFramework/IvyDataTools"/>
<use   name="IvyFramework/IvyAutoMELA"/>
<use   name="root"/>
<use   name="clhep"/>
<use   name="CMS3/Dictionaries"/>
<use   name="CMS3/NtupleMaker"/>


<Flags CXXDEFINES="CMSSW_VERSION=.oOCMSSW_VERSIONOo."/>
<Flags CXXDEFINES="CMSSW_VERSION_MAJOR=.oOCMSSW_VERSION_MAJOROo."/>
<Flags CXXDEFINES="CMSSW_VERSION_MINOR=.oOCMSSW_VERSION_MINOROo."/>
<flags CXXFLAGS="-g -O3 -I$(CMSSW_BASE)/src/IvyFramework/IvyDataTools/interface/ -I$(CMSSW_BASE)/src/IvyFramework/IvyAutoMELA/interface/ -I$(MELA_LIB_PATH)/../../interface/"/>
<flags LDFLAGS="-L$(MELA_LIB_PATH) -L$(MELAANALYTICS_PATH)/EventContainer/lib -L$(MELAANALYTICS_PATH)/GenericMEComputer/lib -L$(MELAANALYTICS_PATH)/CandidateLOCaster/lib -lmcfm_707 -ljhugenmela -lcollier -lJHUGenMELAMELA -lMelaAnalyticsEventContainer -lMelaAnalyticsGenericMEComputer -lMelaAnalyticsCandidateLOCaster" />


<library   file="*.cc" name="CMS3NtupleMakerPlugins">
  <flags EDM_PLUGIN="1"/>
</library>
