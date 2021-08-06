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
<use   name="CMS3/Dictionaries"/>
<use   name="EgammaAnalysis/ElectronTools"/>
<use   name="PhysicsTools/NanoAOD"/>
<use   name="NNKit/FatJetNN"/>
<use   name="JHUGenMELA/MELA"/>
<use   name="CommonLHETools/LHEHandler"/>
<use   name="CMSDataTools/AnalysisTree"/>
<use   name="MelaAnalytics/EventContainer"/>
<use   name="MelaAnalytics/GenericMEComputer"/>
<use   name="MelaAnalytics/CandidateLOCaster"/>
<use   name="root"/>
<use   name="clhep"/>
<use   name="CMS3/NtupleMaker"/>


<library   file="*.cc" name="CMS3NtupleMakerPlugins">
  <Flags CPPDEFINES="CMSSW_VERSION=.oOCMSSW_VERSIONOo."/>
  <Flags CPPDEFINES="CMSSW_VERSION_MAJOR=.oOCMSSW_VERSION_MAJOROo."/>
  <Flags CPPDEFINES="CMSSW_VERSION_MINOR=.oOCMSSW_VERSION_MINOROo."/>
  <flags CPPFLAGS="-I$(CMSSW_BASE)/src/JHUGenMELA/MELA/interface/" />
  <flags CXXFLAGS="-g -O3"/>
  <flags EDM_PLUGIN="1"/>
</library>
