import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

# configure modules via Global Tag
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START53_V19::All'



process.MessageLogger = cms.Service("MessageLogger",
	destinations = cms.untracked.vstring("Log")
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
fileList = cms.untracked.vstring()
#F2F1C38E-79DA-E111-B6E8-003048FFCB6A.root
fileList.extend(['file:/user/scheuch/CMSSW/crab/CMSSW_5_3_12_patch2/src/MuonMatch/MuonMatch/F2F1C38E-79DA-E111-B6E8-003048FFCB6A.root'])
#root://xrootd.unl.edu//store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/34675A33-99F9-E211-A1F4-0025907FD2BA.root

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = fileList
)

process.load("L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff")
process.es_prefer_l1GtTriggerMaskAlgoTrig = cms.ESPrefer("L1GtTriggerMaskAlgoTrigTrivialProducer", "l1GtTriggerMaskAlgoTrig")

process.load('HLTrigger.HLTfilters.hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(False)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('L1_DoubleMu3')

process.load('PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi')
process.muonMatch.src = cms.InputTag("muons")
process.muonMatch.mcPdgId = cms.vint32(13)
process.muonMatch.checkCharge = cms.bool(False)
process.muonMatch.resolveAmbiguities = cms.bool(False)
process.muonMatch.maxDeltaR = cms.double(10.)
process.muonMatch.maxDPtRel = cms.double(100.)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

process.load("L1TriggerConfig.L1GtConfigProducers.l1GtTriggerMenuXml_cfi")
process.l1GtTriggerMenuXml._errorstr = 0;

process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.l1Filter.algorithms = cms.vstring('L1_DoubleMu3')

#process.selectedMuonsGenParticlesMatchNew = cms.EDProducer("MCTruthDeltaRMatcherNew",
#                                              src = cms.InputTag("muons"),
#                                              matched = cms.InputTag("genParticles"),
#                                              distMin = cms.double(0.15),
#                                              matchPDGId = cms.vint32(13)
#)

process.demo = cms.EDAnalyzer('MuonMatch'
)

#process.demo2 = cms.EDAnalyzer('MuonMatchPre'
#)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('AnalysisResult.root')
)

process.p = cms.Path(process.muonMatch*process.demo) #process.l1GtTriggerMenuXml, process.l1Filter*  *process.hltLevel1GTSeed
