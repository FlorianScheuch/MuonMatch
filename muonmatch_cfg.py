import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")




process.MessageLogger = cms.Service("MessageLogger",
	destinations = cms.untracked.vstring("Log")
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
fileList = cms.untracked.vstring()
fileList.extend(['/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0017F55A-EDF8-E211-B97E-003048C692DE.root'])


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = fileList
)

process.load('PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi')
process.muonMatch.src = cms.InputTag("muons")
process.muonMatch.checkCharge = cms.bool(False)
process.muonMatch.resolveAmbiguities = cms.bool(False)
		

#process.selectedMuonsGenParticlesMatchNew = cms.EDProducer("MCTruthDeltaRMatcherNew",
#                                              src = cms.InputTag("muons"),
#                                              matched = cms.InputTag("genParticles"),
#                                              distMin = cms.double(0.15),
#                                              matchPDGId = cms.vint32(13)
#)

process.demo = cms.EDAnalyzer('MuonMatch'
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('AnalysisResult.root')
)

process.p = cms.Path(process.muonMatch*process.demo)
