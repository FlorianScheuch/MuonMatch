import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonCollector")

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:00D4272D-8EF1-E111-BE0F-0017A4770C18.root'
	#'root://xrootd.unl.edu/
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0017F55A-EDF8-E211-B97E-003048C692DE.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/003697A5-03FA-E211-A6C3-0025907DCAA6.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0042F21A-8CF9-E211-BCA5-0025907FD428.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/00794212-CCF9-E211-B7EE-0030487E55C5.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/026499D5-CEF9-E211-B046-0030487E5101.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0445D1E4-2CF9-E211-93B4-0030487EBB2B.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/046A8F09-EEF8-E211-87E5-003048C69408.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/06579C7B-20F9-E211-8D77-0030487E4EB9.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/06CF10CB-34F9-E211-98FF-0030487E5101.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0887DC15-96F9-E211-9B58-0025907DC9BA.root'
    )
)

process.load('PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi')
process.muonMatch.checkCharge = cms.bool(False)
process.muonMatch.resolveAmbiguities = cms.bool(False)

process.load('PhysicsTools.PatAlgos.producersLayer1.muonProducer_cff')

process.out = cms.OutputModule(
	"PoolOutputModule",
	outputCommands = cms.untracked.vstring(
		'drop *',
		'keep *_*muon*_*_*',
		'keep *_*genParticles*_*_*',
		'keep *_*muonMatch*_*_*'),
	fileName = cms.untracked.string('muonCollections.root')
	)

process.outpath = cms.EndPath(process.muonMatch*
	process.makePatMuons*
	process.out)