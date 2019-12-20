import FWCore.ParameterSet.Config as cms

process = cms.Process("vertex")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hlt = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hlt.HLTPaths = ['HLT_PAFullTracks_Multiplicity185*_v*'] # for allphysics
process.hlt.andOr = cms.bool(True)
process.hlt.throw = cms.bool(False)

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)
    
#Reject beam scraping events standard pp configuration
process.NoScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.PAcollisionEventSelection = cms.Sequence(
                                         process.hfCoincFilter *
                                         process.PAprimaryVertexFilter *
                                         process.NoScraping
                                         )

process.eventFilter_HM = cms.Sequence(
#    process.hlt *
#    process.PAcollisionEventSelection
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16DR80/MinBias_TuneCUETP8M1_13TeV-pythia8/AODSIM/NoPU_80X_mcRun2_asymptotic_v14-v1/90000/F6F8F73B-C773-E611-994B-FA163EF5B328.root'
                )
#                                secondaryFileNames = cms.untracked.vstring('')
                            )
process.load("IsoVertexIDAnalysis.IsoVertexIDAnalyzer.vertextree_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('vertextree.root')
                                   )

process.ana = cms.Path(process.eventFilter_HM * process.vertextree_ana)
