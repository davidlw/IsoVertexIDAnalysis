import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000)
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

process.eventFilter = cms.Sequence(
#    process.hlt *
#    process.PAcollisionEventSelection
)
process.eventFilter_step = cms.Path( process.eventFilter )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#'root://cmsxrootd.fnal.gov//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/471/00000/02778130-2ABD-E611-8640-02163E014467.root'
#'root://cmsxrootd.fnal.gov//store/data/Run2016G/ZeroBias/AOD/PromptReco-v1/000/280/330/00000/44EE54CF-8577-E611-BB1E-02163E0143EA.root'
#'root://cmsxrootd.fnal.gov//store/data/Run2016G/Tau/AOD/07Aug17-v1/50002/0CEC9511-31A0-E711-A487-48D539F383F2.root'
'file:/storage1/users/wl33/0CEC9511-31A0-E711-A487-48D539F383F2.root'
                )
#                                secondaryFileNames = cms.untracked.vstring('')
                            )
process.load("IsoVertexIDAnalysis.IsoVertexIDProducer.isolatedOfflinePrimaryVertices_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('isovtxprod.root')
                                   )

process.ana = cms.Path(process.eventFilter * process.isolatedOfflinePrimaryVertices)

#process.output = cms.OutputModule("PoolOutputModule",
#    outputCommands = cms.untracked.vstring('drop *','keep *_isolatedOfflinePrimaryVertices_*_*','keep *_offlinePrimaryVertices_*_*','keep *_generalTracks_*_*'), 
#    fileName = cms.untracked.string('isovtxprod.root'),
#    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_step')),
#    dataset = cms.untracked.PSet(
#     dataTier = cms.untracked.string('AOD')
#    )
#)
#process.output_step = cms.EndPath(process.output)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('vertextree.root')
                                   )

process.load("IsoVertexIDAnalysis.IsoVertexIDAnalyzer.vertextree_cff")
process.vertextree_isolated_ana = process.vertextree_ana.clone()
process.vertextree_isolated_ana.VertexCollection = cms.InputTag('isolatedOfflinePrimaryVertices')

process.ana1 = cms.Path(process.eventFilter * process.vertextree_ana)
process.ana1_iso = cms.Path(process.eventFilter * process.vertextree_isolated_ana)
