import FWCore.ParameterSet.Config as cms

vertextree_ana = cms.EDAnalyzer('IsoVertexIDTreeAnalyzer',

  TrackCollection = cms.InputTag('generalTracks'),
  VertexCollection = cms.InputTag('offlinePrimaryVertices')

)
