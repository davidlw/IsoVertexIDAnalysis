import FWCore.ParameterSet.Config as cms

isolatedOfflinePrimaryVertices = cms.EDProducer('IsoVertexIDProducer',

  TrackCollection = cms.InputTag('generalTracks'),
  VertexCollection = cms.InputTag('offlinePrimaryVertices'),

  sigmaZ = cms.double(0.08),
  nSigmaZ = cms.double(4.5),
  fContaminationMin = cms.double(0.01)
)
