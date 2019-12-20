import FWCore.ParameterSet.Config as cms

isolatedOfflinePrimaryVertices = cms.EDProducer('IsoVertexIDProducer',

  TrackCollection = cms.InputTag('generalTracks'),
  VertexCollection = cms.InputTag('offlinePrimaryVertices'),

  sigmaZ = cms.double(0.1),
  nSigmaZ = cms.double(3),
  fContaminationMin = cms.double(0.01)
)
