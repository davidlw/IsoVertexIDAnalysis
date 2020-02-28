// -*- C++ -*-
//
// Package:    IsoVertexIDAnalysis/IsoVertexIDTreeAnalyzer
// Class:      IsoVertexIDTreeAnalyzer
// 
/**\class IsoVertexIDTreeAnalyzer IsoVertexIDTreeAnalyzer.cc IsoVertexIDAnalysis/IsoVertexIDTreeAnalyzer/plugins/IsoVertexIDTreeAnalyzer.cc

 Description: Vertex tree producer

 Implementation:
     
*/
//
// Original Author:  Wei Li
//         Created:  Thu, 19 Dec 2019 17:42:08 GMT
//
//


// system include files
#include <memory>
#include <exception>
#include <iostream>
#include <memory>
#include <algorithm>
#include <math.h>
#include <vector>
#include <string>
#include <map>

// user include files
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <TTree.h>

//
// class declaration
//

class TTree;

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

#define NMAXVTX 100
#define NMAXTRACKSVTX 300

class IsoVertexIDTreeAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit IsoVertexIDTreeAnalyzer(const edm::ParameterSet&);
      ~IsoVertexIDTreeAnalyzer();

      long trkIdMaker( float pt, float eta, float phi, float vz, float chi2 );

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:

      edm::Service<TFileService> theOutputs;
      edm::EDGetTokenT<reco::TrackCollection> token_tracks;
      edm::EDGetTokenT<edm::View<reco::Track> > token_edmViewTracks;
      edm::EDGetTokenT<reco::VertexCollection> token_vertices;
//    edm::EDGetTokenT<std::vector< PileupSummaryInfo > > token_genPU;   
      edm::EDGetTokenT< std::vector< TrackingVertex> >token_TrackingVtx;
      edm::EDGetTokenT< std::vector< TrackingParticle> >token_TrackingParticle;

      //trking particle association maps
      edm::EDGetTokenT<reco::RecoToSimCollection> associatorMapRTS_;
      edm::EDGetTokenT<reco::SimToRecoCollection> associatorMapSTR_;


      TTree* vertexTree;

      uint runNb;
      uint eventNb;
      uint lsNb;

      bool isSlim_;
      bool isMC_;

      uint nVertices;
      uint nTracks[NMAXVTX];
      uint nTracks_ptGT0p3EtaLT2p4[NMAXVTX];

      float xVtx[NMAXVTX];
      float yVtx[NMAXVTX];
      float zVtx[NMAXVTX];
      float xVtxErr[NMAXVTX];
      float yVtxErr[NMAXVTX];
      float zVtxErr[NMAXVTX];
      float chi2Vtx[NMAXVTX];
      int   ndofVtx[NMAXVTX];
      float minZSepReco[NMAXVTX];
      int   bestMatchGenByWeight[NMAXVTX];
      int   bestMatchGenByWeightMult[NMAXVTX];
      float contaminationFracByWeight[NMAXVTX];
      int   bestMatchGenByTracks[NMAXVTX];
      int   bestMatchGenByTracksMult[NMAXVTX];
      float contaminationFracByTracks[NMAXVTX];
      float trackWeightVtx[NMAXVTX][NMAXTRACKSVTX];
      float trackPtVtx[NMAXVTX][NMAXTRACKSVTX];
      float trackEtaVtx[NMAXVTX][NMAXTRACKSVTX];
      float trackPhiVtx[NMAXVTX][NMAXTRACKSVTX];
      float trackPtErrVtx[NMAXVTX][NMAXTRACKSVTX];
      float trackXVtx[NMAXVTX][NMAXTRACKSVTX];
      float trackYVtx[NMAXVTX][NMAXTRACKSVTX];
      float trackZVtx[NMAXVTX][NMAXTRACKSVTX];
      float trackXYErrVtx[NMAXVTX][NMAXTRACKSVTX];
      float trackZErrVtx[NMAXVTX][NMAXTRACKSVTX];
      bool  trackHPVtx[NMAXVTX][NMAXTRACKSVTX];
      int  trackGenVertex[NMAXVTX][NMAXTRACKSVTX];

      int nVertices_gen;
      std::vector< int > eventIDList;
      std::vector< float > xVtx_gen;
      std::vector< float > yVtx_gen;
      std::vector< float > zVtx_gen;
      std::vector< int > nTracks_gen;
      std::vector< float > minZSepGen;

      void resetArrays();
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
IsoVertexIDTreeAnalyzer::IsoVertexIDTreeAnalyzer(const edm::ParameterSet& iConfig):
isSlim_(iConfig.getParameter<bool>("isSlim")),
isMC_(iConfig.getParameter<bool>("isMC"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

   token_vertices = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VertexCollection"));
   token_tracks = consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("TrackCollection"));

   edm::InputTag TrackingVertexSrc_("mix","MergedTrackTruth");
   edm::InputTag TrackingParticleSrc_("mix","MergedTrackTruth");
   if(isMC_){
     token_TrackingVtx = consumes< std::vector< TrackingVertex > >(TrackingVertexSrc_);
     token_TrackingParticle = consumes< std::vector< TrackingParticle > >(TrackingParticleSrc_);
     token_edmViewTracks = consumes<edm::View<reco::Track>>(iConfig.getParameter<edm::InputTag>("TrackCollection"));
     associatorMapRTS_ = consumes< reco::RecoToSimCollection >(iConfig.getParameter<edm::InputTag>("associatorMap"));
     associatorMapSTR_ = consumes< reco::SimToRecoCollection >(iConfig.getParameter<edm::InputTag>("associatorMap"));
   }

//   edm::InputTag PileupSrc_("addPileupInfo");
//   if(isMC_) token_genPU = consumes< std::vector< PileupSummaryInfo > >(PileupSrc_);

}


IsoVertexIDTreeAnalyzer::~IsoVertexIDTreeAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
IsoVertexIDTreeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   resetArrays();

   using namespace edm;

   edm::Handle< reco::VertexCollection > vertices;
   iEvent.getByToken(token_vertices, vertices);
   if(!vertices->size()) { std::cout<<"Invalid or empty vertex collection!"<<std::endl; vertexTree->Fill(); return; }

   edm::Handle< reco::TrackCollection > tracks;
   iEvent.getByToken(token_tracks, tracks);
   if(!tracks->size()) { std::cout<<"Invalid or empty track collection!"<<std::endl; return; }

   runNb = iEvent.id().run();
   eventNb = iEvent.id().event();
   lsNb = iEvent.luminosityBlock();
 
   //get the tracking Vertex info (truth for tracking purposes)
   if(isMC_){
     edm::Handle< std::vector< TrackingVertex > > trkVtx;
     iEvent.getByToken(token_TrackingVtx, trkVtx);

     edm::Handle< std::vector< TrackingParticle > > trkPart;
     iEvent.getByToken(token_TrackingParticle, trkPart);

     std::vector< TrackingVertex >::const_iterator tv;
     for( tv = trkVtx->begin(); tv != trkVtx->end(); ++tv){
       if( tv->eventId().bunchCrossing()!=0) continue; //skip tracking vertices that are not in the main bunch crossing (reject Out-of-time PU)
       if( tv->nSourceTracks() != 0 ) continue;//Skip vertices that decay from other stuff (not PVs)
       if( fabs( tv->position().z() ) > 30) continue;// skip stuff >30cm in each direction

       //make sure this is from a prompt process (not a g4 decay)
       std::vector<SimVertex>::const_iterator g4Iter;
       int processTotal = 0; 
       for( g4Iter = tv->g4Vertices_begin(); g4Iter != tv->g4Vertices_end(); ++g4Iter){
         processTotal += g4Iter->processType();//processes >0 are not prompt processes, so we don't want them
       }
       if(processTotal > 0) continue;
       
       //std::cout << tv->eventId().event() << " " <<  tv->daughterTracks().size() << std::endl;
       
       //get the gen multiplicity
       int mult = 0;
       TrackingParticleRefVector::iterator daughters;
       for( daughters = tv->daughterTracks_begin(); daughters != tv->daughterTracks_end(); ++daughters){
         if((*daughters)->pt() > 0.3 && TMath::Abs((*daughters)->eta()) < 2.4 && (*daughters)->charge()!=0) mult++;
       }

       //take the first tracking vertex from the event to be it's position
       //Note: this may be slightly inaccurate if a PV splits into multiple tracking vertices, but it looks like the differences are <0.1cm in z on average
       std::vector<int>::iterator it = std::find(eventIDList.begin(), eventIDList.end(), tv->eventId().event());    
       if (it == eventIDList.end()){
         eventIDList.push_back(tv->eventId().event());
         xVtx_gen.push_back(tv->position().x());
         yVtx_gen.push_back(tv->position().y());
         zVtx_gen.push_back(tv->position().z());
         nTracks_gen.push_back(mult);
       } else {
         nTracks_gen.at( std::distance( eventIDList.begin(), it) ) += mult;//add these tracks to a event we already found
       }
     }
     nVertices_gen = eventIDList.size(); 
     for(size_t i = 0; i<zVtx_gen.size(); i++){
       float minSepZ = 999;
       for(size_t j = 0; j<zVtx_gen.size(); j++){
         if(i==j) continue;
         float sep = fabs( zVtx_gen.at(i) - zVtx_gen.at(j) );
         if( sep < minSepZ) minSepZ = sep; 
       }
       minZSepGen.push_back(minSepZ);
     }
   }
   
   //track association maps (need to fix scoping here)
   //This is ugly, but I cannot figure out how to use the RecoToSimCollections below in the vetex loop, so I build my own map...
   std::map<long, int> trk2EvtMap; 
   if(isMC_){ 
     reco::RecoToSimCollection recSimColl;
     reco::SimToRecoCollection simRecColl;
     edm::Handle<reco::SimToRecoCollection > simtorecoCollectionH;
     edm::Handle<reco::RecoToSimCollection > recotosimCollectionH;
     iEvent.getByToken(associatorMapSTR_,simtorecoCollectionH);
     simRecColl= *(simtorecoCollectionH.product());
     iEvent.getByToken(associatorMapRTS_,recotosimCollectionH);
     recSimColl= *(recotosimCollectionH.product()); 

     Handle<edm::View<reco::Track> > tcol;
     iEvent.getByToken(token_edmViewTracks, tcol);
     for(edm::View<reco::Track>::size_type i=0; i<tcol->size(); ++i){
       edm::RefToBase<reco::Track> trk(tcol, i);

       std::vector<std::pair<TrackingParticleRef, double> > tp;
       const TrackingParticle *mtp=0;
       if(recSimColl.find(trk) != recSimColl.end())
       {
         tp = recSimColl[trk];
         mtp = tp.begin()->first.get();
         if(mtp->eventId().bunchCrossing()!=0) continue;
         trk2EvtMap.insert(std::pair<long, int>( trkIdMaker(trk->pt(), trk->eta(), trk->phi(), trk->vz(), trk->chi2() ) , mtp->eventId().event())  );
         //std::cout << mtp->eventId().event() << std::endl;
       }
     }
   }

   /*if(isMC_){
     edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
     iEvent.getByToken(token_genPU, PupInfo);

     std::vector<PileupSummaryInfo>::const_iterator PVI;
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
       if(PVI->getBunchCrossing() == 0){//we only want the in-time bunch crossing
         meanPU_gen = PVI->getTrueNumInteractions();
         nVertices_gen = PVI->getPU_NumInteractions();
         for(size_t i = 0; i<PVI->getPU_zpositions().size() ; i++ ) zVtx_gen.push_back( (PVI->getPU_zpositions()).at(i) );
       }
       //std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
     }
   }*/

   nVertices=0;

//std::cout<<vertices->size()<<std::endl;

   for(unsigned int iv=0; iv<vertices->size(); iv++)
   {
     const reco::Vertex & vtx = (*vertices)[iv];

     if(!vtx.isFake() && vtx.tracksSize()>=2)
     {
       xVtx[nVertices] = vtx.x();
       yVtx[nVertices] = vtx.y();
       zVtx[nVertices] = vtx.z();
       xVtxErr[nVertices] = vtx.xError();
       yVtxErr[nVertices] = vtx.yError();
       zVtxErr[nVertices] = vtx.zError();
       chi2Vtx[nVertices] = vtx.chi2();
       ndofVtx[nVertices] = vtx.ndof();
      
       for(unsigned int iv2=0; iv2<vertices->size(); iv2++)
       {
         if(iv==iv2) continue;
         const reco::Vertex & vtx2 = (*vertices)[iv2];
         if(!vtx2.isFake() && vtx2.tracksSize()>=2)
         {
           float zSep = fabs(vtx.z()-vtx2.z());
           if(zSep < minZSepReco[nVertices]) minZSepReco[nVertices] = zSep;
         }
       }

       float weightNet = 0;
       float trackNet = 0;
       std::vector< float > weightTotals;
       std::vector< float > trackTotals;
       if(isMC_){
         for(size_t i = 0; i<eventIDList.size(); i++) weightTotals.push_back(0);
         for(size_t i = 0; i<eventIDList.size(); i++) trackTotals.push_back(0);
       }

       uint nTracksTmp=0;
       for (reco::Vertex::trackRef_iterator iTrack = vtx.tracks_begin(); iTrack != vtx.tracks_end(); iTrack++) {

          reco::TrackRef track = iTrack->castTo<reco::TrackRef>();

          trackWeightVtx[nVertices][nTracksTmp] = vtx.trackWeight(*iTrack);
          trackPtVtx[nVertices][nTracksTmp] = track->pt();
          trackEtaVtx[nVertices][nTracksTmp] = track->eta();
          trackPhiVtx[nVertices][nTracksTmp] = track->phi();
          trackPtErrVtx[nVertices][nTracksTmp] = track->ptError();
          trackXVtx[nVertices][nTracksTmp] = track->vx();
          trackYVtx[nVertices][nTracksTmp] = track->vy();
          trackZVtx[nVertices][nTracksTmp] = track->vz();
          trackXYErrVtx[nVertices][nTracksTmp] = track->d0Error();
          trackZErrVtx[nVertices][nTracksTmp] = track->dzError();
          trackHPVtx[nVertices][nTracksTmp] = track->quality(reco::TrackBase::highPurity);
//std::cout<<nVertices<<" "<<zVtx[nVertices]<<" "<<nTracksTmp<<" "<<trackPtVtx[nVertices][nTracksTmp]<<" "<<trackZVtx[nVertices][nTracksTmp]<<std::endl;

          if(isMC_){
            int eventNumber = -1;
            long ID = trkIdMaker(track->pt(), track->eta(), track->phi(), track->vz(), track->chi2() );
            std::map<long,int>::iterator isFound = trk2EvtMap.find(ID);
            if(isFound != trk2EvtMap.end() ){
              eventNumber = (isFound)->second;
              
              //std::cout << "Found Number " << eventNumber << std::endl;
              //std::cout << "Distance " << std::distance( eventIDList.begin(), std::find( eventIDList.begin(), eventIDList.end(), eventNumber) ) << std::endl;
              int index = std::distance( eventIDList.begin(), std::find( eventIDList.begin(), eventIDList.end(), eventNumber) );
              trackGenVertex[nVertices][nTracksTmp] = index;
              
              weightNet += vtx.trackWeight(*iTrack);
              trackNet ++; 
              weightTotals.at(index) += vtx.trackWeight(*iTrack);
              trackTotals.at(index)++;
            }
          }

          if(trackHPVtx[nVertices][nTracksTmp] && trackPtVtx[nVertices][nTracksTmp]>0.3 && TMath::Abs(trackEtaVtx[nVertices][nTracksTmp]) < 2.4) nTracks_ptGT0p3EtaLT2p4[nVertices]++;

          if(isSlim_)
          {
            if(!trackHPVtx[nVertices][nTracksTmp]) continue;
            if(trackPtErrVtx[nVertices][nTracksTmp]/trackPtVtx[nVertices][nTracksTmp]>0.1) continue;
            if(fabs(trackZVtx[nVertices][nTracksTmp]-zVtx[nVertices])/pow(trackZErrVtx[nVertices][nTracksTmp]*trackZErrVtx[nVertices][nTracksTmp]+zVtxErr[nVertices]*zVtxErr[nVertices],0.5)>3) continue;
            if(pow(pow(trackXVtx[nVertices][nTracksTmp]-xVtx[nVertices],2)+pow(trackYVtx[nVertices][nTracksTmp]-yVtx[nVertices],2),0.5)/pow(trackXYErrVtx[nVertices][nTracksTmp]*trackXYErrVtx[nVertices][nTracksTmp]+xVtxErr[nVertices]*yVtxErr[nVertices],0.5)>3) continue;             
          }
          nTracksTmp++;
       }

       if(isMC_){
         for(size_t i = 0; i<eventIDList.size(); i++) weightTotals.at(i) = weightTotals.at(i)/weightNet;
         for(size_t i = 0; i<eventIDList.size(); i++) trackTotals.at(i) = trackTotals.at(i)/trackNet;
         std::vector< float >::iterator maxEle = std::max_element(weightTotals.begin(), weightTotals.end());
         bestMatchGenByWeight[nVertices] = std::distance( weightTotals.begin(), maxEle);
         if(bestMatchGenByWeight[nVertices] >=0) bestMatchGenByWeightMult[nVertices] = nTracks_gen[bestMatchGenByWeight[nVertices]];
         contaminationFracByWeight[nVertices] = 1 - *maxEle;
         std::vector< float >::iterator maxEleTrk = std::max_element(trackTotals.begin(), trackTotals.end());
         bestMatchGenByTracks[nVertices] = std::distance( trackTotals.begin(), maxEleTrk);
         if(bestMatchGenByTracks[nVertices] >=0) bestMatchGenByTracksMult[nVertices] = nTracks_gen[bestMatchGenByTracks[nVertices]];
         contaminationFracByTracks[nVertices] = 1 - *maxEleTrk;
       }

       nTracks[nVertices] = nTracksTmp;

       nVertices++;
     }
   }
/*
   for(uint i=0;i<nVertices;i++)
   {
     for(uint j=0;j<nTracks[i];j++)
     {
std::cout<<"After: "<<i<<" "<<zVtx[i]<<" "<<j<<" "<<trackPtVtx[i][j]<<" "<<trackZVtx[i][j]<<std::endl;

     }
   }
*/
   vertexTree->Fill();
}

void
IsoVertexIDTreeAnalyzer::resetArrays()
{
  for(int i=0;i<NMAXVTX;i++)
  {
    xVtx[i] = -999.0;
    yVtx[i] = -999.0;
    zVtx[i] = -999.0;
    xVtxErr[i] = -999.0;
    yVtxErr[i] = -999.0;
    zVtxErr[i] = -999.0;
    chi2Vtx[i] = -999.0;
    ndofVtx[i] = -999.0;
    nTracks[i] = -999;
    nTracks_ptGT0p3EtaLT2p4[i] = 0;
    if(isMC_){
      contaminationFracByWeight[i] = -1;   
      contaminationFracByTracks[i] = -1;   
      bestMatchGenByWeight[i] = -1; 
      bestMatchGenByWeightMult[i] = -1; 
      bestMatchGenByTracks[i] = -1; 
      bestMatchGenByTracksMult[i] = -1; 
    }
    minZSepReco[i] = 999.0;

    for(int j=0;j<NMAXTRACKSVTX;j++)
    {
      trackWeightVtx[i][j] = -999.0;
      trackPtVtx[i][j] = -999.0;
      trackEtaVtx[i][j] = -999.0;
      trackPhiVtx[i][j] = -999.0;
      trackPtErrVtx[i][j] = -999.0;
      trackXVtx[i][j] = -999.0;
      trackYVtx[i][j] = -999.0;
      trackZVtx[i][j] = -999.0;
      trackXYErrVtx[i][j] = -999.0;
      trackZErrVtx[i][j] = -999.0;
      trackHPVtx[i][j] = -999.0;
      if(isMC_) trackGenVertex[i][j] = -1;
    }
  }
  if(isMC_){
    nVertices_gen = -999;
    eventIDList.clear();
    xVtx_gen.clear();
    yVtx_gen.clear();
    zVtx_gen.clear();
    nTracks_gen.clear();
    minZSepGen.clear();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
IsoVertexIDTreeAnalyzer::beginJob()
{
  vertexTree = theOutputs->make<TTree>("vertexTree","vertexTree");
  vertexTree->Branch("RunNb",&runNb,"RunNb/i");
  vertexTree->Branch("LSNb",&lsNb,"LSNb/i");
  vertexTree->Branch("EventNb",&eventNb,"EventNb/i");
  vertexTree->Branch("nVertices",&nVertices,"nVertices/i");
  vertexTree->Branch("xVtx",xVtx,"xVtx[nVertices]/F");
  vertexTree->Branch("yVtx",yVtx,"yVtx[nVertices]/F");
  vertexTree->Branch("zVtx",zVtx,"zVtx[nVertices]/F");
  vertexTree->Branch("nTracks",nTracks,"nTracks[nVertices]/i");
  vertexTree->Branch("nTracks_ptGT0p3EtaLT2p4",nTracks_ptGT0p3EtaLT2p4,"nTracks_ptGT0p3EtaLT2p4[nVertices]/i");
  vertexTree->Branch("trackPtVtx",trackPtVtx,"trackPtVtx[nVertices][300]/F");
  vertexTree->Branch("trackEtaVtx",trackEtaVtx,"trackEtaVtx[nVertices][300]/F");
  vertexTree->Branch("trackPhiVtx",trackPhiVtx,"trackPhiVtx[nVertices][300]/F");
  vertexTree->Branch("minZSepReco",minZSepReco,"minZSepReco[nVertices]/F");

  if(!isSlim_)
  {
    vertexTree->Branch("xVtxErr",xVtxErr,"xVtxErr[nVertices]/F");
    vertexTree->Branch("yVtxErr",yVtxErr,"yVtxErr[nVertices]/F");
    vertexTree->Branch("zVtxErr",zVtxErr,"zVtxErr[nVertices]/F");
    vertexTree->Branch("chi2Vtx",chi2Vtx,"chi2Vtx[nVertices]/F");
    vertexTree->Branch("ndofVtx",ndofVtx,"ndofVtx[nVertices]/i");
    vertexTree->Branch("trackWeightVtx",trackWeightVtx,"trackWeightVtx[nVertices][300]/F");
    vertexTree->Branch("trackPtErrVtx",trackPtErrVtx,"trackPtErrVtx[nVertices][300]/F");
    vertexTree->Branch("trackXVtx",trackXVtx,"trackXVtx[nVertices][300]/F");
    vertexTree->Branch("trackYVtx",trackYVtx,"trackYVtx[nVertices][300]/F");
    vertexTree->Branch("trackZVtx",trackZVtx,"trackZVtx[nVertices][300]/F");
    vertexTree->Branch("trackXYErrVtx",trackYVtx,"trackXYErrVtx[nVertices][300]/F");
    vertexTree->Branch("trackZErrVtx",trackZVtx,"trackZErrVtx[nVertices][300]/F");
    vertexTree->Branch("trackHPVtx",trackHPVtx,"trackHPVtx[nVertices][300]/O");
  }
  if(isMC_){
    vertexTree->Branch("bestMatchGenByWeight",bestMatchGenByWeight,"bestMatchGenByWeight[nVertices]/I");
    vertexTree->Branch("bestMatchGenByWeightMult",bestMatchGenByWeightMult,"bestMatchGenByWeightMult[nVertices]/I");
    vertexTree->Branch("contaminationFracByWeight",contaminationFracByWeight,"contaminationFracByWeight[nVertices]/F");
    vertexTree->Branch("bestMatchGenByTracks",bestMatchGenByTracks,"bestMatchGenByTracks[nVertices]/I");
    vertexTree->Branch("bestMatchGenByTracksMult",bestMatchGenByTracksMult,"bestMatchGenByTracksMult[nVertices]/I");
    vertexTree->Branch("contaminationFracByTracks",contaminationFracByTracks,"contaminationFracByTracks[nVertices]/F");
    vertexTree->Branch("trackGenVertex",trackGenVertex,"trackGenVertex[nVertices][300]/I");
    vertexTree->Branch("nVertices_gen",&nVertices_gen,"nVertices_gen/I");
    vertexTree->Branch("xVtx_gen",&xVtx_gen);
    vertexTree->Branch("yVtx_gen",&yVtx_gen);
    vertexTree->Branch("zVtx_gen",&zVtx_gen);
    vertexTree->Branch("nTracks_ptGT0p3EtaLT2p4_gen",&nTracks_gen);
    vertexTree->Branch("minZSepGen",&minZSepGen);
  }
}

long IsoVertexIDTreeAnalyzer::trkIdMaker( float pt, float eta, float phi, float vz, float chi2 ){
  long ID = (int) (pt*100000) + (int)(fabs(eta)*10000) + (int)(fabs(phi)*1000) + (int)(fabs(vz)*100 + (int)chi2);
  return ID;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
IsoVertexIDTreeAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
IsoVertexIDTreeAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoVertexIDTreeAnalyzer);
