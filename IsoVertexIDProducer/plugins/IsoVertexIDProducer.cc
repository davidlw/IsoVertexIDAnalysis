// -*- C++ -*-
//
// Package:    IsoVertexIDAnalysis/IsoVertexIDProducer
// Class:      IsoVertexIDProducer
// 
/**\class IsoVertexIDProducer IsoVertexIDProducer.cc IsoVertexIDAnalysis/IsoVertexIDProducer/plugins/IsoVertexIDProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wei Li
//         Created:  Thu, 19 Dec 2019 17:41:14 GMT
//
//


// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "Math/ProbFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
//
// class declaration
//

class IsoVertexIDProducer : public edm::stream::EDProducer<> {
   public:
      explicit IsoVertexIDProducer(const edm::ParameterSet&);
      ~IsoVertexIDProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:

      reco::VertexCollection theVtxs;

      float sigmaZ_;
      float nSigmaZ_;
      float fContaminationMin_;

      edm::EDGetTokenT<reco::TrackCollection> token_tracks;
      edm::EDGetTokenT<reco::VertexCollection> token_vertices;

      float fContamination(reco::Vertex vtx2, reco::Vertex vtx);
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

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
IsoVertexIDProducer::IsoVertexIDProducer(const edm::ParameterSet& iConfig) :
sigmaZ_(iConfig.getParameter<double>("sigmaZ")),
nSigmaZ_(iConfig.getParameter<double>("nSigmaZ")),
fContaminationMin_(iConfig.getParameter<double>("fContaminationMin"))
{
  token_vertices = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VertexCollection"));
  token_tracks = consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("TrackCollection"));

  produces< reco::VertexCollection >("");

   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


IsoVertexIDProducer::~IsoVertexIDProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
IsoVertexIDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle< reco::VertexCollection > vertices;
   iEvent.getByToken(token_vertices, vertices);
   if(!vertices->size()) { std::cout<<"Invalid or empty vertex collection!"<<std::endl; return; }
/*
   edm::Handle< reco::TrackCollection > tracks;
   iEvent.getByToken(token_tracks, tracks);
   if(!tracks->size()) { std::cout<<"Invalid or empty track collection!"<<std::endl; return; }
*/
   for(unsigned int iv=0; iv<vertices->size(); ++iv)
   {
     const reco::Vertex & vtx = (*vertices)[iv];
     if(vtx.isFake() || vtx.tracksSize()<2) continue;

     float fcont_sum = 0; 
     for(unsigned int jv=0; jv<vertices->size(); ++jv)
     {
       if(iv == jv) continue;

       const reco::Vertex & vtx2 = (*vertices)[jv];
       if(vtx2.isFake() || vtx2.tracksSize()<2) continue;

       fcont_sum += fContamination(vtx2,vtx); 
     }       
     if(fcont_sum > fContaminationMin_) continue;

     theVtxs.push_back( vtx );
   }

   std::unique_ptr< reco::VertexCollection >
     vertexCandidates( new reco::VertexCollection );
   vertexCandidates->reserve( theVtxs.size() );

   std::copy( theVtxs.begin(),
              theVtxs.end(),
              std::back_inserter(*vertexCandidates) );

   iEvent.put( std::move(vertexCandidates), std::string("") );

   theVtxs.clear();

/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::unique_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(std::move(pOut));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

float
IsoVertexIDProducer::fContamination(reco::Vertex vtx2, reco::Vertex vtx)
{
  float zvtx2 = vtx2.z();
  float ntrk2 = vtx2.tracksSize();

  float zvtx = vtx.z();
  float ntrk = vtx.tracksSize();

  float deltaZ = zvtx2 - zvtx;
  float fcont = ntrk2 * ROOT::Math::gaussian_cdf(-fabs(deltaZ)/sigmaZ_+nSigmaZ_) / ntrk;

  return fcont;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
IsoVertexIDProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
IsoVertexIDProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
IsoVertexIDProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
IsoVertexIDProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
IsoVertexIDProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
IsoVertexIDProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
IsoVertexIDProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoVertexIDProducer);
