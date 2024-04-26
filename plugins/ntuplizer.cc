#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"

float dxy_value(const reco::GenParticle &p, const reco::Vertex &pv){
    float vx = p.vx();
    float vy = p.vy();
    float phi = p.phi();
    float pv_x = pv.x();
    float pv_y = pv.y();
  
    float dxy = -(vx-pv_x)*sin(phi) + (vy-pv_y)*cos(phi);
    return dxy;
}

void print_trackerHits(const edm::Handle<std::vector<PSimHit>> trackerHits, const SimTrack &simTrack) {
  for (size_t i = 0; i < trackerHits->size(); i++) {
    const PSimHit &hit = (*trackerHits)[i];
    if (hit.trackId() == simTrack.trackId()) {
      std::cout << "        -> Hit from particle type: " << hit.particleType() << " || trackId = " << hit.trackId() 
                << "\n                                               localPosition = " << hit.localPosition() 
                << "\n                                               localDirection = " << hit.localDirection() 
                << "\n                                               momentumAtEntry = " << hit.momentumAtEntry() << std::endl;
    }
  }
}

class ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ntuplizer(const edm::ParameterSet&);
      ~ntuplizer();

      edm::ConsumesCollector iC = consumesCollector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::ParameterSet parameters;

      bool isData = true;
      //
      // --- Tokens and Handles
      //

      // trigger bits
      //edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      //edm::Handle<edm::TriggerResults> triggerBits;
    
      // GenParticles
      edm::EDGetTokenT<edm::View<reco::GenParticle> >  theGenParticleCollection;
      edm::Handle<edm::View<reco::GenParticle> > genParticles;

      // SimTracks
      edm::EDGetTokenT<std::vector<SimTrack> >  theSimTrackCollection;
      edm::Handle<std::vector<SimTrack> > simTracks;

      // TrackerHits
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsPixelBarrelHighTofCollection;
      edm::Handle<std::vector<PSimHit>> trackerHitsPixelBarrelHighTof;
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsPixelBarrelLowTofCollection;
      edm::Handle<std::vector<PSimHit>> trackerHitsPixelBarrelLowTof;
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsPixelEndcapHighTofCollection;
      edm::Handle<std::vector<PSimHit>> trackerHitsPixelEndcapHighTof;
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsPixelEndcapLowTofCollection;
      edm::Handle<std::vector<PSimHit>> trackerHitsPixelEndcapLowTof;

      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsTECHighTofCollection; 
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsTECLowTofCollection;
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsTIBHighTofCollection;
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsTIBLowTofCollection;
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsTIDHighTofCollection;
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsTIDLowTofCollection;
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsTOBHighTofCollection;
      edm::EDGetTokenT<std::vector<PSimHit>>  theTrackerHitsTOBLowTofCollection;
      edm::Handle<std::vector<PSimHit>>  trackerHitsTECHighTof; 
      edm::Handle<std::vector<PSimHit>>  trackerHitsTECLowTof;
      edm::Handle<std::vector<PSimHit>>  trackerHitsTIBHighTof;
      edm::Handle<std::vector<PSimHit>>  trackerHitsTIBLowTof;
      edm::Handle<std::vector<PSimHit>>  trackerHitsTIDHighTof;
      edm::Handle<std::vector<PSimHit>>  trackerHitsTIDLowTof;
      edm::Handle<std::vector<PSimHit>>  trackerHitsTOBHighTof;
      edm::Handle<std::vector<PSimHit>>  trackerHitsTOBLowTof;

      //
      // --- Variables
      //
      // Event
      Int_t event = 0;
      Int_t lumiBlock = 0;
      Int_t run = 0;

      //
      // --- Output
      //
      std::string output_filename;
      TH1F *counts;
      TFile *file_out;
      TTree *tree_out;

};

// Constructor
ntuplizer::ntuplizer(const edm::ParameterSet& iConfig) {

   usesResource("TFileService");

   parameters = iConfig;

   counts = new TH1F("counts", "", 1, 0, 1);

   theSimTrackCollection = consumes<std::vector<SimTrack>> (parameters.getParameter<edm::InputTag>("SimTrackCollection"));

   theTrackerHitsPixelBarrelHighTofCollection = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsPixelBarrelHighTofCollection"));
   theTrackerHitsPixelBarrelLowTofCollection = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsPixelBarrelLowTofCollection"));
   theTrackerHitsPixelEndcapHighTofCollection = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsPixelEndcapHighTofCollection"));
   theTrackerHitsPixelEndcapLowTofCollection = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsPixelEndcapLowTofCollection"));
   theTrackerHitsTECHighTofCollection        = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsTECHighTofCollection")); 
   theTrackerHitsTECLowTofCollection         = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsTECLowTofCollection")); 
   theTrackerHitsTIBHighTofCollection        = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsTIBHighTofCollection")); 
   theTrackerHitsTIBLowTofCollection         = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsTIBLowTofCollection")); 
   theTrackerHitsTIDHighTofCollection        = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsTIDHighTofCollection")); 
   theTrackerHitsTIDLowTofCollection         = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsTIDLowTofCollection")); 
   theTrackerHitsTOBHighTofCollection        = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsTOBHighTofCollection")); 
   theTrackerHitsTOBLowTofCollection         = consumes<std::vector<PSimHit>> (parameters.getParameter<edm::InputTag>("TrackerHitsTOBLowTofCollection")); 

   theGenParticleCollection = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("genParticleCollection"));

   //triggerBits_ = consumes<edm::TriggerResults> (parameters.getParameter<edm::InputTag>("bits"));
}


// Destructor
ntuplizer::~ntuplizer() {
}


// beginJob (Before first event)
void ntuplizer::beginJob() {

   std::cout << "Begin Job" << std::endl;

   // Init the file and the TTree
   output_filename = parameters.getParameter<std::string>("nameOfOutput");
   file_out = new TFile(output_filename.c_str(), "RECREATE");
   tree_out = new TTree("Events", "Events");

   // Analyzer parameters
   isData = parameters.getParameter<bool>("isData");

   // Load HLT paths

   // TTree branches
   tree_out->Branch("event", &event, "event/I");
   tree_out->Branch("lumiBlock", &lumiBlock, "lumiBlock/I");
   tree_out->Branch("run", &run, "run/I");

}

// endJob (After event loop has finished)
void ntuplizer::endJob()
{

    std::cout << "End Job" << std::endl;
    file_out->cd();
    tree_out->Write();
    counts->Write();
    file_out->Close();

}


// fillDescriptions
void ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

// Analyze (per event)
void ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   //iEvent.getByToken(muToken, muons);
   iEvent.getByToken(theTrackerHitsPixelBarrelHighTofCollection, trackerHitsPixelBarrelHighTof);
   iEvent.getByToken(theTrackerHitsPixelBarrelLowTofCollection, trackerHitsPixelBarrelLowTof);
   iEvent.getByToken(theTrackerHitsPixelEndcapHighTofCollection, trackerHitsPixelEndcapHighTof);
   iEvent.getByToken(theTrackerHitsPixelEndcapLowTofCollection, trackerHitsPixelEndcapLowTof);
   iEvent.getByToken(theTrackerHitsTECHighTofCollection,trackerHitsTECHighTof);
   iEvent.getByToken(theTrackerHitsTECLowTofCollection, trackerHitsTECLowTof);
   iEvent.getByToken(theTrackerHitsTIBHighTofCollection,trackerHitsTIBHighTof);
   iEvent.getByToken(theTrackerHitsTIBLowTofCollection, trackerHitsTIBLowTof);
   iEvent.getByToken(theTrackerHitsTIDHighTofCollection,trackerHitsTIDHighTof);
   iEvent.getByToken(theTrackerHitsTIDLowTofCollection, trackerHitsTIDLowTof);
   iEvent.getByToken(theTrackerHitsTOBHighTofCollection,trackerHitsTOBHighTof);
   iEvent.getByToken(theTrackerHitsTOBLowTofCollection, trackerHitsTOBLowTof);

   //iEvent.getByToken(triggerBits_, triggerBits);
   iEvent.getByToken(theGenParticleCollection, genParticles);
   iEvent.getByToken(theSimTrackCollection, simTracks);

   // Count number of events read
   counts->Fill(0);

   // -> Event info
   event = iEvent.id().event();
   lumiBlock = iEvent.id().luminosityBlock();
   run = iEvent.id().run();

   std::cout << "------------------------------------------ DEBUG (Event " << iEvent.id().event() << ") -------------------------------------------" << std::endl;

   // Loop over GenParticles
   std::cout << "------ GenParticles" << std::endl;
   for (size_t k = 0; k < genParticles->size(); k++) {
      int index = static_cast<int>(k);
      const reco::GenParticle &gen = (*genParticles)[k];
      if ((abs(gen.pdgId()) == 13 and gen.status() == 1) or (abs(gen.pdgId()) == 1000006 and gen.status() == 106)) {
        float Lxy = sqrt(pow(gen.vx(),2) + pow(gen.vy(),2));
        std::cout << "*** GenParticle " << gen.pdgId() << " [ID:" << index << "] || (pt,eta,phi) = (" << gen.pt() << "," << gen.eta() << "," << gen.phi() 
                                                                             << "), (vx,vy,vz) = (" << gen.vx() << "," << gen.vy() << "," << gen.vz() 
                                                                             << "), (px,py,pz) = (" << gen.px() << "," << gen.py() << "," << gen.pz() 
                                                                             << "), Lxy = " << Lxy 
                                                                             << ", longLived = " << gen.longLived() << std::endl;
      }
    }

    // Vector of PSimHits collection
    //const std::vector<std::vector<PSimHit> > trackerHitCollections = {trackerHitsPixelBarrelHighTof,trackerHitsPixelBarrelLowTof,trackerHitsPixelEndcapHighTof,trackerHitsPixelEndcapLowTof};
    

    // Loop over SimTracks
    std::cout << "------ SimTracks" << std::endl;
    for (size_t j = 0; j < simTracks->size(); j++) {
      const SimTrack &simTrack = (*simTracks)[j];
      if (simTrack.momentum().Pt() < 5) {continue;}
      if (abs(simTrack.type()) == 13) {
        std::cout << "-> SimTrack " << simTrack.type() << " || trackId = " << simTrack.trackId() 
                                                       << ", genpartIndex = " << simTrack.genpartIndex() 
                                                       << ", (pt,eta,phi) = (" << simTrack.momentum().Pt() 
                                                                        << "," << simTrack.momentum().eta() 
                                                                        << "," << simTrack.momentum().phi() << ")" << std::endl; 
    
        // Loop over SimHits
        std::cout << "    ------ PSimHits" << std::endl;
        print_trackerHits(trackerHitsPixelBarrelHighTof, simTrack);
        print_trackerHits(trackerHitsPixelBarrelLowTof, simTrack);
        print_trackerHits(trackerHitsPixelEndcapHighTof, simTrack);
        print_trackerHits(trackerHitsPixelEndcapLowTof, simTrack);
        print_trackerHits(trackerHitsTECHighTof, simTrack);
        print_trackerHits(trackerHitsTECLowTof, simTrack);
        print_trackerHits(trackerHitsTIBHighTof, simTrack);
        print_trackerHits(trackerHitsTIBLowTof, simTrack);
        print_trackerHits(trackerHitsTIDHighTof, simTrack);
        print_trackerHits(trackerHitsTIDLowTof, simTrack);
        print_trackerHits(trackerHitsTOBHighTof, simTrack);
        print_trackerHits(trackerHitsTOBLowTof, simTrack);
        /**for (size_t i = 0; i < trackerHitsPixelBarrelHighTof->size(); i++) {
          const PSimHit &hit = (*trackerHitsPixelBarrelHighTof)[i];
          if (hit.trackId() == simTrack.trackId()) {
            std::cout << "        -> Hit from particle type: " << hit.particleType() << " || trackId = " << hit.trackId() 
                                                                                     << ", localPosition = " << hit.localPosition() 
                                                                                     << ", localDirection = " << hit.localDirection() 
                                                                                     << ", momentumAtEntry = " << hit.momentumAtEntry() << std::endl;
          }
        }
        for (size_t i = 0; i < trackerHitsPixelBarrelLowTof->size(); i++) {
          const PSimHit &hit = (*trackerHitsPixelBarrelLowTof)[i];
          if (hit.trackId() == simTrack.trackId()) {
            std::cout << "        -> Hit from particle type: " << hit.particleType() << " || trackId = " << hit.trackId() 
                                                                                     << ", localPosition = " << hit.localPosition() 
                                                                                     << ", localDirection = " << hit.localDirection() 
                                                                                     << ", momentumAtEntry = " << hit.momentumAtEntry() << std::endl;
          }
        }
        for (size_t i = 0; i < trackerHitsPixelEndcapHighTof->size(); i++) {
          const PSimHit &hit = (*trackerHitsPixelEndcapHighTof)[i];
          if (hit.trackId() == simTrack.trackId()) {
            std::cout << "        -> Hit from particle type: " << hit.particleType() << " || trackId = " << hit.trackId() 
                                                                                     << ", localPosition = " << hit.localPosition() 
                                                                                     << ", localDirection = " << hit.localDirection() 
                                                                                     << ", momentumAtEntry = " << hit.momentumAtEntry() << std::endl;
          }
        }
        for (size_t i = 0; i < trackerHitsPixelEndcapLowTof->size(); i++) {
          const PSimHit &hit = (*trackerHitsPixelEndcapLowTof)[i];
          if (hit.trackId() == simTrack.trackId()) {
            std::cout << "        -> Hit from particle type: " << hit.particleType() << " || trackId = " << hit.trackId() 
                                                                                     << ", localPosition = " << hit.localPosition() 
                                                                                     << ", localDirection = " << hit.localDirection() 
                                                                                     << ", momentumAtEntry = " << hit.momentumAtEntry() << std::endl;
          }
        }**/
      }
    }

   std::cout << "----------------------------------------      END DEBUG       -----------------------------------------" << std::endl;

   //-> Fill tree
   tree_out->Fill();

}

DEFINE_FWK_MODULE(ntuplizer);
