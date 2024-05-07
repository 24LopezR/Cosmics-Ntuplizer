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
#include "DataFormats/DetId/interface/DetId.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/TrackerGeomDet.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"

float deltaR(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1 - eta2;
  float dphi = abs(phi1 - phi2);
  if (dphi > float(M_PI)) {
    dphi -= float(2*M_PI);
  }
  return sqrt(deta*deta + dphi*dphi);
}

float dxy_value(const reco::GenParticle &p, const reco::Vertex &pv){
    float vx = p.vx();
    float vy = p.vy();
    float phi = p.phi();
    float pv_x = pv.x();
    float pv_y = pv.y();
  
    float dxy = -(vx-pv_x)*sin(phi) + (vy-pv_y)*cos(phi);
    return dxy;
}

void print_trackerHit(const PSimHit *hit, const SimTrack &simTrack, const TrackerGeometry *trackerGeom) {
      std::cout << "        -> Hit from particle type: " << hit->particleType() << " || trackId = " << hit->trackId() 
                << "\n                                               localPosition = " << hit->localPosition() 
                << "\n                                               localDirection = " << hit->localDirection() 
                << "\n                                               momentumAtEntry = " << hit->momentumAtEntry()
                << "\n                                               detId = " << hit->detUnitId()
                << "\n                                               tof = " << hit->timeOfFlight() << std::endl;
    
      const DetId id(hit->detUnitId());
      const TrackerGeomDet* trkGeomDet = trackerGeom->idToDet(id);
      const auto& localPos = hit->localPosition();
      const auto& globalPos = trkGeomDet->toGlobal(localPos);
      std::cout << "                                               globalPos = " << globalPos << std::endl;
}

void print_trackerHits(const edm::Handle<std::vector<PSimHit>> trackerHits, const SimTrack &simTrack, const TrackerGeometry *trackerGeom) {
  for (size_t i = 0; i < trackerHits->size(); i++) {
    const PSimHit &hit = (*trackerHits)[i];
    if (hit.trackId() == simTrack.trackId()) {
      std::cout << "        -> Hit from particle type: " << hit.particleType() << " || trackId = " << hit.trackId() 
                << "\n                                               localPosition = " << hit.localPosition() 
                << "\n                                               localDirection = " << hit.localDirection() 
                << "\n                                               momentumAtEntry = " << hit.momentumAtEntry()
                << "\n                                               detId = " << hit.detUnitId()
                << "\n                                               tof = " << hit.timeOfFlight() << std::endl;
    
      const DetId id(hit.detUnitId());
      const TrackerGeomDet* trkGeomDet = trackerGeom->idToDet(id);
      const auto& localPos = hit.localPosition();
      const auto& globalPos = trkGeomDet->toGlobal(localPos);
      std::cout << "                                               globalPos = " << globalPos << std::endl;
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
      virtual void reset();

      edm::ParameterSet parameters;

      bool isData = true;
      bool verbose = true;
      double dRGenSim = 1.0;

      //
      // --- Tokens and Handles
      //

      // trigger bits
      //edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      //edm::Handle<edm::TriggerResults> triggerBits;
    
      const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken;
      const TrackerGeometry *trackerGeom;
      
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


      // GenParticles
      std::vector<int  > gen_pdgId  ; 
      std::vector<int  > gen_status ;
      std::vector<float> gen_pt     ;
      std::vector<float> gen_eta    ;
      std::vector<float> gen_phi    ;
      std::vector<float> gen_vx     ; 
      std::vector<float> gen_vy     ; 
      std::vector<float> gen_vz     ; 
      std::vector<float> gen_px     ; 
      std::vector<float> gen_py     ; 
      std::vector<float> gen_pz     ; 
      std::vector<float> gen_Lxy    ;

      // SimTracks
      std::vector<int  > simTrack_trackId          ;
      std::vector<int  > simTrack_genPartIndex     ;
      std::vector<int  > simTrack_genPartType      ;
      std::vector<float> simTrack_pt               ;
      std::vector<float> simTrack_eta              ;
      std::vector<float> simTrack_phi              ;
 
      // SimHits
      std::vector<int  > simHit_trackId           ;
      std::vector<int  > simHit_detUnitId         ;
      std::vector<int  > simHit_particleType      ;
      std::vector<float> simHit_tof               ;
      std::vector<float> simHit_pabs              ;
      std::vector<float> simHit_momentumAtEntry_x ;
      std::vector<float> simHit_momentumAtEntry_y ;
      std::vector<float> simHit_momentumAtEntry_z ;
      std::vector<float> simHit_globalPosition_x  ;
      std::vector<float> simHit_globalPosition_y  ;
      std::vector<float> simHit_globalPosition_z  ;
      std::vector<float> simHit_genLxy            ;

      //
      // --- Output
      //
      std::string output_filename;
      TH1F *counts;
      TFile *file_out;
      TTree *tree_out;

};

// Constructor
ntuplizer::ntuplizer(const edm::ParameterSet& iConfig) : trackerGeomToken(esConsumes<edm::Transition::Event>()) {

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
   verbose = parameters.getParameter<bool>("verbose");
   dRGenSim = parameters.getParameter<double>("dRGenSim");

   // Load HLT paths

   // TTree branches
   tree_out->Branch("event", &event, "event/I");
   tree_out->Branch("lumiBlock", &lumiBlock, "lumiBlock/I");
   tree_out->Branch("run", &run, "run/I");
      
      // GenParticles
   tree_out->Branch("gen_pdgId",  &gen_pdgId ); 
   tree_out->Branch("gen_status", &gen_status);
   tree_out->Branch("gen_pt",     &gen_pt    );    
   tree_out->Branch("gen_eta",    &gen_eta   );   
   tree_out->Branch("gen_phi",    &gen_phi   );   
   tree_out->Branch("gen_vx",     &gen_vx    );    
   tree_out->Branch("gen_vy",     &gen_vy    );    
   tree_out->Branch("gen_vz",     &gen_vz    );    
   tree_out->Branch("gen_px",     &gen_px    );    
   tree_out->Branch("gen_py",     &gen_py    );    
   tree_out->Branch("gen_pz",     &gen_pz    );    
   tree_out->Branch("gen_Lxy",    &gen_Lxy   );   

   // SimTracks
   tree_out->Branch("simTrack_trackId"      , &simTrack_trackId      ) ;
   tree_out->Branch("simTrack_genPartIndex" , &simTrack_genPartIndex ) ;
   tree_out->Branch("simTrack_genPartType"  , &simTrack_genPartType  ) ;
   tree_out->Branch("simTrack_pt"           , &simTrack_pt           ) ;
   tree_out->Branch("simTrack_eta"          , &simTrack_eta          ) ;
   tree_out->Branch("simTrack_phi"          , &simTrack_phi          ) ;

   // SimHits
   tree_out->Branch("simHit_trackId"          , &simHit_trackId          );
   tree_out->Branch("simHit_detUnitId"        , &simHit_detUnitId        );
   tree_out->Branch("simHit_particleType"     , &simHit_particleType     );
   tree_out->Branch("simHit_tof"              , &simHit_tof              );
   tree_out->Branch("simHit_pabs"             , &simHit_pabs             );
   tree_out->Branch("simHit_momentumAtEntry_x", &simHit_momentumAtEntry_x);
   tree_out->Branch("simHit_momentumAtEntry_y", &simHit_momentumAtEntry_y);
   tree_out->Branch("simHit_momentumAtEntry_z", &simHit_momentumAtEntry_z);
   tree_out->Branch("simHit_globalPosition_x" , &simHit_globalPosition_x );
   tree_out->Branch("simHit_globalPosition_y" , &simHit_globalPosition_y );
   tree_out->Branch("simHit_globalPosition_z" , &simHit_globalPosition_z );
   tree_out->Branch("simHit_genLxy"           , &simHit_genLxy           );

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

   reset();
   trackerGeom = &iSetup.getData(trackerGeomToken);

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

   iEvent.getByToken(theGenParticleCollection, genParticles);
   iEvent.getByToken(theSimTrackCollection, simTracks);

   // Count number of events read
   counts->Fill(0);

   // -> Event info
   event = iEvent.id().event();
   lumiBlock = iEvent.id().luminosityBlock();
   run = iEvent.id().run();

   if (verbose) {std::cout << "------------------------------------------ DEBUG (Event " << iEvent.id().event() << ") -------------------------------------------" << std::endl;}

   // Loop over GenParticles
   if (verbose) {std::cout << "------ GenParticles" << std::endl;}
   for (size_t k = 0; k < genParticles->size(); k++) {
      const reco::GenParticle &gen = (*genParticles)[k];
      float Lxy = sqrt(pow(gen.vx(),2) + pow(gen.vy(),2));
      gen_pdgId.push_back(gen.pdgId());
      gen_status.push_back(gen.status());
      gen_pt.push_back(gen.pt());
      gen_eta.push_back(gen.eta());
      gen_phi.push_back(gen.phi());
      gen_vx.push_back(gen.vx());
      gen_vy.push_back(gen.vy());
      gen_vz.push_back(gen.vz());
      gen_px.push_back(gen.px());
      gen_py.push_back(gen.py());
      gen_pz.push_back(gen.pz());
      gen_Lxy.push_back(Lxy);
    }

    // Loop over SimTracks
    if (verbose) {std::cout << "------ SimTracks" << std::endl;}
    for (size_t j = 0; j < simTracks->size(); j++) {
      const SimTrack &simTrack = (*simTracks)[j];
      if (abs(simTrack.type()) != 13) {continue;}
      simTrack_trackId.push_back(simTrack.trackId());
      simTrack_genPartIndex.push_back(simTrack.genpartIndex());
      simTrack_genPartType.push_back(simTrack.type());
      simTrack_pt.push_back(simTrack.momentum().Pt());
      simTrack_eta.push_back(simTrack.momentum().eta());
      simTrack_phi.push_back(simTrack.momentum().phi());
      if (simTrack.momentum().Pt() < 5) {continue;}
      // Loop over SimHits
      /**if (verbose) {
        std::cout << "    ------ PSimHits" << std::endl;
        print_trackerHits(trackerHitsPixelBarrelHighTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsPixelBarrelLowTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsPixelEndcapHighTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsPixelEndcapLowTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsTECHighTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsTECLowTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsTIBHighTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsTIBLowTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsTIDHighTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsTIDLowTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsTOBHighTof, simTrack, trackerGeom);
        print_trackerHits(trackerHitsTOBLowTof, simTrack, trackerGeom);
      }**/
    }

    // Loop over tracker hit collections
    std::vector<edm::Handle<std::vector<PSimHit>> > trackerHitCollections;
    trackerHitCollections.push_back( trackerHitsPixelBarrelHighTof );
    trackerHitCollections.push_back( trackerHitsPixelBarrelLowTof );
    trackerHitCollections.push_back( trackerHitsPixelEndcapHighTof );
    trackerHitCollections.push_back( trackerHitsPixelEndcapLowTof );
    trackerHitCollections.push_back( trackerHitsTECHighTof);
    trackerHitCollections.push_back( trackerHitsTECLowTof);
    trackerHitCollections.push_back( trackerHitsTIBHighTof);
    trackerHitCollections.push_back( trackerHitsTIBLowTof);
    trackerHitCollections.push_back( trackerHitsTIDHighTof);
    trackerHitCollections.push_back( trackerHitsTIDLowTof);
    trackerHitCollections.push_back( trackerHitsTOBHighTof);
    trackerHitCollections.push_back( trackerHitsTOBLowTof);
    for (auto & trackerHits : trackerHitCollections) {
      for (size_t i = 0; i < trackerHits->size(); i++) {
        const PSimHit *hit = &(*trackerHits)[i];
        if (abs(hit->particleType()) != 13) {continue;}
        // Loop over SimTracks to determine matches
        for (size_t j = 0; j < simTracks->size(); j++) {
          const SimTrack &simTrack = (*simTracks)[j];
          if (abs(simTrack.type()) != 13) {continue;}
          if (simTrack.trackId() != hit->trackId()) {continue;}
          //std::cout << "Hit " << hit->trackId() << " belongs to simTrack " << simTrack.trackId() << std::endl;
          for (size_t k = 0; k < genParticles->size(); k++) {
             const reco::GenParticle &gen = (*genParticles)[k];
             if (abs(gen.pdgId()) != 13 or gen.status() != 1) {continue;}
             if (deltaR(gen.eta(),gen.phi(),simTrack.momentum().eta(),simTrack.momentum().phi()) > dRGenSim) {continue;}
             float Lxy = sqrt(pow(gen.vx(),2) + pow(gen.vy(),2));
             int index = static_cast<int>(k);
             std::cout << "*** GenParticle " << gen.pdgId() << ", status = " << gen.status() << " [ID:" << index << "] || (pt,eta,phi) = (" << gen.pt() << "," << gen.eta() << "," << gen.phi() 
                                                                                  << "), (vx,vy,vz) = (" << gen.vx() << "," << gen.vy() << "," << gen.vz() 
                                                                                  << "), (px,py,pz) = (" << gen.px() << "," << gen.py() << "," << gen.pz() 
                                                                                  << "), Lxy = " << Lxy 
                                                                                  << ", longLived = " << gen.longLived() << std::endl;
             std::cout << "-> SimTrack " << simTrack.type() << " || trackId = " << simTrack.trackId() 
                                                            << ", genpartIndex = " << simTrack.genpartIndex() 
                                                            << ", (pt,eta,phi) = (" << simTrack.momentum().Pt() 
                                                                             << "," << simTrack.momentum().eta() 
                                                                             << "," << simTrack.momentum().phi() << ")" << std::endl; 
             print_trackerHit(hit, simTrack, trackerGeom);
             const DetId id(hit->detUnitId());
             const TrackerGeomDet* trkGeomDet = trackerGeom->idToDet(id);
             const auto& localPos = hit->localPosition();
             const auto& globalPos = trkGeomDet->toGlobal(localPos);
             simHit_trackId.push_back(hit->trackId());
             simHit_detUnitId.push_back(id);
             simHit_particleType.push_back(hit->particleType());
             simHit_tof.push_back(hit->timeOfFlight());
             simHit_pabs.push_back(hit->pabs());
             simHit_momentumAtEntry_x.push_back(hit->momentumAtEntry().x());
             simHit_momentumAtEntry_y.push_back(hit->momentumAtEntry().y());
             simHit_momentumAtEntry_z.push_back(hit->momentumAtEntry().z());
             simHit_globalPosition_x.push_back(globalPos.x());
             simHit_globalPosition_y.push_back(globalPos.y());
             simHit_globalPosition_z.push_back(globalPos.z());

             float genLxy = sqrt(pow(gen.vx(),2) + pow(gen.vy(),2));
             simHit_genLxy.push_back(genLxy);
          }
        }
      }
    }

   if (verbose) {std::cout << "----------------------------------------      END DEBUG       -----------------------------------------" << std::endl;}

   //-> Fill tree
   tree_out->Fill();

}

void ntuplizer::reset() {

      // GenParticles
      gen_pdgId .clear() ; 
      gen_status.clear() ;
      gen_pt    .clear() ;
      gen_eta   .clear() ;
      gen_phi   .clear() ;
      gen_vx    .clear() ; 
      gen_vy    .clear() ; 
      gen_vz    .clear() ; 
      gen_px    .clear() ; 
      gen_py    .clear() ; 
      gen_pz    .clear() ; 
      gen_Lxy   .clear() ;

      // SimTracks
      simTrack_trackId     .clear()     ;
      simTrack_genPartIndex.clear()     ;
      simTrack_genPartType .clear()     ;
      simTrack_pt          .clear()     ;
      simTrack_eta         .clear()     ;
      simTrack_phi         .clear()     ;
 
      // SimHits
      simHit_trackId          .clear() ;
      simHit_detUnitId        .clear() ;
      simHit_particleType     .clear() ;
      simHit_tof              .clear() ;
      simHit_pabs             .clear() ;
      simHit_momentumAtEntry_x.clear() ;
      simHit_momentumAtEntry_y.clear() ;
      simHit_momentumAtEntry_z.clear() ;
      simHit_globalPosition_x .clear() ;
      simHit_globalPosition_y .clear() ;
      simHit_globalPosition_z .clear() ;
      simHit_genLxy           .clear() ;
}

DEFINE_FWK_MODULE(ntuplizer);
