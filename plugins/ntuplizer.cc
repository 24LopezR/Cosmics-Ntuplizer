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
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::Handle<edm::TriggerResults> triggerBits;

      // muons (reco::Muon)
      edm::EDGetTokenT<edm::View<reco::Muon> > muToken;
      edm::Handle<edm::View<reco::Muon> > muons;
    
      // PrimaryVertices
      edm::EDGetTokenT<edm::View<reco::Vertex> > thePrimaryVertexCollection;
      edm::Handle<edm::View<reco::Vertex> > primaryvertices;
      
      // GenParticles
      edm::EDGetTokenT<edm::View<reco::GenParticle> >  theGenParticleCollection;
      edm::Handle<edm::View<reco::GenParticle> > genParticles;

      //
      // --- Variables
      //


      // Trigger tags
      std::vector<std::string> HLTPaths_;
      bool triggerPass[200] = {false};

      // Event
      Int_t event = 0;
      Int_t lumiBlock = 0;
      Int_t run = 0;

      // ----------------------------------
      // muon variables
      // ----------------------------------
      Int_t nmu = 0;
      Int_t mu_isSTA[200] = {0};
      Int_t mu_isGLB[200] = {0};
      Int_t mu_isTRK[200] = {0};
      Int_t mu_isMatchesValid[200] = {0};
      Int_t mu_numberOfMatches[200] = {0};
      Int_t mu_numberOfChambers[200] = {0};
      Int_t mu_numberOfChambersCSCorDT[200] = {0};
      Int_t mu_numberOfMatchedStations[200] = {0};
      Int_t mu_numberOfMatchedRPCLayers[200] = {0};
      Float_t mu_pt[200] = {0.};
      Float_t mu_eta[200] = {0.};
      Float_t mu_phi[200] = {0.};
      Float_t mu_ptError[200] = {0.};
      Float_t mu_dxy[200] = {0.};
      Float_t mu_dz[200] = {0.};
      Float_t mu_normalizedChi2[200] = {0.};
      Float_t mu_charge[200] = {0.};
      Int_t mu_nMuonHits[200] = {0};
      Int_t mu_nValidMuonHits[200] = {0};
      Int_t mu_nValidMuonDTHits[200] = {0};
      Int_t mu_nValidMuonCSCHits[200] = {0};
      Int_t mu_nValidMuonRPCHits[200] = {0};
      Int_t mu_nValidStripHits[200] = {0};
      Int_t mu_nhits[200] = {0};
      Int_t mu_dtStationsWithValidHits[200] = {0};
      Int_t mu_cscStationsWithValidHits[200] = {0};
      Int_t mu_nsegments[200] = {0};

      Float_t mu_iso03_sumPt[200]     = {0.};
      Float_t mu_iso03_emEt[200]      = {0.};
      Float_t mu_iso03_hadEt[200]     = {0.};
      Float_t mu_pfIso03_charged[200] = {0.}; 
      Float_t mu_pfIso03_neutral[200] = {0.};
      Float_t mu_pfIso03_photon [200] = {0.};
      Float_t mu_pfIso03_sumPU  [200] = {0.};
      Float_t mu_pfIso04_charged[200] = {0.};
      Float_t mu_pfIso04_neutral[200] = {0.};
      Float_t mu_pfIso04_photon [200] = {0.};
      Float_t mu_pfIso04_sumPU  [200] = {0.};

      // ----------------------------------
      // PrimaryVertices
      // ----------------------------------
      Int_t nPV;
      Int_t nTruePV;
      Int_t PV_passAcceptance;
      Float_t PV_vx;
      Float_t PV_vy;
      Float_t PV_vz;

      // ----------------------------------
      // GenParticles
      // ----------------------------------
      Int_t nGenMuon;
      Int_t nGenMuon_PFS;
      Int_t nGenMuon_HPFS;
      Int_t nGenMuon_PTDP;
      Int_t nGenMuon_HDP;
      Float_t GenLeptonSel_pt[30];
      Float_t GenLeptonSel_E[30];
      Float_t GenLeptonSel_et[30];
      Float_t GenLeptonSel_eta[30];
      Float_t GenLeptonSel_phi[30];
      Int_t GenLeptonSel_pdgId[30];
      Float_t GenLeptonSel_dxy[30];
      Float_t GenLeptonSel_vx[30];
      Float_t GenLeptonSel_vy[30];
      Float_t GenLeptonSel_vz[30];
      Int_t GenLeptonSel_motherPdgId[30];
      Int_t GenLeptonSel_fromHardProcessFinalState[30];
      Int_t GenLeptonSel_isPromptFinalState[30];
      Int_t GenLeptonSel_isDirectPromptTauDecayProductFinalState[30];
      Int_t GenLeptonSel_isDirectHadronDecayProduct[30];
      
      Int_t nHardProcessParticle;
      Float_t HardProcessParticle_pt[30];
      Float_t HardProcessParticle_E[30];
      Float_t HardProcessParticle_eta[30];
      Float_t HardProcessParticle_phi[30];
      Float_t HardProcessParticle_vx[30];
      Float_t HardProcessParticle_vy[30];
      Float_t HardProcessParticle_vz[30];
      Int_t HardProcessParticle_pdgId[30];

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

   muToken = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("muonCollection"));
   if (!isData) {
     theGenParticleCollection = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("genParticleCollection"));
     thePrimaryVertexCollection = consumes<edm::View<reco::Vertex> >  (parameters.getParameter<edm::InputTag>("PrimaryVertexCollection"));
   }

   triggerBits_ = consumes<edm::TriggerResults> (parameters.getParameter<edm::InputTag>("bits"));
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

   // ----------------------------------
   // Muons
   // ----------------------------------
   tree_out->Branch("nmu", &nmu, "nmu/I");
   tree_out->Branch("mu_isSTA", mu_isSTA, "mu_isSTA[nmu]/I");
   tree_out->Branch("mu_isGLB", mu_isGLB, "mu_isGLB[nmu]/I");
   tree_out->Branch("mu_isTRK", mu_isTRK, "mu_isTRK[nmu]/I");
   tree_out->Branch("mu_isMatchesValid", mu_isMatchesValid, "mu_isMatchesValid[nmu]/I");
   tree_out->Branch("mu_numberOfMatches", mu_numberOfMatches, "mu_numberOfMatches[nmu]/I");
   tree_out->Branch("mu_numberOfChambers", mu_numberOfChambers, "mu_numberOfChambers[nmu]/I");
   tree_out->Branch("mu_numberOfChambersCSCorDT", mu_numberOfChambersCSCorDT, "mu_numberOfChambersCSCorDT[nmu]/I");
   tree_out->Branch("mu_numberOfMatchedStations", mu_numberOfMatchedStations, "mu_numberOfMatchedStations[nmu]/I");
   tree_out->Branch("mu_numberOfMatchedRPCLayers", mu_numberOfMatchedRPCLayers, "mu_numberOfMatchedRPCLayers[nmu]/I");
   tree_out->Branch("mu_pt", mu_pt, "mu_pt[nmu]/F");
   tree_out->Branch("mu_eta", mu_eta, "mu_eta[nmu]/F");
   tree_out->Branch("mu_phi", mu_phi, "mu_phi[nmu]/F");
   tree_out->Branch("mu_ptError", mu_ptError, "mu_ptError[nmu]/F");
   tree_out->Branch("mu_dxy", mu_dxy, "mu_dxy[nmu]/F");
   tree_out->Branch("mu_dz", mu_dz, "mu_dz[nmu]/F");
   tree_out->Branch("mu_normalizedChi2", mu_normalizedChi2, "mu_normalizedChi2[nmu]/F");
   tree_out->Branch("mu_charge", mu_charge, "mu_charge[nmu]/F");
   tree_out->Branch("mu_nMuonHits", mu_nMuonHits, "mu_nMuonHits[nmu]/I");
   tree_out->Branch("mu_nValidMuonHits", mu_nValidMuonHits, "mu_nValidMuonHits[nmu]/I");
   tree_out->Branch("mu_nValidMuonDTHits", mu_nValidMuonDTHits, "mu_nValidMuonDTHits[nmu]/I");
   tree_out->Branch("mu_nValidMuonCSCHits", mu_nValidMuonCSCHits, "mu_nValidMuonCSCHits[nmu]/I");
   tree_out->Branch("mu_nValidMuonRPCHits", mu_nValidMuonRPCHits, "mu_nValidMuonRPCHits[nmu]/I");
   tree_out->Branch("mu_nValidStripHits", mu_nValidStripHits, "mu_nValidStripHits[nmu]/I");
   tree_out->Branch("mu_nhits", mu_nhits, "mu_nhits[nmu]/I");
   tree_out->Branch("mu_dtStationsWithValidHits", mu_dtStationsWithValidHits, "mu_dtStationsWithValidHits[nmu]/I");
   tree_out->Branch("mu_cscStationsWithValidHits", mu_cscStationsWithValidHits, "mu_cscStationsWithValidHits[nmu]/I");
   tree_out->Branch("mu_nsegments", mu_nsegments, "mu_nsegments[nmu]/I");

   tree_out->Branch("mu_iso03_sumPt",     mu_iso03_sumPt,     "mu_iso03_sumPt[nmu]/F"    ); 
   tree_out->Branch("mu_iso03_emEt",      mu_iso03_emEt,      "mu_iso03_emEt[nmu]/F"     ); 
   tree_out->Branch("mu_iso03_hadEt",     mu_iso03_hadEt,     "mu_iso03_hadEt[nmu]/F"    ); 
   tree_out->Branch("mu_pfIso03_charged", mu_pfIso03_charged, "mu_pfIso03_charged[nmu]/F"); 
   tree_out->Branch("mu_pfIso03_neutral", mu_pfIso03_neutral, "mu_pfIso03_neutral[nmu]/F"); 
   tree_out->Branch("mu_pfIso03_photon",  mu_pfIso03_photon,  "mu_pfIso03_photon[nmu]/F" ); 
   tree_out->Branch("mu_pfIso03_sumPU",   mu_pfIso03_sumPU,   "mu_pfIso03_sumPU[nmu]/F"  ); 
   tree_out->Branch("mu_pfIso04_charged", mu_pfIso04_charged, "mu_pfIso04_charged[nmu]/F"); 
   tree_out->Branch("mu_pfIso04_neutral", mu_pfIso04_neutral, "mu_pfIso04_neutral[nmu]/F"); 
   tree_out->Branch("mu_pfIso04_photon",  mu_pfIso04_photon,  "mu_pfIso04_photon[nmu]/F" ); 
   tree_out->Branch("mu_pfIso04_sumPU",   mu_pfIso04_sumPU,   "mu_pfIso04_sumPU[nmu]/F"  ); 

   if (!isData) {
     // ----------------------------------
     // PrimaryVertices
     // ----------------------------------
     tree_out->Branch("nPV", &nPV, "nPV/I");
     tree_out->Branch("nTruePV", &nTruePV, "nTruePV/I");
     tree_out->Branch("PV_passAcceptance", &PV_passAcceptance, "PV_passAcceptance/I");
     tree_out->Branch("PV_vx", &PV_vx, "PV_vx/F");
     tree_out->Branch("PV_vy", &PV_vy, "PV_vy/F");
     tree_out->Branch("PV_vz", &PV_vz, "PV_vz/F");   

     // ----------------------------------
     // GenParticles
     // ----------------------------------
     tree_out->Branch("nGenMuon", &nGenMuon, "nGenMuon/I");
     tree_out->Branch("GenLeptonSel_pt", GenLeptonSel_pt, "GenLeptonSel_pt[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_E", GenLeptonSel_E, "GenLeptonSel_E[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_et", GenLeptonSel_et, "GenLeptonSel_et[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_eta", GenLeptonSel_eta, "GenLeptonSel_eta[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_phi", GenLeptonSel_phi, "GenLeptonSel_phi[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_dxy", GenLeptonSel_dxy, "GenLeptonSel_dxy[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_vx", GenLeptonSel_vx, "GenLeptonSel_vx[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_vy", GenLeptonSel_vy, "GenLeptonSel_vy[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_vz", GenLeptonSel_vz, "GenLeptonSel_vz[nGenMuon]/F");
     tree_out->Branch("GenLeptonSel_pdgId", GenLeptonSel_pdgId, "GenLeptonSel_pdgId[nGenMuon]/I");
     tree_out->Branch("GenLeptonSel_motherPdgId", GenLeptonSel_motherPdgId, "GenLeptonSel_motherPdgId[nGenMuon]/I");
     tree_out->Branch("GenLeptonSel_isPromptFinalState", GenLeptonSel_isPromptFinalState, "GenLeptonSel_isPromptFinalState[nGenMuon]/I");
     tree_out->Branch("GenLeptonSel_fromHardProcessFinalState", GenLeptonSel_fromHardProcessFinalState, "GenLeptonSel_fromHardProcessFinalState[nGenMuon]/I");
     tree_out->Branch("GenLeptonSel_isDirectPromptTauDecayProductFinalState", GenLeptonSel_isDirectPromptTauDecayProductFinalState, "GenLeptonSel_isDirectPromptTauDecayProductFinalState[nGenMuon]/I");
     tree_out->Branch("GenLeptonSel_isDirectHadronDecayProduct", GenLeptonSel_isDirectHadronDecayProduct, "GenLeptonSel_isDirectHadronDecayProduct[nGenMuon]/I");
 
     tree_out->Branch("nHardProcessParticle", &nHardProcessParticle, "nHardProcessParticle/I");
     tree_out->Branch("HardProcessParticle_E", HardProcessParticle_E, "HardProcessParticle_E[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_pt", HardProcessParticle_pt, "HardProcessParticle_pt[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_eta", HardProcessParticle_eta, "HardProcessParticle_eta[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_phi", HardProcessParticle_phi, "HardProcessParticle_phi[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_vx", HardProcessParticle_vx, "HardProcessParticle_vx[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_vy", HardProcessParticle_vy, "HardProcessParticle_vy[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_vz", HardProcessParticle_vz, "HardProcessParticle_vz[nHardProcessParticle]/F");
     tree_out->Branch("HardProcessParticle_pdgId", HardProcessParticle_pdgId, "HardProcessParticle_pdgId[nHardProcessParticle]/I");
   }

   // Trigger branches
   for (unsigned int ihlt = 0; ihlt < HLTPaths_.size(); ihlt++) {
     tree_out->Branch(TString(HLTPaths_[ihlt]), &triggerPass[ihlt]);
   }

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

   iEvent.getByToken(muToken, muons);
   iEvent.getByToken(triggerBits_, triggerBits);
   if (!isData){
      iEvent.getByToken(thePrimaryVertexCollection, primaryvertices);
      iEvent.getByToken(theGenParticleCollection, genParticles);
   }

   // Count number of events read
   counts->Fill(0);


   // -> Event info
   event = iEvent.id().event();
   lumiBlock = iEvent.id().luminosityBlock();
   run = iEvent.id().run();

   // ----------------------------------
   // displacedMuons Collection
   // ----------------------------------
   nmu = 0;;
   for (unsigned int i = 0; i < muons->size(); i++) {
     std::cout << " - - nmu: " << nmu << std::endl;
     const reco::Muon& muon(muons->at(i));
     mu_isDGL[nmu] = muon.isGlobalMuon();
     mu_isDSA[nmu] = muon.isStandAloneMuon();
     mu_isDTK[nmu] = muon.isTrackerMuon();
     mu_isMatchesValid[nmu] = muon.isMatchesValid();
     mu_numberOfMatches[nmu] = muon.numberOfMatches();
     mu_numberOfChambers[nmu] = muon.numberOfChambers();
     mu_numberOfChambersCSCorDT[nmu] = muon.numberOfChambersCSCorDT();
     mu_numberOfMatchedStations[nmu] = muon.numberOfMatchedStations();
     mu_numberOfMatchedRPCLayers[nmu] = muon.numberOfMatchedRPCLayers();
     mu_pt[nmu] = muon.pt();
     mu_eta[nmu] = muon.eta();
     mu_phi[nmu] = muon.phi();
     mu_ptError[nmu] = muon.ptError();
     mu_dxy[nmu] = muon.dxy();
     mu_dz[nmu] = muon.dz();
     mu_normalizedChi2[nmu] = muon.normalizedChi2();
     mu_charge[nmu] = muon.charge();
     mu_nMuonHits[nmu] = muon.track()->hitPattern().numberOfMuonHits();
     mu_nValidMuonHits[nmu] = muon.track()->hitPattern().numberOfValidMuonHits();
     mu_nValidMuonDTHits[nmu] = muon.track()->hitPattern().numberOfValidMuonDTHits();
     mu_nValidMuonCSCHits[nmu] = muon.track()->hitPattern().numberOfValidMuonCSCHits();
     mu_nValidMuonRPCHits[nmu] = muon.track()->hitPattern().numberOfValidMuonRPCHits();
     mu_nValidStripHits[nmu] = muon.track()->hitPattern().numberOfValidStripHits();
     mu_nhits[nmu] = muon.track()->hitPattern().numberOfValidHits();
     mu_dtStationsWithValidHits[nmu] = muon.track()->hitPattern().dtStationsWithValidHits();
     mu_cscStationsWithValidHits[nmu] = muon.track()->hitPattern().cscStationsWithValidHits();
     if ( muon.isStandAloneMuon() ) {
       const reco::Track* outerTrack = (muon.standAloneMuon()).get();
       // Number of DT+CSC segments
       unsigned int nsegments = 0;
       for (trackingRecHit_iterator hit = outerTrack->recHitsBegin(); hit != outerTrack->recHitsEnd(); ++hit) {
         if (!(*hit)->isValid()) continue;
         DetId id = (*hit)->geographicalId();
         if (id.det() != DetId::Muon) continue;
         if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
           nsegments++;
         }
       }
       mu_nsegments[nmu] = nsegments;
     } else {
       mu_nsegments[nmu] = 0;
     }
     
     mu_iso03_sumPt[nmu]     = muon.isolationR03().sumPt;
     mu_iso03_emEt[nmu]      = muon.isolationR03().emEt;
     mu_iso03_hadEt[nmu]     = muon.isolationR03().hadEt;
     mu_pfIso03_charged[nmu] = muon.pfIsolationR03().sumChargedHadronPt; 
     mu_pfIso03_neutral[nmu] = muon.pfIsolationR03().sumNeutralHadronEt;
     mu_pfIso03_photon[nmu]  = muon.pfIsolationR03().sumPhotonEt;
     mu_pfIso03_sumPU[nmu]   = muon.pfIsolationR03().sumPUPt;
     mu_pfIso04_charged[nmu] = muon.pfIsolationR04().sumChargedHadronPt; 
     mu_pfIso04_neutral[nmu] = muon.pfIsolationR04().sumNeutralHadronEt;
     mu_pfIso04_photon[nmu]  = muon.pfIsolationR04().sumPhotonEt;
     mu_pfIso04_sumPU[nmu]   = muon.pfIsolationR04().sumPUPt;

     nmu++;
     std::cout << "End muon" << std::endl;
   }

   if (!isData) {
     // ----------------------------------
     // PrimaryVertices collection
     // ----------------------------------
     nTruePV = 0;
     nPV = primaryvertices->size();

     for (size_t i = 0; i < primaryvertices->size(); i ++){
         const reco::Vertex &current_vertex = (*primaryvertices)[i];
         if(current_vertex.isValid()){ nTruePV++; }
     }

     // The PV information:
     const reco::Vertex &thePrimaryVertex = (*primaryvertices)[0];
     
     PV_vx = thePrimaryVertex.x();
     PV_vy = thePrimaryVertex.y();
     PV_vz = thePrimaryVertex.z();
     PV_passAcceptance = false;
     
     if (!thePrimaryVertex.isFake() && thePrimaryVertex.ndof() > 4 && fabs(thePrimaryVertex.z()) < 25 && thePrimaryVertex.position().rho() <= 2) {
        PV_passAcceptance = true;
     }   
     GlobalPoint _PVpoint(thePrimaryVertex.x(), thePrimaryVertex.y(), thePrimaryVertex.z());

     // ----------------------------------
     // GenParticle collection
     // ----------------------------------
     std::vector<int> iGM; // gen muons count

     reco::GenParticleRef mref;
     reco::GenParticle m;
     if (!isData) {
        // Get the muons that can be reconstructed
        for (size_t i = 0; i < genParticles->size(); i++) {
           const reco::GenParticle &genparticle = (*genParticles)[i];
             if ( abs(genparticle.pdgId()) == 13 && genparticle.status() == 1) {
                 iGM.push_back(i);
             }
         }
         nGenMuon = iGM.size();
         std::cout << "Number of gen muons = " << nGenMuon << std::endl;

         for (size_t i = 0; i < iGM.size(); i++) {
            const reco::GenParticle &genparticle = (*genParticles)[iGM.at(i)];
            
            GenLeptonSel_pt[i] = genparticle.pt();
            GenLeptonSel_E[i] = genparticle.energy();
            GenLeptonSel_et[i] = genparticle.et();
            GenLeptonSel_eta[i] = genparticle.eta();
            GenLeptonSel_phi[i] = genparticle.phi();
            GenLeptonSel_pdgId[i] = genparticle.pdgId();

            // Bottom-up to get the real decaying particle:
            if (genparticle.mother()->pdgId() == genparticle.pdgId()) {

                mref = genparticle.motherRef();
                m = *mref;
                while (m.pdgId() == m.mother()->pdgId()) {
                    mref = m.motherRef();
                    m = *mref;
                }

                GenLeptonSel_vx[i] = m.vx();
                GenLeptonSel_vy[i] = m.vy();
                GenLeptonSel_vz[i] = m.vz();
                GenLeptonSel_dxy[i] = dxy_value(m, thePrimaryVertex); // should be computed here or before?

                if(m.numberOfMothers() != 0){
                    GenLeptonSel_motherPdgId[i] = m.motherRef()->pdgId();
                } else {
                    GenLeptonSel_motherPdgId[i] = 0; 
                }
            } else {

                GenLeptonSel_vx[i] = genparticle.vx();
                GenLeptonSel_vy[i] = genparticle.vy();
                GenLeptonSel_vz[i] = genparticle.vz();
                GenLeptonSel_dxy[i] = dxy_value(genparticle, thePrimaryVertex); // should be computed here or before?

                GenLeptonSel_motherPdgId[i] = genparticle.motherRef()->pdgId();
            }
            
            // Flags
            GenLeptonSel_isPromptFinalState[i] = genparticle.isPromptFinalState();
            GenLeptonSel_fromHardProcessFinalState[i] = genparticle.fromHardProcessFinalState(); // has to be done with the last one
            GenLeptonSel_isDirectPromptTauDecayProductFinalState[i] = genparticle.isDirectPromptTauDecayProductFinalState(); 
            GenLeptonSel_isDirectHadronDecayProduct[i] = genparticle.statusFlags().isDirectHadronDecayProduct(); 
         }
         
         // Counters initialization
         nGenMuon_PFS = 0; 
         nGenMuon_HPFS = 0; 
         nGenMuon_PTDP = 0; 
         nGenMuon_HDP = 0; 
         for (size_t i = 0; i < iGM.size(); i++) {
            if (GenLeptonSel_isPromptFinalState[i]) { nGenMuon_PFS++; }
            if (GenLeptonSel_fromHardProcessFinalState[i]) { nGenMuon_HPFS++; }
            if (GenLeptonSel_isDirectPromptTauDecayProductFinalState[i]) { nGenMuon_PTDP++; }
            if (GenLeptonSel_isDirectHadronDecayProduct[i]) { nGenMuon_HDP++; }
         }

         // Hard Process Collection
         nHardProcessParticle = 0;
         for (size_t i = 0; i < genParticles->size(); i++) {
            const reco::GenParticle &genparticle = (*genParticles)[i];
            if (genparticle.isHardProcess()){
                HardProcessParticle_pt[nHardProcessParticle] = genparticle.pt();
                HardProcessParticle_E[nHardProcessParticle] = genparticle.energy();
                HardProcessParticle_eta[nHardProcessParticle] = genparticle.eta();
                HardProcessParticle_phi[nHardProcessParticle] = genparticle.phi();
                HardProcessParticle_vx[nHardProcessParticle] = genparticle.vx();
                HardProcessParticle_vy[nHardProcessParticle] = genparticle.vy();
                HardProcessParticle_vz[nHardProcessParticle] = genparticle.vz();
                HardProcessParticle_pdgId[nHardProcessParticle] = genparticle.pdgId();
                nHardProcessParticle++;
             }       
         }
     }
   }

   // Check if trigger fired:
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   unsigned int ipath = 0;
   for (auto path : HLTPaths_) {
     std::string path_v = path + "_v";
     // std::cout << path << "\t" << std::endl;
     bool fired = false;
     for (unsigned int itrg = 0; itrg < triggerBits->size(); ++itrg) {
       TString TrigPath = names.triggerName(itrg);
       if (!triggerBits->accept(itrg))
         continue;
       if (!TrigPath.Contains(path_v)){
         continue;
       }
       fired = true;
     }
     triggerPass[ipath] = fired;
     ipath++;
   } 

   //-> Fill tree
   tree_out->Fill();

}

DEFINE_FWK_MODULE(ntuplizer);
