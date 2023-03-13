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

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"



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

      //
      // --- Tokens and Handles
      //

      // trigger bits
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::Handle<edm::TriggerResults> triggerBits;

      // displacedGlobalMuons (reco::Track)
      edm::EDGetTokenT<edm::View<reco::Track> > dglToken;
      edm::Handle<edm::View<reco::Track> > dgls;
      // displacedStandAloneMuons (reco::Track)
      edm::EDGetTokenT<edm::View<reco::Track> > dsaToken;
      edm::Handle<edm::View<reco::Track> > dsas;
      // displacedMuons (reco::Muon // pat::Muon)
      edm::EDGetTokenT<edm::View<reco::Muon> > dmuToken;
      edm::Handle<edm::View<reco::Muon> > dmuons;

      //
      // --- Variables
      //

      bool isData = false;

      // Trigger tags
      std::vector<std::string> HLTPaths_;
      bool triggerPass[200] = {false};

      // Event
      Int_t event = 0;
      Int_t lumiBlock = 0;
      Int_t run = 0;

      // displacedGlobalMuons
      Int_t ndgl = 0;
      Float_t dgl_pt[200] = {0.};
      Float_t dgl_eta[200] = {0.};
      Float_t dgl_phi[200] = {0.};
      Float_t dgl_ptError[200] = {0.};
      Float_t dgl_dxy[200] = {0.};
      Float_t dgl_dz[200] = {0.};
      Float_t dgl_normalizedChi2[200] = {0.};
      Float_t dgl_charge[200] = {0.};
      Int_t dgl_nMuonHits[200] = {0};
      Int_t dgl_nValidMuonHits[200] = {0};
      Int_t dgl_nValidMuonDTHits[200] = {0};
      Int_t dgl_nValidMuonCSCHits[200] = {0};
      Int_t dgl_nValidMuonRPCHits[200] = {0};
      Int_t dgl_nValidStripHits[200] = {0};
      Int_t dgl_nhits[200] = {0};
      Int_t dgl_nLostMuonHits[200] = {0};
      Int_t dgl_nLostMuonDTHits[200] = {0};
      Int_t dgl_nLostMuonCSCHits[200] = {0};
      Int_t dgl_nLostMuonRPCHits[200] = {0};

      // displacedStandAloneMuons
      Int_t ndsa = 0;
      Float_t dsa_pt[200] = {0.};
      Float_t dsa_eta[200] = {0.};
      Float_t dsa_phi[200] = {0.};
      Float_t dsa_ptError[200] = {0.};
      Float_t dsa_dxy[200] = {0.};
      Float_t dsa_dz[200] = {0.};
      Float_t dsa_normalizedChi2[200] = {0.};
      Float_t dsa_charge[200] = {0.};
      Int_t dsa_nMuonHits[200] = {0};
      Int_t dsa_nValidMuonHits[200] = {0};
      Int_t dsa_nValidMuonDTHits[200] = {0};
      Int_t dsa_nValidMuonCSCHits[200] = {0};
      Int_t dsa_nValidMuonRPCHits[200] = {0};
      Int_t dsa_nValidStripHits[200] = {0};
      Int_t dsa_nhits[200] = {0};
      Int_t dsa_nLostMuonHits[200] = {0};
      Int_t dsa_nLostMuonDTHits[200] = {0};
      Int_t dsa_nLostMuonCSCHits[200] = {0};
      Int_t dsa_nLostMuonRPCHits[200] = {0};
      Int_t dsa_dtStationsWithValidHits[200] = {0};
      Int_t dsa_cscStationsWithValidHits[200] = {0};


      // displacedMuons
      Int_t ndmu = 0;
      /*Float_t dmu_pt[200] = {0.};
      Float_t dmu_eta[200] = {0.};
      Float_t dmu_phi[200] = {0.};
      Float_t dmu_normChi2[200] = {0.};
      Float_t dmu_charge[200] = {0.};*/
      Int_t dmu_isDSA[200] = {0};
      Int_t dmu_isDGL[200] = {0};
      Int_t dmu_isDTK[200] = {0};
      Int_t dmu_isMatchesValid[200] = {0};
      Int_t dmu_numberOfMatches[200] = {0};
      Int_t dmu_numberOfChambers[200] = {0};
      Int_t dmu_numberOfChambersCSCorDT[200] = {0};
      Int_t dmu_numberOfMatchedStations[200] = {0};
      Int_t dmu_numberOfMatchedRPCLayers[200] = {0};
      Float_t dmu_dsa_pt[200] = {0.};
      Float_t dmu_dsa_eta[200] = {0.};
      Float_t dmu_dsa_phi[200] = {0.};
      Float_t dmu_dsa_ptError[200] = {0.};
      Float_t dmu_dsa_dxy[200] = {0.};
      Float_t dmu_dsa_dz[200] = {0.};
      Float_t dmu_dsa_normalizedChi2[200] = {0.};
      Float_t dmu_dsa_charge[200] = {0.};
      Int_t dmu_dsa_nMuonHits[200] = {0};
      Int_t dmu_dsa_nValidMuonHits[200] = {0};
      Int_t dmu_dsa_nValidMuonDTHits[200] = {0};
      Int_t dmu_dsa_nValidMuonCSCHits[200] = {0};
      Int_t dmu_dsa_nValidMuonRPCHits[200] = {0};
      Int_t dmu_dsa_nValidStripHits[200] = {0};
      Int_t dmu_dsa_nhits[200] = {0};
      Int_t dmu_dsa_dtStationsWithValidHits[200] = {0};
      Int_t dmu_dsa_cscStationsWithValidHits[200] = {0};
      Int_t dmu_dsa_nsegments[200] = {0};
      Float_t dmu_dgl_pt[200] = {0.};
      Float_t dmu_dgl_eta[200] = {0.};
      Float_t dmu_dgl_phi[200] = {0.};
      Float_t dmu_dgl_ptError[200] = {0.};
      Float_t dmu_dgl_dxy[200] = {0.};
      Float_t dmu_dgl_dz[200] = {0.};
      Float_t dmu_dgl_normalizedChi2[200] = {0.};
      Float_t dmu_dgl_charge[200] = {0.};
      Int_t dmu_dgl_nMuonHits[200] = {0};
      Int_t dmu_dgl_nValidMuonHits[200] = {0};
      Int_t dmu_dgl_nValidMuonDTHits[200] = {0};
      Int_t dmu_dgl_nValidMuonCSCHits[200] = {0};
      Int_t dmu_dgl_nValidMuonRPCHits[200] = {0};
      Int_t dmu_dgl_nValidStripHits[200] = {0};
      Int_t dmu_dgl_nhits[200] = {0};
      Float_t dmu_dtk_pt[200] = {0.};
      Float_t dmu_dtk_eta[200] = {0.};
      Float_t dmu_dtk_phi[200] = {0.};
      Float_t dmu_dtk_ptError[200] = {0.};
      Float_t dmu_dtk_dxy[200] = {0.};
      Float_t dmu_dtk_dz[200] = {0.};
      Float_t dmu_dtk_normalizedChi2[200] = {0.};
      Float_t dmu_dtk_charge[200] = {0.};
      Int_t dmu_dtk_nMuonHits[200] = {0};
      Int_t dmu_dtk_nValidMuonHits[200] = {0};
      Int_t dmu_dtk_nValidMuonDTHits[200] = {0};
      Int_t dmu_dtk_nValidMuonCSCHits[200] = {0};
      Int_t dmu_dtk_nValidMuonRPCHits[200] = {0};
      Int_t dmu_dtk_nValidStripHits[200] = {0};
      Int_t dmu_dtk_nhits[200] = {0};

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

   dglToken = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("displacedGlobalCollection"));
   dsaToken = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("displacedStandAloneCollection"));
   dmuToken = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("displacedMuonCollection"));

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
   HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX3BX");
   HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX");

   // TTree branches
   tree_out->Branch("event", &event, "event/I");
   tree_out->Branch("lumiBlock", &lumiBlock, "lumiBlock/I");
   tree_out->Branch("run", &run, "run/I");

   tree_out->Branch("ndgl", &ndgl, "ndgl/I");
   tree_out->Branch("dgl_pt", dgl_pt, "dgl_pt[ndgl]/F");
   tree_out->Branch("dgl_eta", dgl_eta, "dgl_eta[ndgl]/F");
   tree_out->Branch("dgl_phi", dgl_phi, "dgl_phi[ndgl]/F");
   tree_out->Branch("dgl_ptError", dgl_ptError, "dgl_ptError[ndgl]/F");
   tree_out->Branch("dgl_dxy", dgl_dxy, "dgl_dxy[ndgl]/F");
   tree_out->Branch("dgl_dz", dgl_dz, "dgl_dz[ndgl]/F");
   tree_out->Branch("dgl_normalizedChi2", dgl_normalizedChi2, "dgl_normalizedChi2[ndgl]/F");
   tree_out->Branch("dgl_charge", dgl_charge, "dgl_charge[ndgl]/F");
   tree_out->Branch("dgl_nMuonHits", dgl_nMuonHits, "dgl_nMuonHits[ndgl]/I");
   tree_out->Branch("dgl_nValidMuonHits", dgl_nValidMuonHits, "dgl_nValidMuonHits[ndgl]/I");
   tree_out->Branch("dgl_nValidMuonDTHits", dgl_nValidMuonDTHits, "dgl_nValidMuonDTHits[ndgl]/I");
   tree_out->Branch("dgl_nValidMuonCSCHits", dgl_nValidMuonCSCHits, "dgl_nValidMuonCSCHits[ndgl]/I");
   tree_out->Branch("dgl_nValidMuonRPCHits", dgl_nValidMuonRPCHits, "dgl_nValidMuonRPCHits[ndgl]/I");
   tree_out->Branch("dgl_nValidStripHits", dgl_nValidStripHits, "dgl_nValidStripHits[ndgl]/I");
   tree_out->Branch("dgl_nhits", dgl_nhits, "dgl_nhits[ndgl]/I");
   tree_out->Branch("dgl_nLostMuonHits", dgl_nLostMuonHits, "dgl_nLostMuonHits[ndgl]/I");
   tree_out->Branch("dgl_nLostMuonDTHits", dgl_nLostMuonDTHits, "dgl_nLostMuonDTHits[ndgl]/I");
   tree_out->Branch("dgl_nLostMuonCSCHits", dgl_nLostMuonCSCHits, "dgl_nLostMuonCSCHits[ndgl]/I");
   tree_out->Branch("dgl_nLostMuonRPCHits", dgl_nLostMuonRPCHits, "dgl_nLostMuonRPCHits[ndgl]/I");

   tree_out->Branch("ndsa", &ndsa, "ndsa/I");
   tree_out->Branch("dsa_pt", dsa_pt, "dsa_pt[ndsa]/F");
   tree_out->Branch("dsa_eta", dsa_eta, "dsa_eta[ndsa]/F");
   tree_out->Branch("dsa_phi", dsa_phi, "dsa_phi[ndsa]/F");
   tree_out->Branch("dsa_ptError", dsa_ptError, "dsa_ptError[ndsa]/F");
   tree_out->Branch("dsa_dxy", dsa_dxy, "dsa_dxy[ndsa]/F");
   tree_out->Branch("dsa_dz", dsa_dz, "dsa_dz[ndsa]/F");
   tree_out->Branch("dsa_normalizedChi2", dsa_normalizedChi2, "dsa_normalizedChi2[ndsa]/F");
   tree_out->Branch("dsa_charge", dsa_charge, "dsa_charge[ndsa]/F");
   tree_out->Branch("dsa_nMuonHits", dsa_nMuonHits, "dsa_nMuonHits[ndsa]/I");
   tree_out->Branch("dsa_nValidMuonHits", dsa_nValidMuonHits, "dsa_nValidMuonHits[ndsa]/I");
   tree_out->Branch("dsa_nValidMuonDTHits", dsa_nValidMuonDTHits, "dsa_nValidMuonDTHits[ndsa]/I");
   tree_out->Branch("dsa_nValidMuonCSCHits", dsa_nValidMuonCSCHits, "dsa_nValidMuonCSCHits[ndsa]/I");
   tree_out->Branch("dsa_nValidMuonRPCHits", dsa_nValidMuonRPCHits, "dsa_nValidMuonRPCHits[ndsa]/I");
   tree_out->Branch("dsa_nValidStripHits", dsa_nValidStripHits, "dsa_nValidStripHits[ndsa]/I");
   tree_out->Branch("dsa_nhits", dsa_nhits, "dsa_nhits[ndsa]/I");
   tree_out->Branch("dsa_nLostMuonHits", dsa_nLostMuonHits, "dsa_nLostMuonHits[ndsa]/I");
   tree_out->Branch("dsa_nLostMuonDTHits", dsa_nLostMuonDTHits, "dsa_nLostMuonDTHits[ndsa]/I");
   tree_out->Branch("dsa_nLostMuonCSCHits", dsa_nLostMuonCSCHits, "dsa_nLostMuonCSCHits[ndsa]/I");
   tree_out->Branch("dsa_nLostMuonRPCHits", dsa_nLostMuonRPCHits, "dsa_nLostMuonRPCHits[ndsa]/I");
   tree_out->Branch("dsa_dtStationsWithValidHits", dsa_dtStationsWithValidHits, "dsa_dtStationsWithValidHits[ndsa]/I");
   tree_out->Branch("dsa_cscStationsWithValidHits", dsa_cscStationsWithValidHits, "dsa_cscStationsWithValidHits[ndsa]/I");

   tree_out->Branch("ndmu", &ndmu, "ndmu/I");
   /*tree_out->Branch("dmu_pt", dmu_pt, "dmu_pt[ndmu]/F");
   tree_out->Branch("dmu_eta", dmu_eta, "dmu_eta[ndmu]/F");
   tree_out->Branch("dmu_phi", dmu_phi, "dmu_phi[ndmu]/F");
   tree_out->Branch("dmu_normChi2", dmu_normChi2, "dmu_normChi2[ndmu]/F");
   tree_out->Branch("dmu_charge", dmu_charge, "dmu_charge[ndmu]/F");
   tree_out->Branch("dmu_nMuonHits", dmu_nMuonHits, "dmu_nMuonHits[ndmu]/I");
   tree_out->Branch("dmu_nValidMuonHits", dmu_nValidMuonHits, "dmu_nValidMuonHits[ndmu]/I");
   tree_out->Branch("dmu_nValidMuonDTHits", dmu_nValidMuonDTHits, "dmu_nValidMuonDTHits[ndmu]/I");
   tree_out->Branch("dmu_nValidMuonCSCHits", dmu_nValidMuonCSCHits, "dmu_nValidMuonCSCHits[ndmu]/I");
   tree_out->Branch("dmu_nValidMuonRPCHits", dmu_nValidMuonRPCHits, "dmu_nValidMuonRPCHits[ndmu]/I");
   tree_out->Branch("dmu_nValidStripHits", dmu_nValidStripHits, "dmu_nValidStripHits[ndmu]/I");
   tree_out->Branch("dmu_nhits", dmu_nhits, "dmu_nhits[ndmu]/I");*/
   tree_out->Branch("dmu_isDSA", dmu_isDSA, "dmu_isDSA[ndmu]/I");
   tree_out->Branch("dmu_isDGL", dmu_isDGL, "dmu_isDGL[ndmu]/I");
   tree_out->Branch("dmu_isDTK", dmu_isDTK, "dmu_isDTK[ndmu]/I");
   tree_out->Branch("dmu_isMatchesValid", dmu_isMatchesValid, "dmu_isMatchesValid[ndmu]/I");
   tree_out->Branch("dmu_numberOfMatches", dmu_numberOfMatches, "dmu_numberOfMatches[ndmu]/I");
   tree_out->Branch("dmu_numberOfChambers", dmu_numberOfChambers, "dmu_numberOfChambers[ndmu]/I");
   tree_out->Branch("dmu_numberOfChambersCSCorDT", dmu_numberOfChambersCSCorDT, "dmu_numberOfChambersCSCorDT[ndmu]/I");
   tree_out->Branch("dmu_numberOfMatchedStations", dmu_numberOfMatchedStations, "dmu_numberOfMatchedStations[ndmu]/I");
   tree_out->Branch("dmu_numberOfMatchedRPCLayers", dmu_numberOfMatchedRPCLayers, "dmu_numberOfMatchedRPCLayers[ndmu]/I");
   // dmu_dsa
   tree_out->Branch("dmu_dsa_pt", dmu_dsa_pt, "dmu_dsa_pt[ndmu]/F");
   tree_out->Branch("dmu_dsa_eta", dmu_dsa_eta, "dmu_dsa_eta[ndmu]/F");
   tree_out->Branch("dmu_dsa_phi", dmu_dsa_phi, "dmu_dsa_phi[ndmu]/F");
   tree_out->Branch("dmu_dsa_ptError", dmu_dsa_ptError, "dmu_dsa_ptError[ndmu]/F");
   tree_out->Branch("dmu_dsa_dxy", dmu_dsa_dxy, "dmu_dsa_dxy[ndmu]/F");
   tree_out->Branch("dmu_dsa_dz", dmu_dsa_dz, "dmu_dsa_dz[ndmu]/F");
   tree_out->Branch("dmu_dsa_normalizedChi2", dmu_dsa_normalizedChi2, "dmu_dsa_normalizedChi2[ndmu]/F");
   tree_out->Branch("dmu_dsa_charge", dmu_dsa_charge, "dmu_dsa_charge[ndmu]/F");
   tree_out->Branch("dmu_dsa_nMuonHits", dmu_dsa_nMuonHits, "dmu_dsa_nMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonHits", dmu_dsa_nValidMuonHits, "dmu_dsa_nValidMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonDTHits", dmu_dsa_nValidMuonDTHits, "dmu_dsa_nValidMuonDTHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonCSCHits", dmu_dsa_nValidMuonCSCHits, "dmu_dsa_nValidMuonCSCHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonRPCHits", dmu_dsa_nValidMuonRPCHits, "dmu_dsa_nValidMuonRPCHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidStripHits", dmu_dsa_nValidStripHits, "dmu_dsa_nValidStripHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nhits", dmu_dsa_nhits, "dmu_dsa_nhits[ndmu]/I");
   tree_out->Branch("dmu_dsa_dtStationsWithValidHits", dmu_dsa_dtStationsWithValidHits, "dmu_dsa_dtStationsWithValidHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_cscStationsWithValidHits", dmu_dsa_cscStationsWithValidHits, "dmu_dsa_cscStationsWithValidHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nsegments", dmu_dsa_nsegments, "dmu_dsa_nsegments[ndmu]/I");
   // dmu_dgl
   tree_out->Branch("dmu_dgl_pt", dmu_dgl_pt, "dmu_dgl_pt[ndmu]/F");
   tree_out->Branch("dmu_dgl_eta", dmu_dgl_eta, "dmu_dgl_eta[ndmu]/F");
   tree_out->Branch("dmu_dgl_phi", dmu_dgl_phi, "dmu_dgl_phi[ndmu]/F");
   tree_out->Branch("dmu_dgl_ptError", dmu_dgl_ptError, "dmu_dgl_ptError[ndmu]/F");
   tree_out->Branch("dmu_dgl_dxy", dmu_dgl_dxy, "dmu_dgl_dxy[ndmu]/F");
   tree_out->Branch("dmu_dgl_dz", dmu_dgl_dz, "dmu_dgl_dz[ndmu]/F");
   tree_out->Branch("dmu_dgl_normalizedChi2", dmu_dgl_normalizedChi2, "dmu_dgl_normalizedChi2[ndmu]/F");
   tree_out->Branch("dmu_dgl_charge", dmu_dgl_charge, "dmu_dgl_charge[ndmu]/F");
   tree_out->Branch("dmu_dgl_nMuonHits", dmu_dgl_nMuonHits, "dmu_dgl_nMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonHits", dmu_dgl_nValidMuonHits, "dmu_dgl_nValidMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonDTHits", dmu_dgl_nValidMuonDTHits, "dmu_dgl_nValidMuonDTHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonCSCHits", dmu_dgl_nValidMuonCSCHits, "dmu_dgl_nValidMuonCSCHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonRPCHits", dmu_dgl_nValidMuonRPCHits, "dmu_dgl_nValidMuonRPCHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidStripHits", dmu_dgl_nValidStripHits, "dmu_dgl_nValidStripHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nhits", dmu_dgl_nhits, "dmu_dgl_nhits[ndmu]/I");
   // dmu_dtk
   tree_out->Branch("dmu_dtk_pt", dmu_dtk_pt, "dmu_dtk_pt[ndmu]/F");
   tree_out->Branch("dmu_dtk_eta", dmu_dtk_eta, "dmu_dtk_eta[ndmu]/F");
   tree_out->Branch("dmu_dtk_phi", dmu_dtk_phi, "dmu_dtk_phi[ndmu]/F");
   tree_out->Branch("dmu_dtk_ptError", dmu_dtk_ptError, "dmu_dtk_ptError[ndmu]/F");
   tree_out->Branch("dmu_dtk_dxy", dmu_dtk_dxy, "dmu_dtk_dxy[ndmu]/F");
   tree_out->Branch("dmu_dtk_dz", dmu_dtk_dz, "dmu_dtk_dz[ndmu]/F");
   tree_out->Branch("dmu_dtk_normalizedChi2", dmu_dtk_normalizedChi2, "dmu_dtk_normalizedChi2[ndmu]/F");
   tree_out->Branch("dmu_dtk_charge", dmu_dtk_charge, "dmu_dtk_charge[ndmu]/F");
   tree_out->Branch("dmu_dtk_nMuonHits", dmu_dtk_nMuonHits, "dmu_dtk_nMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nValidMuonHits", dmu_dtk_nValidMuonHits, "dmu_dtk_nValidMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nValidMuonDTHits", dmu_dtk_nValidMuonDTHits, "dmu_dtk_nValidMuonDTHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nValidMuonCSCHits", dmu_dtk_nValidMuonCSCHits, "dmu_dtk_nValidMuonCSCHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nValidMuonRPCHits", dmu_dtk_nValidMuonRPCHits, "dmu_dtk_nValidMuonRPCHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nValidStripHits", dmu_dtk_nValidStripHits, "dmu_dtk_nValidStripHits[ndmu]/I");
   tree_out->Branch("dmu_dtk_nhits", dmu_dtk_nhits, "dmu_dtk_nhits[ndmu]/I");

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

   iEvent.getByToken(dglToken, dgls);
   iEvent.getByToken(dsaToken, dsas);
   iEvent.getByToken(dmuToken, dmuons);
   iEvent.getByToken(triggerBits_, triggerBits);

   // Count number of events read
   counts->Fill(0);


   // -> Event info
   event = iEvent.id().event();
   lumiBlock = iEvent.id().luminosityBlock();
   run = iEvent.id().run();


   // displacedGlobalMuons
   ndgl = 0;
   for (unsigned int i = 0; i < dgls->size(); i++) {
     const reco::Track& dgl(dgls->at(i));
     dgl_pt[ndgl] = dgl.pt();
     dgl_eta[ndgl] = dgl.eta();
     dgl_phi[ndgl] = dgl.phi();
     dgl_ptError[ndgl] = dgl.ptError();
     dgl_dxy[ndgl] = dgl.dxy();
     dgl_dz[ndgl] = dgl.dz();
     dgl_normalizedChi2[ndgl] = dgl.normalizedChi2();
     dgl_charge[ndgl] = dgl.charge();

     dgl_nMuonHits[ndgl] = dgl.hitPattern().numberOfMuonHits();
     dgl_nValidMuonDTHits[ndgl] = dgl.hitPattern().numberOfValidMuonDTHits();
     dgl_nValidMuonCSCHits[ndgl] = dgl.hitPattern().numberOfValidMuonCSCHits();
     dgl_nValidMuonRPCHits[ndgl] = dgl.hitPattern().numberOfValidMuonRPCHits();
     dgl_nValidMuonHits[ndgl] = dgl.hitPattern().numberOfValidMuonHits();
     dgl_nValidStripHits[ndgl] = dgl.hitPattern().numberOfValidStripHits();
     dgl_nhits[ndgl] = dgl.hitPattern().numberOfValidHits();

     dgl_nLostMuonHits[ndgl] = dgl.hitPattern().numberOfLostMuonHits();
     dgl_nLostMuonDTHits[ndgl] = dgl.hitPattern().numberOfLostMuonDTHits();
     dgl_nLostMuonCSCHits[ndgl] = dgl.hitPattern().numberOfLostMuonCSCHits();
     dgl_nLostMuonRPCHits[ndgl] = dgl.hitPattern().numberOfLostMuonRPCHits();
     ndgl++;
   }

   // displacedStandAloneMuons
   ndsa = 0;
   for (unsigned int i = 0; i < dsas->size(); i++) {
     const reco::Track& dsa(dsas->at(i));
     dsa_pt[ndsa] = dsa.pt();
     dsa_eta[ndsa] = dsa.eta();
     dsa_phi[ndsa] = dsa.phi();
     dsa_ptError[ndsa] = dsa.ptError();
     dsa_dxy[ndsa] = dsa.dxy();
     dsa_dz[ndsa] = dsa.dz();
     dsa_normalizedChi2[ndsa] = dsa.normalizedChi2();
     dsa_charge[ndsa] = dsa.charge();

     dsa_nMuonHits[ndsa] = dsa.hitPattern().numberOfMuonHits();
     dsa_nValidMuonHits[ndsa] = dsa.hitPattern().numberOfValidMuonHits();
     dsa_nValidMuonDTHits[ndsa] = dsa.hitPattern().numberOfValidMuonDTHits();
     dsa_nValidMuonCSCHits[ndsa] = dsa.hitPattern().numberOfValidMuonCSCHits();
     dsa_nValidMuonRPCHits[ndsa] = dsa.hitPattern().numberOfValidMuonRPCHits();
     dsa_nValidStripHits[ndsa] = dsa.hitPattern().numberOfValidStripHits();
     dsa_nhits[ndsa] = dsa.hitPattern().numberOfValidHits();

     dsa_nLostMuonHits[ndsa] = dsa.hitPattern().numberOfLostMuonHits();
     dsa_nLostMuonDTHits[ndsa] = dsa.hitPattern().numberOfLostMuonDTHits();
     dsa_nLostMuonCSCHits[ndsa] = dsa.hitPattern().numberOfLostMuonCSCHits();
     dsa_nLostMuonRPCHits[ndsa] = dsa.hitPattern().numberOfLostMuonRPCHits();

     dsa_dtStationsWithValidHits[ndsa] = dsa.hitPattern().dtStationsWithValidHits();
     dsa_cscStationsWithValidHits[ndsa] = dsa.hitPattern().cscStationsWithValidHits();

     ndsa++;
   }

   // displacedMuons
   ndmu = 0;;
   for (unsigned int i = 0; i < dmuons->size(); i++) {
     std::cout << " - - ndmu: " << ndmu << std::endl;
     const reco::Muon& dmuon(dmuons->at(i));
     //dmu_pt[ndmu] = dmuon.pt();
     //dmu_eta[ndmu] = dmuon.eta();
     //dmu_phi[ndmu] = dmuon.phi();
     //dmu_normChi2[ndmu] = dmuon.normChi2();
     //dmu_charge[ndmu] = dmuon.charge();
     dmu_isDGL[ndmu] = dmuon.isGlobalMuon();
     dmu_isDSA[ndmu] = dmuon.isStandAloneMuon();
     dmu_isDTK[ndmu] = dmuon.isTrackerMuon();
     dmu_isMatchesValid[ndmu] = dmuon.isMatchesValid();
     dmu_numberOfMatches[ndmu] = dmuon.numberOfMatches();
     dmu_numberOfChambers[ndmu] = dmuon.numberOfChambers();
     dmu_numberOfChambersCSCorDT[ndmu] = dmuon.numberOfChambersCSCorDT();
     dmu_numberOfMatchedStations[ndmu] = dmuon.numberOfMatchedStations();
     dmu_numberOfMatchedRPCLayers[ndmu] = dmuon.numberOfMatchedRPCLayers();
     /*
     dmu_nMuonHits[ndmu] = dmuon.hitPattern().numberOfMuonHits();
     dmu_nValidMuonHits[ndmu] = dmuon.hitPattern().numberOfValidMuonHits();
     dmu_nValidMuonDTHits[ndmu] = dmuon.hitPattern().numberOfValidMuonDTHits();
     dmu_nValidMuonCSCHits[ndmu] = dmuon.hitPattern().numberOfValidMuonCSCHits();
     dmu_nValidMuonRPCHits[ndmu] = dmuon.hitPattern().numberOfValidMuonRPCHits();
     dmu_nValidStripHits[ndmu] = dmuon.hitPattern().numberOfValidStripHits();
     dmu_nhits[ndmu] = dmuon.hitPattern().numberOfValidHits();
     */

     // Access the DGL track associated to the displacedMuon
     std::cout << "isGlobalMuon: " << dmuon.isGlobalMuon() << std::endl;
     if ( dmuon.isGlobalMuon() ) {
       const reco::Track* globalTrack = (dmuon.combinedMuon()).get();
       dmu_dgl_pt[ndmu] = globalTrack->pt();
       dmu_dgl_eta[ndmu] = globalTrack->eta();
       dmu_dgl_phi[ndmu] = globalTrack->phi();
       dmu_dgl_ptError[ndmu] = globalTrack->ptError();
       dmu_dgl_dxy[ndmu] = globalTrack->dxy();
       dmu_dgl_dz[ndmu] = globalTrack->dz();
       dmu_dgl_normalizedChi2[ndmu] = globalTrack->normalizedChi2();
       dmu_dgl_charge[ndmu] = globalTrack->charge();
       dmu_dgl_nMuonHits[ndmu] = globalTrack->hitPattern().numberOfMuonHits();
       dmu_dgl_nValidMuonHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonHits();
       dmu_dgl_nValidMuonDTHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonDTHits();
       dmu_dgl_nValidMuonCSCHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonCSCHits();
       dmu_dgl_nValidMuonRPCHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonRPCHits();
       dmu_dgl_nValidStripHits[ndmu] = globalTrack->hitPattern().numberOfValidStripHits();
       dmu_dgl_nhits[ndmu] = globalTrack->hitPattern().numberOfValidHits();
     } else {
       dmu_dgl_pt[ndmu] = 0;
       dmu_dgl_eta[ndmu] = 0;
       dmu_dgl_phi[ndmu] = 0;
       dmu_dgl_ptError[ndmu] = 0;
       dmu_dgl_dxy[ndmu] = 0;
       dmu_dgl_dz[ndmu] = 0;
       dmu_dgl_normalizedChi2[ndmu] = 0;
       dmu_dgl_charge[ndmu] = 0;
       dmu_dgl_nMuonHits[ndmu] = 0;
       dmu_dgl_nValidMuonHits[ndmu] = 0;
       dmu_dgl_nValidMuonDTHits[ndmu] = 0;
       dmu_dgl_nValidMuonCSCHits[ndmu] = 0;
       dmu_dgl_nValidMuonRPCHits[ndmu] = 0;
       dmu_dgl_nValidStripHits[ndmu] = 0;
       dmu_dgl_nhits[ndmu] = 0;
     }     

     // Access the DSA track associated to the displacedMuon
     std::cout << "isStandAloneMuon: " << dmuon.isStandAloneMuon() << std::endl;
     if ( dmuon.isStandAloneMuon() ) {
       const reco::Track* outerTrack = (dmuon.standAloneMuon()).get();
       dmu_dsa_pt[ndmu] = outerTrack->pt();
       dmu_dsa_eta[ndmu] = outerTrack->eta();
       dmu_dsa_phi[ndmu] = outerTrack->phi();
       dmu_dsa_ptError[ndmu] = outerTrack->ptError();
       dmu_dsa_dxy[ndmu] = outerTrack->dxy();
       dmu_dsa_dz[ndmu] = outerTrack->dz();
       dmu_dsa_normalizedChi2[ndmu] = outerTrack->normalizedChi2();
       dmu_dsa_charge[ndmu] = outerTrack->charge();
       dmu_dsa_nMuonHits[ndmu] = outerTrack->hitPattern().numberOfMuonHits();
       dmu_dsa_nValidMuonHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonHits();
       dmu_dsa_nValidMuonDTHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonDTHits();
       dmu_dsa_nValidMuonCSCHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonCSCHits();
       dmu_dsa_nValidMuonRPCHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonRPCHits();
       dmu_dsa_nValidStripHits[ndmu] = outerTrack->hitPattern().numberOfValidStripHits();
       dmu_dsa_nhits[ndmu] = outerTrack->hitPattern().numberOfValidHits();
       dmu_dsa_dtStationsWithValidHits[ndmu] = outerTrack->hitPattern().dtStationsWithValidHits();
       dmu_dsa_cscStationsWithValidHits[ndmu] = outerTrack->hitPattern().cscStationsWithValidHits();
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
       dmu_dsa_nsegments[ndmu] = nsegments;
       
     } else {
       dmu_dsa_pt[ndmu] = 0;
       dmu_dsa_eta[ndmu] = 0;
       dmu_dsa_phi[ndmu] = 0;
       dmu_dsa_ptError[ndmu] = 0;
       dmu_dsa_dxy[ndmu] = 0;
       dmu_dsa_dz[ndmu] = 0;
       dmu_dsa_normalizedChi2[ndmu] = 0;
       dmu_dsa_charge[ndmu] = 0;
       dmu_dsa_nMuonHits[ndmu] = 0;
       dmu_dsa_nValidMuonHits[ndmu] = 0;
       dmu_dsa_nValidMuonDTHits[ndmu] = 0;
       dmu_dsa_nValidMuonCSCHits[ndmu] = 0;
       dmu_dsa_nValidMuonRPCHits[ndmu] = 0;
       dmu_dsa_nValidStripHits[ndmu] = 0;
       dmu_dsa_nhits[ndmu] = 0;
       dmu_dsa_dtStationsWithValidHits[ndsa] = 0;
       dmu_dsa_cscStationsWithValidHits[ndsa] = 0;
       dmu_dsa_nsegments[ndmu] = 0;
     }

     // Access the DTK track associated to the displacedMuon
     std::cout << "isTrackerMuon: " << dmuon.isTrackerMuon() << std::endl;
     if ( dmuon.isTrackerMuon() ) {
       const reco::Track* innerTrack = (dmuon.track()).get();
       dmu_dtk_pt[ndmu] = innerTrack->pt();
       dmu_dtk_eta[ndmu] = innerTrack->eta();
       dmu_dtk_phi[ndmu] = innerTrack->phi();
       dmu_dtk_ptError[ndmu] = innerTrack->ptError();
       dmu_dtk_dxy[ndmu] = innerTrack->dxy();
       dmu_dtk_dz[ndmu] = innerTrack->dz();
       dmu_dtk_normalizedChi2[ndmu] = innerTrack->normalizedChi2();
       dmu_dtk_charge[ndmu] = innerTrack->charge();
       dmu_dtk_nMuonHits[ndmu] = innerTrack->hitPattern().numberOfMuonHits();
       dmu_dtk_nValidMuonHits[ndmu] = innerTrack->hitPattern().numberOfValidMuonHits();
       dmu_dtk_nValidMuonDTHits[ndmu] = innerTrack->hitPattern().numberOfValidMuonDTHits();
       dmu_dtk_nValidMuonCSCHits[ndmu] = innerTrack->hitPattern().numberOfValidMuonCSCHits();
       dmu_dtk_nValidMuonRPCHits[ndmu] = innerTrack->hitPattern().numberOfValidMuonRPCHits();
       dmu_dtk_nValidStripHits[ndmu] = innerTrack->hitPattern().numberOfValidStripHits();
       dmu_dtk_nhits[ndmu] = innerTrack->hitPattern().numberOfValidHits();
     } else {
       dmu_dtk_pt[ndmu] = 0;
       dmu_dtk_eta[ndmu] = 0;
       dmu_dtk_phi[ndmu] = 0;
       dmu_dtk_ptError[ndmu] = 0;
       dmu_dtk_dxy[ndmu] = 0;
       dmu_dtk_dz[ndmu] = 0;
       dmu_dtk_normalizedChi2[ndmu] = 0;
       dmu_dtk_charge[ndmu] = 0;
       dmu_dtk_nMuonHits[ndmu] = 0;
       dmu_dtk_nValidMuonHits[ndmu] = 0;
       dmu_dtk_nValidMuonDTHits[ndmu] = 0;
       dmu_dtk_nValidMuonCSCHits[ndmu] = 0;
       dmu_dtk_nValidMuonRPCHits[ndmu] = 0;
       dmu_dtk_nValidStripHits[ndmu] = 0;
       dmu_dtk_nhits[ndmu] = 0;
     }

     ndmu++;
     std::cout << "End muon" << std::endl;
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
       // std::cout << TrigPath << std::endl; 
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
