#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
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
      // trigger objects
      edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone> > triggerObjects_;
      edm::Handle<edm::View<pat::TriggerObjectStandAlone>  > triggerObjects;
      // trigger prescales
      edm::EDGetTokenT<pat::PackedTriggerPrescales>  triggerPrescales_;
      edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

      // displacedGlobalMuons (reco::Track)
      edm::EDGetTokenT<edm::View<reco::Track> > dglToken;
      edm::Handle<edm::View<reco::Track> > dgls;
      // displacedStandAloneMuons (reco::Track)
      edm::EDGetTokenT<edm::View<reco::Track> > dsaToken;
      edm::Handle<edm::View<reco::Track> > dsas;
      // displacedMuons (reco::Muon // pat::Muon)
      edm::EDGetTokenT<edm::View<pat::Muon> > dmuToken;
      edm::Handle<edm::View<pat::Muon> > dmuons;

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


      // displacedMuons
      Int_t ndmu = 0;
      Float_t dmu_pt[200] = {0.};
      Float_t dmu_eta[200] = {0.};
      Float_t dmu_phi[200] = {0.};
      Int_t dmu_isDSA[200] = {0};
      Int_t dmu_isDGL[200] = {0};
      Int_t dmu_isDTK[200] = {0};
      Float_t dmu_dsa_pt[200] = {0.};

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
   dmuToken = consumes<edm::View<pat::Muon> >  (parameters.getParameter<edm::InputTag>("displacedMuonCollection"));

   triggerBits_ = consumes<edm::TriggerResults> (parameters.getParameter<edm::InputTag>("bits"));
   triggerObjects_ = consumes<edm::View<pat::TriggerObjectStandAlone> > (parameters.getParameter<edm::InputTag>("objects"));
   triggerPrescales_ = consumes<pat::PackedTriggerPrescales > (parameters.getParameter<edm::InputTag>("prescales"));
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

   tree_out->Branch("ndmu", &ndmu, "ndmu/I");
   tree_out->Branch("dmu_pt", dmu_pt, "dmu_pt[ndmu]/F");
   tree_out->Branch("dmu_eta", dmu_eta, "dmu_eta[ndmu]/F");
   tree_out->Branch("dmu_phi", dmu_phi, "dmu_phi[ndmu]/F");
   tree_out->Branch("dmu_isDSA", dmu_isDSA, "dmu_isDSA[ndmu]/I");
   tree_out->Branch("dmu_isDGL", dmu_isDGL, "dmu_isDGL[ndmu]/I");
   tree_out->Branch("dmu_isDTK", dmu_isDTK, "dmu_isDTK[ndmu]/I");
   tree_out->Branch("dmu_dsa_pt", dmu_dsa_pt, "dmu_dsa_pt[ndmu]/F");

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
   iEvent.getByToken(triggerObjects_, triggerObjects);

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
     ndsa++;
   }

   // displacedMuons
   ndmu = 0;;
   for (unsigned int i = 0; i < dmuons->size(); i++) {
     const pat::Muon& dmuon(dmuons->at(i));
     dmu_pt[ndmu] = dmuon.pt();
     dmu_eta[ndmu] = dmuon.eta();
     dmu_phi[ndmu] = dmuon.phi();
     dmu_isDGL[ndmu] = dmuon.isGlobalMuon();
     dmu_isDSA[ndmu] = dmuon.isStandAloneMuon();
     dmu_isDTK[ndmu] = dmuon.isTrackerMuon();

     // Access the DSA track associated to the displacedMuon
     if ( dmuon.isStandAloneMuon() ) {
       const reco::Track* outerTrack = (dmuon.standAloneMuon()).get();
       dmu_dsa_pt[ndmu] = outerTrack->pt();
     }

     ndmu++;
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
