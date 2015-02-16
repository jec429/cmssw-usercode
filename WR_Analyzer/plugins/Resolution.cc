#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"


class Resolution : public edm::EDAnalyzer {
public:
  explicit Resolution(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

private:
  const edm::InputTag gen_src;
  const edm::InputTag vertex_src;
  const edm::InputTag gen_jet_src;
  const edm::InputTag electron_src;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronHEEPIdMapToken_;
  const edm::InputTag muon_src;
  const edm::InputTag jet_src;
  const double jet_pt_min;
  const double jet_eta_max;
  const double dilepton_mass_cut;
  const double lljj_mass_cut;

  TH1F* mu_pre0;
  TH1F* mu_pre1;
  TH1F* mu_pre2;
  TH1F* mu_pre3;
  TH1F* mu_pre4;
  TH1F* mu_pre5;
  TH1F* mu_pre6;

  TH1F* global_mu_pt;
  TH1F* best_mu_pt;
  TH1F* picky_mu_pt;
  TH1F* PF_mu_pt;
  TH1F* tpfms_mu_pt;
  TH1F* tuneP_mu_pt;

};

Resolution::Resolution(const edm::ParameterSet& cfg)
  : gen_src(cfg.getParameter<edm::InputTag>("gen_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    gen_jet_src(cfg.getParameter<edm::InputTag>("gen_jet_src")),
    electron_src(cfg.getParameter<edm::InputTag>("electron_src")),
    electronVetoIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("electronVetoIdMap"))),
    electronLooseIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("electronLooseIdMap"))),
    electronMediumIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("electronLooseIdMap"))),
    electronTightIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("electronTightIdMap"))),
    electronHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("electronHEEPIdMap"))),
    muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
    jet_src(cfg.getParameter<edm::InputTag>("jet_src")),
    jet_pt_min(cfg.getParameter<double>("jet_pt_min")),  
    jet_eta_max(cfg.getParameter<double>("jet_eta_max")),
    dilepton_mass_cut(cfg.getParameter<double>("dilepton_mass_cut")),
    lljj_mass_cut(cfg.getParameter<double>("lljj_mass_cut"))
{
  edm::Service<TFileService> fs;
  mu_pre0 = fs->make<TH1F>("mu_pre0", "", 2, 0, 2);
  mu_pre1 = fs->make<TH1F>("mu_pre1", "", 2, 0, 2);
  mu_pre2 = fs->make<TH1F>("mu_pre2", "", 2, 0, 2);
  mu_pre3 = fs->make<TH1F>("mu_pre3", "", 2, 0, 2);
  mu_pre4 = fs->make<TH1F>("mu_pre4", "", 2, 0, 2);
  mu_pre5 = fs->make<TH1F>("mu_pre5", "", 2, 0, 2);
  mu_pre6 = fs->make<TH1F>("mu_pre6", "", 2, 0, 2);

  global_mu_pt = fs->make<TH1F>("global_mu_pt", "", 200, -10, 10);
  best_mu_pt = fs->make<TH1F>("best_mu_pt", "", 200, -10, 10);
  picky_mu_pt = fs->make<TH1F>("picky_mu_pt", "", 200, -10, 10);
  PF_mu_pt = fs->make<TH1F>("PF_mu_pt", "", 200, -10, 10);
  tpfms_mu_pt = fs->make<TH1F>("tpfms_mu_pt", "", 200, -10, 10);
  tuneP_mu_pt = fs->make<TH1F>("tuneP_mu_pt", "", 200, -10, 10);
    
}


void Resolution::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  using namespace std;
  using namespace edm;
  using namespace reco;

  edm::Handle<pat::ElectronCollection> electrons;
  event.getByLabel(electron_src, electrons);

  edm::Handle<pat::MuonCollection> muons;
  event.getByLabel(muon_src, muons);
  edm::Handle<reco::VertexCollection> primary_vertex;
  event.getByLabel(vertex_src, primary_vertex);

  edm::Handle<pat::JetCollection> jets;
  event.getByLabel(jet_src, jets);
  edm::Handle<pat::METCollection> mets;
  event.getByLabel("slimmedMETs", mets);
  edm::Handle<reco::GenParticleCollection> gen_particles;
  event.getByLabel(gen_src, gen_particles);
  edm::Handle<reco::GenJetCollection> gen_jets;
  event.getByLabel(gen_jet_src, gen_jets);

  //std::cout<<"EVENT0"<<std::endl;
  std::vector<pat::Electron> eles;
  std::vector<pat::Muon> mus;
  std::vector<reco::GenParticle> gen_mus;
  std::vector<reco::GenJet> gjets;  
  std::vector<pat::Jet> pjets;

  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel("offlineSlimmedPrimaryVertices", vertices);

  for(auto mu : *muons){
    mu_pre0->Fill(1);
    bool MC_match = false;
    if(mu.pt() > 40){
      mu_pre1->Fill(1);
      if(fabs(mu.eta()) < 2.4){
	mu_pre2->Fill(1);
	
	float dR = 0.1;
	for(auto gen:*gen_particles){
	  if(fabs(gen.pdgId()) == 13 && deltaR(gen.eta(),gen.phi(),mu.eta(),mu.phi()) < dR && gen.status() == 1 && fabs(gen.mother()->pdgId()) != 15) {
	    MC_match = true;
	    dR = deltaR(gen.eta(),gen.phi(),mu.eta(),mu.phi());
	  }
	}
	if(MC_match){
	  mu_pre3->Fill(1);
	if(mu.isHighPtMuon(primary_vertex->at(0)))
	  mu_pre4->Fill(1);	  
	if(mu.isTightMuon(primary_vertex->at(0)))
	  mu_pre5->Fill(1);	
	}
      }
    }
    if((mu.pt() > 40) && (fabs(mu.eta()) < 2.4) && mu.isHighPtMuon(primary_vertex->at(0)) && MC_match)
      mus.push_back(mu);      
  }

  for(auto mu:mus){ 
    //std::cout<<mu.pt()<<std::endl;
    // Find the matched genparticle
    float dR = 0.1;
    std::vector<reco::GenParticle> gmus;
    for(auto gen:*gen_particles){
      if(fabs(gen.pdgId()) == 13 && deltaR(gen.eta(),gen.phi(),mu.eta(),mu.phi()) < dR && gen.status() == 1 && fabs(gen.mother()->pdgId()) != 15) {
	dR = deltaR(gen.eta(),gen.phi(),mu.eta(),mu.phi());
	gmus.push_back(gen);
      }
    }
    global_mu_pt->Fill(mu.globalTrack()->pt() - gmus.back().pt());
    best_mu_pt->Fill(mu.bestTrack()->pt() - gmus.back().pt());
    //picky_mu_pt->Fill(mu.pickyTrack()->pt() - gmus.back().pt());
    //PF_mu_pt->Fill(mu.pfCandidateRef()->pt() - gmus.back().pt());
    //tpfms_mu_pt->Fill(mu.tpfmsTrack()->pt() - gmus.back().pt());
    tuneP_mu_pt->Fill(mu.tunePMuonBestTrack()->pt() - gmus.back().pt());
  }
  

}

DEFINE_FWK_MODULE(Resolution);
