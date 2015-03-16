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


class XsecScale : public edm::EDAnalyzer {
public:
  explicit XsecScale(const edm::ParameterSet&);
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
  const double scale;
  
  TH1F* h_HT;
  TH1F* h_HT_scaled;
  TH1F* h_gen_HT;
  TH1F* h_gen_HT_scaled;

};

XsecScale::XsecScale(const edm::ParameterSet& cfg)
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
    lljj_mass_cut(cfg.getParameter<double>("lljj_mass_cut")),
    scale(cfg.getParameter<double>("scale"))
{
  edm::Service<TFileService> fs;
  h_HT = fs->make<TH1F>("h_HT", "", 100, 0, 4000);
  h_HT_scaled = fs->make<TH1F>("h_HT_scaled", "", 100, 0, 4000);  
  h_gen_HT = fs->make<TH1F>("h_gen_HT", "", 100, 0, 4000);
  h_gen_HT_scaled = fs->make<TH1F>("h_gen_HT_scaled", "", 100, 0, 4000);  
    
}


void XsecScale::analyze(const edm::Event& event, const edm::EventSetup& setup) {
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

  std::vector<pat::Electron> eles;
  std::vector<pat::Muon> mus;
  std::vector<TLorentzVector> best_mus;
  std::vector<reco::GenParticle> gen_mus;
  std::vector<reco::GenParticle> gjets;
  std::vector<reco::GenParticle> gmus_W;  
  std::vector<pat::Jet> pjets;

  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel("offlineSlimmedPrimaryVertices", vertices);

  float ht = 0;
  float ht_gen = 0;

  for(auto jet:*jets)
    ht = ht + jet.pt();
  for(auto gjet:*gen_jets)
    ht_gen = ht_gen + gjet.pt();

  h_HT->Fill(ht);
  h_HT_scaled->Fill(ht,scale);  
  h_gen_HT->Fill(ht_gen);
  h_gen_HT_scaled->Fill(ht_gen,scale);  
  
}

DEFINE_FWK_MODULE(XsecScale);
