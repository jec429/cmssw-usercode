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

  TH1F* global_mu_pt1;
  TH1F* best_mu_pt1;
  TH1F* picky_mu_pt1;
  TH1F* PF_mu_pt1;
  TH1F* tpfms_mu_pt1;
  TH1F* tuneP_mu_pt1;
  TH1F* global_mu_pt2;
  TH1F* best_mu_pt2;
  TH1F* picky_mu_pt2;
  TH1F* PF_mu_pt2;
  TH1F* tpfms_mu_pt2;
  TH1F* tuneP_mu_pt2;
  TH1F* global_mu_pt3;
  TH1F* best_mu_pt3;
  TH1F* picky_mu_pt3;
  TH1F* PF_mu_pt3;
  TH1F* tpfms_mu_pt3;
  TH1F* tuneP_mu_pt3;

  TH1F* global_mu_pt4;
  TH1F* best_mu_pt4;
  TH1F* tuneP_mu_pt4;
  TH1F* global_mu_pt5;
  TH1F* best_mu_pt5;
  TH1F* tuneP_mu_pt5;
  TH1F* global_mu_pt6;
  TH1F* best_mu_pt6;
  TH1F* tuneP_mu_pt6;

  TH1F* iso03;
  TH1F* iso05;
  TH1F* pf_iso03;
  TH1F* pf_iso04;

  TH1F* Mjjll;
  TH1F* Mqqll;
  TH1F* Mjjlglg;
  TH1F* Mqqlglg;
  
  TH2F* typeDiff_pt;
  TH2F* typeDiff;
  TH2F* typeDiff_all;

  //TH2F* global_mu_pt;
  //TH2F* best_mu_pt;
  //TH2F* picky_mu_pt;
  //TH2F* PF_mu_pt;
  //TH2F* tpfms_mu_pt;
  //TH2F* tuneP_mu_pt;

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

  global_mu_pt1 = fs->make<TH1F>("global_mu_pt1", "", 200, -1, 1);
  best_mu_pt1 = fs->make<TH1F>("best_mu_pt1", "", 200, -1, 1);
  picky_mu_pt1 = fs->make<TH1F>("picky_mu_pt1", "", 200, -1, 1);
  PF_mu_pt1 = fs->make<TH1F>("PF_mu_pt1", "", 200, -1, 1);
  tpfms_mu_pt1 = fs->make<TH1F>("tpfms_mu_pt1", "", 200, -1, 1);
  tuneP_mu_pt1 = fs->make<TH1F>("tuneP_mu_pt1", "", 200, -1, 1);
  global_mu_pt2 = fs->make<TH1F>("global_mu_pt2", "", 200, -1, 1);
  best_mu_pt2 = fs->make<TH1F>("best_mu_pt2", "", 200, -1, 1);
  picky_mu_pt2 = fs->make<TH1F>("picky_mu_pt2", "", 200, -1, 1);
  PF_mu_pt2 = fs->make<TH1F>("PF_mu_pt2", "", 200, -1, 1);
  tpfms_mu_pt2 = fs->make<TH1F>("tpfms_mu_pt2", "", 200, -1, 1);
  tuneP_mu_pt2 = fs->make<TH1F>("tuneP_mu_pt2", "", 200, -1, 1);
  global_mu_pt3 = fs->make<TH1F>("global_mu_pt3", "", 200, -1, 1);
  best_mu_pt3 = fs->make<TH1F>("best_mu_pt3", "", 200, -1, 1);
  picky_mu_pt3 = fs->make<TH1F>("picky_mu_pt3", "", 200, -1, 1);
  PF_mu_pt3 = fs->make<TH1F>("PF_mu_pt3", "", 200, -1, 1);
  tpfms_mu_pt3 = fs->make<TH1F>("tpfms_mu_pt3", "", 200, -1, 1);
  tuneP_mu_pt3 = fs->make<TH1F>("tuneP_mu_pt3", "", 200, -1, 1);
  global_mu_pt4 = fs->make<TH1F>("global_mu_pt4", "", 200, -1, 1);
  best_mu_pt4 = fs->make<TH1F>("best_mu_pt4", "", 200, -1, 1);
  tuneP_mu_pt4 = fs->make<TH1F>("tuneP_mu_pt4", "", 200, -1, 1);
  global_mu_pt5 = fs->make<TH1F>("global_mu_pt5", "", 200, -1, 1);
  best_mu_pt5 = fs->make<TH1F>("best_mu_pt5", "", 200, -1, 1);
  tuneP_mu_pt5 = fs->make<TH1F>("tuneP_mu_pt5", "", 200, -1, 1);
  global_mu_pt6 = fs->make<TH1F>("global_mu_pt6", "", 200, -1, 1);
  best_mu_pt6 = fs->make<TH1F>("best_mu_pt6", "", 200, -1, 1);
  tuneP_mu_pt6 = fs->make<TH1F>("tuneP_mu_pt6", "", 200, -1, 1);
  
  iso03 = fs->make<TH1F>("iso03", "", 100, 0, 10);
  iso05 = fs->make<TH1F>("iso05", "", 100, 0, 10);
  pf_iso03 = fs->make<TH1F>("pf_iso03", "", 100, 0, 10);
  pf_iso04 = fs->make<TH1F>("pf_iso04", "", 100, 0, 10);

  Mjjll = fs->make<TH1F>("Mjjll", "", 100, 0, 10000);
  Mqqll = fs->make<TH1F>("Mqqll", "", 100, 0, 10000);
  Mjjlglg = fs->make<TH1F>("Mjjlglg", "", 100, 0, 10000);
  Mqqlglg = fs->make<TH1F>("Mqqlglg", "", 100, 0, 10000);

  typeDiff_all = fs->make<TH2F>("typeDiff_all", "", 6, 0, 6, 6, 0, 6);
  typeDiff = fs->make<TH2F>("typeDiff", "", 6, 0, 6, 6, 0, 6);
  typeDiff_pt = fs->make<TH2F>("typeDiff_pt", "", 200, 0, 1500, 200, 0, 1500);
    
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

  std::vector<pat::Electron> eles;
  std::vector<pat::Muon> mus;
  std::vector<TLorentzVector> best_mus;
  std::vector<reco::GenParticle> gen_mus;
  std::vector<reco::GenParticle> gjets;
  std::vector<reco::GenParticle> gmus_W;  
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

    if((mu.pt() > 40) && (fabs(mu.eta()) < 2.4) && mu.isHighPtMuon(primary_vertex->at(0)) && MC_match){
      mus.push_back(mu);      
      TLorentzVector bmu;
      if(mu.tunePMuonBestTrack().isAvailable())
	bmu.SetPtEtaPhiM(mu.tunePMuonBestTrack()->pt(),mu.tunePMuonBestTrack()->eta(),mu.tunePMuonBestTrack()->phi(),0.1);
      best_mus.push_back(bmu);
    }
  }

  for(const pat::Jet& jet : *jets){  
    if((jet.pt() > jet_pt_min) && (fabs(jet.eta()) < jet_eta_max) && //(jet.getPFConstituents().size() > 1) & 
       (jet.neutralEmEnergyFraction() < 0.99) && (jet.chargedEmEnergyFraction() < 0.99) && 
       (jet.neutralHadronEnergyFraction() < 0.99) && (jet.chargedHadronEnergyFraction() > 0.0) && 
       (jet.muonEnergyFraction() < 0.8 ) && (jet.chargedMultiplicity() > 0)){   
      //jet_pt->Fill(jet.pt());
      pjets.push_back(jet);    
    }
  }

  int n = 0;
  int m = 0;

  for(auto gj:*gen_particles){   
    if(gj.pt() > 40 && gj.mother() != 0){

      if(abs(gj.mother()->pdgId()) == 9900014 && abs(gj.pdgId()) > 9900000)
	n++;

      if(abs(gj.pdgId()) == 13 && gj.status() == 1)
	m++;

      if(abs(gj.mother()->pdgId()) == 9900014 && abs(gj.pdgId()) < 10){
	gjets.push_back(gj);      
      }

      if((abs(gj.mother()->pdgId()) == 9900014 || abs(gj.mother()->pdgId()) == 9900024) && abs(gj.pdgId()) == 13)
	gmus_W.push_back(gj);
      
    }
  }

  if(n > 10){
    cout<<n<<" "<<m<<endl;
    cout<<"Event"<<endl;
    for(auto gj:*gen_particles){
      if(gj.mother() != 0)
	cout<<gj.pdgId()<<" "<<gj.mother()->pdgId()<<" "<<gj.status()<<" "<<gj.pt()<<endl;
    }
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

    gen_mus.push_back(gmus.back());

    typeDiff_all->Fill(mu.muonBestTrackType(),mu.tunePMuonBestTrackType());

    if(mu.muonBestTrackType() != mu.tunePMuonBestTrackType()){
      //std::cout<<"Diff TYPE="<<mu.muonBestTrackType()<<" "<<mu.tunePMuonBestTrackType()<<std::endl;      
      //std::cout<<"TUNEP TYPE="<<mu.tunePMuonBestTrackType()<<std::endl;
      typeDiff->Fill(mu.muonBestTrackType(),mu.tunePMuonBestTrackType());
      if(mu.bestTrack() != 0 && mu.tunePMuonBestTrack().isAvailable())
	typeDiff_pt->Fill(mu.bestTrack()->pt(),mu.tunePMuonBestTrack()->pt());
    }

    if(mu.muonBestTrackType() == 2){
      std::cout<<"StandAlone Muon"<<std::endl;
      if(mu.bestTrack() != 0){
	std::cout<<"Event="<<event.id().event()<<std::endl;
	std::cout<<"Run="<<event.id().run()<<std::endl;
	std::cout<<"Lumi="<<event.id().luminosityBlock()<<std::endl;
	std::cout<<"pT="<< mu.bestTrack()->pt() <<std::endl;
	std::cout<<"eta="<< mu.bestTrack()->eta() <<std::endl;
	std::cout<<"phi="<< mu.bestTrack()->phi() <<std::endl;
      }
      else std::cout<<"No track available"<<std::endl;
    }
      
    iso03->Fill((mu.isolationR03().sumPt+mu.isolationR03().emEt+mu.isolationR03().hadEt)/mu.pt());
    iso05->Fill((mu.isolationR05().sumPt+mu.isolationR05().emEt+mu.isolationR05().hadEt)/mu.pt());
    pf_iso03->Fill((mu.pfIsolationR03().sumChargedHadronPt+max(0.,mu.pfIsolationR03().sumNeutralHadronEt+mu.pfIsolationR03().sumPhotonEt - 0.5*mu.pfIsolationR03().sumPUPt))/mu.pt());
    pf_iso04->Fill((mu.pfIsolationR04().sumChargedHadronPt+max(0.,mu.pfIsolationR04().sumNeutralHadronEt+mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt))/mu.pt());

    if(mu.globalTrack().isAvailable()){
      if(mu.globalTrack()->pt() < 100)
	global_mu_pt1->Fill((mu.globalTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.globalTrack()->pt() < 200 && mu.globalTrack()->pt() > 100)
	global_mu_pt2->Fill((mu.globalTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.globalTrack()->pt() < 300 && mu.globalTrack()->pt() > 200)
	global_mu_pt3->Fill((mu.globalTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.globalTrack()->pt() < 400 && mu.globalTrack()->pt() > 300)
	global_mu_pt4->Fill((mu.globalTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.globalTrack()->pt() < 500 && mu.globalTrack()->pt() > 400)
	global_mu_pt5->Fill((mu.globalTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.globalTrack()->pt() > 500)
	global_mu_pt6->Fill((mu.globalTrack()->pt() - gmus.back().pt())/gmus.back().pt());
    }
    if(mu.bestTrack() != 0){
      if(mu.bestTrack()->pt() < 100)
	best_mu_pt1->Fill((mu.bestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.bestTrack()->pt() < 200 && mu.bestTrack()->pt() > 100)
	best_mu_pt2->Fill((mu.bestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.bestTrack()->pt() < 300 && mu.bestTrack()->pt() > 200)
	best_mu_pt3->Fill((mu.bestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.bestTrack()->pt() < 400 && mu.bestTrack()->pt() > 300)
	best_mu_pt4->Fill((mu.bestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.bestTrack()->pt() < 500 && mu.bestTrack()->pt() > 400)
	best_mu_pt5->Fill((mu.bestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.bestTrack()->pt() > 500)
	best_mu_pt6->Fill((mu.bestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
    }

    //picky_mu_pt->Fill(mu.pickyTrack()->pt() - gmus.back().pt());
    //PF_mu_pt->Fill(mu.pfCandidateRef()->pt() - gmus.back().pt());
    //tpfms_mu_pt->Fill(mu.tpfmsTrack()->pt() - gmus.back().pt());

    if(mu.tunePMuonBestTrack().isAvailable()){
      if(mu.tunePMuonBestTrack()->pt() < 100)
	tuneP_mu_pt1->Fill((mu.tunePMuonBestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.tunePMuonBestTrack()->pt() < 200 && mu.tunePMuonBestTrack()->pt() > 100)
	tuneP_mu_pt2->Fill((mu.tunePMuonBestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.tunePMuonBestTrack()->pt() < 300 && mu.tunePMuonBestTrack()->pt() > 200)
	tuneP_mu_pt3->Fill((mu.tunePMuonBestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.tunePMuonBestTrack()->pt() < 400 && mu.tunePMuonBestTrack()->pt() > 300)
	tuneP_mu_pt4->Fill((mu.tunePMuonBestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.tunePMuonBestTrack()->pt() < 500 && mu.tunePMuonBestTrack()->pt() > 400)
	tuneP_mu_pt5->Fill((mu.tunePMuonBestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
      if(mu.tunePMuonBestTrack()->pt() > 500)
	tuneP_mu_pt6->Fill((mu.tunePMuonBestTrack()->pt() - gmus.back().pt())/gmus.back().pt());
    }
  }
  
  if(best_mus.size() > 1 && pjets.size() > 1)
    Mjjll->Fill((mus[0].p4()+mus[1].p4()+pjets[0].p4()+pjets[1].p4()).M());
  if(best_mus.size() > 1 && gjets.size() > 1)
    Mqqll->Fill((mus[0].p4()+mus[1].p4()+gjets[0].p4()+gjets[1].p4()).M());
  if(gmus_W.size() > 1 && pjets.size() > 1)
    Mjjlglg->Fill((gmus_W[0].p4()+gmus_W[1].p4()+pjets[0].p4()+pjets[1].p4()).M());
  if(gmus_W.size() > 1 && gjets.size() > 1)
    Mqqlglg->Fill((gmus_W[0].p4()+gmus_W[1].p4()+gjets[0].p4()+gjets[1].p4()).M());

}

DEFINE_FWK_MODULE(Resolution);
