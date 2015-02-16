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

class WR_Analyzer : public edm::EDAnalyzer {
public:
  explicit WR_Analyzer(const edm::ParameterSet&);
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

  TH1F* ParticleID;
  TH1F* electronID;
  TH1F* electronID2;
  TH1F* M_ee;
  TH1F* M_emu;
  TH1F* M_mumu;
  TH1F* M_jj;
  TH1F* M_eeqq;
  TH1F* M_eejj;
  TH1F* M_emujj;
  TH1F* M_mumujj;
  TH1F* M_N1;
  TH1F* M_N2;
  TH1F* mMET;
  TH1F* jet_pt;
  TH1F* jet1_pt;
  TH1F* jet2_pt;
  TH1F* jet1_eta;
  TH1F* jet2_eta;
  TH1F* jet1_phi;
  TH1F* jet2_phi;
  TH1F* e1_pt;
  TH1F* e2_pt;
  TH1F* e1_eta;
  TH1F* e2_eta;
  TH1F* e1_phi;
  TH1F* e2_phi;
  TH1F* e1_dz;
  TH1F* e2_dz;
  TH1F* mu1_pt;
  TH1F* mu2_pt;
  TH1F* mu1_eta;
  TH1F* mu2_eta;
  TH1F* mu1_phi;
  TH1F* mu2_phi;

  TH1F* deta_leptons;
  TH1F* dphi_leptons;
  TH1F* deta_jets;
  TH1F* dphi_jets;

  TH1F* N_pv;
  TH1F* N_jets;

  // Electron Plots
  TH1F* fe_pt;
  TH1F* fe_eta;
  TH1F* fe_mass;
  TH1F* fe_pdgId;
  TH1F* fe_idAvailable;
  TH1F* ge_pt;
  TH1F* ge_eta;
  TH1F* ge_mass;
  TH1F* ge_pdgId;
  TH1F* ge_idAvailable;

  TH1F* ele_pre0;
  TH1F* ele_pre1;
  TH1F* ele_pre2;
  TH1F* ele_pre3;
  TH1F* ele_pre4;
  TH1F* ele_pre5;
  TH1F* ele_pre6;
  TH1F* mu_pre0;
  TH1F* mu_pre1;
  TH1F* mu_pre2;
  TH1F* mu_pre3;
  TH1F* mu_pre4;
  TH1F* mu_pre5;
  TH1F* mu_pre6;

};

WR_Analyzer::WR_Analyzer(const edm::ParameterSet& cfg)
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
  ele_pre0 = fs->make<TH1F>("ele_pre0", "", 2, 0, 2);
  ele_pre1 = fs->make<TH1F>("ele_pre1", "", 2, 0, 2);
  ele_pre2 = fs->make<TH1F>("ele_pre2", "", 2, 0, 2);
  ele_pre3 = fs->make<TH1F>("ele_pre3", "", 2, 0, 2);
  ele_pre4 = fs->make<TH1F>("ele_pre4", "", 2, 0, 2);
  ele_pre5 = fs->make<TH1F>("ele_pre5", "", 2, 0, 2);
  ele_pre6 = fs->make<TH1F>("ele_pre6", "", 2, 0, 2);
  mu_pre0 = fs->make<TH1F>("mu_pre0", "", 2, 0, 2);
  mu_pre1 = fs->make<TH1F>("mu_pre1", "", 2, 0, 2);
  mu_pre2 = fs->make<TH1F>("mu_pre2", "", 2, 0, 2);
  mu_pre3 = fs->make<TH1F>("mu_pre3", "", 2, 0, 2);
  mu_pre4 = fs->make<TH1F>("mu_pre4", "", 2, 0, 2);
  mu_pre5 = fs->make<TH1F>("mu_pre5", "", 2, 0, 2);
  mu_pre6 = fs->make<TH1F>("mu_pre6", "", 2, 0, 2);
  ParticleID = fs->make<TH1F>("ParticleID", "", 40, -9910000, 9910000);
  electronID = fs->make<TH1F>("electronID", "", 10, 0, 10);
  electronID2 = fs->make<TH1F>("electronID2", "", 10, 0, 10);
  M_ee = fs->make<TH1F>("M_ee", "", 100, 0, 3000);
  M_emu = fs->make<TH1F>("M_emu", "", 100, 0, 3000);
  M_mumu = fs->make<TH1F>("M_mumu", "", 100, 0, 3000);
  M_jj = fs->make<TH1F>("M_jj", "", 100, 0, 3000);
  M_eeqq = fs->make<TH1F>("M_eeqq", "", 100, 0, 6000);
  M_eejj = fs->make<TH1F>("M_eejj", "", 100, 0, 6000);
  M_emujj = fs->make<TH1F>("M_emujj", "", 100, 0, 6000);
  M_mumujj = fs->make<TH1F>("M_mumujj", "", 100, 0, 6000);
  M_N1 = fs->make<TH1F>("M_N1", "", 100, 0, 5000);
  M_N2 = fs->make<TH1F>("M_N2", "", 100, 0, 5000);
  mMET = fs->make<TH1F>("mMET", "", 100, 0, 1000);
  jet_pt = fs->make<TH1F>("jet_pt", "", 100, 0, 2000);
  jet1_pt = fs->make<TH1F>("jet1_pt", "", 100, 0, 2000);
  jet2_pt = fs->make<TH1F>("jet2_pt", "", 100, 0, 2000);
  jet1_eta = fs->make<TH1F>("jet1_eta", "", 100, -3, 3);
  jet2_eta = fs->make<TH1F>("jet2_eta", "", 100, -3, 3);
  jet1_phi = fs->make<TH1F>("jet1_phi", "", 100, -3.15, 3.15);
  jet2_phi = fs->make<TH1F>("jet2_phi", "", 100, -3.15, 3.15);
  e1_pt = fs->make<TH1F>("e1_pt", "", 100, 0, 1000);
  e2_pt = fs->make<TH1F>("e2_pt", "", 100, 0, 1000);
  e1_eta = fs->make<TH1F>("e1_eta", "", 100, -3, 3);
  e2_eta = fs->make<TH1F>("e2_eta", "", 100, -3, 3);
  e1_phi = fs->make<TH1F>("e1_phi", "", 100, -3.15, 3.15);
  e2_phi = fs->make<TH1F>("e2_phi", "", 100, -3.15, 3.15);
  e1_dz = fs->make<TH1F>("e1_dz", "", 100, 0, 30);
  e2_dz = fs->make<TH1F>("e2_dz", "", 100, 0, 30);
  mu1_pt = fs->make<TH1F>("mu1_pt", "", 100, 0, 2000);
  mu2_pt = fs->make<TH1F>("mu2_pt", "", 100, 0, 2000);
  mu1_eta = fs->make<TH1F>("mu1_eta", "", 100, -3, 3);
  mu2_eta = fs->make<TH1F>("mu2_eta", "", 100, -3, 3);
  mu1_phi = fs->make<TH1F>("mu1_phi", "", 100, -3.15, 3.15);
  mu2_phi = fs->make<TH1F>("mu2_phi", "", 100, -3.15, 3.15);

  deta_leptons = fs->make<TH1F>("deta_leptons", "", 100, -5, 5);
  dphi_leptons = fs->make<TH1F>("dphi_leptons", "", 100, -3.15, 3.15);
  deta_jets = fs->make<TH1F>("deta_jets", "", 100, -5, 5);
  dphi_jets = fs->make<TH1F>("dphi_jets", "", 100, -3.15, 3.15);

  N_pv = fs->make<TH1F>("N_pv", "", 100, 0, 30);
  N_jets = fs->make<TH1F>("N_jets", "", 40, 0, 40);

  // Electron Plots
  fe_pt = fs->make<TH1F>("fe_pt", "", 100, 0, 1000);
  fe_eta = fs->make<TH1F>("fe_eta", "", 100, -3, 3);
  fe_mass = fs->make<TH1F>("fe_mass", "", 100, -1, 1);
  fe_pdgId = fs->make<TH1F>("fe_pdgId", "", 100, -50, 50);
  fe_idAvailable = fs->make<TH1F>("fe_idAvailable", "", 2, 0, 2);
  ge_pt = fs->make<TH1F>("ge_pt", "", 100, 0, 1000);
  ge_eta = fs->make<TH1F>("ge_eta", "", 100, -3, 3);
  ge_mass = fs->make<TH1F>("ge_mass", "", 100, -1, 1);
  ge_pdgId = fs->make<TH1F>("ge_pdgId", "", 100, -50, 50);
  ge_idAvailable = fs->make<TH1F>("ge_idAvailable", "", 2, 0, 2);

  
}


void WR_Analyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  using namespace std;
  using namespace edm;
  using namespace reco;  

  edm::Handle<pat::ElectronCollection> electrons;
  event.getByLabel(electron_src, electrons);
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
  event.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  event.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  event.getByToken(electronMediumIdMapToken_,medium_id_decisions);
  event.getByToken(electronTightIdMapToken_,tight_id_decisions);
  event.getByToken(electronHEEPIdMapToken_,heep_id_decisions);

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
  std::vector<reco::GenParticle> WR_daughters;
  std::vector<reco::GenJet> gjets;  
  std::vector<pat::Jet> pjets;

  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel("offlineSlimmedPrimaryVertices", vertices);

  for(const reco::GenParticle& gen : *gen_particles){
    ParticleID->Fill(gen.pdgId());
  }

  for(std::vector<pat::Electron>::const_iterator el = electrons->begin();
       el != electrons->end(); el++){

    // Look up the ID decision for this electron in 
    // the ValueMap object and store it. We need a Ptr object as the key.
    const Ptr<pat::Electron> elPtr(electrons, el - electrons->begin() );   

    ele_pre0->Fill(1);
    if(el->pt() > 40){
      bool isPassVeto = false;
      bool isPassHEEP = false;
      ele_pre1->Fill(1);
      if(fabs(el->superCluster()->eta()) < 1.442 ||  (fabs(el->superCluster()->eta()) > 1.56 && fabs(el->superCluster()->eta()) < 2.5)){
	ele_pre2->Fill(1);
	if(el->passConversionVeto()){
	  ele_pre3->Fill(1);
	  e1_dz->Fill(fabs(vertices->at(0).z() - el->vz()));
	  if(fabs(vertices->at(0).z() - el->vz()) < 1){
	    ele_pre4->Fill(1);	  
	    if((*heep_id_decisions)[ elPtr ])
	      ele_pre5->Fill(1);
	    if((*veto_id_decisions)[ elPtr ])
	      ele_pre6->Fill(1);

	    isPassVeto  = (*veto_id_decisions)[ elPtr ];
	    //isPassLoose = (*loose_id_decisions)[ elPtr ];
	    //isPassMedium = (*medium_id_decisions)[ elPtr ];
	    //isPassTight = (*tight_id_decisions)[ elPtr ];
	    isPassHEEP = (*heep_id_decisions)[ elPtr ];
	  }
	}
      }
      
      if(isPassHEEP){
	eles.push_back(*el);
	ge_pt ->Fill(el->pt());
	ge_eta ->Fill(el->eta());
	ge_mass ->Fill(el->mass());
	ge_pdgId ->Fill(el->pdgId());
	//ge_idAvailable ->Fill(el->isElectronIDAvailable("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto"));
      }
      else {
	fe_pt ->Fill(el->pt());
	fe_eta ->Fill(el->eta());
	fe_mass ->Fill(el->mass());
	fe_pdgId ->Fill(el->pdgId());
	//fe_idAvailable ->Fill(el->isElectronIDAvailable("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto"));
      }
      electronID->Fill(isPassHEEP);
      electronID2->Fill(isPassVeto);
    }
  }

  for(const pat::Muon& mu : *muons){
    mu_pre0->Fill(1);
    if(mu.pt() > 40){
      mu_pre1->Fill(1);
      if(fabs(mu.eta()) < 2.4){
	mu_pre2->Fill(1);
	bool MC_match = false;
	float dR = 0.1;
	for(const reco::GenParticle& gen : *gen_particles){
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
    if((mu.pt() > 40) && (fabs(mu.eta()) < 2.4) && mu.isHighPtMuon(primary_vertex->at(0)))
      mus.push_back(mu);    
  }

  if(gen_jets->size() > 1){
    gjets.push_back(gen_jets->at(0));
    gjets.push_back(gen_jets->at(1));
  }
  
  for(const pat::Jet& jet : *jets){  
    if((jet.pt() > jet_pt_min) && (fabs(jet.eta()) < jet_eta_max) && //(jet.getPFConstituents().size() > 1) & 
       (jet.neutralEmEnergyFraction() < 0.99) && (jet.chargedEmEnergyFraction() < 0.99) && 
       (jet.neutralHadronEnergyFraction() < 0.99) && (jet.chargedHadronEnergyFraction() > 0.0) && 
       (jet.muonEnergyFraction() < 0.8 ) && (jet.chargedMultiplicity() > 0)){   
      jet_pt->Fill(jet.pt());
      pjets.push_back(jet);    
    }
  }

  pat::MET met = mets->at(0);
  double Mee = 0.0;
  double Meejj = 0.0;

  if(eles.size() > 1){
    Mee = (eles[0].p4()+eles[1].p4()).M();
    if(pjets.size() > 1)
      Meejj = (pjets[0].p4()+pjets[1].p4()+eles[0].p4()+eles[1].p4()).M();
  }
  
  if(Mee > dilepton_mass_cut && Meejj > lljj_mass_cut){
    e1_pt->Fill(eles[0].pt());
    e1_eta->Fill(eles[0].eta());
    e1_phi->Fill(eles[0].phi());
    e2_pt->Fill(eles[1].pt());
    e2_eta->Fill(eles[1].eta());
    e2_phi->Fill(eles[1].phi());
    deta_leptons->Fill(eles[0].eta()-eles[1].eta());
    dphi_leptons->Fill(deltaPhi(eles[0].phi(),eles[1].phi()));
    if(pjets.size() > 1){
      M_N1->Fill((pjets[0].p4()+pjets[1].p4()+eles[0].p4()).M());
      M_N2->Fill((pjets[0].p4()+pjets[1].p4()+eles[1].p4()).M());    
      M_ee->Fill(Mee);
      M_eejj->Fill(Meejj);      
    }    
  }

  mMET->Fill(met.pt());
  N_pv->Fill(1);
  N_jets->Fill(jets->size());
  if(pjets.size() > 1){
    jet1_pt->Fill(pjets[0].pt());
    jet1_eta->Fill(pjets[0].eta());
    jet1_phi->Fill(pjets[0].phi());  
    jet2_pt->Fill(pjets[1].pt());
    jet2_eta->Fill(pjets[1].eta());
    jet2_phi->Fill(pjets[1].phi());  
    deta_jets->Fill(pjets[0].eta()-pjets[1].eta());
    dphi_jets->Fill(deltaPhi(pjets[0].phi(),pjets[1].phi()));
    M_jj->Fill((pjets[0].p4()+pjets[1].p4()).M());
  }

  double Mmumu = 0.0;
  double Mmumujj = 0.0;
  if(mus.size() > 1){
    Mmumu = (mus[0].p4()+mus[1].p4()).M();  
    if(pjets.size() > 1) 
      Mmumujj = (pjets[0].p4()+pjets[1].p4()+mus[0].p4()+mus[1].p4()).M();
  }
  if(Mmumu > dilepton_mass_cut && Mmumujj > lljj_mass_cut){
    mu1_pt->Fill(mus[0].p4().pt());
    mu1_eta->Fill(mus[0].p4().eta());
    mu1_phi->Fill(mus[0].p4().phi());
    mu2_pt->Fill(mus[1].p4().pt());
    mu2_eta->Fill(mus[1].p4().eta());
    mu2_phi->Fill(mus[1].p4().phi());
    M_mumu->Fill(Mmumu);
    if(pjets.size() > 1){
      M_mumujj->Fill(Mmumujj);
    }
  }

  if((mus.size() > 0) && (eles.size() > 0) ){
    M_emu->Fill((eles[0].p4()+mus[0].p4()).M());
    if(pjets.size() > 1) 
      M_emujj->Fill((pjets[0].p4()+pjets[1].p4()+eles[0].p4()+mus[0].p4()).M());    
  }
  if(WR_daughters.size()>3)
    M_eeqq->Fill((WR_daughters[0].p4()+WR_daughters[1].p4()+WR_daughters[2].p4()+WR_daughters[3].p4()).M());  

}

DEFINE_FWK_MODULE(WR_Analyzer);
