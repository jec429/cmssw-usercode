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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

class Comparison : public edm::EDAnalyzer {
public:
  explicit Comparison(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

private:
  edm::EDGetTokenT<edm::ValueMap<bool> > electronHEEPIdMapToken_;

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
  TH1F* g_e1_pt;
  TH1F* g_e2_pt;
  TH1F* g_e1_eta;
  TH1F* g_e2_eta;
  TH1F* g_e1_phi;
  TH1F* g_e2_phi;
  TH1F* r_e1_pt;
  TH1F* r_e2_pt;
  TH1F* r_e1_eta;
  TH1F* r_e2_eta;
  TH1F* r_e1_phi;
  TH1F* r_e2_phi;
  TH1F* e_id;
  TH1F* d_gr_e_pt;
  TH1F* d_gr_e_eta;
  TH1F* d_gr_e_phi;
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

  TH1F* deta_jjl1;
  TH1F* deta_jjl2;
  TH1F* dphi_jjl1;
  TH1F* dphi_jjl2;

  TH1F* N_pv;
  TH1F* N_jets;

  TH1F* PV_z;
};

Comparison::Comparison(const edm::ParameterSet& cfg):
  electronHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("electronHEEPIdMap")))

{
  edm::Service<TFileService> fs;

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
  g_e1_pt = fs->make<TH1F>("g_e1_pt", "", 100, 0, 1000);
  g_e2_pt = fs->make<TH1F>("g_e2_pt", "", 100, 0, 1000);
  g_e1_eta = fs->make<TH1F>("g_e1_eta", "", 100, -3, 3);
  g_e2_eta = fs->make<TH1F>("g_e2_eta", "", 100, -3, 3);
  g_e1_phi = fs->make<TH1F>("g_e1_phi", "", 100, -3.15, 3.15);
  g_e2_phi = fs->make<TH1F>("g_e2_phi", "", 100, -3.15, 3.15);
  r_e1_pt = fs->make<TH1F>("r_e1_pt", "", 200, 0, 1000);
  r_e2_pt = fs->make<TH1F>("r_e2_pt", "", 100, 0, 1000);
  r_e1_eta = fs->make<TH1F>("r_e1_eta", "", 100, -3, 3);
  r_e2_eta = fs->make<TH1F>("r_e2_eta", "", 100, -3, 3);
  r_e1_phi = fs->make<TH1F>("r_e1_phi", "", 100, -3.15, 3.15);
  r_e2_phi = fs->make<TH1F>("r_e2_phi", "", 100, -3.15, 3.15);
  e_id = fs->make<TH1F>("e_id","",2,0,2);
  d_gr_e_pt = fs->make<TH1F>("d_gr_e_pt", "", 100, -500, 500);
  d_gr_e_eta = fs->make<TH1F>("d_gr_e_eta", "", 100, -5, 5);
  d_gr_e_phi = fs->make<TH1F>("d_gr_e_phi", "", 100, -6.3, 6.3);
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

  deta_jjl1 = fs->make<TH1F>("deta_jjl1", "", 100, -5, 5);
  dphi_jjl1 = fs->make<TH1F>("dphi_jjl1", "", 100, -3.15, 3.15);
  deta_jjl2 = fs->make<TH1F>("deta_jjl2", "", 100, -5, 5);
  dphi_jjl2 = fs->make<TH1F>("dphi_jjl2", "", 100, -3.15, 3.15);

  N_pv = fs->make<TH1F>("N_pv", "", 100, 0, 30);
  N_jets = fs->make<TH1F>("N_jets", "", 40, 0, 40);
  //NumLeptons->SetTitle(";number of leptons from top decays;events");
  PV_z = fs->make<TH1F>("PV_z", "", 100, 0, 30);
}


void Comparison::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  using namespace std;
  using namespace edm;
  using namespace reco;

  edm::Handle<reco::GenParticleCollection> gen_particles;
  event.getByLabel("prunedGenParticles", gen_particles);
  edm::Handle<reco::GenJetCollection> gen_jets;
  event.getByLabel("slimmedGenJets", gen_jets);
  edm::Handle<pat::ElectronCollection> electrons;
  event.getByLabel("slimmedElectrons", electrons);

  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
  event.getByToken(electronHEEPIdMapToken_,heep_id_decisions);


  //std::cout<<"EVENT"<<std::endl;
  std::vector<reco::GenParticle> eles;
  std::vector<reco::GenParticle> WR_daughters;
  std::vector<reco::GenJet> jets;  

  for(std::vector<pat::Electron>::const_iterator el = electrons->begin(); el != electrons->end(); el++){
    if(el->pt() > 20){// && el->pt() < 50){
      if(fabs(el->superCluster()->eta()) < 1.4442){
	if(el->passConversionVeto()){
	  //if(fabs(vertices->at(0).z() - el->vz()) < 1){
	  bool MC_match = false;
	  float dR = 0.1;
	  for(const reco::GenParticle& gen : *gen_particles){
	    if(fabs(gen.pdgId()) == 11 && deltaR(gen.eta(),gen.phi(),el->eta(),el->phi()) < dR && gen.status() == 1 && fabs(gen.mother()->pdgId()) != 15) {
	      //std::cout<<"MOM="<<gen.motherRef()->pdgId()<<std::endl;
	      MC_match = true;
	      dR = deltaR(gen.eta(),gen.phi(),el->eta(),el->phi());
	      eles.push_back(gen);
	    }		
	  }  
	  if(MC_match){
	    //std::cout<<el->ecalDrivenMomentum().x() - el->p4().x()<<std::endl;
	    const Ptr<pat::Electron> elPtr(electrons, el - electrons->begin() );
	    e_id->Fill(el->ecalDriven());
	    g_e1_pt->Fill(eles.back().pt());
	    g_e1_eta->Fill(eles.back().eta());
	    g_e1_phi->Fill(eles.back().phi());
	    if((*heep_id_decisions)[elPtr])
	      r_e1_pt->Fill(el->pt());
	    r_e1_eta->Fill(el->eta());
	    r_e1_phi->Fill(el->phi());
	    d_gr_e_pt->Fill(eles.back().pt()-el->pt());
	    d_gr_e_eta->Fill(eles.back().eta()-el->eta());
	    d_gr_e_phi->Fill(eles.back().phi()-el->phi());
	  }
	    //}
	}
      }
    }    
  }

}

DEFINE_FWK_MODULE(Comparison);
