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

class GenHistos : public edm::EDAnalyzer {
public:
  explicit GenHistos(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

private:
  const edm::InputTag gen_src;
  const edm::InputTag gen_jet_src;
  const double jet_pt_min;
  const double jet_eta_max;
  const double dilepton_mass_cut;
  const double lljj_mass_cut;

  TH1F* ParticleID;
  TH1F* M_ee;
  TH1F* M_emu;
  TH1F* M_mumu;
  TH1F* M_jj;
  TH1F* M_eeqq;
  TH1F* M_eejj;
  TH1F* M_emujj;
  TH1F* M_mumujj;
  TH1F* M_N1_ee;
  TH1F* M_N2_ee;
  TH1F* M_N1_mumu;
  TH1F* M_N2_mumu;
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

  TH1F* deta_ee;
  TH1F* dphi_ee;
  TH1F* deta_mumu;
  TH1F* dphi_mumu;
  TH1F* deta_jets;
  TH1F* dphi_jets;

  TH1F* deta_jje1;
  TH1F* deta_jje2;
  TH1F* dphi_jje1;
  TH1F* dphi_jje2;
  TH1F* deta_jjmu1;
  TH1F* deta_jjmu2;
  TH1F* dphi_jjmu1;
  TH1F* dphi_jjmu2;

  TH1F* dR_e1j1;
  TH1F* dR_e1j2;
  TH1F* dR_e2j1;
  TH1F* dR_e2j2;
  TH1F* dR_mu1j1;
  TH1F* dR_mu1j2;
  TH1F* dR_mu2j1;
  TH1F* dR_mu2j2;

  TH1F* N_pv;
  TH1F* N_jets;

  TH1F* PV_z;
};

GenHistos::GenHistos(const edm::ParameterSet& cfg)
  : gen_src(cfg.getParameter<edm::InputTag>("gen_src")),
    gen_jet_src(cfg.getParameter<edm::InputTag>("gen_jet_src")),
    jet_pt_min(cfg.getParameter<double>("jet_pt_min")),  
    jet_eta_max(cfg.getParameter<double>("jet_eta_max")),
    dilepton_mass_cut(cfg.getParameter<double>("dilepton_mass_cut")),
    lljj_mass_cut(cfg.getParameter<double>("lljj_mass_cut"))
{
  edm::Service<TFileService> fs;

  ParticleID = fs->make<TH1F>("ParticleID", "", 40, -9910000, 9910000);
  M_ee = fs->make<TH1F>("M_ee", "", 100, 0, 3000);
  M_emu = fs->make<TH1F>("M_emu", "", 100, 0, 3000);
  M_mumu = fs->make<TH1F>("M_mumu", "", 100, 0, 3000);
  M_jj = fs->make<TH1F>("M_jj", "", 100, 0, 3000);
  M_eeqq = fs->make<TH1F>("M_eeqq", "", 100, 0, 6000);
  M_eejj = fs->make<TH1F>("M_eejj", "", 100, 0, 6000);
  M_emujj = fs->make<TH1F>("M_emujj", "", 100, 0, 6000);
  M_mumujj = fs->make<TH1F>("M_mumujj", "", 100, 0, 6000);
  M_N1_ee = fs->make<TH1F>("M_N1_ee", "", 100, 0, 5000);
  M_N2_ee = fs->make<TH1F>("M_N2_ee", "", 100, 0, 5000);
  M_N1_mumu = fs->make<TH1F>("M_N1_mumu", "", 100, 0, 5000);
  M_N2_mumu = fs->make<TH1F>("M_N2_mumu", "", 100, 0, 5000);
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

  deta_ee = fs->make<TH1F>("deta_ee", "", 100, -5, 5);
  dphi_ee = fs->make<TH1F>("dphi_ee", "", 100, -3.15, 3.15);
  deta_mumu = fs->make<TH1F>("deta_mumu", "", 100, -5, 5);
  dphi_mumu = fs->make<TH1F>("dphi_mumu", "", 100, -3.15, 3.15);
  deta_jets = fs->make<TH1F>("deta_jets", "", 100, -5, 5);
  dphi_jets = fs->make<TH1F>("dphi_jets", "", 100, -3.15, 3.15);

  deta_jje1 = fs->make<TH1F>("deta_jje1", "", 100, -5, 5);
  dphi_jje1 = fs->make<TH1F>("dphi_jje1", "", 100, -3.15, 3.15);
  deta_jje2 = fs->make<TH1F>("deta_jje2", "", 100, -5, 5);
  dphi_jje2 = fs->make<TH1F>("dphi_jje2", "", 100, -3.15, 3.15);
  deta_jjmu1 = fs->make<TH1F>("deta_jjmu1", "", 100, -5, 5);
  dphi_jjmu1 = fs->make<TH1F>("dphi_jjmu1", "", 100, -3.15, 3.15);
  deta_jjmu2 = fs->make<TH1F>("deta_jjmu2", "", 100, -5, 5);
  dphi_jjmu2 = fs->make<TH1F>("dphi_jjmu2", "", 100, -3.15, 3.15);

  dR_e1j1 = fs->make<TH1F>("dR_e1j1", "", 100, 0, 1);
  dR_e1j2 = fs->make<TH1F>("dR_e1j2", "", 100, 0, 1);
  dR_e2j1 = fs->make<TH1F>("dR_e2j1", "", 100, 0, 1);
  dR_e2j2 = fs->make<TH1F>("dR_e2j2", "", 100, 0, 1);

  dR_mu1j1 = fs->make<TH1F>("dR_mu1j1", "", 100, 0, 1);
  dR_mu1j2 = fs->make<TH1F>("dR_mu1j2", "", 100, 0, 1);
  dR_mu2j1 = fs->make<TH1F>("dR_mu2j1", "", 100, 0, 1);
  dR_mu2j2 = fs->make<TH1F>("dR_mu2j2", "", 100, 0, 1);

  N_pv = fs->make<TH1F>("N_pv", "", 100, 0, 30);
  N_jets = fs->make<TH1F>("N_jets", "", 40, 0, 40);
  //NumLeptons->SetTitle(";number of leptons from top decays;events");
  PV_z = fs->make<TH1F>("PV_z", "", 100, 0, 30);
}


void GenHistos::analyze(const edm::Event& event, const edm::EventSetup& setup) {

  edm::Handle<reco::GenParticleCollection> gen_particles;
  event.getByLabel(gen_src, gen_particles);
  //std::cout<<"EVENT"<<std::endl;
  std::vector<reco::GenParticle> eles;
  std::vector<reco::GenParticle> mus;
  std::vector<reco::GenParticle> WR_daughters;
  std::vector<reco::GenJet> jets;  

  for(const reco::GenParticle& gen : *gen_particles){    
    if(gen.mother() != 0){
      if((fabs(gen.mother()->pdgId()) == 9900024) || (fabs(gen.mother()->pdgId()) == 9900012) ){
	WR_daughters.push_back(gen);
	//std::cout<<gen.pdgId()<<" "<<gen.status()<<std::endl;
      }
    }
    if(fabs(gen.status()) != 1 )       
	 continue;
    //std::cout<<"EVENT1="<<gen.pdgId()<<std::endl;
    ParticleID->Fill(gen.pdgId());
    //std::cout<<"EVENT2"<<std::endl;
    if(fabs(gen.pdgId()) == 11) {
      //std::cout<<"status="<<gen.status()<<std::endl;
      if(gen.mother() != 0){
	if((fabs(gen.mother()->pdgId()) == 11) || (fabs(gen.mother()->pdgId()) == 24) || (fabs(gen.mother()->pdgId()) > 9900000))
	  eles.push_back(gen);
      }
    }
    if(fabs(gen.pdgId()) == 13) {
      //std::cout<<"status="<<gen.status()<<std::endl;
      if(gen.mother() != 0){
	if((fabs(gen.mother()->pdgId()) == 13) || (fabs(gen.mother()->pdgId()) == 24) || (fabs(gen.mother()->pdgId()) > 9900000))
	  mus.push_back(gen);
      }
    }


  }

  PV_z->Fill(gen_particles->at(2).vz());

  edm::Handle<reco::GenJetCollection> gen_jets;
  event.getByLabel(gen_jet_src, gen_jets);
  //std::cout<<"Number jets="<<gen_jets->size()<<std::endl;

  jets.push_back(gen_jets->at(0));
  jets.push_back(gen_jets->at(1));

  jet1_pt->Fill(gen_jets->at(0).pt());
  jet2_pt->Fill(gen_jets->at(1).pt());
  jet1_eta->Fill(gen_jets->at(0).eta());
  jet2_eta->Fill(gen_jets->at(1).eta());
  jet1_phi->Fill(gen_jets->at(0).phi());
  jet2_phi->Fill(gen_jets->at(1).phi());

  deta_jets->Fill(gen_jets->at(0).eta()-gen_jets->at(1).eta());
  dphi_jets->Fill(gen_jets->at(0).phi()-gen_jets->at(1).phi());

  for(const reco::GenJet& jet : *gen_jets){
    if((jet.pt() > 40) & (fabs(jet.eta()) < 2.5) )
      jet_pt->Fill(jet.pt());
  }

  double Mee = 0;
  double Meejj = 0;

  if(eles.size() > 1){
    Mee = (eles[0].p4()+eles[1].p4()).M();
    Meejj = (jets[0].p4()+jets[1].p4()+eles[0].p4()+eles[1].p4()).M();
  }
  double Mmumu = 0;
  double Mmumujj = 0;

  if(mus.size() > 1){
    Mmumu = (mus[0].p4()+mus[1].p4()).M();
    Mmumujj = (jets[0].p4()+jets[1].p4()+mus[0].p4()+mus[1].p4()).M();
  }
  //std::cout<<"EVENT2"<<std::endl;
  if(Mee > dilepton_mass_cut && Meejj > lljj_mass_cut){
    e1_pt->Fill(eles[0].pt());
    e1_eta->Fill(eles[0].eta());
    e1_phi->Fill(eles[0].phi());
    //std::cout<<"EVENT1"<<std::endl;
    e1_dz->Fill(fabs(eles[0].vz() - gen_particles->at(2).vz()));
    e2_pt->Fill(eles[1].pt());
    e2_eta->Fill(eles[1].eta());
    e2_phi->Fill(eles[1].phi());
    e2_dz->Fill(fabs(eles[1].vz() - gen_particles->at(2).vz()));

    deta_ee->Fill(eles[0].eta()-eles[1].eta());
    dphi_ee->Fill(deltaPhi(eles[0].phi(),eles[1].phi()));
    M_N1_ee->Fill((jets[0].p4()+jets[1].p4()+eles[0].p4()).M());
    M_N2_ee->Fill((jets[0].p4()+jets[1].p4()+eles[1].p4()).M());    
    M_ee->Fill(Mee);
    M_eejj->Fill(Meejj);    
    deta_jje1->Fill((jets[0].p4()+jets[1].p4()+eles[0].p4()).eta());
    dphi_jje1->Fill((jets[0].p4()+jets[1].p4()+eles[0].p4()).phi());
    deta_jje2->Fill((jets[0].p4()+jets[1].p4()+eles[1].p4()).eta());
    dphi_jje2->Fill((jets[0].p4()+jets[1].p4()+eles[1].p4()).phi());
    dR_e1j1->Fill(deltaR(eles[0],jets[0]));    
    dR_e1j2->Fill(deltaR(eles[0],jets[1]));    
    dR_e2j1->Fill(deltaR(eles[1],jets[0]));    
    dR_e2j2->Fill(deltaR(eles[1],jets[1])); 
  }

  jet1_pt->Fill(jets[0].pt());
  jet1_eta->Fill(jets[0].eta());
  jet1_phi->Fill(jets[0].phi());  
  jet2_pt->Fill(jets[1].pt());
  jet2_eta->Fill(jets[1].eta());
  jet2_phi->Fill(jets[1].phi());  
  M_jj->Fill((jets[0].p4()+jets[1].p4()).M());
  deta_jets->Fill(jets[0].eta()-jets[1].eta());
  dphi_jets->Fill(deltaPhi(jets[0].phi(),jets[1].phi()));

  if(WR_daughters.size()>3)
    M_eeqq->Fill((WR_daughters[0].p4()+WR_daughters[1].p4()+WR_daughters[2].p4()+WR_daughters[3].p4()).M());  
  
  if(Mmumu > dilepton_mass_cut && Mmumujj > lljj_mass_cut){
    mu1_pt->Fill(mus[0].pt());
    mu1_eta->Fill(mus[0].eta());
    mu1_phi->Fill(mus[0].phi());
    //std::cout<<"EVENT1"<<std::endl;
    mu2_pt->Fill(mus[1].pt());
    mu2_eta->Fill(mus[1].eta());
    mu2_phi->Fill(mus[1].phi());
    M_mumu->Fill(Mmumu);
    M_mumujj->Fill(Mmumujj);
    deta_mumu->Fill(mus[0].eta()-mus[1].eta());
    dphi_mumu->Fill(deltaPhi(mus[0].phi(),mus[1].phi()));
    M_N1_mumu->Fill((jets[0].p4()+jets[1].p4()+mus[0].p4()).M());
    M_N2_mumu->Fill((jets[0].p4()+jets[1].p4()+mus[1].p4()).M());
    deta_jjmu1->Fill((jets[0].p4()+jets[1].p4()+mus[0].p4()).eta());
    dphi_jjmu1->Fill((jets[0].p4()+jets[1].p4()+mus[0].p4()).phi());
    deta_jjmu2->Fill((jets[0].p4()+jets[1].p4()+mus[1].p4()).eta());
    dphi_jjmu2->Fill((jets[0].p4()+jets[1].p4()+mus[1].p4()).phi());
    dR_mu1j1->Fill(deltaR(mus[0],jets[0]));    
    dR_mu1j2->Fill(deltaR(mus[0],jets[1]));    
    dR_mu2j1->Fill(deltaR(mus[1],jets[0]));    
    dR_mu2j2->Fill(deltaR(mus[1],jets[1]));    
  }

  N_pv->Fill(1);
  N_jets->Fill(jets.size());


}

DEFINE_FWK_MODULE(GenHistos);
