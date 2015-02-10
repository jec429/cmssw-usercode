#include <memory>
#include <vector>

#include "TLorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "Math/VectorUtil.h"

class EventSelector : public edm::EDFilter {
public:
  explicit EventSelector(const edm::ParameterSet&);

private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);  
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

};

EventSelector::EventSelector(const edm::ParameterSet& cfg) 
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
}

namespace {
  template <typename T>
  T mag(T x, T y) {
    return sqrt(x*x + y*y);
  }

  template <typename T>
  T mag(T x, T y, T z) {
    return sqrt(x*x + y*y + z*z);
  }
}

bool EventSelector::filter(edm::Event& event, const edm::EventSetup&) {
  using namespace std;
  using namespace edm;
  using namespace reco;

  edm::Handle<reco::GenParticleCollection> gen_particles;
  event.getByLabel(gen_src, gen_particles);
  edm::Handle<reco::VertexCollection> primary_vertex;
  event.getByLabel(vertex_src, primary_vertex);

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

  edm::Handle<pat::JetCollection> jets;
  event.getByLabel(jet_src, jets);

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  // Jet Selection
  std::vector<pat::Jet> good_jets;
  for(const pat::Jet& jet : *jets){
    if((jet.pt() > jet_pt_min) && (fabs(jet.eta()) < jet_eta_max) && //(jet.getPFConstituents().size() > 1) & 
       (jet.neutralEmEnergyFraction() < 0.99) && (jet.chargedEmEnergyFraction() < 0.99) && 
       (jet.neutralHadronEnergyFraction() < 0.99) && (jet.chargedHadronEnergyFraction() > 0.0) && 
       (jet.muonEnergyFraction() < 0.8 ) && (jet.chargedMultiplicity() > 0)) 
      good_jets.push_back(jet);
  }
  if(good_jets.size() < 2) 
    return false;
        
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  // Lepton Selection
  
  std::vector<Int_t>   passVetoId_;     
  std::vector<Int_t>   passLooseId_; 
  std::vector<Int_t>   passMediumId_; 
  std::vector<Int_t>   passTightId_; 
  std::vector<Int_t>   passHEEPId_; 
  
  passVetoId_.clear();     
  passTightId_.clear();  
  passHEEPId_.clear();  
  
  std::vector<pat::Electron> good_electrons;
  std::vector<pat::Muon> good_muons;

  for( std::vector<pat::Electron>::const_iterator el = electrons->begin();
       el != electrons->end(); el++){
    // Look up the ID decision for this electron in 
    // the ValueMap object and store it. We need a Ptr object as the key.
    const Ptr<pat::Electron> elPtr(electrons, el - electrons->begin() );
    bool isPassVeto  = (*veto_id_decisions)[ elPtr ];
    bool isPassLoose = (*loose_id_decisions)[ elPtr ];
    bool isPassMedium = (*medium_id_decisions)[ elPtr ];
    bool isPassTight = (*tight_id_decisions)[ elPtr ];
    bool isPassHEEP = (*heep_id_decisions)[ elPtr ];
    passVetoId_.push_back( isPassVeto );
    passLooseId_.push_back( isPassLoose );
    passMediumId_.push_back( isPassMedium );
    passTightId_.push_back( isPassTight );
    passHEEPId_.push_back( isPassHEEP );
    //std::cout<<"HEEP="<<isPassHEEP<<std::endl;
    
    if(isPassHEEP)
      good_electrons.push_back(*el);
   
  }  
  
  for(const pat::Muon& mu : *muons){
    if((mu.pt() > 40) && (fabs(mu.eta()) < 2.4) && (mu.isHighPtMuon(primary_vertex->at(0))))
      good_muons.push_back(mu);    
  }

  //if((good_electrons.size() < 2) && (good_muons.size() < 2))
  //return false;

  double Mee = 0.0;
  double Meejj = 0.0;
  if(good_electrons.size() > 1){
    Mee = (good_electrons[0].p4()+good_electrons[1].p4()).M();
    if(good_jets.size() > 1)
      Meejj = (good_jets[0].p4()+good_jets[1].p4()+good_electrons[0].p4()+good_electrons[1].p4()).M();
  }
  
  if(Mee < dilepton_mass_cut || Meejj < lljj_mass_cut)
    return false;

  return true;
}

DEFINE_FWK_MODULE(EventSelector);
