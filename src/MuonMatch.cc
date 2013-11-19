// -*- C++ -*-
//
// Package:    MuonMatch
// Class:      MuonMatch
// 
/**\class MuonMatch MuonMatch.cc L1Analysis/MuonMatch/src/MuonMatch.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Florian Scheuch
//         Created:  Thu Oct 17 10:45:56 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/PhysObjectMatcher.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "TMath.h"
#include "TH1D.h"

//Include Muon Propagator
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
//
// class declaration
//

class MuonMatch : public edm::EDAnalyzer {
   public:
      explicit MuonMatch(const edm::ParameterSet&);
      ~MuonMatch();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual bool equalsGenParticles(reco::GenParticle gp1, reco::GenParticle gp2);
      virtual void defineNumberOfMuonGhosts(edm::Handle<edm::View<reco::Muon>> muonHandle, edm::Handle<reco::MuonCollection> muonCollectionHandle, edm::Handle<reco::GenParticleCollection> genParticleHandle, edm::Handle<reco::GenParticleMatch> genParticleMatchHandle, double lower_energy_threshold, TH1D* histogram, TH1D* histogram_barrel, TH1D* histogram_endcap, const edm::Event& iEvent);
      
    //  PropagateToMuon propagatorToMuon;
      // ----------member data ---------------------------
	//TH1D* counter;
	TH1D* counter_collection; //No of RecoMuons/GenMuons
	TH1D* phi_ghosts; //phi direction of ghosts
	TH1D* eta_ghosts; //eta direction of ghosts

	TH1D* phi_gen; //phi direction of generated mu
	TH1D* eta_gen; //eta direction of generated mu

	TH1D* gen_ghost_energy; // Energy of the ghost producing particle
	TH1D* muonHandleSize;
	TH1D* eventNoCatcher;
	TH1D* number_of_ghosts[10];
	TH1D* number_of_ghosts_barrel[10];
	TH1D* number_of_ghosts_endcap[10];
	TH1D* gen_mu_energy;
	//Event counter
	int event_count;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonMatch::MuonMatch(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


MuonMatch::~MuonMatch()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void MuonMatch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	event_count += 1;
	if(event_count%100000 == 0) std::cout << "Analyzing event number " << event_count << std::endl; 
/*	edm::Handle<std::vector<reco::>> ak5GenJets; //genParticleHandle
	iEvent.getByLabel("ak5GenJets", ak5GenJets);
	
	edm::Handle<edm::View<reco::PFJet> > ak5PFJets; //muonHandle
	iEvent.getByLabel("ak5PFJetsJEC", ak5PFJets);
	
	edm::Handle<edm::Association<reco::GenJetCollection>> genJetMatch; //genParticleMatch
	iEvent.getByLabel("JetGenJetMatch", genJetMatch);
*/	
	//MuonCollection
	edm::Handle<edm::View<reco::Muon>> muonHandle;
	iEvent.getByLabel("muons",muonHandle);
	
	edm::Handle<reco::MuonCollection> muonCollectionHandle;
	iEvent.getByLabel("muons", muonCollectionHandle);
	//get gen particle collection
	edm::Handle<reco::GenParticleCollection> genParticleHandle;
	iEvent.getByLabel("genParticles",genParticleHandle);

	//Get association with genParticle matches
	edm::Handle<reco::GenParticleMatch> genParticleMatchHandle;
	iEvent.getByLabel("muonMatch",genParticleMatchHandle);
	//EndMuonCollection
	
	//edm::Handle<reco::CandMatchMap> match; //new
	//iEvent.getByLabel("muonMatch", match); //new
	
	
	for(unsigned int energy = 0; energy < 10; energy ++){
		defineNumberOfMuonGhosts(muonHandle, muonCollectionHandle, genParticleHandle, genParticleMatchHandle, energy*5, number_of_ghosts[energy], number_of_ghosts_barrel[energy], number_of_ghosts_endcap[energy], iEvent);
	}
	


	

		
		
		//if (mcMatch.isNonnull() && mcMatch.isAvailable()) {
		//	reco::GenParticleCollection genparticlecollection = *(mcMatch.get());
		//	counter_collection->Fill(genparticlecollection->size());
		//Analyse Hier
		//}
	
	

	int genMuonCounter = 0;

	//std::cout << genParticleHandle->size() << " Groesse" << std::endl;
	for(unsigned int idx = 0; idx < genParticleHandle->size(); idx++){
	//	std::cout << "Durchlauf " << idx << std::endl;
		if (genParticleHandle->at(idx).pdgId() == 13 || genParticleHandle->at(idx).pdgId() == -13) {
			genMuonCounter = genMuonCounter+1;
			gen_mu_energy->Fill(genParticleHandle->at(idx).energy());
			phi_gen->Fill(genParticleHandle->at(idx).phi());
			eta_gen->Fill(genParticleHandle->at(idx).eta());
		}
	}
	//std::cout << "GenMuonCounter: " << genMuonCounter << std::endl;
	counter_collection->Fill(1.*muonHandle->size()/genMuonCounter);
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonMatch::beginJob() {
	edm::Service<TFileService> fs;
	
	TFileDirectory dir_number_of_ghosts = fs->mkdir("Number_of_ghosts_sum");
	TFileDirectory dir_number_of_ghosts_barrel = fs->mkdir("Number_of_ghosts_barrel");
	TFileDirectory dir_number_of_ghosts_endcap = fs->mkdir("Number_of_ghosts_endcap");

	const double PI  = 3.141592653589793238462;

	//counter = fs->make<TH1D>("counter","counter", 10, 0, 9);
	counter_collection = fs->make<TH1D>("counter_collection", "counter_collection", 100, 0, 2);
	eta_ghosts = fs->make<TH1D>("Eta distribution of ghosts", "Eta distribution of ghosts", 2000, 0, 2.5);
	phi_ghosts = fs->make<TH1D>("Phi distribution of ghosts", "Phi distribution of ghosts", 192, -1*PI, PI);
	gen_ghost_energy = fs->make<TH1D>("Energy of the ghost produciong particle", "Energy of the ghost producing particle", 5000, 0, 4999);
	gen_mu_energy = fs->make<TH1D>("Energy of the generated mu", "Energy of the generated mu", 5000, 0, 4999);
	eta_gen = fs->make<TH1D>("Eta distribution of generated mu", "Eta distribution of generated mu", 2000, 0, 2.5);
	phi_gen = fs->make<TH1D>("Phi distribution of generated mu", "Phi distribution of generated mu", 192, -1*PI, PI);
	
	number_of_ghosts[0] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 0 GeV","Number of ghosts, pt_min = 0 GeV", 4, -0.5, 3.5);
	number_of_ghosts[1] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 5 GeV","Number of ghosts, pt_min = 5 GeV", 4, -0.5, 3.5);
	number_of_ghosts[2] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 10 GeV","Number of ghosts, pt_min = 10 GeV", 4, -0.5, 3.5);
	number_of_ghosts[3] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 15 GeV","Number of ghosts, pt_min = 15 GeV", 4, -0.5, 3.5);
	number_of_ghosts[4] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 20 GeV","Number of ghosts, pt_min = 20 GeV", 4, -0.5, 3.5);
	number_of_ghosts[5] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 25 GeV","Number of ghosts, pt_min = 25 GeV", 4, -0.5, 3.5);
	number_of_ghosts[6] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 30 GeV","Number of ghosts, pt_min = 30 GeV", 4, -0.5, 3.5);
	number_of_ghosts[7] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 35 GeV","Number of ghosts, pt_min = 35 GeV", 4, -0.5, 3.5);
	number_of_ghosts[8] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 40 GeV","Number of ghosts, pt_min = 40 GeV", 4, -0.5, 3.5);
	number_of_ghosts[9] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 45 GeV","Number of ghosts, pt_min = 45 GeV", 4, -0.5, 3.5);

	number_of_ghosts_barrel[0] = dir_number_of_ghosts_barrel.make<TH1D>("Number of ghosts in barrel region, pt_min = 0 GeV","Number of ghosts (barrel), pt_min = 0 GeV", 4, -0.5, 3.5);
	number_of_ghosts_barrel[1] = dir_number_of_ghosts_barrel.make<TH1D>("Number of ghosts in barrel region, pt_min = 5 GeV","Number of ghosts (barrel), pt_min = 5 GeV", 4, -0.5, 3.5);
	number_of_ghosts_barrel[2] = dir_number_of_ghosts_barrel.make<TH1D>("Number of ghosts in barrel region, pt_min = 10 GeV","Number of ghosts (barrel), pt_min = 10 GeV", 4, -0.5, 3.5);
	number_of_ghosts_barrel[3] = dir_number_of_ghosts_barrel.make<TH1D>("Number of ghosts in barrel region, pt_min = 15 GeV","Number of ghosts (barrel), pt_min = 15 GeV", 4, -0.5, 3.5);
	number_of_ghosts_barrel[4] = dir_number_of_ghosts_barrel.make<TH1D>("Number of ghosts in barrel region, pt_min = 20 GeV","Number of ghosts (barrel), pt_min = 20 GeV", 4, -0.5, 3.5);
	number_of_ghosts_barrel[5] = dir_number_of_ghosts_barrel.make<TH1D>("Number of ghosts in barrel region, pt_min = 25 GeV","Number of ghosts (barrel), pt_min = 25 GeV", 4, -0.5, 3.5);
	number_of_ghosts_barrel[6] = dir_number_of_ghosts_barrel.make<TH1D>("Number of ghosts in barrel region, pt_min = 30 GeV","Number of ghosts (barrel), pt_min = 30 GeV", 4, -0.5, 3.5);
	number_of_ghosts_barrel[7] = dir_number_of_ghosts_barrel.make<TH1D>("Number of ghosts in barrel region, pt_min = 35 GeV","Number of ghosts (barrel), pt_min = 35 GeV", 4, -0.5, 3.5);
	number_of_ghosts_barrel[8] = dir_number_of_ghosts_barrel.make<TH1D>("Number of ghosts in barrel region, pt_min = 40 GeV","Number of ghosts (barrel), pt_min = 40 GeV", 4, -0.5, 3.5);
	number_of_ghosts_barrel[9] = dir_number_of_ghosts_barrel.make<TH1D>("Number of ghosts in barrel region, pt_min = 45 GeV","Number of ghosts (barrel), pt_min = 45 GeV", 4, -0.5, 3.5);

	number_of_ghosts_endcap[0] = dir_number_of_ghosts_endcap.make<TH1D>("Number of ghosts in endcap region, pt_min = 0 GeV","Number of ghosts (endcap), pt_min = 0 GeV", 4, -0.5, 3.5);
	number_of_ghosts_endcap[1] = dir_number_of_ghosts_endcap.make<TH1D>("Number of ghosts in endcap region, pt_min = 5 GeV","Number of ghosts (endcap), pt_min = 5 GeV", 4, -0.5, 3.5);	
	number_of_ghosts_endcap[2] = dir_number_of_ghosts_endcap.make<TH1D>("Number of ghosts in endcap region, pt_min = 10 GeV","Number of ghosts (endcap), pt_min = 10 GeV", 4, -0.5, 3.5);
	number_of_ghosts_endcap[3] = dir_number_of_ghosts_endcap.make<TH1D>("Number of ghosts in endcap region, pt_min = 15 GeV","Number of ghosts (endcap), pt_min = 15 GeV", 4, -0.5, 3.5);
	number_of_ghosts_endcap[4] = dir_number_of_ghosts_endcap.make<TH1D>("Number of ghosts in endcap region, pt_min = 20 GeV","Number of ghosts (endcap), pt_min = 20 GeV", 4, -0.5, 3.5);
	number_of_ghosts_endcap[5] = dir_number_of_ghosts_endcap.make<TH1D>("Number of ghosts in endcap region, pt_min = 25 GeV","Number of ghosts (endcap), pt_min = 25 GeV", 4, -0.5, 3.5);
	number_of_ghosts_endcap[6] = dir_number_of_ghosts_endcap.make<TH1D>("Number of ghosts in endcap region, pt_min = 30 GeV","Number of ghosts (endcap), pt_min = 30 GeV", 4, -0.5, 3.5);
	number_of_ghosts_endcap[7] = dir_number_of_ghosts_endcap.make<TH1D>("Number of ghosts in endcap region, pt_min = 35 GeV","Number of ghosts (endcap), pt_min = 35 GeV", 4, -0.5, 3.5);
	number_of_ghosts_endcap[8] = dir_number_of_ghosts_endcap.make<TH1D>("Number of ghosts in endcap region, pt_min = 40 GeV","Number of ghosts (endcap), pt_min = 40 GeV", 4, -0.5, 3.5);
	number_of_ghosts_endcap[9] = dir_number_of_ghosts_endcap.make<TH1D>("Number of ghosts in endcap region, pt_min = 45 GeV","Number of ghosts (endcap), pt_min = 45 GeV", 4, -0.5, 3.5);
	
	event_count = 0;
}

//--- determines the number of muon ghosts in all events with a certain energy cut
void
MuonMatch::defineNumberOfMuonGhosts(edm::Handle<edm::View<reco::Muon>> muonHandle, edm::Handle<reco::MuonCollection> muonCollectionHandle, edm::Handle<reco::GenParticleCollection> genParticleHandle, edm::Handle<reco::GenParticleMatch> genParticleMatchHandle, double lower_energy_threshold, TH1D* histogram, TH1D* histogram_barrel, TH1D* histogram_endcap, const edm::Event& iEvent){
	std::vector<reco::GenParticle> processed_particles;
	int count_barrel = 0; // counts the number of ghosts in the barrel region
	int count_endcap = 0; // counts the number of ghosts in the endcap region
	for (unsigned int idx = 0; idx < muonHandle->size(); idx++) { // loop over all reco muons and find identical corresponding genParticles
		
		edm::RefToBase<reco::Muon> muonRef = muonHandle->refAt(idx);
		reco::Muon amuon(*(muonRef.get())); //RecoMuon in Variable amuon
		
		reco::GenParticleRef genparticleref = (*genParticleMatchHandle)[muonRef];
		reco::MuonRef ref_muon = muonRef.castTo<reco::MuonRef>();
		if (genparticleref.isNonnull() && genparticleref.isAvailable() && (amuon.isStandAloneMuon() || amuon.isGlobalMuon())) {
			const reco::GenParticle genparticle = *(genparticleref.get());

			if(genparticle.pdgId() != amuon.pdgId()){
				edm::LogInfo ("PdgId not fitting") << "Gen ID: " << genparticle.pdgId() << " Reco ID: " << amuon.pdgId() << " With energy > " << lower_energy_threshold << " GeV in LumiBlock: " << iEvent.id().luminosityBlock() << " Event No: " << iEvent.id().event();
			}
			
			count_barrel = 0;
			count_endcap = 0;
			for(unsigned int jdx = 0; jdx < processed_particles.size(); jdx++){
				reco::GenParticle gp = processed_particles.at(jdx);
				if(equalsGenParticles(gp, genparticle)){
					if(amuon.pt() > lower_energy_threshold){
						if(amuon.eta()<1.24 && amuon.eta()>-1.24) count_barrel++;
						else count_endcap++;
					}
					if(lower_energy_threshold==0 && amuon.pt() > lower_energy_threshold){ //nur bei dem ersten durchgang auffüllen
						phi_ghosts->Fill(amuon.phi());
						eta_ghosts->Fill(amuon.eta());
						gen_ghost_energy->Fill(genparticle.energy());
					}
		//			std::cout << "Particles are equal: " << count << std::endl;
				}
			}
			processed_particles.push_back(genparticle);
		}
	}
	if(count_barrel > 0){
					//std::cout << "Found " << count << " Ghosts with energy > " << lower_energy_threshold << " GeV in LumiBlock: " << iEvent.id().luminosityBlock() << " Event No: " << iEvent.id().event() << std::endl;
				edm::LogInfo ("Ghost_found_barrel") << "Found " << count_barrel << " Ghosts in barrel region with energy > " << lower_energy_threshold << " GeV in LumiBlock: " << iEvent.id().luminosityBlock() << " Event No: " << iEvent.id().event();
	}
	if(count_endcap > 0){
					//std::cout << "Found " << count << " Ghosts with energy > " << lower_energy_threshold << " GeV in LumiBlock: " << iEvent.id().luminosityBlock() << " Event No: " << iEvent.id().event() << std::endl;
				edm::LogInfo ("Ghost_found_endcap") << "Found " << count_endcap << " Ghosts in endcap region with energy > " << lower_energy_threshold << " GeV in LumiBlock: " << iEvent.id().luminosityBlock() << " Event No: " << iEvent.id().event();
	}
	histogram->Fill(count_barrel+count_endcap);
	histogram_barrel->Fill(count_barrel);
	histogram_endcap->Fill(count_endcap);
}	










			//edm::OwnVector<reco::Candidate> ownVector(); //new
			//ownVector.pushBack(&amuon); //new
			//edm::Ref<edm::OwnVector<reco::Candidate>> refOwnVector(ownVector, muonRef.key()); //new
			
			//reco::CandMatchMap::const_iterator i = match->find(refOwnVector); //new
			//counter_collection->Fill(match->numberOfAssociations(muonRef.castTo())); //new
				//reco::CandidateRef mcMatch = (*match)[muonRef]; //new




// ------------ method called once each job just after ending the event loop  ------------
void 
MuonMatch::endJob() 
{
}

//Checks if the GenParticles are equal in condition of charge, p_t
bool
MuonMatch::equalsGenParticles(reco::GenParticle gp1, reco::GenParticle gp2){
	if(gp1.pt() != gp2.pt()){
		//counter->Fill(gp1.pt()/gp2.pt());
	//	std::cout << "Pt is different" << std::endl;
		return false;
	}
	if(gp1.charge() != gp2.charge()){
	//	std::cout << "Charge is different" << std::endl;
		return false;
	}
	return true;
}

// ------------ method called when starting to processes a run  ------------
void 
MuonMatch::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
	//propagatorToMuon.init(iSetup);
}

// ------------ method called when ending the processing of a run  ------------
void 
MuonMatch::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MuonMatch::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MuonMatch::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonMatch::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonMatch);


/*
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/003697A5-03FA-E211-A6C3-0025907DCAA6.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0042F21A-8CF9-E211-BCA5-0025907FD428.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/00794212-CCF9-E211-B7EE-0030487E55C5.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/026499D5-CEF9-E211-B046-0030487E5101.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0445D1E4-2CF9-E211-93B4-0030487EBB2B.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/046A8F09-EEF8-E211-87E5-003048C69408.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/06579C7B-20F9-E211-8D77-0030487E4EB9.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/06CF10CB-34F9-E211-98FF-0030487E5101.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0887DC15-96F9-E211-9B58-0025907DC9BA.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0A0CE369-14F9-E211-AC48-0030487D83B9.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0A107581-07F9-E211-8EF0-0030487D5D4D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0C06E963-CCF9-E211-B8FB-0030487D5EAF.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0C84489F-6DF9-E211-BBF6-002590A2CCF2.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0E2E9642-3AF9-E211-83D8-003048D43656.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/0E5356CD-92F9-E211-872C-00266CF2AACC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/1009B78F-07FA-E211-A817-00266CFFA048.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/1010EE5A-4FF9-E211-B3BB-0025901D493A.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/12D57B94-9DF9-E211-8404-00266CF32E2C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/12D9980C-BEF9-E211-B057-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/14B1A4D5-0FF9-E211-9660-0030487D70FD.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/1694007F-93F9-E211-8A6B-003048C66184.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/16F68992-3AF9-E211-B192-0030487D5EBD.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/1805F48C-88F9-E211-A022-0025907FD442.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/1C87CE7D-32F9-E211-9241-003048C692DE.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/1C9F3005-CCF9-E211-9E79-0030487E0A2D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/1CCD75C6-79F9-E211-9E38-0025907DCA64.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/1ED74858-45F9-E211-A0BC-0030487DE7C1.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/2067C8F5-31F9-E211-AF14-0030487D5E5F.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/207FDC5A-CDF9-E211-952F-0030487D5E9D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/20F323AA-54F9-E211-94FC-002481E945C4.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/2276C9EC-A0F9-E211-9166-0025904B1452.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/22CDD520-8CF9-E211-88D1-002590A2CDA8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/24533375-57F9-E211-ABA6-003048D43986.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/2463AEE2-4FF9-E211-BB25-0025901D493A.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/263D2C4B-8DF9-E211-92F1-003048D439AC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/264B3ED0-8DF9-E211-A2A2-0025907FD2D8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/28AA0B9B-9EF9-E211-9D0C-00266CF332B4.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/28D533BD-D3F9-E211-9678-0030487D5EA9.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/2C26FB99-BCF9-E211-A6A0-0030487D5E9D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/2CCFCCB3-BEF9-E211-BC93-0030487D5E9D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/2CF1ADB4-92F9-E211-842F-003048D43700.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/2E06EE53-BFF9-E211-9894-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/2E3BFA95-87F9-E211-84FC-002590A2CD44.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/301356A9-02F9-E211-9E70-0030487E4EC5.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/30552712-CCF9-E211-B52D-0030487E55BB.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/305AC52D-3BF9-E211-AF93-003048D43656.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/30FD03E7-9FF9-E211-B0F3-0025901D493A.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/321C07D5-6AF9-E211-B7C4-0025907FD242.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/321CCA6F-80F9-E211-87BE-00266CFFA7BC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/32251453-9AF9-E211-B6B9-0025907DC9F8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/32B72F7D-C2F9-E211-AEF9-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/3409B563-99F9-E211-9D77-00266CFFA6F8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/34675A33-99F9-E211-A1F4-0025907FD2BA.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/363B0839-10F9-E211-BDE2-0030487D70FD.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/384CDB6A-E3F9-E211-96AD-0030487D5059.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/3A4F464D-BEF9-E211-B4DA-0030487D5E9D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/3C4E603C-41F9-E211-84A2-0030487D5EB5.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/3C90CF63-9AF9-E211-950B-0025901D4844.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/3CB3D79F-05FA-E211-A6DF-0025907DC9C6.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/3CC206B8-00F9-E211-9483-0030487D5D89.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/3E00E132-8EF9-E211-A095-003048F0E82C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/400A8052-00F9-E211-B0AC-0030487D5D89.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/409D41EB-8AF9-E211-A0D4-0025907FD428.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/42AF1AFE-8AF9-E211-B73F-002590A2CD44.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4420CB95-34F9-E211-9097-0030487E52A3.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/44C87EBB-12F9-E211-B9A2-0030487D70FF.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/44EA31EF-45F9-E211-AE68-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/46164723-74F9-E211-B0DE-00266CF9C1AC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/465D376D-43F9-E211-A555-003048D43730.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4662085E-D1F9-E211-A66E-0030487E4EB9.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/46CEC32C-9FF9-E211-A8C7-00266CF332B4.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/46F2F104-20F9-E211-ABD5-0030487E4EB9.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4825DFBA-6FF9-E211-8B7A-0025907DC9C6.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4850B90F-93F9-E211-8D6B-003048C66184.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/48675DB0-96F9-E211-9EB6-00266CF9C1AC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4876E14E-97F9-E211-BCF2-00266CFFA6A8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/48AFE15B-06F9-E211-AC75-0030487D8633.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/48D58CE4-8DF9-E211-96E7-0025907FD2D6.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4A1A4701-98F9-E211-A280-003048D4399C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4A2D28C8-5AF9-E211-A52E-002590A2CD68.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4AC04E72-4CF9-E211-A84E-003048C66184.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4ACDE8D4-D3F9-E211-B1DF-0030487D5059.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4C021C5A-94F9-E211-858B-0025907DCA0C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4E7E7F3A-91F9-E211-8BD3-003048C68A80.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/4EF4B0E8-9FF9-E211-8F83-0025904B1452.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/5072BA30-05FA-E211-838B-00266CFFA148.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/5231C3F5-99F9-E211-A3C5-0025907DC9F8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/52A54094-4BF9-E211-B28A-0030487D5E53.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/52C2D6D7-1CFA-E211-8817-0030487D5059.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/54724834-6BF9-E211-A70E-0025907FD242.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/5472990E-9FF9-E211-A0B0-0025904B1452.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/54D9C2BD-D2F9-E211-A60B-0030487D5EA9.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/588EACA0-F8F8-E211-A431-0030487D5E49.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/58F7C7BA-45F9-E211-AA38-0030487DE7C1.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/5A66EA0F-9EF9-E211-A5CD-00266CF332B4.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/5C87EB17-90F9-E211-9DB2-0025907FD280.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/602FE384-7AF9-E211-80C3-0025901D42C0.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/60328A76-06FA-E211-B3FC-00266CFFA7BC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/607B264B-89F9-E211-9F34-002590A2CD06.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/609BB7F2-BCF9-E211-BC38-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/622476DC-DEF9-E211-88D6-0030487E4EC5.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/6401353C-91F9-E211-A305-0030487D5D95.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/644019F3-43F9-E211-A540-003048D439C6.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/662D8F59-2DF9-E211-A049-0030487EBB2B.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/66888134-9BF9-E211-B222-002481E10CC2.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/669559BE-9AF9-E211-BCE2-002481E10CC2.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/6AC5A9B3-94F9-E211-B5EA-0025907DC9D2.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/6C995596-53F9-E211-A895-00266CFFA6F8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/6E0EF95A-A1F9-E211-8C9C-0025904B1452.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/6E5490F0-E2F9-E211-AFD2-00266CFFA604.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/6E890C54-8DF9-E211-B253-003048D47912.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/6EA55A15-99F9-E211-9883-0025901D4B4C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/70BAC087-9FF9-E211-89BE-0025904B1452.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/7203DC71-CCF9-E211-B436-0030487D5059.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/724D6099-BDF9-E211-BD96-0030487D5E9D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/725E62DD-9DF9-E211-A5D3-003048C66BBE.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/72981B11-CCF9-E211-AB90-0030487D5D89.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/747DD363-BCF9-E211-A583-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/74F5D506-95F9-E211-ABC7-003048C69412.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/765CEA06-9FF9-E211-91E2-00266CF32E2C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/7811B559-CFF9-E211-82DF-0030487D7B79.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/7C24CA09-BCF9-E211-B6B4-0030487D5E9D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/7C988E2E-0AF9-E211-B940-0030487E5179.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/7E05AC50-88F9-E211-892C-0025907FD2B2.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/7E60E516-92F9-E211-9791-0030487D70FF.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/7E8DD0A0-95F9-E211-99CE-0025907DC9BA.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/7EAE2F78-0BFA-E211-80E9-0030487DE7C5.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/7EE6B395-41F9-E211-9D18-003048C662D4.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/80D315DC-3BF9-E211-A383-003048D43656.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/842C7027-70F9-E211-AF11-0025907DC9C6.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/847F4CBB-93F9-E211-98D3-0025907DCA9C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/84809E30-35F9-E211-A868-003048C69408.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/84B2378C-18F9-E211-937E-0030487E54B5.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/865084EF-90F9-E211-8C57-003048D439A6.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/866A0013-DCF9-E211-9505-0030487D5E53.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/8679FA73-CCF9-E211-AAF0-0030487D7109.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/88230625-F5F9-E211-A28C-0025907FD430.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/886F33D1-48F9-E211-A88F-003048D4DFA4.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/889959A7-04F9-E211-9BB3-0030487D5D89.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/88AA117C-93F9-E211-AB44-0025907FD4A0.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/88F167A1-8FF9-E211-811B-003048C692E2.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/8C6E099D-0DF9-E211-B217-003048C693C8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/8CFA28FA-13F9-E211-9C74-0030487D5EBB.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/8EB1E1D4-D5F9-E211-8F35-0030487D5D95.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/902A0299-8FF9-E211-914D-003048D4DFBA.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/90B4D98A-82F9-E211-A8A5-0030487D5D91.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/90E0A535-8FF9-E211-848B-0025907DCA72.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/92863F66-81F9-E211-A69A-00266CFFA7BC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/92A15391-9FF9-E211-A25E-00266CF32E2C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/92A9C91D-81F9-E211-9BE5-002590A2CDA8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/9406B30C-77F9-E211-B1FE-002590494CC2.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/965F6941-8BF9-E211-A79C-0025907FD2BA.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/96762E1A-31F9-E211-B1B9-0030487E52A3.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/9A209C7C-99F9-E211-8555-0025907DC9F8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/9A74C541-3EF9-E211-AB16-0030487E52A3.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/9AC8905F-BFF9-E211-9B95-0030487D5E9D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/9E84BC0E-02F9-E211-9615-003048D439AC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/A00E038F-1CFA-E211-9ACB-0030487D8121.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/A022D4B5-8EF9-E211-AF1B-0025901D4844.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/A0C10A11-CCF9-E211-B000-0030487D7109.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/A2492CE2-D6F9-E211-A679-0030487D605B.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/A2FED4E5-F8F8-E211-9954-0030487E52A3.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/A686DD55-94F9-E211-8868-003048C69412.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/A826DA66-99F9-E211-96C6-00266CFFA1AC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/A8CC23D3-A1F9-E211-9967-0025904B1452.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/A8E54EA4-FBF8-E211-B67F-003048C6903A.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/AA397DE5-80F9-E211-BD5B-00266CFFA7BC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/AA7C4132-9BF9-E211-AEC5-0025907DC9F8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/AAC79F9D-BDF9-E211-BECC-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/AC16AE72-BBF9-E211-AB1C-0030487D5E9D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/AC4EAEB9-09FA-E211-BCE7-0025901D4C18.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/AC8330A5-9EF9-E211-85D6-0025904B1452.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/AC965214-96F9-E211-83DF-003048C69412.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/AE3EEDB8-BBF9-E211-B20F-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/B0371249-93F9-E211-BE7D-00266CF2AACC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/B0418E99-9EF9-E211-84EC-00266CF32E2C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/B0B8936B-8CF9-E211-89BF-0025907DC9F2.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/B0E6C022-A2F9-E211-A768-0025904B1026.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/B0F3C7DB-07F9-E211-88C1-0030487E4B99.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/B2005E01-38F9-E211-882A-003048D462DC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/B496283B-FBF8-E211-9158-0030487E4EC5.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/B6A37517-CCF9-E211-B87E-0030487E52A1.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/B8E713D1-BFF9-E211-A7AF-0030487D5E9D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/BE29643A-37F9-E211-B1BD-0030487D5EAF.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/BEAA21A8-8AF9-E211-BE6B-0025907FD3CC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/BEDAAEC5-96F9-E211-BB0F-00266CF33118.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/C2186044-BBF9-E211-AB8D-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/C23DD4F1-3DF9-E211-8B37-003048D37570.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/C2FCAB31-9EF9-E211-8567-00266CF32E2C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/C420CD2A-BDF9-E211-971D-0030487D5E9D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/C42FB9C2-84F9-E211-A616-002481E94B26.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/C46ED3DF-08FA-E211-A5EF-002590494E34.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/C63274EB-91F9-E211-8F37-0030487D5D95.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/C6AB2B44-D6F9-E211-90FB-0030487D5D95.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/CA3D6E53-EEF9-E211-9EDB-0030487D5E45.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/CA4F6F29-97F9-E211-96B9-00266CF9C1AC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/CA938222-00FA-E211-A721-00266CFFA278.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/CC0E767C-86F9-E211-AA0D-0025907FD424.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/CCB92A11-0BF9-E211-A851-0030487D5D8F.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/CCDEDBA5-D1F9-E211-B288-0030487D5D95.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/CE646CB4-73F9-E211-89DC-00266CF9C1AC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D03796A1-D2F9-E211-831C-0030487D8121.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D08E5926-7AF9-E211-8716-0025907DCA64.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D0F3938B-F3F8-E211-9146-0030487D5DBD.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D23DA220-13F9-E211-9209-0030487D70FF.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D29AEC63-0AF9-E211-90B5-0030487D5E73.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D2B1DB5B-8EF9-E211-BF89-0025907FD2D8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D434B37F-95F9-E211-87B4-003048C69412.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D47BE41D-07F9-E211-8576-0030487D5D4D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D601F2C6-95F9-E211-88ED-003048C68A8E.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D6067BD4-9AF9-E211-BE2E-0025907DC9F8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D66E8390-91F9-E211-B669-003048D439A6.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D6B3C140-9FF9-E211-868F-0025901D493A.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D6DD6004-63F9-E211-85AF-00266CF332B4.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D83F3EC5-33F9-E211-8061-0030487E55C5.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D84149D7-CDF9-E211-A0BA-0030487E52A3.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/D84E8F44-95F9-E211-A6F1-0025907DC9D2.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/DAD0F92F-84F9-E211-AF51-002481E94B26.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/DC2AABE3-DBF9-E211-B25A-0030487D83B9.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/DCAA20D6-99F9-E211-8550-00266CFFA6F8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/DE0446A3-9BF9-E211-B231-002590A2CD06.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/DE1052E5-BEF9-E211-B751-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/DEACB55E-F6F8-E211-9EFD-0030487E52A3.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/E01AD8C3-99F9-E211-8ACA-0025901D4B4C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/E06FB584-0BF9-E211-8CC4-0030487D5EAF.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/E4EE388C-98F9-E211-9D51-0025901D4B4C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/E64AA55F-96F9-E211-B707-00266CF33118.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/E8B6A50B-C2F9-E211-9367-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/EA536CD2-86F9-E211-8961-0025907FD424.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/EE1A8589-A0F9-E211-AA95-0025904B1452.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/EEA409FC-97F9-E211-A627-0025901D4B4C.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/EECB2C87-8BF9-E211-A290-0025907FD428.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/F00E4ADC-FCF8-E211-A716-0030487D8581.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/F0538580-BEF9-E211-B291-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/F2FDBDD0-8EF9-E211-B144-0030487E4ED5.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/F4BA47B2-05F9-E211-996A-0030487D8563.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/F64EEA3E-90F9-E211-B3FC-003048D439A6.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/F6E07EA6-3EF9-E211-B1A6-0030487E52A3.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/F811A35F-97F9-E211-8C1D-00266CF33118.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/FA040145-3CF9-E211-986D-003048D43656.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/FA9B8E07-05F9-E211-9CFD-0030487D7109.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/FC8D72EF-62F9-E211-BFC0-002590494E36.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/FCF8B31A-7CF9-E211-B41A-003048D436EA.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/FEC7B91D-DEF9-E211-957A-0030487E0A29.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/10000/FEEAD48B-5FF9-E211-B10A-0025907FD2BA.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/06E3B99B-BAF9-E211-97AB-0030487EBB25.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/16E2E61E-BAF9-E211-8ECA-0030487D43E1.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/26718C34-BAF9-E211-9063-0030487D5DAF.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/28348A80-BAF9-E211-99D8-0030487D5E49.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/307B2A88-96F9-E211-AA48-0025904B12A8.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/343C9658-C2F9-E211-9787-0030487D5059.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/38F0169B-BAF9-E211-8496-0030487E4B8D.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/447D955D-C2F9-E211-9383-0030487D8581.root'])
fileList.extend(['/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/4CA3E040-7FF9-E211-89CC-003048F0E5A2.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/56E4D6AF-BAF9-E211-BD7A-0030487D5D89.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/5CE03BA8-BAF9-E211-95DB-0030487D5E49.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/62BFC496-95F9-E211-95B9-003048F0E3BC.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/6ED00299-BAF9-E211-A056-0030487E0A29.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/78CE402B-BAF9-E211-95BB-0030487EBB2B.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/9407E293-BAF9-E211-B4F6-0030487D43E1.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/9E76461B-BAF9-E211-A891-0030487D5D89.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/BE13E393-BAF9-E211-9492-0030487EBB2B.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/CA267885-BAF9-E211-9BAC-0030487D811F.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/CC136FA9-BAF9-E211-9BAA-0030487D5DAF.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/E8A34392-BAF9-E211-97F0-0030487D5E53.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/F29DB49F-BAF9-E211-82AC-0030487E55BB.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/F46F22C2-BAF9-E211-8A01-0030487EBB21.root',
	'/store/mc/Summer12_DR53X/ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S10_START53_V19-v1/20000/F6564719-BAF9-E211-8CAB-0030487D811F.root'])
*/