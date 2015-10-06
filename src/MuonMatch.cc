/*
 *  MuonMatch.cc
 *  
 *
 *  Created by Florian Scheuch on 09.04.14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 *
 *	WARNING: THE COUNTING OF THE DISTANCES BETWEEN RECO MUONS AND GEN PARTICLE IS ONLY VALID IF #EVENTS WITH 1 GHOST >> #EVENTS WITH n>1 GHOSTS!
 *	
 *	Ghost definition?: Ghost is deltaR > 0.1 away from generating particle
 *	TODO: Have a look at the events with ghosts
 *
 */

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
#include "DataFormats/VertexReco/interface/Vertex.h"

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
#include "TH2D.h"

//Include Muon Propagator
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"


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
	
	TH1D* recoMuR;
	TH1D* deltaVzGenReco;
	TH1D* recoGenMuR;
	TH1D* recoGenMuRNoGhost;
	TH1D* pt_ratio_ghost_gen;
	TH1D* pt_ratio_not_ghost_gen;
	TH2D* pt_ratio_2d_gen_reco;
	TH2D* pt_ratio_2d_ghost_gen_reco;
	TH2D* pt_ratio_2d_corr_gen_reco;
	TH1D* deltaRGhostReco;
	TH1D* deltaPtGhostReco;
	TH2D* pt_ratio_deltaR_pt_2d_gen_reco;
	TH2D* pt_ratio_deltaR_pt_2d_gen_reco_all_particles;
	TH1D* calo_likelihood;
	TH1D* calo_likelihood_ghost;
	TH1D* calo_likelihood_not_ghost;
	TH1D* calo_likelihood_ghost_event;
	TH1D* ho_energy;
	TH1D* ho_energy_ghost_event;
	TH1D* ho_energy_ghost;
	TH1D* ho_energy_no_ghost;
	TH1D* ho_energy_3x3;
	TH1D* ho_energy_ghost_event_3x3;
	TH1D* ho_energy_ghost_3x3;
	TH1D* ho_energy_no_ghost_3x3;
};

MuonMatch::MuonMatch(const edm::ParameterSet& iConfig) {
	//now do what ever initialization is needed
}

MuonMatch::~MuonMatch() {
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}

void MuonMatch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	edm::Handle<edm::View<reco::Muon>> muonHandle;
	iEvent.getByLabel("muons",muonHandle);
	
	edm::Handle<reco::MuonCollection> muonCollectionHandle;
	iEvent.getByLabel("muons", muonCollectionHandle);
	
	// Get GenParticleCollection
	edm::Handle<reco::GenParticleCollection> genParticleHandle;
	iEvent.getByLabel("genParticles",genParticleHandle);
	
	// Get association with genParticle matches
	edm::Handle<reco::GenParticleMatch> genParticleMatchHandle;
	iEvent.getByLabel("muonMatch",genParticleMatchHandle);
	
	
	std::vector<reco::GenParticle> processed_particles; // Container for GenParticles that already have been associated
	std::vector<edm::RefToBase<reco::Muon>> processed_particles_ref;
	
	// Search through all reco mu
	for(unsigned int idx = 0; idx < muonHandle->size(); idx++){		
		edm::RefToBase<reco::Muon> muonRef = muonHandle->refAt(idx);
		reco::Muon amuon(*(muonRef.get()));
		
		// Check if muon in barrel region
		if(amuon.eta()>.8 && amuon.eta()<-.8) continue;
		
		// Just look for global muons
		if(!amuon.isGlobalMuon()) continue;
		
		// Get associated gen particle
		reco::GenParticleRef genparticleref = (*genParticleMatchHandle)[muonRef];
		reco::MuonRef ref_muon = muonRef.castTo<reco::MuonRef>();
		
		// Check if genParticle is present
		if (!genparticleref.isNonnull() || !genparticleref.isAvailable()){
			continue;
		}

		//CaloTowerRef caloRef = amuon.caloTower();
		//const CaloTower* caloTow = caloRef.get();
		//std::cout << caloTow->energyInHO() << " GeV in HO" << std::endl;
		ho_energy -> Fill(amuon.calEnergy().ho);
		ho_energy_3x3 -> Fill(amuon.calEnergy().hoS9);
		
		const reco::GenParticle genparticle = *(genparticleref.get());
		
		deltaVzGenReco -> Fill(genparticle.vz() - amuon.vz());
	//	std::cout << "Vz diff Gen Reco: " << fabs(genparticle.vz() - amuon.vz()) << " vz1: " << genparticle.vz() << " vz2: " << amuon.vz() << std::endl;
		
		// Check if genparticle and reco mu are from the same vertex
		if(fabs(genparticle.vz() - amuon.vz()) > 0.01) continue;
		if(fabs(genparticle.vy() - amuon.vy()) > 0.01) continue;
		if(fabs(genparticle.vx() - amuon.vx()) > 0.01) continue;
		
		double del_R_reco_gen = sqrt(
									 (genparticle.phi()-amuon.phi()) * (genparticle.phi()-amuon.phi()) +
									 (genparticle.eta()-amuon.eta()) * (genparticle.eta()-amuon.eta())
									 ); // Distance of reco mu and matching gen mu
		del_R_reco_gen = fmod(del_R_reco_gen, TMath::Pi()); //Normalize from 0 to Pi
		recoGenMuRNoGhost -> Fill(del_R_reco_gen);
		calo_likelihood -> Fill(amuon.caloCompatibility());
		

		// Go through handled genParticles, if genParticle was handled before -> this reco muon or the first muon is a muon ghost
		for(unsigned int jdx = 0; jdx < processed_particles.size(); jdx++){
			reco::GenParticle gp = processed_particles.at(jdx);
			if(equalsGenParticles(gp, genparticle)){
				reco::Muon firstmuon(*(processed_particles_ref.at(jdx).get()));

				//Fill 2D Plot gen pt reco pt
				pt_ratio_2d_gen_reco -> Fill(genparticle.pt(), firstmuon.pt());
				pt_ratio_2d_gen_reco -> Fill(genparticle.pt(), amuon.pt());

				double del_R = sqrt(
									(firstmuon.phi()-amuon.phi()) * (firstmuon.phi()-amuon.phi()) +
									(firstmuon.eta()-amuon.eta()) * (firstmuon.eta()-amuon.eta())
				); // Distance of both reco Mu
				del_R = fmod(del_R, TMath::Pi()); //Normalize from 0 to Pi
				double del_R_first_reco_gen = sqrt(
											 (genparticle.phi()-firstmuon.phi()) * (genparticle.phi()-firstmuon.phi()) +
											 (genparticle.eta()-firstmuon.eta()) * (genparticle.eta()-firstmuon.eta())
											 ); // Distance of frist reco mu and matching gen mu
				del_R_first_reco_gen = fmod(del_R_first_reco_gen, TMath::Pi()); //Normalize from 0 to Pi
				recoMuR -> Fill(del_R);
				recoGenMuR -> Fill(del_R_reco_gen);
				recoGenMuR -> Fill(del_R_first_reco_gen);
				
				//Fill 2D Plot for all mu in ghost event
				pt_ratio_deltaR_pt_2d_gen_reco_all_particles -> Fill(del_R_reco_gen, amuon.pt());
				pt_ratio_deltaR_pt_2d_gen_reco_all_particles -> Fill(del_R_first_reco_gen, firstmuon.pt());
				calo_likelihood_ghost_event -> Fill(amuon.caloCompatibility());
				calo_likelihood_ghost_event -> Fill(firstmuon.caloCompatibility());
				ho_energy_ghost_event -> Fill(amuon.calEnergy().ho);
				ho_energy_ghost_event -> Fill(firstmuon.calEnergy().ho);
				ho_energy_ghost_event_3x3 -> Fill(amuon.calEnergy().hoS9);
				ho_energy_ghost_event_3x3 -> Fill(firstmuon.calEnergy().hoS9);

				//Filter for the ghost (see definition at top of file)
				if(del_R_first_reco_gen > 0.1){ //If ghost -> Fill pt ratio
					pt_ratio_ghost_gen -> Fill(firstmuon.pt() / genparticle.pt());
					//Fill correlation plot for ghosts only
					pt_ratio_2d_ghost_gen_reco -> Fill(genparticle.pt(), firstmuon.pt());
					pt_ratio_deltaR_pt_2d_gen_reco -> Fill(del_R_first_reco_gen, firstmuon.pt());
					calo_likelihood_ghost -> Fill(firstmuon.caloCompatibility());
					ho_energy_ghost -> Fill(firstmuon.calEnergy().ho);
					ho_energy_ghost_3x3 -> Fill(firstmuon.calEnergy().hoS9);
				} else{ //If reco mu seems not to be a ghost -> Fill pt ratio of other distribution
					pt_ratio_not_ghost_gen -> Fill(firstmuon.pt() / genparticle.pt());
					//fill gen and reco pt for not ghosts
					pt_ratio_2d_corr_gen_reco -> Fill(genparticle.pt(), firstmuon.pt());
					calo_likelihood_not_ghost -> Fill(firstmuon.caloCompatibility());
					ho_energy_no_ghost -> Fill(firstmuon.calEnergy().ho);
					ho_energy_no_ghost_3x3 -> Fill(firstmuon.calEnergy().hoS9);
				}
				if(del_R_reco_gen > 0.1){ //If ghost -> Fill pt ratio
					pt_ratio_ghost_gen -> Fill(amuon.pt() / genparticle.pt());
					edm::LogInfo ("GhostFound") << "LumiBlock: " << iEvent.id().luminosityBlock() << " Event No: " << iEvent.id().event() << " Pt: " << amuon.pt() << " Pt: " << firstmuon.pt() << " GenPt: " << genparticle.pt() << " DelRecoGen: " << del_R_reco_gen;
				std::cout << "LumiBlock: " << iEvent.id().luminosityBlock() << " Event No: " << iEvent.id().event() << " Pt: " << amuon.pt() << " Pt: " << firstmuon.pt() << " GenPt: " << genparticle.pt() << " vz: " << amuon.vz() << " vz: " << firstmuon.vz() << " vzGen: " << genparticle.vz() << "DelRecoGen: " << del_R_reco_gen << std::endl;
					//Fill correlation plot for ghosts only
					pt_ratio_2d_ghost_gen_reco -> Fill(genparticle.pt(), amuon.pt());
					pt_ratio_deltaR_pt_2d_gen_reco -> Fill(del_R_reco_gen, amuon.pt());
					calo_likelihood_ghost -> Fill(amuon.caloCompatibility());
					ho_energy_ghost -> Fill(amuon.calEnergy().ho);
					ho_energy_ghost_3x3 -> Fill(amuon.calEnergy().hoS9);
				} else { //If reco mu seems not to be a ghost -> Fill pt ratio of other distribution
					pt_ratio_not_ghost_gen -> Fill(amuon.pt() / genparticle.pt());
					//fill gen and reco pt for not ghosts
					pt_ratio_2d_corr_gen_reco -> Fill(genparticle.pt(), amuon.pt());
					calo_likelihood_not_ghost -> Fill(amuon.caloCompatibility()); //Calo likelihood
					ho_energy_no_ghost -> Fill(amuon.calEnergy().ho); //HO deposit
					ho_energy_no_ghost_3x3 -> Fill(amuon.calEnergy().hoS9); //HO deposit
				}
			}
		}

		// Fill the processed particles in the vector
		processed_particles.push_back(genparticle);
		processed_particles_ref.push_back(muonRef);
	}
}

void MuonMatch::beginJob() {
	const double PI  = 3.141592653589793238462;
	
	edm::Service<TFileService> fs;
	//TFileDirectory dir_number_of_ghosts = fs->mkdir("Number_of_ghosts_sum");
	
	recoMuR = fs -> make<TH1D>("Delta R between reco mu within an ghost event", "Delta R between reco mu within ghost event", 1000, 0, 10);
	recoGenMuR = fs -> make<TH1D>("Delta R between reco mu ghosts and gen mu", "Delta R between reco mu ghosts and gen mu", 1000, 0, 10);
	recoGenMuRNoGhost = fs -> make<TH1D>("Delta R between reco mu and matching gen mu", "Delta R between reco mu and matching gen mu", 1000, 0, 10);
	deltaVzGenReco = fs -> make<TH1D>("Distance of the vertices in z direction of reco mu and associated gen mu", "Distance of the vertices in z direction of reco mu and associated gen mu", 2000, -1, 1);
	pt_ratio_ghost_gen = fs -> make<TH1D>("Pt ratio between generated particle and ghost", "Pt ratio between generated particle and ghost", 1000, 0, 10);
	pt_ratio_not_ghost_gen = fs -> make<TH1D>("Pt ratio between generated particle and reco mu in a ghost event", "Pt ratio between generated particle and reco mu in a ghost event", 1000, 0, 10);
	pt_ratio_2d_gen_reco = fs -> make<TH2D>("2D Pt ratio gen reco", "2D Pt ratio gen reco", 200, 0, 100, 250, 0, 100);
	pt_ratio_2d_ghost_gen_reco = fs -> make<TH2D>("2D Pt ratio ghost gen reco", "2D Pt ratio ghost gen reco", 200, 0, 500, 100, 0, 100);
	pt_ratio_2d_corr_gen_reco = fs -> make<TH2D>("2D Pt ratio corr gen reco", "2D Pt ratio corr gen reco", 200, 0, 500, 100, 0, 100);
	deltaRGhostReco = fs -> make<TH1D>("Delta R between both mu matched to one gen", "Delta R between both mu matched to one gen", 1000, 0, 4);
	deltaPtGhostReco = fs -> make<TH1D>("Delta Pt between both mu matched to one gen", "Delta Pt between both mu matched to one gen", 1000, 0, 4);
	pt_ratio_deltaR_pt_2d_gen_reco = fs -> make<TH2D>("2D Delta R vs- pt gen reco", "2D Delta R vs- pt gen reco", 200, 0, 4, 250, 0, 100);
	pt_ratio_deltaR_pt_2d_gen_reco_all_particles = fs -> make<TH2D>("2D Delta R vs- pt gen reco All particles", "2D Delta R vs- pt gen reco All particles", 200, 0, 4, 250, 0, 100);
	calo_likelihood = fs -> make<TH1D>("Likelihood for muon in calo system", "Likelihood for muon in calo system", 1000, 0, 1);
	calo_likelihood_ghost = fs -> make<TH1D>("Likelihood for muon in calo system for ghosts", "Likelihood for muon in calo system for ghosts", 1000, 0, 1);
	calo_likelihood_not_ghost = fs -> make<TH1D>("Likelihood for muon in calo system for not ghosts", "Likelihood for muon in calo system for not ghosts", 1000, 0, 1);
	calo_likelihood_ghost_event = fs -> make<TH1D>("Likelihood for muon in calo system in ghost event", "Likelihood for muon in calo system in ghost event", 1000, 0, 1);
	ho_energy = fs -> make<TH1D>("Energy deposited in HO System", "Energy deposited in HO System", 1000, 0, 10);
	ho_energy_ghost_event = fs -> make<TH1D>("Energy deposited in HO System for ghost events", "Energy deposited in HO System for ghost events", 1000, 0, 10);
	ho_energy_ghost = fs -> make<TH1D>("Energy deposited in HO System for ghosts", "Energy deposited in HO System for ghosts", 1000, 0, 10);
	ho_energy_no_ghost = fs -> make<TH1D>("Energy deposited in HO System for non ghosts", "Energy deposited in HO System for non ghosts", 1000, 0, 10);
	ho_energy_3x3 = fs -> make<TH1D>("Energy deposited in HO System 3x3", "Energy deposited in HO System 3x3", 1000, 0, 10);
	ho_energy_ghost_event_3x3 = fs -> make<TH1D>("Energy deposited in HO System for ghost events 3x3", "Energy deposited in HO System for ghost events 3x3", 1000, 0, 10);
	ho_energy_ghost_3x3 = fs -> make<TH1D>("Energy deposited in HO System for ghosts 3x3", "Energy deposited in HO System for ghosts 3x3", 1000, 0, 10);
	ho_energy_no_ghost_3x3 = fs -> make<TH1D>("Energy deposited in HO System for non ghosts 3x3", "Energy deposited in HO System for non ghosts 3x3", 1000, 0, 10);
	//number_of_ghosts[0] = dir_number_of_ghosts.make<TH1D>("Number of ghosts, pt_min = 0 GeV","Number of ghosts, pt_min = 0 GeV", 4, -0.5, 3.5);
}

void MuonMatch::endJob(){
}

//Checks if the GenParticles are equal in condition of charge, p_t
bool MuonMatch::equalsGenParticles(reco::GenParticle gp1, reco::GenParticle gp2){
	if(gp1.pt() != gp2.pt()){
		return false;
	}
	if(gp1.charge() != gp2.charge()){
		return false;
	}
	if(gp1.px() != gp2.px()){
		return false;
	}
	if(gp1.py() != gp2.py()){
		return false;
	}
	if(gp1.vx() != gp2.vx()){
		return false;
	}
	return true;
}

// ------------ method called when starting to processes a run  ------------
void MuonMatch::beginRun(edm::Run const&, edm::EventSetup const& iSetup){
	//propagatorToMuon.init(iSetup);
}

// ------------ method called when ending the processing of a run  ------------
void MuonMatch::endRun(edm::Run const&, edm::EventSetup const&){
}

// ------------ method called when starting to processes a luminosity block  ------------
void MuonMatch::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}

// ------------ method called when ending the processing of a luminosity block  ------------
void MuonMatch::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonMatch::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonMatch);

