// -*- C++ -*-
//
// Package:    Z2mumu
// Class:      Z2mumu
// 
/**\class Z2mumu Z2mumu.cc Gen/Z2mumu/src/Z2mumu.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Aleena Rafique
//         Created:  Tue Dec 11 18:55:26 PKT 2012
// $Id$
//
//


// system include files
#include <memory>
#include <vector>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "DataFormats/Common/interface/Handle.h"
#include "TStyle.h"
#include "math.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFile.h" // ROOT, for saving file
#include "TVirtualPad.h" // ROOT, for interactive graphics
#include "TApplication.h" // ROOT , for interactive graphics
//
// class declaration
//

class Z2mumu : public edm::EDAnalyzer {
   public:
      explicit Z2mumu(const edm::ParameterSet&);
      ~Z2mumu();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        TH1D* charge_mult_mu;
        TH1D* charge_mu;
        TH1D* px_mu;
        TH1D* py_mu;
        TH1D* pz_mu;
        TH1D* e_mu;
        TH1D* mass_mu;
        TH1D* pT_mu;
        TH1D* eta_mu;
        TH1D* phi_mu;
        TH1D* pT_mu_from_px_py;
        TH1D* inv_mass_mu;
        TH2D* eta_vs_energy;
        TH2D* eta_vs_pT;
        TH1D* delta_eta;
        TH1D* delta_phi;
        TH1D* inv_Z_mass;
        TH1D* pT_Z;
        TH2D* twoD_eta;
        TH2D* twoD_phi;
        TH2D* twoD_pT;

      // ----------member data ---------------------------
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
Z2mumu::Z2mumu(const edm::ParameterSet& iConfig):charge_mu(0),charge_mult_mu(0),px_mu(0),py_mu(0),pz_mu(0),e_mu(0),mass_mu(0),pT_mu(0),eta_mu(0),phi_mu(0),pT_mu_from_px_py(0),inv_mass_mu(0),eta_vs_energy(0),eta_vs_pT(0),delta_eta(0),delta_phi(0),inv_Z_mass(0),pT_Z(0),twoD_eta(0),twoD_phi(0),twoD_pT(0)

{
   //now do what ever initialization is needed

}


Z2mumu::~Z2mumu()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Z2mumu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace HepMC;
   using namespace reco;
   using namespace std;
//   using namespace RooFit ;

Handle< std::vector<reco::GenParticle> > EvtHandle ;
iEvent.getByLabel(edm::InputTag("genParticles"),EvtHandle);
std::vector<reco::GenParticle>::const_iterator i;
        int nChg = 0;
        double eta1=0,eta2=0,phi1=0,phi2=0,e1=0,e2=0,px1=0,px2=0,py1=0,py2=0,pz1=0,pz2=0,pT1=0,pT2=0;
        for(i=EvtHandle->begin();i!=EvtHandle->end();++i)
        {
                const Candidate * d = (*i).mother();
                if ((*i).pdgId()==13 || (*i).pdgId()==-13)
                if (d->pdgId()==23)
//              if ((*i).pt()>20 && fabs((*i).eta())<2.1)
                {
                        ++nChg;
                        double charge_muons=((*i).charge());
                        double px_muons=((*i).px());
                        double py_muons=((*i).py());
                        double pz_muons=((*i).pz());
                        double e_muons=((*i).energy());
                        double mass_muons=((*i).mass());
                        double pT_muons=((*i).pt());
                        double eta_muons=((*i).eta());
                        double phi_muons=((*i).phi());
                        double sq_pT_mu = ((px_muons*px_muons)+(py_muons*py_muons));
                        double pT_from_px_py = sqrt(sq_pT_mu);
                        double sq_mass_mu = ((e_muons*e_muons)-((px_muons*px_muons)+(py_muons*py_muons)+(pz_muons*pz_muons)));
                        double invMass = sqrt(sq_mass_mu);

                        charge_mu->Fill(charge_muons);
                        px_mu->Fill(px_muons);
                        py_mu->Fill(py_muons);
                        pz_mu->Fill(pz_muons);
                        e_mu->Fill(e_muons);
                        mass_mu->Fill(mass_muons);
                        pT_mu->Fill(pT_muons);
                        eta_mu->Fill(eta_muons);
                        phi_mu->Fill(phi_muons);
                        pT_mu_from_px_py->Fill(pT_from_px_py);
                        inv_mass_mu->Fill(invMass);
                        eta_vs_energy->Fill(eta_muons,e_muons);
                        eta_vs_pT->Fill(eta_muons,pT_muons);

                }
		if ((*i).pdgId()==13 && d->pdgId()==23)
		{
                        eta1 = ((*i).eta());
                        phi1 = ((*i).phi());
                        e1 = ((*i).energy());
                        px1 = ((*i).px());
                        py1 = ((*i).py());
                        pz1 = ((*i).pz());
                        pT1 = ((*i).pt());
		}
                 if ((*i).pdgId()==-13 && d->pdgId()==23)
                {
                        eta2 = ((*i).eta());
                        phi2 = ((*i).phi());
                        e2 = ((*i).energy());
                        px2 = ((*i).px());
                        py2 = ((*i).py());
                        pz2 = ((*i).pz());
                        pT2 = ((*i).pt());
                }
        }
        double del_eta = (eta1 - eta2);
        double del_phi = (phi1 - phi2);
        double sq_Z_mass = (((e1+e2)*(e1+e2))-(((px1+px2)*(px1+px2))+((py1+py2)*(py1+py2))+((pz1+pz2)*(pz1+pz2))));
        double invZmass = sqrt(sq_Z_mass);
        double sq_pT_Z = (((px1+px2)*(px1+px2))+((py1+py2)*(py1+py2)));
        double Z_pT = sqrt(sq_pT_Z);
        charge_mult_mu->Fill(nChg);
        delta_eta->Fill(del_eta);
        delta_phi->Fill(del_phi);
        if(del_phi > 1.)
        {
        if(invZmass >= 78.0 && invZmass <102.0)
        inv_Z_mass->Fill(invZmass);
        pT_Z      ->Fill(Z_pT);
        }
        twoD_eta  ->Fill(eta1,eta2);
        twoD_phi  ->Fill(phi1,phi2);
        twoD_pT   ->Fill(pT1,pT2);



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
Z2mumu::beginJob()
{
        edm::Service<TFileService> fs;
        charge_mu = fs->make<TH1D> ("charge_mu","charge of muons",100,-2.,2.);
        charge_mult_mu = fs->make<TH1D> ("charge_mult_mu","muon charged multiplicity",100,0.,10.);
        px_mu = fs->make<TH1D> ("px_mu","px of muons",100,-100.,100.);
        py_mu = fs->make<TH1D> ("py_mu","py of muons",100,-100.,100.);
        pz_mu = fs->make<TH1D> ("pz_mu","pz of muons",100,-100.,100.);
        e_mu = fs->make<TH1D> ("e_mu","energy of muons",100,0.,200.);
        mass_mu = fs->make<TH1D> ("mass_mu","mass of muons",100,0.05,0.15);
        pT_mu = fs->make<TH1D> ("pT_mu","pT of muons",100,0.,100.);
        eta_mu = fs->make<TH1D> ("eta_mu","pseudorapidity plot for muons",100,-2.5,2.5);
        phi_mu = fs->make<TH1D> ("phi_mu","Azimuthal angle plot for muons",100,-3.14,3.14);
        pT_mu_from_px_py = fs->make<TH1D> ("pT_mu_using_px_py","pT plot for muons from px and py",100,0.,100.);
        inv_mass_mu = fs->make<TH1D> ("inv_mass_mu","invariant muon mass using 4-momentum",100,0.05,0.15);
        eta_vs_energy = fs->make<TH2D> ("eta_vs_energy","2D plot eta vs energy of muons",100,-2.5,2.5,100,0.,200.);
        eta_vs_pT = fs->make<TH2D> ("eta_vs_pTmu","2D plot eta vs pT of muons",100,-2.5,2.5,100,0.,100.);
         delta_eta = fs->make<TH1D> ("delta_eta","pseudorapidity difference plot for muons",100,-2.5,2.5);
         delta_phi = fs->make<TH1D> ("delta_phi","delta phi plot for muons",100,-3.14,3.14);
         inv_Z_mass = fs->make<TH1D> ("inv_Z_mass","invariant mass of Z boson plot",100,60.,120.);
         pT_Z = fs->make<TH1D> ("pT_Z","pT Z plot",100,0.,100.);
         twoD_eta = fs->make<TH2D> ("twoD_eta","2D plot for eta1 and eta2",100,-2.5,2.5,100,-2.5,2.5);
         twoD_phi = fs->make<TH2D> ("twoD_phi","2D plot for phi1 and phi2",100,-3.14,3.14,100,-3.14,3.14);
         twoD_pT = fs->make<TH2D> ("twoD_pT","2D plot for pT1 and pT2",100,-20.,100.,100,-20.,100.);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Z2mumu::endJob() 
{
        charge_mult_mu->GetXaxis()->SetTitle("charge multiplicity of muons");
        charge_mult_mu->GetYaxis()->SetTitle("Entries/bin");
        charge_mu->GetXaxis()->SetTitle("charge of muons");
        charge_mu->GetYaxis()->SetTitle("Entries/bin");
        px_mu->GetXaxis()->SetTitle("px of muons[Gev/c]");
        px_mu->GetYaxis()->SetTitle("Entries/bin");
        py_mu->GetXaxis()->SetTitle("py of muons[Gev/c]");
        py_mu->GetYaxis()->SetTitle("Entries/bin");
        pz_mu->GetXaxis()->SetTitle("pz of muons[Gev/c]");
        pz_mu->GetYaxis()->SetTitle("Entries/bin");
        e_mu->GetXaxis()->SetTitle("energy of muons[Gev]");
        e_mu->GetYaxis()->SetTitle("Entries/bin");
        mass_mu->GetXaxis()->SetTitle("mass of muons[Gev/c^2]");
        mass_mu->GetYaxis()->SetTitle("Entries/bin");
        pT_mu->GetXaxis()->SetTitle("pT of muons[Gev/c]");
        pT_mu->GetYaxis()->SetTitle("Entries/bin");
        eta_mu->GetXaxis()->SetTitle("eta of muons");
        eta_mu->GetYaxis()->SetTitle("Entries/bin");
        phi_mu->GetXaxis()->SetTitle("phi angle for muons");
        phi_mu->GetYaxis()->SetTitle("Entries/bin");
        pT_mu_from_px_py->GetXaxis()->SetTitle("pT of muons[Gev/c]");
        pT_mu_from_px_py->GetYaxis()->SetTitle("Entries/bin");
        inv_mass_mu->GetXaxis()->SetTitle("mass of muons[Gev/c^2]");
        inv_mass_mu->GetYaxis()->SetTitle("Entries/bin");
        eta_vs_energy->GetXaxis()->SetTitle("eta of muons");
        eta_vs_energy->GetYaxis()->SetTitle("energy of muons[Gev]");
        eta_vs_pT->GetXaxis()->SetTitle("eta of muons");
        eta_vs_pT->GetYaxis()->SetTitle("pT of muons[Gev/c]");
        delta_eta->GetXaxis()->SetTitle("delta eta for muons");
        delta_eta->GetYaxis()->SetTitle("Entries/bin");
        delta_phi->GetXaxis()->SetTitle("delta phi for muons");
        delta_phi->GetYaxis()->SetTitle("Entries/bin");
        inv_Z_mass->GetXaxis()->SetTitle("invariant mass of Z boson[Gev/c^2]");
        inv_Z_mass->GetYaxis()->SetTitle("Entries/bin");
        pT_Z->GetXaxis()->SetTitle("pT of Z boson[Gev/c]");
        pT_Z->GetYaxis()->SetTitle("Entries/bin");
        twoD_eta->GetXaxis()->SetTitle("eta for muons");
        twoD_eta->GetYaxis()->SetTitle("eta for antimuons");
        twoD_phi->GetXaxis()->SetTitle("phi angle for muons");
        twoD_phi->GetYaxis()->SetTitle("phi angle for antimuons");
        twoD_pT->GetXaxis()->SetTitle("pT for muons[Gev/c]");
        twoD_pT->GetYaxis()->SetTitle("pT for antimuons[Gev/c]");
}

// ------------ method called when starting to processes a run  ------------
void 
Z2mumu::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Z2mumu::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Z2mumu::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Z2mumu::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Z2mumu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Z2mumu);
