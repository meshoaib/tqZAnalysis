// -*- C++ -*-
//
// Package:    Z2ep
// Class:      Z2ep
// 
/**\class Z2ep Z2ep.cc EHEP_Start/Z2ep/src/Z2ep.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Muhammad Shoaib
//         Created:  Wed Dec 11 10:37:26 PKT 2013
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
//-------------------------------
#include "TFile.h"
#include "TH1F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//
// class declaration
//

class Z2ep : public edm::EDAnalyzer {
   public:
      explicit Z2ep(const edm::ParameterSet&);
      ~Z2ep();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze (const edm::Event&, const edm::EventSetup&);
      virtual void endJob  () ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun  (edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock  (edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      // class TH1D: public TH1, public TArrayD
      // TH1D : histograms with one double per channel. Maximum precision 14 digits
      // ---------------------------------------------------------------------------
      // following are pointers of type TH1D
      //      
        TH1D* charge_mult_elec;
        TH1D* charge_elec;
        TH1D* px_elec;
        TH1D* py_elec;
        TH1D* pz_elec;
        TH1D* e_elec;
        TH1D* mass_elec;
        TH1D* pT_elec;
        TH1D* eta_elec;
        TH1D* phi_elec;
        TH1D* pT_elec_from_px_py;
        TH1D* inv_mass_elec;
        TH2D* eta_vs_energy;
        TH2D* eta_vs_pT;
        TH1D* delta_eta;
        TH1D* delta_phi;
        TH1D* inv_Z_mass;
        TH1D* pT_Z;
        TH2D* twoD_eta;
        TH2D* twoD_phi;
        TH2D* twoD_pT;
      // -------------------------
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
Z2ep::Z2ep(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


Z2ep::~Z2ep()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Z2ep::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace HepMC;
    using namespace reco;
    using namespace std;
//  using namespace RooFit ;

    Handle< std::vector<reco::GenParticle> > EvtHandle ;
    iEvent.getByLabel(edm::InputTag("genParticles"),EvtHandle);
    std::vector<reco::GenParticle>::const_iterator j;
        int    nChg =0;
        double  eta1=0,
                eta2=0,
                phi1=0,
                phi2=0,
                e1  =0,
                e2  =0,
                px1 =0,
                px2 =0,
                py1 =0,
                py2 =0,
                pz1 =0,
                pz2 =0,
                pT1 =0,
                pT2 =0;
        for(j=EvtHandle->begin();j!=EvtHandle->end();++j)
        {
                const Candidate * d = (*j).mother();
                if ((*j).pdgId()==11 || (*j).pdgId()==-11)
                if (d->pdgId()==23)
//              if ((*j).pt()>20 && fabs((*j).eta())<2.1)
                {
                ++nChg;
                double charge_electrons = ((*j).charge());
                double px_electrons     = ((*j).px());
                double py_electrons     = ((*j).py());
                double pz_electrons     = ((*j).pz());
                double e_electrons      = ((*j).energy());
                double mass_electrons   = ((*j).mass());
                double pT_electrons     = ((*j).pt());
                double eta_electrons    = ((*j).eta());
                double phi_electrons    = ((*j).phi());
                double sq_pT_elec       = ((px_electrons*px_electrons)
                                         + (py_electrons*py_electrons));
                double pT_from_px_py    = sqrt(sq_pT_elec);
                double sq_mass_elec     = ((e_electrons*e_electrons)
                                          -((px_electrons*px_electrons)
                                          + (py_electrons*py_electrons)
                                          + (pz_electrons*pz_electrons)));
                double invMass          = sqrt(sq_mass_elec);

                charge_elec        ->Fill(charge_electrons);
                px_elec            ->Fill(px_electrons);
                py_elec            ->Fill(py_electrons);
                pz_elec            ->Fill(pz_electrons);
                e_elec             ->Fill(e_electrons);
                mass_elec          ->Fill(mass_electrons);
                pT_elec            ->Fill(pT_electrons);
                eta_elec           ->Fill(eta_electrons);
                phi_elec           ->Fill(phi_electrons);
                pT_elec_from_px_py ->Fill(pT_from_px_py);
                inv_mass_elec      ->Fill(invMass);
                eta_vs_energy      ->Fill(eta_electrons,e_electrons);
                eta_vs_pT          ->Fill(eta_electrons,pT_electrons);

                }
                    if ((*j).pdgId()==11 && d->pdgId()==23)
                    {
                        eta1 = ((*j).eta());
                        phi1 = ((*j).phi());
                        e1   = ((*j).energy());
                        px1  = ((*j).px());
                        py1  = ((*j).py());
                        pz1  = ((*j).pz());
                        pT1  = ((*j).pt());
                    }
                    if ((*j).pdgId()==-11 && d->pdgId()==23)
                    {
                        eta2 = ((*j).eta());
                        phi2 = ((*j).phi());
                        e2   = ((*j).energy());
                        px2  = ((*j).px());
                        py2  = ((*j).py());
                        pz2  = ((*j).pz());
                        pT2  = ((*j).pt());
                    }
        }
        double del_eta   = (eta1 - eta2);
        double del_phi   = (phi1 - phi2);
        double sq_Z_mass = (((e1+e2)*(e1+e2))-(((px1+px2)*(px1+px2))+((py1+py2)*(py1+py2))+((pz1+pz2)*(pz1+pz2))));
        double invZmass  = sqrt(sq_Z_mass);
        double sq_pT_Z   = (((px1+px2)*(px1+px2))+((py1+py2)*(py1+py2)));
        double Z_pT      = sqrt(sq_pT_Z);
        charge_mult_elec -> Fill(nChg);
        delta_eta        -> Fill(del_eta);
        delta_phi        -> Fill(del_phi);
        if(del_phi > 1.)
        {
        if(invZmass >= 78.0 && invZmass <102.0)
        inv_Z_mass       -> Fill(invZmass);
        pT_Z             -> Fill(Z_pT);
        }
        twoD_eta         -> Fill(eta1,eta2);
        twoD_phi         -> Fill(phi1,phi2);
        twoD_pT          -> Fill(pT1,pT2);




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
Z2ep::beginJob()
{
edm::Service<TFileService> fs;

    charge_elec        = fs->make<TH1D> ("charge_elec","charge of electrons",100,-2.,2.);
    charge_mult_elec   = fs->make<TH1D> ("charge_mult_elec","electron charged multiplicity",100,0.,10.);
    px_elec            = fs->make<TH1D> ("px_elec","px of electrons",100,0.,100.);
    py_elec            = fs->make<TH1D> ("py_elec","py of electrons",100,0.,100.);
    pz_elec            = fs->make<TH1D> ("pz_elec","pz of electrons",100,0.,100.);
    e_elec             = fs->make<TH1D> ("e_elec","energy of electrons",100,0.,100.);
    mass_elec          = fs->make<TH1D> ("mass_elec","mass of electrons",100,0.0001,0.001);
    pT_elec            = fs->make<TH1D> ("pT_elec","pT of electrons",100,0.,100.);
    eta_elec           = fs->make<TH1D> ("eta_elec","pseudorapidity plot for electrons",100,-2.5,2.5);
    phi_elec           = fs->make<TH1D> ("phi_elec","Azimuthal angle plot for electrons",100,-3.14,3.14);
    pT_elec_from_px_py = fs->make<TH1D> ("pT_elec_using_px_py","pT plot for electrons from px and py",100,0.,100.);
    inv_mass_elec      = fs->make<TH1D> ("inv_mass_elec","invariant electron mass using 4-momentum",100,0.0001,0.001);
    eta_vs_energy      = fs->make<TH2D> ("eta_vs_energy","2D plot eta vs energy of electrons",100,-2.5,2.5,100,0.,100.);
    eta_vs_pT          = fs->make<TH2D> ("eta_vs_pTelec","2D plot eta vs pT of electrons",100,-2.5,2.5,100,0.,100.);
    delta_eta          = fs->make<TH1D> ("delta_eta","pseudorapidity difference plot for electrons",100,-2.5,2.5);
    delta_phi          = fs->make<TH1D> ("delta_phi","delta phi plot for electrons",100,-3.14,3.14);
    inv_Z_mass         = fs->make<TH1D> ("inv_Z_mass","invariant mass of Z boson plot",100,60.,120.);
    pT_Z               = fs->make<TH1D> ("pT_Z","pT Z plot",100,0.,100.);
    twoD_eta           = fs->make<TH2D> ("twoD_eta","2D plot for eta1 and eta2",100,-2.5,2.5,100,-2.5,2.5);
    twoD_phi           = fs->make<TH2D> ("twoD_phi","2D plot for phi1 and phi2",100,-3.14,3.14,100,-3.14,3.14);
    twoD_pT            = fs->make<TH2D> ("twoD_pT","2D plot for pT1 and pT2",100,-2.,100.,100,-2.,100.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Z2ep::endJob() 
{
    charge_mult_elec   ->GetXaxis()->SetTitle("charge multiplicity of electrons");
    charge_mult_elec   ->GetYaxis()->SetTitle("Entries/bin");
    charge_elec        ->GetXaxis()->SetTitle("charge of electrons");
    charge_elec        ->GetYaxis()->SetTitle("Entries/bin");
    px_elec            ->GetXaxis()->SetTitle("px of electrons[Gev/c]");
    px_elec            ->GetYaxis()->SetTitle("Entries/bin");
    py_elec            ->GetXaxis()->SetTitle("py of electrons[Gev/c]");
    py_elec            ->GetYaxis()->SetTitle("Entries/bin");
    pz_elec            ->GetXaxis()->SetTitle("pz of electrons[Gev/c]");
    pz_elec            ->GetYaxis()->SetTitle("Entries/bin");
    e_elec             ->GetXaxis()->SetTitle("energy of electrons[Gev]");
    e_elec             ->GetYaxis()->SetTitle("Entries/bin");
    mass_elec          ->GetXaxis()->SetTitle("mass of electrons[Gev/c^2]");
    mass_elec          ->GetYaxis()->SetTitle("Entries/bin");
    pT_elec            ->GetXaxis()->SetTitle("pT of electrons[Gev/c]");
    pT_elec            ->GetYaxis()->SetTitle("Entries/bin");
    eta_elec           ->GetXaxis()->SetTitle("eta of electrons");
    eta_elec           ->GetYaxis()->SetTitle("Entries/bin");
    phi_elec           ->GetXaxis()->SetTitle("phi angle for electrons");
    phi_elec           ->GetYaxis()->SetTitle("Entries/bin");
    pT_elec_from_px_py ->GetXaxis()->SetTitle("pT of electrons[Gev/c]");
    pT_elec_from_px_py ->GetYaxis()->SetTitle("Entries/bin");
    inv_mass_elec      ->GetXaxis()->SetTitle("mass of electrons[Gev/c^2]");
    inv_mass_elec      ->GetYaxis()->SetTitle("Entries/bin");
    eta_vs_energy      ->GetXaxis()->SetTitle("eta of electrons");
    eta_vs_energy      ->GetYaxis()->SetTitle("energy of electrons[Gev]");
    eta_vs_pT          ->GetXaxis()->SetTitle("eta of electrons");
    eta_vs_pT          ->GetYaxis()->SetTitle("pT of electrons[Gev/c]");
    delta_eta          ->GetXaxis()->SetTitle("delta eta for electrons");
    delta_phi          ->GetXaxis()->SetTitle("delta phi for electrons");
    delta_phi          ->GetYaxis()->SetTitle("Entries/bin");
    inv_Z_mass         ->GetXaxis()->SetTitle("invariant mass of Z boson[Gev/c^2]");
    inv_Z_mass         ->GetYaxis()->SetTitle("Entries/bin");
    pT_Z               ->GetXaxis()->SetTitle("pT of Z boson[Gev/c]");
    pT_Z               ->GetYaxis()->SetTitle("Entries/bin");
    twoD_eta           ->GetXaxis()->SetTitle("eta for electrons");
    twoD_eta           ->GetYaxis()->SetTitle("eta for positrons");
    twoD_phi           ->GetXaxis()->SetTitle("phi angle for electrons");
    twoD_phi           ->GetYaxis()->SetTitle("phi angle for positrons");
    twoD_pT            ->GetXaxis()->SetTitle("pT for electrons[Gev/c]");
    twoD_pT            ->GetYaxis()->SetTitle("pT for positrons[Gev/c]");

}

// ------------ method called when starting to processes a run  ------------
void 
Z2ep::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Z2ep::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Z2ep::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Z2ep::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Z2ep::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Z2ep);
