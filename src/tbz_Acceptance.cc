// -*- C++ -*-
//
// Package:    tbz_Acceptance
// Class:      tbz_Acceptance
//
/**\class tbz_Acceptance tbz_Acceptance.cc tbz_Acceptance/tbz_Acceptance/src/tbz_Acceptance.cc

 Description: [one line class summary]
 1 lepton and b-quark for w reconstruction (only lepton decay of W is to be studied)

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Muhammad Shoaib 
//      Start date:  Sunday, 29/01/2014
//         Created:  _____________________________
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/BTauReco/interface/JetTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "CLHEP/Vector/EulerAngles.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "math.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TLorentzVector.h"

//
// class declaration
//

class tbz_Acceptance : public edm::EDAnalyzer {
   public:
      explicit tbz_Acceptance(const edm::ParameterSet&);
     ~tbz_Acceptance();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      
      //---- histograms------
      TH1D* H1_pt_elec_Z    ;
      TH1D* H1_eta_elec_Z   ;
      TH1D* H1_pt_positrn_Z ;
      TH1D* H1_eta_positrn_Z;
      TH1D* H1_pt_muon1_Z   ;
      TH1D* H1_eta_muon1_Z ;
      TH1D* H1_pt_muon2_Z   ;
      TH1D* H1_eta_muon2_Z ;
      TH1D* H1_wmu_pt_W     ;
      TH1D* H1_wmu_eta_W    ;
      TH1D* H1_wNu_pt       ;
      TH1D* H1_wNu_eta      ;
      TH1D* H1_wElec_pt     ;
      TH1D* H1_wElec_eta    ;
      TH1D* H1_wElec_Nu_pt  ;
      TH1D* H1_wElec_Nu_eta  ;
      
      TH1D* H1_welec_pt      ;
      TH1D* H1_welec_eta     ;
      
      TH1D* H1_wu_pt         ;
      TH1D* H1_wu_eta        ;
      
      TH1D* H1_zuu_pt        ;
      TH1D* H1_zuu_eta       ;
      
      TH1D* H1_zee_pt        ;
      TH1D* H1_zee_eta       ;
      //-
      TH1D* H1_AllZmuons_pt   ;
      TH1D* H1_AllZmuons_eta  ;
      TH1D* H1_AllZelec_pt    ;
      TH1D* H1_AllZelec_eta   ;
      
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
tbz_Acceptance::tbz_Acceptance(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
}


tbz_Acceptance::~tbz_Acceptance()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
//
// member functions
//

// ------------ method called for each event  ------------
void
tbz_Acceptance::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{

      using namespace  std   ;
      using namespace  edm   ;
      using namespace  pat   ;
      using namespace  reco  ;
      
      //Gen particles loop
      Handle<GenParticleCollection> genParticle      ;             
      iEvent.getByLabel("genParticles",genParticle)  ;
      
      // bool isZ=false;
      // bool isW=false;
      
      //----------
      TLorentzVector elec_from_Z       (0,0,0,0)    ;  
      TLorentzVector positrn_from_Z    (0,0,0,0)    ;
      TLorentzVector elec_from_W       (0,0,0,0)    ;
      TLorentzVector positrn_from_W    (0,0,0,0)    ;
      TLorentzVector wElec_nuerino     (0,0,0,0)    ;
      TLorentzVector wElec             (0,0,0,0)    ;
      TLorentzVector wMuons_nuetrino   (0,0,0,0)    ;
      TLorentzVector wMuons            (0,0,0,0)    ;
      TLorentzVector muon1             (0,0,0,0)    ;
      TLorentzVector muon2             (0,0,0,0)    ;
      TLorentzVector Zee_Inv_Mass      (0,0,0,0)    ;
      TLorentzVector Zee_cand          (0,0,0,0)    ;
      TLorentzVector Zuu_cand          (0,0,0,0)    ;
      TLorentzVector Wu_cand           (0,0,0,0)    ;
      TLorentzVector Welec_cand        (0,0,0,0)    ;
      
      

   for(size_t i = 0; i < genParticle->size(); ++ i)
   {
     if( genParticle->size() < 2.) continue ;
       const GenParticle & p = (*genParticle)[i]    ; 
       // int id = p.pdgId()                        ;      
       if(abs(p.pdgId()) ==11 && p.mother()->pdgId()==23)
       {
          double AllElec_pt                   ;
          double AllElec_eta                  ;
          AllElec_pt     = p.pt()             ;
          H1_AllZelec_pt ->Fill(AllElec_pt)   ;
          AllElec_eta    = p.eta()            ;
          H1_AllZelec_eta ->Fill(AllElec_eta) ;      
       
       
         if(p.pdgId()==11 && p.mother()->pdgId()==23) // electrons for Z
         {
           double pt_elec =   0.            ;      
           double eta_elec = -100.          ;
           pt_elec  = p.pt()                ;
           H1_pt_elec_Z->Fill(pt_elec)      ;
           eta_elec = p.eta()               ;
           H1_eta_elec_Z->Fill(eta_elec)    ;
           elec_from_Z.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy()) ;
         }
         if(p.pdgId()== -11 && p.mother()->pdgId()==23) //positrons
         {
         
          double pt_positrn  = 0.                   ;
          double eta_positrn = 0.                   ;
          pt_positrn         = p.pt()               ;
          H1_pt_positrn_Z       ->Fill(pt_positrn)  ;
          eta_positrn        = p.eta()              ;
          H1_eta_positrn_Z      ->Fill(eta_positrn) ;
          
          positrn_from_Z.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy());
          
         }
       }
       
       
       if(abs(p.pdgId()) ==13 && p.mother()->pdgId() ==23 )
       {
         double AllMuons_pt                   ;
         double AllMuons_eta                  ;         
         AllMuons_pt     = p.pt()             ;
         H1_AllZmuons_pt ->Fill(AllMuons_pt)  ;
         AllMuons_eta    = p.eta()            ;
         H1_AllZmuons_eta ->Fill(AllMuons_eta);
         
         if(p.pdgId() ==13 && p.mother()->pdgId() ==23) // Muons for Z
         {         
          double pt_muon1  =  0.                                 ;
          double eta_muon1 = -100                                ;
          pt_muon1         = p.pt()                              ;
          H1_pt_muon1_Z       ->Fill(pt_muon1)                   ;
          eta_muon1        = p.eta()                             ;
          H1_eta_muon1_Z      ->Fill(eta_muon1)                  ;
                                                                 
          muon1.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy())   ;
         }
         if( p.pdgId() == -13 && p.mother()->pdgId() ==23)
         {
          double pt_mu2    =   0.                                ;
          double eta_mu2   = -100.                               ; 
          pt_mu2           = p.pt()                              ;
          H1_pt_muon2_Z       ->Fill(pt_mu2)                     ;
          eta_mu2          = p.eta()                             ;
          H1_eta_muon2_Z      ->Fill(eta_mu2)                    ;          
          muon2.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy())   ;          
         }
       } 
       
         if(abs(p.pdgId()) ==13 || abs(p.pdgId()) ==14 && abs(p.mother()->pdgId()) ==24 ) // Muons for W
         {
            if(abs(p.pdgId()) ==13 && abs(p.mother()->pdgId()) ==24 )
             {
             
             double wmu_pt  =  0.                                    ;
             double wmu_eta   = -100.                                ;
                                                                     
             wmu_pt         = p.pt()                                 ;
             H1_wmu_pt_W       ->Fill(wmu_pt)                        ;
             wmu_eta          = p.eta()                              ;
             H1_wmu_eta_W      ->Fill(wmu_eta)                       ;
                                                                     
             wMuons.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy())   ;
             
             }
             if(abs(p.pdgId()) ==14 && abs(p.mother()->pdgId()) ==24 )
             {
             
             double wNu_pt      =  0.                                        ;
             double wNu_eta     = -100.                                      ;                                                                             
             wNu_pt             = p.pt()                                     ;
             H1_wNu_pt          ->Fill(wNu_pt)                               ;
             wNu_eta            = p.eta()                                    ;
             H1_wNu_eta          ->Fill(wNu_eta)                             ;                                                                             
             wMuons_nuetrino.SetPxPyPzE( p.px(), p.py(), p.pz(), p.energy() )  ;             
             }
         }
         
         if(abs(p.pdgId()) ==11 || abs(p.pdgId()) ==12 && abs(p.mother()->pdgId()) ==24) 
         {
           if( abs(p.pdgId()) ==11 && abs(p.mother()->pdgId()) ==24 )
             {             
             double wElec_pt       =  0.            ;
             double wElec_eta      = -100.          ;
             
             wElec_pt              = p.pt()         ;
             H1_wElec_pt           ->Fill(wElec_pt );
             wElec_eta             = p.eta()        ;             
             H1_wElec_eta          ->Fill(wElec_eta);
             
             wElec.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy());             
             }
            
            if( abs(p.pdgId()) == 12 && abs(p.mother()->pdgId()) == 24 )
             {             
             double wElec_Nu_pt    = 0.                  ;
             double wElec_Nu_eta   = -100.               ;                                                         
             wElec_Nu_pt           = p.pt()              ;
             H1_wElec_Nu_pt        ->Fill(wElec_Nu_pt )  ;
             wElec_Nu_eta          = p.eta()             ;             
             H1_wElec_Nu_eta       ->Fill(wElec_Nu_eta)  ;   
             
             wElec_nuerino.SetPxPyPzE( p.px(), p.py(), p.pz(), p.energy() );
             
             }
         }                 
   }//genParticle loop ends 
   double zee_pt                             ;
   double zee_eta                            ;
   double zuu_pt                             ;
   double zuu_eta                            ;
   double wu_pt                              ;
   double wu_eta                             ;
   double welec_pt                           ;
   double welec_eta                          ;
   
   Zee_cand    = elec_from_Z + positrn_from_Z ;
   
   zee_pt      = Zee_cand.Pt()                ;
   if(zee_pt > 0.)
   {
   H1_zee_pt   ->Fill(zee_pt )                ;
   zee_eta     = Zee_cand.Eta()               ;
   H1_zee_eta  ->Fill(zee_eta)                ;
   }
 
   Zuu_cand    = muon1 + muon2                ;
   zuu_pt      = Zuu_cand.Pt()                ;
   if(zuu_pt > 0.)
   {
   zuu_eta     = Zuu_cand.Eta()               ;
   H1_zuu_pt   ->Fill( zuu_pt )               ; 
   H1_zuu_eta  ->Fill( zuu_eta)               ; 
   }
   Wu_cand     = wMuons + wMuons_nuetrino     ;
   wu_pt       = Wu_cand.Pt()                 ;
   if(wu_pt > 0.)
   {
   wu_eta      = Wu_cand.Eta()                ;
   H1_wu_pt    ->Fill( wu_pt  )               ;
   H1_wu_eta   ->Fill( wu_eta )               ;
   }
   
   Welec_cand  = wElec + wElec_nuerino        ;
   welec_pt    = Welec_cand.Pt()              ;
   if(welec_pt > 0.)
   {
   welec_eta   = Welec_cand.Eta()             ;
   H1_welec_pt   ->Fill(welec_pt  )           ;
   H1_welec_eta  ->Fill(welec_eta )           ;
   }
   
}


   #ifdef THIS_IS_AN_EVENT_EXAMPLE
      Handle<ExampleData> pIn;
      iEvent.getByLabel("example",pIn);
   #endif

   #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
      ESHandle<SetupData> pSetup;
      iSetup.get<SetupRecord>().get(pSetup);
   #endif



// ------------ method called once each job just before starting event loop  ------------
void
tbz_Acceptance::beginJob()
{

edm::Service<TFileService> fs;
H1_pt_elec_Z     = fs->make<TH1D>("H1_pt_elec_Z","H1_pt_elec_Z", 100, 0., 500.);
H1_eta_elec_Z    = fs->make<TH1D>("H1_eta_elec_Z", "H1_eta_elec_Z", 100, -6., 6.);
H1_pt_positrn_Z  = fs->make<TH1D>("H1_pt_positrn_Z", "H1_pt_positrn_Z", 100, 0., 500.);
H1_eta_positrn_Z = fs->make<TH1D>("H1_eta_positrn_Z","H1_eta_positrn_Z", 100, -6., 6.);
H1_pt_muon1_Z    = fs->make<TH1D>("H1_pt_muon1_Z", "H1_pt_muon1_Z", 100, 0., 500.);
H1_eta_muon1_Z   = fs->make<TH1D>("H1_etat_muon1_Z", "H1_etat_muon1_Z", 100, -6., 6.);
H1_pt_muon2_Z    = fs->make<TH1D>("H1_pt_muon2_Z", "H1_pt_muon2_Z", 100, 0., 500.);
H1_eta_muon2_Z   = fs->make<TH1D>("H1_eta_muon2_Z", "H1_eta_muon2_Z", 100, -6., 6.);
H1_wmu_pt_W      = fs->make<TH1D>("H1_wmu_pt_W", "H1_wmu_pt_W", 100, 0., 500.);
H1_wmu_eta_W     = fs->make<TH1D>("H1_wmu_eta_W", "H1_wmu_eta_W", 100, -6., 6.);
H1_wNu_pt        = fs->make<TH1D>("H1_wNu_pt   ", "H1_wNu_pt   ", 100, 0., 500.)    ;
H1_wNu_eta       = fs->make<TH1D>("H1_wNu_pt   ", "H1_wNu_pt   ", 100, -6., 6.)     ;
H1_wElec_pt      = fs->make<TH1D>("H1_wElec_pt ", "H1_wElec_pt ", 100, 0., 500.)    ;
H1_wElec_eta     = fs->make<TH1D>("H1_wElec_eta", "H1_wElec_eta", 100, -6., 6.)     ;
H1_wElec_Nu_pt   = fs->make<TH1D>("H1_wElec_Nu_pt ","H1_wElec_Nu_pt ", 100, 0., 500);
H1_wElec_Nu_eta  = fs->make<TH1D>("H1_wElec_Nu_eta","H1_wElec_Nu_eta", 100, -6., 6.);

H1_welec_pt     = fs->make<TH1D>(" H1_welec_pt ", "H1_welec_pt ", 100, 0., 500.)    ; 
H1_welec_eta    = fs->make<TH1D>(" H1_welec_eta", "H1_welec_eta", 100, -6, 6.)      ; 
                                                                                    
H1_wu_pt        = fs->make<TH1D>(" H1_wu_pt    ", "H1_wu_pt    ", 100, 0., 500.)    ; 
H1_wu_eta       = fs->make<TH1D>(" H1_wu_eta   ", "H1_wu_eta   ", 100, -6., 6.)     ; 
                                                                                    
H1_zuu_pt       = fs->make<TH1D>(" H1_zuu_pt   ", "H1_zuu_pt   ", 100, 0., 500.)    ; 
H1_zuu_eta      = fs->make<TH1D>(" H1_zuu_eta  ", "H1_zuu_eta  ", 100, -6., 6.)     ; 
                                                                                    
H1_zee_pt       = fs->make<TH1D>(" H1_zee_pt   ", "H1_zee_pt   ", 100, 0., 500.)    ; 
H1_zee_eta      = fs->make<TH1D>(" H1_zee_eta  ", "H1_zee_eta  ", 100, -6., 6.)     ; 

//---ALL---

H1_AllZmuons_pt  = fs->make<TH1D>("AllZmuons_pt","AllZmuons_pt", 100, 0., 500.)     ;
H1_AllZmuons_eta = fs->make<TH1D>("AllZmuons_eta","AllZmuons_eta",100, -6, 6)       ;
H1_AllZelec_pt   = fs->make<TH1D>("AllZelec_pt", "AllZelec_pt", 100, 0., 500.)      ;
H1_AllZelec_eta  = fs->make<TH1D>("AllZelec_eta","AllZelec_eta",100, -6, 6)         ;

    
}

// ------------ method called once each job just after ending the event loop  ------------
void
tbz_Acceptance::endJob()
{

}

// ------------ method called when starting to processes a run  ------------
void
tbz_Acceptance::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
tbz_Acceptance::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
tbz_Acceptance::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
tbz_Acceptance::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
tbz_Acceptance::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(tbz_Acceptance);
