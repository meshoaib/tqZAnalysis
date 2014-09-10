// -*- C++ -*-
//
// Package:    tbZ_Final
// Class:      tbZ_Final
//
/**\class tbZ_Final tbZ_Final.cc tbZ_Final/tbZ_Final/src/tbZ_Final.cc

 Description: [one line class summary]
 1 lepton and b-quark for w reconstruction (only lepton decay of W is to be studied)

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Muhammad Shoaib and Abid Farooq
//      Start date:  Sunday, 22/12/2013
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

class tbZ_Final : public edm::EDAnalyzer {
   public:
      explicit tbZ_Final(const edm::ParameterSet&);
     ~tbZ_Final();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // class TH1D: public TH1, public TArrayD
      // TH1D : histograms with one double per channel. Maximum precision 14 digits
      // ---------------------------------------------------------------------------
      // following are pointers of type TH1D
      //
      // for Z-boson
      TH1D* H1_DEL_Electron_Eta;
      TH1D* H1_DEL_Electron_Phi;
      TH1D* H1_Z_Mass_from_Electron_Inv_Mass;
      TH1D* H1_Z_Pt_from_Electron_Pt;
      TH2D* H2_Z_Electron_vs_Positron_Eta;
      TH2D* H2_Z_Electron_vs_Positron_Phi;
      TH2D* H2_Z_Electron_vs_Positron_Pt;

      TH1D* H1_DEL_Muon_Eta;
      TH1D* H1_DEL_Muon_Phi;
      TH1D* H1_Z_Mass_from_Muon_Inv_Mass;
      TH1D* H1_Z_Pt_from_Muon_Pt;
      TH2D* H2_Z_Muon_min_vs_Muon_pls_Eta;
      TH2D* H2_Z_Muon_min_vs_Muon_pls_Phi;
      TH2D* H2_Z_Muon_min_vs_Muon_pls_Pt;

      /* for top quark */

      TH1D* H1_for_t_b_Mass;

      TH1D* H1_for_t_W_from_e_Mass;
      TH1D* H1_for_t_W_from_Mu_Mass;
      TH1D* H1_for_t_W_from_e_Eta;
      TH1D* H1_for_t_W_from_Mu_Eta;
      TH1D* H1_for_t_W_from_e_Pt;
      TH1D* H1_for_t_W_from_Mu_Pt;

      TH1D* H1_t_from_W_e_Mass;
      TH1D* H1_t_from_W_Mu_Mass;
      TH1D* H1_t_from_W_e_Eta;
      TH1D* H1_t_from_W_Mu_Eta;
      TH1D* H1_t_from_W_e_Pt;
      TH1D* H1_t_from_W_Mu_Pt;
      
      
      /*  Deltas */
      
      TH1D* H1_b_pt; 
      TH1D* H1_b_eta;
      
      TH1D* H1_pt_for_W_from_e;   
      TH1D* H1_Eta_for_W_from_e;  
      TH1D* H1_del_pt_b_W_from_e; 
      TH1D* H1_del_eta_b_W_from_e;
      
      TH1D* H1_pt_for_W_from_Mu;
      TH1D* H1_Eta_for_W_from_Mu;  
      TH1D* H1_del_pt_b_W_from_Mu; 
      TH1D* H1_del_eta_b_W_from_Mu;
      
      
      
      //  free b qrk
      TH1D* H1_free_b_pt  ;
      TH1D* H1_free_b_eta ;
      TH1D* H1_free_b_Mass;
      // comparison
      
      TH1D* H1_del_t_Z_eee_Pt   ;
      TH1D* H1_del_t_Z_eeMu_Pt  ;
      TH1D* H1_del_t_Z_MuMue_Pt ;
      TH1D* H1_del_t_Z_MuMuMu_Pt;

      
/*       --------------------------------------------------------*/

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
tbZ_Final::tbZ_Final(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
}


tbZ_Final::~tbZ_Final()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
//
// member functions
//

// ------------ method called for each event  ------------
void
tbZ_Final::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

using namespace  std;
using namespace  edm;
using namespace  pat;
using namespace reco;

//Gen particles loop
   Handle<GenParticleCollection> genParticle; // The default generator particle collection
                                              // is reco::GenParticleCollection, which is
                                              // a typedef for std::vector<reco::GenParticle>.
   iEvent.getByLabel("genParticles",genParticle);
std::vector<int> ParticleAt;
std::vector<int> ParticleID;

bool isZ=false;
bool isb=false;
bool ist=false;
bool ifb=false;

   for(size_t i = 0; i < genParticle->size(); ++ i)
   {
      const GenParticle & p = (*genParticle)[i]; // p is const GenParticle reference
                                                 // initialized to (*genParticle)[i],
                                                 // and (*genParticle)[i] is ith element
                                                 // of genParticle vector

      int id = p.pdgId();   // pgdId of ith element of genParticle vector

          if ((abs(p.pdgId())==11 || abs(p.pdgId())==13
             || abs(p.pdgId())==12 || abs(p.pdgId())==14) && /* either e or Mu         */
               abs(p.mother()->pdgId())==24              && /* Mother must be W       */
               abs(p.mother()->pdgId())!=23              && /* Mother must not be Z   */
               (p.mother()->mother()->pdgId())==6)          /* Grand-Mother must be t */
         //if((p.mother()->mother()->pdgId())==6)
         {
         ParticleAt.push_back(i);
         ParticleID.push_back(abs(p.pdgId()));
         ist=true;
         }

         if((p.pdgId())==5 && (p.mother()->pdgId())==6) // b with mother t
         {
         ParticleAt.push_back(i);
         ParticleID.push_back(abs(p.pdgId()));
         isb=true;
         }

         if(abs(p.pdgId())==11 && abs(p.mother()->pdgId())==23) // electrons for Z
         {
         ParticleAt.push_back(i);
         ParticleID.push_back(abs(p.pdgId()));
         isZ=true;
         }

         if(abs(p.pdgId())==13 && abs(p.mother()->pdgId())==23) // Muons for Z
         {
         ParticleAt.push_back(i);
         ParticleID.push_back(abs(p.pdgId()));
         isZ=true;
         }
         if(p.pdgId()==5 && p.mother()->pdgId()!=6 && p.mother()->pdgId()!=2212/*&& p.mother()->mother()->pdgId()!=6*/) // free b
         {
         ParticleAt.push_back(i);
         ParticleID.push_back(abs(p.pdgId()));
         ifb=true;
         }
         
   }//genParticle loop ends for Z-boson
if(ist && isb && isZ && ifb)
{
cout<<"==========================================================================="<<endl;
   for(Size_t k=0; k<ParticleAt.size();k++)
   {
   const GenParticle & p = (*genParticle)[ParticleAt[k]];
   int id = p.pdgId();
   // if(abs(p.mother()->mother()->pdgId())==6)
   // {cout<<"Particle at : "<<setw(3)<<ParticleAt[k]<<" have pdgId() : "<<setw(3)<<id<<
        // "  has mother of pgdId() :"<<setw(3)<<p.mother()->pdgId()<<" Grand Mother : "<<setw(3)<<p.mother()->mother()->pdgId()<<endl;
   // }
   // else
   // {
   // cout<<"Particle at : "<<setw(3)<<ParticleAt[k]<<" have pdgId() : "<<setw(3)<<id<<
        // "  has mother of pgdId() :"<<setw(3)<<p.mother()->pdgId()<<endl;
   // }
   cout<<" p.pdgId() : "<<p.pdgId()<<"   Mother  "<<p.mother()->pdgId()<<endl; //" Grand Mother : "<<p.mother()->mother()->pdgId()<<endl;
   }
   
cout<<"==========================================================================="<<endl;
}


if(ist && isb && isZ && ifb)
{
cout<<"==========================================================================="<<endl;

double Z_Electron_pls_Charge=0.0
      ,Z_Electron_pls_Px=0.0
      ,Z_Electron_pls_Py=0.0
      ,Z_Electron_pls_Pz=0.0
      ,Z_Electron_pls_Energy=0.0
      ,Z_Electron_pls_Mass=0.0
      ,Z_Electron_pls_Pt=0.0
      ,Z_Electron_pls_Eta=0.0
      ,Z_Electron_pls_Phi=0.0;

double Z_Electron_min_Charge  =0.0
      ,Z_Electron_min_Px      =0.0
      ,Z_Electron_min_Py      =0.0
      ,Z_Electron_min_Pz      =0.0
      ,Z_Electron_min_Energy  =0.0
      ,Z_Electron_min_Mass    =0.0
      ,Z_Electron_min_Pt      =0.0
      ,Z_Electron_min_Eta     =0.0
      ,Z_Electron_min_Phi=0.0;

double Z_Muon_pls_Charge=0.0 
      ,Z_Muon_pls_Px    =0.0
      ,Z_Muon_pls_Py    =0.0
      ,Z_Muon_pls_Pz    =0.0
      ,Z_Muon_pls_Energy=0.0
      ,Z_Muon_pls_Mass  =0.0
      ,Z_Muon_pls_Pt    =0.0
      ,Z_Muon_pls_Eta   =0.0
      ,Z_Muon_pls_Phi=0.0;

double Z_Muon_min_Charge =0.0
      ,Z_Muon_min_Px     =0.0
      ,Z_Muon_min_Py     =0.0
      ,Z_Muon_min_Pz     =0.0
      ,Z_Muon_min_Energy =0.0
      ,Z_Muon_min_Mass   =0.0
      ,Z_Muon_min_Pt     =0.0
      ,Z_Muon_min_Eta    =0.0
      ,Z_Muon_min_Phi=0.0;
//SetPtEtaPhiE(pt,eta,phi,e)
TLorentzVector   b_qrk_p4       (0,0,0,0)
                ,Elec_pls_p4(0,0,0,0)
                ,Elec_nu_p4 (0,0,0,0)
                ,Muon_pls_p4(0,0,0,0)
                ,Muon_nu_p4 (0,0,0,0);
                
TLorentzVector   free_b_qrk_p4       (0,0,0,0)
                ,free_b_qrk_pt_eta_phi_e_p4(0,0,0,0);
                
TLorentzVector   b_qrk_pt_eta_phi_e_p4   (0,0,0,0)
                ,Elec_pls_pt_eta_phi_e_p4(0,0,0,0)
                ,Elec_nu_pt_eta_phi_e_p4 (0,0,0,0)
                ,Muon_pls_pt_eta_phi_e_p4(0,0,0,0)
                ,Muon_nu_pt_eta_phi_e_p4 (0,0,0,0);

TLorentzVector   W_boson_p4  (0,0,0,0)
                ,W_from_e_p4(0,0,0,0)
                ,W_from_Mu_p4 (0,0,0,0)
                ,t_from_e_p4 (0,0,0,0)
                ,t_from_Mu_p4 (0,0,0,0);
                
TLorentzVector     W_boson_pt_eta_phi_e_p4  (0,0,0,0)
                , W_from_e_pt_eta_phi_e_p4(0,0,0,0)
                ,W_from_Mu_pt_eta_phi_e_p4 (0,0,0,0)
                , t_from_e_pt_eta_phi_e_p4 (0,0,0,0)
                ,t_from_Mu_pt_eta_phi_e_p4 (0,0,0,0);
                
   bool is_W_e=false;
   bool is_W_Mu=false;
   bool is_W_e_nu=false;
   bool is_W_Mu_nu=false;
   bool is_b=false;
   bool is_fb=false;
   
   double    b_pt=0.0
            ,W_from_e_pt=0.0
            ,W_from_Mu_pt=0.0
            ,b_eta=0.0
            ,W_from_e_eta=0.0
            ,W_from_Mu_eta=0.0;
   double    del_pt_b_W_from_e=0.0
            ,del_eta_b_W_from_e=0.0
            ,del_pt_b_W_from_Mu=0.0
            ,del_eta_b_W_from_Mu=0.0;
            
   double   free_b_pt =0.0
           ,free_b_eta=0.0;
double t_from_W_Mu_Pt=0.0;
double t_from_W_Mu_Eta=0.0;            
double t_from_W_e_Pt=0.0   ;         
double t_from_W_e_Eta=0.0   ;         
            
 bool is_z_ee=false;
 bool is_z_MuMu=false;
   for(Size_t k=0; k<ParticleAt.size();k++)
   {
   const GenParticle & p = (*genParticle)[ParticleAt[k]];
   int id = p.pdgId();

   /*     ================================================
                  Z-Boson re-generation loop
          -------------------------------------------------         */

   if((p.mother()->pdgId())==23) // Daughters of Z loop
   {
      if(abs(p.pdgId())== 11) // Electrons of Z
      {
      is_z_ee=true;
         if((p.pdgId())==11) // e- of Z
         {
            Z_Electron_min_Eta    = (p.eta());
            Z_Electron_min_Phi    = (p.phi());
            Z_Electron_min_Energy = (p.energy());
            Z_Electron_min_Px     = (p.px());
            Z_Electron_min_Py     = (p.py());
            Z_Electron_min_Pz     = (p.pz());
            Z_Electron_min_Pt     = (p.pt());
         } // e- of Z ends
         if((p.pdgId())==-11) // e+ of Z
         {
            Z_Electron_pls_Eta    = (p.eta());
            Z_Electron_pls_Phi    = (p.phi());
            Z_Electron_pls_Energy = (p.energy());
            Z_Electron_pls_Px     = (p.px());
            Z_Electron_pls_Py     = (p.py());
            Z_Electron_pls_Pz     = (p.pz());
            Z_Electron_pls_Pt     = (p.pt());
         } // e+ of Z ends
      } // Electrons of Z ends

      if(abs(p.pdgId())==13) // Muons of Z
      {
      is_z_MuMu=true;
         if((p.pdgId())==13) // Mu- of Z
         {
            Z_Muon_min_Eta    = (p.eta());
            Z_Muon_min_Phi    = (p.phi());
            Z_Muon_min_Energy = (p.energy());
            Z_Muon_min_Px     = (p.px());
            Z_Muon_min_Py     = (p.py());
            Z_Muon_min_Pz     = (p.pz());
            Z_Muon_min_Pt     = (p.pt());
         } // Mu- of Z ends
         if((p.pdgId())==-13) // Mu+ of Z
         {
            Z_Muon_pls_Eta    = (p.eta());
            Z_Muon_pls_Phi    = (p.phi());
            Z_Muon_pls_Energy = (p.energy());
            Z_Muon_pls_Px     = (p.px());
            Z_Muon_pls_Py     = (p.py());
            Z_Muon_pls_Pz     = (p.pz());
            Z_Muon_pls_Pt     = (p.pt());
         } // Mu+ of Z ends
      } // Muons of Z ends

   } // Daughters of Z loop ends

   /*     ================================================
                           b-quark loop
          -------------------------------------------------         */
   if((p.pdgId())==5 && (p.mother()->pdgId())==6 ) // b quark having mother t
   { is_b=true;
      b_qrk_p4.SetPx(p.px());
      b_qrk_p4.SetPy(p.py());
      b_qrk_p4.SetPz(p.pz());
      b_qrk_p4.SetE (p.energy());
      b_qrk_pt_eta_phi_e_p4.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
      b_pt=b_qrk_p4.Pt();
      b_eta=b_qrk_p4.Eta();
      H1_b_pt  -> Fill(b_pt);
      H1_b_eta -> Fill(b_eta);
      H1_for_t_b_Mass         ->Fill(b_qrk_p4.M());
   } // b quark if check ends
   
   /*     ================================================
                        free b-quark loop
          -------------------------------------------------         */
   if(p.pdgId()==5 && p.mother()->pdgId()!=6 && p.mother()->pdgId()!=2212 ) // free b quark having mother t
   { is_fb=true;
      free_b_qrk_p4.SetPx(p.px());
      free_b_qrk_p4.SetPy(p.py());
      free_b_qrk_p4.SetPz(p.pz());
      free_b_qrk_p4.SetE (p.energy());
      free_b_qrk_pt_eta_phi_e_p4.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
      free_b_pt      = free_b_qrk_p4.Pt();
      free_b_eta     =free_b_qrk_p4.Eta();
      H1_free_b_pt   -> Fill(free_b_pt);
      H1_free_b_eta  -> Fill(free_b_eta);
      H1_free_b_Mass ->Fill(free_b_qrk_p4.M());
   } // free b quark if check ends




   /*     ============================================================================
           Leptons ( e+ and Mu+ ) and neutrino (e_nu and Mu_nu) to top (t) quark loop
          ----------------------------------------------------------------------------         */

   if((p.mother()->mother()->pdgId())==6 && p.pdgId()!=5 )  // Grand Daughters of top
   {
   if(p.pdgId()==-11) // e+ grand daughter of t
   {
      is_W_e=true;
      Elec_pls_p4.SetPx(p.px());
      Elec_pls_p4.SetPy(p.py());
      Elec_pls_p4.SetPz(p.pz());
      Elec_pls_p4.SetE (p.energy());
      Elec_pls_pt_eta_phi_e_p4.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
   } // e+ grand daughter of t ends
   if(p.pdgId()==12) // e_nu grand daughter of t
   {
      is_W_e_nu=true;
      Elec_nu_p4.SetPx(p.px());
      Elec_nu_p4.SetPy(p.py());
      Elec_nu_p4.SetPz(p.pz());
      Elec_nu_p4.SetE (p.energy());
      Elec_nu_pt_eta_phi_e_p4.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
   } // e_nu grand daughter of t ends
   if(p.pdgId()==-13) // Mu+ grand daughter of t
   {
      is_W_Mu=true;
      Muon_pls_p4.SetPx(p.px());
      Muon_pls_p4.SetPy(p.py());
      Muon_pls_p4.SetPz(p.pz());
      Muon_pls_p4.SetE (p.energy());
      Muon_pls_pt_eta_phi_e_p4.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
   } // Mu+ grand daughter of t ends
   if(p.pdgId()==14) // Mu_nu grand daughter of t
   {
      is_W_Mu_nu=true;
      Muon_nu_p4.SetPx(p.px());
      Muon_nu_p4.SetPy(p.py());
      Muon_nu_p4.SetPz(p.pz());
      Muon_nu_p4.SetE (p.energy());
      Muon_nu_pt_eta_phi_e_p4.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
   } // Mu_nu grand daughter of t ends

   } // Grand Daughters of top ends

   if(is_W_e && is_W_e_nu && is_b) // top from W+ from e+
   {
   W_from_e_p4 = Elec_pls_p4 + Elec_nu_p4;
   t_from_e_p4 = W_from_e_p4 + b_qrk_p4; 
   
   W_from_e_pt=W_from_e_p4.Pt();
   W_from_e_eta=W_from_e_p4.Eta();
   
   del_pt_b_W_from_e=b_pt-W_from_e_pt;
   del_eta_b_W_from_e=b_eta-W_from_e_eta;
   
   H1_pt_for_W_from_e    ->Fill(W_from_e_pt);
   H1_Eta_for_W_from_e   ->Fill(W_from_e_eta);
   H1_del_pt_b_W_from_e  ->Fill(del_pt_b_W_from_e);
   H1_del_eta_b_W_from_e ->Fill(del_eta_b_W_from_e);
   
   W_from_e_pt_eta_phi_e_p4.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
   t_from_e_pt_eta_phi_e_p4.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());

   H1_for_t_W_from_e_Mass  ->Fill(W_from_e_p4.M());
   H1_for_t_W_from_e_Pt    ->Fill(W_from_e_p4.Pt());
   H1_for_t_W_from_e_Eta   ->Fill(W_from_e_p4.Eta());
   H1_t_from_W_e_Mass      ->Fill(t_from_e_p4.M());
   H1_t_from_W_e_Eta       ->Fill(t_from_e_p4.Eta());
   H1_t_from_W_e_Pt        ->Fill(t_from_e_p4.Pt());
   
   t_from_W_e_Pt=t_from_e_p4.Pt();
   t_from_W_e_Eta=t_from_e_p4.Eta();
   
   }
   if(is_W_Mu && is_W_Mu_nu && is_b) // top from W+ from Mu+
   {
   W_from_Mu_p4= Muon_pls_p4 + Muon_nu_p4;
   t_from_Mu_p4= W_from_Mu_p4+ b_qrk_p4;
   
   W_from_Mu_pt=W_from_Mu_p4.Pt();
   W_from_Mu_eta=W_from_Mu_p4.Eta();
   
   del_pt_b_W_from_Mu=b_pt-W_from_Mu_pt;
   del_eta_b_W_from_Mu=b_eta-W_from_Mu_eta;
   
   H1_pt_for_W_from_Mu    ->Fill(W_from_Mu_pt);
   H1_Eta_for_W_from_Mu   ->Fill(W_from_Mu_eta);
   H1_del_pt_b_W_from_Mu  ->Fill(del_pt_b_W_from_Mu);
   H1_del_eta_b_W_from_Mu ->Fill(del_eta_b_W_from_Mu);
   
   W_from_Mu_pt_eta_phi_e_p4.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
   t_from_Mu_pt_eta_phi_e_p4.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
   
   H1_for_t_W_from_Mu_Mass ->Fill(W_from_Mu_p4.M());
   H1_for_t_W_from_Mu_Eta  ->Fill(W_from_Mu_p4.Eta());
   H1_for_t_W_from_Mu_Pt   ->Fill(W_from_Mu_p4.Pt());
   H1_t_from_W_Mu_Mass     ->Fill(t_from_Mu_p4.M());
   H1_t_from_W_Mu_Eta      ->Fill(t_from_Mu_p4.Eta());
   H1_t_from_W_Mu_Pt       ->Fill(t_from_Mu_p4.Pt());
   
   t_from_W_Mu_Pt=t_from_Mu_p4.Pt();
   t_from_W_Mu_Eta=t_from_Mu_p4.Eta();
   
   }

   
} // Loop over event particles
 


   // ==================================================================
 //            Z generation
              double Z_Pt_from_Electron_Pt;
              double Z_Pt_from_Muon_Pt;
 //          from electron and positron
        if(is_z_ee)
        {
        double DEL_Electron_Eta    = (Z_Electron_min_Eta - Z_Electron_pls_Eta);
        double DEL_Electron_Phi    = (Z_Electron_min_Phi - Z_Electron_pls_Phi);
        double Z_Mass_from_Electron_Inv_Mass_SQUARE  = (  ((Z_Electron_min_Energy+Z_Electron_pls_Energy)
                               *(Z_Electron_min_Energy+Z_Electron_pls_Energy))
                            -( ((Z_Electron_min_Px+Z_Electron_pls_Px)
                               *(Z_Electron_min_Px+Z_Electron_pls_Px))
                              +((Z_Electron_min_Py+Z_Electron_pls_Py)
                               *(Z_Electron_min_Py+Z_Electron_pls_Py))
                              +((Z_Electron_min_Pz+Z_Electron_pls_Pz)
                               *(Z_Electron_min_Pz+Z_Electron_pls_Pz)) )  ); // e^2-(Px^2 + Py^2 + Pz^2)
        double Z_Mass_from_Electron_Inv_Mass   = sqrt(Z_Mass_from_Electron_Inv_Mass_SQUARE);               // sqrt( e^2-(Px^2 + Py^2 + Pz^2) )
        double Z_Pt_from_Electron_Pt_SQUARE    = ( ((Z_Electron_min_Px+Z_Electron_pls_Px)
                              *(Z_Electron_min_Px+Z_Electron_pls_Px))
                             +((Z_Electron_min_Py+Z_Electron_pls_Py)
                              *(Z_Electron_min_Py+Z_Electron_pls_Py)) );     // Px^2 + Py^2
        Z_Pt_from_Electron_Pt       = sqrt(Z_Pt_from_Electron_Pt_SQUARE);                 // sqrt( Px^2 + Py^2)
        //charge_mult_elec  -> Fill(nChg);
        H1_DEL_Electron_Eta                -> Fill(DEL_Electron_Eta);
        H1_DEL_Electron_Phi                -> Fill(DEL_Electron_Phi);
        H1_Z_Mass_from_Electron_Inv_Mass   -> Fill(Z_Mass_from_Electron_Inv_Mass);
        H1_Z_Pt_from_Electron_Pt           -> Fill(Z_Pt_from_Electron_Pt);
        H2_Z_Electron_vs_Positron_Eta      -> Fill(Z_Electron_min_Eta,Z_Electron_pls_Eta);
        H2_Z_Electron_vs_Positron_Phi      -> Fill(Z_Electron_min_Phi,Z_Electron_pls_Phi);
        H2_Z_Electron_vs_Positron_Pt       -> Fill(Z_Electron_min_Pt,Z_Electron_pls_Pt);
        }

  //          from Mu- and Mu+
        if(is_z_MuMu)
        {
        double DEL_Muon_Eta    = (Z_Muon_min_Eta - Z_Muon_pls_Eta);
        double DEL_Muon_Phi    = (Z_Muon_min_Phi - Z_Muon_pls_Phi);
        double Z_Mass_from_Muon_Inv_Mass_SQUARE  = (  ((Z_Muon_min_Energy+Z_Muon_pls_Energy)
                               *(Z_Muon_min_Energy+Z_Muon_pls_Energy))
                            -( ((Z_Muon_min_Px+Z_Muon_pls_Px)
                               *(Z_Muon_min_Px+Z_Muon_pls_Px))
                              +((Z_Muon_min_Py+Z_Muon_pls_Py)
                               *(Z_Muon_min_Py+Z_Muon_pls_Py))
                              +((Z_Muon_min_Pz+Z_Muon_pls_Pz)
                               *(Z_Muon_min_Pz+Z_Muon_pls_Pz)) )  ); // e^2-(Px^2 + Py^2 + Pz^2)
        double Z_Mass_from_Muon_Inv_Mass   = sqrt(Z_Mass_from_Muon_Inv_Mass_SQUARE);               // sqrt( e^2-(Px^2 + Py^2 + Pz^2) )
        double Z_Pt_from_Muon_Pt_SQUARE    = ( ((Z_Muon_min_Px+Z_Muon_pls_Px)
                              *(Z_Muon_min_Px+Z_Muon_pls_Px))
                             +((Z_Muon_min_Py+Z_Muon_pls_Py)
                              *(Z_Muon_min_Py+Z_Muon_pls_Py)) );     // Px^2 + Py^2
        Z_Pt_from_Muon_Pt       = sqrt(Z_Pt_from_Muon_Pt_SQUARE);                 // sqrt( Px^2 + Py^2)
        //charge_mult_elec  -> Fill(nChg);
        H1_DEL_Muon_Eta                 -> Fill(DEL_Muon_Eta);
        H1_DEL_Muon_Phi                 -> Fill(DEL_Muon_Phi);
        H1_Z_Mass_from_Muon_Inv_Mass    -> Fill(Z_Mass_from_Muon_Inv_Mass);
        H1_Z_Pt_from_Muon_Pt            -> Fill(Z_Pt_from_Muon_Pt);
        H2_Z_Muon_min_vs_Muon_pls_Eta   -> Fill(Z_Muon_min_Eta,Z_Muon_pls_Eta);
        H2_Z_Muon_min_vs_Muon_pls_Phi   -> Fill(Z_Muon_min_Phi,Z_Muon_pls_Phi);
        H2_Z_Muon_min_vs_Muon_pls_Pt    -> Fill(Z_Muon_min_Pt,Z_Muon_pls_Pt);
        }
        
        /*   tbZ  properties   */
        if(ist && isb)
        {
        double del_t_Z_eee_Pt   =0.0 
              ,del_t_Z_eeMu_Pt  =0.0
              ,del_t_Z_MuMue_Pt =0.0
              ,del_t_Z_MuMuMu_Pt=0.0;
              
        double del_t_Z_eee_Eta   =0.0
              ,del_t_Z_eeMu_Eta  =0.0
              ,del_t_Z_MuMue_Eta =0.0
              ,del_t_Z_MuMuMu_Eta=0.0;
              
              del_t_Z_eee_Pt   = Z_Pt_from_Electron_Pt-t_from_W_e_Pt;
              del_t_Z_eeMu_Pt  = Z_Pt_from_Electron_Pt-t_from_W_Mu_Pt;
              del_t_Z_MuMue_Pt = Z_Pt_from_Muon_Pt-t_from_W_e_Pt;
              del_t_Z_MuMuMu_Pt= Z_Pt_from_Muon_Pt-t_from_W_Mu_Pt;
              
              H1_del_t_Z_eee_Pt   -> Fill(del_t_Z_eee_Pt   );
              H1_del_t_Z_eeMu_Pt  -> Fill(del_t_Z_eeMu_Pt  );
              H1_del_t_Z_MuMue_Pt -> Fill(del_t_Z_MuMue_Pt );
              H1_del_t_Z_MuMuMu_Pt-> Fill(del_t_Z_MuMuMu_Pt);
              
              
              // del_t_Z_eee_Eta   = Z_Pt_from_Electron_Pt-t_from_W_e_Pt;
              // del_t_Z_eeMu_Eta  = Z_Pt_from_Electron_Pt-t_from_W_Mu_Pt;
              // del_t_Z_MuMue_Eta = Z_Pt_from_Muon_Pt-t_from_W_e_Pt;
              // del_t_Z_MuMuMu_Eta= Z_Pt_from_Muon_Pt-t_from_W_Mu_Pt;
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
}


// ------------ method called once each job just before starting event loop  ------------
void
tbZ_Final::beginJob()
{
edm::Service<TFileService> fs;

    H1_DEL_Electron_Eta              = fs->make<TH1D> ("Delta_Eta of e+ and e-",  "pseudorapidity difference plot for electrons",   100,-2.5,2.5);
    H1_DEL_Electron_Phi              = fs->make<TH1D> ("Delta_Phi of e+ and e-",  "delta phi plot for electrons",                   100,-3.14,3.14);
    H1_Z_Mass_from_Electron_Inv_Mass = fs->make<TH1D> ("Mass of Z from e+ and e-","invariant mass of Z boson from Electrons plot",                100,60.,120.);
    H1_Z_Pt_from_Electron_Pt         = fs->make<TH1D> ("Pt   of Z from e+ and e-","pT Z plot from Electrons",                                           100,0.,100.);
    H2_Z_Electron_vs_Positron_Eta    = fs->make<TH2D> ("e+ Eta vs e- Eta",        "2D plot for e- and e+",                       100,-2.5,2.5,100,-2.5,2.5);
    H2_Z_Electron_vs_Positron_Phi    = fs->make<TH2D> ("e+ Phi vs e- Phi",        "2D plot for e- and e+",                       100,-3.14,3.14,100,-3.14,3.14);
    H2_Z_Electron_vs_Positron_Pt     = fs->make<TH2D> ("e+ Pt vs e- Pt",          "2D plot for e- and e+",                          100,-2.,100.,100,-2.,100.);

    H1_DEL_Muon_Eta                  = fs->make<TH1D> ("Delta_Eta of Mu+ and Mu-",  "pseudorapidity difference plot for Muons",   100,-2.5,2.5);
    H1_DEL_Muon_Phi                  = fs->make<TH1D> ("Delta_Phi of Mu+ and Mu-",  "delta phi plot for Muon",                   100,-3.14,3.14);
    H1_Z_Mass_from_Muon_Inv_Mass     = fs->make<TH1D> ("Mass of Z from Mu+ and Mu-","invariant mass of Z boson from Muons plot",                100,60.,120.);
    H1_Z_Pt_from_Muon_Pt             = fs->make<TH1D> ("Pt   of Z from Mu+ and Mu-","pT Z plot from Muons ",                                           100,0.,100.);
    H2_Z_Muon_min_vs_Muon_pls_Eta    = fs->make<TH2D> ("Mu+ Eta vs Mu- Eta",        "2D plot for Mu- Eta and Mu+ Eta",                       100,-2.5,2.5,100,-2.5,2.5);
    H2_Z_Muon_min_vs_Muon_pls_Phi    = fs->make<TH2D> ("Mu+ Phi vs Mu- Phi",        "2D plot for Mu- Eta and Mu+ Eta",                       100,-3.14,3.14,100,-3.14,3.14);
    H2_Z_Muon_min_vs_Muon_pls_Pt     = fs->make<TH2D> ("Mu+ Pt vs Mu- Pt",          "2D plot for Mu- Eta and Mu+ Eta",                          100,-2.,100.,100,-2.,100.);

    /*  for top quark */
    H1_for_t_b_Mass                 =fs->make<TH1D>("H1_for_t_b_Mass","H1_for_t_b_Mass"                 ,100,0.,20.);
    
    H1_for_t_W_from_e_Mass          =fs->make<TH1D>("H1_for_t_W_from_e_Mass ","H1_for_t_W_from_e_Mass " ,100,0.,200.);
    H1_for_t_W_from_Mu_Mass         =fs->make<TH1D>("H1_for_t_W_from_Mu_Mass","H1_for_t_W_from_Mu_Mass" ,100,0.,200.);
    H1_for_t_W_from_e_Eta           =fs->make<TH1D>("H1_for_t_W_from_e_Eta  ","H1_for_t_W_from_e_Eta  " ,200,-10.,10.);
    H1_for_t_W_from_Mu_Eta          =fs->make<TH1D>("H1_for_t_W_from_Mu_Eta ","H1_for_t_W_from_Mu_Eta " ,200,-10.,10.);
    H1_for_t_W_from_e_Pt            =fs->make<TH1D>("H1_for_t_W_from_e_Pt   ","H1_for_t_W_from_e_Pt   " ,100,0.,100.);
    H1_for_t_W_from_Mu_Pt           =fs->make<TH1D>("H1_for_t_W_from_Mu_Pt  ","H1_for_t_W_from_Mu_Pt  " ,100,0.,100.);
                                  
    H1_t_from_W_e_Mass              =fs->make<TH1D>("H1_t_from_W_e_Mass     ","H1_t_from_W_e_Mass     " ,100,100.,250.);
    H1_t_from_W_Mu_Mass             =fs->make<TH1D>("H1_t_from_W_Mu_Mass    ","H1_t_from_W_Mu_Mass    " ,100,100.,250.);
    H1_t_from_W_e_Eta               =fs->make<TH1D>("H1_t_from_W_e_Eta      ","H1_t_from_W_e_Eta      " ,100,-10.,10.);
    H1_t_from_W_Mu_Eta              =fs->make<TH1D>("H1_t_from_W_Mu_Eta     ","H1_t_from_W_Mu_Eta     " ,105,-10.,10.);
    H1_t_from_W_e_Pt                =fs->make<TH1D>("H1_t_from_W_e_Pt       ","H1_t_from_W_e_Pt       " ,100,0.,100.);
    H1_t_from_W_Mu_Pt               =fs->make<TH1D>("H1_t_from_W_Mu_Pt      ","H1_t_from_W_Mu_Pt      " ,100,0.,100.);
    
    
    H1_b_pt                        =fs->make<TH1D>("H1_b_pt               ","H1_b_pt               ",100,0.,100.);
    H1_b_eta                       =fs->make<TH1D>("H1_b_eta              ","H1_b_eta              ",100,-10.,10.);
                                                                                                     
    H1_pt_for_W_from_e             =fs->make<TH1D>("H1_pt_for_W_from_e    ","H1_pt_for_W_from_e    ",100,0.,100.);
    H1_Eta_for_W_from_e            =fs->make<TH1D>("H1_Eta_for_W_from_e   ","H1_Eta_for_W_from_e   ",100,-10.,10.);
    H1_del_pt_b_W_from_e           =fs->make<TH1D>("H1_del_pt_b_W_from_e  ","H1_del_pt_b_W_from_e  ",100,-10.,10.);
    H1_del_eta_b_W_from_e          =fs->make<TH1D>("H1_del_eta_b_W_from_e ","H1_del_eta_b_W_from_e ",100,-10.,10.);
                                                                                                     
    H1_pt_for_W_from_Mu            =fs->make<TH1D>("H1_pt_for_W_from_Mu   ","H1_pt_for_W_from_Mu   ",100,0.,100.);
    H1_Eta_for_W_from_Mu           =fs->make<TH1D>("H1_Eta_for_W_from_Mu  ","H1_Eta_for_W_from_Mu  ",100,-10.,10.);
    H1_del_pt_b_W_from_Mu          =fs->make<TH1D>("H1_del_pt_b_W_from_Mu ","H1_del_pt_b_W_from_Mu ",100,-10.,10.);
    H1_del_eta_b_W_from_Mu         =fs->make<TH1D>("H1_del_eta_b_W_from_Mu","H1_del_eta_b_W_from_Mu",100,-10.,10.);
    
    // comparison
    
    H1_del_t_Z_eee_Pt             =fs->make<TH1D>("H1_del_t_Z_eee_Pt   ","H1_del_t_Z_eee_Pt   ",100,-10.0,10.0);
    H1_del_t_Z_eeMu_Pt            =fs->make<TH1D>("H1_del_t_Z_eeMu_Pt  ","H1_del_t_Z_eeMu_Pt  ",100,-10.0,10.0);
    H1_del_t_Z_MuMue_Pt           =fs->make<TH1D>("H1_del_t_Z_MuMue_Pt ","H1_del_t_Z_MuMue_Pt ",100,-10.0,10.0);
    H1_del_t_Z_MuMuMu_Pt          =fs->make<TH1D>("H1_del_t_Z_MuMuMu_Pt","H1_del_t_Z_MuMuMu_Pt",100,-10.0,10.0);
    
    // free b qrk
    H1_free_b_pt                 =fs->make<TH1D>("H1_free_b_pt   "," H1_free_b_pt  ",100,0.,100.); 
    H1_free_b_eta                =fs->make<TH1D>("H1_free_b_eta  "," H1_free_b_eta ",100,-10.,10.); 
    H1_free_b_Mass               =fs->make<TH1D>("H1_free_b_Mass "," H1_free_b_Mass",100,0.,20.); 
                                 
    
                                   

}

// ------------ method called once each job just after ending the event loop  ------------
void
tbZ_Final::endJob()
{

}

// ------------ method called when starting to processes a run  ------------
void
tbZ_Final::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
tbZ_Final::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
tbZ_Final::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
tbZ_Final::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
tbZ_Final::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(tbZ_Final);
