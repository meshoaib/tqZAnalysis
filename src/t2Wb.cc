// -*- C++ -*-
//
// Package:    t2Wb
// Class:      t2Wb
//
/**\class t2Wb t2Wb.cc t2Wb/t2Wb/src/t2Wb.cc

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
#include "CLHEP/Vector/EulerAngles.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "math.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TLorentzVector.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
//
// class declaration
//

class t2Wb : public edm::EDAnalyzer {
   public:
      explicit t2Wb(const edm::ParameterSet&);
     ~t2Wb();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

//GENERATOR
        TH1D * GEN_LEP_W;
        TH1D * GEN_LEP_TOP;
        TH1D * GEN_LEP_W_ETA;
        TH1D * GEN_LEP_W_PT;
        TH1D * GEN_LEP_TOP_ETA;
        TH1D * GEN_LEP_TOP_PT;
        TH1D * GEN_muon_pt;
        TH1D * GEN_electron_pt;
        TH1D * GEN_b_M;
//Eff   
        TH1D * HIST_GEN_MUON;
        TH1D * HIST_GEN_ELECTRON;
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
t2Wb::t2Wb(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs; // TFileService is a CMSSW Framework Service that 
                                  // allows one to create ROOT objects in multiple 
                                  // modules and store those histograms in the same
                                  // ROOT file.
   
   

//Generator Level
   GEN_LEP_W         =fs->make<TH1D>("GEN_LEP_W","GEN_LEP_W"            ,50,60.,120.);
   GEN_LEP_TOP       =fs->make<TH1D>("GEN_LEP_TOP","GEN_LEP_TOP"        ,200,140.,200.);
   GEN_LEP_W_ETA     =fs->make<TH1D>("GEN_LEP_W_ETA","GEN_LEP_W_ETA"    ,50,-5.,+5.);
   GEN_LEP_W_PT      =fs->make<TH1D>("GEN_LEP_W_PT","GEN_LEP_W_PT"      ,50,0.,400.);
   GEN_LEP_TOP_ETA   =fs->make<TH1D>("GEN_LEP_TOP_ETA","GEN_LEP_TOP_ETA",50,-5.,+5.);
   GEN_LEP_TOP_PT    =fs->make<TH1D>("GEN_LEP_TOP_PT","GEN_LEP_TOP_PT"  ,25,0.,450.);
   GEN_muon_pt       =fs->make<TH1D>("GEN_muon_pt","GEN_muon_PT"        ,50,0.,400.);
   GEN_electron_pt   =fs->make<TH1D>("GEN_electron_pt","GEN_electron_pt",50,0.,400.);
   GEN_b_M           =fs->make<TH1D>("GEN_b_M","GEN_b_mass"             ,50,0.,20.);
//Effeciency
   HIST_GEN_MUON     =fs->make<TH1D>("HIST_GEN_MUON","Gen_general_muons_PT"       ,10,0.,250.);
   HIST_GEN_ELECTRON =fs->make<TH1D>("HIST_GEN_ELECTRON","Gen_general_electron_pt",10,0.,250.);
}


t2Wb::~t2Wb()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
//
// member functions
//

// ------------ method called for each event  ------------
void
t2Wb::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

using namespace std;
using namespace edm;
using namespace pat;
using namespace reco;

//Gen particles loop
   Handle<GenParticleCollection> genParticle; // The default generator particle collection 
                                              // is reco::GenParticleCollection, which is 
                                              // a typedef for std::vector<reco::GenParticle>.
   iEvent.getByLabel("genParticles",genParticle);
   //std::vector<reco::genParticle>::const_iterator j;
   
   
   // TLorentzVector is a general four-vector class, 
   // which can be used either for the description 
   // of position and time (x,y,z,t) or momentum and 
   // energy (px,py,pz,E).   
   TLorentzVector    lep_Muon_p4      (0,0,0,0)
                    ,lep_Mu_plus_p4   (0,0,0,0)
                    ,lep_Mu_minus_p4  (0,0,0,0)
                    ,nu_Muon_p4       (0,0,0,0)
                    ,nu_Mu_plus_p4    (0,0,0,0)
                    ,nu_Mu_minus_p4   (0,0,0,0);
                    
   TLorentzVector    lepton_p4(0,0,0,0)
                    ,nutrino_p4(0,0,0,0);
                    
   TLorentzVector    lep_electron_p4(0,0,0,0)
                    ,lep_e_p4       (0,0,0,0)
                    ,lep_p_p4       (0,0,0,0)
                    ,nu_e_p4        (0,0,0,0)
                    ,nu_p_p4        (0,0,0,0)
                    ,nu_electron_p4 (0,0,0,0);
                    
   TLorentzVector    lep_p4 (0,0,0,0)
                   , nu_p4  (0,0,0,0)
                   , b_p4   (0,0,0,0)
                   , bbar_p4(0,0,0,0);
                   
   TLorentzVector    topplus (0,0,0,0)
                   , topminus(0,0,0,0)
                   , Wplus   (0,0,0,0)
                   , Wminus  (0,0,0,0);
                   
   TLorentzVector top_lep(0,0,0,0);
   
   TLorentzVector     w_lep(0,0,0,0);
   
   TLorentzVector  Wplus_lep   (0,0,0,0)
                 , Wminus_lep  (0,0,0,0)
                 , topPlus_lep (0,0,0,0)
                 , topminus_lep(0,0,0,0);
                 
   double  a, b, eta, eta1;
   //for(j=EvtHandle->begin();j!=EvtHandle->end();++j)
   for(size_t i = 0; i < genParticle->size(); ++ i)
   {
      const GenParticle & p = (*genParticle)[i]; // p is const GenParticle reference
                                                 // initialized to (*genParticle)[i],
                                                 // and (*genParticle)[i] is ith element
                                                 // of genParticle vector
      //const Candidate * p = (*j);
                                                     
      int id = p.pdgId();   // pgdId of ith element of genParticle vector
      
      //int st = p.status();
      //bool chka=false;        bool chkb=false;
      if (abs(p.pdgId())==13 && abs(p.mother()->pdgId())==24 && abs(p.mother()->mother()->pdgId())==6)
         {
         // daughter Muon (Mu+ or Mu-) is detected
         // its histogram is filled
         HIST_GEN_MUON->Fill(p.pt());
         }
      if (abs(p.pdgId())==11 && abs(p.mother()->pdgId())==24 && abs(p.mother()->mother()->pdgId())==6)
         {
         // daughter electron (e+ or e-) is detected
         // its histogram is filled
         HIST_GEN_ELECTRON->Fill(p.pt());
         }

      //         top
      if (abs(id)==6 )
      {
         unsigned int n = p.numberOfDaughters(); // reco::CompositePtrCandidate
                                                 // virtual size_t numberOfDaughters () const
                                                 // number of daughters
         for(size_t j=0;j<n;++j)
         {
            const Candidate * d = p.daughter( j ); // virtual Candidate * daughter (size_type)
                                                   // return daughter at a given position, 
                                                   // i = 0, ... numberOfDaughters() - 1
            int dauId=d->pdgId();
            //           (W+)         bottom (b)
            if(abs(dauId)!=24  && abs(dauId)!=5 )
               continue;

            // for b quark
            //          (b)
            if (dauId == 5)
            {
               b_p4.SetPx(d->px());
               b_p4.SetPy(d->py());
               b_p4.SetPz(d->pz());
               b_p4.SetE(d->energy());
            }
            //      anti_(b)
            if (dauId == -5)
            {
               bbar_p4.SetPx(d->px());
               bbar_p4.SetPy(d->py());
               bbar_p4.SetPz(d->pz());
               bbar_p4.SetE(d->energy());
            }
            
            bool chk1=false;  // will be true if W+ will be generated from Mu+
            bool chk2=false;  // will be true if W+ will be generated from  e+
            bool chk3=false;  // will be true if W- will be generated from Mu-
            bool chk4=false;  // will be true if W- will be generated from  e-

            // if it is W+ (pdgId()=24) then enter the following if
            if (dauId==24 )
            {
               unsigned int m=d->numberOfDaughters();
               for(size_t k=0; k<m;++k) // loop over number of daughters of W+
               {
                  const Candidate * w_dau = d->daughter(k);
                  int w_dauId=w_dau->pdgId();
                  //          W+             tau+   tau_neutrino
                  if(w_dauId==24 || w_dauId==-15 || w_dauId==16)
                     continue;
                  //cout<<"Selected W+ daughter Id is :"<<w_dauId<<endl;
                  //        (Mu+)
                  if(w_dauId==-13 )
                  {
                     chk1=true;
                     lep_Mu_plus_p4.SetPx(w_dau->px());
                     lep_Mu_plus_p4.SetPy(w_dau->py());
                     lep_Mu_plus_p4.SetPz(w_dau->pz());
                     lep_Mu_plus_p4.SetE (w_dau->energy());
                  }
                  //Mu_(neutrino)
                  if(w_dauId== 14 )
                  {
                     nu_Mu_minus_p4.SetPx(w_dau->px());
                     nu_Mu_minus_p4.SetPy(w_dau->py());
                     nu_Mu_minus_p4.SetPz(w_dau->pz());
                     nu_Mu_minus_p4.SetE (w_dau->energy());
                  }
                  //         (e+)
                  if(w_dauId==-11 )
                  {
                     chk2=true;
                     lep_p_p4.SetPx(w_dau->px());
                     lep_p_p4.SetPy(w_dau->py());
                     lep_p_p4.SetPz(w_dau->pz());
                     lep_p_p4.SetE (w_dau->energy());
                  }
                  //         (nu)
                  if(w_dauId==12  )
                  {
                     nu_e_p4.SetPx(w_dau->px());
                     nu_e_p4.SetPy(w_dau->py());
                     nu_e_p4.SetPz(w_dau->pz());
                     nu_e_p4.SetE (w_dau->energy());
                  }
               }//lepton+hadron loop*/
            }
            // if it is W- (pdgId()=-24) then enter the following if
            if (dauId==-24 )
            {
               unsigned int m=d->numberOfDaughters();
               for(size_t k=0; k<m;++k)  // loop over number of daughters of W-
               {
                  const Candidate * w_dau = d->daughter(k);
                  int w_dauId=w_dau->pdgId();
                  //          W-             tau    anti_tau_nutrino
                  if(w_dauId==-24 || w_dauId==15 ||  w_dauId==-16)
                     continue;

                  //        (Mu-)
                  if(w_dauId==13 )
                  {
                     chk3=true;
                     lep_Mu_minus_p4.SetPx(w_dau->px());
                     lep_Mu_minus_p4.SetPy(w_dau->py());
                     lep_Mu_minus_p4.SetPz(w_dau->pz());
                     lep_Mu_minus_p4.SetE (w_dau->energy());
                  }
                  // anti_Mu_nutrino
                  if(w_dauId==-14)
                  {
                     nu_Mu_plus_p4.SetPx(w_dau->px());
                     nu_Mu_plus_p4.SetPy(w_dau->py());
                     nu_Mu_plus_p4.SetPz(w_dau->pz());
                     nu_Mu_plus_p4.SetE (w_dau->energy());
                  }
                  //         (e-)
                  if(w_dauId==11  )
                  {
                     chk4=true;
                     lep_e_p4.SetPx(w_dau->px());
                     lep_e_p4.SetPy(w_dau->py());
                     lep_e_p4.SetPz(w_dau->pz());
                     lep_e_p4.SetE (w_dau->energy());
                  }
                  //     anti_(nu)
                  if(w_dauId==-12  )
                  {
                     nu_p_p4.SetPx(w_dau->px());
                     nu_p_p4.SetPy(w_dau->py());
                     nu_p_p4.SetPz(w_dau->pz());
                     nu_p_p4.SetE (w_dau->energy());
                  }
               }//lepton+hadron loop*/
            }//W- loop
            
            lep_Muon_p4    =lep_Mu_minus_p4+lep_Mu_plus_p4;
            lep_electron_p4=lep_e_p4+lep_p_p4;
            
            nu_Muon_p4    =nu_Mu_minus_p4+nu_Mu_plus_p4;
            nu_electron_p4=nu_e_p4+nu_p_p4;
            
            lepton_p4 =lep_Muon_p4+lep_electron_p4;
            nutrino_p4=nu_Muon_p4+nu_electron_p4;
            
            GEN_muon_pt    ->Fill(lep_Muon_p4.Pt());
            GEN_electron_pt->Fill(lep_electron_p4.Pt());
            
            w_lep=lepton_p4+nutrino_p4;
            
            if(chk1)      // Mu+
            {
               lep_Muon_p4=lep_Mu_plus_p4;
               w_lep=lep_Muon_p4+nutrino_p4;
               top_lep=w_lep+b_p4;
            }
            else if(chk2) //  e+
            {
               lep_electron_p4=lep_p_p4;
               w_lep=lep_electron_p4+nutrino_p4;
               top_lep=w_lep+bbar_p4;
            }
            else if(chk3) // Mu-
            {
               lep_Muon_p4=lep_Mu_minus_p4;
               w_lep=lep_Muon_p4+nutrino_p4;
               top_lep=w_lep+b_p4;
            }
            else if(chk4) //  e-
            {
               lep_electron_p4=lep_e_p4;
               w_lep=lep_electron_p4+nutrino_p4;
               top_lep=w_lep+bbar_p4;
            }
               
            //cout<<"top leptonic mass is"<<top_lep.M()<<endl;
            a=b_p4.M();
            GEN_b_M        ->Fill(a);
            a=top_lep.M();
            GEN_LEP_TOP    ->Fill(a);
            //eta=top_lep.Eta();
           // GEN_LEP_TOP_ETA->Fill(eta);
            b=top_lep.Pt();
            GEN_LEP_TOP_PT ->Fill(b);
            a=w_lep.M();
            GEN_LEP_W      ->Fill(a);
            //eta=w_lep.Eta();
            //GEN_LEP_W_ETA  ->Fill(eta);
            b=w_lep.Pt();
            GEN_LEP_W_PT   ->Fill(b);
            
         }//top daughter loop ends
      }//top if block ends
   }//genParticle loop ends
   
   
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
t2Wb::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
t2Wb::endJob()
{
 
}

// ------------ method called when starting to processes a run  ------------
void
t2Wb::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
t2Wb::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
t2Wb::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
t2Wb::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
t2Wb::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(t2Wb);
