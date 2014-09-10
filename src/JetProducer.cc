// -*- C++ -*-
//
// Package:    JetProducer
/**\class JetProducer JetProducer.cc JetProducer/JetProducer/src/JetProducer.cc

 Description: [one line class summary]

*/
// Original Author:  Muhammad Shoaib
//         Created:  Wed Jul 31 12:12:45 PKT 2013
// $Id$

// system include files
#include <memory>
#include <vector>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//-----------------------------------------------------
#include "DataFormats/JetReco/interface/Jet.h"
//#include "RecoJets/JetAnalyzers/interface/JetPlotsExample.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//for b-tagging 
#include "DataFormats/BTauReco/interface/JetTag.h"
//-----------------------------------------------------
#include "TH1.h"
#include "TH2.h"
#include "DataFormats/Common/interface/Handle.h"
#include "TStyle.h"
#include "math.h"
#include "TAxis.h"
#include "vector"
#include "TLegend.h"
#include "TFile.h" // ROOT, for saving file
#include "TVirtualPad.h" // ROOT, for interactive graphics
#include "TApplication.h" // ROOT , for interactive graphics
#include "TProfile.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/defs.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//-----------------------------------------------------
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
// class declaration
//
using namespace std;
using namespace edm;
using namespace reco;

class JetProducer : public edm::EDProducer {
   public:
      explicit JetProducer(const edm::ParameterSet&);
      ~JetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      double JetsPtCut_          ;
      double JetsEtaCut_         ;
    //-------------
    TH1D* H_JetsSize_Producer    ;
    TH1D* H_JetsPt               ;
      
};
// constructors and destructor
//
JetProducer::JetProducer(const edm::ParameterSet& iConfig)
{

    // produces vector ofJets Collection 
    
    // produces<std::vector<reco::CaloJet> >()                 ;
    produces<std::vector<reco::PFJet> >()                      ;
    produces<std::vector<reco::JetTag> >()                     ;     
    JetsPtCut_  = iConfig.getParameter<double>("JetsPtCut")    ;
    JetsEtaCut_ = iConfig.getParameter<double>("JetsEtaCut")   ;
    
}


JetProducer::~JetProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// member functions
// ------------ method called to produce the data  ------------
void
JetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     
     // Handle<CaloJetCollection > jets;
     edm::Handle<reco::PFJetCollection> jets                                                 ;
     // iEvent.getByLabel(string("ak5CaloJets"), jets)  ;ak5PFJets                            
     iEvent.getByLabel(string("ak5PFJets"), jets)                                             ;
     cout<<"PRINt JETS"<<endl                                                                 ;     
     // std::vector<reco::CaloJet> * newjets = new std::vector<reco::CaloJet>()               ;
     std::vector<reco::PFJet> * newjets = new std::vector<reco::PFJet>()                      ;      
     unsigned int ijet   = 0                                                                  ;
     // bool passed         = false                                                           ;                                                                                              
     double PFJetNHEF /*PFJetNHEF1*/,PFJetCHEF,PFJetNEMF,PFJetCEMF,NeutralHadIso,photonIso    ;     
     // for (CaloJetCollection::const_iterator JetsProd = jets->begin(); JetsProd !=jets->end();++JetsProd)  
     
     for (PFJetCollection::const_iterator JetsProd = jets->begin(); JetsProd !=jets->end();++JetsProd)
       {
       
         // if(JetsProd->pt() < 25.)                  continue     ;
         if(JetsProd->pt() < JetsPtCut_ )continue                  ;
         H_JetsPt->Fill(JetsProd->pt())                            ;
         cout<<"Pt_Jets_in_producer: "<<JetsProd->pt()<<endl       ;
         if( fabs(JetsProd->eta()) > JetsEtaCut_ )    continue     ;
         
         PFJetNHEF = JetsProd ->neutralHadronEnergyFraction()      ;
         PFJetCHEF = JetsProd ->chargedHadronEnergyFraction()      ;
         PFJetNEMF = JetsProd ->neutralEmEnergyFraction()          ;
         PFJetCEMF = JetsProd ->chargedEmEnergyFraction()          ;   
         
         cout << "neutralHadronEnergyFraction: = "<< PFJetNHEF  << ", chargedHadronEnergyFraction: = " << PFJetCHEF          ;
         cout << ", neutralEmEnergyFraction =: "    << PFJetNEMF <<  ", chargedEmEnergyFraction = "     << PFJetCEMF <<endl  ;
      
         if(PFJetNHEF > 0.90 )  continue                 ;
         if(PFJetCHEF <= 0. )   continue                 ;
         if(PFJetNEMF > 0.90 )  continue                 ;
         if(PFJetCEMF > 0.99 )  continue                 ;
         
         cout << "neutralHadronEnergyFraction1: = "<< PFJetNHEF  << ", chargedHadronEnergyFraction1: = " << PFJetCHEF           ;
         cout << ", neutralEmEnergyFraction1 =: "    << PFJetNEMF <<  ", chargedEmEnergyFraction1: = "     << PFJetCEMF <<endl  ;         
         
         ijet++                                          ; 
         H_JetsSize_Producer -> Fill(ijet)               ;
         newjets->push_back(* JetsProd)                  ;
         if( newjets->size()<=0)              continue   ;         
      }
      
      if(newjets) cout<< "size of my jet container "<< newjets->size()  <<endl        ;
      // std::auto_ptr<std::vector<reco::CaloJet> > ptr(newjets)                      ;
      std::auto_ptr<std::vector<reco::PFJet> > ptr(newjets)                           ;
      iEvent.put(ptr)                                                                 ;
      if(newjets)cout<<"jets registered in this event"<<endl                          ;      
      if(newjets==0) cout<<" WARNING: Selected jet container is empty!!!!!"<<endl     ;
      
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetProducer::beginJob()
{

edm::Service<TFileService> fs1                                                                 ;
H_JetsSize_Producer = fs1->make<TH1D> ("H_JetsSize_Producer","H_JetsSize_Producer",100, 0, 20) ;
H_JetsPt            = fs1->make<TH1D> ("H_JetsPt","H_JetsPt_Producer", 40, 0., 400.)           ;


}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
JetProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
JetProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetProducer);
