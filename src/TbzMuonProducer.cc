#include "MyAnalysis/TbZ/interface/tbzMuonProducer.h"

 using namespace edm  ;
 using namespace std  ;
 using namespace reco ;

class PtGreater_Muon {
public:
        template <typename T> bool operator () (const T& i, const T& j) {
        return (i.pt() > j.pt());
        }
};
// constructors and destructor
TbzMuonProducer::TbzMuonProducer(const edm::ParameterSet& iConfig)
{
       //now do what ever other initialization is needed
      // produces vector of muons
      // produces<std::vector<reco::Muon> >()                                  ; //original
      // reco::PFCandidate
      produces<std::vector<reco::Muon> >()                                  ;
      // produces<std::vector<reco::PFCandidateCollection> >()              ;
      muonPtCut_         = iConfig.getParameter<double>("muonPtCut")        ;
      muonEtaCut_        = iConfig.getParameter<double>("muonEtaCut")       ;
      muonTag_           = iConfig.getParameter<edm::InputTag>("muonTag")   ;
   
}
TbzMuonProducer::~TbzMuonProducer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}
// member functions
// ------------ method called to produce the data  ------------
void
TbzMuonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
      bool result =false                ;
      std::vector<reco::Muon> * wmuons = new std::vector<reco::Muon>()      ;
      std::auto_ptr<std::vector<Muon> > ptr(wmuons)                         ; 
      
      // std::vector<reco::PFCandidate> * wmuons = new std::vector<reco::Muon>() ;//orginal
      // std::auto_ptr<std::vector<Muon> > ptr(wmuons)                           ;//original
      
      Handle<View<MET> > metCollection                                      ;
      iEvent.getByLabel(edm::InputTag("htMetAK5"), metCollection)           ;
      edm::Ptr<reco::MET> met( metCollection,0)                             ;
      ///all events                                                         
      m_muonCutFlow->Fill(0)                                                ; 

      //============== vertex information ==============================    
      edm::Handle< std::vector<reco::Vertex> > vertices                   ;
      iEvent.getByLabel("offlinePrimaryVertices", vertices)               ;
     // pick the first (i.e. highest sum pt) vertex                     
     if (!vertices.isValid())
      {
       std::cout<<"Didja hear the one about the empty vertex collection?\n"   ;
       iEvent.put(ptr)                                                        ;
       return                                                                 ;
      }
      // require in the event that there is at least one reconstructed vertex
      if(vertices->size()<=0) {
      iEvent.put(ptr)  ;
      return           ;
      }
      
      // pick the first (i.e. highest sum pt) vertex
      const reco::Vertex* theVertex=&(vertices->front());
      
      // require that the vertex meets certain criteria
      if(theVertex->ndof()<5   ||  fabs(theVertex->z())>24.0   
                               ||  fabs(theVertex->position().rho())>2.0 )
       {
         iEvent.put(ptr)        ;
         return                 ;      
      }
     ////after good vertex                                              
       m_muonCutFlow->Fill(1)   ;
      //================================================================
       
       //trigger 
       if (m_triggerSelector and m_triggerCache.setEvent(iEvent,iSetup))
       {
          // if the L1 or HLT configurations have changed, (re)initialize the filters (including during the first event)
          if (m_triggerCache.configurationUpdated())
             {
                m_triggerSelector ->init(m_triggerCache)      ;
                result = (*m_triggerSelector)(m_triggerCache) ;
             }
       }
       //===============================================================       
       if (result ==1)
       {      
          std::vector<reco::Vertex>::const_iterator itv  ;
          int NVtx = 0                                   ;
          // now, count vertices                         
          for (itv = vertices->begin(); itv != vertices->end(); ++itv)
             {
              // require that the vertex meets certain criteria
              if(itv->ndof()<5) continue                       ;
              if(fabs(itv->z())>50.0) continue                 ;
              if(fabs(itv->position().rho())>2.0) continue     ;
              ++NVtx                                           ;
             }
          NpVertices -> Fill(float(NVtx),1)                    ;
       }//if(==1)

      if(!result) {
      
                  }
        m_muonCutFlow->Fill(2)                        ; 
        
       edm::Handle<reco::MuonCollection> Wleptons     ;
       iEvent.getByLabel(muonTag_,Wleptons)           ; 
       reco::MuonCollection myMuon_new( Wleptons->begin(), Wleptons->end() )   ;//original
       // int nmuon = myMuon_new.size()                                        ;
       std::sort(myMuon_new.begin(), myMuon_new.end(), PtGreater_Muon())       ;
       
       // Jets
       // Handle<reco::CaloJetCollection > jets                   ;
       // edm::Handle < std::vector <reco::PFJet> > jets          ;
       // edm::Handle<edm::View<reco::PFJet> > jets                ;
       
       // edm::Handle<reco::PFJetCollection> jets               ;
       // iEvent.getByLabel("myJetProdLabel", jets)             ;
       // // reco::CaloJetCollection::const_iterator JetsProd          ;
       // reco::PFJetCollection::const_iterator JetsProd        ;
       
       // // reco::CaloJetCollection myJets(jets->begin(), jets->end()) ;
          // reco::PFJetCollection myJets(jets->begin(), jets->end()) ;
       // sort(myJets.begin(), myJets.end(), PtGreater())            ;
       
       // Get PFCandidate collection
       edm::Handle<edm::View<reco::PFCandidate> > pfCollection   ;
       iEvent.getByLabel("particleFlow",pfCollection)            ; 
       unsigned int pfCollectionSize = pfCollection->size()      ;
      // Ask for PfMuon consistency if requested
      
             bool muonFound = false;
       
     //instantiate helper class
     tbz::TbZUtility tbzHelper  ;
     //fill selected muons     
     int tt=0                   ;
     bool passPt=false          ;
     bool passEta=false         ;
     bool isTight=false         ;
     double TransMom_Muon       ;
     double ChargedHadronPt     ;
     double sumNeutralHadronEt  ;
     double sumPhotonEt         ;
     float sumPUPt              ;
     double Iso                 ;
         
       if(Wleptons->size() )m_muonCutFlow->Fill(3)      ; 
       int nmu= 0                                       ;
       for(reco::MuonCollection::const_iterator muon  = myMuon_new.begin()
                                               ; muon != myMuon_new.end(); ++muon)
       {
       
       edm::Ptr<reco::Muon> myMuon(&myMuon_new,tt)           ;
       reco::TrackRef tk = myMuon ->innerTrack()             ;
       
       bool isTightMuon   =   false                          ;
       bool isLooseMuon   =   false                          ;
       bool passLeadingPt =   false                          ;
       
       //muon.pfIsolationR04().sumChargedHadronPt;
       ChargedHadronPt    =   myMuon->pfIsolationR04().sumChargedHadronPt                                ; 
       sumNeutralHadronEt =   myMuon->pfIsolationR04().sumNeutralHadronEt                                ;
       sumPhotonEt        =   myMuon->pfIsolationR04().sumPhotonEt                                       ;
       sumPUPt            =   myMuon->pfIsolationR04().sumPUPt                                           ;

       // Proposed Configuration (DeltaBeta Correction) at cone of 0.4: 
       Iso = (ChargedHadronPt +  max(0. , sumNeutralHadronEt + sumPhotonEt - 0.5*sumPUPt))/myMuon->pt()  ; 
       
       cout<< "Muon_Isolation: "<< Iso <<endl     ;
       if(Iso > 0.12) continue                    ;
       cout<<"Muon_Iso_afterCut: "<< Iso <<endl   ;
       
       isTightMuon        =   tbzHelper.isTightMuon (  myMuon,  theVertex )            ;
       isLooseMuon        =   tbzHelper.isLooseMuon (  myMuon,  theVertex )            ; 
       
       if( nmu == 0 && muon->pt() > 20.) passLeadingPt = true                          ;
       cout << "Leading Muon True: " << passLeadingPt<<endl                            ;
       
       if( nmu == 0 && passLeadingPt ) pT_1st_LeadingMuon->Fill(muon->pt())            ;
       if( nmu == 1 && muon->pt()>= muonPtCut_ ) pT_2nd_LeadingMuon->Fill(muon->pt())  ;
       if ( nmu == 2 && muon->pt()>= muonPtCut_) pT_3rd_LeadingMuon->Fill(muon->pt())  ;
       
       if( muon->pt()<= muonPtCut_)             continue    ;
       TransMom_Muon = muon->pt()                           ;
       TransMom_Muon_all->Fill(TransMom_Muon)               ;
       passPt=true                                          ;
       
       // double Muon_Eta = muon->eta()                     ;
       // double Muon_Phi = muon->phi()                     ;
       
       if ( fabs(muon->eta())>= muonEtaCut_)    continue    ;
       // -------15-07-14---------------------------------------
       
       // Loose Electron count--- 150714 ----
         int NLooseMu = 0 ;
         if(isLooseMuon)
         {
         NLooseMu ++ ;
         tbzHelper.setNLooseElec(NLooseMu);
         }
         //------------------------------------
       
       
       
       
       for (unsigned int j=0; j<pfCollectionSize; j++) {
       
                   edm::RefToBase<reco::PFCandidate> pf = pfCollection->refAt(j);
 
                   // Check the muon is in the PF collection when required
                   bool thisIsTheMuon = false;
                   if (tk.isNonnull() && pf->trackRef()==tk) {
                         thisIsTheMuon = true;
                         muonFound = true    ;
                   }  
                   
}
             // Do not take this muon (under explicit request) if it is not a PfMuon
             if (!muonFound) continue;
        //-------------------------------------------------------
        
       //--------05-06-14----------------------------------------
       
        // double DR_jet_Muon =1000.;
        // for(JetsProd=jets->begin(); JetsProd != jets->end(); ++JetsProd)
         // {      
         
            // double PFJetNHEF = JetsProd ->neutralHadronEnergyFraction()  ;
            // double PFJetCHEF = JetsProd ->chargedHadronEnergyFraction()  ;
            // double PFJetNEMF = JetsProd ->neutralEmEnergyFraction()      ;
            // double PFJetCEMF = JetsProd ->chargedEmEnergyFraction()      ;     
     
            // if(PFJetNHEF > 0.99 )  continue                      ;
            // if(PFJetCHEF <= 0. )   continue                      ;
            // if(PFJetNEMF > 0.99 )  continue                      ;
            // if(PFJetCEMF > 0.99 )  continue                      ;
            
            // double   jets_eta       =   JetsProd->eta()          ;
            // double   jets_phi       =   JetsProd->phi()          ;               
            // double   DEta_jet_Muon  =   jets_eta - Muon_Eta      ;
            // double   DPhi_jet_Muon  =   jets_phi - Muon_Phi      ;
            
            // DR_jet_Muon             =  sqrt((DEta_jet_Muon)*(DEta_jet_Muon)
                                       // +(DPhi_jet_Muon)*(DPhi_jet_Muon))                       ;
                                       
            // cout<<"deltaR_JetMu_Before_**break**_statement_without_ISO: "<<DR_jet_Muon<<endl   ;
            // H_deltaR_JetMuon1 ->Fill(DR_jet_Muon)                                              ;
            // // if(DR_jet_Muon < 0.5)  break                                                    ;
            // cout<<"deltaR_JetMu_After_**break**_statement_without_ISO: "<<DR_jet_Muon<<endl    ;
               
            // }// end of jet loop-for deltR
            
             // // if(DR_jet_Muon < 0.5)  continue                                 ;                     
             // cout<<"After_*DR_JetMu(continue)*_statement_without_ISO" <<endl    ;
             // H_deltaR_JetMuon2 ->Fill(DR_jet_Muon)                              ;
            //-------------------------------------------------------------------------------------------
               //pfIsolationR04()
             double Isolation = muon->isolationR03().sumPt             ;
             Isolation_mu->Fill(Isolation)                             ; 
             if(muon->isolationR03().sumPt >= .2)    continue          ; //2.5
             passEta=true                                              ;
             if(isTightMuon ==false)                  continue         ; // we should use this from very start ??
             
             isTight=true                                              ;
             
             //--- number of tight muons ------
             nmu++                                                     ;       
             cout <<"numberMuon: " << nmu << endl                      ;
             Number_of_tightMuons->Fill(nmu)                           ;
             // //--------------------------------
             
             cout<<"isTight muon "<< isTightMuon<<endl                 ;
             wmuons->push_back(*muon)                                  ;
             if( wmuons->size()<=0.)                  continue         ;
             tt++                                                      ;                   
             
       }
       
       //========================================================
       
       if(passPt)m_muonCutFlow->Fill(4)                               ;
       if(passEta && passPt)m_muonCutFlow->Fill(5)                    ;
       if(isTight& passEta && passPt)m_muonCutFlow->Fill(6)           ;       
       //has one lepton                                               
       if(wmuons->size() ==1 )m_muonCutFlow->Fill(7)                  ;
       if(wmuons->size() ==2 )m_muonCutFlow->Fill(8)                  ;
       if(wmuons->size() ==3 )m_muonCutFlow->Fill(9)                  ;
       if(wmuons->size() >3 )m_muonCutFlow->Fill(10)                  ; 
       
       //==========================================================
       cout<<"size of selected muons : "<<wmuons->size()<<endl        ;
       iEvent.put(ptr)                                                ;
       //==========================================================
}
// ---- method called once each job just before starting event loop
void 
TbzMuonProducer::beginJob()
{

   edm::Service<TFileService> fs                                                                  ;
   m_muonCutFlow      = fs->make<TH1D> ("muonsProd_cutFlow","muons cut flow",15,0.,15.)           ;
   inv_Z_mass         = fs->make<TH1D> ("inv_Z_mass","invariant mass of Z boson plot",60,60.,120.);
   pT_Z               = fs->make<TH1D> ("pT_Z","pT Z plot",50,0.,100.)                            ;
   w_mT               = fs->make<TH1D> ("trans_mass","W transverse  mass plot",100,20.,120.)      ;
   met                = fs->make<TH1D> ("met","MET plot",100,0.,100.)                             ;
   acop               = fs->make<TH1D> ("acop","acop plot",100,-3.14,3.14)                        ;
   pT_1st_LeadingMuon = fs->make<TH1D> ("pT_1st_LeadingMuon","pT 1st Leading Mu",100,0.,200.)     ;
   pT_2nd_LeadingMuon = fs->make<TH1D> ("pT_2nd_LeadingMuon","pT 2nd Leading Mu",100,0.,200.)     ;
   pT_3rd_LeadingMuon = fs->make<TH1D> ("pT_3rd_LeadingMuon","pT 3rd Leading Mu",100,0.,200.)     ;
   //----------                                                                             
   TransMom_Muon_all  = fs->make<TH1D> ("TransMom_Muon_all", "Sum_allMuon_pt", 100, 0. , 200)     ;
   //---                                                                                           
   Isolation_mu       = fs->make<TH1D> ("Isolation_mu", "Isolation_muon", 100, 0., 10.)           ;
   
   H_deltaR_JetMuon1  = fs->make<TH1D> ("H_deltaR_JetMuon1","deltaR_JetMuon1_befroe_break", 100, -10., 10. );
   H_deltaR_JetMuon2  = fs->make<TH1D> ("H_deltaR_JetMuon2","deltaR_JetMuon1_after_continue", 100, -10., 10. );
   //----16-06-2014---
   Number_of_tightMuons = fs->make<TH1D> ("Number_of_tightMuons","Number_of_tightMuons", 10, 0, 10);   
   //---
   NpVertices          = fs->make<TH1D> ("NpVertices","Primary_Vertices",10,0,10);
   //---
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TbzMuonProducer::endJob()
{
  std::vector<std::string > cuts                          ;
  cuts.push_back("All events........");
  cuts.push_back("# events after good vertex........")    ;  
  cuts.push_back("# events after muon trigger........")   ;
  cuts.push_back("# events with muons.........")          ;
  cuts.push_back("# events with muon pT >25GeV..........");
  cuts.push_back("# events with eta <2.1 ..........")     ;
  cuts.push_back("# events with Tight muon cut..........");
  cuts.push_back("# events with 1 muons..........")       ;
  cuts.push_back("# events with 2 muons..........")       ;
  cuts.push_back("# events with 3 muons..........")       ;
  cuts.push_back("# events with >3 muons..........")      ;


  for(unsigned  i=0; i<cuts.size(); i++)
    {          
    std::cout<< "i : "<< i<<" ......... "<<cuts.at(i)<<"......................" << m_muonCutFlow->GetBinContent(i+1)<<std::endl;
    }
  std::cout<< ".........................End of muon CUT Flow in PRODUCER......................" << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void 
TbzMuonProducer::beginRun(edm::Run&, edm::EventSetup const&)
{

}

// ------------ method called when ending the processing of a run  ------------
void 
TbzMuonProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TbzMuonProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TbzMuonProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TbzMuonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TbzMuonProducer);
