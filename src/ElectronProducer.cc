#include "MyAnalysis/TbZ/interface/ElectronProducer.h"
#include "MyAnalysis/TbZ/interface/TbZUtility.h"

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps    ;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals    ;

// constructors and destructor
ElectronProducer::ElectronProducer(const edm::ParameterSet& iConfig)
{
	// get input parameters
    electronsInputTag_      = iConfig.getParameter<edm::InputTag>("electronsInputTag")              ;
    conversionsInputTag_    = iConfig.getParameter<edm::InputTag>("conversionsInputTag")            ;
    beamSpotInputTag_       = iConfig.getParameter<edm::InputTag>("beamSpotInputTag")               ;
    rhoIsoInputTag          = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag")                 ;
    primaryVertexInputTag_  = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag")          ;
    isoValInputTags_        = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags")  ;
    electronPtCut_          = iConfig.getParameter<double>("electronPtCut")                         ;
    electronEtaCut_         = iConfig.getParameter<double>("electronEtaCut")                        ;
    //----
    muonsInputTag_          = iConfig.getParameter<edm::InputTag>("muonsInputTag")                  ;
    //----
     // debug
     printDebug_         = iConfig.getParameter<bool>("printDebug")                                 ;
     // produces vector of electrons                                                                
     produces<std::vector<reco::GsfElectron> >()                                                    ;

}

ElectronProducer::~ElectronProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// member functions
// ------------ method called to produce the data  -----------

void
ElectronProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    edm::Handle<reco::PFJetCollection> jets                  ;
    iEvent.getByLabel("myJetProdLabel", jets)                ;
    
    reco::PFJetCollection::const_iterator JetsProd           ;
    reco::PFJetCollection myJets(jets->begin(), jets->end()) ;
   // sort(myJets.begin(), myJets.end(), PtGreater())        ;
   
    edm::Handle<reco::GsfElectronCollection> els_h           ;// original
    iEvent.getByLabel(electronsInputTag_, els_h)             ;
    
    
    //---- Muons ---
     edm::Handle<vector <reco::Muon> > muonColl               ; //original reco::PFCandidate
    iEvent.getByLabel(muonsInputTag_, muonColl)              ;
    
    // conversions
    edm::Handle<reco::ConversionCollection> conversions_h    ;
    iEvent.getByLabel(conversionsInputTag_, conversions_h)   ;

    // iso deposits
    IsoDepositVals isoVals(isoValInputTags_.size());
    
    for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
        iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
    }

    // beam spot
    edm::Handle<reco::BeamSpot> beamspot_h                     ;
    iEvent.getByLabel(beamSpotInputTag_, beamspot_h)           ;
    const reco::BeamSpot &beamSpot = *(beamspot_h.product())   ;

    // vertices
    edm::Handle<reco::VertexCollection> vtx_h                  ;
    iEvent.getByLabel(primaryVertexInputTag_, vtx_h)           ;

    // rho for isolation
    edm::Handle<double> rhoIso_h                               ;
    iEvent.getByLabel(rhoIsoInputTag, rhoIso_h)                ;
    double rhoIso = *(rhoIso_h.product())                      ;

    // this will be the new object collection // reco::PFCandidate //reco::PFCandidateCollection
     std::vector<reco::GsfElectron> * electronprod = new std::vector<reco::GsfElectron>()  ; //original
    std::auto_ptr<std::vector<reco::GsfElectron> > ptr(electronprod)                      ; // original
    
    //total events
        electronCutFlow_producer->Fill(0)                                                 ;
        // loop on electrons                                                              
        int nelec = 0                                                                     ;
        unsigned int n = els_h->size()                                                    ;
        bool    passPt            = false                                                 ;
        bool    passEta           = false                                                 ;
        bool    isTight           = false                                                 ;
        bool    passLeadingElec   = false                                                 ;
        double  TransMom_Elec                                                             ;
        //----                                                                            
        //vector<reco::Muon> muon_preSel_producer                                         ;    
        //----                                                                            
        std::cout<<"Electron Producer: size of electron without cuts "<< els_h->size()<<std::endl            ;
        
         int nmuons           = 0. ;
         int Number_Muons     = 0. ;
         int NLooseEle        = 0 ;
         
         for(unsigned int mi=0; mi < muonColl->size(); ++mi)
         {          
            
            nmuons++                                          ;
            H1_no_OfMuon->Fill(nmuons)                        ; 
            
            Number_Muons               = muonColl->size()     ;
            Number_tightMuons_Producer ->Fill(Number_Muons)   ;  
           
         }  
         
        for(unsigned int i = 0; i < n; ++i) {
        
        //total electrons
        electronCutFlow_producer->Fill(1)       ;
        
        // get reference to electron            
        reco::GsfElectronRef ele(els_h, i)      ;
        
        // get particle flow isolation
        double iso_ch =  (*(isoVals)[0])[ele]   ;
        double iso_em = (*(isoVals)[1])[ele]    ;
        double iso_nh = (*(isoVals)[2])[ele]    ;
        // test ID

// working points
bool veto       = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::VETO,ele, conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
bool loose      = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, ele, conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
bool medium     = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, ele, conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
bool tight      = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, ele, conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);

// eop/fbrem cuts for extra tight ID
bool fbremeopin = EgammaCutBasedEleId::PassEoverPCuts(ele);
// cuts to match tight trigger requirements
bool trigtight = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERTIGHT, ele)  ;
// for gg2011 WP70 trigger                                                                     
bool trigwp70 = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERWP70, ele)    ;

         if(nelec == 0 && ele->pt() > 20.) passLeadingElec =true                               ;
         if (nelec == 0 && passLeadingElec) Leading_Elec_pt_1st->Fill(ele->pt())               ;
         if (nelec == 1 && ele->pt() >= electronPtCut_ ) Leading_Elec_pt_2nd->Fill(ele->pt())  ;
         if (nelec == 2 && ele->pt() >= electronPtCut_ ) Leading_Elec_pt_3rd->Fill(ele->pt())  ;
         nelec++                                                                               ;        
        
        double elec_eta = ele->eta()   ; 
        double elec_phi = ele->phi()   ;
        double elec_pt  = ele->pt()    ;
        //--------------------------------------deltaR between electron-muon-----
        
         double MuonEta_Prodcr      = 1000                     ;
         double MuonPhi_Prodcr      = 1000                     ;
         double deltaPhi_ElecMu     = 1000                     ;
         double deltaEta_ElecMu     = 1000                     ;
         double deltaR_ElecMu       = 1000                     ;
         
         
         // Loose Electron count--- 150714 ----
         
         tbz::TbZUtility tbzHelper     ;
         if(loose)
         {
         NLooseEle ++ ;
         }
         cout<<"Loose_electrons_in_producer: " <<  NLooseEle <<endl;
         //------------------------------------
         
         // bool   MuonCollectionSize  = false                     ;
         // bool   ElecCollectionSize  = false                     ;
         
         // if (els_h->size() > 0)     ElecCollectionSize  = true  ;
         // if (muonColl->size() > 0)  MuonCollectionSize  = true  ;
        // if( MuonCollectionSize == true && ElecCollectionSize  == true )
        // {
        for(unsigned int mi=0; mi < muonColl->size(); ++mi)
            {          
         
           edm::Ptr<reco::Muon> myMuon(muonColl,mi)                                  ;  //reco::PFCandidate //reco::PFCandidateCollection
           MuonEta_Prodcr = myMuon->eta()                                               ;
           MuonPhi_Prodcr = myMuon->phi()                                               ;
           deltaEta_ElecMu = elec_eta - MuonEta_Prodcr                                  ;
           deltaPhi_ElecMu = elec_phi - MuonPhi_Prodcr                                  ;
           
           // if (fabs(deltaPhi_ElecMu) > 3.14)                                            
           // deltaPhi_ElecMu = 6.28 - fabs(deltaPhi_ElecMu)                            ;
           // float DELTAR = sqrt(pow(dphi,2)+pow(deta,2))*1.0;                        
           
           // deltaR_ElecMu   = sqrt(pow(deltaPhi_ElecMu,2) + pow(deltaEta_ElecMu,2))*1.0  ;
           
          deltaR_ElecMu    =  sqrt(( deltaEta_ElecMu)*(deltaEta_ElecMu)
                                             +( deltaPhi_ElecMu)*(deltaPhi_ElecMu))     ;
                                             
           H1_deltaR_ElecMu ->Fill(deltaR_ElecMu)                                       ;
           
           if(deltaR_ElecMu < 0.1) break                                                ;  //0.3 before
           
               }         
               
            if(deltaR_ElecMu < 0.3) continue             ;            //we were testing with 0.1
            deltaR_ElecMu_Cut ->Fill(deltaR_ElecMu)      ;      
           
        //--------ISO-----------------------------------------
        // effective area for isolation
        float AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, elec_eta, ElectronEffectiveArea::kEleEAData2012);

         // apply to neutrals
         double rhoPrime = std::max(rhoIso, 0.0);
         double iso_n = std::max(iso_nh + iso_em - rhoPrime * AEff, 0.0);

         // compute final isolation
         double iso = (iso_n + iso_ch) / elec_pt            ;
         //if (iso > 2.0) continue;                         
         Electron_ISO->Fill(iso)                            ;  // filling histogram before applying isolation cut on electrons 
         if (iso > 2.0) continue                            ; // Isolation cut applied here
         
       //-----------------------------------------------------------------------------  
       
         //pt cut                                           
         if(ele->pt()<= electronPtCut_ )continue            ;
         //if(ele->pt()<=25.)continue;                      
         TransMom_Elec = ele->pt()                          ;
         Trans_Momt_Ele->Fill(TransMom_Elec)                ;
         passPt=true                                        ;      
         //if ( fabs(ele->eta())>=2.1) continue;
         if ( fabs(ele->eta())>=electronEtaCut_) continue   ;
         passEta=true                                       ;
         //eta cut
      
      if(!tight) continue;
      //if(!loose) continue;
      //------------------- 2nd------------------------------------------------------------------
         double var0          =  ele -> pfIsolationVariables().chargedHadronIso        ;
         double NeutralHadIso =  ele -> pfIsolationVariables().neutralHadronIso        ;
         double photonIso     =  ele -> pfIsolationVariables().photonIso               ;
         double ele_pt_new    =  ele -> pt()                                           ;
         //---
         double ElecEnergy    =  ele->energy() ;
         cout << "Elec_Producer ---- Electron_Energy: " << ElecEnergy <<endl;
         //---
         double relIso_elec    = (var0 + NeutralHadIso + photonIso) /ele_pt_new        ;
         cout<<"Elec_Producer ---- Isolation_Before_Electron_Cut_With_No_DeltaRcut: "<<relIso_elec <<endl ;
         
         H1_ElectronIso ->Fill(relIso_elec)                                            ;
         if(relIso_elec > 0.10)                                      continue          ;//0.2
         
         cout<<"Elec_Producer ---- Isolation_After_Electron_Cut_With_No_DeltaRcut: "<<relIso_elec <<endl  ;
            H1_ElectronIso_afterCUT ->Fill(relIso_elec);
            
      //-----------------------------------------------------------------------------------------     
          float  Electron_IsopT = ele-> dr03TkSumPt()                      ;
         float Elec_EcalRechit =  ele-> dr03EcalRecHitSumEt()              ;
         float Elec_HcalSumEt =  ele-> dr03HcalTowerSumEt()               ;
        // float Elec_HcalSumEt2 =  ele-> dr03HcalDepth2TowerSumEt()         ;

         double relIso_elec1 = (Electron_IsopT + Elec_EcalRechit + Elec_HcalSumEt )/ele -> pt() ;
         cout << "Elec_Producer ---- NewIso: "<<relIso_elec1<<endl ;
         // if(relIso_elec1 > 0.2) continue     ;
         cout << "Elec_Producer ---- NewIso_after_cut: "<<relIso_elec1<<endl;
      
         double DR_jet_elec = 1000.;
         
         for(JetsProd=jets->begin(); JetsProd != jets->end(); ++JetsProd)
         {      
         
            double PFJetNHEF = JetsProd ->neutralHadronEnergyFraction()  ;
            double PFJetCHEF = JetsProd ->chargedHadronEnergyFraction()  ;
            double PFJetNEMF = JetsProd ->neutralEmEnergyFraction()      ;
            double PFJetCEMF = JetsProd ->chargedEmEnergyFraction()      ;     
     
            if(PFJetNHEF > 0.99 )  continue                              ;
            if(PFJetCHEF <= 0. )   continue                              ;
            if(PFJetNEMF > 0.99 )  continue                              ;
            if(PFJetCEMF > 0.99 )  continue                              ;
     
            
            double   jets_eta1             =  JetsProd->eta()            ;
            double   jets_phi1             =  JetsProd->phi()            ;               
            double   DEta_jet_elec         =  jets_eta1 - ele->eta()     ;
            double   DPhi_jet_elec         =  jets_phi1 - ele->phi()     ;
            double   JetsEnergy            =  JetsProd->energy()         ;
            cout << "Elec_Producer ---- Jets_Energy: " << JetsEnergy <<endl                 ;
            
            DR_jet_elec                    =  sqrt(( DEta_jet_elec)*(DEta_jet_elec)
                                             +( DPhi_jet_elec)*(DPhi_jet_elec))          ;
                                             
              cout<<"Elec_Producer ---- deltaR_Before_**break**_statement_without_ISO: "<<DR_jet_elec<<endl ;
              H_deltaR_jetElec_bforeBREAK ->Fill(DR_jet_elec)                            ;
               
               if(DR_jet_elec < 0.5)
               {
               
               cout << "Elec_Producer ---- Electron_Energy1: " << ElecEnergy <<endl            ;
               cout << "Elec_Producer ---- Jets_Energy1: "     << JetsEnergy <<endl            ;
               H2_ElescJet_Energy           ->Fill(JetsEnergy,ElecEnergy)   ;
               H2_ElecEta_JetsEta           ->Fill(jets_eta1,elec_eta)      ;
               H2_ElecPhi_JetsPhi           ->Fill(jets_phi1,elec_phi)      ;
               
               }
               
               // if(DR_jet_elec == 0) cout<<"delta_Eta: " << DEta_jet_elec <<"delta_Phi: "<<DPhi_jet_elec<<endl;
               if(DR_jet_elec < 0.5) H_delta05R_jetElec ->Fill(DR_jet_elec);
               
               // if(DR_jet_elec < 0.5)  break ;            
              
            }// end of jet loop-for deltR
            
               if(DR_jet_elec > 0.5) deltaR_grter05_jetElec ->Fill(DR_jet_elec);
               //---------------------            
            // if(DR_jet_elec < 0.5)  continue                              ;  
            cout<<"Elec_Producer ---- After_DR_jet_elec(continue): " << DR_jet_elec <<endl  ;
            H_deltaR_jetElec_aftrBREAK->Fill(DR_jet_elec)                ;
            //isTight electron
                              
            isTight=true ;
            // print decisions

        if (printDebug_) {
        
            printf("%u %u %u : ",       iEvent.id().run(), iEvent.luminosityBlock(), iEvent.id().event());
            printf("veto(%i), ",        veto)         ;
            printf("loose(%i), ",       loose)        ;
            printf("medium(%i), ",      medium)       ;
            printf("tight(%i), ",       tight)        ;
            printf("trigtight(%i), ",   trigtight)    ;
            printf("trigwp70(%i), ",    trigwp70)     ;
            printf("fbremeopin(%i)\n",  fbremeopin)   ;
            printf("eta (%f)\n",  ele->eta())         ;
            printf("pt (%f)\n",  ele->pt() )          ;            
        
    }

      std::cout<<"Electron_Producer: electron size before push_back: "<< electronprod->size()<<std::endl;
      
      electronprod->push_back(*ele)            ;
      if( electronprod->size()<=0.) continue       ;
      std::cout<<"Electron_Producer: electron size after push_back: "<< electronprod->size()<<std::endl;
      
      
      
      }
      
            if(passPt)electronCutFlow_producer->Fill(2)                         ;
            if(passEta && passPt)electronCutFlow_producer->Fill(3)              ;
            if(isTight && passEta && passPt)electronCutFlow_producer->Fill(4)   ;
            //has one electron                                                  
            if(electronprod->size() ==1 )electronCutFlow_producer->Fill(5)      ;
            if(electronprod->size() ==2 )electronCutFlow_producer->Fill(6)      ;
            if(electronprod->size() ==3 )electronCutFlow_producer->Fill(7)      ;
            if(electronprod->size() >3 )electronCutFlow_producer->Fill(8)       ;

            iEvent.put(ptr);
            std::cout<<"Electron_Producer: electron registered in this event size : "<< electronprod->size() << std::endl;

}

// ------------ method called once each job just before starting event loop  ------------
void 
ElectronProducer::beginJob()
{
   edm::Service<TFileService> fs;
   electronCutFlow_producer = fs->make<TH1D> ("electronsProd_cutFlow","electrons cut flow",15,0.,15.);
   Leading_Elec_pt_1st = fs->make<TH1D> ("Leading_elec_pt_1st", "1st_elec_pt", 100, 0, 200);
   Leading_Elec_pt_2nd = fs->make<TH1D> ("Leading_elec_pt_2nd", "2nd_elec_pt", 100, 0, 200);
   Leading_Elec_pt_3rd = fs->make<TH1D> ("Leading_elec_pt_3rd", "3rd_elec_pt", 100, 0, 200);
   //----------------
   Trans_Momt_Ele = fs->make<TH1D> ("Trans_Momt_Ele","Sum_All_elePt",100, 0., 200.);
   //----------------------------------------------------------------------------------
   Electron_ISO = fs->make<TH1D>  ("Electron_ISO", "Isolation of electrons", 100, 0.0, 10.0);
   
   H1_ElectronIso              = fs->make<TH1D>  ("H1_ElectronIso", "H1_ElectronIso", 100, -10., 10.0);
	H1_ElectronIso_afterCUT     = fs->make<TH1D>  ("H1_ElectronIso_afterCUT", "H1_ElectronIso_afterCUT", 100, -10., 10.0);
   H_deltaR_jetElec_bforeBREAK = fs->make<TH1D>  ("H_deltaR_jetElec_bforeBREAK", "H_deltaR_jetElec_bforeBREAK", 10, -1., 9.0);
   H_deltaR_jetElec_aftrBREAK  = fs->make<TH1D>  ("H_deltaR_jetElec_aftrBREAK", "H_deltaR_jetElec_aftr_Continue", 10, -1., 9);
   //---
   H2_ElescJet_Energy          = fs->make<TH2D>  ("H2_ElescJet_Energy","H2_ElescJet_Energy",100,0,400, 100, 0, 400);
   //---
   H2_ElecEta_JetsEta          = fs->make<TH2D> ("H2_ElecEta_JetsEta", "H2_ElecEta_JetsEta",100, -3, 3, 100, -3,3);
   H2_ElecPhi_JetsPhi          = fs->make<TH2D> ("H2_ElecPhi_JetsPhi", "H2_ElecPhi_JetsPhi",100, -3, 3, 100, -3,3);
   //----
   Number_tightMuons_Producer  = fs->make<TH1D> ("Number_tightMuons_Producer","Number_tightMuons_Producer",10,0,10);
   //----
   H1_deltaR_ElecMu            = fs->make<TH1D> ("H1_deltaR_ElecMu","H1_deltaR_ElecMu",100, -10, 10) ;
   H1_no_OfMuon                = fs->make<TH1D> ("H1_no_OfMuon","H1_no_OfMuon", 10, 0, 10)           ;
   deltaR_ElecMu_Cut           = fs->make<TH1D> ("deltaR_ElecMu_Cut","deltaR_ElecMu_Cut",100,-10,10) ;
   //
   H_delta05R_jetElec         = fs->make<TH1D> ("delta05R_jetElec","delta05R_jetElec",3,-1,1);
   deltaR_grter05_jetElec     = fs->make<TH1D>  ("deltaR_grter05_jetElec","deltaR_grter05_jetElec",10,0,10);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronProducer::endJob() {


	std::vector<std::string > ecuts;
	ecuts.push_back("All events........");
	ecuts.push_back("# events with electrons.........");
	ecuts.push_back("# events with electron pT >20GeV..........");
	ecuts.push_back("# events with eta <2.1 ..........");
	ecuts.push_back("# events with Tight electron cut..........");
	ecuts.push_back("# events with 1 electron..........");
	ecuts.push_back("# events with 2 electron..........");
	ecuts.push_back("# events with 3 electron..........");
	ecuts.push_back("# events with >3 electron..........");


	for(unsigned  a=0; a<ecuts.size(); a++){

	//if( i<cuts.size() )

      cout<< "a : "<< a<<" ......... "<<ecuts.at(a)<<"......." << electronCutFlow_producer->GetBinContent(a+1)<<endl;

    }

	cout<< ".........................End of Electron CUT Flow in PRODUCER......................" << endl;


}

// ------------ method called when starting to processes a run  ------------
void 
ElectronProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ElectronProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ElectronProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ElectronProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronProducer);
