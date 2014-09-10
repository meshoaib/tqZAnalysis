#ifndef TbZUtility_h 
#define TbZUtility_h

/**_________________________________________________________________
   class:   utility class to implement all common methods
By AA in July, 2013
__________________________________________________________**/


#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

using namespace edm;
using namespace std;
using namespace reco;

namespace tbz{

class TbZUtility {

public:
  TbZUtility()
  {
  nLooseMu =   0.;
 nLooseElec = 0.;
  };
  ~TbZUtility(){};
  //---150714---
  void setNLooseMu(int&)    ;
  int  getNLooseMu()        ;
  void setNLooseElec(int&)  ;
  int getNLooseElec()       ;
  //----
  
 // bool isTightMuon (  edm::Ptr<reco::PFCandidateCollection>,  const reco::Vertex* theVertex );
 bool isTightMuon (  edm::Ptr<reco::Muon>,  const reco::Vertex* theVertex );
 bool isMediumMuon( edm::Ptr<reco::Muon>,  const reco::Vertex* theVertex  );
 bool isLooseMuon (  edm::Ptr<reco::Muon>,  const reco::Vertex* theVertex );


 bool isTightElectron (  edm::Ptr<reco::GsfElectron> );
 bool isMediumElectron( edm::Ptr<reco::Electron>  );
 bool isLooseElectron (  edm::Ptr<reco::Electron> );


 bool sfosMassInRange(Double_t , Double_t, std::vector<double> &);
 bool esfosMassInRange(Double_t , Double_t, std::vector<double> &);

 // void makePairs(const std::vector<reco::PFCandidateCollection>&, std::vector<double>&, std::pair<int, int>&, std::vector< reco::NamedCompositeCandidate >& );
  void makePairs(const std::vector<reco::Muon>&, std::vector<double>&, std::pair<int, int>&, std::vector< reco::NamedCompositeCandidate >& );
 void makeEPairs(const std::vector<reco::GsfElectron>&, std::vector<double>&, std::pair<int, int>&, std::vector< reco::NamedCompositeCandidate >& );
 ///overload for edm handle
 // void makePairs( std::vector<reco::Muon>*, std::vector<double>&, std::pair<int, int>&, std::vector< reco::NamedCompositeCandidate >& );
//std::vector<reco::Muon>
//edm::Ptr<reco::Muon>
//if(doTruthMatch_)
//{
 void isTrue(  edm::Ptr<reco::GsfElectron> &
              ,edm::Handle<std::vector<reco::GenParticle> > &
              ,std::vector<reco::GenParticle> &
              ,double &, double &, double &, double &);
              
 void isTrueWenu(  edm::Ptr<reco::GsfElectron> &
              ,edm::Handle<std::vector<reco::GenParticle> > &
              ,std::vector<reco::GenParticle> &
              ,double &, double &, double &, double &);              
//const std::vector<reco::GsfElectron> &e_in
//const std::vector<reco::GsfElectron> &const std::vector<GenParticle> &genParticles
//reco::PFCandidateCollection
 void doTruthMatching( edm::Ptr<reco::Muon>&
                       ,edm::Handle< std::vector<reco::GenParticle> >&
                       ,std::vector<reco::GenParticle> &
                       ,double &, double &, double &, double &);
                       
                  
                  
void doWTruthMatching(   edm::Ptr<reco::Muon>&
                       ,edm::Handle< std::vector<reco::GenParticle> >&
                       ,std::vector<reco::GenParticle> &
                       ,double &, double &, double &, double &);
                  
                  
//}
 //void  DeltaR(const std::vector<reco::GsfElectron> &ev,  const std::vector<reco::CaloJetCollection > &JetsProd, double &, double &, double &, double &); 

private:

int nLooseElec ;
int nLooseMu   ;
    // edm::Ptr<Muon> muon_;
 // edm::Ptr<PFCandidateCollection> muon_;
    // edm::Ptr<MET>  neutrino_;  
 // edm::Ptr<PFMET>  neutrino_;
};

//typedef std::vector<reco::TbZUtitlity> TbZUtitlityCollection;
}
#endif
