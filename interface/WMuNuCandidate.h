

#ifndef TbZ_WMuNuCandidate_h 
#define TbZ_WMuNuCandidate_h

#include <map>
#include <memory>

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "CommonTools/CandUtils/interface/CandCombiner.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
//----------------------------------------------------
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
//----------------------------------------------------
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

         using namespace edm    ;
         using namespace std    ;
         using namespace reco   ;
         
namespace reco{

class WMuNuCandidate : public reco::CompositeCandidate {         
public:
  WMuNuCandidate();
  // WMuNuCandidate(edm::Ptr<reco::PFCandidateCollection>,edm::Ptr<reco::PFMET>); //reco::PFMET//reco::MET
  WMuNuCandidate(edm::Ptr<reco::Muon>,edm::Ptr<reco::MET>);
  //WMuNuCandidate(edm::Ptr<reco::GsfElectron>,edm::Ptr<reco::MET>);

  ~WMuNuCandidate();
 //WARNING: W Candidates combine the information from a Muon with the (px,py) information of the MET as the Neutrino
      // --> There is no Pz information!!!!
      // Be very careful when using the default Candidate functions (.mass, .mt, .et, etc). They may not be what you are looking for :-).

      // Example: Candidates have a mt() function which computes the tranverse mass from E & pz.
      // As MET does not have pz information... WMuNuCandidates have an alternative function used in the WMuNu Inclusive Analysis
      // to compute mt just from px, py:
      // Transverse Mass from px, py:
      double massT() const;
      // Transverse Energy from px, py:
      double eT() const;

      // Acoplanarity between the muon and the MET
      double acop() const;

      const reco::Muon & getMuon()     const {return *muon_;}                   ;
         // const reco::PFCandidateCollection & getMuon()     const {return *muon_;}  ;
         const reco::MET  & getNeutrino() const {return *neutrino_;}               ;//original
         // const reco::PFMET  & getNeutrino() const {return *neutrino_;}               ;


private:


      edm::Ptr<Muon> muon_; // original
      // edm::Ptr<PFCandidateCollection> muon_;
      edm::Ptr<MET>  neutrino_; // original
      // edm::Ptr<PFMET>  neutrino_; 
      
};


typedef std::vector<reco::WMuNuCandidate> WMuNuCandidateCollection;
}


#endif

