#include <iostream>
//#include "FirstAnalysis/TBZAnalysis/interface/TbZUtility.h"
#include "MyAnalysis/TbZ/interface/TbZUtility.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
//-----
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//-----
//---
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//--------

// TbZUtility::TbZUtility(){}
// TbZUtility::~TbZUtility(){}
    using namespace  std    ;
    using namespace  edm    ;
    using namespace  HepMC  ;
    using namespace  reco   ;
 //constructor
 // tbz::TbZUtility::TbZUtility()
 // {
 // nLooseMu =   0.;
 // nLooseElec = 0.;
// }
bool tbz::TbZUtility::isTightMuon(edm::Ptr<reco::Muon> muon, const reco::Vertex* theVertex)// original
// bool tbz::TbZUtility::isTightMuon(edm::Ptr<reco::PFCandidateCollection> muon, const reco::Vertex* theVertex)
{
    bool isTight =false; 
    if(!muon || !theVertex->isValid()) 
    {
    std::cout<<" empty muon or theVertex pointer inside TbZUtility"<<std::endl; 
    return false;
    } 
    if(muon->isGlobalMuon() && muon->isPFMuon()) 
    {
    reco::TrackRef track = muon->globalTrack();
    if (!track )
    {
    std::cout<<"WARNING NO  ASSOCIATED TRACK FOR  MUONS!!!!"<<std::endl;
    return false;
    }

    if(muon->globalTrack()->normalizedChi2() < 10.)
      {
      if(muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0) 
      {
      if(muon->numberOfMatchedStations() > 1) 
        {
        if(fabs(muon->muonBestTrack()->dxy(theVertex->position())) < 0.01)//0.2
          {
          if(fabs(muon->muonBestTrack()->dz(theVertex->position())) < 0.1)// 0.5 
            {
            if(muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0) 
             {
             if(muon->track()->hitPattern().trackerLayersWithMeasurement() > 8) //5
               {
               if(muon->isolationR03().sumPt < .2 ) //2.5
                 {
                  isTight = true;     
                 }//isolation
               }//tracker layers 
             }//Pixel hits
           }//dz
         }//dxy
       }//matched station
     }//hitPattern
   }//global track
 }//global or PFMuon
return isTight;

}

bool tbz::TbZUtility::isMediumMuon(edm::Ptr<reco::Muon>, const reco::Vertex* theVertex )
{
     bool isMedium=true;
     std::cout<<"WARNING MEIDUM MUON CUTS NOT IMPLEMENTED YET!!!!"<<std::endl;
     return isMedium;
}
bool tbz::TbZUtility::isLooseMuon(edm::Ptr<reco::Muon>  muon, const reco::Vertex* theVertex )
{
    bool isLoose =false;
    //std::cout<<"WARNING LOOSE  MUON CUTS NOT IMPLEMENTED YET!!!!"<<std::endl;
    bool isTight =false;
         if(!muon || !theVertex->isValid())
      {
        std::cout<<" empty muon or theVertex pointer inside TbZUtility"<<std::endl;
        return false;
      }

        //tight muon selection cuts
        //if(muon->pt()<=25 || fabs(muon->eta())>=2.1) continue;

        if(muon->isGlobalMuon() && muon->isPFMuon())
	  {
	    /*

	    reco::TrackRef track = muon->globalTrack();
	    if (!track )
	      {
		std::cout<<"WARNING NO  ASSOCIATED TRACK FOR  MUONS!!!!"<<std::endl;
		return false;
	      }

	    if(muon->globalTrack()->normalizedChi2() < 10.)
	      {
	      
	      if(muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0)
		  {
		    if(muon->numberOfMatchedStations() > 1)
		      {
			if(fabs(muon->muonBestTrack()->dxy(theVertex->position())) < 0.2)
			  {
			    if(fabs(muon->muonBestTrack()->dz(theVertex->position())) < 0.5)
			      {
				if(muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0)
				  {
				    if(muon->track()->hitPattern().trackerLayersWithMeasurement() > 5)
				      {
		*/
				     // if(muon->isolationR03().sumPt < 2.5)
					 // {

					    isLoose = true;
					    //}//isolation
					    /*
				      }//tracker layers
				  }//Pixel hits
			      }//dz
			  }//dxy
		      }//matched station
		  }//hitPattern
								       

	      }//global track

*/

	  }//global or PFMuon


 return isLoose;
}
bool tbz::TbZUtility::isTightElectron( edm::Ptr<reco::GsfElectron> elec)
{
    bool isTight= false;
    // std::cout<<"WARNING TIGHT  ELECTRON CUTS NOT IMPLEMENTED YET!!!!"<<std::endl;
    // if(elec-> isElectron())
    // {
    // if(fabs(elec->gsfTrack()->dxy(vertex_->position())) < 0.04) //fabs(electron.gsfTrack()->dxy(vertex_->position()))
      // {
      // if(elec -> passConversionVeto() = true)
         // {
         // if(elec -> electronID("mvaTrigV0") > 0.5 )
            // {
            // if(elec -> gsfTrack()->trackerExpectedHitsInner().numberOfHits() <=0)
            // {
         // isTight = true;
           // }
         // }
      // }
   // } 
  // }   
  return isTight;
}
bool tbz::TbZUtility::isMediumElectron(edm::Ptr<reco::Electron> )
{
    bool isMedium=true;
    std::cout<<"WARNING MEDIUM  ELECTRON CUTS NOT IMPLEMENTED YET!!!!"<<std::endl;
    return isMedium;
}

bool tbz::TbZUtility::isLooseElectron(edm::Ptr<reco::Electron> )
{
    bool isLoose=true;
    std::cout<<"WARNING LOOSE  ELECTRON CUTS NOT IMPLEMENTED YET!!!!"<<std::endl;
    return isLoose;
}

bool tbz::TbZUtility::sfosMassInRange(Double_t low, Double_t high, std::vector<double> &_SFOS_masses)
{
	for(unsigned int i = 0; i < _SFOS_masses.size(); ++i)
     {
     if(_SFOS_masses.at(i) >= low && _SFOS_masses.at(i) <= high) return true;
     }
    return false;
     }//end of Muon sfosMassInRange
      //.............................brackets are OK up to this........

void  tbz::TbZUtility::makePairs( const std::vector<reco::Muon> &mv, 
                                        std::vector<double>& sfos_masses,
                                        std::pair < int,int>& minM_indexPair,
                                        std::vector<reco::NamedCompositeCandidate >& DY)
{
	sfos_masses.clear();
	double zmass = 91.1876;
	double min=1000.;

	int j=0;
	int jj=0;

	int m_i =100;
	int m_ii=100;
	double hard_pt =0.;

	for(std::vector<reco::Muon>::const_iterator muItr = mv.begin(); muItr != mv.end(); muItr++)
	{
	jj=0;
	for(std::vector<reco::Muon>::const_iterator muItr2 = mv.begin(); muItr2 != mv.end(); muItr2++)
	{
	std::cout<<"charge1 "<<(*muItr).charge()<<"  charge2 "<<(*muItr2).charge() <<std::endl;

	if(muItr2 > muItr && ( (*muItr).charge()*(*muItr2).charge() ) <0 )
		{
	double dphi = (*muItr).phi() - (*muItr2).phi();
	if( fabs(dphi)< 0.5) continue; // 0.1 previous
	double mass  =  ( (*muItr).p4() + (*muItr2).p4() ).mass();
	reco::NamedCompositeCandidate DY_tmp;
	DY_tmp.setP4( (*muItr).p4() + (*muItr2).p4() );
	//DY.push_back(DY_tmp);
	std::cout<<"charge1 "<<(*muItr).charge()<<"  charge2 "<<(*muItr2).charge()<< " Z mass "<< mass << " dphi "<< dphi<<std::endl;
	//sfos_masses.push_back ( mass );

	if( (*muItr).pt() > (*muItr2).pt() )
	hard_pt =  (*muItr).pt();
	else  {
	hard_pt =  (*muItr2).pt();
	}

	double m =  zmass - mass;
	if(fabs(m) <min  && hard_pt > 20.) { 
	  min = fabs(m); m_i = j; m_ii = jj; 

	  if(DY.size()) DY.clear();
	  DY.push_back(DY_tmp);
	  sfos_masses.push_back ( mass );

	}

		}//end of SFOS anti-matching  

	jj++;
	}//end of second loop
	
	j++;
	
	}//END OF MV LOOP

	minM_indexPair= std::make_pair(m_i, m_ii);

	std::cout<<"closest Z mass "<<min<<"  with muon indices : ( "<<m_i<<",   "<<m_ii<<" )"<< " cont size : "<<  mv.size()<<std::endl; 

}//end of rePair

//-------------------------------------------------------------- 
void  tbz::TbZUtility::makeEPairs(const std::vector<reco::GsfElectron> &ev, std::vector<double>& el_sfos_masses,  std::pair <int, int>& minE_indexPair,  std::vector<reco::NamedCompositeCandidate >& DY1 )
{
        el_sfos_masses.clear();
        double zmass1 = 91.1876;
        double min1=1000.;


        int e=0;
        int ee=0;

        int m_e =100;
        int m_ee=100;
	double hard_pt =0.;

        for(std::vector<reco::GsfElectron>::const_iterator electItr = ev.begin();
                                                           electItr != ev.end() ;
                                                           electItr++            )
        {
        ee=0;
        for(std::vector<reco::GsfElectron>::const_iterator electItr2 = ev.begin();
                                                           electItr2 != ev.end() ;
                                                           electItr2++            )
        {

	//if(electItr2 <= electItr)continue;

	double dphi = (*electItr).phi() - (*electItr2).phi();

	std::cout<<"charge1 "<<(*electItr).charge()<<"  charge2 "<<(*electItr2).charge()<< " dphi "<< dphi<<std::endl;


        if(electItr2 > electItr && ( (*electItr).charge()*(*electItr2).charge() ) <0 )
                {
        //double dphi = (*electItr).phi() - (*electItr2).phi();
        if( fabs(dphi)< 0.5) continue; //previous 1.
	double mass1  =  ( (*electItr).p4() + (*electItr2).p4() ).mass();
        reco::NamedCompositeCandidate DY_tmp1;
        DY_tmp1.setP4( (*electItr).p4() + (*electItr2).p4() );
	//      DY1.push_back(DY_tmp1);
        std::cout<< " Z mass "<< mass1 << " dphi "<< dphi<<std::endl;
        //el_sfos_masses.push_back ( mass1 );
        
	if( (*electItr).pt() > (*electItr2).pt() )
	  hard_pt =  (*electItr).pt();
        else {
	  hard_pt =  (*electItr2).pt();
        }


        double m1 =  zmass1 - mass1;
        if(fabs(m1) <min1  && hard_pt >20.) { 
	  min1 = fabs(m1); m_e = e; m_ee = ee; 

	  if(DY1.size())DY1.clear();
	  DY1.push_back(DY_tmp1);
	  el_sfos_masses.push_back ( mass1 );
	}

                }//end of SFOS anti-matching        

	ee++;
        }//end of second loop

        e++;

        }//END OF MV LOOP    


	minE_indexPair= std::make_pair(m_e, m_ee);        
  cout<<"closest Z mass "<<min1<<"  with Electron indices : ( ";
  cout<<m_e<<",   "<<m_ee<<" )"<< " cont size : "<<  ev.size()<<endl;  

}//end of electron rePair


//------- MC Truth Matching ----- from here -----
 // if(doTruthMatch_)
// {
void  tbz::TbZUtility::isTrue(edm::Ptr<reco::GsfElectron> &e_in
                             ,edm::Handle< std::vector<reco::GenParticle> > &genParticles
                             ,std::vector<reco::GenParticle> &GenElectron
                             ,double &elec_z_mass, double &elec_eta, double& elec_phi
                             ,double &elec_pt_ratio)
{        
    using namespace edm;
    using namespace HepMC;
    using namespace reco;
    using namespace std;  
    
    std::vector<reco::GenParticle>::const_iterator    j;     
    for(j=genParticles->begin(); j != genParticles->end(); ++j)
     {
        elec_pt_ratio =0.;
        elec_z_mass=0.;
        elec_eta =-1000.;
        elec_phi= -1000.;        
        const Candidate * d1 = (*j).mother();
        if ((*j).pdgId() ==11 || (*j).pdgId() == -11 && d1->pdgId()==23)
        //if(d1->pdgId()==23)
        {
         double deta_elec = (*j).eta() - e_in->eta();
         double dphi_elec = (*j).phi() - e_in->phi();
         double dR_elec = sqrt(deta_elec*deta_elec + dphi_elec*dphi_elec);
         cout<< " dR1_elec @@@@@@@@@@@@@@@@@ *******************: "<<dR_elec<<endl;
         GenElectron.push_back((*j));
        } 
     }
}
//===========================================================================================
//--------------25-01-14
void  tbz::TbZUtility::isTrueWenu(edm::Ptr<reco::GsfElectron> &e_in
                             ,edm::Handle< std::vector<reco::GenParticle> > &genParticles
                             ,std::vector<reco::GenParticle> &GenElectron_W
                             /*,std::vector<reco::GenParticle> &GenNuetrino_W*/
                             ,double &w_enu_mass, double &wenu_eta, double& wenu_phi
                             ,double &pt_ratio_elecW)
{
  using namespace edm;
  using namespace HepMC;
  using namespace reco;
  using namespace std;
  
  std::vector<reco::GenParticle>::const_iterator k;
  for(k=genParticles->begin();k!=genParticles->end();++k)
    {
        pt_ratio_elecW =0.                    ;
        w_enu_mass=0.                         ;
        wenu_eta =-1000.                      ;
        wenu_phi= -1000.                      ;
        const Candidate * p = (*k).mother()   ;
        
        if(abs((*k).pdgId()) == 11 || abs((*k).pdgId()) == 12   /*either MU or MU-Neutrino*/
                && abs(p->pdgId() ) == 24                      /*Mother must be W*/
                &&(p->mother()->pdgId()) == 6 )               /*Grand-Mother must be t*/
        {
        double deta_elec_W = (*k).eta() - e_in->eta()                                 ;
        double dphi_elec_W = (*k).phi() - e_in->phi()                                 ;
        double dR_elecW   = sqrt(deta_elec_W*deta_elec_W + dphi_elec_W*dphi_elec_W)   ;
        cout<< " dR1_mu_W @@@@@@@@@@@ **********: "<<dR_elecW<<endl                   ;
        GenElectron_W.push_back((*k))                                                 ;     
        }
       if(abs((*k).pdgId()) == 5 && abs(p->pdgId() ) == 6)
       {
       GenElectron_W.push_back((*k));
       }
    }
}
     
//=========================================================================================== 
   
void tbz::TbZUtility::doTruthMatching( edm::Ptr<reco::Muon> &muon  
                                      ,edm::Handle< std::vector<reco::GenParticle> > &EvtHandle 
                                      ,std::vector<reco::GenParticle> &MuonGenPart
                                      /*,std::vector<reco::GenParticle> &MuonGenPart_W*/
                                      ,double &mu_z_mass, double &eta, double& phi, double &pt_ratio)
{
  using namespace edm    ;
  using namespace HepMC  ;
  using namespace reco   ;
  using namespace std    ;
  
  std::vector<reco::GenParticle>::const_iterator i;
  for(i=EvtHandle->begin();i!=EvtHandle->end();++i)
    {
        pt_ratio  =  0.                      ;
        mu_z_mass =  0.                      ;
        eta       = -1000.                   ;
        phi       = -1000.                   ;
        const Candidate * d = (*i).mother()  ;
        
        if ((*i).pdgId()==13 || (*i).pdgId()==-13)
        if (d->pdgId()==23 ) 
         {
        double deta = (*i).eta() - muon->eta()              ;
        double dphi = (*i).phi() - muon->phi()              ;
        double dR = sqrt(deta*deta+dphi*dphi )              ;
        cout<< " dR1 @@@@@@@@@@@@@@***********: "<<dR<<endl ;
        MuonGenPart.push_back((*i))                         ;   
        
       }             
    }
}
//--------------24-01-14
void tbz::TbZUtility::doWTruthMatching( edm::Ptr<reco::Muon> &muon  
                                      ,edm::Handle< std::vector<reco::GenParticle> > &EvtHandle 
                                      ,std::vector<reco::GenParticle> &MuonGenPart_W
                                      ,double &mu_z_mass, double &eta, double& phi, double &pt_ratio)
{
  using namespace edm       ;
  using namespace HepMC     ;
  using namespace reco      ;
  using namespace std       ;
  
  std::vector<reco::GenParticle>::const_iterator J;
  for(J=EvtHandle->begin();J!=EvtHandle->end();++J)
    {
        pt_ratio   =  0.    ;
        mu_z_mass  =  0.    ;
        eta        = -1000. ;
        phi        = -1000. ;
        
        const Candidate * d = (*J).mother();
        if( abs((*J).pdgId()) == 13 || abs((*J).pdgId()) == 14 /*either MU or MU-Neutrino*/
                && d->pdgId() == 24                           /*Mother must be W*/
                && (d->mother()->pdgId())==6 )               /*Grand-Mother must be t*/
           {
        double deta_mu_W = (*J).eta() - muon->eta();
        double dphi_mu_W = (*J).phi() - muon->phi();
        double dR_muW    = sqrt(deta_mu_W*deta_mu_W + dphi_mu_W*dphi_mu_W);
        cout<< " dR1_mu_W @@@@@@@@@@@ **********: "<<dR_muW<<endl;
        MuonGenPart_W.push_back((*J));          
          }
        if(abs((*J).pdgId()) == 5 && abs(d->pdgId()) == 6 )
        {
        MuonGenPart_W.push_back((*J));
        }
    }
}
//--150714--
  using namespace edm       ;
  using namespace HepMC     ;
  using namespace reco      ;
  using namespace std       ;

 void tbz::TbZUtility::setNLooseMu(int &nLM)
 {
 nLooseMu = nLM;
 }
 int tbz::TbZUtility::getNLooseMu()
 {
 return nLooseMu;
 }
 void tbz::TbZUtility::setNLooseElec(int &nLE)
 {
 nLooseElec = nLE;
 }
 int tbz::TbZUtility::getNLooseElec()
{
return nLooseElec;
} 
 
 //----



