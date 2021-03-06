import FWCore.ParameterSet.Config as cms

process = cms.Process("DEMO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('RecoBTag/Configuration/RecoBTag_cff')
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
#JetMETCorrections
#process.load("RecoMET.METFilters.metFilters_cff")
#process.load("JetMETCorrections.Type1MET.caloMETCorrections_cff")
#process.load("JetMETCorrections.Type1MET.correctionTermsCaloMet_cff")
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff")
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0RecoTrack_cff")
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")
#process.load("JetMETCorrections.Type1MET.correctedMet_cff")

#---------12-05-14-------
#process.load("CommonTools.ParticleFlow.pfIsolation_cfg")
#from CommonTools.ParticleFlow import pfIsolation_cfg

########################################
# Summer12 with CMSSW_5_3_X and GT START53_V7A
# DoubleMuon trigger 
# HLT_Mu17_Mu8_v17 or HLT_Mu17_TkMu8_v10 
# DoubleElectron trigger 
# HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17
########################################

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #DATA2012A6
    #DATA2012A7
    #DATA2012A8
    #DATA2012A9
    #Data2101A

    #SignalMC
    #TbZSignalMC
    #TbZSignalMC
    #Data2012A
    #ZZ_MC
    #WW_2012MC8TeV
    #TTbarMC8TeV
    #ZMC_2012
    #WZ_2012MC
    #WJETS_MC
   
    #SignalMC
   'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/A27919CE-486A-E311-8EDD-002590D0B0C8.root'

   #2012_EE
   # 'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/58C97ECA-45DA-E111-B7B9-008CFA001B18.root' 
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/90985F3F-11DA-E111-8278-848F69FD294F.root'

    #2012_DataSet
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/F4FBA42A-A6CF-E111-BAAB-003048678A7E.root'

    #DoublMu
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/AC00CD43-B9D2-E111-97FD-001A92971B20.root'
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/527B7FC8-19CF-E111-B116-00304867BEC0.root'

    #MuEG
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/0A47AE1F-BFD9-E111-8AB4-20CF3019DF0C.root'
    
    #WZ_MC
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/38576F66-61DD-E111-A13B-002481E0DC4E.root'	
    ##'file:/afs/cern.ch/user/m/mwaqas/TBZ14/TBZMET/CMSSW_5_3_15/src/MyAnalysis/TbZ/001C0B68-536A-E311-B25F-002590D0B066.root'
    #'file:/afs/cern.ch/work/m/meshoaib/analysis/CMSSW_5_3_15/src/MyAnalysis/TbZ/001C0B68-536A-E311-B25F-002590D0B066.root'
    #'file:/afs/cern.ch/work/m/meshoaib/Analysis_TBZ/CMSSW_5_3_7/src/FirstAnalysis/TBZAnalysis/A00B400A-B1CF-E111-BCEF-0026189438C2.root'
    #'file:/afs/cern.ch/work/m/meshoaib/Analysis_TBZ/CMSSW_5_3_7/src/FirstAnalysis/TBZAnalysis/441B933C-3846-E211-9928-00266CFFA768.root'
    #'file:/ehep/ehep_data/datasets/DoubleMu_2012RunA/13Jul2012-v1/AC00CD43-B9D2-E111-97FD-001A92971B20.root'    
    #'file:/ehep/ehep_data/datasets/MC_2012/TT/0C8518DA-42F1-E111-BD2F-003048D3CB3C.root'

    )
)

#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')## To delete duplicate events, this included checkAllFilesOpened

process.source.duplicateCheckMode = cms.untracked.string('checkAllFilesOpened')

process.TFileService = cms.Service("TFileService",
        fileName= cms.string("DoubleElectronB.root")
)

print "before" 

process.myProducerLabel = cms.EDProducer('TbzMuonProducer',
                                         muonTag = cms.InputTag("muons"),
                                         muonPtCut = cms.double(20.)    ,
                                         muonEtaCut = cms.double(2.1)                                         
)                                                                                  

process.myLoseMuons = cms.EDProducer('TbzLoseMuonProducer',
                                         muonTag = cms.InputTag("muons"),
                                         muonPtCut = cms.double(10.)    ,
                                         muonEtaCut = cms.double(2.1)                                         
)                                                                 

process.myJetProdLabel = cms.EDProducer('JetProducer',
                                         JetsPtCut  = cms.double(30.), #as suggested by jets contact person
                                         JetsEtaCut = cms.double(3.0)
                                        )

from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

# particle flow isolation
#

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso

process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons') 
#print "NewIsolation"
#process.tightElecSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence )
#process.Newpfiso = cms.Sequence(process.tightElecSequence )

process.myElectronProdLabel = cms.EDProducer('ElectronProducer'  ,
            electronsInputTag = cms.InputTag("gsfElectrons")     ,
            muonsInputTag = cms.InputTag("myLoseMuons")      ,
            conversionsInputTag = cms.InputTag("allConversions") ,
            beamSpotInputTag  = cms.InputTag("offlineBeamSpot")  ,
            rhoIsoInputTag   = cms.InputTag("kt6PFJetsForIsolation","rho"),
            primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),

            isoValInputTags =  cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                             cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                              cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),	

                                             electronPtCut = cms.double(20.),
                                             electronEtaCut = cms.double(2.1),
                                             printDebug              = cms.bool(True)
)

print "after" 


print "NewIsolation"

#from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
#process.tightElecSequence = setupPFElectronIso(process, 'gsfElectrons') # reco::GsfElectron
#print "Hello Iso"
process.tightElecSequence = setupPFElectronIso(process, 'myElectronProdLabel') #electronprod
print "Hello Iso"
process.Newpfiso = cms.Sequence(process.tightElecSequence )

print "After New ISolation"


#process.metFilters
process.topAna = cms.EDAnalyzer('TbZTopAnalyzer',

MutriggerSelection   = cms.string('(HLT_Mu17_Mu8_* OR HLT_Mu17_TkMu8_*)'),
MuEGtriggerSelection = cms.string('(HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_* OR HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),
EtriggerSelection    = cms.string('(HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_*)'),

triggerConfiguration =  cms.PSet(
  hltResults = cms.InputTag('TriggerResults','','HLT'),
  l1tResults = cms.InputTag(''),
  daqPartitions = cms.uint32(1),
  l1tIgnoreMask = cms.bool( False ),
  l1techIgnorePrescales = cms.bool( False ),
  throw  = cms.bool( True )
),

                                electronsInputTag     = cms.InputTag("myElectronProdLabel")                        ,
                                #electronsInputTag    = cms.InputTag("gsfElectrons")                               ,
                                muonsInputTag         = cms.InputTag("myProducerLabel")                            , 
                             #  muonsInputTag         = cms.InputTag("muons")                                      ,
                                jetsInputTag          = cms.InputTag("ak5PFJets")                                  ,#ak5PFJetsak5CaloJets
                                bjetInputTag          = cms.InputTag("combinedSecondaryVertexBJetTags")            ,
                                metInputTag           = cms.InputTag("pfMet")                                      ,
                                conversionsInputTag   = cms.InputTag("allConversions")                             ,
                                beamSpotInputTag      = cms.InputTag("offlineBeamSpot")                            ,
                                rhoIsoInputTag        = cms.InputTag("kt6PFJetsForIsolation","rho")                ,

				#rhoIsoInputTagAna  = cms.InputTag("kt6PFJetsForIsolation","rho")                   ,

                                primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices")                     ,
                               
				NewisoValInputTags =  cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso') ,
                               					 cms.InputTag('elPFIsoValueGamma03PFIdPFIso')      ,
            					                 cms.InputTag('elPFIsoValueNeutral03PFIdPFIso'))   ,                               
                                
                                metPtCut      = cms.double(25.)                                    ,
				vertexSrc     = cms.InputTag("offlinePrimaryVertices")             ,
                                BtagEtaCut    = cms.double(2.4)                                    ,
                                BtagPtCut     = cms.double(25)                                     ,
                                BtagDiscrCut  = cms.double(0.679)                                  ,
                                DPHiENue      = cms.double(1.0)                                    ,
                                DPHiMuNue     = cms.double(1.0)                                    ,
                                MaxZMass      = cms.double(102.0)                                  ,
                                MinZMAss      = cms.double(78.0)                                   ,
                                JetsPtCut     = cms.double(20.)                                    ,
                                ElecPtCut     = cms.double(20.)                                    ,
                                doTruthMatch  = cms.bool(False)                                    ,
				electronPtCut = cms.double(20.)                                    ,
                                electronEtaCut = cms.double(2.1)                                   ,

                                #realdata      = cms.bool(True)                                    ,
				realdata      = cms.bool(False)                                    ,
				doPileup      = cms.bool(False)                                    ,	
                                printDebug    = cms.bool(True)

                              )
print "after 1" 
process.dump=cms.EDAnalyzer('EventContentAnalyzer')
print "after 2" 

process.p = cms.Path(process.myProducerLabel*process.myLoseMuons*process.myJetProdLabel*process.kt6PFJetsForIsolation*process.pfiso*process.myElectronProdLabel*process.Newpfiso*process.topAna)


print "after 3" 
