import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("runRivetAnalysis")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring('file:/portal/ekpbms1/home/mschnepf/QCD/NP_corrections/runtime/CMSSW_10_6_2/src/RivetAnalyses/ZplusJet/MCSim_1.root') )
#process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring('file:/storage/gridka-nrg/mschnepf/MC_Production_official/CreateAOD/AODSIM_0.root') )
#process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring('file:/portal/ekpbms1/home/mschnepf/QCD/NP_corrections/runtime/CMSSW_10_6_0/src/Rivet/NP_Correction/6ECFB7E5-B9EB-E611-B128-0CC47A4D7662.root') )
#process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/5C962DEB-DEEB-E611-A148-FA163E4700A9.root') )
process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring(os.environ['FILE_NAMES']) )

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("genParticles"),
    genEventInfo = cms.InputTag("generator"),
    signalParticlePdgIds = cms.vint32(),
)
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.AnalysisNames = cms.vstring('ZplusJet_Partonic')
process.rivetAnalyzer.OutputFile ='out.yoda'

process.p = cms.Path(process.generator*process.rivetAnalyzer)
