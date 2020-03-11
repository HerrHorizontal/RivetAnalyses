import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("runRivetAnalysis")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

#process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring('file:/portal/ekpbms1/home/mschnepf/QCD/NP_corrections/runtime/CMSSW_10_6_2/src/RivetAnalyses/ZplusJet/MCSim_1.root') )
#process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring('file:/storage/gridka-nrg/mschnepf/MC_Production_official/CreateAOD/AODSIM_0.root') )
#process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring('file:/portal/ekpbms1/home/mschnepf/QCD/NP_corrections/runtime/CMSSW_10_6_0/src/Rivet/NP_Correction/6ECFB7E5-B9EB-E611-B128-0CC47A4D7662.root') )
process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring(os.environ['FILE_NAMES']) )

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("mergedGenParticles"),
    genEventInfo = cms.InputTag("generator"),
    signalParticlePdgIds = cms.vint32(),
)

process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi")
process.particleLevel.HepMCCollection = cms.InputTag("genParticles2HepMC:unsmeared")

process.path = cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel)


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
#process.rivetAnalyzer.AnalysisNames = cms.vstring('NPCorrections_ZplusJet')
process.rivetAnalyzer.AnalysisNames = cms.vstring('ZplusJet')
process.rivetAnalyzer.OutputFile ='out.yoda'

process.p = cms.Path(process.generator*process.rivetAnalyzer)
