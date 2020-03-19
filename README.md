# RivetAnalyses
Repo for private and WORKINPROGRESS Rivetanalyses 
The official CMS proposals on implementing a Rivet analysis based on your analysis can be found [here](https://twiki.cern.ch/twiki/bin/viewauth/CMS/Rivet "TWIKI:Rivet").

## Usage
### Standalone Rivet
For standalone usage of the Rivet analyses on hepMC files see [Getting Started](https://rivet.hepforge.org/trac/wiki/GettingStarted "Getting Started with Rivet") on the official Rivet website.
### Running with CMSSW
It is alternatively possible to run Rivet analyses on CMS GENSIM, FULLSIM, AODSIM, MINIAODSIM and NONOAODSIM using the [GenParticle2HepMCConverter](https://twiki.cern.ch/twiki/bin/view/CMS/GenParticles2HepMCConverter "TWIKI: GenParticle2HepMCConverter"). Be aware that the available information in the different SIM tiers differs.

To run a Rivet analysis on GENSIM, FULLSIM and AODSIM follow the instructions given [here](https://twiki.cern.ch/twiki/bin/view/CMS/RivetontoAODSIM "TWIKI:RivetontoAODSIM") and run CMSSW with your adapted configuration file. The compatibility has been tested for CMSSSW_10_6_0 and Rivet 2.7.0 running on GENSIM and AODSIM files.

For MINIAODSIM and NANOAODSIM additional EDProducers have to be interposed as described [here](https://twiki.cern.ch/twiki/bin/view/CMS/ParticleLevelProducer "TWIKI:ParticleLevelProducer").

You can find example configurations for the presented analyses in the `CMSSWconfigs` directory. 
