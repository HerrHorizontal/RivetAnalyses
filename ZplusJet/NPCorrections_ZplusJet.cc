// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class NPCorrections_ZplusJet : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(NPCorrections_ZplusJet);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState ffs(Cuts::abseta < 5 && Cuts::pT > 100*MeV);
      const ChargedFinalState cfs(ffs);
      declare(ffs, "FS");
      declare(cfs, "CFS");
      // leptons
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      declare(bare_leps, "bare_leptons");
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 10*GeV);
      declare(dressed_leps, "dressed_leptons");
      // jet collection
      FastJets jets(ffs, FastJets::ANTIKT, 0.4, JetAlg::Invisibles::NONE)
      declare(jets, "Jets");

      // out of acceptance particles treat as invisible
      VetoedFinalState fs_onlyinacc(ffs, (Cuts::abspid::MUON && Cuts::abseta > 2.4) || 
                                    (Cuts::abspid::PHOTON && Cuts::abseta > 3.0) || 
                                    (Cuts::abspid::ELECTRON && Cuts::abseta > 3.0));
      declare(MissingMomentum(fs_onlyinacc), "invisibles");

      // Book histograms
      //// Book histograms with equidistant bins
      _hist_NMus = bookHisto1D("NMuons", 10, 0, 10);
      _hist_MuPlusPt = bookHisto1D("MuPlusPt", 11, 25, 300);
      _hist_MuPlusEta = bookHisto1D("MuPlusEta", 24, -2.4, 2.4);
      _hist_MuPlusPhi = bookHisto1D("MuPlusPhi", 30, -3.14, 3.14);
      _hist_MuMinusPt = bookHisto1D("MuMinusPt", 11, 25, 300);
      _hist_MuMinusEta = bookHisto1D("MuMinusEta", 24, -2.4, 2.4);
      _hist_MuMinusPhi = bookHisto1D("MuMinusPhi", 30, -3.14, 3.14);
      _hist_ZY = bookHisto1D("ZY", 12, -2.4, 2.4);
      _hist_ZM = bookHisto1D("ZM", 20, 71, 111);
      _hist_ZPhi = bookHisto1D("ZPhi", 30, -3.14159, 3.14159);
      _hist_NJets = bookHisto1D("NJets", 10, 0, 10);
      _hist_Jet1Y = bookHisto1D("Jet1Y", 12, -2.4, 2.4);
      _hist_Jet1Phi = bookHisto1D("Jet1Phi", 31, -3.142, 3.142);
      _hist_Jet2Y = bookHisto1D("Jet2Y", 30, -2.4, 2.4);
      _hist_Jet2Phi = bookHisto1D("Jet2Phi", 31, -3.142, 3.142);
      _hist_Jet3Y = bookHisto1D("Jet3Y", 30, -2.4, 2.4);
      _hist_Jet3Phi = bookHisto1D("Jet1Phi", 31, -3.142, 3.142);
      _hist_MET = bookHisto1D("MET", 70, 0, 350);

      //// Book histograms with variable bin size
      ////// Define bin edges
      std::vector<double> binedges_ZPt = {25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 220, 250, 400, 1000};

      std::vector<double> binedges_Jet1Pt = {20, 25, 30, 40, 50, 75, 125, 175, 225, 300, 400, 500, 700, 1000};
      std::vector<double> binedges_Jet2Pt = {5, 10, 20, 30, 40, 50, 75, 125, 175, 250, 400};
      std::vector<double> binedges_Jet3Pt = {5, 10, 20, 30, 40, 50, 75, 125, 175, 250, 400};
      std::vector<double> binedges_JetAvePt = {5, 10, 15, 20, 25, 30, 40, 50, 75, 125, 175, 225, 300, 400};
      std::vector<double> binedges_Jet1Eta = {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -1.930,
                                             -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0.000, 0.261, 0.522, 0.783, 1.044, 
                                             1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 
                                             5.191};
      ////// Book Histograms with vector
      _hist_ZPt = bookHisto1D("ZPt", binedges_ZPt);
      _hist_Jet1Pt = bookHisto1D("ZPt", binedges_Jet1Pt);
      _hist_Jet2Pt = bookHisto1D("ZPt", binedges_Jet2Pt);
      _hist_Jet3Pt = bookHisto1D("ZPt", binedges_Jet3Pt);
      _hist_JetAvePt = bookHisto1D("ZPt", binedges_JetAvePt);
      _hist_Jet1Eta = bookHisto1D("ZPt", binedges_Jet1Eta);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here
      //// Find muons and jets
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "dressed_leptons").dressedLeptons(Cuts::abspid::MUON && Cuts::abseta < 2.4 && Cuts::pT > 25*GeV);
      Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::absrap < 2.4 && Cuts::pT > 0.1*GeV);

      //// Require at least two opposite sign muons compatible with Z-boson mass
      if (muons.size() < 2) vetoEvent;
      for () {
        //// Require at least one opposite charge pair
        if (muons[i].pid() != -muons[j].pid()) vetoEvent;
        //// Require invariant mass of both muons to be compatibale with Z-mass
        if (!(71.1876*GeV <= (muons[0].mom()+muons[1].mom()).mass() <= 111.1876*GeV)) vetoEvent;
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_YYYY); // normalize to unity
      scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _hist_NMus, _hist_MuPlusPt, _hist_MuPlusEta, _hist_MuPlusPhi, _hist_MuMinusPt, _hist_MuMinusEta, _hist_MuMinusPhi,
              _hist_ZY, _hist_ZM, _hist_ZPhi, _hist_NJets, _hist_Jet1Y, _hist_Jet1Phi, _hist_Jet2Y, _hist_Jet2Phi, _hist_Jet3Y, _hist_Jet3Phi, _hist_MET;
    Histo1DPtr _hist_ZPt, _hist_Jet1Pt, _hist_Jet2Pt, _hist_Jet3Pt, _hist_JetAvePt, _hist_Jet1Eta;

    Profile1DPtr _p_AAAA;
    CounterPtr _c_BBBB;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(NPCorrections_ZplusJet);


}
