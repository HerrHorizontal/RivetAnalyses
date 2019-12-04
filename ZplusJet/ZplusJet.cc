// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Projections/MissingMomentum.hh"

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

      /// Initialise and register projections
      const FinalState ffs(Cuts::abseta < 5 && Cuts::pT > 100*MeV);
      const ChargedFinalState cfs(ffs);
      declare(ffs, "FS");
      declare(cfs, "CFS");
      /// leptons
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      declare(bare_leps, "bare_leptons");
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 25*GeV);
      declare(dressed_leps, "dressed_leptons");
      /// jet collection
      FastJets jets(ffs, FastJets::ANTIKT, 0.4);
      declare(jets, "Jets");

      /// out of acceptance particles treat as invisible
      VetoedFinalState fs_onlyinacc(ffs, (Cuts::abspid ==  PID::MUON && Cuts::abseta > 2.4) || 
                                    (Cuts::abspid == PID::PHOTON && Cuts::abseta > 3.0) || 
                                    (Cuts::abspid == PID::ELECTRON && Cuts::abseta > 3.0));
      declare(MissingMomentum(fs_onlyinacc), "invisibles");

      /// Book histograms
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
//      book(_h["NMus"], "N_Muons", 10, 0, 10);
//      book(_h["MuPlusPt"], "pT_Mu+", 11, 25, 300);
//      book(_h["MuPlusEta"], "eta_Mu+", 24, -2.4, 2.4);
//      book(_h["MuPlusPhi"], "phi_Mu+", 30, -3.14, 3.14);
//      book(_h("MuMinusPt"), "pT_Mu-", 11, 25, 300);
//      book(_h["MuMinusEta"], "eta_Mu-", 24, -2.4, 2.4);
//      book(_h["MuMinusPhi"], "phi_Mu-", 30, -3.14, 3.14);
//      book(_h["ZY", "y_Z", 12, -2.4, 2.4);
//      book(_h["ZM"], "m_Z", 20, 71, 111);
//      book(_h["ZPhi"], "phi_Z", 30, -3.14159, 3.14159);
//      book(_h["NJets"], "N_Jets", 10, 0, 10);
//      book(_h["Jet1Y"], "y_Jet1", 12, -2.4, 2.4);
//      book(_h["Jet1Phi"], "phi_Jet1", 31, -3.142, 3.142);
//      book(_h["Jet2Y"], "y_Jet2", 30, -2.4, 2.4);
//      book(_h["Jet2Phi"], "phi_Jet2", 31, -3.142, 3.142);
//      book(_h["Jet3Y"], "y_Jet3", 30, -2.4, 2.4);
//      book(_h["Jet3Phi"], "phi_Jet3", 31, -3.142, 3.142);
//      book(_h["MET"], "MET", 70, 0, 350);

      //// Book histograms with variable bin size
      ///// Define bin edges
      vector<double> binedges_ZPt = {25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 220, 250, 400, 1000};

      vector<double> binedges_Jet1Pt = {20, 25, 30, 40, 50, 75, 125, 175, 225, 300, 400, 500, 700, 1000};
      vector<double> binedges_Jet2Pt = {5, 10, 20, 30, 40, 50, 75, 125, 175, 250, 400};
      vector<double> binedges_Jet3Pt = {5, 10, 20, 30, 40, 50, 75, 125, 175, 250, 400};
      vector<double> binedges_JetAvePt = {5, 10, 15, 20, 25, 30, 40, 50, 75, 125, 175, 225, 300, 400};
      vector<double> binedges_Jet1Eta = {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -1.930,
                                             -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0.000, 0.261, 0.522, 0.783, 1.044, 
                                             1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 
                                             5.191};
      vector<double> binedges_PhiStarEta = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50};
      ///// Book Histograms with vector
      _hist_ZPt = bookHisto1D("ZPt", binedges_ZPt);
      _hist_Jet1Pt = bookHisto1D("Jet1Pt", binedges_Jet1Pt);
      _hist_Jet2Pt = bookHisto1D("Jet2Pt", binedges_Jet2Pt);
      _hist_Jet3Pt = bookHisto1D("Jet3Pt", binedges_Jet3Pt);
      _hist_JetAvePt = bookHisto1D("JetAvePt", binedges_JetAvePt);
      _hist_Jet1Eta = bookHisto1D("Jet1Eta", binedges_Jet1Eta);
//      book(_h["ZPt"], "pT_Z", binedges_ZPt);
//      book(_h["Jet1Pt"], "pZ_Jet1", binedges_Jet1Pt);
//      book(_h["Jet2Pt"], "pT_Jet2", binedges_Jet2Pt);
//      book(_h["Jet3Pt"], "pT_Jet3", binedges_Jet3Pt);
//      book(_h["JetAvePt"], "<pT>_Jets", binedges_JetAvePt);
//      book(_h["Jet1Eta"], "eta_Jet1", binedges_Jet1Eta);
//      book(_h["PhiStarEta"], "Phi*_eta", binedges_PhiStarEta);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// Find muons and jets
      //vector<DressedLepton> muons = apply<DressedLeptons>(event, "dressed_leptons").dressedLeptons(Cuts::abspid == PID::MUON && Cuts::abseta < 2.4 && Cuts::pT > 25*GeV);
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "dressed_leptons").dressedLeptons();

      Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::absrap < 2.4 && Cuts::pT > 0.1*GeV);

      /// Remove jet-muon overlap
      idiscardIfAnyDeltaRLess(jets, muons, 0.3);

      /// Require at least one hard jet
      if (jets.size()<1) vetoEvent;      
      if (jets.at(0).pT() <= 20*GeV) vetoEvent;

      /// Require at least two opposite sign muons compatible with Z-boson mass and keep the pair closest to Zboson mass
      bool _bosoncandidateexists = false;
      double _massdiff = 20*GeV;
      DressedLepton _muon = muons.at(0);
      DressedLepton _antimuon = muons.at(0);

      if (muons.size() < 2) vetoEvent;

      for (unsigned int it = 1; it < muons.size(); ++it) {
        for (unsigned int jt = 0; jt < it; ++jt) {
          double _candidatemass = (muons.at(it).mom() + muons.at(jt).mom()).mass();
          if (muons.at(it).pid() == -muons.at(jt).pid() && abs(_candidatemass - 91.1876*GeV) < _massdiff) {
            _bosoncandidateexists = true;
            _massdiff = abs(_candidatemass - 91.1876*GeV);
            if (muons.at(it).pid() > 0) {
              _muon = muons.at(it);
              _antimuon = muons.at(jt);
            }
            else {
              _muon = muons.at(jt);
              _antimuon = muons.at(it);
            }
          }
          else continue;
        }
      }

      if (!(_bosoncandidateexists)) vetoEvent;

      /// Fill muon related histograms
      _hist_NMus -> fill(muons.size());
      //_h["NMus"] -> fill(muons.size());

      const double phi_Mu = _muon.phi();
      _hist_MuMinusPhi -> fill(phi_Mu);
      //_h["MuMinusPhi"] -> fill(phi_Mu);
      const double phi_AntiMu = _antimuon.phi();
      _hist_MuPlusPhi -> fill(phi_AntiMu);
      //_h["MuPlusPhi"] -> fill(phi_AntiMu);

      const double eta_Mu = _muon.eta();
      _hist_MuMinusEta -> fill(eta_Mu);
      //_h["MuMinusEta"] -> fill(eta_Mu);
      const double eta_AntiMu = _antimuon.eta();
      _hist_MuPlusEta -> fill(eta_AntiMu);
      //_h["MuPlusEta"] -> fill(eta_AntiMu);

      const double pT_Mu = _muon.pT()/GeV;
      _hist_MuMinusPt -> fill(pT_Mu);
      //_h["MuMinusPt"] -> fill(pT_Mu);
      const double pT_AntiMu = _antimuon.pT()/GeV;
      _hist_MuPlusPt -> fill(pT_AntiMu);
      //_h["MuPlusPt"] -> fill(pT_AntiMu);

      /// Fill Zboson related histograms
      const double m_Z = (_muon.mom() + _antimuon.mom()).mass()/GeV;
      _hist_ZM -> fill(m_Z);
      //_h["ZM"] -> fill(m_Z);

      const double phi_Z = (_muon.mom() + _antimuon.mom()).phi();
      _hist_ZPhi -> fill(phi_Z);
      //_h["ZPhi"] -> fill(phi_Z);

      const double rap_Z = (_muon.mom() + _antimuon.mom()).rap();
      _hist_ZY -> fill(rap_Z);
      //_h["ZY"] -> fill(rap_Z);

      const double pT_Z = (_muon.mom() + _antimuon.mom()).pT()/GeV;
      _hist_ZPt -> fill(pT_Z);
      //_h["ZPt"] -> fill(pT_Z);

      //const double pi = 3.14159265358979323846;
      const double thetastar = acos(tanh((_antimuon.mom().eta() - _muon.mom().eta())/2));
      const double phistareta = tan(HALFPI - (_antimuon.mom().phi() - _muon.mom().phi())/2)*sin(thetastar);
      _hist_PhiStarEta -> fill(phistareta);
      //_h["PhiStarEta"] -> fill(phistareta);

      /// Fill jet related histograms
      _hist_NJets -> fill(jets.size());
      //_h["NJets"] -> fill(jets.size());

      _hist_Jet1Phi -> fill(jets.at(0).phi());
      _hist_Jet2Phi -> fill(jets.at(1).phi());
      _hist_Jet3Phi -> fill(jets.at(2).phi());
      //_h["Jet1Phi"] -> fill(jets.at(0).phi());
      //_h["Jet2Phi"] -> fill(jets.at(1).phi());
      //_h["Jet3Phi"] -> fill(jets.at(2).phi());

      _hist_Jet1Eta -> fill(jets.at(0).eta());
      //_h["Jet1Eta"] -> fill(jets.at(0).eta());

      _hist_Jet1Y -> fill(jets.at(0).rap());
      _hist_Jet2Y -> fill(jets.at(1).rap());
      _hist_Jet3Y -> fill(jets.at(3).rap());
      //_h["Jet1Y"] -> fill(jets.at(0).rap());
      //_h["Jet2Y"] -> fill(jets.at(1).rap());
      //_h["Jet3Y"] -> fill(jets.at(3).rap());

      const double pT_aveJet = sum(jets, pT, 0)/GeV/jets.size();
      _hist_JetAvePt -> fill(pT_aveJet);
      //_h["JetAvePt"] -> fill(pT_aveJet);

      _hist_Jet1Pt -> fill(jets.at(0).pT()/GeV);
      _hist_Jet2Pt -> fill(jets.at(1).pT()/GeV);
      _hist_Jet3Pt -> fill(jets.at(2).pT()/GeV);
      //_h["Jet1Pt"] -> fill(jets.at(0).pT()/GeV);
      //_h["Jet2Pt"] -> fill(jets.at(1).pT()/GeV);
      //_h["Jet3Pt"] -> fill(jets.at(2).pT()/GeV);

      

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h_YYYY); // normalize to unity
      const double sf = crossSection() / picobarn / sumOfWeights();
      for (auto hist : _h) {
        scale(hist.second, sf); // norm to cross section
      }

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _hist_NMus, _hist_MuPlusPt, _hist_MuPlusEta, _hist_MuPlusPhi, _hist_MuMinusPt, _hist_MuMinusEta, _hist_MuMinusPhi,
              _hist_ZY, _hist_ZM, _hist_ZPhi, _hist_NJets, _hist_Jet1Y, _hist_Jet1Phi, _hist_Jet2Y, _hist_Jet2Phi, _hist_Jet3Y, _hist_Jet3Phi, _hist_MET;
    Histo1DPtr _hist_ZPt, _hist_Jet1Pt, _hist_Jet2Pt, _hist_Jet3Pt, _hist_JetAvePt, _hist_Jet1Eta, _hist_PhiStarEta;

    map<string, Histo1DPtr> _h;
    //Profile1DPtr _p_AAAA;
    //CounterPtr _c_BBBB;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(NPCorrections_ZplusJet);


}
