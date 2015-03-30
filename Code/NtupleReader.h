//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 30 10:49:51 2015 by ROOT version 5.32/00
// from TChain t/
//////////////////////////////////////////////////////////

#ifndef NtupleReader_h
#define NtupleReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtupleReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          DataMC_Type;
   vector<double>  *beamspot_par;
   vector<double>  *beamspot_cov;
   Double_t        beamspot_emittanceX;
   Double_t        beamspot_emittanceY;
   Double_t        beamspot_betaStar;
   vector<double>  *Vtx_chi2;
   vector<unsigned int> *Vtx_nTrk;
   vector<float>   *Vtx_ndof;
   vector<double>  *Vtx_x;
   vector<double>  *Vtx_y;
   vector<double>  *Vtx_z;
   vector<vector<vector<double> > > *Vtx_Cov;
   vector<vector<int> > *Vtx_Track_idx;
   vector<vector<float> > *Vtx_Track_Weights;
   vector<bool>    *Vtx_isFake;
   vector<vector<vector<double> > > *Vtx_TracksP4;
   Bool_t          isPatMuon;
   vector<vector<double> > *Muon_p4;
   vector<vector<double> > *Muon_Poca;
   vector<bool>    *Muon_isGlobalMuon;
   vector<bool>    *Muon_isStandAloneMuon;
   vector<bool>    *Muon_isTrackerMuon;
   vector<bool>    *Muon_isCaloMuon;
   vector<bool>    *Muon_isIsolationValid;
   vector<bool>    *Muon_isQualityValid;
   vector<bool>    *Muon_isTimeValid;
   vector<float>   *Muon_emEt03;
   vector<float>   *Muon_emVetoEt03;
   vector<float>   *Muon_hadEt03;
   vector<float>   *Muon_hadVetoEt03;
   vector<int>     *Muon_nJets03;
   vector<int>     *Muon_nTracks03;
   vector<float>   *Muon_sumPt03;
   vector<float>   *Muon_trackerVetoPt03;
   vector<float>   *Muon_emEt05;
   vector<float>   *Muon_emVetoEt05;
   vector<float>   *Muon_hadEt05;
   vector<float>   *Muon_hadVetoEt05;
   vector<int>     *Muon_nJets05;
   vector<int>     *Muon_nTracks05;
   vector<float>   *Muon_sumPt05;
   vector<float>   *Muon_trackerVetoPt05;
   vector<float>   *Muon_sumChargedHadronPt03;
   vector<float>   *Muon_sumChargedParticlePt03;
   vector<float>   *Muon_sumNeutralHadronEt03;
   vector<float>   *Muon_sumNeutralHadronEtHighThreshold03;
   vector<float>   *Muon_sumPhotonEt03;
   vector<float>   *Muon_sumPhotonEtHighThreshold03;
   vector<float>   *Muon_sumPUPt03;
   vector<float>   *Muon_sumChargedHadronPt04;
   vector<float>   *Muon_sumChargedParticlePt04;
   vector<float>   *Muon_sumNeutralHadronEt04;
   vector<float>   *Muon_sumNeutralHadronEtHighThreshold04;
   vector<float>   *Muon_sumPhotonEt04;
   vector<float>   *Muon_sumPhotonEtHighThreshold04;
   vector<float>   *Muon_sumPUPt04;
   vector<unsigned int> *Muon_Track_idx;
   vector<int>     *Muon_hitPattern_pixelLayerwithMeas;
   vector<int>     *Muon_numberOfMatchedStations;
   vector<float>   *Muon_normChi2;
   vector<int>     *Muon_hitPattern_numberOfValidMuonHits;
   vector<int>     *Muon_innerTrack_numberofValidHits;
   vector<int>     *Muon_numberOfMatches;
   vector<int>     *Muon_numberOfChambers;
   vector<bool>    *Muon_isPFMuon;
   vector<int>     *Muon_numberofValidPixelHits;
   vector<int>     *Muon_trackerLayersWithMeasurement;
   vector<int>     *Muon_charge;
   vector<int>     *Muon_trackCharge;
   vector<int>     *Muon_pdgid;
   vector<double>  *Muon_B;
   vector<double>  *Muon_M;
   vector<vector<double> > *Muon_par;
   vector<vector<double> > *Muon_cov;
   Bool_t          isPatElectron;
   vector<vector<double> > *Electron_p4;
   vector<vector<double> > *Electron_Poca;
   vector<float>   *Electron_Gsf_deltaEtaEleClusterTrackAtCalo;
   vector<float>   *Electron_Gsf_deltaEtaSeedClusterTrackAtCalo;
   vector<float>   *Electron_Gsf_deltaEtaSuperClusterTrackAtVtx;
   vector<float>   *Electron_Gsf_deltaPhiEleClusterTrackAtCalo;
   vector<float>   *Electron_Gsf_deltaPhiSeedClusterTrackAtCalo;
   vector<float>   *Electron_Gsf_deltaPhiSuperClusterTrackAtVtx;
   vector<float>   *Electron_Gsf_dr03EcalRecHitSumE;
   vector<float>   *Electron_Gsf_dr03HcalDepth1TowerSumEt;
   vector<float>   *Electron_Gsf_dr03HcalDepth1TowerSumEtBc;
   vector<float>   *Electron_Gsf_dr03HcalDepth2TowerSumEt;
   vector<float>   *Electron_Gsf_dr03HcalDepth2TowerSumEtBc;
   vector<float>   *Electron_Gsf_dr03HcalTowerSumEt;
   vector<float>   *Electron_Gsf_dr03HcalTowerSumEtBc;
   vector<float>   *Electron_Gsf_dr03TkSumPt;
   vector<bool>    *Electron_Gsf_passingCutBasedPreselection;
   vector<bool>    *Electron_Gsf_passingMvaPreselection;
   vector<int>     *Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits;
   vector<double>  *Electron_supercluster_e;
   vector<double>  *Electron_supercluster_phi;
   vector<double>  *Electron_supercluster_eta;
   vector<float>   *Electron_supercluster_centroid_x;
   vector<float>   *Electron_supercluster_centroid_y;
   vector<float>   *Electron_supercluster_centroid_z;
   vector<unsigned int> *Electron_Track_idx;
   vector<float>   *Electron_ecalRecHitSumEt03;
   vector<float>   *Electron_hcalDepth1TowerSumEt03;
   vector<float>   *Electron_hcalDepth1TowerSumEtBc03;
   vector<float>   *Electron_hcalDepth2TowerSumEt03;
   vector<float>   *Electron_hcalDepth2TowerSumEtBc03;
   vector<float>   *Electron_tkSumPt03;
   vector<float>   *Electron_ecalRecHitSumEt04;
   vector<float>   *Electron_hcalDepth1TowerSumEt04;
   vector<float>   *Electron_hcalDepth1TowerSumEtBc04;
   vector<float>   *Electron_hcalDepth2TowerSumEt04;
   vector<float>   *Electron_hcalDepth2TowerSumEtBc04;
   vector<float>   *Electron_tkSumPt04;
   vector<float>   *Electron_chargedHadronIso;
   vector<float>   *Electron_neutralHadronIso;
   vector<float>   *Electron_photonIso;
   vector<double>  *Electron_isoDeposits_chargedHadronIso04;
   vector<double>  *Electron_isoDeposits_neutralHadronIso04;
   vector<double>  *Electron_isoDeposits_photonIso04;
   vector<double>  *Electron_isoDeposits_chargedHadronIso03;
   vector<double>  *Electron_isoDeposits_neutralHadronIso03;
   vector<double>  *Electron_isoDeposits_photonIso03;
   vector<float>   *Electron_sigmaIetaIeta;
   vector<float>   *Electron_hadronicOverEm;
   vector<float>   *Electron_fbrem;
   vector<float>   *Electron_eSuperClusterOverP;
   vector<float>   *Electron_ecalEnergy;
   vector<float>   *Electron_trackMomentumAtVtx;
   vector<int>     *Electron_numberOfMissedHits;
   vector<bool>    *Electron_HasMatchedConversions;
   Double_t        RhoIsolationAllInputTags;
   vector<int>     *Electron_charge;
   vector<int>     *Electron_trackCharge;
   vector<int>     *Electron_pdgid;
   vector<double>  *Electron_B;
   vector<double>  *Electron_M;
   vector<vector<double> > *Electron_par;
   vector<vector<double> > *Electron_cov;
   vector<double>  *Electron_RegEnergy;
   vector<double>  *Electron_RegEnergyError;
   vector<double>  *Electron_Rho_kt6PFJets;
   vector<double>  *Electron_MVA_TrigNoIP_discriminator;
   vector<double>  *Electron_MVA_NonTrig_discriminator;
   vector<double>  *Electron_MVA_Trig_discriminator;
   vector<vector<double> > *PFTau_p4;
   vector<vector<double> > *PFTau_Poca;
   vector<bool>    *PFTau_isTightIsolation;
   vector<bool>    *PFTau_isMediumIsolation;
   vector<bool>    *PFTau_isLooseIsolation;
   vector<bool>    *PFTau_isTightIsolationDBSumPtCorr;
   vector<bool>    *PFTau_isMediumIsolationDBSumPtCorr;
   vector<bool>    *PFTau_isLooseIsolationDBSumPtCorr;
   vector<bool>    *PFTau_isVLooseIsolationDBSumPtCorr;
   vector<bool>    *PFTau_isHPSAgainstElectronsLoose;
   vector<bool>    *PFTau_isHPSAgainstElectronsMedium;
   vector<bool>    *PFTau_isHPSAgainstElectronsTight;
   vector<bool>    *PFTau_isHPSAgainstMuonLoose;
   vector<bool>    *PFTau_isHPSAgainstMuonMedium;
   vector<bool>    *PFTau_isHPSAgainstMuonTight;
   vector<bool>    *PFTau_isHPSAgainstMuonLoose2;
   vector<bool>    *PFTau_isHPSAgainstMuonMedium2;
   vector<bool>    *PFTau_isHPSAgainstMuonTight2;
   vector<bool>    *PFTau_isHPSByDecayModeFinding;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits;
   vector<float>   *PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByLooseIsolationMVA;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByMediumIsolationMVA;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByTightIsolationMVA;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2;
   vector<bool>    *PFTau_HPSPFTauDiscriminationByTightIsolationMVA2;
   vector<int>     *PFTau_hpsDecayMode;
   vector<int>     *PFTau_Charge;
   vector<vector<int> > *PFTau_Track_idx;
   vector<vector<double> > *PFTau_TIP_primaryVertex_pos;
   vector<vector<double> > *PFTau_TIP_primaryVertex_cov;
   vector<vector<double> > *PFTau_TIP_secondaryVertex_pos;
   vector<vector<double> > *PFTau_TIP_secondaryVertex_cov;
   vector<vector<double> > *PFTau_TIP_secondaryVertex_vtxchi2;
   vector<vector<double> > *PFTau_TIP_secondaryVertex_vtxndof;
   vector<vector<double> > *PFTau_TIP_primaryVertex_vtxchi2;
   vector<vector<double> > *PFTau_TIP_primaryVertex_vtxndof;
   vector<vector<double> > *PFTau_TIP_flightLength;
   vector<vector<double> > *PFTau_TIP_flightLengthSig;
   vector<vector<double> > *PFTau_a1_lvp;
   vector<vector<double> > *PFTau_a1_cov;
   vector<vector<int> > *PFTau_a1_charge;
   vector<vector<int> > *PFTau_a1_pdgid;
   vector<vector<double> > *PFTau_a1_B;
   vector<vector<double> > *PFTau_a1_M;
   vector<vector<vector<double> > > *PFTau_daughterTracks;
   vector<vector<vector<double> > > *PFTau_daughterTracks_cov;
   vector<vector<int> > *PFTau_daughterTracks_charge;
   vector<vector<int> > *PFTau_daughterTracks_pdgid;
   vector<vector<double> > *PFTau_daughterTracks_B;
   vector<vector<double> > *PFTau_daughterTracks_M;
   vector<vector<vector<double> > > *PFTau_daughterTracks_poca;
   vector<vector<double> > *PFTau_3PS_A1_LV;
   vector<vector<double> > *PFTau_3PS_M_A1;
   vector<vector<double> > *PFTau_3PS_M_12;
   vector<vector<double> > *PFTau_3PS_M_13;
   vector<vector<double> > *PFTau_3PS_M_23;
   vector<vector<int> > *PFTau_3PS_Tau_Charge;
   vector<vector<float> > *PFTau_3PS_LCchi2;
   vector<vector<int> > *PFTau_3PS_has3ProngSolution;
   vector<vector<vector<double> > > *PFTau_3PS_Tau_LV;
   vector<vector<vector<double> > > *PFTau_PionsP4;
   vector<vector<int> > *PFTau_PionsCharge;
   vector<vector<vector<double> > > *PFTau_PiZeroP4;
   vector<vector<int> > *PFTau_PiZeroNumOfPhotons;
   vector<vector<int> > *PFTau_PiZeroNumOfElectrons;
   vector<vector<vector<double> > > *PFTau_ChargedHadronsP4;
   vector<vector<vector<int> > > *PFTau_ChargedHadronsCharge;
   vector<vector<vector<double> > > *PFTau_GammaP4;
   vector<vector<vector<double> > > *PFTau_Photons_p4_inDR05;
   vector<vector<double> > *PFTau_MatchedPFJetP4;
   vector<vector<vector<double> > > *PFTau_MatchedPFJetGammasP4;
   vector<vector<vector<double> > > *PFTau_MatchedPFJetSCVariables;
   vector<vector<vector<double> > > *PFTau_MatchedPFJetPhotonVariables;
   vector<double>  *PFTau_PhotonEnergyFraction;
   vector<vector<int> > *PFTau_photon_hasPixelSeed;
   vector<vector<float> > *PFTau_photon_hadronicOverEm;
   vector<vector<float> > *PFTau_photon_sigmaIetaIeta;
   vector<vector<float> > *PFTau_photon_trkSumPtHollowConeDR04;
   vector<vector<float> > *PFTau_photon_ecalRecHitSumEtConeDR04;
   vector<vector<float> > *PFTau_photon_hcalTowerSumEtConeDR04;
   vector<float>   *PFTau_photon_rho;
   Bool_t          isPatJet;
   vector<vector<double> > *PFJet_p4;
   vector<float>   *PFJet_chargedEmEnergy;
   vector<float>   *PFJet_chargedHadronEnergy;
   vector<int>     *PFJet_chargedHadronMultiplicity;
   vector<float>   *PFJet_chargedMuEnergy;
   vector<int>     *PFJet_chargedMultiplicity;
   vector<float>   *PFJet_electronEnergy;
   vector<int>     *PFJet_electronMultiplicity;
   vector<float>   *PFJet_HFEMEnergy;
   vector<int>     *PFJet_HFEMMultiplicity;
   vector<float>   *PFJet_HFHadronEnergy;
   vector<int>     *PFJet_HFHadronMultiplicity;
   vector<float>   *PFJet_muonEnergy;
   vector<int>     *PFJet_muonMultiplicity;
   vector<float>   *PFJet_neutralEmEnergy;
   vector<float>   *PFJet_neutralHadronEnergy;
   vector<int>     *PFJet_neutralHadronMultiplicity;
   vector<float>   *PFJet_photonEnergy;
   vector<int>     *PFJet_photonMultiplicity;
   vector<float>   *PFJet_jetArea;
   vector<float>   *PFJet_maxDistance;
   vector<int>     *PFJet_nConstituents;
   vector<float>   *PFJet_pileup;
   vector<float>   *PFJet_etaetaMoment;
   vector<float>   *PFJet_etaphiMoment;
   vector<vector<int> > *PFJet_Track_idx;
   vector<int>     *PFJet_MatchedHPS_idx;
   vector<int>     *PFJet_numberOfDaughters;
   vector<float>   *PFJet_chargedEmEnergyFraction;
   vector<float>   *PFJet_chargedHadronEnergyFraction;
   vector<float>   *PFJet_neutralHadronEnergyFraction;
   vector<float>   *PFJet_neutralEmEnergyFraction;
   vector<float>   *PFJet_PUJetID_discr;
   vector<bool>    *PFJet_PUJetID_looseWP;
   vector<bool>    *PFJet_PUJetID_mediumWP;
   vector<bool>    *PFJet_PUJetID_tightWP;
   vector<int>     *PFJet_partonFlavour;
   vector<float>   *PFJet_bDiscriminator;
   vector<vector<vector<double> > > *PFJet_TracksP4;
   vector<int>     *PFJet_nTrk;
   vector<float>   *PFJet_JECuncertainty;
   vector<vector<double> > *PFJet_GenJet_p4;
   vector<vector<vector<double> > > *PFJet_GenJet_Constituents_p4;
   vector<vector<double> > *PFJet_GenJetNoNu_p4;
   vector<vector<vector<double> > > *PFJet_GenJetNoNu_Constituents_p4;
   Bool_t          isPatMET;
   Double_t        MET_Uncorr_et;
   Float_t         MET_Uncorr_pt;
   Double_t        MET_Uncorr_phi;
   Float_t         MET_Uncorr_sumET;
   Double_t        MET_Uncorr_significance;
   Double_t        MET_Uncorr_significance_xx;
   Double_t        MET_Uncorr_significance_xy;
   Double_t        MET_Uncorr_significance_yy;
   Float_t         MET_Uncorr_MuonEtFraction;
   Float_t         MET_Uncorr_NeutralEMFraction;
   Float_t         MET_Uncorr_NeutralHadEtFraction;
   Float_t         MET_Uncorr_Type6EtFraction;
   Float_t         MET_Uncorr_Type7EtFraction;
   Double_t        MET_CorrT0rt_et;
   Float_t         MET_CorrT0rt_pt;
   Double_t        MET_CorrT0rt_phi;
   Float_t         MET_CorrT0rt_sumET;
   Float_t         MET_CorrT0rt_MuonEtFraction;
   Float_t         MET_CorrT0rt_NeutralEMFraction;
   Float_t         MET_CorrT0rt_NeutralHadEtFraction;
   Float_t         MET_CorrT0rt_Type6EtFraction;
   Float_t         MET_CorrT0rt_Type7EtFraction;
   Double_t        MET_CorrT0rtT1_et;
   Float_t         MET_CorrT0rtT1_pt;
   Double_t        MET_CorrT0rtT1_phi;
   Float_t         MET_CorrT0rtT1_sumET;
   Float_t         MET_CorrT0rtT1_MuonEtFraction;
   Float_t         MET_CorrT0rtT1_NeutralEMFraction;
   Float_t         MET_CorrT0rtT1_NeutralHadEtFraction;
   Float_t         MET_CorrT0rtT1_Type6EtFraction;
   Float_t         MET_CorrT0rtT1_Type7EtFraction;
   Double_t        MET_CorrT0pc_et;
   Float_t         MET_CorrT0pc_pt;
   Double_t        MET_CorrT0pc_phi;
   Float_t         MET_CorrT0pc_sumET;
   Float_t         MET_CorrT0pc_MuonEtFraction;
   Float_t         MET_CorrT0pc_NeutralEMFraction;
   Float_t         MET_CorrT0pc_NeutralHadEtFraction;
   Float_t         MET_CorrT0pc_Type6EtFraction;
   Float_t         MET_CorrT0pc_Type7EtFraction;
   Double_t        MET_CorrT0pcT1_et;
   Float_t         MET_CorrT0pcT1_pt;
   Double_t        MET_CorrT0pcT1_phi;
   Float_t         MET_CorrT0pcT1_sumET;
   Float_t         MET_CorrT0pcT1_MuonEtFraction;
   Float_t         MET_CorrT0pcT1_NeutralEMFraction;
   Float_t         MET_CorrT0pcT1_NeutralHadEtFraction;
   Float_t         MET_CorrT0pcT1_Type6EtFraction;
   Float_t         MET_CorrT0pcT1_Type7EtFraction;
   Double_t        MET_CorrT0rtTxy_et;
   Float_t         MET_CorrT0rtTxy_pt;
   Double_t        MET_CorrT0rtTxy_phi;
   Float_t         MET_CorrT0rtTxy_sumET;
   Float_t         MET_CorrT0rtTxy_MuonEtFraction;
   Float_t         MET_CorrT0rtTxy_NeutralEMFraction;
   Float_t         MET_CorrT0rtTxy_NeutralHadEtFraction;
   Float_t         MET_CorrT0rtTxy_Type6EtFraction;
   Float_t         MET_CorrT0rtTxy_Type7EtFraction;
   Double_t        MET_CorrT0rtT1Txy_et;
   Float_t         MET_CorrT0rtT1Txy_pt;
   Double_t        MET_CorrT0rtT1Txy_phi;
   Float_t         MET_CorrT0rtT1Txy_sumET;
   Float_t         MET_CorrT0rtT1Txy_MuonEtFraction;
   Float_t         MET_CorrT0rtT1Txy_NeutralEMFraction;
   Float_t         MET_CorrT0rtT1Txy_NeutralHadEtFraction;
   Float_t         MET_CorrT0rtT1Txy_Type6EtFraction;
   Float_t         MET_CorrT0rtT1Txy_Type7EtFraction;
   Double_t        MET_CorrT0pcTxy_et;
   Float_t         MET_CorrT0pcTxy_pt;
   Double_t        MET_CorrT0pcTxy_phi;
   Float_t         MET_CorrT0pcTxy_sumET;
   Float_t         MET_CorrT0pcTxy_MuonEtFraction;
   Float_t         MET_CorrT0pcTxy_NeutralEMFraction;
   Float_t         MET_CorrT0pcTxy_NeutralHadEtFraction;
   Float_t         MET_CorrT0pcTxy_Type6EtFraction;
   Float_t         MET_CorrT0pcTxy_Type7EtFraction;
   Double_t        MET_CorrT0pcT1Txy_et;
   Float_t         MET_CorrT0pcT1Txy_pt;
   Double_t        MET_CorrT0pcT1Txy_phi;
   Float_t         MET_CorrT0pcT1Txy_sumET;
   Float_t         MET_CorrT0pcT1Txy_MuonEtFraction;
   Float_t         MET_CorrT0pcT1Txy_NeutralEMFraction;
   Float_t         MET_CorrT0pcT1Txy_NeutralHadEtFraction;
   Float_t         MET_CorrT0pcT1Txy_Type6EtFraction;
   Float_t         MET_CorrT0pcT1Txy_Type7EtFraction;
   Double_t        MET_CorrT1_et;
   Float_t         MET_CorrT1_pt;
   Double_t        MET_CorrT1_phi;
   Float_t         MET_CorrT1_sumET;
   Float_t         MET_CorrT1_MuonEtFraction;
   Float_t         MET_CorrT1_NeutralEMFraction;
   Float_t         MET_CorrT1_NeutralHadEtFraction;
   Float_t         MET_CorrT1_Type6EtFraction;
   Float_t         MET_CorrT1_Type7EtFraction;
   Double_t        MET_CorrT1Txy_et;
   Float_t         MET_CorrT1Txy_pt;
   Double_t        MET_CorrT1Txy_phi;
   Float_t         MET_CorrT1Txy_sumET;
   Float_t         MET_CorrT1Txy_MuonEtFraction;
   Float_t         MET_CorrT1Txy_NeutralEMFraction;
   Float_t         MET_CorrT1Txy_NeutralHadEtFraction;
   Float_t         MET_CorrT1Txy_Type6EtFraction;
   Float_t         MET_CorrT1Txy_Type7EtFraction;
   Double_t        MET_CorrCaloT1_et;
   Float_t         MET_CorrCaloT1_pt;
   Double_t        MET_CorrCaloT1_phi;
   Float_t         MET_CorrCaloT1_sumET;
   Double_t        MET_CorrCaloT1T2_et;
   Float_t         MET_CorrCaloT1T2_pt;
   Double_t        MET_CorrCaloT1T2_phi;
   Float_t         MET_CorrCaloT1T2_sumET;
   Double_t        MET_CorrMVA_et;
   Float_t         MET_CorrMVA_pt;
   Double_t        MET_CorrMVA_phi;
   Float_t         MET_CorrMVA_sumET;
   Double_t        MET_CorrMVA_significance;
   Double_t        MET_CorrMVA_significance_xx;
   Double_t        MET_CorrMVA_significance_xy;
   Double_t        MET_CorrMVA_significance_yy;
   Float_t         MET_CorrMVA_MuonEtFraction;
   Float_t         MET_CorrMVA_NeutralEMFraction;
   Float_t         MET_CorrMVA_NeutralHadEtFraction;
   Float_t         MET_CorrMVA_Type6EtFraction;
   Float_t         MET_CorrMVA_Type7EtFraction;
   vector<vector<double> > *MET_CorrMVA_srcMuon_p4;
   vector<vector<double> > *MET_CorrMVA_srcElectron_p4;
   vector<vector<double> > *MET_CorrMVA_srcTau_p4;
   Double_t        MET_CorrMVAMuTau_et;
   Float_t         MET_CorrMVAMuTau_pt;
   Double_t        MET_CorrMVAMuTau_phi;
   Float_t         MET_CorrMVAMuTau_sumET;
   Double_t        MET_CorrMVAMuTau_significance;
   Double_t        MET_CorrMVAMuTau_significance_xx;
   Double_t        MET_CorrMVAMuTau_significance_xy;
   Double_t        MET_CorrMVAMuTau_significance_yy;
   Float_t         MET_CorrMVAMuTau_MuonEtFraction;
   Float_t         MET_CorrMVAMuTau_NeutralEMFraction;
   Float_t         MET_CorrMVAMuTau_NeutralHadEtFraction;
   Float_t         MET_CorrMVAMuTau_Type6EtFraction;
   Float_t         MET_CorrMVAMuTau_Type7EtFraction;
   vector<vector<double> > *MET_CorrMVAMuTau_srcMuon_p4;
   vector<vector<double> > *MET_CorrMVAMuTau_srcTau_p4;
   Double_t        MET_Type1Corr_et;
   Float_t         MET_Type1Corr_pt;
   Double_t        MET_Type1Corr_phi;
   Float_t         MET_Type1Corr_sumET;
   Float_t         MET_Type1Corr_MuonEtFraction;
   Float_t         MET_Type1Corr_NeutralEMFraction;
   Float_t         MET_Type1Corr_NeutralHadEtFraction;
   Float_t         MET_Type1Corr_Type6EtFraction;
   Float_t         MET_Type1Corr_Type7EtFraction;
   Double_t        MET_Type1p2Corr_et;
   Float_t         MET_Type1p2Corr_pt;
   Double_t        MET_Type1p2Corr_phi;
   Float_t         MET_Type1p2Corr_sumET;
   Float_t         MET_Type1p2Corr_MuonEtFraction;
   Float_t         MET_Type1p2Corr_NeutralEMFraction;
   Float_t         MET_Type1p2Corr_NeutralHadEtFraction;
   Float_t         MET_Type1p2Corr_Type6EtFraction;
   Float_t         MET_Type1p2Corr_Type7EtFraction;
   Double_t        MET_Type1CorrElectronUp_et;
   Double_t        MET_Type1CorrElectronDown_et;
   Double_t        MET_Type1CorrMuonUp_et;
   Double_t        MET_Type1CorrMuonDown_et;
   Double_t        MET_Type1CorrTauUp_et;
   Double_t        MET_Type1CorrTauDown_et;
   Double_t        MET_Type1CorrJetResUp_et;
   Double_t        MET_Type1CorrJetResDown_et;
   Double_t        MET_Type1CorrJetEnUp_et;
   Double_t        MET_Type1CorrJetEnDown_et;
   Double_t        MET_Type1CorrUnclusteredUp_et;
   Double_t        MET_Type1CorrUnclusteredDown_et;
   Double_t        MET_Type1p2CorrElectronUp_et;
   Double_t        MET_Type1p2CorrElectronDown_et;
   Double_t        MET_Type1p2CorrMuonUp_et;
   Double_t        MET_Type1p2CorrMuonDown_et;
   Double_t        MET_Type1p2CorrTauUp_et;
   Double_t        MET_Type1p2CorrTauDown_et;
   Double_t        MET_Type1p2CorrJetResUp_et;
   Double_t        MET_Type1p2CorrJetResDown_et;
   Double_t        MET_Type1p2CorrJetEnUp_et;
   Double_t        MET_Type1p2CorrJetEnDown_et;
   Double_t        MET_Type1p2CorrUnclusteredUp_et;
   Double_t        MET_Type1p2CorrUnclusteredDown_et;
   UInt_t          Event_EventNumber;
   UInt_t          Event_RunNumber;
   Int_t           Event_bunchCrossing;
   Int_t           Event_orbitNumber;
   UInt_t          Event_luminosityBlock;
   Bool_t          Event_isRealData;
   Float_t         PileupInfo_TrueNumInteractions_nm1;
   Float_t         PileupInfo_TrueNumInteractions_n0;
   Float_t         PileupInfo_TrueNumInteractions_np1;
   Float_t         PUWeight;
   Float_t         PUWeight_p5;
   Float_t         PUWeight_m5;
   Float_t         PUWeight3D;
   Float_t         PUWeight3D_p5;
   Float_t         PUWeight3D_m5;
   Float_t         PUWeightFineBins;
   Float_t         TauSpinnerWeight;
   Float_t         SelEffWeight;
   Float_t         MinVisPtFilter;
   Float_t         KinWeightPt;
   Float_t         KinWeightEta;
   Float_t         KinWeightMassPt;
   Float_t         EmbeddedWeight;
   vector<vector<double> > *Track_p4;
   vector<vector<double> > *Track_Poca;
   vector<double>  *Track_chi2;
   vector<double>  *Track_ndof;
   vector<unsigned short> *Track_numberOfLostHits;
   vector<unsigned short> *Track_numberOfValidHits;
   vector<unsigned int> *Track_qualityMask;
   vector<int>     *Track_charge;
   vector<int>     *Track_pdgid;
   vector<double>  *Track_B;
   vector<double>  *Track_M;
   vector<vector<double> > *Track_par;
   vector<vector<double> > *Track_cov;
   UInt_t          GenEventInfoProduct_signalProcessID;
   Float_t         GenEventInfoProduct_weight;
   vector<double>  *GenEventInfoProduct_weights;
   Float_t         GenEventInfoProduct_qScale;
   Float_t         GenEventInfoProduct_alphaQED;
   Float_t         GenEventInfoProduct_alphaQCD;
   Int_t           GenEventInfoProduct_id1;
   Int_t           GenEventInfoProduct_id2;
   Double_t        GenEventInfoProduct_x1;
   Double_t        GenEventInfoProduct_x2;
   Double_t        GenEventInfoProduct_scalePDF;
   vector<vector<float> > *MC_p4;
   vector<int>     *MC_pdgid;
   vector<int>     *MC_charge;
   vector<int> *MC_midx;
   vector<vector<int> > *MC_childpdgid;
   vector<vector<int> > *MC_childidx;
   vector<int>     *MC_status;
   vector<vector<double> > *MCSignalParticle_p4;
   vector<int>     *MCSignalParticle_pdgid;
   vector<vector<int> > *MCSignalParticle_childpdgid;
   vector<int>     *MCSignalParticle_charge;
   vector<vector<double> > *MCSignalParticle_Poca;
   vector<vector<unsigned int> > *MCSignalParticle_Tauidx;
   vector<vector<vector<double> > > *MCTauandProd_p4;
   vector<vector<vector<double> > > *MCTauandProd_Vertex;
   vector<vector<int> > *MCTauandProd_pdgid;
   vector<vector<unsigned int> > *MCTauandProd_midx;
   vector<vector<int> > *MCTauandProd_charge;
   vector<unsigned int> *MCTau_JAK;
   vector<unsigned int> *MCTau_DecayBitMask;
   vector<string>  *HTLTriggerName;
   vector<bool>    *TriggerAccept;
   vector<bool>    *TriggerError;
   vector<bool>    *TriggerWasRun;
   vector<unsigned int> *HLTPrescale;
   vector<unsigned int> *NHLTL1GTSeeds;
   vector<unsigned int> *L1SEEDPrescale;
   vector<bool>    *L1SEEDInvalidPrescale;
   vector<bool>    *L1SEEDisTechBit;
   vector<vector<float> > *TriggerMatchMuon;
   vector<vector<float> > *TriggerMatchJet;
   vector<vector<float> > *TriggerMatchTau;
   vector<vector<float> > *HLTTrigger_objs_Pt;
   vector<vector<float> > *HLTTrigger_objs_Eta;
   vector<vector<float> > *HLTTrigger_objs_Phi;
   vector<vector<float> > *HLTTrigger_objs_E;
   vector<vector<int> > *HLTTrigger_objs_Id;
   vector<string>  *HLTTrigger_objs_trigger;
   vector<string>  *L1TriggerName;
   vector<bool>    *L1TriggerDecision;
   vector<int>     *L1ErrorCode;
   vector<unsigned int> *L1Prescale;

   // List of branches
   TBranch        *b_DataMC_Type;   //!
   TBranch        *b_beamspot_par;   //!
   TBranch        *b_beamspot_cov;   //!
   TBranch        *b_beamspot_emittanceX;   //!
   TBranch        *b_beamspot_emittanceY;   //!
   TBranch        *b_beamspot_betaStar;   //!
   TBranch        *b_Vtx_chi2;   //!
   TBranch        *b_Vtx_nTrk;   //!
   TBranch        *b_Vtx_ndof;   //!
   TBranch        *b_Vtx_x;   //!
   TBranch        *b_Vtx_y;   //!
   TBranch        *b_Vtx_z;   //!
   TBranch        *b_Vtx_Cov;   //!
   TBranch        *b_Vtx_Track_idx;   //!
   TBranch        *b_Vtx_Track_Weights;   //!
   TBranch        *b_Vtx_isFake;   //!
   TBranch        *b_Vtx_TracksP4;   //!
   TBranch        *b_isPatMuon;   //!
   TBranch        *b_Muon_p4;   //!
   TBranch        *b_Muon_Poca;   //!
   TBranch        *b_Muon_isGlobalMuon;   //!
   TBranch        *b_Muon_isStandAloneMuon;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_Muon_isCaloMuon;   //!
   TBranch        *b_Muon_isIsolationValid;   //!
   TBranch        *b_Muon_isQualityValid;   //!
   TBranch        *b_Muon_isTimeValid;   //!
   TBranch        *b_Muon_emEt03;   //!
   TBranch        *b_Muon_emVetoEt03;   //!
   TBranch        *b_Muon_hadEt03;   //!
   TBranch        *b_Muon_hadVetoEt03;   //!
   TBranch        *b_Muon_nJets03;   //!
   TBranch        *b_Muon_nTracks03;   //!
   TBranch        *b_Muon_sumPt03;   //!
   TBranch        *b_Muon_trackerVetoPt03;   //!
   TBranch        *b_Muon_emEt05;   //!
   TBranch        *b_Muon_emVetoEt05;   //!
   TBranch        *b_Muon_hadEt05;   //!
   TBranch        *b_Muon_hadVetoEt05;   //!
   TBranch        *b_Muon_nJets05;   //!
   TBranch        *b_Muon_nTracks05;   //!
   TBranch        *b_Muon_sumPt05;   //!
   TBranch        *b_Muon_trackerVetoPt05;   //!
   TBranch        *b_Muon_sumChargedHadronPt03;   //!
   TBranch        *b_Muon_sumChargedParticlePt03;   //!
   TBranch        *b_Muon_sumNeutralHadronEt03;   //!
   TBranch        *b_Muon_sumNeutralHadronEtHighThreshold03;   //!
   TBranch        *b_Muon_sumPhotonEt03;   //!
   TBranch        *b_Muon_sumPhotonEtHighThreshold03;   //!
   TBranch        *b_Muon_sumPUPt03;   //!
   TBranch        *b_Muon_sumChargedHadronPt04;   //!
   TBranch        *b_Muon_sumChargedParticlePt04;   //!
   TBranch        *b_Muon_sumNeutralHadronEt04;   //!
   TBranch        *b_Muon_sumNeutralHadronEtHighThreshold04;   //!
   TBranch        *b_Muon_sumPhotonEt04;   //!
   TBranch        *b_Muon_sumPhotonEtHighThreshold04;   //!
   TBranch        *b_Muon_sumPUPt04;   //!
   TBranch        *b_Muon_Track_idx;   //!
   TBranch        *b_Muon_hitPattern_pixelLayerwithMeas;   //!
   TBranch        *b_Muon_numberOfMatchedStations;   //!
   TBranch        *b_Muon_normChi2;   //!
   TBranch        *b_Muon_hitPattern_numberOfValidMuonHits;   //!
   TBranch        *b_Muon_innerTrack_numberofValidHits;   //!
   TBranch        *b_Muon_numberOfMatches;   //!
   TBranch        *b_Muon_numberOfChambers;   //!
   TBranch        *b_Muon_isPFMuon;   //!
   TBranch        *b_Muon_numberofValidPixelHits;   //!
   TBranch        *b_Muon_trackerLayersWithMeasurement;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_trackCharge;   //!
   TBranch        *b_Muon_pdgid;   //!
   TBranch        *b_Muon_B;   //!
   TBranch        *b_Muon_M;   //!
   TBranch        *b_Muon_par;   //!
   TBranch        *b_Muon_cov;   //!
   TBranch        *b_isPatElectron;   //!
   TBranch        *b_Electron_p4;   //!
   TBranch        *b_Electron_Poca;   //!
   TBranch        *b_Electron_Gsf_deltaEtaEleClusterTrackAtCalo;   //!
   TBranch        *b_Electron_Gsf_deltaEtaSeedClusterTrackAtCalo;   //!
   TBranch        *b_Electron_Gsf_deltaEtaSuperClusterTrackAtVtx;   //!
   TBranch        *b_Electron_Gsf_deltaPhiEleClusterTrackAtCalo;   //!
   TBranch        *b_Electron_Gsf_deltaPhiSeedClusterTrackAtCalo;   //!
   TBranch        *b_Electron_Gsf_deltaPhiSuperClusterTrackAtVtx;   //!
   TBranch        *b_Electron_Gsf_dr03EcalRecHitSumE;   //!
   TBranch        *b_Electron_Gsf_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electron_Gsf_dr03HcalDepth1TowerSumEtBc;   //!
   TBranch        *b_Electron_Gsf_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_Electron_Gsf_dr03HcalDepth2TowerSumEtBc;   //!
   TBranch        *b_Electron_Gsf_dr03HcalTowerSumEt;   //!
   TBranch        *b_Electron_Gsf_dr03HcalTowerSumEtBc;   //!
   TBranch        *b_Electron_Gsf_dr03TkSumPt;   //!
   TBranch        *b_Electron_Gsf_passingCutBasedPreselection;   //!
   TBranch        *b_Electron_Gsf_passingMvaPreselection;   //!
   TBranch        *b_Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits;   //!
   TBranch        *b_Electron_supercluster_e;   //!
   TBranch        *b_Electron_supercluster_phi;   //!
   TBranch        *b_Electron_supercluster_eta;   //!
   TBranch        *b_Electron_supercluster_centroid_x;   //!
   TBranch        *b_Electron_supercluster_centroid_y;   //!
   TBranch        *b_Electron_supercluster_centroid_z;   //!
   TBranch        *b_Electron_Track_idx;   //!
   TBranch        *b_Electron_ecalRecHitSumEt03;   //!
   TBranch        *b_Electron_hcalDepth1TowerSumEt03;   //!
   TBranch        *b_Electron_hcalDepth1TowerSumEtBc03;   //!
   TBranch        *b_Electron_hcalDepth2TowerSumEt03;   //!
   TBranch        *b_Electron_hcalDepth2TowerSumEtBc03;   //!
   TBranch        *b_Electron_tkSumPt03;   //!
   TBranch        *b_Electron_ecalRecHitSumEt04;   //!
   TBranch        *b_Electron_hcalDepth1TowerSumEt04;   //!
   TBranch        *b_Electron_hcalDepth1TowerSumEtBc04;   //!
   TBranch        *b_Electron_hcalDepth2TowerSumEt04;   //!
   TBranch        *b_Electron_hcalDepth2TowerSumEtBc04;   //!
   TBranch        *b_Electron_tkSumPt04;   //!
   TBranch        *b_Electron_chargedHadronIso;   //!
   TBranch        *b_Electron_neutralHadronIso;   //!
   TBranch        *b_Electron_photonIso;   //!
   TBranch        *b_Electron_isoDeposits_chargedHadronIso04;   //!
   TBranch        *b_Electron_isoDeposits_neutralHadronIso04;   //!
   TBranch        *b_Electron_isoDeposits_photonIso04;   //!
   TBranch        *b_Electron_isoDeposits_chargedHadronIso03;   //!
   TBranch        *b_Electron_isoDeposits_neutralHadronIso03;   //!
   TBranch        *b_Electron_isoDeposits_photonIso03;   //!
   TBranch        *b_Electron_sigmaIetaIeta;   //!
   TBranch        *b_Electron_hadronicOverEm;   //!
   TBranch        *b_Electron_fbrem;   //!
   TBranch        *b_Electron_eSuperClusterOverP;   //!
   TBranch        *b_Electron_ecalEnergy;   //!
   TBranch        *b_Electron_trackMomentumAtVtx;   //!
   TBranch        *b_Electron_numberOfMissedHits;   //!
   TBranch        *b_Electron_HasMatchedConversions;   //!
   TBranch        *b_RhoIsolationAllInputTags;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_trackCharge;   //!
   TBranch        *b_Electron_pdgid;   //!
   TBranch        *b_Electron_B;   //!
   TBranch        *b_Electron_M;   //!
   TBranch        *b_Electron_par;   //!
   TBranch        *b_Electron_cov;   //!
   TBranch        *b_Electron_RegEnergy;   //!
   TBranch        *b_Electron_RegEnergyError;   //!
   TBranch        *b_Electron_Rho_kt6PFJets;   //!
   TBranch        *b_Electron_MVA_TrigNoIP_discriminator;   //!
   TBranch        *b_Electron_MVA_NonTrig_discriminator;   //!
   TBranch        *b_Electron_MVA_Trig_discriminator;   //!
   TBranch        *b_PFTau_p4;   //!
   TBranch        *b_PFTau_Poca;   //!
   TBranch        *b_PFTau_isTightIsolation;   //!
   TBranch        *b_PFTau_isMediumIsolation;   //!
   TBranch        *b_PFTau_isLooseIsolation;   //!
   TBranch        *b_PFTau_isTightIsolationDBSumPtCorr;   //!
   TBranch        *b_PFTau_isMediumIsolationDBSumPtCorr;   //!
   TBranch        *b_PFTau_isLooseIsolationDBSumPtCorr;   //!
   TBranch        *b_PFTau_isVLooseIsolationDBSumPtCorr;   //!
   TBranch        *b_PFTau_isHPSAgainstElectronsLoose;   //!
   TBranch        *b_PFTau_isHPSAgainstElectronsMedium;   //!
   TBranch        *b_PFTau_isHPSAgainstElectronsTight;   //!
   TBranch        *b_PFTau_isHPSAgainstMuonLoose;   //!
   TBranch        *b_PFTau_isHPSAgainstMuonMedium;   //!
   TBranch        *b_PFTau_isHPSAgainstMuonTight;   //!
   TBranch        *b_PFTau_isHPSAgainstMuonLoose2;   //!
   TBranch        *b_PFTau_isHPSAgainstMuonMedium2;   //!
   TBranch        *b_PFTau_isHPSAgainstMuonTight2;   //!
   TBranch        *b_PFTau_isHPSByDecayModeFinding;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByLooseIsolationMVA;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByMediumIsolationMVA;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByTightIsolationMVA;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2;   //!
   TBranch        *b_PFTau_HPSPFTauDiscriminationByTightIsolationMVA2;   //!
   TBranch        *b_PFTau_hpsDecayMode;   //!
   TBranch        *b_PFTau_Charge;   //!
   TBranch        *b_PFTau_Track_idx;   //!
   TBranch        *b_PFTau_TIP_primaryVertex_pos;   //!
   TBranch        *b_PFTau_TIP_primaryVertex_cov;   //!
   TBranch        *b_PFTau_TIP_secondaryVertex_pos;   //!
   TBranch        *b_PFTau_TIP_secondaryVertex_cov;   //!
   TBranch        *b_PFTau_TIP_secondaryVertex_vtxchi2;   //!
   TBranch        *b_PFTau_TIP_secondaryVertex_vtxndof;   //!
   TBranch        *b_PFTau_TIP_primaryVertex_vtxchi2;   //!
   TBranch        *b_PFTau_TIP_primaryVertex_vtxndof;   //!
   TBranch        *b_PFTau_TIP_flightLength;   //!
   TBranch        *b_PFTau_TIP_flightLengthSig;   //!
   TBranch        *b_PFTau_a1_lvp;   //!
   TBranch        *b_PFTau_a1_cov;   //!
   TBranch        *b_PFTau_a1_charge;   //!
   TBranch        *b_PFTau_a1_pdgid;   //!
   TBranch        *b_PFTau_a1_B;   //!
   TBranch        *b_PFTau_a1_M;   //!
   TBranch        *b_PFTau_daughterTracks;   //!
   TBranch        *b_PFTau_daughterTracks_cov;   //!
   TBranch        *b_PFTau_daughterTracks_charge;   //!
   TBranch        *b_PFTau_daughterTracks_pdgid;   //!
   TBranch        *b_PFTau_daughterTracks_B;   //!
   TBranch        *b_PFTau_daughterTracks_M;   //!
   TBranch        *b_PFTau_daughterTracks_poca;   //!
   TBranch        *b_PFTau_3PS_A1_LV;   //!
   TBranch        *b_PFTau_3PS_M_A1;   //!
   TBranch        *b_PFTau_3PS_M_12;   //!
   TBranch        *b_PFTau_3PS_M_13;   //!
   TBranch        *b_PFTau_3PS_M_23;   //!
   TBranch        *b_PFTau_3PS_Tau_Charge;   //!
   TBranch        *b_PFTau_3PS_LCchi2;   //!
   TBranch        *b_PFTau_3PS_has3ProngSolution;   //!
   TBranch        *b_PFTau_3PS_Tau_LV;   //!
   TBranch        *b_PFTau_PionsP4;   //!
   TBranch        *b_PFTau_PionsCharge;   //!
   TBranch        *b_PFTau_PiZeroP4;   //!
   TBranch        *b_PFTau_PiZeroNumOfPhotons;   //!
   TBranch        *b_PFTau_PiZeroNumOfElectrons;   //!
   TBranch        *b_PFTau_ChargedHadronsP4;   //!
   TBranch        *b_PFTau_ChargedHadronsCharge;   //!
   TBranch        *b_PFTau_GammaP4;   //!
   TBranch        *b_PFTau_Photons_p4_inDR05;   //!
   TBranch        *b_PFTau_MatchedPFJetP4;   //!
   TBranch        *b_PFTau_MatchedPFJetGammasP4;   //!
   TBranch        *b_PFTau_MatchedPFJetSCVariables;   //!
   TBranch        *b_PFTau_MatchedPFJetPhotonVariables;   //!
   TBranch        *b_PFTau_PhotonEnergyFraction;   //!
   TBranch        *b_PFTau_photon_hasPixelSeed;   //!
   TBranch        *b_PFTau_photon_hadronicOverEm;   //!
   TBranch        *b_PFTau_photon_sigmaIetaIeta;   //!
   TBranch        *b_PFTau_photon_trkSumPtHollowConeDR04;   //!
   TBranch        *b_PFTau_photon_ecalRecHitSumEtConeDR04;   //!
   TBranch        *b_PFTau_photon_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_PFTau_photon_rho;   //!
   TBranch        *b_isPatJet;   //!
   TBranch        *b_PFJet_p4;   //!
   TBranch        *b_PFJet_chargedEmEnergy;   //!
   TBranch        *b_PFJet_chargedHadronEnergy;   //!
   TBranch        *b_PFJet_chargedHadronMultiplicity;   //!
   TBranch        *b_PFJet_chargedMuEnergy;   //!
   TBranch        *b_PFJet_chargedMultiplicity;   //!
   TBranch        *b_PFJet_electronEnergy;   //!
   TBranch        *b_PFJet_electronMultiplicity;   //!
   TBranch        *b_PFJet_HFEMEnergy;   //!
   TBranch        *b_PFJet_HFEMMultiplicity;   //!
   TBranch        *b_PFJet_HFHadronEnergy;   //!
   TBranch        *b_PFJet_HFHadronMultiplicity;   //!
   TBranch        *b_PFJet_muonEnergy;   //!
   TBranch        *b_PFJet_muonMultiplicity;   //!
   TBranch        *b_PFJet_neutralEmEnergy;   //!
   TBranch        *b_PFJet_neutralHadronEnergy;   //!
   TBranch        *b_PFJet_neutralHadronMultiplicity;   //!
   TBranch        *b_PFJet_photonEnergy;   //!
   TBranch        *b_PFJet_photonMultiplicity;   //!
   TBranch        *b_PFJet_jetArea;   //!
   TBranch        *b_PFJet_maxDistance;   //!
   TBranch        *b_PFJet_nConstituents;   //!
   TBranch        *b_PFJet_pileup;   //!
   TBranch        *b_PFJet_etaetaMoment;   //!
   TBranch        *b_PFJet_etaphiMoment;   //!
   TBranch        *b_PFJet_Track_idx;   //!
   TBranch        *b_PFJet_MatchedHPS_idx;   //!
   TBranch        *b_PFJet_numberOfDaughters;   //!
   TBranch        *b_PFJet_chargedEmEnergyFraction;   //!
   TBranch        *b_PFJet_chargedHadronEnergyFraction;   //!
   TBranch        *b_PFJet_neutralHadronEnergyFraction;   //!
   TBranch        *b_PFJet_neutralEmEnergyFraction;   //!
   TBranch        *b_PFJet_PUJetID_discr;   //!
   TBranch        *b_PFJet_PUJetID_looseWP;   //!
   TBranch        *b_PFJet_PUJetID_mediumWP;   //!
   TBranch        *b_PFJet_PUJetID_tightWP;   //!
   TBranch        *b_PFJet_partonFlavour;   //!
   TBranch        *b_PFJet_bDiscriminator;   //!
   TBranch        *b_PFJet_TracksP4;   //!
   TBranch        *b_PFJet_nTrk;   //!
   TBranch        *b_PFJet_JECuncertainty;   //!
   TBranch        *b_PFJet_GenJet_p4;   //!
   TBranch        *b_PFJet_GenJet_Constituents_p4;   //!
   TBranch        *b_PFJet_GenJetNoNu_p4;   //!
   TBranch        *b_PFJet_GenJetNoNu_Constituents_p4;   //!
   TBranch        *b_isPatMET;   //!
   TBranch        *b_MET_Uncorr_et;   //!
   TBranch        *b_MET_Uncorr_pt;   //!
   TBranch        *b_MET_Uncorr_phi;   //!
   TBranch        *b_MET_Uncorr_sumET;   //!
   TBranch        *b_MET_Uncorr_significance;   //!
   TBranch        *b_MET_Uncorr_significance_xx;   //!
   TBranch        *b_MET_Uncorr_significance_xy;   //!
   TBranch        *b_MET_Uncorr_significance_yy;   //!
   TBranch        *b_MET_Uncorr_MuonEtFraction;   //!
   TBranch        *b_MET_Uncorr_NeutralEMFraction;   //!
   TBranch        *b_MET_Uncorr_NeutralHadEtFraction;   //!
   TBranch        *b_MET_Uncorr_Type6EtFraction;   //!
   TBranch        *b_MET_Uncorr_Type7EtFraction;   //!
   TBranch        *b_MET_CorrT0rt_et;   //!
   TBranch        *b_MET_CorrT0rt_pt;   //!
   TBranch        *b_MET_CorrT0rt_phi;   //!
   TBranch        *b_MET_CorrT0rt_sumET;   //!
   TBranch        *b_MET_CorrT0rt_MuonEtFraction;   //!
   TBranch        *b_MET_CorrT0rt_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrT0rt_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrT0rt_Type6EtFraction;   //!
   TBranch        *b_MET_CorrT0rt_Type7EtFraction;   //!
   TBranch        *b_MET_CorrT0rtT1_et;   //!
   TBranch        *b_MET_CorrT0rtT1_pt;   //!
   TBranch        *b_MET_CorrT0rtT1_phi;   //!
   TBranch        *b_MET_CorrT0rtT1_sumET;   //!
   TBranch        *b_MET_CorrT0rtT1_MuonEtFraction;   //!
   TBranch        *b_MET_CorrT0rtT1_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrT0rtT1_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrT0rtT1_Type6EtFraction;   //!
   TBranch        *b_MET_CorrT0rtT1_Type7EtFraction;   //!
   TBranch        *b_MET_CorrT0pc_et;   //!
   TBranch        *b_MET_CorrT0pc_pt;   //!
   TBranch        *b_MET_CorrT0pc_phi;   //!
   TBranch        *b_MET_CorrT0pc_sumET;   //!
   TBranch        *b_MET_CorrT0pc_MuonEtFraction;   //!
   TBranch        *b_MET_CorrT0pc_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrT0pc_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrT0pc_Type6EtFraction;   //!
   TBranch        *b_MET_CorrT0pc_Type7EtFraction;   //!
   TBranch        *b_MET_CorrT0pcT1_et;   //!
   TBranch        *b_MET_CorrT0pcT1_pt;   //!
   TBranch        *b_MET_CorrT0pcT1_phi;   //!
   TBranch        *b_MET_CorrT0pcT1_sumET;   //!
   TBranch        *b_MET_CorrT0pcT1_MuonEtFraction;   //!
   TBranch        *b_MET_CorrT0pcT1_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrT0pcT1_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrT0pcT1_Type6EtFraction;   //!
   TBranch        *b_MET_CorrT0pcT1_Type7EtFraction;   //!
   TBranch        *b_MET_CorrT0rtTxy_et;   //!
   TBranch        *b_MET_CorrT0rtTxy_pt;   //!
   TBranch        *b_MET_CorrT0rtTxy_phi;   //!
   TBranch        *b_MET_CorrT0rtTxy_sumET;   //!
   TBranch        *b_MET_CorrT0rtTxy_MuonEtFraction;   //!
   TBranch        *b_MET_CorrT0rtTxy_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrT0rtTxy_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrT0rtTxy_Type6EtFraction;   //!
   TBranch        *b_MET_CorrT0rtTxy_Type7EtFraction;   //!
   TBranch        *b_MET_CorrT0rtT1Txy_et;   //!
   TBranch        *b_MET_CorrT0rtT1Txy_pt;   //!
   TBranch        *b_MET_CorrT0rtT1Txy_phi;   //!
   TBranch        *b_MET_CorrT0rtT1Txy_sumET;   //!
   TBranch        *b_MET_CorrT0rtT1Txy_MuonEtFraction;   //!
   TBranch        *b_MET_CorrT0rtT1Txy_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrT0rtT1Txy_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrT0rtT1Txy_Type6EtFraction;   //!
   TBranch        *b_MET_CorrT0rtT1Txy_Type7EtFraction;   //!
   TBranch        *b_MET_CorrT0pcTxy_et;   //!
   TBranch        *b_MET_CorrT0pcTxy_pt;   //!
   TBranch        *b_MET_CorrT0pcTxy_phi;   //!
   TBranch        *b_MET_CorrT0pcTxy_sumET;   //!
   TBranch        *b_MET_CorrT0pcTxy_MuonEtFraction;   //!
   TBranch        *b_MET_CorrT0pcTxy_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrT0pcTxy_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrT0pcTxy_Type6EtFraction;   //!
   TBranch        *b_MET_CorrT0pcTxy_Type7EtFraction;   //!
   TBranch        *b_MET_CorrT0pcT1Txy_et;   //!
   TBranch        *b_MET_CorrT0pcT1Txy_pt;   //!
   TBranch        *b_MET_CorrT0pcT1Txy_phi;   //!
   TBranch        *b_MET_CorrT0pcT1Txy_sumET;   //!
   TBranch        *b_MET_CorrT0pcT1Txy_MuonEtFraction;   //!
   TBranch        *b_MET_CorrT0pcT1Txy_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrT0pcT1Txy_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrT0pcT1Txy_Type6EtFraction;   //!
   TBranch        *b_MET_CorrT0pcT1Txy_Type7EtFraction;   //!
   TBranch        *b_MET_CorrT1_et;   //!
   TBranch        *b_MET_CorrT1_pt;   //!
   TBranch        *b_MET_CorrT1_phi;   //!
   TBranch        *b_MET_CorrT1_sumET;   //!
   TBranch        *b_MET_CorrT1_MuonEtFraction;   //!
   TBranch        *b_MET_CorrT1_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrT1_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrT1_Type6EtFraction;   //!
   TBranch        *b_MET_CorrT1_Type7EtFraction;   //!
   TBranch        *b_MET_CorrT1Txy_et;   //!
   TBranch        *b_MET_CorrT1Txy_pt;   //!
   TBranch        *b_MET_CorrT1Txy_phi;   //!
   TBranch        *b_MET_CorrT1Txy_sumET;   //!
   TBranch        *b_MET_CorrT1Txy_MuonEtFraction;   //!
   TBranch        *b_MET_CorrT1Txy_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrT1Txy_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrT1Txy_Type6EtFraction;   //!
   TBranch        *b_MET_CorrT1Txy_Type7EtFraction;   //!
   TBranch        *b_MET_CorrCaloT1_et;   //!
   TBranch        *b_MET_CorrCaloT1_pt;   //!
   TBranch        *b_MET_CorrCaloT1_phi;   //!
   TBranch        *b_MET_CorrCaloT1_sumET;   //!
   TBranch        *b_MET_CorrCaloT1T2_et;   //!
   TBranch        *b_MET_CorrCaloT1T2_pt;   //!
   TBranch        *b_MET_CorrCaloT1T2_phi;   //!
   TBranch        *b_MET_CorrCaloT1T2_sumET;   //!
   TBranch        *b_MET_CorrMVA_et;   //!
   TBranch        *b_MET_CorrMVA_pt;   //!
   TBranch        *b_MET_CorrMVA_phi;   //!
   TBranch        *b_MET_CorrMVA_sumET;   //!
   TBranch        *b_MET_CorrMVA_significance;   //!
   TBranch        *b_MET_CorrMVA_significance_xx;   //!
   TBranch        *b_MET_CorrMVA_significance_xy;   //!
   TBranch        *b_MET_CorrMVA_significance_yy;   //!
   TBranch        *b_MET_CorrMVA_MuonEtFraction;   //!
   TBranch        *b_MET_CorrMVA_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrMVA_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrMVA_Type6EtFraction;   //!
   TBranch        *b_MET_CorrMVA_Type7EtFraction;   //!
   TBranch        *b_MET_CorrMVA_srcMuon_p4;   //!
   TBranch        *b_MET_CorrMVA_srcElectron_p4;   //!
   TBranch        *b_MET_CorrMVA_srcTau_p4;   //!
   TBranch        *b_MET_CorrMVAMuTau_et;   //!
   TBranch        *b_MET_CorrMVAMuTau_pt;   //!
   TBranch        *b_MET_CorrMVAMuTau_phi;   //!
   TBranch        *b_MET_CorrMVAMuTau_sumET;   //!
   TBranch        *b_MET_CorrMVAMuTau_significance;   //!
   TBranch        *b_MET_CorrMVAMuTau_significance_xx;   //!
   TBranch        *b_MET_CorrMVAMuTau_significance_xy;   //!
   TBranch        *b_MET_CorrMVAMuTau_significance_yy;   //!
   TBranch        *b_MET_CorrMVAMuTau_MuonEtFraction;   //!
   TBranch        *b_MET_CorrMVAMuTau_NeutralEMFraction;   //!
   TBranch        *b_MET_CorrMVAMuTau_NeutralHadEtFraction;   //!
   TBranch        *b_MET_CorrMVAMuTau_Type6EtFraction;   //!
   TBranch        *b_MET_CorrMVAMuTau_Type7EtFraction;   //!
   TBranch        *b_MET_CorrMVAMuTau_srcMuon_p4;   //!
   TBranch        *b_MET_CorrMVAMuTau_srcTau_p4;   //!
   TBranch        *b_MET_Type1Corr_et;   //!
   TBranch        *b_MET_Type1Corr_pt;   //!
   TBranch        *b_MET_Type1Corr_phi;   //!
   TBranch        *b_MET_Type1Corr_sumET;   //!
   TBranch        *b_MET_Type1Corr_MuonEtFraction;   //!
   TBranch        *b_MET_Type1Corr_NeutralEMFraction;   //!
   TBranch        *b_MET_Type1Corr_NeutralHadEtFraction;   //!
   TBranch        *b_MET_Type1Corr_Type6EtFraction;   //!
   TBranch        *b_MET_Type1Corr_Type7EtFraction;   //!
   TBranch        *b_MET_Type1p2Corr_et;   //!
   TBranch        *b_MET_Type1p2Corr_pt;   //!
   TBranch        *b_MET_Type1p2Corr_phi;   //!
   TBranch        *b_MET_Type1p2Corr_sumET;   //!
   TBranch        *b_MET_Type1p2Corr_MuonEtFraction;   //!
   TBranch        *b_MET_Type1p2Corr_NeutralEMFraction;   //!
   TBranch        *b_MET_Type1p2Corr_NeutralHadEtFraction;   //!
   TBranch        *b_MET_Type1p2Corr_Type6EtFraction;   //!
   TBranch        *b_MET_Type1p2Corr_Type7EtFraction;   //!
   TBranch        *b_MET_Type1CorrElectronUp_et;   //!
   TBranch        *b_MET_Type1CorrElectronDown_et;   //!
   TBranch        *b_MET_Type1CorrMuonUp_et;   //!
   TBranch        *b_MET_Type1CorrMuonDown_et;   //!
   TBranch        *b_MET_Type1CorrTauUp_et;   //!
   TBranch        *b_MET_Type1CorrTauDown_et;   //!
   TBranch        *b_MET_Type1CorrJetResUp_et;   //!
   TBranch        *b_MET_Type1CorrJetResDown_et;   //!
   TBranch        *b_MET_Type1CorrJetEnUp_et;   //!
   TBranch        *b_MET_Type1CorrJetEnDown_et;   //!
   TBranch        *b_MET_Type1CorrUnclusteredUp_et;   //!
   TBranch        *b_MET_Type1CorrUnclusteredDown_et;   //!
   TBranch        *b_MET_Type1p2CorrElectronUp_et;   //!
   TBranch        *b_MET_Type1p2CorrElectronDown_et;   //!
   TBranch        *b_MET_Type1p2CorrMuonUp_et;   //!
   TBranch        *b_MET_Type1p2CorrMuonDown_et;   //!
   TBranch        *b_MET_Type1p2CorrTauUp_et;   //!
   TBranch        *b_MET_Type1p2CorrTauDown_et;   //!
   TBranch        *b_MET_Type1p2CorrJetResUp_et;   //!
   TBranch        *b_MET_Type1p2CorrJetResDown_et;   //!
   TBranch        *b_MET_Type1p2CorrJetEnUp_et;   //!
   TBranch        *b_MET_Type1p2CorrJetEnDown_et;   //!
   TBranch        *b_MET_Type1p2CorrUnclusteredUp_et;   //!
   TBranch        *b_MET_Type1p2CorrUnclusteredDown_et;   //!
   TBranch        *b_Event_EventNumber;   //!
   TBranch        *b_Event_RunNumber;   //!
   TBranch        *b_Event_bunchCrossing;   //!
   TBranch        *b_Event_orbitNumber;   //!
   TBranch        *b_Event_luminosityBlock;   //!
   TBranch        *b_Event_isRealData;   //!
   TBranch        *b_PileupInfo_TrueNumInteractions_nm1;   //!
   TBranch        *b_PileupInfo_TrueNumInteractions_n0;   //!
   TBranch        *b_PileupInfo_TrueNumInteractions_np1;   //!
   TBranch        *b_PUWeight;   //!
   TBranch        *b_PUWeight_p5;   //!
   TBranch        *b_PUWeight_m5;   //!
   TBranch        *b_PUWeight3D;   //!
   TBranch        *b_PUWeight3D_p5;   //!
   TBranch        *b_PUWeight3D_m5;   //!
   TBranch        *b_PUWeightFineBins;   //!
   TBranch        *b_TauSpinnerWeight;   //!
   TBranch        *b_SelEffWeight;   //!
   TBranch        *b_MinVisPtFilter;   //!
   TBranch        *b_KinWeightPt;   //!
   TBranch        *b_KinWeightEta;   //!
   TBranch        *b_KinWeightMassPt;   //!
   TBranch        *b_EmbeddedWeight;   //!
   TBranch        *b_Track_p4;   //!
   TBranch        *b_Track_Poca;   //!
   TBranch        *b_Track_chi2;   //!
   TBranch        *b_Track_ndof;   //!
   TBranch        *b_Track_numberOfLostHits;   //!
   TBranch        *b_Track_numberOfValidHits;   //!
   TBranch        *b_Track_qualityMask;   //!
   TBranch        *b_Track_charge;   //!
   TBranch        *b_Track_pdgid;   //!
   TBranch        *b_Track_B;   //!
   TBranch        *b_Track_M;   //!
   TBranch        *b_Track_par;   //!
   TBranch        *b_Track_cov;   //!
   TBranch        *b_GenEventInfoProduct_signalProcessID;   //!
   TBranch        *b_GenEventInfoProduct_weight;   //!
   TBranch        *b_GenEventInfoProduct_weights;   //!
   TBranch        *b_GenEventInfoProduct_qScale;   //!
   TBranch        *b_GenEventInfoProduct_alphaQED;   //!
   TBranch        *b_GenEventInfoProduct_alphaQCD;   //!
   TBranch        *b_GenEventInfoProduct_id1;   //!
   TBranch        *b_GenEventInfoProduct_id2;   //!
   TBranch        *b_GenEventInfoProduct_x1;   //!
   TBranch        *b_GenEventInfoProduct_x2;   //!
   TBranch        *b_GenEventInfoProduct_scalePDF;   //!
   TBranch        *b_MC_p4;   //!
   TBranch        *b_MC_pdgid;   //!
   TBranch        *b_MC_charge;   //!
   TBranch        *b_MC_midx;   //!
   TBranch        *b_MC_childpdgid;   //!
   TBranch        *b_MC_childidx;   //!
   TBranch        *b_MC_status;   //!
   TBranch        *b_MCSignalParticle_p4;   //!
   TBranch        *b_MCSignalParticle_pdgid;   //!
   TBranch        *b_MCSignalParticle_childpdgid;   //!
   TBranch        *b_MCSignalParticle_charge;   //!
   TBranch        *b_MCSignalParticle_Poca;   //!
   TBranch        *b_MCSignalParticle_Tauidx;   //!
   TBranch        *b_MCTauandProd_p4;   //!
   TBranch        *b_MCTauandProd_Vertex;   //!
   TBranch        *b_MCTauandProd_pdgid;   //!
   TBranch        *b_MCTauandProd_midx;   //!
   TBranch        *b_MCTauandProd_charge;   //!
   TBranch        *b_MCTau_JAK;   //!
   TBranch        *b_MCTau_DecayBitMask;   //!
   TBranch        *b_HTLTriggerName;   //!
   TBranch        *b_TriggerAccept;   //!
   TBranch        *b_TriggerError;   //!
   TBranch        *b_TriggerWasRun;   //!
   TBranch        *b_HLTPrescale;   //!
   TBranch        *b_NHLTL1GTSeeds;   //!
   TBranch        *b_L1SEEDPrescale;   //!
   TBranch        *b_L1SEEDInvalidPrescale;   //!
   TBranch        *b_L1SEEDisTechBit;   //!
   TBranch        *b_TriggerMatchMuon;   //!
   TBranch        *b_TriggerMatchJet;   //!
   TBranch        *b_TriggerMatchTau;   //!
   TBranch        *b_HLTTrigger_objs_Pt;   //!
   TBranch        *b_HLTTrigger_objs_Eta;   //!
   TBranch        *b_HLTTrigger_objs_Phi;   //!
   TBranch        *b_HLTTrigger_objs_E;   //!
   TBranch        *b_HLTTrigger_objs_Id;   //!
   TBranch        *b_HLTTrigger_objs_trigger;   //!
   TBranch        *b_L1TriggerName;   //!
   TBranch        *b_L1TriggerDecision;   //!
   TBranch        *b_L1ErrorCode;   //!
   TBranch        *b_L1Prescale;   //!

   NtupleReader(TTree *tree=0);
   virtual ~NtupleReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual double   PFTau_Mass(unsigned int i);
};

#endif

#ifdef NtupleReader_cxx
NtupleReader::NtupleReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("t",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("t","");
      chain->Add("/net/scratch_cms/institut_3b/kargoll/TestPrecision/TauNtuple_5_3_22_patch1-Jan_28_2015/CMSSW_5_3_22_patch1/src/SkimProduction/CRAB/dy_tautau0/TauNtuple.root/t");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

NtupleReader::~NtupleReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NtupleReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NtupleReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NtupleReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   beamspot_par = 0;
   beamspot_cov = 0;
   Vtx_chi2 = 0;
   Vtx_nTrk = 0;
   Vtx_ndof = 0;
   Vtx_x = 0;
   Vtx_y = 0;
   Vtx_z = 0;
   Vtx_Cov = 0;
   Vtx_Track_idx = 0;
   Vtx_Track_Weights = 0;
   Vtx_isFake = 0;
   Vtx_TracksP4 = 0;
   Muon_p4 = 0;
   Muon_Poca = 0;
   Muon_isGlobalMuon = 0;
   Muon_isStandAloneMuon = 0;
   Muon_isTrackerMuon = 0;
   Muon_isCaloMuon = 0;
   Muon_isIsolationValid = 0;
   Muon_isQualityValid = 0;
   Muon_isTimeValid = 0;
   Muon_emEt03 = 0;
   Muon_emVetoEt03 = 0;
   Muon_hadEt03 = 0;
   Muon_hadVetoEt03 = 0;
   Muon_nJets03 = 0;
   Muon_nTracks03 = 0;
   Muon_sumPt03 = 0;
   Muon_trackerVetoPt03 = 0;
   Muon_emEt05 = 0;
   Muon_emVetoEt05 = 0;
   Muon_hadEt05 = 0;
   Muon_hadVetoEt05 = 0;
   Muon_nJets05 = 0;
   Muon_nTracks05 = 0;
   Muon_sumPt05 = 0;
   Muon_trackerVetoPt05 = 0;
   Muon_sumChargedHadronPt03 = 0;
   Muon_sumChargedParticlePt03 = 0;
   Muon_sumNeutralHadronEt03 = 0;
   Muon_sumNeutralHadronEtHighThreshold03 = 0;
   Muon_sumPhotonEt03 = 0;
   Muon_sumPhotonEtHighThreshold03 = 0;
   Muon_sumPUPt03 = 0;
   Muon_sumChargedHadronPt04 = 0;
   Muon_sumChargedParticlePt04 = 0;
   Muon_sumNeutralHadronEt04 = 0;
   Muon_sumNeutralHadronEtHighThreshold04 = 0;
   Muon_sumPhotonEt04 = 0;
   Muon_sumPhotonEtHighThreshold04 = 0;
   Muon_sumPUPt04 = 0;
   Muon_Track_idx = 0;
   Muon_hitPattern_pixelLayerwithMeas = 0;
   Muon_numberOfMatchedStations = 0;
   Muon_normChi2 = 0;
   Muon_hitPattern_numberOfValidMuonHits = 0;
   Muon_innerTrack_numberofValidHits = 0;
   Muon_numberOfMatches = 0;
   Muon_numberOfChambers = 0;
   Muon_isPFMuon = 0;
   Muon_numberofValidPixelHits = 0;
   Muon_trackerLayersWithMeasurement = 0;
   Muon_charge = 0;
   Muon_trackCharge = 0;
   Muon_pdgid = 0;
   Muon_B = 0;
   Muon_M = 0;
   Muon_par = 0;
   Muon_cov = 0;
   Electron_p4 = 0;
   Electron_Poca = 0;
   Electron_Gsf_deltaEtaEleClusterTrackAtCalo = 0;
   Electron_Gsf_deltaEtaSeedClusterTrackAtCalo = 0;
   Electron_Gsf_deltaEtaSuperClusterTrackAtVtx = 0;
   Electron_Gsf_deltaPhiEleClusterTrackAtCalo = 0;
   Electron_Gsf_deltaPhiSeedClusterTrackAtCalo = 0;
   Electron_Gsf_deltaPhiSuperClusterTrackAtVtx = 0;
   Electron_Gsf_dr03EcalRecHitSumE = 0;
   Electron_Gsf_dr03HcalDepth1TowerSumEt = 0;
   Electron_Gsf_dr03HcalDepth1TowerSumEtBc = 0;
   Electron_Gsf_dr03HcalDepth2TowerSumEt = 0;
   Electron_Gsf_dr03HcalDepth2TowerSumEtBc = 0;
   Electron_Gsf_dr03HcalTowerSumEt = 0;
   Electron_Gsf_dr03HcalTowerSumEtBc = 0;
   Electron_Gsf_dr03TkSumPt = 0;
   Electron_Gsf_passingCutBasedPreselection = 0;
   Electron_Gsf_passingMvaPreselection = 0;
   Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits = 0;
   Electron_supercluster_e = 0;
   Electron_supercluster_phi = 0;
   Electron_supercluster_eta = 0;
   Electron_supercluster_centroid_x = 0;
   Electron_supercluster_centroid_y = 0;
   Electron_supercluster_centroid_z = 0;
   Electron_Track_idx = 0;
   Electron_ecalRecHitSumEt03 = 0;
   Electron_hcalDepth1TowerSumEt03 = 0;
   Electron_hcalDepth1TowerSumEtBc03 = 0;
   Electron_hcalDepth2TowerSumEt03 = 0;
   Electron_hcalDepth2TowerSumEtBc03 = 0;
   Electron_tkSumPt03 = 0;
   Electron_ecalRecHitSumEt04 = 0;
   Electron_hcalDepth1TowerSumEt04 = 0;
   Electron_hcalDepth1TowerSumEtBc04 = 0;
   Electron_hcalDepth2TowerSumEt04 = 0;
   Electron_hcalDepth2TowerSumEtBc04 = 0;
   Electron_tkSumPt04 = 0;
   Electron_chargedHadronIso = 0;
   Electron_neutralHadronIso = 0;
   Electron_photonIso = 0;
   Electron_isoDeposits_chargedHadronIso04 = 0;
   Electron_isoDeposits_neutralHadronIso04 = 0;
   Electron_isoDeposits_photonIso04 = 0;
   Electron_isoDeposits_chargedHadronIso03 = 0;
   Electron_isoDeposits_neutralHadronIso03 = 0;
   Electron_isoDeposits_photonIso03 = 0;
   Electron_sigmaIetaIeta = 0;
   Electron_hadronicOverEm = 0;
   Electron_fbrem = 0;
   Electron_eSuperClusterOverP = 0;
   Electron_ecalEnergy = 0;
   Electron_trackMomentumAtVtx = 0;
   Electron_numberOfMissedHits = 0;
   Electron_HasMatchedConversions = 0;
   Electron_charge = 0;
   Electron_trackCharge = 0;
   Electron_pdgid = 0;
   Electron_B = 0;
   Electron_M = 0;
   Electron_par = 0;
   Electron_cov = 0;
   Electron_RegEnergy = 0;
   Electron_RegEnergyError = 0;
   Electron_Rho_kt6PFJets = 0;
   Electron_MVA_TrigNoIP_discriminator = 0;
   Electron_MVA_NonTrig_discriminator = 0;
   Electron_MVA_Trig_discriminator = 0;
   PFTau_p4 = 0;
   PFTau_Poca = 0;
   PFTau_isTightIsolation = 0;
   PFTau_isMediumIsolation = 0;
   PFTau_isLooseIsolation = 0;
   PFTau_isTightIsolationDBSumPtCorr = 0;
   PFTau_isMediumIsolationDBSumPtCorr = 0;
   PFTau_isLooseIsolationDBSumPtCorr = 0;
   PFTau_isVLooseIsolationDBSumPtCorr = 0;
   PFTau_isHPSAgainstElectronsLoose = 0;
   PFTau_isHPSAgainstElectronsMedium = 0;
   PFTau_isHPSAgainstElectronsTight = 0;
   PFTau_isHPSAgainstMuonLoose = 0;
   PFTau_isHPSAgainstMuonMedium = 0;
   PFTau_isHPSAgainstMuonTight = 0;
   PFTau_isHPSAgainstMuonLoose2 = 0;
   PFTau_isHPSAgainstMuonMedium2 = 0;
   PFTau_isHPSAgainstMuonTight2 = 0;
   PFTau_isHPSByDecayModeFinding = 0;
   PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection = 0;
   PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection = 0;
   PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection = 0;
   PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection = 0;
   PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits = 0;
   PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits = 0;
   PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits = 0;
   PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits = 0;
   PFTau_HPSPFTauDiscriminationByLooseIsolationMVA = 0;
   PFTau_HPSPFTauDiscriminationByMediumIsolationMVA = 0;
   PFTau_HPSPFTauDiscriminationByTightIsolationMVA = 0;
   PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2 = 0;
   PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2 = 0;
   PFTau_HPSPFTauDiscriminationByTightIsolationMVA2 = 0;
   PFTau_hpsDecayMode = 0;
   PFTau_Charge = 0;
   PFTau_Track_idx = 0;
   PFTau_TIP_primaryVertex_pos = 0;
   PFTau_TIP_primaryVertex_cov = 0;
   PFTau_TIP_secondaryVertex_pos = 0;
   PFTau_TIP_secondaryVertex_cov = 0;
   PFTau_TIP_secondaryVertex_vtxchi2 = 0;
   PFTau_TIP_secondaryVertex_vtxndof = 0;
   PFTau_TIP_primaryVertex_vtxchi2 = 0;
   PFTau_TIP_primaryVertex_vtxndof = 0;
   PFTau_TIP_flightLength = 0;
   PFTau_TIP_flightLengthSig = 0;
   PFTau_a1_lvp = 0;
   PFTau_a1_cov = 0;
   PFTau_a1_charge = 0;
   PFTau_a1_pdgid = 0;
   PFTau_a1_B = 0;
   PFTau_a1_M = 0;
   PFTau_daughterTracks = 0;
   PFTau_daughterTracks_cov = 0;
   PFTau_daughterTracks_charge = 0;
   PFTau_daughterTracks_pdgid = 0;
   PFTau_daughterTracks_B = 0;
   PFTau_daughterTracks_M = 0;
   PFTau_daughterTracks_poca = 0;
   PFTau_3PS_A1_LV = 0;
   PFTau_3PS_M_A1 = 0;
   PFTau_3PS_M_12 = 0;
   PFTau_3PS_M_13 = 0;
   PFTau_3PS_M_23 = 0;
   PFTau_3PS_Tau_Charge = 0;
   PFTau_3PS_LCchi2 = 0;
   PFTau_3PS_has3ProngSolution = 0;
   PFTau_3PS_Tau_LV = 0;
   PFTau_PionsP4 = 0;
   PFTau_PionsCharge = 0;
   PFTau_PiZeroP4 = 0;
   PFTau_PiZeroNumOfPhotons = 0;
   PFTau_PiZeroNumOfElectrons = 0;
   PFTau_ChargedHadronsP4 = 0;
   PFTau_ChargedHadronsCharge = 0;
   PFTau_GammaP4 = 0;
   PFTau_Photons_p4_inDR05 = 0;
   PFTau_MatchedPFJetP4 = 0;
   PFTau_MatchedPFJetGammasP4 = 0;
   PFTau_MatchedPFJetSCVariables = 0;
   PFTau_MatchedPFJetPhotonVariables = 0;
   PFTau_PhotonEnergyFraction = 0;
   PFTau_photon_hasPixelSeed = 0;
   PFTau_photon_hadronicOverEm = 0;
   PFTau_photon_sigmaIetaIeta = 0;
   PFTau_photon_trkSumPtHollowConeDR04 = 0;
   PFTau_photon_ecalRecHitSumEtConeDR04 = 0;
   PFTau_photon_hcalTowerSumEtConeDR04 = 0;
   PFTau_photon_rho = 0;
   PFJet_p4 = 0;
   PFJet_chargedEmEnergy = 0;
   PFJet_chargedHadronEnergy = 0;
   PFJet_chargedHadronMultiplicity = 0;
   PFJet_chargedMuEnergy = 0;
   PFJet_chargedMultiplicity = 0;
   PFJet_electronEnergy = 0;
   PFJet_electronMultiplicity = 0;
   PFJet_HFEMEnergy = 0;
   PFJet_HFEMMultiplicity = 0;
   PFJet_HFHadronEnergy = 0;
   PFJet_HFHadronMultiplicity = 0;
   PFJet_muonEnergy = 0;
   PFJet_muonMultiplicity = 0;
   PFJet_neutralEmEnergy = 0;
   PFJet_neutralHadronEnergy = 0;
   PFJet_neutralHadronMultiplicity = 0;
   PFJet_photonEnergy = 0;
   PFJet_photonMultiplicity = 0;
   PFJet_jetArea = 0;
   PFJet_maxDistance = 0;
   PFJet_nConstituents = 0;
   PFJet_pileup = 0;
   PFJet_etaetaMoment = 0;
   PFJet_etaphiMoment = 0;
   PFJet_Track_idx = 0;
   PFJet_MatchedHPS_idx = 0;
   PFJet_numberOfDaughters = 0;
   PFJet_chargedEmEnergyFraction = 0;
   PFJet_chargedHadronEnergyFraction = 0;
   PFJet_neutralHadronEnergyFraction = 0;
   PFJet_neutralEmEnergyFraction = 0;
   PFJet_PUJetID_discr = 0;
   PFJet_PUJetID_looseWP = 0;
   PFJet_PUJetID_mediumWP = 0;
   PFJet_PUJetID_tightWP = 0;
   PFJet_partonFlavour = 0;
   PFJet_bDiscriminator = 0;
   PFJet_TracksP4 = 0;
   PFJet_nTrk = 0;
   PFJet_JECuncertainty = 0;
   PFJet_GenJet_p4 = 0;
   PFJet_GenJet_Constituents_p4 = 0;
   PFJet_GenJetNoNu_p4 = 0;
   PFJet_GenJetNoNu_Constituents_p4 = 0;
   MET_CorrMVA_srcMuon_p4 = 0;
   MET_CorrMVA_srcElectron_p4 = 0;
   MET_CorrMVA_srcTau_p4 = 0;
   MET_CorrMVAMuTau_srcMuon_p4 = 0;
   MET_CorrMVAMuTau_srcTau_p4 = 0;
   Track_p4 = 0;
   Track_Poca = 0;
   Track_chi2 = 0;
   Track_ndof = 0;
   Track_numberOfLostHits = 0;
   Track_numberOfValidHits = 0;
   Track_qualityMask = 0;
   Track_charge = 0;
   Track_pdgid = 0;
   Track_B = 0;
   Track_M = 0;
   Track_par = 0;
   Track_cov = 0;
   GenEventInfoProduct_weights = 0;
   MC_p4 = 0;
   MC_pdgid = 0;
   MC_charge = 0;
   MC_midx = 0;
   MC_childpdgid = 0;
   MC_childidx = 0;
   MC_status = 0;
   MCSignalParticle_p4 = 0;
   MCSignalParticle_pdgid = 0;
   MCSignalParticle_childpdgid = 0;
   MCSignalParticle_charge = 0;
   MCSignalParticle_Poca = 0;
   MCSignalParticle_Tauidx = 0;
   MCTauandProd_p4 = 0;
   MCTauandProd_Vertex = 0;
   MCTauandProd_pdgid = 0;
   MCTauandProd_midx = 0;
   MCTauandProd_charge = 0;
   MCTau_JAK = 0;
   MCTau_DecayBitMask = 0;
   HTLTriggerName = 0;
   TriggerAccept = 0;
   TriggerError = 0;
   TriggerWasRun = 0;
   HLTPrescale = 0;
   NHLTL1GTSeeds = 0;
   L1SEEDPrescale = 0;
   L1SEEDInvalidPrescale = 0;
   L1SEEDisTechBit = 0;
   TriggerMatchMuon = 0;
   TriggerMatchJet = 0;
   TriggerMatchTau = 0;
   HLTTrigger_objs_Pt = 0;
   HLTTrigger_objs_Eta = 0;
   HLTTrigger_objs_Phi = 0;
   HLTTrigger_objs_E = 0;
   HLTTrigger_objs_Id = 0;
   HLTTrigger_objs_trigger = 0;
   L1TriggerName = 0;
   L1TriggerDecision = 0;
   L1ErrorCode = 0;
   L1Prescale = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("DataMC_Type", &DataMC_Type, &b_DataMC_Type);
   fChain->SetBranchAddress("beamspot_par", &beamspot_par, &b_beamspot_par);
   fChain->SetBranchAddress("beamspot_cov", &beamspot_cov, &b_beamspot_cov);
   fChain->SetBranchAddress("beamspot_emittanceX", &beamspot_emittanceX, &b_beamspot_emittanceX);
   fChain->SetBranchAddress("beamspot_emittanceY", &beamspot_emittanceY, &b_beamspot_emittanceY);
   fChain->SetBranchAddress("beamspot_betaStar", &beamspot_betaStar, &b_beamspot_betaStar);
   fChain->SetBranchAddress("Vtx_chi2", &Vtx_chi2, &b_Vtx_chi2);
   fChain->SetBranchAddress("Vtx_nTrk", &Vtx_nTrk, &b_Vtx_nTrk);
   fChain->SetBranchAddress("Vtx_ndof", &Vtx_ndof, &b_Vtx_ndof);
   fChain->SetBranchAddress("Vtx_x", &Vtx_x, &b_Vtx_x);
   fChain->SetBranchAddress("Vtx_y", &Vtx_y, &b_Vtx_y);
   fChain->SetBranchAddress("Vtx_z", &Vtx_z, &b_Vtx_z);
   fChain->SetBranchAddress("Vtx_Cov", &Vtx_Cov, &b_Vtx_Cov);
   fChain->SetBranchAddress("Vtx_Track_idx", &Vtx_Track_idx, &b_Vtx_Track_idx);
   fChain->SetBranchAddress("Vtx_Track_Weights", &Vtx_Track_Weights, &b_Vtx_Track_Weights);
   fChain->SetBranchAddress("Vtx_isFake", &Vtx_isFake, &b_Vtx_isFake);
   fChain->SetBranchAddress("Vtx_TracksP4", &Vtx_TracksP4, &b_Vtx_TracksP4);
   fChain->SetBranchAddress("isPatMuon", &isPatMuon, &b_isPatMuon);
   fChain->SetBranchAddress("Muon_p4", &Muon_p4, &b_Muon_p4);
   fChain->SetBranchAddress("Muon_Poca", &Muon_Poca, &b_Muon_Poca);
   fChain->SetBranchAddress("Muon_isGlobalMuon", &Muon_isGlobalMuon, &b_Muon_isGlobalMuon);
   fChain->SetBranchAddress("Muon_isStandAloneMuon", &Muon_isStandAloneMuon, &b_Muon_isStandAloneMuon);
   fChain->SetBranchAddress("Muon_isTrackerMuon", &Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
   fChain->SetBranchAddress("Muon_isCaloMuon", &Muon_isCaloMuon, &b_Muon_isCaloMuon);
   fChain->SetBranchAddress("Muon_isIsolationValid", &Muon_isIsolationValid, &b_Muon_isIsolationValid);
   fChain->SetBranchAddress("Muon_isQualityValid", &Muon_isQualityValid, &b_Muon_isQualityValid);
   fChain->SetBranchAddress("Muon_isTimeValid", &Muon_isTimeValid, &b_Muon_isTimeValid);
   fChain->SetBranchAddress("Muon_emEt03", &Muon_emEt03, &b_Muon_emEt03);
   fChain->SetBranchAddress("Muon_emVetoEt03", &Muon_emVetoEt03, &b_Muon_emVetoEt03);
   fChain->SetBranchAddress("Muon_hadEt03", &Muon_hadEt03, &b_Muon_hadEt03);
   fChain->SetBranchAddress("Muon_hadVetoEt03", &Muon_hadVetoEt03, &b_Muon_hadVetoEt03);
   fChain->SetBranchAddress("Muon_nJets03", &Muon_nJets03, &b_Muon_nJets03);
   fChain->SetBranchAddress("Muon_nTracks03", &Muon_nTracks03, &b_Muon_nTracks03);
   fChain->SetBranchAddress("Muon_sumPt03", &Muon_sumPt03, &b_Muon_sumPt03);
   fChain->SetBranchAddress("Muon_trackerVetoPt03", &Muon_trackerVetoPt03, &b_Muon_trackerVetoPt03);
   fChain->SetBranchAddress("Muon_emEt05", &Muon_emEt05, &b_Muon_emEt05);
   fChain->SetBranchAddress("Muon_emVetoEt05", &Muon_emVetoEt05, &b_Muon_emVetoEt05);
   fChain->SetBranchAddress("Muon_hadEt05", &Muon_hadEt05, &b_Muon_hadEt05);
   fChain->SetBranchAddress("Muon_hadVetoEt05", &Muon_hadVetoEt05, &b_Muon_hadVetoEt05);
   fChain->SetBranchAddress("Muon_nJets05", &Muon_nJets05, &b_Muon_nJets05);
   fChain->SetBranchAddress("Muon_nTracks05", &Muon_nTracks05, &b_Muon_nTracks05);
   fChain->SetBranchAddress("Muon_sumPt05", &Muon_sumPt05, &b_Muon_sumPt05);
   fChain->SetBranchAddress("Muon_trackerVetoPt05", &Muon_trackerVetoPt05, &b_Muon_trackerVetoPt05);
   fChain->SetBranchAddress("Muon_sumChargedHadronPt03", &Muon_sumChargedHadronPt03, &b_Muon_sumChargedHadronPt03);
   fChain->SetBranchAddress("Muon_sumChargedParticlePt03", &Muon_sumChargedParticlePt03, &b_Muon_sumChargedParticlePt03);
   fChain->SetBranchAddress("Muon_sumNeutralHadronEt03", &Muon_sumNeutralHadronEt03, &b_Muon_sumNeutralHadronEt03);
   fChain->SetBranchAddress("Muon_sumNeutralHadronEtHighThreshold03", &Muon_sumNeutralHadronEtHighThreshold03, &b_Muon_sumNeutralHadronEtHighThreshold03);
   fChain->SetBranchAddress("Muon_sumPhotonEt03", &Muon_sumPhotonEt03, &b_Muon_sumPhotonEt03);
   fChain->SetBranchAddress("Muon_sumPhotonEtHighThreshold03", &Muon_sumPhotonEtHighThreshold03, &b_Muon_sumPhotonEtHighThreshold03);
   fChain->SetBranchAddress("Muon_sumPUPt03", &Muon_sumPUPt03, &b_Muon_sumPUPt03);
   fChain->SetBranchAddress("Muon_sumChargedHadronPt04", &Muon_sumChargedHadronPt04, &b_Muon_sumChargedHadronPt04);
   fChain->SetBranchAddress("Muon_sumChargedParticlePt04", &Muon_sumChargedParticlePt04, &b_Muon_sumChargedParticlePt04);
   fChain->SetBranchAddress("Muon_sumNeutralHadronEt04", &Muon_sumNeutralHadronEt04, &b_Muon_sumNeutralHadronEt04);
   fChain->SetBranchAddress("Muon_sumNeutralHadronEtHighThreshold04", &Muon_sumNeutralHadronEtHighThreshold04, &b_Muon_sumNeutralHadronEtHighThreshold04);
   fChain->SetBranchAddress("Muon_sumPhotonEt04", &Muon_sumPhotonEt04, &b_Muon_sumPhotonEt04);
   fChain->SetBranchAddress("Muon_sumPhotonEtHighThreshold04", &Muon_sumPhotonEtHighThreshold04, &b_Muon_sumPhotonEtHighThreshold04);
   fChain->SetBranchAddress("Muon_sumPUPt04", &Muon_sumPUPt04, &b_Muon_sumPUPt04);
   fChain->SetBranchAddress("Muon_Track_idx", &Muon_Track_idx, &b_Muon_Track_idx);
   fChain->SetBranchAddress("Muon_hitPattern_pixelLayerwithMeas", &Muon_hitPattern_pixelLayerwithMeas, &b_Muon_hitPattern_pixelLayerwithMeas);
   fChain->SetBranchAddress("Muon_numberOfMatchedStations", &Muon_numberOfMatchedStations, &b_Muon_numberOfMatchedStations);
   fChain->SetBranchAddress("Muon_normChi2", &Muon_normChi2, &b_Muon_normChi2);
   fChain->SetBranchAddress("Muon_hitPattern_numberOfValidMuonHits", &Muon_hitPattern_numberOfValidMuonHits, &b_Muon_hitPattern_numberOfValidMuonHits);
   fChain->SetBranchAddress("Muon_innerTrack_numberofValidHits", &Muon_innerTrack_numberofValidHits, &b_Muon_innerTrack_numberofValidHits);
   fChain->SetBranchAddress("Muon_numberOfMatches", &Muon_numberOfMatches, &b_Muon_numberOfMatches);
   fChain->SetBranchAddress("Muon_numberOfChambers", &Muon_numberOfChambers, &b_Muon_numberOfChambers);
   fChain->SetBranchAddress("Muon_isPFMuon", &Muon_isPFMuon, &b_Muon_isPFMuon);
   fChain->SetBranchAddress("Muon_numberofValidPixelHits", &Muon_numberofValidPixelHits, &b_Muon_numberofValidPixelHits);
   fChain->SetBranchAddress("Muon_trackerLayersWithMeasurement", &Muon_trackerLayersWithMeasurement, &b_Muon_trackerLayersWithMeasurement);
   fChain->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_trackCharge", &Muon_trackCharge, &b_Muon_trackCharge);
   fChain->SetBranchAddress("Muon_pdgid", &Muon_pdgid, &b_Muon_pdgid);
   fChain->SetBranchAddress("Muon_B", &Muon_B, &b_Muon_B);
   fChain->SetBranchAddress("Muon_M", &Muon_M, &b_Muon_M);
   fChain->SetBranchAddress("Muon_par", &Muon_par, &b_Muon_par);
   fChain->SetBranchAddress("Muon_cov", &Muon_cov, &b_Muon_cov);
   fChain->SetBranchAddress("isPatElectron", &isPatElectron, &b_isPatElectron);
   fChain->SetBranchAddress("Electron_p4", &Electron_p4, &b_Electron_p4);
   fChain->SetBranchAddress("Electron_Poca", &Electron_Poca, &b_Electron_Poca);
   fChain->SetBranchAddress("Electron_Gsf_deltaEtaEleClusterTrackAtCalo", &Electron_Gsf_deltaEtaEleClusterTrackAtCalo, &b_Electron_Gsf_deltaEtaEleClusterTrackAtCalo);
   fChain->SetBranchAddress("Electron_Gsf_deltaEtaSeedClusterTrackAtCalo", &Electron_Gsf_deltaEtaSeedClusterTrackAtCalo, &b_Electron_Gsf_deltaEtaSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("Electron_Gsf_deltaEtaSuperClusterTrackAtVtx", &Electron_Gsf_deltaEtaSuperClusterTrackAtVtx, &b_Electron_Gsf_deltaEtaSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("Electron_Gsf_deltaPhiEleClusterTrackAtCalo", &Electron_Gsf_deltaPhiEleClusterTrackAtCalo, &b_Electron_Gsf_deltaPhiEleClusterTrackAtCalo);
   fChain->SetBranchAddress("Electron_Gsf_deltaPhiSeedClusterTrackAtCalo", &Electron_Gsf_deltaPhiSeedClusterTrackAtCalo, &b_Electron_Gsf_deltaPhiSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("Electron_Gsf_deltaPhiSuperClusterTrackAtVtx", &Electron_Gsf_deltaPhiSuperClusterTrackAtVtx, &b_Electron_Gsf_deltaPhiSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("Electron_Gsf_dr03EcalRecHitSumE", &Electron_Gsf_dr03EcalRecHitSumE, &b_Electron_Gsf_dr03EcalRecHitSumE);
   fChain->SetBranchAddress("Electron_Gsf_dr03HcalDepth1TowerSumEt", &Electron_Gsf_dr03HcalDepth1TowerSumEt, &b_Electron_Gsf_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_Gsf_dr03HcalDepth1TowerSumEtBc", &Electron_Gsf_dr03HcalDepth1TowerSumEtBc, &b_Electron_Gsf_dr03HcalDepth1TowerSumEtBc);
   fChain->SetBranchAddress("Electron_Gsf_dr03HcalDepth2TowerSumEt", &Electron_Gsf_dr03HcalDepth2TowerSumEt, &b_Electron_Gsf_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("Electron_Gsf_dr03HcalDepth2TowerSumEtBc", &Electron_Gsf_dr03HcalDepth2TowerSumEtBc, &b_Electron_Gsf_dr03HcalDepth2TowerSumEtBc);
   fChain->SetBranchAddress("Electron_Gsf_dr03HcalTowerSumEt", &Electron_Gsf_dr03HcalTowerSumEt, &b_Electron_Gsf_dr03HcalTowerSumEt);
   fChain->SetBranchAddress("Electron_Gsf_dr03HcalTowerSumEtBc", &Electron_Gsf_dr03HcalTowerSumEtBc, &b_Electron_Gsf_dr03HcalTowerSumEtBc);
   fChain->SetBranchAddress("Electron_Gsf_dr03TkSumPt", &Electron_Gsf_dr03TkSumPt, &b_Electron_Gsf_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_Gsf_passingCutBasedPreselection", &Electron_Gsf_passingCutBasedPreselection, &b_Electron_Gsf_passingCutBasedPreselection);
   fChain->SetBranchAddress("Electron_Gsf_passingMvaPreselection", &Electron_Gsf_passingMvaPreselection, &b_Electron_Gsf_passingMvaPreselection);
   fChain->SetBranchAddress("Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits", &Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits, &b_Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits);
   fChain->SetBranchAddress("Electron_supercluster_e", &Electron_supercluster_e, &b_Electron_supercluster_e);
   fChain->SetBranchAddress("Electron_supercluster_phi", &Electron_supercluster_phi, &b_Electron_supercluster_phi);
   fChain->SetBranchAddress("Electron_supercluster_eta", &Electron_supercluster_eta, &b_Electron_supercluster_eta);
   fChain->SetBranchAddress("Electron_supercluster_centroid_x", &Electron_supercluster_centroid_x, &b_Electron_supercluster_centroid_x);
   fChain->SetBranchAddress("Electron_supercluster_centroid_y", &Electron_supercluster_centroid_y, &b_Electron_supercluster_centroid_y);
   fChain->SetBranchAddress("Electron_supercluster_centroid_z", &Electron_supercluster_centroid_z, &b_Electron_supercluster_centroid_z);
   fChain->SetBranchAddress("Electron_Track_idx", &Electron_Track_idx, &b_Electron_Track_idx);
   fChain->SetBranchAddress("Electron_ecalRecHitSumEt03", &Electron_ecalRecHitSumEt03, &b_Electron_ecalRecHitSumEt03);
   fChain->SetBranchAddress("Electron_hcalDepth1TowerSumEt03", &Electron_hcalDepth1TowerSumEt03, &b_Electron_hcalDepth1TowerSumEt03);
   fChain->SetBranchAddress("Electron_hcalDepth1TowerSumEtBc03", &Electron_hcalDepth1TowerSumEtBc03, &b_Electron_hcalDepth1TowerSumEtBc03);
   fChain->SetBranchAddress("Electron_hcalDepth2TowerSumEt03", &Electron_hcalDepth2TowerSumEt03, &b_Electron_hcalDepth2TowerSumEt03);
   fChain->SetBranchAddress("Electron_hcalDepth2TowerSumEtBc03", &Electron_hcalDepth2TowerSumEtBc03, &b_Electron_hcalDepth2TowerSumEtBc03);
   fChain->SetBranchAddress("Electron_tkSumPt03", &Electron_tkSumPt03, &b_Electron_tkSumPt03);
   fChain->SetBranchAddress("Electron_ecalRecHitSumEt04", &Electron_ecalRecHitSumEt04, &b_Electron_ecalRecHitSumEt04);
   fChain->SetBranchAddress("Electron_hcalDepth1TowerSumEt04", &Electron_hcalDepth1TowerSumEt04, &b_Electron_hcalDepth1TowerSumEt04);
   fChain->SetBranchAddress("Electron_hcalDepth1TowerSumEtBc04", &Electron_hcalDepth1TowerSumEtBc04, &b_Electron_hcalDepth1TowerSumEtBc04);
   fChain->SetBranchAddress("Electron_hcalDepth2TowerSumEt04", &Electron_hcalDepth2TowerSumEt04, &b_Electron_hcalDepth2TowerSumEt04);
   fChain->SetBranchAddress("Electron_hcalDepth2TowerSumEtBc04", &Electron_hcalDepth2TowerSumEtBc04, &b_Electron_hcalDepth2TowerSumEtBc04);
   fChain->SetBranchAddress("Electron_tkSumPt04", &Electron_tkSumPt04, &b_Electron_tkSumPt04);
   fChain->SetBranchAddress("Electron_chargedHadronIso", &Electron_chargedHadronIso, &b_Electron_chargedHadronIso);
   fChain->SetBranchAddress("Electron_neutralHadronIso", &Electron_neutralHadronIso, &b_Electron_neutralHadronIso);
   fChain->SetBranchAddress("Electron_photonIso", &Electron_photonIso, &b_Electron_photonIso);
   fChain->SetBranchAddress("Electron_isoDeposits_chargedHadronIso04", &Electron_isoDeposits_chargedHadronIso04, &b_Electron_isoDeposits_chargedHadronIso04);
   fChain->SetBranchAddress("Electron_isoDeposits_neutralHadronIso04", &Electron_isoDeposits_neutralHadronIso04, &b_Electron_isoDeposits_neutralHadronIso04);
   fChain->SetBranchAddress("Electron_isoDeposits_photonIso04", &Electron_isoDeposits_photonIso04, &b_Electron_isoDeposits_photonIso04);
   fChain->SetBranchAddress("Electron_isoDeposits_chargedHadronIso03", &Electron_isoDeposits_chargedHadronIso03, &b_Electron_isoDeposits_chargedHadronIso03);
   fChain->SetBranchAddress("Electron_isoDeposits_neutralHadronIso03", &Electron_isoDeposits_neutralHadronIso03, &b_Electron_isoDeposits_neutralHadronIso03);
   fChain->SetBranchAddress("Electron_isoDeposits_photonIso03", &Electron_isoDeposits_photonIso03, &b_Electron_isoDeposits_photonIso03);
   fChain->SetBranchAddress("Electron_sigmaIetaIeta", &Electron_sigmaIetaIeta, &b_Electron_sigmaIetaIeta);
   fChain->SetBranchAddress("Electron_hadronicOverEm", &Electron_hadronicOverEm, &b_Electron_hadronicOverEm);
   fChain->SetBranchAddress("Electron_fbrem", &Electron_fbrem, &b_Electron_fbrem);
   fChain->SetBranchAddress("Electron_eSuperClusterOverP", &Electron_eSuperClusterOverP, &b_Electron_eSuperClusterOverP);
   fChain->SetBranchAddress("Electron_ecalEnergy", &Electron_ecalEnergy, &b_Electron_ecalEnergy);
   fChain->SetBranchAddress("Electron_trackMomentumAtVtx", &Electron_trackMomentumAtVtx, &b_Electron_trackMomentumAtVtx);
   fChain->SetBranchAddress("Electron_numberOfMissedHits", &Electron_numberOfMissedHits, &b_Electron_numberOfMissedHits);
   fChain->SetBranchAddress("Electron_HasMatchedConversions", &Electron_HasMatchedConversions, &b_Electron_HasMatchedConversions);
   fChain->SetBranchAddress("RhoIsolationAllInputTags", &RhoIsolationAllInputTags, &b_RhoIsolationAllInputTags);
   fChain->SetBranchAddress("Electron_charge", &Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_trackCharge", &Electron_trackCharge, &b_Electron_trackCharge);
   fChain->SetBranchAddress("Electron_pdgid", &Electron_pdgid, &b_Electron_pdgid);
   fChain->SetBranchAddress("Electron_B", &Electron_B, &b_Electron_B);
   fChain->SetBranchAddress("Electron_M", &Electron_M, &b_Electron_M);
   fChain->SetBranchAddress("Electron_par", &Electron_par, &b_Electron_par);
   fChain->SetBranchAddress("Electron_cov", &Electron_cov, &b_Electron_cov);
   fChain->SetBranchAddress("Electron_RegEnergy", &Electron_RegEnergy, &b_Electron_RegEnergy);
   fChain->SetBranchAddress("Electron_RegEnergyError", &Electron_RegEnergyError, &b_Electron_RegEnergyError);
   fChain->SetBranchAddress("Electron_Rho_kt6PFJets", &Electron_Rho_kt6PFJets, &b_Electron_Rho_kt6PFJets);
   fChain->SetBranchAddress("Electron_MVA_TrigNoIP_discriminator", &Electron_MVA_TrigNoIP_discriminator, &b_Electron_MVA_TrigNoIP_discriminator);
   fChain->SetBranchAddress("Electron_MVA_NonTrig_discriminator", &Electron_MVA_NonTrig_discriminator, &b_Electron_MVA_NonTrig_discriminator);
   fChain->SetBranchAddress("Electron_MVA_Trig_discriminator", &Electron_MVA_Trig_discriminator, &b_Electron_MVA_Trig_discriminator);
   fChain->SetBranchAddress("PFTau_p4", &PFTau_p4, &b_PFTau_p4);
   fChain->SetBranchAddress("PFTau_Poca", &PFTau_Poca, &b_PFTau_Poca);
   fChain->SetBranchAddress("PFTau_isTightIsolation", &PFTau_isTightIsolation, &b_PFTau_isTightIsolation);
   fChain->SetBranchAddress("PFTau_isMediumIsolation", &PFTau_isMediumIsolation, &b_PFTau_isMediumIsolation);
   fChain->SetBranchAddress("PFTau_isLooseIsolation", &PFTau_isLooseIsolation, &b_PFTau_isLooseIsolation);
   fChain->SetBranchAddress("PFTau_isTightIsolationDBSumPtCorr", &PFTau_isTightIsolationDBSumPtCorr, &b_PFTau_isTightIsolationDBSumPtCorr);
   fChain->SetBranchAddress("PFTau_isMediumIsolationDBSumPtCorr", &PFTau_isMediumIsolationDBSumPtCorr, &b_PFTau_isMediumIsolationDBSumPtCorr);
   fChain->SetBranchAddress("PFTau_isLooseIsolationDBSumPtCorr", &PFTau_isLooseIsolationDBSumPtCorr, &b_PFTau_isLooseIsolationDBSumPtCorr);
   fChain->SetBranchAddress("PFTau_isVLooseIsolationDBSumPtCorr", &PFTau_isVLooseIsolationDBSumPtCorr, &b_PFTau_isVLooseIsolationDBSumPtCorr);
   fChain->SetBranchAddress("PFTau_isHPSAgainstElectronsLoose", &PFTau_isHPSAgainstElectronsLoose, &b_PFTau_isHPSAgainstElectronsLoose);
   fChain->SetBranchAddress("PFTau_isHPSAgainstElectronsMedium", &PFTau_isHPSAgainstElectronsMedium, &b_PFTau_isHPSAgainstElectronsMedium);
   fChain->SetBranchAddress("PFTau_isHPSAgainstElectronsTight", &PFTau_isHPSAgainstElectronsTight, &b_PFTau_isHPSAgainstElectronsTight);
   fChain->SetBranchAddress("PFTau_isHPSAgainstMuonLoose", &PFTau_isHPSAgainstMuonLoose, &b_PFTau_isHPSAgainstMuonLoose);
   fChain->SetBranchAddress("PFTau_isHPSAgainstMuonMedium", &PFTau_isHPSAgainstMuonMedium, &b_PFTau_isHPSAgainstMuonMedium);
   fChain->SetBranchAddress("PFTau_isHPSAgainstMuonTight", &PFTau_isHPSAgainstMuonTight, &b_PFTau_isHPSAgainstMuonTight);
   fChain->SetBranchAddress("PFTau_isHPSAgainstMuonLoose2", &PFTau_isHPSAgainstMuonLoose2, &b_PFTau_isHPSAgainstMuonLoose2);
   fChain->SetBranchAddress("PFTau_isHPSAgainstMuonMedium2", &PFTau_isHPSAgainstMuonMedium2, &b_PFTau_isHPSAgainstMuonMedium2);
   fChain->SetBranchAddress("PFTau_isHPSAgainstMuonTight2", &PFTau_isHPSAgainstMuonTight2, &b_PFTau_isHPSAgainstMuonTight2);
   fChain->SetBranchAddress("PFTau_isHPSByDecayModeFinding", &PFTau_isHPSByDecayModeFinding, &b_PFTau_isHPSByDecayModeFinding);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection", &PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection, &b_PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection", &PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection, &b_PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection", &PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection, &b_PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection", &PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection, &b_PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits", &PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits, &b_PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits", &PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits, &b_PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", &PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits, &b_PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits", &PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits, &b_PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByLooseIsolationMVA", &PFTau_HPSPFTauDiscriminationByLooseIsolationMVA, &b_PFTau_HPSPFTauDiscriminationByLooseIsolationMVA);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByMediumIsolationMVA", &PFTau_HPSPFTauDiscriminationByMediumIsolationMVA, &b_PFTau_HPSPFTauDiscriminationByMediumIsolationMVA);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByTightIsolationMVA", &PFTau_HPSPFTauDiscriminationByTightIsolationMVA, &b_PFTau_HPSPFTauDiscriminationByTightIsolationMVA);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2", &PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2, &b_PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2", &PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2, &b_PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2);
   fChain->SetBranchAddress("PFTau_HPSPFTauDiscriminationByTightIsolationMVA2", &PFTau_HPSPFTauDiscriminationByTightIsolationMVA2, &b_PFTau_HPSPFTauDiscriminationByTightIsolationMVA2);
   fChain->SetBranchAddress("PFTau_hpsDecayMode", &PFTau_hpsDecayMode, &b_PFTau_hpsDecayMode);
   fChain->SetBranchAddress("PFTau_Charge", &PFTau_Charge, &b_PFTau_Charge);
   fChain->SetBranchAddress("PFTau_Track_idx", &PFTau_Track_idx, &b_PFTau_Track_idx);
   fChain->SetBranchAddress("PFTau_TIP_primaryVertex_pos", &PFTau_TIP_primaryVertex_pos, &b_PFTau_TIP_primaryVertex_pos);
   fChain->SetBranchAddress("PFTau_TIP_primaryVertex_cov", &PFTau_TIP_primaryVertex_cov, &b_PFTau_TIP_primaryVertex_cov);
   fChain->SetBranchAddress("PFTau_TIP_secondaryVertex_pos", &PFTau_TIP_secondaryVertex_pos, &b_PFTau_TIP_secondaryVertex_pos);
   fChain->SetBranchAddress("PFTau_TIP_secondaryVertex_cov", &PFTau_TIP_secondaryVertex_cov, &b_PFTau_TIP_secondaryVertex_cov);
   fChain->SetBranchAddress("PFTau_TIP_secondaryVertex_vtxchi2", &PFTau_TIP_secondaryVertex_vtxchi2, &b_PFTau_TIP_secondaryVertex_vtxchi2);
   fChain->SetBranchAddress("PFTau_TIP_secondaryVertex_vtxndof", &PFTau_TIP_secondaryVertex_vtxndof, &b_PFTau_TIP_secondaryVertex_vtxndof);
   fChain->SetBranchAddress("PFTau_TIP_primaryVertex_vtxchi2", &PFTau_TIP_primaryVertex_vtxchi2, &b_PFTau_TIP_primaryVertex_vtxchi2);
   fChain->SetBranchAddress("PFTau_TIP_primaryVertex_vtxndof", &PFTau_TIP_primaryVertex_vtxndof, &b_PFTau_TIP_primaryVertex_vtxndof);
   fChain->SetBranchAddress("PFTau_TIP_flightLength", &PFTau_TIP_flightLength, &b_PFTau_TIP_flightLength);
   fChain->SetBranchAddress("PFTau_TIP_flightLengthSig", &PFTau_TIP_flightLengthSig, &b_PFTau_TIP_flightLengthSig);
   fChain->SetBranchAddress("PFTau_a1_lvp", &PFTau_a1_lvp, &b_PFTau_a1_lvp);
   fChain->SetBranchAddress("PFTau_a1_cov", &PFTau_a1_cov, &b_PFTau_a1_cov);
   fChain->SetBranchAddress("PFTau_a1_charge", &PFTau_a1_charge, &b_PFTau_a1_charge);
   fChain->SetBranchAddress("PFTau_a1_pdgid", &PFTau_a1_pdgid, &b_PFTau_a1_pdgid);
   fChain->SetBranchAddress("PFTau_a1_B", &PFTau_a1_B, &b_PFTau_a1_B);
   fChain->SetBranchAddress("PFTau_a1_M", &PFTau_a1_M, &b_PFTau_a1_M);
   fChain->SetBranchAddress("PFTau_daughterTracks", &PFTau_daughterTracks, &b_PFTau_daughterTracks);
   fChain->SetBranchAddress("PFTau_daughterTracks_cov", &PFTau_daughterTracks_cov, &b_PFTau_daughterTracks_cov);
   fChain->SetBranchAddress("PFTau_daughterTracks_charge", &PFTau_daughterTracks_charge, &b_PFTau_daughterTracks_charge);
   fChain->SetBranchAddress("PFTau_daughterTracks_pdgid", &PFTau_daughterTracks_pdgid, &b_PFTau_daughterTracks_pdgid);
   fChain->SetBranchAddress("PFTau_daughterTracks_B", &PFTau_daughterTracks_B, &b_PFTau_daughterTracks_B);
   fChain->SetBranchAddress("PFTau_daughterTracks_M", &PFTau_daughterTracks_M, &b_PFTau_daughterTracks_M);
   fChain->SetBranchAddress("PFTau_daughterTracks_poca", &PFTau_daughterTracks_poca, &b_PFTau_daughterTracks_poca);
   fChain->SetBranchAddress("PFTau_3PS_A1_LV", &PFTau_3PS_A1_LV, &b_PFTau_3PS_A1_LV);
   fChain->SetBranchAddress("PFTau_3PS_M_A1", &PFTau_3PS_M_A1, &b_PFTau_3PS_M_A1);
   fChain->SetBranchAddress("PFTau_3PS_M_12", &PFTau_3PS_M_12, &b_PFTau_3PS_M_12);
   fChain->SetBranchAddress("PFTau_3PS_M_13", &PFTau_3PS_M_13, &b_PFTau_3PS_M_13);
   fChain->SetBranchAddress("PFTau_3PS_M_23", &PFTau_3PS_M_23, &b_PFTau_3PS_M_23);
   fChain->SetBranchAddress("PFTau_3PS_Tau_Charge", &PFTau_3PS_Tau_Charge, &b_PFTau_3PS_Tau_Charge);
   fChain->SetBranchAddress("PFTau_3PS_LCchi2", &PFTau_3PS_LCchi2, &b_PFTau_3PS_LCchi2);
   fChain->SetBranchAddress("PFTau_3PS_has3ProngSolution", &PFTau_3PS_has3ProngSolution, &b_PFTau_3PS_has3ProngSolution);
   fChain->SetBranchAddress("PFTau_3PS_Tau_LV", &PFTau_3PS_Tau_LV, &b_PFTau_3PS_Tau_LV);
   fChain->SetBranchAddress("PFTau_PionsP4", &PFTau_PionsP4, &b_PFTau_PionsP4);
   fChain->SetBranchAddress("PFTau_PionsCharge", &PFTau_PionsCharge, &b_PFTau_PionsCharge);
   fChain->SetBranchAddress("PFTau_PiZeroP4", &PFTau_PiZeroP4, &b_PFTau_PiZeroP4);
   fChain->SetBranchAddress("PFTau_PiZeroNumOfPhotons", &PFTau_PiZeroNumOfPhotons, &b_PFTau_PiZeroNumOfPhotons);
   fChain->SetBranchAddress("PFTau_PiZeroNumOfElectrons", &PFTau_PiZeroNumOfElectrons, &b_PFTau_PiZeroNumOfElectrons);
   fChain->SetBranchAddress("PFTau_ChargedHadronsP4", &PFTau_ChargedHadronsP4, &b_PFTau_ChargedHadronsP4);
   fChain->SetBranchAddress("PFTau_ChargedHadronsCharge", &PFTau_ChargedHadronsCharge, &b_PFTau_ChargedHadronsCharge);
   fChain->SetBranchAddress("PFTau_GammaP4", &PFTau_GammaP4, &b_PFTau_GammaP4);
   fChain->SetBranchAddress("PFTau_Photons_p4_inDR05", &PFTau_Photons_p4_inDR05, &b_PFTau_Photons_p4_inDR05);
   fChain->SetBranchAddress("PFTau_MatchedPFJetP4", &PFTau_MatchedPFJetP4, &b_PFTau_MatchedPFJetP4);
   fChain->SetBranchAddress("PFTau_MatchedPFJetGammasP4", &PFTau_MatchedPFJetGammasP4, &b_PFTau_MatchedPFJetGammasP4);
   fChain->SetBranchAddress("PFTau_MatchedPFJetSCVariables", &PFTau_MatchedPFJetSCVariables, &b_PFTau_MatchedPFJetSCVariables);
   fChain->SetBranchAddress("PFTau_MatchedPFJetPhotonVariables", &PFTau_MatchedPFJetPhotonVariables, &b_PFTau_MatchedPFJetPhotonVariables);
   fChain->SetBranchAddress("PFTau_PhotonEnergyFraction", &PFTau_PhotonEnergyFraction, &b_PFTau_PhotonEnergyFraction);
   fChain->SetBranchAddress("PFTau_photon_hasPixelSeed", &PFTau_photon_hasPixelSeed, &b_PFTau_photon_hasPixelSeed);
   fChain->SetBranchAddress("PFTau_photon_hadronicOverEm", &PFTau_photon_hadronicOverEm, &b_PFTau_photon_hadronicOverEm);
   fChain->SetBranchAddress("PFTau_photon_sigmaIetaIeta", &PFTau_photon_sigmaIetaIeta, &b_PFTau_photon_sigmaIetaIeta);
   fChain->SetBranchAddress("PFTau_photon_trkSumPtHollowConeDR04", &PFTau_photon_trkSumPtHollowConeDR04, &b_PFTau_photon_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("PFTau_photon_ecalRecHitSumEtConeDR04", &PFTau_photon_ecalRecHitSumEtConeDR04, &b_PFTau_photon_ecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("PFTau_photon_hcalTowerSumEtConeDR04", &PFTau_photon_hcalTowerSumEtConeDR04, &b_PFTau_photon_hcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("PFTau_photon_rho", &PFTau_photon_rho, &b_PFTau_photon_rho);
   fChain->SetBranchAddress("isPatJet", &isPatJet, &b_isPatJet);
   fChain->SetBranchAddress("PFJet_p4", &PFJet_p4, &b_PFJet_p4);
   fChain->SetBranchAddress("PFJet_chargedEmEnergy", &PFJet_chargedEmEnergy, &b_PFJet_chargedEmEnergy);
   fChain->SetBranchAddress("PFJet_chargedHadronEnergy", &PFJet_chargedHadronEnergy, &b_PFJet_chargedHadronEnergy);
   fChain->SetBranchAddress("PFJet_chargedHadronMultiplicity", &PFJet_chargedHadronMultiplicity, &b_PFJet_chargedHadronMultiplicity);
   fChain->SetBranchAddress("PFJet_chargedMuEnergy", &PFJet_chargedMuEnergy, &b_PFJet_chargedMuEnergy);
   fChain->SetBranchAddress("PFJet_chargedMultiplicity", &PFJet_chargedMultiplicity, &b_PFJet_chargedMultiplicity);
   fChain->SetBranchAddress("PFJet_electronEnergy", &PFJet_electronEnergy, &b_PFJet_electronEnergy);
   fChain->SetBranchAddress("PFJet_electronMultiplicity", &PFJet_electronMultiplicity, &b_PFJet_electronMultiplicity);
   fChain->SetBranchAddress("PFJet_HFEMEnergy", &PFJet_HFEMEnergy, &b_PFJet_HFEMEnergy);
   fChain->SetBranchAddress("PFJet_HFEMMultiplicity", &PFJet_HFEMMultiplicity, &b_PFJet_HFEMMultiplicity);
   fChain->SetBranchAddress("PFJet_HFHadronEnergy", &PFJet_HFHadronEnergy, &b_PFJet_HFHadronEnergy);
   fChain->SetBranchAddress("PFJet_HFHadronMultiplicity", &PFJet_HFHadronMultiplicity, &b_PFJet_HFHadronMultiplicity);
   fChain->SetBranchAddress("PFJet_muonEnergy", &PFJet_muonEnergy, &b_PFJet_muonEnergy);
   fChain->SetBranchAddress("PFJet_muonMultiplicity", &PFJet_muonMultiplicity, &b_PFJet_muonMultiplicity);
   fChain->SetBranchAddress("PFJet_neutralEmEnergy", &PFJet_neutralEmEnergy, &b_PFJet_neutralEmEnergy);
   fChain->SetBranchAddress("PFJet_neutralHadronEnergy", &PFJet_neutralHadronEnergy, &b_PFJet_neutralHadronEnergy);
   fChain->SetBranchAddress("PFJet_neutralHadronMultiplicity", &PFJet_neutralHadronMultiplicity, &b_PFJet_neutralHadronMultiplicity);
   fChain->SetBranchAddress("PFJet_photonEnergy", &PFJet_photonEnergy, &b_PFJet_photonEnergy);
   fChain->SetBranchAddress("PFJet_photonMultiplicity", &PFJet_photonMultiplicity, &b_PFJet_photonMultiplicity);
   fChain->SetBranchAddress("PFJet_jetArea", &PFJet_jetArea, &b_PFJet_jetArea);
   fChain->SetBranchAddress("PFJet_maxDistance", &PFJet_maxDistance, &b_PFJet_maxDistance);
   fChain->SetBranchAddress("PFJet_nConstituents", &PFJet_nConstituents, &b_PFJet_nConstituents);
   fChain->SetBranchAddress("PFJet_pileup", &PFJet_pileup, &b_PFJet_pileup);
   fChain->SetBranchAddress("PFJet_etaetaMoment", &PFJet_etaetaMoment, &b_PFJet_etaetaMoment);
   fChain->SetBranchAddress("PFJet_etaphiMoment", &PFJet_etaphiMoment, &b_PFJet_etaphiMoment);
   fChain->SetBranchAddress("PFJet_Track_idx", &PFJet_Track_idx, &b_PFJet_Track_idx);
   fChain->SetBranchAddress("PFJet_MatchedHPS_idx", &PFJet_MatchedHPS_idx, &b_PFJet_MatchedHPS_idx);
   fChain->SetBranchAddress("PFJet_numberOfDaughters", &PFJet_numberOfDaughters, &b_PFJet_numberOfDaughters);
   fChain->SetBranchAddress("PFJet_chargedEmEnergyFraction", &PFJet_chargedEmEnergyFraction, &b_PFJet_chargedEmEnergyFraction);
   fChain->SetBranchAddress("PFJet_chargedHadronEnergyFraction", &PFJet_chargedHadronEnergyFraction, &b_PFJet_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("PFJet_neutralHadronEnergyFraction", &PFJet_neutralHadronEnergyFraction, &b_PFJet_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("PFJet_neutralEmEnergyFraction", &PFJet_neutralEmEnergyFraction, &b_PFJet_neutralEmEnergyFraction);
   fChain->SetBranchAddress("PFJet_PUJetID_discr", &PFJet_PUJetID_discr, &b_PFJet_PUJetID_discr);
   fChain->SetBranchAddress("PFJet_PUJetID_looseWP", &PFJet_PUJetID_looseWP, &b_PFJet_PUJetID_looseWP);
   fChain->SetBranchAddress("PFJet_PUJetID_mediumWP", &PFJet_PUJetID_mediumWP, &b_PFJet_PUJetID_mediumWP);
   fChain->SetBranchAddress("PFJet_PUJetID_tightWP", &PFJet_PUJetID_tightWP, &b_PFJet_PUJetID_tightWP);
   fChain->SetBranchAddress("PFJet_partonFlavour", &PFJet_partonFlavour, &b_PFJet_partonFlavour);
   fChain->SetBranchAddress("PFJet_bDiscriminator", &PFJet_bDiscriminator, &b_PFJet_bDiscriminator);
   fChain->SetBranchAddress("PFJet_TracksP4", &PFJet_TracksP4, &b_PFJet_TracksP4);
   fChain->SetBranchAddress("PFJet_nTrk", &PFJet_nTrk, &b_PFJet_nTrk);
   fChain->SetBranchAddress("PFJet_JECuncertainty", &PFJet_JECuncertainty, &b_PFJet_JECuncertainty);
   fChain->SetBranchAddress("PFJet_GenJet_p4", &PFJet_GenJet_p4, &b_PFJet_GenJet_p4);
   fChain->SetBranchAddress("PFJet_GenJet_Constituents_p4", &PFJet_GenJet_Constituents_p4, &b_PFJet_GenJet_Constituents_p4);
   fChain->SetBranchAddress("PFJet_GenJetNoNu_p4", &PFJet_GenJetNoNu_p4, &b_PFJet_GenJetNoNu_p4);
   fChain->SetBranchAddress("PFJet_GenJetNoNu_Constituents_p4", &PFJet_GenJetNoNu_Constituents_p4, &b_PFJet_GenJetNoNu_Constituents_p4);
   fChain->SetBranchAddress("isPatMET", &isPatMET, &b_isPatMET);
   fChain->SetBranchAddress("MET_Uncorr_et", &MET_Uncorr_et, &b_MET_Uncorr_et);
   fChain->SetBranchAddress("MET_Uncorr_pt", &MET_Uncorr_pt, &b_MET_Uncorr_pt);
   fChain->SetBranchAddress("MET_Uncorr_phi", &MET_Uncorr_phi, &b_MET_Uncorr_phi);
   fChain->SetBranchAddress("MET_Uncorr_sumET", &MET_Uncorr_sumET, &b_MET_Uncorr_sumET);
   fChain->SetBranchAddress("MET_Uncorr_significance", &MET_Uncorr_significance, &b_MET_Uncorr_significance);
   fChain->SetBranchAddress("MET_Uncorr_significance_xx", &MET_Uncorr_significance_xx, &b_MET_Uncorr_significance_xx);
   fChain->SetBranchAddress("MET_Uncorr_significance_xy", &MET_Uncorr_significance_xy, &b_MET_Uncorr_significance_xy);
   fChain->SetBranchAddress("MET_Uncorr_significance_yy", &MET_Uncorr_significance_yy, &b_MET_Uncorr_significance_yy);
   fChain->SetBranchAddress("MET_Uncorr_MuonEtFraction", &MET_Uncorr_MuonEtFraction, &b_MET_Uncorr_MuonEtFraction);
   fChain->SetBranchAddress("MET_Uncorr_NeutralEMFraction", &MET_Uncorr_NeutralEMFraction, &b_MET_Uncorr_NeutralEMFraction);
   fChain->SetBranchAddress("MET_Uncorr_NeutralHadEtFraction", &MET_Uncorr_NeutralHadEtFraction, &b_MET_Uncorr_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_Uncorr_Type6EtFraction", &MET_Uncorr_Type6EtFraction, &b_MET_Uncorr_Type6EtFraction);
   fChain->SetBranchAddress("MET_Uncorr_Type7EtFraction", &MET_Uncorr_Type7EtFraction, &b_MET_Uncorr_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrT0rt_et", &MET_CorrT0rt_et, &b_MET_CorrT0rt_et);
   fChain->SetBranchAddress("MET_CorrT0rt_pt", &MET_CorrT0rt_pt, &b_MET_CorrT0rt_pt);
   fChain->SetBranchAddress("MET_CorrT0rt_phi", &MET_CorrT0rt_phi, &b_MET_CorrT0rt_phi);
   fChain->SetBranchAddress("MET_CorrT0rt_sumET", &MET_CorrT0rt_sumET, &b_MET_CorrT0rt_sumET);
   fChain->SetBranchAddress("MET_CorrT0rt_MuonEtFraction", &MET_CorrT0rt_MuonEtFraction, &b_MET_CorrT0rt_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrT0rt_NeutralEMFraction", &MET_CorrT0rt_NeutralEMFraction, &b_MET_CorrT0rt_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrT0rt_NeutralHadEtFraction", &MET_CorrT0rt_NeutralHadEtFraction, &b_MET_CorrT0rt_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrT0rt_Type6EtFraction", &MET_CorrT0rt_Type6EtFraction, &b_MET_CorrT0rt_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrT0rt_Type7EtFraction", &MET_CorrT0rt_Type7EtFraction, &b_MET_CorrT0rt_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtT1_et", &MET_CorrT0rtT1_et, &b_MET_CorrT0rtT1_et);
   fChain->SetBranchAddress("MET_CorrT0rtT1_pt", &MET_CorrT0rtT1_pt, &b_MET_CorrT0rtT1_pt);
   fChain->SetBranchAddress("MET_CorrT0rtT1_phi", &MET_CorrT0rtT1_phi, &b_MET_CorrT0rtT1_phi);
   fChain->SetBranchAddress("MET_CorrT0rtT1_sumET", &MET_CorrT0rtT1_sumET, &b_MET_CorrT0rtT1_sumET);
   fChain->SetBranchAddress("MET_CorrT0rtT1_MuonEtFraction", &MET_CorrT0rtT1_MuonEtFraction, &b_MET_CorrT0rtT1_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtT1_NeutralEMFraction", &MET_CorrT0rtT1_NeutralEMFraction, &b_MET_CorrT0rtT1_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrT0rtT1_NeutralHadEtFraction", &MET_CorrT0rtT1_NeutralHadEtFraction, &b_MET_CorrT0rtT1_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtT1_Type6EtFraction", &MET_CorrT0rtT1_Type6EtFraction, &b_MET_CorrT0rtT1_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtT1_Type7EtFraction", &MET_CorrT0rtT1_Type7EtFraction, &b_MET_CorrT0rtT1_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrT0pc_et", &MET_CorrT0pc_et, &b_MET_CorrT0pc_et);
   fChain->SetBranchAddress("MET_CorrT0pc_pt", &MET_CorrT0pc_pt, &b_MET_CorrT0pc_pt);
   fChain->SetBranchAddress("MET_CorrT0pc_phi", &MET_CorrT0pc_phi, &b_MET_CorrT0pc_phi);
   fChain->SetBranchAddress("MET_CorrT0pc_sumET", &MET_CorrT0pc_sumET, &b_MET_CorrT0pc_sumET);
   fChain->SetBranchAddress("MET_CorrT0pc_MuonEtFraction", &MET_CorrT0pc_MuonEtFraction, &b_MET_CorrT0pc_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrT0pc_NeutralEMFraction", &MET_CorrT0pc_NeutralEMFraction, &b_MET_CorrT0pc_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrT0pc_NeutralHadEtFraction", &MET_CorrT0pc_NeutralHadEtFraction, &b_MET_CorrT0pc_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrT0pc_Type6EtFraction", &MET_CorrT0pc_Type6EtFraction, &b_MET_CorrT0pc_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrT0pc_Type7EtFraction", &MET_CorrT0pc_Type7EtFraction, &b_MET_CorrT0pc_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcT1_et", &MET_CorrT0pcT1_et, &b_MET_CorrT0pcT1_et);
   fChain->SetBranchAddress("MET_CorrT0pcT1_pt", &MET_CorrT0pcT1_pt, &b_MET_CorrT0pcT1_pt);
   fChain->SetBranchAddress("MET_CorrT0pcT1_phi", &MET_CorrT0pcT1_phi, &b_MET_CorrT0pcT1_phi);
   fChain->SetBranchAddress("MET_CorrT0pcT1_sumET", &MET_CorrT0pcT1_sumET, &b_MET_CorrT0pcT1_sumET);
   fChain->SetBranchAddress("MET_CorrT0pcT1_MuonEtFraction", &MET_CorrT0pcT1_MuonEtFraction, &b_MET_CorrT0pcT1_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcT1_NeutralEMFraction", &MET_CorrT0pcT1_NeutralEMFraction, &b_MET_CorrT0pcT1_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrT0pcT1_NeutralHadEtFraction", &MET_CorrT0pcT1_NeutralHadEtFraction, &b_MET_CorrT0pcT1_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcT1_Type6EtFraction", &MET_CorrT0pcT1_Type6EtFraction, &b_MET_CorrT0pcT1_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcT1_Type7EtFraction", &MET_CorrT0pcT1_Type7EtFraction, &b_MET_CorrT0pcT1_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtTxy_et", &MET_CorrT0rtTxy_et, &b_MET_CorrT0rtTxy_et);
   fChain->SetBranchAddress("MET_CorrT0rtTxy_pt", &MET_CorrT0rtTxy_pt, &b_MET_CorrT0rtTxy_pt);
   fChain->SetBranchAddress("MET_CorrT0rtTxy_phi", &MET_CorrT0rtTxy_phi, &b_MET_CorrT0rtTxy_phi);
   fChain->SetBranchAddress("MET_CorrT0rtTxy_sumET", &MET_CorrT0rtTxy_sumET, &b_MET_CorrT0rtTxy_sumET);
   fChain->SetBranchAddress("MET_CorrT0rtTxy_MuonEtFraction", &MET_CorrT0rtTxy_MuonEtFraction, &b_MET_CorrT0rtTxy_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtTxy_NeutralEMFraction", &MET_CorrT0rtTxy_NeutralEMFraction, &b_MET_CorrT0rtTxy_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrT0rtTxy_NeutralHadEtFraction", &MET_CorrT0rtTxy_NeutralHadEtFraction, &b_MET_CorrT0rtTxy_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtTxy_Type6EtFraction", &MET_CorrT0rtTxy_Type6EtFraction, &b_MET_CorrT0rtTxy_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtTxy_Type7EtFraction", &MET_CorrT0rtTxy_Type7EtFraction, &b_MET_CorrT0rtTxy_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtT1Txy_et", &MET_CorrT0rtT1Txy_et, &b_MET_CorrT0rtT1Txy_et);
   fChain->SetBranchAddress("MET_CorrT0rtT1Txy_pt", &MET_CorrT0rtT1Txy_pt, &b_MET_CorrT0rtT1Txy_pt);
   fChain->SetBranchAddress("MET_CorrT0rtT1Txy_phi", &MET_CorrT0rtT1Txy_phi, &b_MET_CorrT0rtT1Txy_phi);
   fChain->SetBranchAddress("MET_CorrT0rtT1Txy_sumET", &MET_CorrT0rtT1Txy_sumET, &b_MET_CorrT0rtT1Txy_sumET);
   fChain->SetBranchAddress("MET_CorrT0rtT1Txy_MuonEtFraction", &MET_CorrT0rtT1Txy_MuonEtFraction, &b_MET_CorrT0rtT1Txy_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtT1Txy_NeutralEMFraction", &MET_CorrT0rtT1Txy_NeutralEMFraction, &b_MET_CorrT0rtT1Txy_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrT0rtT1Txy_NeutralHadEtFraction", &MET_CorrT0rtT1Txy_NeutralHadEtFraction, &b_MET_CorrT0rtT1Txy_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtT1Txy_Type6EtFraction", &MET_CorrT0rtT1Txy_Type6EtFraction, &b_MET_CorrT0rtT1Txy_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrT0rtT1Txy_Type7EtFraction", &MET_CorrT0rtT1Txy_Type7EtFraction, &b_MET_CorrT0rtT1Txy_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcTxy_et", &MET_CorrT0pcTxy_et, &b_MET_CorrT0pcTxy_et);
   fChain->SetBranchAddress("MET_CorrT0pcTxy_pt", &MET_CorrT0pcTxy_pt, &b_MET_CorrT0pcTxy_pt);
   fChain->SetBranchAddress("MET_CorrT0pcTxy_phi", &MET_CorrT0pcTxy_phi, &b_MET_CorrT0pcTxy_phi);
   fChain->SetBranchAddress("MET_CorrT0pcTxy_sumET", &MET_CorrT0pcTxy_sumET, &b_MET_CorrT0pcTxy_sumET);
   fChain->SetBranchAddress("MET_CorrT0pcTxy_MuonEtFraction", &MET_CorrT0pcTxy_MuonEtFraction, &b_MET_CorrT0pcTxy_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcTxy_NeutralEMFraction", &MET_CorrT0pcTxy_NeutralEMFraction, &b_MET_CorrT0pcTxy_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrT0pcTxy_NeutralHadEtFraction", &MET_CorrT0pcTxy_NeutralHadEtFraction, &b_MET_CorrT0pcTxy_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcTxy_Type6EtFraction", &MET_CorrT0pcTxy_Type6EtFraction, &b_MET_CorrT0pcTxy_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcTxy_Type7EtFraction", &MET_CorrT0pcTxy_Type7EtFraction, &b_MET_CorrT0pcTxy_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcT1Txy_et", &MET_CorrT0pcT1Txy_et, &b_MET_CorrT0pcT1Txy_et);
   fChain->SetBranchAddress("MET_CorrT0pcT1Txy_pt", &MET_CorrT0pcT1Txy_pt, &b_MET_CorrT0pcT1Txy_pt);
   fChain->SetBranchAddress("MET_CorrT0pcT1Txy_phi", &MET_CorrT0pcT1Txy_phi, &b_MET_CorrT0pcT1Txy_phi);
   fChain->SetBranchAddress("MET_CorrT0pcT1Txy_sumET", &MET_CorrT0pcT1Txy_sumET, &b_MET_CorrT0pcT1Txy_sumET);
   fChain->SetBranchAddress("MET_CorrT0pcT1Txy_MuonEtFraction", &MET_CorrT0pcT1Txy_MuonEtFraction, &b_MET_CorrT0pcT1Txy_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcT1Txy_NeutralEMFraction", &MET_CorrT0pcT1Txy_NeutralEMFraction, &b_MET_CorrT0pcT1Txy_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrT0pcT1Txy_NeutralHadEtFraction", &MET_CorrT0pcT1Txy_NeutralHadEtFraction, &b_MET_CorrT0pcT1Txy_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcT1Txy_Type6EtFraction", &MET_CorrT0pcT1Txy_Type6EtFraction, &b_MET_CorrT0pcT1Txy_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrT0pcT1Txy_Type7EtFraction", &MET_CorrT0pcT1Txy_Type7EtFraction, &b_MET_CorrT0pcT1Txy_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrT1_et", &MET_CorrT1_et, &b_MET_CorrT1_et);
   fChain->SetBranchAddress("MET_CorrT1_pt", &MET_CorrT1_pt, &b_MET_CorrT1_pt);
   fChain->SetBranchAddress("MET_CorrT1_phi", &MET_CorrT1_phi, &b_MET_CorrT1_phi);
   fChain->SetBranchAddress("MET_CorrT1_sumET", &MET_CorrT1_sumET, &b_MET_CorrT1_sumET);
   fChain->SetBranchAddress("MET_CorrT1_MuonEtFraction", &MET_CorrT1_MuonEtFraction, &b_MET_CorrT1_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrT1_NeutralEMFraction", &MET_CorrT1_NeutralEMFraction, &b_MET_CorrT1_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrT1_NeutralHadEtFraction", &MET_CorrT1_NeutralHadEtFraction, &b_MET_CorrT1_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrT1_Type6EtFraction", &MET_CorrT1_Type6EtFraction, &b_MET_CorrT1_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrT1_Type7EtFraction", &MET_CorrT1_Type7EtFraction, &b_MET_CorrT1_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrT1Txy_et", &MET_CorrT1Txy_et, &b_MET_CorrT1Txy_et);
   fChain->SetBranchAddress("MET_CorrT1Txy_pt", &MET_CorrT1Txy_pt, &b_MET_CorrT1Txy_pt);
   fChain->SetBranchAddress("MET_CorrT1Txy_phi", &MET_CorrT1Txy_phi, &b_MET_CorrT1Txy_phi);
   fChain->SetBranchAddress("MET_CorrT1Txy_sumET", &MET_CorrT1Txy_sumET, &b_MET_CorrT1Txy_sumET);
   fChain->SetBranchAddress("MET_CorrT1Txy_MuonEtFraction", &MET_CorrT1Txy_MuonEtFraction, &b_MET_CorrT1Txy_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrT1Txy_NeutralEMFraction", &MET_CorrT1Txy_NeutralEMFraction, &b_MET_CorrT1Txy_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrT1Txy_NeutralHadEtFraction", &MET_CorrT1Txy_NeutralHadEtFraction, &b_MET_CorrT1Txy_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrT1Txy_Type6EtFraction", &MET_CorrT1Txy_Type6EtFraction, &b_MET_CorrT1Txy_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrT1Txy_Type7EtFraction", &MET_CorrT1Txy_Type7EtFraction, &b_MET_CorrT1Txy_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrCaloT1_et", &MET_CorrCaloT1_et, &b_MET_CorrCaloT1_et);
   fChain->SetBranchAddress("MET_CorrCaloT1_pt", &MET_CorrCaloT1_pt, &b_MET_CorrCaloT1_pt);
   fChain->SetBranchAddress("MET_CorrCaloT1_phi", &MET_CorrCaloT1_phi, &b_MET_CorrCaloT1_phi);
   fChain->SetBranchAddress("MET_CorrCaloT1_sumET", &MET_CorrCaloT1_sumET, &b_MET_CorrCaloT1_sumET);
   fChain->SetBranchAddress("MET_CorrCaloT1T2_et", &MET_CorrCaloT1T2_et, &b_MET_CorrCaloT1T2_et);
   fChain->SetBranchAddress("MET_CorrCaloT1T2_pt", &MET_CorrCaloT1T2_pt, &b_MET_CorrCaloT1T2_pt);
   fChain->SetBranchAddress("MET_CorrCaloT1T2_phi", &MET_CorrCaloT1T2_phi, &b_MET_CorrCaloT1T2_phi);
   fChain->SetBranchAddress("MET_CorrCaloT1T2_sumET", &MET_CorrCaloT1T2_sumET, &b_MET_CorrCaloT1T2_sumET);
   fChain->SetBranchAddress("MET_CorrMVA_et", &MET_CorrMVA_et, &b_MET_CorrMVA_et);
   fChain->SetBranchAddress("MET_CorrMVA_pt", &MET_CorrMVA_pt, &b_MET_CorrMVA_pt);
   fChain->SetBranchAddress("MET_CorrMVA_phi", &MET_CorrMVA_phi, &b_MET_CorrMVA_phi);
   fChain->SetBranchAddress("MET_CorrMVA_sumET", &MET_CorrMVA_sumET, &b_MET_CorrMVA_sumET);
   fChain->SetBranchAddress("MET_CorrMVA_significance", &MET_CorrMVA_significance, &b_MET_CorrMVA_significance);
   fChain->SetBranchAddress("MET_CorrMVA_significance_xx", &MET_CorrMVA_significance_xx, &b_MET_CorrMVA_significance_xx);
   fChain->SetBranchAddress("MET_CorrMVA_significance_xy", &MET_CorrMVA_significance_xy, &b_MET_CorrMVA_significance_xy);
   fChain->SetBranchAddress("MET_CorrMVA_significance_yy", &MET_CorrMVA_significance_yy, &b_MET_CorrMVA_significance_yy);
   fChain->SetBranchAddress("MET_CorrMVA_MuonEtFraction", &MET_CorrMVA_MuonEtFraction, &b_MET_CorrMVA_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrMVA_NeutralEMFraction", &MET_CorrMVA_NeutralEMFraction, &b_MET_CorrMVA_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrMVA_NeutralHadEtFraction", &MET_CorrMVA_NeutralHadEtFraction, &b_MET_CorrMVA_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrMVA_Type6EtFraction", &MET_CorrMVA_Type6EtFraction, &b_MET_CorrMVA_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrMVA_Type7EtFraction", &MET_CorrMVA_Type7EtFraction, &b_MET_CorrMVA_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrMVA_srcMuon_p4", &MET_CorrMVA_srcMuon_p4, &b_MET_CorrMVA_srcMuon_p4);
   fChain->SetBranchAddress("MET_CorrMVA_srcElectron_p4", &MET_CorrMVA_srcElectron_p4, &b_MET_CorrMVA_srcElectron_p4);
   fChain->SetBranchAddress("MET_CorrMVA_srcTau_p4", &MET_CorrMVA_srcTau_p4, &b_MET_CorrMVA_srcTau_p4);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_et", &MET_CorrMVAMuTau_et, &b_MET_CorrMVAMuTau_et);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_pt", &MET_CorrMVAMuTau_pt, &b_MET_CorrMVAMuTau_pt);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_phi", &MET_CorrMVAMuTau_phi, &b_MET_CorrMVAMuTau_phi);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_sumET", &MET_CorrMVAMuTau_sumET, &b_MET_CorrMVAMuTau_sumET);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_significance", &MET_CorrMVAMuTau_significance, &b_MET_CorrMVAMuTau_significance);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_significance_xx", &MET_CorrMVAMuTau_significance_xx, &b_MET_CorrMVAMuTau_significance_xx);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_significance_xy", &MET_CorrMVAMuTau_significance_xy, &b_MET_CorrMVAMuTau_significance_xy);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_significance_yy", &MET_CorrMVAMuTau_significance_yy, &b_MET_CorrMVAMuTau_significance_yy);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_MuonEtFraction", &MET_CorrMVAMuTau_MuonEtFraction, &b_MET_CorrMVAMuTau_MuonEtFraction);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_NeutralEMFraction", &MET_CorrMVAMuTau_NeutralEMFraction, &b_MET_CorrMVAMuTau_NeutralEMFraction);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_NeutralHadEtFraction", &MET_CorrMVAMuTau_NeutralHadEtFraction, &b_MET_CorrMVAMuTau_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_Type6EtFraction", &MET_CorrMVAMuTau_Type6EtFraction, &b_MET_CorrMVAMuTau_Type6EtFraction);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_Type7EtFraction", &MET_CorrMVAMuTau_Type7EtFraction, &b_MET_CorrMVAMuTau_Type7EtFraction);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_srcMuon_p4", &MET_CorrMVAMuTau_srcMuon_p4, &b_MET_CorrMVAMuTau_srcMuon_p4);
   fChain->SetBranchAddress("MET_CorrMVAMuTau_srcTau_p4", &MET_CorrMVAMuTau_srcTau_p4, &b_MET_CorrMVAMuTau_srcTau_p4);
   fChain->SetBranchAddress("MET_Type1Corr_et", &MET_Type1Corr_et, &b_MET_Type1Corr_et);
   fChain->SetBranchAddress("MET_Type1Corr_pt", &MET_Type1Corr_pt, &b_MET_Type1Corr_pt);
   fChain->SetBranchAddress("MET_Type1Corr_phi", &MET_Type1Corr_phi, &b_MET_Type1Corr_phi);
   fChain->SetBranchAddress("MET_Type1Corr_sumET", &MET_Type1Corr_sumET, &b_MET_Type1Corr_sumET);
   fChain->SetBranchAddress("MET_Type1Corr_MuonEtFraction", &MET_Type1Corr_MuonEtFraction, &b_MET_Type1Corr_MuonEtFraction);
   fChain->SetBranchAddress("MET_Type1Corr_NeutralEMFraction", &MET_Type1Corr_NeutralEMFraction, &b_MET_Type1Corr_NeutralEMFraction);
   fChain->SetBranchAddress("MET_Type1Corr_NeutralHadEtFraction", &MET_Type1Corr_NeutralHadEtFraction, &b_MET_Type1Corr_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_Type1Corr_Type6EtFraction", &MET_Type1Corr_Type6EtFraction, &b_MET_Type1Corr_Type6EtFraction);
   fChain->SetBranchAddress("MET_Type1Corr_Type7EtFraction", &MET_Type1Corr_Type7EtFraction, &b_MET_Type1Corr_Type7EtFraction);
   fChain->SetBranchAddress("MET_Type1p2Corr_et", &MET_Type1p2Corr_et, &b_MET_Type1p2Corr_et);
   fChain->SetBranchAddress("MET_Type1p2Corr_pt", &MET_Type1p2Corr_pt, &b_MET_Type1p2Corr_pt);
   fChain->SetBranchAddress("MET_Type1p2Corr_phi", &MET_Type1p2Corr_phi, &b_MET_Type1p2Corr_phi);
   fChain->SetBranchAddress("MET_Type1p2Corr_sumET", &MET_Type1p2Corr_sumET, &b_MET_Type1p2Corr_sumET);
   fChain->SetBranchAddress("MET_Type1p2Corr_MuonEtFraction", &MET_Type1p2Corr_MuonEtFraction, &b_MET_Type1p2Corr_MuonEtFraction);
   fChain->SetBranchAddress("MET_Type1p2Corr_NeutralEMFraction", &MET_Type1p2Corr_NeutralEMFraction, &b_MET_Type1p2Corr_NeutralEMFraction);
   fChain->SetBranchAddress("MET_Type1p2Corr_NeutralHadEtFraction", &MET_Type1p2Corr_NeutralHadEtFraction, &b_MET_Type1p2Corr_NeutralHadEtFraction);
   fChain->SetBranchAddress("MET_Type1p2Corr_Type6EtFraction", &MET_Type1p2Corr_Type6EtFraction, &b_MET_Type1p2Corr_Type6EtFraction);
   fChain->SetBranchAddress("MET_Type1p2Corr_Type7EtFraction", &MET_Type1p2Corr_Type7EtFraction, &b_MET_Type1p2Corr_Type7EtFraction);
   fChain->SetBranchAddress("MET_Type1CorrElectronUp_et", &MET_Type1CorrElectronUp_et, &b_MET_Type1CorrElectronUp_et);
   fChain->SetBranchAddress("MET_Type1CorrElectronDown_et", &MET_Type1CorrElectronDown_et, &b_MET_Type1CorrElectronDown_et);
   fChain->SetBranchAddress("MET_Type1CorrMuonUp_et", &MET_Type1CorrMuonUp_et, &b_MET_Type1CorrMuonUp_et);
   fChain->SetBranchAddress("MET_Type1CorrMuonDown_et", &MET_Type1CorrMuonDown_et, &b_MET_Type1CorrMuonDown_et);
   fChain->SetBranchAddress("MET_Type1CorrTauUp_et", &MET_Type1CorrTauUp_et, &b_MET_Type1CorrTauUp_et);
   fChain->SetBranchAddress("MET_Type1CorrTauDown_et", &MET_Type1CorrTauDown_et, &b_MET_Type1CorrTauDown_et);
   fChain->SetBranchAddress("MET_Type1CorrJetResUp_et", &MET_Type1CorrJetResUp_et, &b_MET_Type1CorrJetResUp_et);
   fChain->SetBranchAddress("MET_Type1CorrJetResDown_et", &MET_Type1CorrJetResDown_et, &b_MET_Type1CorrJetResDown_et);
   fChain->SetBranchAddress("MET_Type1CorrJetEnUp_et", &MET_Type1CorrJetEnUp_et, &b_MET_Type1CorrJetEnUp_et);
   fChain->SetBranchAddress("MET_Type1CorrJetEnDown_et", &MET_Type1CorrJetEnDown_et, &b_MET_Type1CorrJetEnDown_et);
   fChain->SetBranchAddress("MET_Type1CorrUnclusteredUp_et", &MET_Type1CorrUnclusteredUp_et, &b_MET_Type1CorrUnclusteredUp_et);
   fChain->SetBranchAddress("MET_Type1CorrUnclusteredDown_et", &MET_Type1CorrUnclusteredDown_et, &b_MET_Type1CorrUnclusteredDown_et);
   fChain->SetBranchAddress("MET_Type1p2CorrElectronUp_et", &MET_Type1p2CorrElectronUp_et, &b_MET_Type1p2CorrElectronUp_et);
   fChain->SetBranchAddress("MET_Type1p2CorrElectronDown_et", &MET_Type1p2CorrElectronDown_et, &b_MET_Type1p2CorrElectronDown_et);
   fChain->SetBranchAddress("MET_Type1p2CorrMuonUp_et", &MET_Type1p2CorrMuonUp_et, &b_MET_Type1p2CorrMuonUp_et);
   fChain->SetBranchAddress("MET_Type1p2CorrMuonDown_et", &MET_Type1p2CorrMuonDown_et, &b_MET_Type1p2CorrMuonDown_et);
   fChain->SetBranchAddress("MET_Type1p2CorrTauUp_et", &MET_Type1p2CorrTauUp_et, &b_MET_Type1p2CorrTauUp_et);
   fChain->SetBranchAddress("MET_Type1p2CorrTauDown_et", &MET_Type1p2CorrTauDown_et, &b_MET_Type1p2CorrTauDown_et);
   fChain->SetBranchAddress("MET_Type1p2CorrJetResUp_et", &MET_Type1p2CorrJetResUp_et, &b_MET_Type1p2CorrJetResUp_et);
   fChain->SetBranchAddress("MET_Type1p2CorrJetResDown_et", &MET_Type1p2CorrJetResDown_et, &b_MET_Type1p2CorrJetResDown_et);
   fChain->SetBranchAddress("MET_Type1p2CorrJetEnUp_et", &MET_Type1p2CorrJetEnUp_et, &b_MET_Type1p2CorrJetEnUp_et);
   fChain->SetBranchAddress("MET_Type1p2CorrJetEnDown_et", &MET_Type1p2CorrJetEnDown_et, &b_MET_Type1p2CorrJetEnDown_et);
   fChain->SetBranchAddress("MET_Type1p2CorrUnclusteredUp_et", &MET_Type1p2CorrUnclusteredUp_et, &b_MET_Type1p2CorrUnclusteredUp_et);
   fChain->SetBranchAddress("MET_Type1p2CorrUnclusteredDown_et", &MET_Type1p2CorrUnclusteredDown_et, &b_MET_Type1p2CorrUnclusteredDown_et);
   fChain->SetBranchAddress("Event_EventNumber", &Event_EventNumber, &b_Event_EventNumber);
   fChain->SetBranchAddress("Event_RunNumber", &Event_RunNumber, &b_Event_RunNumber);
   fChain->SetBranchAddress("Event_bunchCrossing", &Event_bunchCrossing, &b_Event_bunchCrossing);
   fChain->SetBranchAddress("Event_orbitNumber", &Event_orbitNumber, &b_Event_orbitNumber);
   fChain->SetBranchAddress("Event_luminosityBlock", &Event_luminosityBlock, &b_Event_luminosityBlock);
   fChain->SetBranchAddress("Event_isRealData", &Event_isRealData, &b_Event_isRealData);
   fChain->SetBranchAddress("PileupInfo_TrueNumInteractions_nm1", &PileupInfo_TrueNumInteractions_nm1, &b_PileupInfo_TrueNumInteractions_nm1);
   fChain->SetBranchAddress("PileupInfo_TrueNumInteractions_n0", &PileupInfo_TrueNumInteractions_n0, &b_PileupInfo_TrueNumInteractions_n0);
   fChain->SetBranchAddress("PileupInfo_TrueNumInteractions_np1", &PileupInfo_TrueNumInteractions_np1, &b_PileupInfo_TrueNumInteractions_np1);
   fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
   fChain->SetBranchAddress("PUWeight_p5", &PUWeight_p5, &b_PUWeight_p5);
   fChain->SetBranchAddress("PUWeight_m5", &PUWeight_m5, &b_PUWeight_m5);
   fChain->SetBranchAddress("PUWeight3D", &PUWeight3D, &b_PUWeight3D);
   fChain->SetBranchAddress("PUWeight3D_p5", &PUWeight3D_p5, &b_PUWeight3D_p5);
   fChain->SetBranchAddress("PUWeight3D_m5", &PUWeight3D_m5, &b_PUWeight3D_m5);
   fChain->SetBranchAddress("PUWeightFineBins", &PUWeightFineBins, &b_PUWeightFineBins);
   fChain->SetBranchAddress("TauSpinnerWeight", &TauSpinnerWeight, &b_TauSpinnerWeight);
   fChain->SetBranchAddress("SelEffWeight", &SelEffWeight, &b_SelEffWeight);
   fChain->SetBranchAddress("MinVisPtFilter", &MinVisPtFilter, &b_MinVisPtFilter);
   fChain->SetBranchAddress("KinWeightPt", &KinWeightPt, &b_KinWeightPt);
   fChain->SetBranchAddress("KinWeightEta", &KinWeightEta, &b_KinWeightEta);
   fChain->SetBranchAddress("KinWeightMassPt", &KinWeightMassPt, &b_KinWeightMassPt);
   fChain->SetBranchAddress("EmbeddedWeight", &EmbeddedWeight, &b_EmbeddedWeight);
   fChain->SetBranchAddress("Track_p4", &Track_p4, &b_Track_p4);
   fChain->SetBranchAddress("Track_Poca", &Track_Poca, &b_Track_Poca);
   fChain->SetBranchAddress("Track_chi2", &Track_chi2, &b_Track_chi2);
   fChain->SetBranchAddress("Track_ndof", &Track_ndof, &b_Track_ndof);
   fChain->SetBranchAddress("Track_numberOfLostHits", &Track_numberOfLostHits, &b_Track_numberOfLostHits);
   fChain->SetBranchAddress("Track_numberOfValidHits", &Track_numberOfValidHits, &b_Track_numberOfValidHits);
   fChain->SetBranchAddress("Track_qualityMask", &Track_qualityMask, &b_Track_qualityMask);
   fChain->SetBranchAddress("Track_charge", &Track_charge, &b_Track_charge);
   fChain->SetBranchAddress("Track_pdgid", &Track_pdgid, &b_Track_pdgid);
   fChain->SetBranchAddress("Track_B", &Track_B, &b_Track_B);
   fChain->SetBranchAddress("Track_M", &Track_M, &b_Track_M);
   fChain->SetBranchAddress("Track_par", &Track_par, &b_Track_par);
   fChain->SetBranchAddress("Track_cov", &Track_cov, &b_Track_cov);
   fChain->SetBranchAddress("GenEventInfoProduct_signalProcessID", &GenEventInfoProduct_signalProcessID, &b_GenEventInfoProduct_signalProcessID);
   fChain->SetBranchAddress("GenEventInfoProduct_weight", &GenEventInfoProduct_weight, &b_GenEventInfoProduct_weight);
   fChain->SetBranchAddress("GenEventInfoProduct_weights", &GenEventInfoProduct_weights, &b_GenEventInfoProduct_weights);
   fChain->SetBranchAddress("GenEventInfoProduct_qScale", &GenEventInfoProduct_qScale, &b_GenEventInfoProduct_qScale);
   fChain->SetBranchAddress("GenEventInfoProduct_alphaQED", &GenEventInfoProduct_alphaQED, &b_GenEventInfoProduct_alphaQED);
   fChain->SetBranchAddress("GenEventInfoProduct_alphaQCD", &GenEventInfoProduct_alphaQCD, &b_GenEventInfoProduct_alphaQCD);
   fChain->SetBranchAddress("GenEventInfoProduct_id1", &GenEventInfoProduct_id1, &b_GenEventInfoProduct_id1);
   fChain->SetBranchAddress("GenEventInfoProduct_id2", &GenEventInfoProduct_id2, &b_GenEventInfoProduct_id2);
   fChain->SetBranchAddress("GenEventInfoProduct_x1", &GenEventInfoProduct_x1, &b_GenEventInfoProduct_x1);
   fChain->SetBranchAddress("GenEventInfoProduct_x2", &GenEventInfoProduct_x2, &b_GenEventInfoProduct_x2);
   fChain->SetBranchAddress("GenEventInfoProduct_scalePDF", &GenEventInfoProduct_scalePDF, &b_GenEventInfoProduct_scalePDF);
   fChain->SetBranchAddress("MC_p4", &MC_p4, &b_MC_p4);
   fChain->SetBranchAddress("MC_pdgid", &MC_pdgid, &b_MC_pdgid);
   fChain->SetBranchAddress("MC_charge", &MC_charge, &b_MC_charge);
   fChain->SetBranchAddress("MC_midx", &MC_midx, &b_MC_midx);
   fChain->SetBranchAddress("MC_childpdgid", &MC_childpdgid, &b_MC_childpdgid);
   fChain->SetBranchAddress("MC_childidx", &MC_childidx, &b_MC_childidx);
   fChain->SetBranchAddress("MC_status", &MC_status, &b_MC_status);
   fChain->SetBranchAddress("MCSignalParticle_p4", &MCSignalParticle_p4, &b_MCSignalParticle_p4);
   fChain->SetBranchAddress("MCSignalParticle_pdgid", &MCSignalParticle_pdgid, &b_MCSignalParticle_pdgid);
   fChain->SetBranchAddress("MCSignalParticle_childpdgid", &MCSignalParticle_childpdgid, &b_MCSignalParticle_childpdgid);
   fChain->SetBranchAddress("MCSignalParticle_charge", &MCSignalParticle_charge, &b_MCSignalParticle_charge);
   fChain->SetBranchAddress("MCSignalParticle_Poca", &MCSignalParticle_Poca, &b_MCSignalParticle_Poca);
   fChain->SetBranchAddress("MCSignalParticle_Tauidx", &MCSignalParticle_Tauidx, &b_MCSignalParticle_Tauidx);
   fChain->SetBranchAddress("MCTauandProd_p4", &MCTauandProd_p4, &b_MCTauandProd_p4);
   fChain->SetBranchAddress("MCTauandProd_Vertex", &MCTauandProd_Vertex, &b_MCTauandProd_Vertex);
   fChain->SetBranchAddress("MCTauandProd_pdgid", &MCTauandProd_pdgid, &b_MCTauandProd_pdgid);
   fChain->SetBranchAddress("MCTauandProd_midx", &MCTauandProd_midx, &b_MCTauandProd_midx);
   fChain->SetBranchAddress("MCTauandProd_charge", &MCTauandProd_charge, &b_MCTauandProd_charge);
   fChain->SetBranchAddress("MCTau_JAK", &MCTau_JAK, &b_MCTau_JAK);
   fChain->SetBranchAddress("MCTau_DecayBitMask", &MCTau_DecayBitMask, &b_MCTau_DecayBitMask);
   fChain->SetBranchAddress("HTLTriggerName", &HTLTriggerName, &b_HTLTriggerName);
   fChain->SetBranchAddress("TriggerAccept", &TriggerAccept, &b_TriggerAccept);
   fChain->SetBranchAddress("TriggerError", &TriggerError, &b_TriggerError);
   fChain->SetBranchAddress("TriggerWasRun", &TriggerWasRun, &b_TriggerWasRun);
   fChain->SetBranchAddress("HLTPrescale", &HLTPrescale, &b_HLTPrescale);
   fChain->SetBranchAddress("NHLTL1GTSeeds", &NHLTL1GTSeeds, &b_NHLTL1GTSeeds);
   fChain->SetBranchAddress("L1SEEDPrescale", &L1SEEDPrescale, &b_L1SEEDPrescale);
   fChain->SetBranchAddress("L1SEEDInvalidPrescale", &L1SEEDInvalidPrescale, &b_L1SEEDInvalidPrescale);
   fChain->SetBranchAddress("L1SEEDisTechBit", &L1SEEDisTechBit, &b_L1SEEDisTechBit);
   fChain->SetBranchAddress("TriggerMatchMuon", &TriggerMatchMuon, &b_TriggerMatchMuon);
   fChain->SetBranchAddress("TriggerMatchJet", &TriggerMatchJet, &b_TriggerMatchJet);
   fChain->SetBranchAddress("TriggerMatchTau", &TriggerMatchTau, &b_TriggerMatchTau);
   fChain->SetBranchAddress("HLTTrigger_objs_Pt", &HLTTrigger_objs_Pt, &b_HLTTrigger_objs_Pt);
   fChain->SetBranchAddress("HLTTrigger_objs_Eta", &HLTTrigger_objs_Eta, &b_HLTTrigger_objs_Eta);
   fChain->SetBranchAddress("HLTTrigger_objs_Phi", &HLTTrigger_objs_Phi, &b_HLTTrigger_objs_Phi);
   fChain->SetBranchAddress("HLTTrigger_objs_E", &HLTTrigger_objs_E, &b_HLTTrigger_objs_E);
   fChain->SetBranchAddress("HLTTrigger_objs_Id", &HLTTrigger_objs_Id, &b_HLTTrigger_objs_Id);
   fChain->SetBranchAddress("HLTTrigger_objs_trigger", &HLTTrigger_objs_trigger, &b_HLTTrigger_objs_trigger);
   fChain->SetBranchAddress("L1TriggerName", &L1TriggerName, &b_L1TriggerName);
   fChain->SetBranchAddress("L1TriggerDecision", &L1TriggerDecision, &b_L1TriggerDecision);
   fChain->SetBranchAddress("L1ErrorCode", &L1ErrorCode, &b_L1ErrorCode);
   fChain->SetBranchAddress("L1Prescale", &L1Prescale, &b_L1Prescale);
   Notify();
}

Bool_t NtupleReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NtupleReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NtupleReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
double NtupleReader::PFTau_Mass(unsigned int i)
{
  double m(-999);
  if(PFTau_isHPSByDecayModeFinding->size() >0 && PFTau_isHPSByDecayModeFinding->at(i)){
      double E = PFTau_p4->at(i).at(0);
      double px = PFTau_p4->at(i).at(1);
      double py = PFTau_p4->at(i).at(2);
      double pz = PFTau_p4->at(i).at(3);

      double p2 = px*px + py*py + pz*pz;

      m = (E*E > p2) ? sqrt(E*E - p2) : -sqrt(p2 - E*E);
  }
  return m;
}

#endif // #ifdef NtupleReader_cxx
