#ifndef ZToTaumuTauh_h
#define ZToTaumuTauh_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "ReferenceScaleFactors.h"

class ZToTaumuTauh : public Selection {

 public:
  ZToTaumuTauh(TString Name_, TString id_);
  virtual ~ZToTaumuTauh();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,
	  PrimeVtx,
	  NMuId,
	  NMuKin,
	  NMuIso,
	  NTauId,
	  NTauKin,
	  NTauIso,
	  ChargeSum,
	  TauDecayMode,
	  TauFLSigma,
	  MT_MuMET,
	  NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

  ReferenceScaleFactors *RSF;

  // cut values
  double cMu_dxy, cMu_dz, cMu_relIso, cMu_pt, cMu_eta, cMu_dRHltMatch;
  double cTau_pt, cTau_eta, cMuTau_dR, cTau_IsoRaw, cTau_dRHltMatch;
  std::vector<TString> cTriggerNames;

  double OneProngNoPiWeight;
  //dummy values
  int TriggerOkDummy, selVertexDummy, selMuonDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;
  double SB_lowerLimit, SB_upperLimit;
  bool Scaleby_Counting;
  bool verbose, Use_Embedded;

  // Histograms

  std::vector<std::vector<TH1D>* > Extradist1d_OS, Extradist1d_SS;

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH1D> NMtTauMET;
  std::vector<TH1D> NMvis, Mvis_SignalOnly, Mvis_SignalOnly_genMu, Mvis_SignalOnly_genA1, Mvis_SignalOnly_genTaumu, Mvis_SignalOnly_genTauh;
  std::vector<TH1D> NSignal_SB_WJets;
  std::vector<TH1D> NSB_Data;
  std::vector<TH1D> NSB;
  std::vector<TH1D> Mu_pt, Mu_eta, Mu_phi, Tau_pt, Tau_eta, Tau_phi, MET_phi;
  std::vector<TH1D> TauFL, TauFLSigned, TauFLSigmaSigned, TauFLSigmaUnsigned;
  std::vector<TH1D> A1mass, A1mass10GeV;

  std::vector<TH1D> Mvis3Prong, Mvis1Prong, MvisIncl;
  std::vector<TH1D> MTMuMET3Prong, MTMuMET1Prong, MTMuMETIncl;

  std::vector<TH1D> dR_selTauh_genTauh, dR_selMu_genMu;
  std::vector<TH1D> POCAPV_Mag;
  std::vector<TH1D> Phi_SVPV, Phi_genTauh, Theta_SVPV, Theta_genTauh, dPhi_SVPV_genTauh, dTheta_SVPV_genTauh, Angle_SVPV_genTauh;
  std::vector<TH1D> Phi_POCAPV, Phi_genTaumu, Theta_POCAPV, Theta_genTaumu, dPhi_POCAPV_genTaumu, dTheta_POCAPV_genTaumu;
  std::vector<TH1D> dPhi_MinusSVPV_genTaumu, dTheta_MinusSVPV_genTaumu, Angle_MinusSVPV_genTaumu;
  std::vector<TH1D> GJ_Tauh, GJ_Taumu;
  std::vector<TH1D> dPhi_DiTauGen, Pt_DiTauGen, Pt_ZGen, M_ZGen, M_DiTauPtBalance, dM_DiTau, dPt_GenTaumuPtBalance, dP_GenTaumuPtBalance, dP_GenTauh;
  std::vector<TH2D> dP_GenTauMuPtBalance_vs_dPTauh, Pt_vs_dPhi_DiTauGen;
  std::vector<TH2D> TauFLSigmaCut_vs_Res, TauFLSigma_vs_Res;

  std::vector<TH1D> TPTF_TauA1_pRes_BestSolution, TPTF_TauA1_pRes_FitSolution;
  std::vector<TH2D> TPTF_TauA1_BestSolution_vs_FitSolution, TPTF_TauA1_pRes_vs_GenA1Mass_BestSolution, TPTF_TauA1_pRes_vs_GenA1Mass_FitSolution, TPTF_TauA1_pRes_vs_GenGJAngle_FitSolution, TPTF_TauA1_pRes_vs_GenGJAngle_BestSolution, TPTF_TauA1_pRes_vs_RecoChi2_FitSolution;
  std::vector<TH2D> TPTF_A1_pRes_vs_GenGJAngle;
  std::vector<TH2D> TPTF_TauA1_pRes_vs_RecoA1Mass_FitSolution, TPTF_TauA1_pRes_vs_RecoGJAngle_FitSolution, TPTF_A1_pRes_vs_RecoGJAngle;

  std::vector<TH1D> Reco_ZMass, Reco_ZMass_UnboostedGenZ, Reco_EventFit_Solution, Reco_A1Fit_Solution, Reco_Chi2, Reco_Chi2_FitSolutionOnly, Reco_Chi2_FitSolutionOnlyLargeScale, Reco_ConstrainedDeltaSum, Reco_ConstrainedDeltaMass, Reco_ConstrainedDeltaPt, Reco_NIter;
  std::vector<TH1D> GenReco_ZMass, GenReco_EventFit_Solution, GenReco_A1Fit_Solution, GenReco_Chi2, GenReco_Chi2_FitSolutionOnly, GenReco_ConstrainedDeltaSum, GenReco_NIter;
  std::vector<TH1D> Reco_PtRes_TauA1, Reco_PtRes_TauMu;
  std::vector<TH1D> Reco_ZMass_MassScan, Reco_ZMasswithProbWeight_MassScan, Reco_ProbStack_MassScan, Reco_ZMass_PDF;
  std::vector<TH1D> GenZ_Pt_Unboosted, RecoZ_Pt_Unboosted;

  std::vector<TH1D> NQCD;
  std::vector<TH1D> QCD_MT_MuMET_A, QCD_MT_MuMET_B, QCD_MT_MuMET_C, QCD_MT_MuMET_D;
  std::vector<TH1D> QCD_MET_A, QCD_MET_B, QCD_MET_C, QCD_MET_D;

  bool selectMuon_Id(unsigned i, unsigned vertex);
  bool selectMuon_Kinematics(unsigned i);
  bool selectMuon_Isolation(unsigned i);
  bool selectMuon_AntiIsolation(unsigned i);

  bool selectPFTau_Id(unsigned i);
  bool selectPFTau_Id(unsigned i, std::vector<int> muonCollection);
  bool selectPFTau_Isolation(unsigned i);
  bool selectPFTau_Kinematics(unsigned i);
  double Reconstruct_hadronicTauEnergy(unsigned i);

 private:

};
#endif
