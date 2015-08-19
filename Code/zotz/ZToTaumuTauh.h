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
	  DiMuonVeto,
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
  double cDiMuVeto_pt, cDiMuVeto_eta, cDiMuVeto_dz, cDiMuVeto_dBetaCombIso, cDiMuVeto_dR;
  double cTau_pt, cTau_eta, cMuTau_dR, cTau_IsoRaw, cTau_dRHltMatch;
  std::vector<TString> cTriggerNames;

  double OneProngNoPiWeight;
  //dummy values
  int TriggerOkDummy, selVertexDummy, selMuonDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;
  double SB_lowerLimit, SB_upperLimit;
  bool Scaleby_Counting;
  bool selection_verbose, Use_Embedded;

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
  std::vector<TH1D> Mu_pt, Mu_eta, Mu_phi, Tau_pt, Tau_eta, Tau_phi, Tau_Mass_Inclusive, Tau_Mass_sq_Inclusive, Tau_Mass_Inclusive_NoTLV, Tau_Mass_Inclusive_UnFitTracks, Tau_Mass_Inclusive_ReFitTracks, Tau_Mass_Difference_PFTau_UnFitTracks_3PS, Tau_Mass_Difference_PFTau_ReFitTracks_3PS;
  std::vector<TH1D> MET_phi;
  std::vector<TH1D> TauFL_NoTauFLSigmaCut, TauFLSigned_NoTauFLSigmaCut, TauFLSigmaSigned, TauFLSigmaUnsigned;
  std::vector<TH1D> A1mass, A1mass10GeV, A1massRefit, dA1mass_PFTau_Refit;
  std::vector<TH1D> A1_Phi_Res, A1_Theta_Res;

  std::vector<TH1D> Mvis3Prong, Mvis1Prong, MvisIncl;
  std::vector<TH1D> MTMuMET3Prong, MTMuMET1Prong, MTMuMETIncl;

  std::vector<TH1D> dR_selTauh_genTauh, dR_selMu_genMu;
  std::vector<TH1D> POCAPV_Mag;
  std::vector<TH1D> Phi_SVPV, Phi_genTauh, Theta_SVPV, Theta_genTauh, dPhi_SVPV_genTauh, dTheta_SVPV_genTauh, Angle_SVPV_genTauh;
  std::vector<TH2D> dPhi_SVPV_genTauh_vs_TauFL, dPhi_SVPV_genTauhPlus_vs_TauFL, dPhi_SVPV_genTauhMinus_vs_TauFL;
  std::vector<TH1D> Phi_POCAPV, Phi_genTaumu, Theta_POCAPV, Theta_genTaumu, dPhi_POCAPV_genTaumu, dTheta_POCAPV_genTaumu;
  std::vector<TH1D> dPhi_MinusSVPV_genTaumu, dTheta_MinusSVPV_genTaumu, Angle_MinusSVPV_genTaumu;
  std::vector<TH1D> dPhi_GenTauMu_GenMu, dPhi_GenTauMu_RecoMu, dTheta_GenTauMu_GenMu, dTheta_GenTauMu_RecoMu;

  std::vector<TH1D> GJAngle_Over_GJAngleMax_StraightTau, GJAngle_Over_GJAngleMax_HelixTau;
  std::vector<TH1D> dGJAngle_GJAngleMAX_StraightTau, dGJAngle_GJAngleMAX_HelixTau, Angle_HelixTau_StraightTau, dGJAngle_HelixTau_StraightTau, dGJAngle_HelixTau_StraightTauOverGJAngle, TauA1_Reco_Solution_StraightTau, TauA1_Reco_Solution_HelixTau;
  std::vector<TH1D> NUnphysical_StraightTau_HelixTau;

  std::vector<TH1D> Gen_TauA1_GJ, Gen_TauMu_GJ;
  std::vector<TH1D> Gen_TauA1_P, Gen_TauA1_Pt, Gen_TauA1_Px, Gen_TauA1_Py, Gen_TauA1_Pz, Gen_TauA1_Phi, Gen_TauA1_Eta;
  std::vector<TH1D> Gen_TauMu_P, Gen_TauMu_Pt, Gen_TauMu_Px, Gen_TauMu_Py, Gen_TauMu_Pz, Gen_TauMu_Phi, Gen_TauMu_Eta;
  std::vector<TH1D> Gen_TauA1_P_noSel, Gen_TauA1_Pt_noSel, Gen_TauA1_Px_noSel, Gen_TauA1_Py_noSel, Gen_TauA1_Pz_noSel, Gen_TauA1_Phi_noSel, Gen_TauA1_Eta_noSel;
  std::vector<TH1D> Gen_TauMu_P_noSel, Gen_TauMu_Pt_noSel, Gen_TauMu_Px_noSel, Gen_TauMu_Py_noSel, Gen_TauMu_Pz_noSel, Gen_TauMu_Phi_noSel, Gen_TauMu_Eta_noSel;
  std::vector<TH1D> Gen_Z_Pt, Gen_Z_M, Gen_Z_Phi, Gen_Z_Eta;
  std::vector<TH1D> Gen_Z_Pt_noSel, Gen_Z_Phi_noSel, Gen_Z_Eta_noSel;
  std::vector<TH1D> Gen_DiTau_dPhi, Gen_DiTau_Pt, Gen_DiTau_PtBalance_M, Gen_DiTau_PtBalance_dM, Gen_TauMu_PtBalance_Pt, Gen_TauMu_PtBalance_dP, Gen_TauA1_dP;
  std::vector<TH1D> Gen_TPTF_TauA1_Solution_NoSelection, Gen_TPTF_TauA1_Solution_WithSelection;
  std::vector<TH2D> Gen_Z_Pt_vs_MET;
  std::vector<TH2D> Gen_Z_Pt_vs_VtxTracksPt;
  std::vector<TH2D> Gen_Z_Phi_vs_VtxTracksPhi;
  std::vector<TH1D> VtxTracksPtRes, VtxTracksPhiCorrectedRes;
  std::vector<TH2D> dP_GenTauMuPtBalance_vs_dPTauh, Pt_vs_dPhi_DiTauGen;
  std::vector<TH1D> TauFL_WithTauFLSigmaCut;
  std::vector<TH2D> TauFLSigmaCut_vs_Res, TauFLSigma_vs_Res, TauFLSigma_vs_UnphysicalAll, TauFL_vs_UnphysicalAll;
  std::vector<TH1D> TauFLSigma_vs_UnphysicalProb, TauFL_vs_UnphysicalProb;

  std::vector<TH1D> TPTF_TauA1_pRes_BestSolution, TPTF_TauA1_pRes_FitSolution;
  std::vector<TH2D> TPTF_TauA1_BestSolution_vs_FitSolution, TPTF_TauA1_RightSolution_vs_FitSolution, TPTF_TauA1_pRes_vs_GenA1Mass_BestSolution, TPTF_TauA1_pRes_vs_GenA1Mass_FitSolution, TPTF_TauA1_pRes_vs_GenGJAngle_FitSolution, TPTF_TauA1_pRes_vs_GenGJAngle_BestSolution, TPTF_TauA1_pRes_vs_RecoChi2_FitSolution;
  std::vector<TH2D> TPTF_A1_pRes_vs_GenGJAngle;
  std::vector<TH2D> TPTF_TauA1_pRes_vs_RecoA1Mass_FitSolution, TPTF_TauA1_pRes_vs_RecoGJAngle_FitSolution, TPTF_TauA1_p_iRes_sq_vs_RecoGJAngle_FitSolution, TPTF_A1_pRes_vs_RecoGJAngle;
  std::vector<TH2D> TPTF_TauA1_pxRes_vs_RecoGJAngle_FitSolution, TPTF_TauA1_pyRes_vs_RecoGJAngle_FitSolution, TPTF_TauA1_pzRes_vs_RecoGJAngle_FitSolution, TPTF_TauA1_ptRes_vs_RecoGJAngle_FitSolution;
  std::vector<TH2D> TPTF_TauA1_pxsqRes_vs_RecoGJAngle_FitSolution, TPTF_TauA1_pysqRes_vs_RecoGJAngle_FitSolution, TPTF_TauA1_pzsqRes_vs_RecoGJAngle_FitSolution;
  std::vector<TH2D> TPTF_TauA1_p_orthoRes_vs_RecoGJAngle_FitSolution, TPTF_TauA1_p_paralRes_vs_RecoGJAngle_FitSolution;
  std::vector<TH1D> TPTF_TauA1_pxsq_Reco, TPTF_TauA1_pysq_Reco, TPTF_TauA1_pzsq_Reco, TPTF_TauA1_pxsq_Gen, TPTF_TauA1_pysq_Gen, TPTF_TauA1_pzsq_Gen;

  std::vector<TH2D> TPTF_TauA1_ptRes_vs_ptGen, TPTF_TauA1_ptRes_vs_ptReco;

  std::vector<TH1D> TPTF_Neutrino_UnFitTracks_Mass;
  std::vector<TH2D> TPTF_Neutrino_UnFitTracks_Mass_vs_TauFL;
  std::vector<TH1D> TPTF_Neutrino_ReFitTracks_Mass;
  std::vector<TH2D> TPTF_Neutrino_ReFitTracks_Mass_vs_TauFL;
  std::vector<TH1D> TPTF_Neutrino_PFTau_Mass;
  std::vector<TH2D> TPTF_Neutrino_PFTau_Mass_vs_TauFL;
  std::vector<TH1D> TPTF_Neutrino_RefitPFTau_HelixGenTau_Mass;
  std::vector<TH1D> TPTF_Neutrino_GenA1_StraightGenTau_Mass;
  std::vector<TH1D> TPTF_Neutrino_GenA1_HelixGenTau_Mass;

  std::vector<TH1D> TransTrk_Failure_withSelection, TransTrk_Failure_noSelection;

  std::vector<TH1D> Est_Z_Pt_wTruth, Est_Z_PtRes_wTruth;
  std::vector<TH2D> Est_Z_Pt_wTruth_vs_GenZ_Pt;
  std::vector<TH1D> Est_Z_Pt_alwaysMinus, Est_Z_PtRes_alwaysMinus;
  std::vector<TH2D> Est_Z_Pt_alwaysMinus_vs_GenZ_Pt;
  std::vector<TH1D> Est_Z_Energy_wTruth, Est_Z_Energy_alwaysMinus;
  std::vector<TH1D> Est_Z_EnergyRes_wTruth, Est_Z_EnergyRes_alwaysMinus;
  std::vector<TH1D> Est_Z_EnergyRes_wTruth2, Est_Z_EnergyRes_alwaysMinus2;
  std::vector<TH1D> Est_TauMu_PtRes_wTruth, Est_TauMu_PtRes_wTruth2;
  std::vector<TH1D> Est_TauMu_PtRes_wTruth_NoZMass, Est_Z_M_wTruth_NoZMass, Est_Z_EnergyRes_wTruth_NoZMass;
  std::vector<TH1D> Est_Z_M_wTruth2;
  std::vector<TH1D> Est_TauMu_wMET_PtRes, Est_TauMu_wMET_PhiRes, Est_TauMu_wMET_EtaRes, Est_Z_wMET_PtRes, Est_Z_wMET_PhiRes;

  std::vector<TH1D> Reco_TauA1_P, Reco_TauA1_Pt, Reco_TauA1_Px, Reco_TauA1_Py, Reco_TauA1_Pz, Reco_TauA1_Phi, Reco_TauA1_Eta;
  std::vector<TH1D> Reco_TauMu_P, Reco_TauMu_Pt, Reco_TauMu_Px, Reco_TauMu_Py, Reco_TauMu_Pz, Reco_TauMu_Phi, Reco_TauMu_Eta;
  std::vector<TH1D> Tau_pt_wo_FLSigmaCut, Tau_phi_wo_FLSigmaCut, Tau_eta_wo_FLSigmaCut;
  std::vector<TH1D> Reco_ZMass, Reco_ZMass_UnboostedGenZ, Reco_EventFit_Solution, Reco_A1Fit_Solution, Reco_Chi2_FitSolutionOnlyLargeScale, Reco_ConstrainedDeltaSum, Reco_ConstrainedDeltaMass, Reco_ConstrainedDeltaPt;
  std::vector<TH1D> Reco_Z_Energy_Res, RecoZ_Pt;
  std::vector<TH1D> GenReco_ZMass, GenReco_EventFit_Solution, GenReco_A1Fit_Solution, GenReco_Chi2, GenReco_Chi2_FitSolutionOnly, GenReco_ConstrainedDeltaSum, GenReco_NIter;
  std::vector<TH1D> Reco_PtRes_TauA1, Reco_PtRes_TauA1_AmbPoint0, Reco_PtRes_TauA1_AmbPoint12, Reco_PtRes_TauA1_AmbPoint1, Reco_PtRes_TauMu, Reco_PtRes_TauMu_AmbPoint0, Reco_PtRes_TauMu_AmbPoint12, Reco_PtRes_TauMu_AmbPoint1;
  std::vector<TH1D> Reco_PtRes_TauA1_LowZPt, Reco_PtRes_TauMu_LowZPt;
  std::vector<TH1D> Reco_PtRes_TauA1_NoFit, Reco_PtRes_TauA1_AmbPoint0_NoFit, Reco_PtRes_TauA1_AmbPoint12_NoFit, Reco_PtRes_TauMu_NoFit, Reco_PtRes_TauMu_AmbPoint0_NoFit, Reco_PtRes_TauMu_AmbPoint12_NoFit;
  std::vector<TH1D> Reco_PhiRes_TauMu_PreFit, Reco_ThetaRes_TauMu_PreFit;
  std::vector<TH1D> Reco_TauMu_DeltaPX_FitImpact, Reco_TauMu_DeltaPY_FitImpact, Reco_TauMu_DeltaPZ_FitImpact;
  std::vector<TH1D> Reco_TauA1_DeltaPX_FitImpact, Reco_TauA1_DeltaPY_FitImpact, Reco_TauA1_DeltaPZ_FitImpact;
  std::vector<TH1D> Reco_TauMu_ResCosTheta, Reco_TauMu_DeltaPhi_FitImpact;
  std::vector<TH1D> Reco_dPhi_TauMuTauA1_AfterFit, Reco_dPhi_TauMuTauA1_BeforeFit;
  std::vector<TH1D> Reco_dPhi_TauMuTauA1_AfterFit_lowBoost, Reco_dPhi_TauMuTauA1_BeforeFit_lowBoost;
  std::vector<TH1D> Reco_ZMass_MassScan, Reco_ZMass_MassScanUnboosted, Reco_ZMasswithProbWeight_MassScan, Reco_ProbStack_MassScan, Reco_ZMass_PDF;
  std::vector<TH1D> GenZ_Pt_Unboosted, RecoZ_Pt_Unboosted;
  std::vector<TH1D> Reco_Z_PhiRes, Reco_Z_PhiRes_noAmb, Reco_Z_PhiRes_wAmb, Reco_Z_PhiRes_pickedrightAmb, Reco_Z_PhiRes_pickedwrongAmb;
  std::vector<TH1D> Reco_Z_EtaRes, Reco_Z_EtaRes_noAmb, Reco_Z_EtaRes_wAmb, Reco_Z_EtaRes_pickedrightAmb, Reco_Z_EtaRes_pickedwrongAmb;
  std::vector<TH1D> Reco_Z_PRes, Reco_Z_PRes_noAmb, Reco_Z_PRes_wAmb, Reco_Z_PRes_pickedrightAmb, Reco_Z_PRes_pickedwrongAmb;
  std::vector<TH1D> Reco_NIter, Reco_NIter_noAmb, Reco_NIter_wAmb, Reco_NIter_pickedrightAmb, Reco_NIter_pickedwrongAmb;
  std::vector<TH1D> Reco_Chi2, Reco_Chi2_noAmb, Reco_Chi2_wAmb, Reco_Chi2_rightAmb, Reco_Chi2_wrongAmb, Reco_Chi2_pickedrightAmb, Reco_Chi2_pickedwrongAmb;
  std::vector<TH1D> Reco_Chi2_orig, Reco_Chi2_orig_noAmb, Reco_Chi2_orig_wAmb, Reco_Chi2_orig_rightAmb, Reco_Chi2_orig_wrongAmb, Reco_Chi2_orig_pickedrightAmb, Reco_Chi2_orig_pickedwrongAmb;
  std::vector<TH1D> Reco_Chi2_SC, Reco_Chi2_SC_noAmb, Reco_Chi2_SC_wAmb, Reco_Chi2_SC_rightAmb, Reco_Chi2_SC_wrongAmb, Reco_Chi2_SC_pickedrightAmb, Reco_Chi2_SC_pickedwrongAmb;
  std::vector<TH1D> Reco_Chi2_HC, Reco_Chi2_HC_noAmb, Reco_Chi2_HC_wAmb, Reco_Chi2_HC_rightAmb, Reco_Chi2_HC_wrongAmb, Reco_Chi2_HC_pickedrightAmb, Reco_Chi2_HC_pickedwrongAmb;

  std::vector<TH1D> Reco_Z_Px, Reco_Z_Py, Reco_Z_Pz, Reco_Z_Pt, Reco_Z_PtRes, Reco_Z_Phi, Reco_Z_Eta;

  std::vector<TH1D> Reco_Chi2_diff, Reco_Chi2_orig_diff, Reco_Chi2_HC_diff, Reco_Chi2_SC_diff, Reco_Chi2_origplusSC_diff;
  std::vector<TH2D> Reco_Chi2_diff_vs_correctAssignment;

  std::vector<TH2D> Reco_FitSolution_byChi2_Full_vs_RightSolution;
  std::vector<TH2D> Reco_FitSolution_byChi2_orig_vs_RightSolution, Reco_FitSolution_byChi2_SC_vs_RightSolution, Reco_FitSolution_byChi2_HC_vs_RightSolution;
  std::vector<TH2D> Reco_FitSolution_byChi2_origplusSC_vs_RightSolution;

  std::vector<TH1D> Reco_TauA1_P_wRecoil, Reco_TauA1_Pt_wRecoil, Reco_TauA1_Px_wRecoil, Reco_TauA1_Py_wRecoil, Reco_TauA1_Pz_wRecoil, Reco_TauA1_Phi_wRecoil, Reco_TauA1_Eta_wRecoil;
  std::vector<TH1D> Reco_TauMu_P_wRecoil, Reco_TauMu_Pt_wRecoil, Reco_TauMu_Px_wRecoil, Reco_TauMu_Py_wRecoil, Reco_TauMu_Pz_wRecoil, Reco_TauMu_Phi_wRecoil, Reco_TauMu_Eta_wRecoil;
  std::vector<TH1D> Reco_EventFit_Solution_wRecoil;
  std::vector<TH1D> Reco_PtRes_TauA1_wRecoil, Reco_PtRes_TauMu_wRecoil, Reco_dPhi_TauMuTauA1_AfterFit_wRecoil, Reco_dPhi_TauMuTauA1_BeforeFit_wRecoil;
  std::vector<TH1D> Reco_Chi2_FitSolutionOnly_wRecoil, Reco_Chi2_Full_wRecoil, Reco_Chi2_Orig_wRecoil, Reco_Chi2_SC_wRecoil, Reco_Chi2_HC_wRecoil, Reco_Chi2_OrigProb_wRecoil;
  std::vector<TH1D> Reco_Chi2_diff_wRecoil, Reco_Chi2_orig_diff_wRecoil;
  std::vector<TH1D> Reco_ZMass_wRecoil, Reco_NIter_wRecoil, Reco_Z_Energy_Res_wRecoil;
  std::vector<TH1D> Reco_PtRes_TauA1_wRecoil_PreFit, Reco_PtRes_TauMu_wRecoil_PreFit;
  std::vector<TH1D> Reco_PtRes_TauA1_wRecoil_AmbZero, Reco_PtRes_TauA1_wRecoil_wAmb;
  std::vector<TH1D> Reco_PtRes_TauMu_wRecoil_AmbZero, Reco_PtRes_TauMu_wRecoil_wAmb;
  std::vector<TH1D> Reco_Z_Px_wRecoil, Reco_Z_Py_wRecoil, Reco_Z_Pz_wRecoil, Reco_Z_Pt_wRecoil, Reco_Z_PtRes_wRecoil, Reco_Z_Phi_wRecoil, Reco_Z_PhiRes_wRecoil, Reco_Z_Eta_wRecoil, Reco_Z_EtaRes_wRecoil;

  std::vector<TH2D> Reco_TauA1_ptRes_vs_ptGen_wRecoil, Reco_TauA1_ptRes_vs_ptReco_wRecoil, Reco_TauMu_ptRes_vs_ptGen_wRecoil, Reco_TauMu_ptRes_vs_ptReco_wRecoil;
  std::vector<TH2D> Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr, Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr;
  std::vector<TH2D> Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_wAmb, Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_AmbZero, Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_wAmb, Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_AmbZero;
  std::vector<TH2D> Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit, Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit, Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit;
  std::vector<TH2D> Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit, Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit, Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit;

  std::vector<TH2D> Reco_FitSolution_byChi2_Full_vs_RightSolution_wRecoil;
  std::vector<TH2D> Reco_FitSolution_byChi2_orig_vs_RightSolution_wRecoil, Reco_FitSolution_byChi2_SC_vs_RightSolution_wRecoil;
  std::vector<TH2D> Reco_FitSolution_byChi2_HC_vs_RightSolution_wRecoil, Reco_FitSolution_byChi2_origplusSC_vs_RightSolution_wRecoil;

  std::vector<TH1D> TauMu_Start_Collinear_PtRes, TauMu_Start_MET_PtRes, TauMu_Start_PtBalance_PtRes, TauMu_Start_EventRecoil_PtRes;
  std::vector<TH2D> TauMu_Start_Collinear_PtReco_vs_PtGen, TauMu_Start_MET_PtReco_vs_PtGen, TauMu_Start_PtBalance_PtReco_vs_PtGen, TauMu_Start_EventRecoil_PtReco_vs_PtGen;

  std::vector<TH1D> TauMu_Start_MET_PtRes_AfterMC, TauMu_Start_PtBalance_PtRes_AfterMC, TauMu_Start_EventRecoil_PtRes_AfterMC;
  std::vector<TH1D> TauMu_Start_dPhi_TauMuTauH_MET, TauMu_Start_dPhi_TauMuTauH_EventRecoil, TauMu_Start_dPhi_TauMuTauH_PtBalance;
  std::vector<TH1D> Z_Start_MET_PtRes, Z_Start_PtBalance_PtRes, Z_Start_EventRecoil_PtRes, Z_Start_MET_PhiRes, Z_Start_PtBalance_PhiRes, Z_Start_EventRecoil_PhiRes;

  std::vector<TH1D> Z_Start_MET_noNu_PtRes, Z_Start_MET_noNu_PhiRes, Z_Start_MET_noNu_M, Z_Start_MET_noNu_dE;
  std::vector<TH1D> Z_Start_MET_noNu_PtRes_withRot, Z_Start_MET_noNu_PhiRes_withRot, Z_Start_MET_noNu_M_withRot;
  std::vector<TH1D> Z_Start_B2B_M, Z_Start_B2B_dE;

  std::vector<TH1D> Mu_TP_phi0, Mu_TP_lambda, Mu_TP_dxy, Mu_TP_dz, Mu_TP_kappa, Mu_TP_POCA_quadrant, Mu_TP_POCA_quadrantVlad, Mu_TP_POCA_quadrantby_dxyphi0;
  std::vector<TH1D> Mu_TP_Poca_quadrantData;
  std::vector<TH1D> Mu_TP_Poca_quadrantMCDY;
  std::vector<TH2D> Mu_TP_Poca_xy, Mu_TP_Vertex_xy, Mu_TP_RefitVertex_xy, Mu_TP_BeamSpot_xy, Mu_TP_NTP_Poca_xy;

  std::vector<TH1D> MVAMET_metobject_XX, MVAMET_ptobject_XX;
  std::vector<TH1D> A1_Par_Px, A1_Par_Py, A1_Par_Pz, A1_Par_M, A1_Cov_Pxx, A1_Cov_Pyy, A1_Cov_Pzz, A1_Cov_Pxy, A1_Cov_Pxz, A1_Cov_Pyz,A1_Pull_Px, A1_Pull_Py, A1_Pull_Pz, A1_Pull_M;
  std::vector<TH1D> A1_Pull_Px_AmbZero, A1_Pull_Py_AmbZero, A1_Pull_Pz_AmbZero, A1_Pull_M_AmbZero, A1_Pull_Px_wAmb, A1_Pull_Py_wAmb, A1_Pull_Pz_wAmb, A1_Pull_M_wAmb;
  std::vector<TH1D> SV_Par_x, SV_Par_y, SV_Par_z, SV_Cov_xx, SV_Cov_yy, SV_Cov_zz, SV_Cov_xy, SV_Cov_xz, SV_Cov_yz, SV_Pull_Px, SV_Pull_Py, SV_Pull_Pz;
  std::vector<TH1D> SV_Pull_Px_AmbZero, SV_Pull_Py_AmbZero, SV_Pull_Pz_AmbZero, SV_Pull_Px_wAmb, SV_Pull_Py_wAmb, SV_Pull_Pz_wAmb;
  std::vector<TH1D> TauA1_Par_Px_AmbZero, TauA1_Par_Py_AmbZero, TauA1_Par_Pz_AmbZero, TauA1_Cov_Pxx_AmbZero, TauA1_Cov_Pyy_AmbZero, TauA1_Cov_Pzz_AmbZero, TauA1_Cov_Pxy_AmbZero, TauA1_Cov_Pxz_AmbZero, TauA1_Cov_Pyz_AmbZero, TauA1_Pull_Px_AmbZero, TauA1_Pull_Py_AmbZero, TauA1_Pull_Pz_AmbZero;
  std::vector<TH1D> TauA1_Par_Px_CorrectAmb, TauA1_Par_Py_CorrectAmb, TauA1_Par_Pz_CorrectAmb, TauA1_Cov_Pxx_CorrectAmb, TauA1_Cov_Pyy_CorrectAmb, TauA1_Cov_Pzz_CorrectAmb, TauA1_Cov_Pxy_CorrectAmb, TauA1_Cov_Pxz_CorrectAmb, TauA1_Cov_Pyz_CorrectAmb, TauA1_Pull_Px_CorrectAmb, TauA1_Pull_Py_CorrectAmb, TauA1_Pull_Pz_CorrectAmb;
  std::vector<TH1D> TauA1_Par_Px_WrongAmb, TauA1_Par_Py_WrongAmb, TauA1_Par_Pz_WrongAmb, TauA1_Cov_Pxx_WrongAmb, TauA1_Cov_Pyy_WrongAmb, TauA1_Cov_Pzz_WrongAmb, TauA1_Cov_Pxy_WrongAmb, TauA1_Cov_Pxz_WrongAmb, TauA1_Cov_Pyz_WrongAmb, TauA1_Pull_Px_WrongAmb, TauA1_Pull_Py_WrongAmb, TauA1_Pull_Pz_WrongAmb;

  std::vector<TH1D> TauFLSigmaVlad, TauFLSigmaVlad_PhiA1, TauFLSigmaVlad_PhiTau, TauFLSigmaVlad_PhiTauNoCorr, TauFLSigmaVlad_PhiZnoCorr, TauFLSigmaVlad_PhiZwCorr;
  std::vector<TH1D> TauFLSigmaAlex, TauFLSigmaAlex_PhiA1, TauFLSigmaAlex_PhiTau, TauFLSigmaAlex_PhiTauNoCorr, TauFLSigmaAlex_PhiZnoCorr, TauFLSigmaAlex_PhiZwCorr;

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
  bool selectMuon_DiMuonVeto(unsigned i, unsigned i_vtx);

  double Reconstruct_hadronicTauEnergy(unsigned i);
  double GJAngleMax(TLorentzVector A1);
  LorentzVectorParticle CorrectRecoTauMomentumBias(LorentzVectorParticle RecoTau, TLorentzVector RecoA1, std::vector<double> BiasInGJAngleBins);
  TLorentzVector BoostToRestFrame(TLorentzVector TLV1, TLorentzVector TLV2);
  TLorentzVector TauHelixP4AtSV(unsigned int selTau, TLorentzVector Tau);
  TLorentzVector GenTauHelixP4AtSV(unsigned int selTau, TLorentzVector Tau);
  TVector2 ZPtCollinearTauMuEstimator(TrackParticle Muon, TLorentzVector Tauh, double PhiRecoil);
  TLorentzVector TauMuEstimator(TLorentzVector Tauh, TLorentzVector Muon);
  TLorentzVector TauMuEstimator2(TrackParticle Muon, TLorentzVector Tauh, double PhiRecoil);
  TLorentzVector TauMuEstimatorNoZMass(TrackParticle Muon, TLorentzVector Tauh, double PhiRecoil);
  TLorentzVector TauMuFullEstimate(TVector3 PV, TrackParticle Muon, LorentzVectorParticle Tauh, TVector2 TauMuPt, TVector3 &Intersection);
  TLorentzVector ZEstimatorWithMET(bool rotate, TVector3 PV, TrackParticle Muon, LorentzVectorParticle A1, LorentzVectorParticle Tauh, TVector2 MET);
  TLorentzVector ZEstimatorB2B(bool rotate, TVector3 PV, TrackParticle Muon, LorentzVectorParticle A1, LorentzVectorParticle Tauh);

 private:

};
#endif
