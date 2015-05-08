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
  std::vector<TH1D> Gen_DiTau_dPhi, Gen_DiTau_Pt, Gen_Z_Pt, Gen_Z_M, Gen_DiTau_PtBalance_M, Gen_DiTau_PtBalance_dM, Gen_TauMu_PtBalance_Pt, Gen_TauMu_PtBalance_dP, Gen_TauA1_dP;
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
  std::vector<TH1D> TPTF_TauA1_p_Reco, TPTF_TauA1_pt_Reco, TPTF_TauA1_px_Reco, TPTF_TauA1_py_Reco, TPTF_TauA1_pz_Reco, TPTF_TauA1_p_Gen, TPTF_TauA1_pt_Gen, TPTF_TauA1_px_Gen, TPTF_TauA1_py_Gen, TPTF_TauA1_pz_Gen;
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

  std::vector<TH1D> Reco_ZMass, Reco_ZMass_UnboostedGenZ, Reco_EventFit_Solution, Reco_A1Fit_Solution, Reco_Chi2, Reco_Chi2_FitSolutionOnly, Reco_Chi2_FitSolutionOnlyLargeScale, Reco_ConstrainedDeltaSum, Reco_ConstrainedDeltaMass, Reco_ConstrainedDeltaPt, Reco_NIter;
  std::vector<TH1D> Reco_Z_Energy_Res, RecoZ_Pt;
  std::vector<TH1D> GenReco_ZMass, GenReco_EventFit_Solution, GenReco_A1Fit_Solution, GenReco_Chi2, GenReco_Chi2_FitSolutionOnly, GenReco_ConstrainedDeltaSum, GenReco_NIter;
  std::vector<TH1D> Reco_PtRes_TauA1, Reco_PtRes_TauA1_AmbPoint0, Reco_PtRes_TauA1_AmbPoint12, Reco_PtRes_TauA1_AmbPoint1, Reco_PtRes_TauMu, Reco_PtRes_TauMu_AmbPoint0, Reco_PtRes_TauMu_AmbPoint12, Reco_PtRes_TauMu_AmbPoint1;
  std::vector<TH1D> Reco_PtRes_TauA1_LowZPt, Reco_PtRes_TauMu_LowZPt;
  std::vector<TH1D> Reco_PtRes_TauA1_NoFit, Reco_PtRes_TauA1_AmbPoint0_NoFit, Reco_PtRes_TauA1_AmbPoint12_NoFit, Reco_PtRes_TauMu_NoFit, Reco_PtRes_TauMu_AmbPoint0_NoFit, Reco_PtRes_TauMu_AmbPoint12_NoFit;
  std::vector<TH1D> Reco_TauMu_DeltaPX_FitImpact, Reco_TauMu_DeltaPY_FitImpact, Reco_TauMu_DeltaPZ_FitImpact;
  std::vector<TH1D> Reco_TauA1_DeltaPX_FitImpact, Reco_TauA1_DeltaPY_FitImpact, Reco_TauA1_DeltaPZ_FitImpact;
  std::vector<TH1D> Reco_TauMu_ResCosTheta, Reco_TauMu_ResPhi;
  std::vector<TH1D> Reco_dPhi_TauMuTauA1_AfterFit, Reco_dPhi_TauMuTauA1_BeforeFit;
  std::vector<TH1D> Reco_ZMass_MassScan, Reco_ZMass_MassScanUnboosted, Reco_ZMasswithProbWeight_MassScan, Reco_ProbStack_MassScan, Reco_ZMass_PDF;
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
 private:

};
#endif
