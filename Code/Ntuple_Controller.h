//Ntuple_Controller.h HEADER FILE

#ifndef Ntuple_Controller_h
#define Ntuple_Controller_h


// Root include files
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TVectorT.h"
#include "TSystem.h"

// Include files (C & C++ libraries)
#include<iostream>
#include <vector>
#include <string.h>

#include "NtupleReader.h"

#include "HistoConfig.h"
#ifdef USE_TauSpinner
#include "TauSpinerInterface.h"
#endif
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay.h"

#ifdef USE_SVfit
#include "DataFormats/SVFitObject.h"
#include "SVFitStorage.h"
#include "SVfitProvider.h"
#endif


#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"

// Rochester muon momentum correction
#include "CommonFiles/rochcor2012jan22.h"

// small struct needed to allow sorting indices by some value
struct sortIdxByValue {
    bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
        return left.second > right.second;
    }
};

///////////////////////////////////////////////////////////////////////////////
//*****************************************************************************
//*
//*   Class: Ntuple_Controller
//*   
//*   Purpose: The purpose of this class is to provide a interface to the
//*            Ntuple
//*
//*   Designed by: Ian Nugent
//*
//*
//*****************************************************************************
///////////////////////////////////////////////////////////////////////////////


class Ntuple_Controller{
 private:
  NtupleReader *Ntp;
  TFile *newfile;
  TTree *SkimmedTree;
  int nbytes;
  int jentry;
  int nb;
  bool copyTree;

  int currentEvent;

  bool cannotObtainHiggsMass; // avoid repeated printing of warning when running locally

  // Ntuple Access Functions
  virtual void Branch_Setup(TString B_Name, int type);
  virtual void Branch_Setup(){}

  // Functions to configure objects
  virtual void ConfigureObjects(); 
  void doElectrons();
  void doPhotons();
  void doJets();
  void doMuons();
  void doTaus();
  void doMET();
  unsigned int ObjEvent;

  // helper functions for internal calculations
  void printMCDecayChain(unsigned int par, unsigned int level = 0, bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false);

  // Object Variables
  std::vector<TLorentzVector> electrons_default;
  std::vector<TLorentzVector> photons_default;
  std::vector<TLorentzVector> jets_default;
  std::vector<TLorentzVector> muons_default;
  std::vector<TLorentzVector> taus_default;
  TLorentzVector              met_default;
  std::vector<TLorentzVector> electrons;
  std::vector<TLorentzVector> photons;
  std::vector<TLorentzVector> jets;
  std::vector<TLorentzVector> muons;
  std::vector<TLorentzVector> taus;
  TLorentzVector              met;

  // TString flags for object corrections
  TString tauCorrection;
  TString muonCorrection;
  TString elecCorrection;
  TString jetCorrection;

  // Systematic controls variables
  int theSys;
  HistoConfig HConfig;

  // Interfaces
#ifdef USE_TauSpinner  
  TauSpinerInterface TauSpinerInt;
#endif
  HistoConfig HistoC;

  // Fit Variables
  LorentzVectorParticle               theTau;
  std::vector<LorentzVectorParticle>  daughter;
  double                              LC_chi2;
  double                              ndof;
  bool                                fitStatus;
  bool                                isInit;

  // muon correction related objects
  rochcor2012*   rmcor;
  std::vector<TLorentzVector> Muon_corrected_p4;
  void           CorrectMuonP4();
  bool           Muon_isCorrected;

  // helpers for SVFit
#ifdef USE_SVfit
  // create SVFitObject from standard muon and standard tau_h
  void runAndSaveSVFit_MuTauh(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned muIdx, unsigned tauIdx, double scaleMu, double scaleTau, bool save = true);
  // create SVFitObject from standard muon and fully reconstructed 3prong tau
  void runAndSaveSVFit_MuTau3p(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned muIdx, TLorentzVector tauLV, LorentzVectorParticle neutrino, double scaleMu, double scaleTau, bool save = true);
#endif

 public:
  // Constructor
  Ntuple_Controller(std::vector<TString> RootFiles);

  // Destructor
  virtual ~Ntuple_Controller() ;

  // Event initializer
  void InitEvent();

  //TauSpiner function
  double TauSpinerGet(int SpinType);
  void TauSpinerSetSignal(int signalcharge){
#ifdef USE_TauSpinner
TauSpinerInt.SetTauSignalCharge(signalcharge);
#endif
}
   enum beamspot{BS_x0,BS_y0,BS_z0,BS_sigmaZ,BS_dxdz,BS_dydz,BS_BeamWidthX,NBS_par};
   enum TrackQuality {
     undefQuality = -1, loose = 0, tight = 1, highPurity = 2,
     confirmed = 3, goodIterative = 4, looseSetWithPV = 5, highPuritySetWithPV = 6,
     qualitySize = 7
   };
  enum TrackPar{i_qoverp = 0, i_lambda, i_phi, i_dxy,i_dsz};

  // access to SVFit
  #ifdef USE_SVfit
  SVFitObject* getSVFitResult_MuTauh(SVFitStorage& svFitStor, TString metType, unsigned muIdx, unsigned tauIdx, unsigned rerunEvery = 5000, TString suffix = "", double scaleMu = 1 , double scaleTau = 1);
  SVFitObject* getSVFitResult_MuTau3p(SVFitStorage& svFitStor, TString metType, unsigned muIdx, TLorentzVector tauLV, LorentzVectorParticle neutrino, TString suffix = "", double scaleMu = 1, double scaleTau = 1);
  #endif


  // Ntuple Access Functions 
  virtual Int_t Get_Entries();
  virtual void Get_Event(int _jentry);
  virtual Int_t Get_EventIndex();
  virtual TString Get_File_Name();

  //Ntuple Cloning Functions
  virtual void CloneTree(TString n);
  virtual void SaveCloneTree();
  inline void AddEventToCloneTree(){if(copyTree)SkimmedTree->Fill();}

  // Systematic controls
  enum    Systematic {Default=0,NSystematics};

  int     SetupSystematics(TString sys_);
  void    SetSysID(int sysid){theSys=sysid;}


  // Data/MC switch and thin
  bool isData(){return Ntp->Event_isRealData;}
  void ThinTree();

  // Set object corrections to be applied
  void SetTauCorrections(TString tauCorr){tauCorrection = tauCorr;}
  void SetMuonCorrections(TString muonCorr){muonCorrection = muonCorr;}
  void SetElecCorrections(TString elecCorr){elecCorrection = elecCorr;}
  void SetJetCorrections(TString jetCorr){jetCorrection = jetCorr;}
  // corresponding getters
  const TString& GetTauCorrections() const {return tauCorrection;}
  const TString& GetMuonCorrections() const {return muonCorrection;}
  const TString& GetElecCorrections() const {return elecCorrection;}
  const TString& GetJetCorrections() const {return jetCorrection;}

  // Information from input Ntuple path name
  TString GetInputNtuplePath();
  TString GetInputDatasetName();
  TString GetInputPublishDataName();
  int getSampleHiggsMass();
  int readHiggsMassFromString(TString input);

  // resonance mass
  int getHiggsSampleMassFromGenInfo();
  double getResonanceMassFromGenInfo(bool useZ0 = true, bool useHiggs0 = true, bool useW = true);

  // Physics Variable Get Functions
  // Event Variables
  int64_t GetMCID();
  int GetStrippedMCID();
  unsigned int RunNumber(){return Ntp->Event_RunNumber;}
  unsigned int EventNumber(){ return Ntp->Event_EventNumber;}
  int BunchCrossing(){ return Ntp->Event_bunchCrossing;}
  int OrbitNumber(){ return Ntp->Event_orbitNumber;}
  unsigned int LuminosityBlock(){return Ntp->Event_luminosityBlock;}
  float           PileupInfo_TrueNumInteractions_nm1(){return Ntp->PileupInfo_TrueNumInteractions_nm1;}
  float           PileupInfo_TrueNumInteractions_n0(){return Ntp->PileupInfo_TrueNumInteractions_n0;}
  float           PileupInfo_TrueNumInteractions_np1(){return Ntp->PileupInfo_TrueNumInteractions_np1;}
  float        PUWeight(){return Ntp->PUWeight;}
  float        PUWeight_p5(){return Ntp->PUWeight_p5;}
  float        PUWeight3D_m5(){return Ntp->PUWeight3D_m5;}
  float        PUWeight3D(){return Ntp->PUWeight3D;}
  float        PUWeight3D_p5(){return Ntp->PUWeight3D_p5;}
  float        PUWeight_m5(){return Ntp->PUWeight_m5;}
  float		   PUWeightFineBins(){return Ntp->PUWeightFineBins;}

  // embedding
  float		Embedding_TauSpinnerWeight(){return Ntp->TauSpinnerWeight;};
  float 	Embedding_SelEffWeight(){return Ntp->SelEffWeight;}
  float 	Embedding_MinVisPtFilter(){return Ntp->MinVisPtFilter;}
  float 	Embedding_KinWeightPt(){return Ntp->KinWeightPt;}
  float 	Embedding_KinWeightEta(){return Ntp->KinWeightEta;}
  float 	Embedding_KinWeightMassPt(){return Ntp->KinWeightMassPt;}
  // don't use this combined weight blindly. Check which weights you should apply.
  float		EmbeddedWeight(){
	  return Ntp->TauSpinnerWeight * Ntp->SelEffWeight * Ntp->MinVisPtFilter * Ntp->KinWeightPt * Ntp->KinWeightEta * Ntp->KinWeightMassPt;
  }

  TVectorT<double>      beamspot_par(){TVectorT<double> BS(NBS_par);for(unsigned int i=0;i<NBS_par;i++)BS(i)=Ntp->beamspot_par->at(i);return BS;}

  TMatrixTSym<double>   beamspot_cov(){
    TMatrixTSym<double> BS_cov(NBS_par);
    unsigned int l=0;
    for(unsigned int i=0;i<NBS_par;i++){
      for(unsigned int j=i;j<NBS_par;j++){
	BS_cov(i,j)=Ntp->beamspot_cov->at(l);
      }
    }
    return BS_cov;
  }
  
  double  beamspot_emittanceX(){return Ntp->beamspot_emittanceX;}
  double  beamspot_emittanceY(){return Ntp->beamspot_emittanceY;}
  double  beamspot_betaStar(){return Ntp->beamspot_betaStar;}

  // Vertex Information
  unsigned int NVtx(){return Ntp->Vtx_ndof->size();}
  TVector3     Vtx(unsigned int i){return TVector3(Ntp->Vtx_x->at(i),Ntp->Vtx_y->at(i),Ntp->Vtx_z->at(i));}
  double       Vtx_chi2(unsigned int i){return Ntp->Vtx_chi2->at(i);}
  unsigned     Vtx_nTrk(unsigned int i){return Ntp->Vtx_nTrk->at(i);}
  float        Vtx_ndof(unsigned int i){return Ntp->Vtx_ndof->at(i);}
  TMatrixF     Vtx_Cov(unsigned int i);
  std::vector<int>  Vtx_Track_idx(unsigned int i){return Ntp->Vtx_Track_idx->at(i);}
  bool Vtx_isFake(unsigned int i){return Ntp->Vtx_isFake->at(i);}
  TLorentzVector Vtx_TracksP4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->Vtx_TracksP4->at(i).at(j).at(1),Ntp->Vtx_TracksP4->at(i).at(j).at(2),Ntp->Vtx_TracksP4->at(i).at(j).at(3),Ntp->Vtx_TracksP4->at(i).at(j).at(0));}

  bool isVtxGood(unsigned int i);
  bool isGoodVtx(unsigned int i);

  // Muon information
  unsigned int   NMuons(){return Ntp->Muon_p4->size();}
  TLorentzVector Muon_p4(unsigned int i, TString corr = "default");
  TVector3       Muon_Poca(unsigned int i){return TVector3(Ntp->Muon_Poca->at(i).at(0),Ntp->Muon_Poca->at(i).at(1),Ntp->Muon_Poca->at(i).at(2));}
  bool           Muon_isGlobalMuon(unsigned int i){return Ntp->Muon_isGlobalMuon->at(i);}
  bool           Muon_isStandAloneMuon(unsigned int i){return Ntp->Muon_isStandAloneMuon->at(i);}
  bool           Muon_isTrackerMuon(unsigned int i){return Ntp->Muon_isTrackerMuon->at(i);}
  bool           Muon_isCaloMuon(unsigned int i){return Ntp->Muon_isCaloMuon->at(i);}
  bool           Muon_isIsolationValid(unsigned int i){return Ntp->Muon_isIsolationValid->at(i);}
  bool           Muon_isQualityValid(unsigned int i){return Ntp->Muon_isQualityValid->at(i);}
  bool           Muon_isTimeValid(unsigned int i){return Ntp->Muon_isTimeValid->at(i);}
  float          Muon_emEt03(unsigned int i){return Ntp->Muon_emEt03->at(i);}
  float          Muon_emVetoEt03(unsigned int i){return Ntp->Muon_emVetoEt03->at(i);}
  float          Muon_hadEt03(unsigned int i){return Ntp->Muon_hadEt03->at(i);}
  float          Muon_hadVetoEt03(unsigned int i){return Ntp->Muon_hadVetoEt03->at(i);}
  int 	         Muon_nJets03(unsigned int i){return Ntp->Muon_nJets03->at(i);}
  int            Muon_nTracks03(unsigned int i){return Ntp->Muon_nTracks03->at(i);}
  float          Muon_sumPt03(unsigned int i){return Ntp->Muon_sumPt03->at(i);}
  float          Muon_trackerVetoPt03(unsigned int i){return Ntp->Muon_trackerVetoPt03->at(i);}
  float          Muon_emEt05(unsigned int i){return Ntp->Muon_emEt05->at(i);}
  float          Muon_emVetoEt05(unsigned int i){return Ntp->Muon_emVetoEt05->at(i);}
  float          Muon_hadEt05(unsigned int i){return Ntp->Muon_hadEt05->at(i);}
  float          Muon_hadVetoEt05(unsigned int i){return Ntp->Muon_hadVetoEt05->at(i);}
  int            Muon_nJets05(unsigned int i){return Ntp->Muon_nJets05->at(i);}
  int            Muon_nTracks05(unsigned int i){return Ntp->Muon_nTracks05->at(i);}
  float          Muon_sumPt05(unsigned int i){return Ntp->Muon_sumPt05->at(i);}
  float          Muon_trackerVetoPt05(unsigned int i){return Ntp->Muon_trackerVetoPt05->at(i);}
  unsigned int   Muon_Track_idx(unsigned int i){return Ntp->Muon_Track_idx->at(i);}
  int            Muon_hitPattern_pixelLayerwithMeas(unsigned int i){return Ntp->Muon_hitPattern_pixelLayerwithMeas->at(i);}
  int            Muon_numberOfMatchedStations(unsigned int i){return Ntp->Muon_numberOfMatchedStations->at(i);}
  float          Muon_normChi2(unsigned int i){return Ntp->Muon_normChi2->at(i);}
  int            Muon_hitPattern_numberOfValidMuonHits(unsigned int i){return Ntp->Muon_hitPattern_numberOfValidMuonHits->at(i);}
  int            Muon_innerTrack_numberofValidHits(unsigned int i){return Ntp->Muon_innerTrack_numberofValidHits->at(i);}
  int            Muon_numberOfMatches(unsigned int i){return Ntp->Muon_numberOfMatches->at(i);}
  int            Muon_numberOfChambers(unsigned int i){return Ntp->Muon_numberOfChambers->at(i);}
  int            Muon_Charge(unsigned int i){return Ntp->Muon_charge->at(i);}
  int            Muon_trackCharge(unsigned int i){return Ntp->Muon_trackCharge->at(i);}

  bool           Muon_isPFMuon(unsigned int i){return Ntp->Muon_isPFMuon->at(i);}                                                     
  float          Muon_sumChargedHadronPt03(unsigned int i){return Ntp->Muon_sumChargedHadronPt03->at(i);}                             // sum-pt of charged Hadron					                               
  float          Muon_sumChargedParticlePt03(unsigned int i){return Ntp->Muon_sumChargedParticlePt03->at(i);}			      // sum-pt of charged Particles(inludes e/mu)			    
  float          Muon_sumNeutralHadronEt03(unsigned int i){return Ntp->Muon_sumNeutralHadronEt03->at(i);}			      // sum pt of neutral hadrons					    
  float          Muon_sumNeutralHadronEtHighThreshold03(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold03->at(i);}   // sum pt of neutral hadrons with a higher threshold		    
  float          Muon_sumPhotonEt03(unsigned int i){return Ntp->Muon_sumPhotonEt03->at(i);}					      // sum pt of PF photons					    
  float          Muon_sumPhotonEtHighThreshold03(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold03->at(i);}		      // sum pt of PF photons with a higher threshold		    
  float          Muon_sumPUPt03(unsigned int i){return Ntp->Muon_sumPUPt03->at(i);}						      // sum pt of charged Particles not from PV (for Pu corrections)  
																      								    
  float          Muon_sumChargedHadronPt04(unsigned int i){return Ntp->Muon_sumChargedHadronPt04->at(i);}			      // sum-pt of charged Hadron					    
  float          Muon_sumChargedParticlePt04(unsigned int i){return Ntp->Muon_sumChargedParticlePt04->at(i);}			      // sum-pt of charged Particles(inludes e/mu)			    
  float          Muon_sumNeutralHadronEt04(unsigned int i){return Ntp->Muon_sumNeutralHadronEt04->at(i);}			      // sum pt of neutral hadrons					    
  float          Muon_sumNeutralHadronEtHighThreshold04(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold04->at(i);}   // sum pt of neutral hadrons with a higher threshold		    
  float          Muon_sumPhotonEt04(unsigned int i){return Ntp->Muon_sumPhotonEt04->at(i);}					      // sum pt of PF photons					    
  float          Muon_sumPhotonEtHighThreshold04(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold04->at(i);}		      // sum pt of PF photons with a higher threshold		    
  float          Muon_sumPUPt04(unsigned int i){return Ntp->Muon_sumPUPt04->at(i);}						      // sum pt of charged Particles not from PV (for Pu corrections)  

  int            Muon_numberofValidPixelHits(unsigned int i){return Ntp->Muon_numberofValidPixelHits->at(i);}
  int            Muon_trackerLayersWithMeasurement(unsigned int i){return Ntp->Muon_trackerLayersWithMeasurement->at(i);}


  TrackParticle Muon_TrackParticle(unsigned int i){
    TMatrixT<double>    mu_par(TrackParticle::NHelixPar,1);
    TMatrixTSym<double> mu_cov(TrackParticle::NHelixPar);
    unsigned int l=0;
    for(int k=0; k<TrackParticle::NHelixPar; k++){
      mu_par(k,0)=Ntp->Muon_par->at(i).at(k);
      for(int j=k; j<TrackParticle::NHelixPar; j++){
        mu_cov(k,j)=Ntp->Muon_cov->at(i).at(l);
        l++;
      }
    }
    return TrackParticle(mu_par,mu_cov,Ntp->Muon_pdgid->at(i),Ntp->Muon_M->at(i),Ntp->Muon_charge->at(i),Ntp->Muon_B->at(i));
  }

  bool           isGoodMuon(unsigned int i);
  bool           isGoodMuon_nooverlapremoval(unsigned int i);

  bool			 isTightMuon(unsigned int i);
  bool			 isTightMuon(unsigned int i, unsigned int j, TString corr = "default");
  bool           isSelectedMuon(unsigned int i, unsigned int j, double impact_xy, double impact_z, TString corr = "default");
  bool			 isLooseMuon(unsigned int i);
  float          Muon_RelIso(unsigned int i, TString corr = "default");

  //Base Tau Information (PF)
   unsigned int      NPFTaus(){return Ntp->PFTau_p4->size();}
   TLorentzVector	 PFTau_p4(unsigned int i, TString corr = "default");
   TVector3          PFTau_Poca(unsigned int i){return TVector3(Ntp->PFTau_Poca->at(i).at(0),Ntp->PFTau_Poca->at(i).at(1),Ntp->PFTau_Poca->at(i).at(2));}
   bool PFTau_isTightIsolation(unsigned int i){return Ntp->PFTau_isTightIsolation->at(i);}
   bool PFTau_isMediumIsolation(unsigned int i){return  Ntp->PFTau_isMediumIsolation->at(i);}
   bool PFTau_isLooseIsolation(unsigned int i){return  Ntp->PFTau_isLooseIsolation->at(i);}
   bool PFTau_isTightIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isTightIsolationDBSumPtCorr->at(i);}
   bool PFTau_isMediumIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isMediumIsolationDBSumPtCorr->at(i);}
   bool PFTau_isLooseIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isLooseIsolationDBSumPtCorr->at(i);}
   bool PFTau_isVLooseIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isVLooseIsolationDBSumPtCorr->at(i);}
   bool PFTau_isHPSAgainstElectronsLoose(unsigned int i){return  Ntp->PFTau_isHPSAgainstElectronsLoose->at(i);}
   bool PFTau_isHPSAgainstElectronsMedium(unsigned int i){return  Ntp->PFTau_isHPSAgainstElectronsMedium->at(i);}
   bool PFTau_isHPSAgainstElectronsTight(unsigned int i){return  Ntp->PFTau_isHPSAgainstElectronsTight->at(i);}
   bool PFTau_isHPSAgainstMuonLoose(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonLoose->at(i);}
   bool PFTau_isHPSAgainstMuonMedium(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonMedium->at(i);}
   bool PFTau_isHPSAgainstMuonTight(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonTight->at(i);}
   bool PFTau_isHPSAgainstMuonLoose2(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonLoose2->at(i);}
   bool PFTau_isHPSAgainstMuonMedium2(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonMedium2->at(i);}
   bool PFTau_isHPSAgainstMuonTight2(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonTight2->at(i);}
   bool PFTau_isHPSByDecayModeFinding(unsigned int i){return  Ntp->PFTau_isHPSByDecayModeFinding->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection->at(i);}
   bool PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits->at(i);}
   bool PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits->at(i);}
   float PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits->at(i);}
   bool PFTau_HPSPFTauDiscriminationByLooseIsolationMVA(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMediumIsolationMVA(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA->at(i);}
   bool PFTau_HPSPFTauDiscriminationByTightIsolationMVA(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA->at(i);}
   bool PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2->at(i);}
   bool PFTau_HPSPFTauDiscriminationByTightIsolationMVA2(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA2->at(i);}
   int PFTau_hpsDecayMode(unsigned int i){return  Ntp->PFTau_hpsDecayMode->at(i);}
   int PFTau_Charge(unsigned int i){return  Ntp->PFTau_Charge->at(i);}
   std::vector<int> PFTau_Track_idx(unsigned int i){return  Ntp->PFTau_Track_idx->at(i);}
   TVector3 PFTau_TIP_primaryVertex_pos(unsigned int i){return  TVector3(Ntp->PFTau_TIP_primaryVertex_pos->at(i).at(0),Ntp->PFTau_TIP_primaryVertex_pos->at(i).at(1),Ntp->PFTau_TIP_primaryVertex_pos->at(i).at(2));}
   TMatrixTSym<double> PFTau_TIP_primaryVertex_cov(unsigned int i);
   bool PFTau_TIP_hassecondaryVertex(unsigned int i){if(Ntp->PFTau_TIP_secondaryVertex_pos->at(i).size()==3)return true; return false;}
   TVector3 PFTau_TIP_secondaryVertex_pos(unsigned int i){return  TVector3(Ntp->PFTau_TIP_secondaryVertex_pos->at(i).at(0),Ntp->PFTau_TIP_secondaryVertex_pos->at(i).at(1),Ntp->PFTau_TIP_secondaryVertex_pos->at(i).at(2));}
   TMatrixTSym<double> PFTau_TIP_secondaryVertex_cov(unsigned int i);
   double PFTau_TIP_secondaryVertex_vtxchi2(unsigned int i){if(Ntp->PFTau_TIP_secondaryVertex_vtxchi2->at(i).size()==1) return  Ntp->PFTau_TIP_secondaryVertex_vtxchi2->at(i).at(0); return 0;}
   double PFTau_TIP_secondaryVertex_vtxndof(unsigned int i){if(Ntp->PFTau_TIP_secondaryVertex_vtxndof->at(i).size()==1) return  Ntp->PFTau_TIP_secondaryVertex_vtxndof->at(i).at(0);  return 0;}
   bool PFTau_TIP_hasA1Momentum(unsigned int i){if(Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar)return true; return false;}
   LorentzVectorParticle PFTau_a1_lvp(unsigned int i);
   TLorentzVector PFTau_3PS_A1_LV(unsigned int i){return PFTau_a1_lvp(i).LV();}
   int PFTau_NdaughterTracks(unsigned int i){return Ntp->PFTau_daughterTracks_M->at(i).size();}
   std::vector<TrackParticle> PFTau_daughterTracks(unsigned int i);
   int PFTau_NdaughtersReFitTracks_p4(unsigned int i) {return Ntp->PFTau_PionsP4->at(i).size();};
   std::vector<TLorentzVector> PFTau_daughterReFitTracks_p4(unsigned int i);
   std::vector<TVector3> PFTau_daughterTracks_poca(unsigned int i);   
   TMatrixTSym<double> PFTau_FlightLength3d_cov(unsigned int i){return  PFTau_TIP_secondaryVertex_cov(i)+PFTau_TIP_primaryVertex_cov(i);}
   TVector3 PFTau_FlightLength3d(unsigned int i){return PFTau_TIP_secondaryVertex_pos(i)-PFTau_TIP_primaryVertex_pos(i);}
   TMatrixTSym<double> PF_Tau_FlightLegth3d_TauFrame_cov(unsigned int i);
   TVector3 PF_Tau_FlightLegth3d_TauFrame(unsigned int i);
   double   PFTau_FlightLength_significance(unsigned int i){float e=PFTau_FlightLength_error(i); if(e>0) return PFTau_FlightLength(i)/e; return 0;}
   double	PFTau_FlightLenght_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov );
   double   PFTau_FlightLength_error(unsigned int i){return sqrt(PF_Tau_FlightLegth3d_TauFrame_cov(i)(LorentzVectorParticle::vz,LorentzVectorParticle::vz));}
   double   PFTau_FlightLength(unsigned int i){return PFTau_FlightLength3d(i).Mag();}
   
   double   PFTau_Mass(unsigned int i);

   bool ThreeProngTauFit(unsigned int i, unsigned int j,LorentzVectorParticle &theTau,std::vector<LorentzVectorParticle> &daughter,double &LC_chi2, double &phisign){
	   if(Ntp->PFTau_TIP_secondaryVertex_vtxchi2->at(i).size()==1 && Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar){
		   // Tau Solver
		   TVector3 pv=PFTau_TIP_primaryVertex_pos(i);
		   TMatrixTSym<double> pvcov=PFTau_TIP_primaryVertex_cov(i);
		   LorentzVectorParticle a1=PFTau_a1_lvp(i);
		   TauA1NuConstrainedFitter TauA1NU(j,a1,pv,pvcov);
		   bool fitStatus= TauA1NU.Fit();
		   if(fitStatus){
			   phisign =TauA1NU.GetTauRotationSignificance();
			   theTau=TauA1NU.GetMother();
			   //	 std::cout<<" theTau.Px()   "<< theTau.LV().Px() <<"   " <<" theTau.Py()   "<< theTau.LV().Py() <<"   " <<" theTau.Pz()   "<< theTau.LV().Pz() <<"   " <<std::endl;
			   daughter=TauA1NU.GetReFitDaughters();
			   return true;
		   }
	   }
	   return false;
   }

/*    bool ThreeProngTauFit(unsigned int i, unsigned int j,LorentzVectorParticle &theTau,std::vector<LorentzVectorParticle> &daughter,double &LC_chi2){ */
/*      if(Ntp->PFTau_TIP_secondaryVertex_vtxchi2->at(i).size()==1 && Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar){ */
/*        // Tau Solver */
/*        TVector3 pv=PFTau_TIP_primaryVertex_pos(i); */
/*        TMatrixTSym<double> pvcov=PFTau_TIP_primaryVertex_cov(i); */
/*        LorentzVectorParticle a1=PFTau_a1_lvp(i); */
/*        TauA1NuConstrainedFitter TauA1NU(j,a1,pv,pvcov); */
/*        bool fitStatus= TauA1NU.Fit(); */
/*        if(fitStatus){ */
/* 	 theTau=TauA1NU.GetMother(); */
/* 	 daughter=TauA1NU.GetReFitDaughters(); */
/* 	 return  true; */
/*        } */
/*      } */
/*      return false; */
/*    } */
   bool EventFit(unsigned int i, unsigned int MuonIndex, LorentzVectorParticle TauA1, LorentzVectorParticle &theZ,std::vector<LorentzVectorParticle> &RefitDaughters, std::vector<LorentzVectorParticle> &InitDaughters,double &LC_chi2, int &Niterat, double &csum, double MassConstraint){//, TVectorD &par_0, TVectorD &par){
	   TrackParticle Muon = Muon_TrackParticle(MuonIndex);
	   //std::cout << "Muon Track: Kappa:" << Muon.Parameter(TrackParticle::kappa) << std::endl;
	   //std::cout << "Muon Track: Lambda:" << Muon.Parameter(TrackParticle::lambda) << std::endl;
	   //std::cout << "Muon Track: Phi:" << Muon.Parameter(TrackParticle::phi) << std::endl;
	   //std::cout << "Muon Track: dxy:" << Muon.Parameter(TrackParticle::dxy) << std::endl;
	   //std::cout << "Muon Track: dz:" << Muon.Parameter(TrackParticle::dz) << std::endl;
	   //std::cout << "Muon Track: P:" << Muon_p4(MuonIndex).P() << std::endl;
	   //std::cout << "Muon Track: Pt:" << Muon_p4(MuonIndex).Pt() << std::endl;

	   TLorentzVector MuonLV = TLorentzVector(Ntp->Muon_p4->at(MuonIndex).at(1),Ntp->Muon_p4->at(MuonIndex).at(2),Ntp->Muon_p4->at(MuonIndex).at(3),Ntp->Muon_p4->at(MuonIndex).at(0));

	   if(Ntp->PFTau_TIP_secondaryVertex_vtxchi2->at(i).size()==1 &&
		   Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar){
		   // Tau Solver
		   TVector3 pv=PFTau_TIP_primaryVertex_pos(i);
		   TMatrixTSym<double> pvcov=PFTau_TIP_primaryVertex_cov(i);
		   //std::cout << "Print pvcov" << std::endl;
		   //pvcov.Print();
		   TVector3 MuonPoca(Ntp->Muon_Poca->at(MuonIndex).at(0),Ntp->Muon_Poca->at(MuonIndex).at(1),Ntp->Muon_Poca->at(MuonIndex).at(2));
		   DiTauConstrainedFitter Z2Tau(TauA1, Muon, pv, pvcov);
		   Z2Tau.SetMassConstraint(MassConstraint);
		   Z2Tau.SetMaxDelta(1.0);
		   Z2Tau.SetNIterMax(100);
		   Z2Tau.SetEpsilon(0.01);
		   //std::cout << "DEBUG 6" << std::endl;
		   bool fitStatus= Z2Tau.Fit();
		   //std::cout << "DEBUG 7" << std::endl;
		   if(fitStatus && Z2Tau.isConverged()){
			   RefitDaughters = Z2Tau.GetReFitDaughters();
			   InitDaughters = Z2Tau.GetInitialDaughters();
			   //std::cout << "InitDaughters.at(0): "; InitDaughters.at(0).LV().Print();
			   //std::cout << "InitDaughters.at(1): "; InitDaughters.at(1).LV().Print();
			   //std::cout << "RefitDaughters.at(0): "; RefitDaughters.at(0).LV().Print();
			   //std::cout << "RefitDaughters.at(1): "; RefitDaughters.at(1).LV().Print();

			   LorentzVectorParticle TZ=Z2Tau.GetMother();
			   theZ = Z2Tau.GetMother();
			   //std::cout << "Mother: "; theZ.LV().Print();
			   //std::cout << "Mother pt: " << theZ.LV().Pt() << std::endl;
			   LorentzVectorParticle TauMuInitialEstimate =Z2Tau.GetTauMuEstimate();
			   //std::cout << "TauMuInitialEstimate: "; TauMuInitialEstimate.LV().Print();
			   //par_0 = Z2Tau.Getpar_0();
			   //par = Z2Tau.Getpar();
			   daughter=Z2Tau.GetReFitDaughters();
			   LC_chi2=Z2Tau.ChiSquare();
			   Niterat = Z2Tau.NIter();
			   csum = Z2Tau.CSum();

			   return  true;
		   }
	   }
	   return false;
   }

   bool AmbiguitySolver(std::vector<bool> A1Fit, std::vector<bool> EventFit, std::vector<double> Probs, int &IndexToReturn, bool &AmbuguityPoint);
   bool AmbiguitySolverByChi2(std::vector<bool> A1Fit, std::vector<bool> EventFit, std::vector<double> Chi2s, int &IndexToReturn, bool &AmbuguityPoint);

   ////////////////////////////////////////////////
   // wrapper for backwards compatibility to KFit do not use in new code!!!
   int NKFTau(){return NPFTaus();}
   int KFTau_discriminatorByQC(unsigned int i,unsigned int j=0){
     TLorentzVector a1=PFTau_a1_lvp(i).LV();
     if(0.8<a1.M() && a1.M()<1.5)return true;
     return false;
   }

/*    bool KFTau_discriminatorByKFit(unsigned int i, unsigned int j=0){ */
/*      LorentzVectorParticle theTau; */
/*      std::vector<LorentzVectorParticle> daughter; */
/*      double LC_chi2; */
/*      return ThreeProngTauFit(i,j,theTau,daughter,LC_chi2); */
/*    } */

/*    TLorentzVector KFTau_TauFit_p4(unsigned int i, unsigned int j=0){ */
/*      LorentzVectorParticle theTau; */
/*      std::vector<LorentzVectorParticle> daughter; */
/*      double LC_chi2; */
/*      if(ThreeProngTauFit(i,j,theTau,daughter,LC_chi2))return theTau.LV();return TLorentzVector(0,0,0,0); */
/*    } */
   int KFTau_MatchedHPS_idx(unsigned int i){return i;}
/*    bool isGoodKFTau(unsigned int i,unsigned int j=0){ */
/*      LorentzVectorParticle theTau; */
/*      std::vector<LorentzVectorParticle> daughter; */
/*      double LC_chi2; */
/*      if(ThreeProngTauFit(i,j,theTau,daughter,LC_chi2)){ */
/*        if(KFTau_discriminatorByQC(i,j)){ */
/* 	 if(PFTau_hpsDecayMode(KFTau_MatchedHPS_idx(i)) == 10 && PFTau_isHPSByDecayModeFinding(i)){ */
/* 	   if(PFTau_isVLooseIsolationDBSumPtCorr(KFTau_MatchedHPS_idx(i))){ */
/* 	     return true; */
/* 	   } */
/* 	 } */
/*        } */
/*      } */
/*      return false; */
/*    } */
   
/*    TLorentzVector   KFTau_TauVis_p4(unsigned int i,unsigned int j=0){return PFTau_a1_lvp(i).LV();} */
/*    TLorentzVector   KFTau_Neutrino_p4(unsigned int i,unsigned int j=0){return (KFTau_TauFit_p4(i,j)-PFTau_a1_lvp(i).LV());} */
   int      KFTau_nKinTaus(){return NKFTau();}
   int      KFTau_indexOfFitInfo(unsigned int i){return i;}
   TVector3 KFTau_Fit_TauPrimVtx(unsigned int i){return PFTau_TIP_primaryVertex_pos(i);}
   float    KFTau_Fit_ndf(unsigned int i,unsigned int j=0){return PFTau_TIP_secondaryVertex_vtxndof(i);}
/*    float    KFTau_Fit_chi2(unsigned int i,unsigned int j=0){ */
/*      LorentzVectorParticle theTau; */
/*      std::vector<LorentzVectorParticle> daughter; */
/*      double LC_chi2; */
/*      if(ThreeProngTauFit(i,j,theTau,daughter,LC_chi2)){ */
/*        return PFTau_TIP_secondaryVertex_vtxchi2(i)+LC_chi2; */
/*      }  */
/*      return 999; */
/*    } */
//   float    KFTau_Fit_Chi2Prob(unsigned int i,unsigned int j=0){return TMath::Prob(KFTau_Fit_chi2(i,j),(int)KFTau_Fit_ndf(i,j));}
   int      KFTau_Fit_charge(unsigned int i){if(Ntp->PFTau_a1_charge->at(i).size()>0) return Ntp->PFTau_a1_charge->at(i).at(0); return 0;}
   int      KFTau_Fit_csum(unsigned int i,unsigned int j=0){return 0;}
   int      KFTau_Fit_iterations(unsigned int i, unsigned int j=0){return 0;}
   double   KFTau_RotatedVtx(unsigned int i, unsigned int j=0){return 0;}
   
   double   KFTau_Fit_RefitVisibleMass(unsigned int i, unsigned int j=0){return PFTau_a1_lvp(i).LV().M();}
   double   KFTau_Fit_SV_PV_significance(unsigned int i, unsigned int j=0){return PFTau_FlightLength_significance(i);}
   double   KFTau_Fit_SV_PV_error(unsigned int i, unsigned int j=0){return PFTau_FlightLength_error(i);}
   double   KFTau_Fit_SV_PV_distance(unsigned int i, unsigned int j=0){return PFTau_FlightLength(i);}

   /////////////////////////////////////////////////////////////////////////////////


      // Jet Information
   unsigned int       NPFJets(){return Ntp->PFJet_p4->size();}
   TLorentzVector     PFJet_p4(unsigned int i, TString corr = "default");
   float              PFJet_chargedEmEnergy(unsigned int i){return Ntp->PFJet_chargedEmEnergy->at(i);}
   float              PFJet_chargedHadronEnergy(unsigned int i){return Ntp->PFJet_chargedHadronEnergy->at(i);}
   int	              PFJet_chargedHadronMultiplicity(unsigned int i){return Ntp->PFJet_chargedHadronMultiplicity->at(i);}
   float              PFJet_chargedMuEnergy(unsigned int i){return Ntp->PFJet_chargedMuEnergy->at(i);}
   int	              PFJet_chargedMultiplicity(unsigned int i){return Ntp->PFJet_chargedMultiplicity->at(i);}
   float              PFJet_electronEnergy(unsigned int i){return Ntp->PFJet_electronEnergy->at(i);}
   int	              PFJet_electronMultiplicity(unsigned int i){return Ntp->PFJet_electronMultiplicity->at(i);}
   float              PFJet_HFEMEnergy(unsigned int i){return Ntp->PFJet_HFEMEnergy->at(i);}
   int	              PFJet_HFEMMultiplicity(unsigned int i){return Ntp->PFJet_HFEMMultiplicity->at(i);}
   float              PFJet_HFHadronEnergy(unsigned int i){return Ntp->PFJet_HFHadronEnergy->at(i);}
   int	              PFJet_HFHadronMultiplicity(unsigned int i){return Ntp->PFJet_HFHadronMultiplicity->at(i);}
   float              PFJet_muonEnergy(unsigned int i){return Ntp->PFJet_muonEnergy->at(i);}
   int	              PFJet_muonMultiplicity(unsigned int i){return Ntp->PFJet_muonMultiplicity->at(i);}
   float              PFJet_neutralEmEnergy(unsigned int i){return Ntp->PFJet_neutralEmEnergy->at(i);}
   float              PFJet_neutralHadronEnergy(unsigned int i){return Ntp->PFJet_neutralHadronEnergy->at(i);}
   int	              PFJet_neutralHadronMultiplicity(unsigned int i){return Ntp->PFJet_neutralHadronMultiplicity->at(i);}
   float              PFJet_photonEnergy(unsigned int i){return Ntp->PFJet_photonEnergy->at(i);}
   int	              PFJet_photonMultiplicity(unsigned int i){return Ntp->PFJet_photonMultiplicity->at(i);}
   float              PFJet_jetArea(unsigned int i){return Ntp->PFJet_jetArea->at(i);}
   float              PFJet_maxDistance(unsigned int i){return Ntp->PFJet_maxDistance->at(i);}
   int                PFJet_nConstituents(unsigned int i){return Ntp->PFJet_nConstituents->at(i);}
   float              PFJet_pileup(unsigned int i){return Ntp->PFJet_pileup->at(i);}
   float              PFJet_etaetaMoment(unsigned int i){return Ntp->PFJet_etaetaMoment->at(i);}
   float              PFJet_etaphiMoment(unsigned int i){return Ntp->PFJet_etaphiMoment->at(i);}
   std::vector<int>   PFJet_Track_idx(unsigned int i){return Ntp->PFJet_Track_idx->at(i);}
   int                PFJet_MatchedHPS_idx(unsigned int i){return Ntp->PFJet_MatchedHPS_idx->at(i);}
   int                PFJet_numberOfDaughters(unsigned int i){return Ntp->PFJet_numberOfDaughters->at(i);}
   float              PFJet_chargedEmEnergyFraction(unsigned int i){return Ntp->PFJet_chargedEmEnergyFraction->at(i);}
   float              PFJet_chargedHadronEnergyFraction(unsigned int i){return Ntp->PFJet_chargedHadronEnergyFraction->at(i);}
   float              PFJet_neutralHadronEnergyFraction(unsigned int i){return Ntp->PFJet_neutralHadronEnergyFraction->at(i);}
   float              PFJet_neutralEmEnergyFraction(unsigned int i){return Ntp->PFJet_neutralEmEnergyFraction->at(i);}
   bool               isGoodJet(unsigned int i);
   bool               isGoodJet_nooverlapremoval(unsigned int i);
   bool               isJetID(unsigned int i, TString corr = "default");
   int	              PFJet_nTrk(unsigned int i){return Ntp->PFJet_nTrk->at(i);}
   TLorentzVector     PFJet_TracksP4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->PFJet_TracksP4->at(i).at(j).at(1),Ntp->PFJet_TracksP4->at(i).at(j).at(2),Ntp->PFJet_TracksP4->at(i).at(j).at(3),Ntp->PFJet_TracksP4->at(i).at(j).at(0));}
   int                PFJet_nTracks(unsigned int i){return Ntp->PFJet_TracksP4->at(i).size();}
   float			  PFJet_JECuncertainty(unsigned int i){return Ntp->PFJet_JECuncertainty->at(i);}
   unsigned int       PFJet_NGenJets(){return Ntp->PFJet_GenJet_p4->size();}
   unsigned int       PFJet_NGenJetsNoNu(){return Ntp->PFJet_GenJetNoNu_p4->size();}
   TLorentzVector     PFJet_GenJet_p4(unsigned int i){return TLorentzVector(Ntp->PFJet_GenJet_p4->at(i).at(1),Ntp->PFJet_GenJet_p4->at(i).at(2),Ntp->PFJet_GenJet_p4->at(i).at(3),Ntp->PFJet_GenJet_p4->at(i).at(0));}
   TLorentzVector     PFJet_GenJet_Constituents_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->PFJet_GenJet_Constituents_p4->at(i).at(j).at(1),Ntp->PFJet_GenJet_Constituents_p4->at(i).at(j).at(2),Ntp->PFJet_GenJet_Constituents_p4->at(i).at(j).at(3),Ntp->PFJet_GenJet_Constituents_p4->at(i).at(j).at(0));}
   TLorentzVector     PFJet_GenJetNoNu_p4(unsigned int i){return TLorentzVector(Ntp->PFJet_GenJetNoNu_p4->at(i).at(1),Ntp->PFJet_GenJetNoNu_p4->at(i).at(2),Ntp->PFJet_GenJetNoNu_p4->at(i).at(3),Ntp->PFJet_GenJetNoNu_p4->at(i).at(0));}
   TLorentzVector     PFJet_GenJetNoNu_Constituents_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->PFJet_GenJetNoNu_Constituents_p4->at(i).at(j).at(1),Ntp->PFJet_GenJetNoNu_Constituents_p4->at(i).at(j).at(2),Ntp->PFJet_GenJetNoNu_Constituents_p4->at(i).at(j).at(3),Ntp->PFJet_GenJetNoNu_Constituents_p4->at(i).at(j).at(0));}

   float              PFJet_PUJetID_discr(unsigned int i){return Ntp->PFJet_PUJetID_discr->at(i);}
   bool	              PFJet_PUJetID_looseWP(unsigned int i){return Ntp->PFJet_PUJetID_looseWP->at(i);}
   bool	              PFJet_PUJetID_mediumWP(unsigned int i){return Ntp->PFJet_PUJetID_mediumWP->at(i);}
   bool	              PFJet_PUJetID_tightWP(unsigned int i){return Ntp->PFJet_PUJetID_tightWP->at(i);}

   int	              PFJet_partonFlavour(unsigned int i){return Ntp->PFJet_partonFlavour->at(i);}
   float              PFJet_bDiscriminator(unsigned int i){return Ntp->PFJet_bDiscriminator->at(i);}
   //float              PFJet_BTagWeight(unsigned int i){return Ntp->PFJet_BTagWeight->at(i);} // not implemented at the moment

   double 			  rundependentJetPtCorrection(double jeteta, int runnumber);
   double             JERCorrection(TLorentzVector jet, double dr=0.25, TString corr = "default"); // dr=0.25 from AN2013_416_v4
   TLorentzVector     PFJet_matchGenJet(TLorentzVector jet, double dr);
   double             JetEnergyResolutionCorr(double jeteta);
   double             JetEnergyResolutionCorrErr(double jeteta);

   //MET information
   double             MET_Uncorr_et(){return Ntp->MET_Uncorr_et;}
   double             MET_Uncorr_phi(){return Ntp->MET_Uncorr_phi;}
   double             MET_Uncorr_ex(){return Ntp->MET_Uncorr_et*cos(Ntp->MET_Uncorr_phi);}
   double             MET_Uncorr_ey(){return Ntp->MET_Uncorr_et*sin(Ntp->MET_Uncorr_phi);}
   double             MET_Uncorr_significance(){return Ntp->MET_Uncorr_significance;}
   double             MET_Uncorr_significance_xx(){return Ntp->MET_Uncorr_significance_xx;}
   double             MET_Uncorr_significance_xy(){return Ntp->MET_Uncorr_significance_xy;}
   double             MET_Uncorr_significance_yy(){return Ntp->MET_Uncorr_significance_yy;}

   double             MET_CorrT0rt_et(){return Ntp->MET_CorrT0rt_et;}
   double             MET_CorrT0rt_phi(){return Ntp->MET_CorrT0rt_phi;}
   double             MET_CorrT0rt_ex(){return Ntp->MET_CorrT0rt_et*cos(Ntp->MET_CorrT0rt_phi);}
   double             MET_CorrT0rt_ey(){return Ntp->MET_CorrT0rt_et*sin(Ntp->MET_CorrT0rt_phi);}
   double             MET_CorrT0rtT1_et(){return Ntp->MET_CorrT0rtT1_et;}
   double             MET_CorrT0rtT1_phi(){return Ntp->MET_CorrT0rtT1_phi;}
   double             MET_CorrT0rtT1_ex(){return Ntp->MET_CorrT0rtT1_et*cos(Ntp->MET_CorrT0rtT1_phi);}
   double             MET_CorrT0rtT1_ey(){return Ntp->MET_CorrT0rtT1_et*sin(Ntp->MET_CorrT0rtT1_phi);}
   double             MET_CorrT0pc_et(){return Ntp->MET_CorrT0pc_et;}
   double             MET_CorrT0pc_phi(){return Ntp->MET_CorrT0pc_phi;}
   double             MET_CorrT0pc_ex(){return Ntp->MET_CorrT0pc_et*cos(Ntp->MET_CorrT0pc_phi);}
   double             MET_CorrT0pc_ey(){return Ntp->MET_CorrT0pc_et*sin(Ntp->MET_CorrT0pc_phi);}
   double             MET_CorrT0pcT1_et(){return Ntp->MET_CorrT0pcT1_et;}
   double             MET_CorrT0pcT1_phi(){return Ntp->MET_CorrT0pcT1_phi;}
   double             MET_CorrT0pcT1_ex(){return Ntp->MET_CorrT0pcT1_et*cos(Ntp->MET_CorrT0pcT1_phi);}
   double             MET_CorrT0pcT1_ey(){return Ntp->MET_CorrT0pcT1_et*sin(Ntp->MET_CorrT0pcT1_phi);}
   double             MET_CorrT0rtTxy_et(){return Ntp->MET_CorrT0rtTxy_et;}
   double             MET_CorrT0rtTxy_phi(){return Ntp->MET_CorrT0rtTxy_phi;}
   double             MET_CorrT0rtTxy_ex(){return Ntp->MET_CorrT0rtTxy_et*cos(Ntp->MET_CorrT0rtTxy_phi);}
   double             MET_CorrT0rtTxy_ey(){return Ntp->MET_CorrT0rtTxy_et*sin(Ntp->MET_CorrT0rtTxy_phi);}
   double             MET_CorrT0rtT1Txy_et(){return Ntp->MET_CorrT0rtT1Txy_et;}
   double             MET_CorrT0rtT1Txy_phi(){return Ntp->MET_CorrT0rtT1Txy_phi;}
   double             MET_CorrT0rtT1Txy_ex(){return Ntp->MET_CorrT0rtT1Txy_et*cos(Ntp->MET_CorrT0rtT1Txy_phi);}
   double             MET_CorrT0rtT1Txy_ey(){return Ntp->MET_CorrT0rtT1Txy_et*sin(Ntp->MET_CorrT0rtT1Txy_phi);}
   double             MET_CorrT0pcTxy_et(){return Ntp->MET_CorrT0pcTxy_et;}
   double             MET_CorrT0pcTxy_phi(){return Ntp->MET_CorrT0pcTxy_phi;}
   double             MET_CorrT0pcTxy_ex(){return Ntp->MET_CorrT0pcTxy_et*cos(Ntp->MET_CorrT0pcTxy_phi);}
   double             MET_CorrT0pcTxy_ey(){return Ntp->MET_CorrT0pcTxy_et*sin(Ntp->MET_CorrT0pcTxy_phi);}
   double             MET_CorrT0pcT1Txy_et(){return Ntp->MET_CorrT0pcT1Txy_et;}
   double             MET_CorrT0pcT1Txy_phi(){return Ntp->MET_CorrT0pcT1Txy_phi;}
   double             MET_CorrT0pcT1Txy_ex(){return Ntp->MET_CorrT0pcT1Txy_et*cos(Ntp->MET_CorrT0pcT1Txy_phi);}
   double             MET_CorrT0pcT1Txy_ey(){return Ntp->MET_CorrT0pcT1Txy_et*sin(Ntp->MET_CorrT0pcT1Txy_phi);}
   double             MET_CorrT1_et(){return Ntp->MET_CorrT1_et;}
   double             MET_CorrT1_phi(){return Ntp->MET_CorrT1_phi;}
   double             MET_CorrT1_ex(){return Ntp->MET_CorrT1_et*cos(Ntp->MET_CorrT1_phi);}
   double             MET_CorrT1_ey(){return Ntp->MET_CorrT1_et*sin(Ntp->MET_CorrT1_phi);}
   double             MET_CorrT1Txy_et(){return Ntp->MET_CorrT1Txy_et;}
   double             MET_CorrT1Txy_phi(){return Ntp->MET_CorrT1Txy_phi;}
   double             MET_CorrT1Txy_ex(){return Ntp->MET_CorrT1Txy_et*cos(Ntp->MET_CorrT1Txy_phi);}
   double             MET_CorrT1Txy_ey(){return Ntp->MET_CorrT1Txy_et*sin(Ntp->MET_CorrT1Txy_phi);}
   double             MET_CorrCaloT1_et(){return Ntp->MET_CorrCaloT1_et;}
   double             MET_CorrCaloT1_phi(){return Ntp->MET_CorrCaloT1_phi;}
   double             MET_CorrCaloT1_ex(){return Ntp->MET_CorrCaloT1_et*cos(Ntp->MET_CorrCaloT1_phi);}
   double             MET_CorrCaloT1_ey(){return Ntp->MET_CorrCaloT1_et*sin(Ntp->MET_CorrCaloT1_phi);}
   double             MET_CorrCaloT1T2_et(){return Ntp->MET_CorrCaloT1T2_et;}
   double             MET_CorrCaloT1T2_phi(){return Ntp->MET_CorrCaloT1T2_phi;}
   double             MET_CorrCaloT1T2_ex(){return Ntp->MET_CorrCaloT1T2_et*cos(Ntp->MET_CorrCaloT1T2_phi);}
   double             MET_CorrCaloT1T2_ey(){return Ntp->MET_CorrCaloT1T2_et*sin(Ntp->MET_CorrCaloT1T2_phi);}

   double             MET_CorrMVA_et(){return Ntp->MET_CorrMVA_et;}
   double             MET_CorrMVA_phi(){return Ntp->MET_CorrMVA_phi;}
   double             MET_CorrMVA_ex(){return Ntp->MET_CorrMVA_et*cos(Ntp->MET_CorrMVA_phi);}
   double             MET_CorrMVA_ey(){return Ntp->MET_CorrMVA_et*sin(Ntp->MET_CorrMVA_phi);}
   double             MET_CorrMVA_significance(){return Ntp->MET_CorrMVA_significance;}
   double             MET_CorrMVA_significance_xx(){return Ntp->MET_CorrMVA_significance_xx;}
   double             MET_CorrMVA_significance_xy(){return Ntp->MET_CorrMVA_significance_xy;}
   double             MET_CorrMVA_significance_yy(){return Ntp->MET_CorrMVA_significance_yy;}
   unsigned int   	  NMET_CorrMVA_srcMuons(){return Ntp->MET_CorrMVA_srcMuon_p4->size();}
   TLorentzVector     MET_CorrMVA_srcMuon_p4(unsigned int i){return TLorentzVector(Ntp->MET_CorrMVA_srcMuon_p4->at(i).at(1),Ntp->MET_CorrMVA_srcMuon_p4->at(i).at(2),Ntp->MET_CorrMVA_srcMuon_p4->at(i).at(3),Ntp->MET_CorrMVA_srcMuon_p4->at(i).at(0));}
   bool				  findCorrMVASrcMuon(unsigned int muon_idx, int &mvaSrcMuon_idx, float &dR );
   unsigned int   	  NMET_CorrMVA_srcElectrons(){return Ntp->MET_CorrMVA_srcElectron_p4->size();}
   TLorentzVector     MET_CorrMVA_srcElectron_p4(unsigned int i){return TLorentzVector(Ntp->MET_CorrMVA_srcElectron_p4->at(i).at(1),Ntp->MET_CorrMVA_srcElectron_p4->at(i).at(2),Ntp->MET_CorrMVA_srcElectron_p4->at(i).at(3),Ntp->MET_CorrMVA_srcElectron_p4->at(i).at(0));}
   bool				  findCorrMVASrcElectron(unsigned int elec_idx, int &mvaSrcElectron_idx, float &dR );
   unsigned int   	  NMET_CorrMVA_srcTaus(){return Ntp->MET_CorrMVA_srcTau_p4->size();}
   TLorentzVector     MET_CorrMVA_srcTau_p4(unsigned int i){return TLorentzVector(Ntp->MET_CorrMVA_srcTau_p4->at(i).at(1),Ntp->MET_CorrMVA_srcTau_p4->at(i).at(2),Ntp->MET_CorrMVA_srcTau_p4->at(i).at(3),Ntp->MET_CorrMVA_srcTau_p4->at(i).at(0));}
   bool				  findCorrMVASrcTau(unsigned int tau_idx, int &mvaSrcTau_idx, float &dR );

   double             MET_CorrMVAMuTau_et(){return Ntp->MET_CorrMVAMuTau_et;}
   double             MET_CorrMVAMuTau_phi(){return Ntp->MET_CorrMVAMuTau_phi;}
   double             MET_CorrMVAMuTau_ex(){return Ntp->MET_CorrMVAMuTau_et*cos(Ntp->MET_CorrMVAMuTau_phi);}
   double             MET_CorrMVAMuTau_ey(){return Ntp->MET_CorrMVAMuTau_et*sin(Ntp->MET_CorrMVAMuTau_phi);}
   double             MET_CorrMVAMuTau_significance(){return Ntp->MET_CorrMVAMuTau_significance;}
   double             MET_CorrMVAMuTau_significance_xx(){return Ntp->MET_CorrMVAMuTau_significance_xx;}
   double             MET_CorrMVAMuTau_significance_xy(){return Ntp->MET_CorrMVAMuTau_significance_xy;}
   double             MET_CorrMVAMuTau_significance_yy(){return Ntp->MET_CorrMVAMuTau_significance_yy;}
   unsigned int   	  NMET_CorrMVAMuTau_srcMuons(){return Ntp->MET_CorrMVAMuTau_srcMuon_p4->size();}
   TLorentzVector     MET_CorrMVAMuTau_srcMuon_p4(unsigned int i){return TLorentzVector(Ntp->MET_CorrMVAMuTau_srcMuon_p4->at(i).at(1),Ntp->MET_CorrMVAMuTau_srcMuon_p4->at(i).at(2),Ntp->MET_CorrMVAMuTau_srcMuon_p4->at(i).at(3),Ntp->MET_CorrMVAMuTau_srcMuon_p4->at(i).at(0));}
   bool				  findCorrMVAMuTauSrcMuon(unsigned int muon_idx, int &mvaMuTauSrcMuon_idx, float &dR );
   unsigned int   	  NMET_CorrMVAMuTau_srcTaus(){return Ntp->MET_CorrMVAMuTau_srcTau_p4->size();}
   TLorentzVector     MET_CorrMVAMuTau_srcTau_p4(unsigned int i){return TLorentzVector(Ntp->MET_CorrMVAMuTau_srcTau_p4->at(i).at(1),Ntp->MET_CorrMVAMuTau_srcTau_p4->at(i).at(2),Ntp->MET_CorrMVAMuTau_srcTau_p4->at(i).at(3),Ntp->MET_CorrMVAMuTau_srcTau_p4->at(i).at(0));}
   bool				  findCorrMVAMuTauSrcTau(unsigned int tau_idx, int &mvaMuTauSrcTau_idx, float &dR );

   ///////////////////////////////
   //
   // MET uncertainties
   //
   // Type1      = T0 + T1 + Txy (recommendation by JetMET POG and Christian Veelken)
   // Type1Type2 = T0 + T1 + Txy + calibration for unclustered energy (better MET response, worse MET resolution)
   //
   double             MET_Type1CorrElectronUp_et(){ return Ntp->MET_Type1CorrElectronUp_et; }
   double             MET_Type1CorrElectronDown_et(){ return Ntp->MET_Type1CorrElectronDown_et; }
   double             MET_Type1CorrMuonUp_et(){ return Ntp->MET_Type1CorrMuonUp_et; }
   double             MET_Type1CorrMuonDown_et(){ return Ntp->MET_Type1CorrMuonDown_et; }
   double             MET_Type1CorrTauUp_et(){ return Ntp->MET_Type1CorrTauUp_et; }
   double             MET_Type1CorrTauDown_et(){ return Ntp->MET_Type1CorrTauDown_et; }
   double             MET_Type1CorrJetResUp_et(){ return Ntp->MET_Type1CorrJetResUp_et; }
   double             MET_Type1CorrJetResDown_et(){ return Ntp->MET_Type1CorrJetResDown_et; }
   double             MET_Type1CorrJetEnUp_et(){ return Ntp->MET_Type1CorrJetEnUp_et; }
   double             MET_Type1CorrJetEnDown_et(){ return Ntp->MET_Type1CorrJetEnDown_et; }
   double             MET_Type1CorrUnClusteredUp_et(){ return Ntp->MET_Type1CorrUnclusteredUp_et; }
   double             MET_Type1CorrUnClusteredDown_et(){ return Ntp->MET_Type1CorrUnclusteredDown_et; }
   double             MET_Type1p2CorrElectronUp_et(){ return Ntp->MET_Type1p2CorrElectronUp_et; }
   double             MET_Type1p2CorrElectronDown_et(){ return Ntp->MET_Type1p2CorrElectronDown_et; }
   double             MET_Type1p2CorrMuonUp_et(){ return Ntp->MET_Type1p2CorrMuonUp_et; }
   double             MET_Type1p2CorrMuonDown_et(){ return Ntp->MET_Type1p2CorrMuonDown_et; }
   double             MET_Type1p2CorrTauUp_et(){ return Ntp->MET_Type1p2CorrTauUp_et; }
   double             MET_Type1p2CorrTauDown_et(){ return Ntp->MET_Type1p2CorrTauDown_et; }
   double             MET_Type1p2CorrJetResUp_et(){ return Ntp->MET_Type1p2CorrJetResUp_et; }
   double             MET_Type1p2CorrJetResDown_et(){ return Ntp->MET_Type1p2CorrJetResDown_et; }
   double             MET_Type1p2CorrJetEnUp_et(){ return Ntp->MET_Type1p2CorrJetEnUp_et; }
   double             MET_Type1p2CorrJetEnDown_et(){ return Ntp->MET_Type1p2CorrJetEnDown_et; }
   double             MET_Type1p2CorrUnclusteredUp_et(){ return Ntp->MET_Type1p2CorrUnclusteredUp_et; }
   double             MET_Type1p2CorrUnclusteredDown_et(){ return Ntp->MET_Type1p2CorrUnclusteredDown_et; }

   //Track Information
   unsigned int      NTracks(){return Ntp->Track_p4->size();}
   TLorentzVector    Track_p4(unsigned int i){return TLorentzVector(Ntp->Track_p4->at(i).at(1),Ntp->Track_p4->at(i).at(2),Ntp->Track_p4->at(i).at(3),Ntp->Track_p4->at(i).at(0));}
   TVector3          Track_Poca(unsigned int i){return TVector3(Ntp->Track_Poca->at(i).at(0),Ntp->Track_Poca->at(i).at(1),Ntp->Track_Poca->at(i).at(2));}
   int               Track_charge(unsigned int i){return Ntp->Track_charge->at(i);}
   double            Track_chi2(unsigned int i){return Ntp->Track_chi2->at(i);}
   double            Track_ndof(unsigned int i){return Ntp->Track_ndof->at(i);}
   unsigned short    Track_numberOfLostHits(unsigned int i){return Ntp->Track_numberOfLostHits->at(i);}
   unsigned short    Track_numberOfValidHits(unsigned int i){return Ntp->Track_numberOfValidHits->at(i);}
   unsigned int      Track_qualityMask(unsigned int i){return Ntp->Track_qualityMask->at(i);}

   TrackParticle Track_TrackParticle(unsigned int i){
     TMatrixT<double>    track_par(TrackParticle::NHelixPar,1);
     TMatrixTSym<double> track_cov(TrackParticle::NHelixPar);
     unsigned int l=0;
     for(int k=0; k<TrackParticle::NHelixPar; k++){
       track_par(k,0)=Ntp->Track_par->at(i).at(k);
       for(int j=k; j<TrackParticle::NHelixPar; j++){
	 track_cov(k,j)=Ntp->Track_cov->at(i).at(l);
	 l++;
       }
     }
     return TrackParticle(track_par,track_cov,Ntp->Track_pdgid->at(i),Ntp->Track_M->at(i),Ntp->Track_charge->at(i),Ntp->Track_B->at(i));
   }

   // MC Information
   // Signal particles (Z0,W+/-,H0,H+/-)
   unsigned int               NMCSignalParticles(){return Ntp->MCSignalParticle_p4->size();}
   TLorentzVector             MCSignalParticle_p4(unsigned int i){return TLorentzVector(Ntp->MCSignalParticle_p4->at(i).at(1),Ntp->MCSignalParticle_p4->at(i).at(2),Ntp->MCSignalParticle_p4->at(i).at(3),Ntp->MCSignalParticle_p4->at(i).at(0));}
   int                        MCSignalParticle_pdgid(unsigned int i){return Ntp->MCSignalParticle_pdgid->at(i);}
   int                        MCSignalParticle_charge(unsigned int i){return Ntp->MCSignalParticle_charge->at(i);}
   TVector3                   MCSignalParticle_Poca(unsigned int i){return TVector3(Ntp->MCSignalParticle_Poca->at(i).at(0),Ntp->MCSignalParticle_Poca->at(i).at(1),Ntp->MCSignalParticle_Poca->at(i).at(2));}
   std::vector<unsigned int>  MCSignalParticle_Tauidx(unsigned int i){return Ntp->MCSignalParticle_Tauidx->at(i);}
   // full MC chain
   unsigned int               NMCParticles(){return Ntp->MC_p4->size();}
   TLorentzVector             MCParticle_p4(unsigned int i){return TLorentzVector(Ntp->MC_p4->at(i).at(1),Ntp->MC_p4->at(i).at(2),Ntp->MC_p4->at(i).at(3),Ntp->MC_p4->at(i).at(0));}
   int                        MCParticle_pdgid(unsigned int i){return Ntp->MC_pdgid->at(i);}
   int                        MCParticle_charge(unsigned int i){return Ntp->MC_charge->at(i);}
   int              		  MCParticle_midx(unsigned int i){return Ntp->MC_midx->at(i);}
   std::vector<int>           MCParticle_childpdgid(unsigned int i){return Ntp->MC_childpdgid->at(i);}
   std::vector<int>           MCParticle_childidx(unsigned int i){return Ntp->MC_childidx->at(i);}
   int						  MCParticle_status(unsigned int i){return Ntp->MC_status->at(i);}
   int 						  getMatchTruthIndex(TLorentzVector tvector, int pid, double dr);
   int						  matchTruth(TLorentzVector tvector);
   bool						  matchTruth(TLorentzVector tvector, int pid, double dr);
   // decay tree functionality
   bool						  MCParticle_hasMother(unsigned int i){return Ntp->MC_midx->at(i) >= 0;}
   void						  printMCDecayChainOfMother(unsigned int i, bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false); // decay chain of object i
   void						  printMCDecayChainOfEvent(bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false); // full event decay chain
   std::string				  MCParticleToString(unsigned int par, bool printStatus = false, bool printPt = false, bool printEtaPhi = false);


   // Tau decays (Tau is first element of vector)
   int NMCTaus(){return Ntp->MCTauandProd_p4->size();}
   TLorentzVector MCTau_p4(unsigned int i){return MCTauandProd_p4(i,0);}
   int MCTau_pdgid(unsigned int i){return MCTauandProd_pdgid(i,0);}
   int MCTau_charge(unsigned int i){return MCTauandProd_charge(i,0);}
   unsigned int MCTau_JAK(unsigned int i){return Ntp->MCTau_JAK->at(i);}
   unsigned int MCTau_DecayBitMask(unsigned int i){return Ntp->MCTau_DecayBitMask->at(i);}
   int MCTau_getDaughterOfType(unsigned int i_mcTau, int daughter_pdgid, bool ignoreCharge = true);
   int MCTau_true3prongAmbiguity(unsigned int i);
   int matchTauTruth(unsigned int i_hpsTau, bool onlyHadrDecays = false);
   TLorentzVector BoostToRestFrame(TLorentzVector TLV1, TLorentzVector TLV2);
   TLorentzVector MCTau_invisiblePart(unsigned int i);
   TLorentzVector MCTau_visiblePart(unsigned int i);

   //Tau and decay products
   int NMCTauDecayProducts(unsigned int i){if(0<=i && i<(unsigned int)NMCTaus()) return Ntp->MCTauandProd_p4->at(i).size(); return 0;}
   TLorentzVector MCTauandProd_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->MCTauandProd_p4->at(i).at(j).at(1),Ntp->MCTauandProd_p4->at(i).at(j).at(2),Ntp->MCTauandProd_p4->at(i).at(j).at(3),Ntp->MCTauandProd_p4->at(i).at(j).at(0));}
   int MCTauandProd_pdgid(unsigned int i, unsigned int j){return Ntp->MCTauandProd_pdgid->at(i).at(j);}
   unsigned int MCTauandProd_midx(unsigned int i, unsigned int j){return Ntp->MCTauandProd_midx->at(i).at(j);}
   int MCTauandProd_charge(unsigned int i, unsigned int j){return Ntp->MCTauandProd_charge->at(i).at(j);}
   TVector3 MCTauandProd_Vertex(unsigned int i, unsigned int j){
     return TVector3(Ntp->MCTauandProd_Vertex->at(i).at(j).at(0),Ntp->MCTauandProd_Vertex->at(i).at(j).at(1),Ntp->MCTauandProd_Vertex->at(i).at(j).at(2));
   }
   bool hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &idx);
   bool hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,unsigned int &tau1_idx, unsigned int &tau2_idx);


   // overlap of objects
   bool jethasMuonOverlap(unsigned int jet_idx,unsigned int &muon_idx);
   bool muonhasJetOverlap(unsigned int muon_idx,unsigned int &jet_idx);
   bool muonhasJetMatch(unsigned int muon_idx,unsigned int &jet_idx);


   // Electrons
   unsigned int       NElectrons(){return Ntp->Electron_p4->size();}
   TLorentzVector	  Electron_p4(unsigned int i, TString corr = "default");
   TVector3           Electron_Poca(unsigned int i){return TVector3(Ntp->Electron_Poca->at(i).at(0),Ntp->Electron_Poca->at(i).at(1),Ntp->Electron_Poca->at(i).at(2));}
   int   Electron_Charge(unsigned int i){return Ntp->Electron_charge->at(i);}
   float   Electron_Gsf_deltaEtaEleClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaEtaEleClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaEtaSeedClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaEtaSeedClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(unsigned int i){return Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx->at(i);}
   float   Electron_Gsf_deltaPhiEleClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaPhiEleClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaPhiSeedClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaPhiSeedClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(unsigned int i){return Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx->at(i);}
   float   Electron_Gsf_dr03EcalRecHitSumE(unsigned int i){return Ntp->Electron_Gsf_dr03EcalRecHitSumE->at(i);}
   float   Electron_Gsf_dr03HcalDepth1TowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth1TowerSumEt->at(i);}
   float   Electron_Gsf_dr03HcalDepth1TowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth1TowerSumEtBc->at(i);}
   float   Electron_Gsf_dr03HcalDepth2TowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth2TowerSumEt->at(i);}
   float   Electron_Gsf_dr03HcalDepth2TowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth2TowerSumEtBc->at(i);}
   float   Electron_Gsf_dr03HcalTowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalTowerSumEt->at(i);}
   float   Electron_Gsf_dr03HcalTowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalTowerSumEtBc->at(i);}
   float   Electron_Gsf_dr03TkSumPt(unsigned int i){return Ntp->Electron_Gsf_dr03TkSumPt->at(i);}
   bool    Electron_Gsf_passingCutBasedPreselection(unsigned int i){return Ntp->Electron_Gsf_passingCutBasedPreselection->at(i);}
   bool    Electron_Gsf_passingMvaPreselection(unsigned int i){return Ntp->Electron_Gsf_passingMvaPreselection->at(i);}
   int     Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits(unsigned int i){return Ntp->Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits->at(i);}
   double  Electron_supercluster_e(unsigned int i){return Ntp->Electron_supercluster_e->at(i);}
   double  Electron_supercluster_phi(unsigned int i){return Ntp->Electron_supercluster_phi->at(i);}
   double  Electron_supercluster_eta(unsigned int i){return Ntp->Electron_supercluster_eta->at(i);}
   float   Electron_supercluster_centroid_x(unsigned int i){return Ntp->Electron_supercluster_centroid_x->at(i);}
   float   Electron_supercluster_centroid_y(unsigned int i){return Ntp->Electron_supercluster_centroid_y->at(i);}
   float   Electron_supercluster_centroid_z(unsigned int i){return Ntp->Electron_supercluster_centroid_z->at(i);}
   unsigned int Electron_Track_idx(unsigned int i){return Ntp->Electron_Track_idx->at(i);}

   float    Electron_ecalRecHitSumEt03(unsigned int i){return Ntp->Electron_ecalRecHitSumEt03->at(i);}
   float    Electron_hcalDepth1TowerSumEt03(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEt03->at(i);}
   float    Electron_hcalDepth1TowerSumEtBc03(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEtBc03->at(i);}
   float    Electron_hcalDepth2TowerSumEt03(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEt03->at(i);}
   float    Electron_hcalDepth2TowerSumEtBc03(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEtBc03->at(i);}
   float    Electron_tkSumPt03(unsigned int i){return Ntp->Electron_tkSumPt03->at(i);}
   float    Electron_ecalRecHitSumEt04(unsigned int i){return Ntp->Electron_ecalRecHitSumEt04->at(i);}
   float    Electron_hcalDepth1TowerSumEt04(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEt04->at(i);}
   float    Electron_hcalDepth1TowerSumEtBc04(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEtBc04->at(i);}
   float    Electron_hcalDepth2TowerSumEt04(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEt04->at(i);}
   float    Electron_hcalDepth2TowerSumEtBc04(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEtBc04->at(i);}
   float    Electron_tkSumPt04(unsigned int i){return Ntp->Electron_tkSumPt04->at(i);}
   
   float    Electron_chargedHadronIso(unsigned int i){return Ntp->Electron_chargedHadronIso->at(i);}
   float    Electron_neutralHadronIso(unsigned int i){return Ntp->Electron_neutralHadronIso->at(i);}
   float    Electron_photonIso(unsigned int i){return Ntp->Electron_photonIso->at(i);}
   
   double   Electron_isoDeposits_chargedHadronIso04(unsigned int i){return Ntp->Electron_isoDeposits_chargedHadronIso04->at(i);}
   double   Electron_isoDeposits_neutralHadronIso04(unsigned int i){return Ntp->Electron_isoDeposits_neutralHadronIso04->at(i);}
   double   Electron_isoDeposits_photonIso04(unsigned int i){return Ntp->Electron_isoDeposits_photonIso04->at(i);}
   double   Electron_isoDeposits_chargedHadronIso03(unsigned int i){return Ntp->Electron_isoDeposits_chargedHadronIso03->at(i);}
   double   Electron_isoDeposits_neutralHadronIso03(unsigned int i){return Ntp->Electron_isoDeposits_neutralHadronIso03->at(i);}
   double   Electron_isoDeposits_photonIso03(unsigned int i){return Ntp->Electron_isoDeposits_photonIso03->at(i);}

   float    Electron_sigmaIetaIeta(unsigned int i){return Ntp->Electron_sigmaIetaIeta->at(i);}
   float    Electron_hadronicOverEm(unsigned int i){return Ntp->Electron_hadronicOverEm->at(i);}
   float    Electron_fbrem(unsigned int i){return Ntp->Electron_fbrem->at(i);}
   float    Electron_eSuperClusterOverP(unsigned int i){return Ntp->Electron_eSuperClusterOverP->at(i);}
   float    Electron_ecalEnergy(unsigned int i){return Ntp->Electron_ecalEnergy->at(i);}
   float    Electron_trackMomentumAtVtx(unsigned int i){return Ntp->Electron_trackMomentumAtVtx->at(i);}
   int      Electron_numberOfMissedHits(unsigned int i){return Ntp->Electron_numberOfMissedHits->at(i);}
   bool     Electron_HasMatchedConversions(unsigned int i){return Ntp->Electron_HasMatchedConversions->at(i);}

   double   Electron_MVA_Trig_discriminator(unsigned int i){return Ntp->Electron_MVA_Trig_discriminator->at(i);}
   double   Electron_MVA_TrigNoIP_discriminator(unsigned int i){return Ntp->Electron_MVA_TrigNoIP_discriminator->at(i);}
   double   Electron_MVA_NonTrig_discriminator(unsigned int i){return Ntp->Electron_MVA_NonTrig_discriminator->at(i);}
   double   RhoIsolationAllInputTags(){return Ntp->RhoIsolationAllInputTags;}

   double   Electron_RegEnergy(unsigned int i){return Ntp->Electron_RegEnergy->at(i);}
   double   Electron_RegEnergyError(unsigned int i){return Ntp->Electron_RegEnergyError->at(i);}

   TrackParticle Electron_TrackParticle(unsigned int i){
     TMatrixT<double>    e_par(TrackParticle::NHelixPar,1);
     TMatrixTSym<double> e_cov(TrackParticle::NHelixPar);
     unsigned int l=0;
     for(int k=0; k<TrackParticle::NHelixPar; k++){
       e_par(k,0)=Ntp->Electron_par->at(i).at(k);
       for(int j=k; j<TrackParticle::NHelixPar; j++){
	 e_cov(k,j)=Ntp->Electron_cov->at(i).at(l);
	 l++;
       }
     }
     return TrackParticle(e_par,e_cov,Ntp->Electron_pdgid->at(i),Ntp->Electron_M->at(i),Ntp->Electron_charge->at(i),Ntp->Electron_B->at(i));
   }

   bool isTrigPreselElectron(unsigned int i);
   bool isTrigNoIPPreselElectron(unsigned int i);
   bool isMVATrigElectron(unsigned int i, TString corr = "default");
   bool isMVATrigNoIPElectron(unsigned int i, TString corr = "default");
   bool isMVANonTrigElectron(unsigned int i, unsigned int j, TString corr = "default");
   bool isTightElectron(unsigned int i, TString corr = "default");
   bool isTightElectron(unsigned int i, unsigned int j, TString corr = "default");
   float Electron_RelIso03(unsigned int i, TString corr = "default");
   double Electron_RelIsoDep03(unsigned int i, TString corr = "default");//todo: better name
   float Electron_RelIso04(unsigned int i, TString corr = "default");
   double Electron_RelIsoDep04(unsigned int i, TString corr = "default");//todo: better name
   double Electron_Aeff_R04(double Eta);
   double Electron_Aeff_R03(double Eta);
   bool isSelectedElectron(unsigned int i, unsigned int j, double impact_xy, double impact_z, TString corr = "default");

   // Trigger
   bool         TriggerAccept(TString n);
   unsigned int HLTPrescale(TString n);
   unsigned int L1SEEDPrescale(TString n);
   bool         GetTriggerIndex(TString n, unsigned int &i);
   double 		matchTrigger(TLorentzVector obj, std::vector<TString> trigger, std::string objectType);
   bool 		matchTrigger(TLorentzVector obj, double dr_cut, std::vector<TString> trigger, std::string objectType);
   bool			matchTrigger(TLorentzVector obj, double dr_cut, TString trigger, std::string objectType);
   unsigned int NHLTTriggers(){return Ntp->HTLTriggerName->size();}
   std::string  HTLTriggerName(unsigned int i){return Ntp->HTLTriggerName->at(i);}
   bool         TriggerAccept(unsigned int i){return Ntp->TriggerAccept->at(i);}
   bool         TriggerError(unsigned int i){return Ntp->TriggerError->at(i);}
   bool         TriggerWasRun(unsigned int i){return Ntp->TriggerWasRun->at(i);}
   unsigned int HLTPrescale(unsigned int i){return Ntp->HLTPrescale->at(i);}
   unsigned int NHLTL1GTSeeds(unsigned int i){return Ntp->NHLTL1GTSeeds->at(i);}
   unsigned int L1SEEDPrescale(unsigned int i){return Ntp->L1SEEDPrescale->at(i);}
   bool         L1SEEDInvalidPrescale(unsigned int i){return Ntp->L1SEEDInvalidPrescale->at(i);}
   unsigned int NHLTTriggerObject(unsigned int i){return Ntp->HLTTrigger_objs_Eta->at(i).size();}
   TLorentzVector HLTTriggerObject_p4(unsigned int i, unsigned int j){
     TLorentzVector L(0,0,0,0); 
     if(j<Ntp->HLTTrigger_objs_Eta->at(i).size())L.SetPtEtaPhiM(Ntp->HLTTrigger_objs_Pt->at(i).at(j),Ntp->HLTTrigger_objs_Eta->at(i).at(j), Ntp->HLTTrigger_objs_Phi->at(i).at(j),0.0);
     return L;
   }
   int          NHLTTrigger_objs(){return Ntp->HLTTrigger_objs_Pt->size();}
   int          NHLTTrigger_objs(unsigned int i){return Ntp->HLTTrigger_objs_Pt->at(i).size();}
   float        HLTTrigger_objs_Pt(unsigned int i, unsigned int j){return Ntp->HLTTrigger_objs_Pt->at(i).at(j);}
   float        HLTTrigger_objs_Eta(unsigned int i, unsigned int j){return Ntp->HLTTrigger_objs_Eta->at(i).at(j);}
   float        HLTTrigger_objs_Phi(unsigned int i, unsigned int j){return Ntp->HLTTrigger_objs_Phi->at(i).at(j);}
   float        HLTTrigger_objs_E(unsigned int i,unsigned int j){return Ntp->HLTTrigger_objs_E->at(i).at(j);}
   int          HLTTrigger_objs_Id(unsigned int i,unsigned int j){return Ntp->HLTTrigger_objs_Id->at(i).at(j);}
   std::string  HLTTrigger_objs_trigger(unsigned int i){return Ntp->HLTTrigger_objs_trigger->at(i);}

   // helper functions
   double       dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
   double		dxySigned(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
   double       dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
   double       dzSigned(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
   double       vertexSignificance(TVector3 vec, unsigned int vertex);
   double		transverseMass(double pt1, double phi1, double pt2, double phi2){return sqrt(2 * pt1 * pt2 * (1 - cos(phi1 - phi2)));}
   std::vector<int> sortObjects(std::vector<int> indices, std::vector<double> values);
   std::vector<int> sortPFJetsByPt();
   std::vector<int> sortDefaultObjectsByPt(TString objectType);

   int          GenEventInfoProduct_id1(){return Ntp->GenEventInfoProduct_id1;}
   int          GenEventInfoProduct_id2(){return Ntp->GenEventInfoProduct_id2;}
   double       GenEventInfoProduct_x1(){return Ntp->GenEventInfoProduct_x1;}
   double       GenEventInfoProduct_x2(){return Ntp->GenEventInfoProduct_x2;}
   double       GenEventInfoProduct_scalePDF(){return Ntp->GenEventInfoProduct_scalePDF;}
};

#endif

