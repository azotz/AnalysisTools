#include "ZToTaumuTauh.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SkimConfig.h"
#include "PDG_Var.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"

ZToTaumuTauh::ZToTaumuTauh(TString Name_, TString id_):
  Selection(Name_,id_),
  cMu_dxy(0.045),// 100 is dummy value
  cMu_dz(0.2),
  cMu_relIso(0.1),
  cMu_pt(20.),
  cMu_eta(2.1),
  cMu_dRHltMatch(0.5),
  cTau_pt(20.),//20: Recommended by Tau POG
  cTau_eta(2.3),
  cMuTau_dR(0.3),
  cTau_IsoRaw(1.5),
  cTau_dRHltMatch(0.5)
{
	TString trigNames[] = {"HLT_IsoMu18_eta2p1_LooseIsoPFTau20","HLT_IsoMu17_eta2p1_LooseIsoPFTau20"};
	std::vector<TString> temp (trigNames, trigNames + sizeof(trigNames) / sizeof(TString) );
	cTriggerNames = temp;

	//https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013#TauES_and_decay_mode_scale_facto
	OneProngNoPiWeight = 0.88;

	//dummy values

	TriggerOkDummy = -1;
	selVertexDummy = -1;
	selMuon_IsoDummy = -1;
	selMuon_AntiIsoDummy = -1;
	selTauDummy = -1;
	ChargeSumDummy = -999;
	MTDummy = -999;
	MvisDummy = -999;
	TauFLSigmaDummy = -999;

	//WJets BG Method
	SB_lowerLimit = 70; //in GeV, Region > SB_lowerLimit dominated by W+Jets
	SB_upperLimit = 140; //in GeV, Region < SB_upperLimit
	Scaleby_Counting = true; // = false --> Scale by Integral

	tau_corr = "";

	//Set verbose boolean
	verbose = false;
	Use_Embedded = false;
}

ZToTaumuTauh::~ZToTaumuTauh(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZToTaumuTauh::~ZToTaumuTauh Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZToTaumuTauh::~ZToTaumuTauh()" << std::endl;
}

void  ZToTaumuTauh::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)    cut.at(TriggerOk)=0;
    if(i==PrimeVtx)     cut.at(PrimeVtx)=1;
    if(i==NMuId)		cut.at(NMuId)=1;
    if(i==NMuKin)		cut.at(NMuKin)=1;
    if(i==NMuIso)		cut.at(NMuIso)=1;
    if(i==NTauId)		cut.at(NTauId)=1;
    if(i==NTauKin)		cut.at(NTauKin)=1;
    if(i==NTauIso)		cut.at(NTauIso)=1;
    if(i==ChargeSum)	cut.at(ChargeSum)=0;
    if(i==MT_MuMET)		cut.at(MT_MuMET)=30;
    if(i==TauDecayMode)	cut.at(TauDecayMode)=10;//10
    if(i==TauFLSigma)	cut.at(TauFLSigma)=3;//3
  }

  TString hlabel;
  TString htitle;
  for(unsigned int i_cut=0; i_cut<NCuts; i_cut++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i_cut;
  
    if(i_cut==PrimeVtx){
      title.at(i_cut)="Number of Prime Vertices $(N>$";
      title.at(i_cut)+=cut.at(PrimeVtx);
      title.at(i_cut)+=")";
      htitle=title.at(i_cut);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,41,-0.5,40.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,41,-0.5,40.5,hlabel,"Events"));
    }
    else if(i_cut==TriggerOk){
      title.at(i_cut)="Trigger ";
      hlabel="Trigger ";

      std::vector<TH1D> Nm1Temp = HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,cTriggerNames.size()+2,-1.5,cTriggerNames.size()+0.5,hlabel,"Events");
      std::vector<TH1D> Nm0Temp = HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,cTriggerNames.size()+2,-1.5,cTriggerNames.size()+0.5,hlabel,"Events");
      for (unsigned i_hist = 0; i_hist < Nm1Temp.size(); i_hist++){
    	  Nm1Temp.at(i_hist).GetXaxis()->SetBinLabel(1,"not fired");
    	  Nm0Temp.at(i_hist).GetXaxis()->SetBinLabel(1,"not fired");
    	  for (unsigned i_bin = 2; i_bin < cTriggerNames.size()+2; i_bin++){
    		  Nm1Temp.at(i_hist).GetXaxis()->SetBinLabel(i_bin,cTriggerNames.at(i_bin-2));
    		  Nm0Temp.at(i_hist).GetXaxis()->SetBinLabel(i_bin,cTriggerNames.at(i_bin-2));
    	  }
    	  Nm1Temp.at(i_hist).GetXaxis()->SetBinLabel(cTriggerNames.size()+2,"multiple fired");
    	  Nm0Temp.at(i_hist).GetXaxis()->SetBinLabel(cTriggerNames.size()+2,"multiple fired");
      }
      Nminus1.push_back(Nm1Temp);
      Nminus0.push_back(Nm0Temp);
    }
    else if(i_cut==NMuId){
    	title.at(i_cut)="Number $\\mu_{ID} >=$";
    	title.at(i_cut)+=cut.at(NMuId);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{ID}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==NMuKin){
    	title.at(i_cut)="Number $\\mu_{kin} >=$";
    	title.at(i_cut)+=cut.at(NMuKin);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{Kin}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuKin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuKin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i_cut==NMuIso){
    	title.at(i_cut)="Number $\\mu_{Iso} >=$";
    	title.at(i_cut)+=cut.at(NMuIso);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{Iso}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuIso_",htitle,6,-0.5,5.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuIso_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i_cut==NTauId){
    	title.at(i_cut)="Number $\\tau_{ID} >=$";
    	title.at(i_cut)+=cut.at(NTauId);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{ID}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==NTauKin){
    	title.at(i_cut)="Number $\\tau_{Kin} >=$";
    	title.at(i_cut)+=cut.at(NTauKin);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{Kin}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauKin_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauKin_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==NTauIso){
    	title.at(i_cut)="Number $\\tau_{Iso} >=$";
    	title.at(i_cut)+=cut.at(NTauIso);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{Iso}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauIso_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauIso_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==ChargeSum){
    	title.at(i_cut)="Sum of Charges = ";
    	title.at(i_cut)+=cut.at(ChargeSum);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Sum of Charges";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ChargeSum_",htitle,5,-2.5,2.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ChargeSum_",htitle,5,-2.5,2.5,hlabel,"Events"));
    }
    else if(i_cut==MT_MuMET){
    	title.at(i_cut)="$m_{T}(\\mu,MET) <$";
    	title.at(i_cut)+=cut.at(MT_MuMET);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="m_{T}(#mu,MET)";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MT_MuMET_",htitle,75,0,150.,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MT_MuMET_",htitle,75,0,150.,hlabel,"Events"));
    }
    else if(i_cut==TauDecayMode){
    	title.at(i_cut)="Tau Decay Mode = ";
    	title.at(i_cut)+=cut.at(TauDecayMode);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="TauDecayMode";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauDecayMode_",htitle,15,0,15,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauDecayMode_",htitle,15,0,15,hlabel,"Events"));
    }
    else if(i_cut==TauFLSigma){
    	title.at(i_cut)="TauFLSigma = ";
    	title.at(i_cut)+=cut.at(TauFLSigma);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="TauFLSigma";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauFLSigma_",htitle,80,-10,30,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauFLSigma_",htitle,80,-10,30,hlabel,"Events"));
    }
  }
  // Setup NPassed Histograms
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",41,-0.5,40.5,"Number of Vertices","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of good Vertices","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  NMtTauMET=HConfig.GetTH1D(Name+"_NMtTauMET","NMtTauMET",75,0,150,"m_{T}(#tau,MET)","Events");

  NMvis=HConfig.GetTH1D(Name+"_NMvis","NMvis",75,-0.5,149.5,"m_{vis}(#mu,#tau)","Events");

  NSB=HConfig.GetTH1D(Name+"_NSB","NSB",2,0.5,2.5,"SB Region (all Samples)","Events");

  Mu_pt=HConfig.GetTH1D(Name+"_Mu_pt","Mu_pt",50,0,100,"Muon p_{t}","Events");
  Mu_phi=HConfig.GetTH1D(Name+"_Mu_phi","Mu_phi",30,-3.14159265359,3.14159265359,"Muon #phi","Events");
  Mu_eta=HConfig.GetTH1D(Name+"_Mu_eta","Mu_eta",50,-2.5,2.5,"Muon #eta","Events");

  Tau_pt=HConfig.GetTH1D(Name+"_Tau_pt","Tau_pt",50,0,100,"Tau p_{t}","Events");
  Tau_phi=HConfig.GetTH1D(Name+"_Tau_phi","Tau_phi",30,-3.14159265359,3.14159265359,"Tau #phi","Events");
  Tau_eta=HConfig.GetTH1D(Name+"_Tau_eta","Tau_eta",50,-2.5,2.5,"Tau #eta","Events");

  MET_phi=HConfig.GetTH1D(Name+"_MET_phi","MET_phi",30,-3.14159265359,3.14159265359,"MET #phi","Events");

  TauFL=HConfig.GetTH1D(Name+"_TauFL","TauFL",100,-2,4,"tau flight length","Events");
  TauFLSigned=HConfig.GetTH1D(Name+"_TauFLSigned","TauFLSigned",100,-2,4,"tau flight length Signed","Events");

  TauFLSigmaSigned=HConfig.GetTH1D(Name+"_TauFLSigmaSigned","TauFLSigma Signed",80,-10,30,"TauFLSigma Signed","Events");
  TauFLSigmaUnsigned=HConfig.GetTH1D(Name+"_TauFLSigmaUnsigned","TauFLSigma Unsigned",80,-10,30,"TauFLSigma Unsigned","Events");

  A1mass=HConfig.GetTH1D(Name+"_A1mass","A1mass",100,0,2,"A1mass","Events");
  A1mass10GeV=HConfig.GetTH1D(Name+"_A1massGeV","A1massGeV",100,0,10,"A1mass","Events");

  Mvis3Prong=HConfig.GetTH1D(Name+"_Mvis3Prong","Mvis3Prong",75,-0.5,149.5,"m_{vis}(#mu,#tau) 3 prong","Events");
  Mvis1Prong=HConfig.GetTH1D(Name+"_Mvis1Prong","Mvis1Prong",75,-0.5,149.5,"m_{vis}(#mu,#tau) 1 prong","Events");
  MvisIncl=HConfig.GetTH1D(Name+"_MvisAll","MvisAll",75,-0.5,149.5,"m_{vis}(#mu,#tau) incl.","Events");

  MTMuMET3Prong=HConfig.GetTH1D(Name+"_MTMuMET3Prong","MTMuMET3Prong",75,0,150,"m_{T}(#mu,MET) 3 prong","Events");
  MTMuMET1Prong=HConfig.GetTH1D(Name+"_MTMuMET1Prong","MTMuMET1Prong",75,0,150,"m_{T}(#mu,MET) 1 prong","Events");
  MTMuMETIncl=HConfig.GetTH1D(Name+"_MTMuMETIncl","MTMuMETIncl",75,0,150,"m_{T}(#mu,MET) incl","Events");

  //Gen Studies
  Mvis_SignalOnly=HConfig.GetTH1D(Name+"_Mvis_SignalOnly","Mvis_SignalOnly",75,-0.5,149.5,"m_{vis}(#mu,#tau)","Events");
  Mvis_SignalOnly_genMu=HConfig.GetTH1D(Name+"_Mvis_SignalOnly_genMu","Mvis_SignalOnly_genMu",75,-0.5,149.5,"m_{vis}(#mu_{gen},#tau)","Events");
  Mvis_SignalOnly_genA1=HConfig.GetTH1D(Name+"_Mvis_SignalOnly_genA1","Mvis_SignalOnly_genA1",75,-0.5,149.5,"m_{vis}(#mu,a_{1,gen})","Events");
  Mvis_SignalOnly_genTaumu=HConfig.GetTH1D(Name+"_Mvis_SignalOnly_genTaumu","Mvis_SignalOnly_genTaumu",75,-0.5,149.5,"m_{vis}(#tau_{#mu,gen},a_{1})","Events");
  Mvis_SignalOnly_genTauh=HConfig.GetTH1D(Name+"_Mvis_SignalOnly_genTauh","Mvis_SignalOnly_genTauh",75,-0.5,149.5,"m_{vis}(#mu,#tau_{h,gen})","Events");

  dR_selTauh_genTauh=HConfig.GetTH1D(Name+"_dR_selTauh_genTauh","dR_selTauh_genTauh",1,0,50,"dR_selTauh_genTauh","Events");
  dR_selMu_genMu=HConfig.GetTH1D(Name+"_dR_selMu_genMu","dR_selMu_genMu",1,0,50,"dR_selMu_genMu","Events");

  POCAPV_Mag=HConfig.GetTH1D(Name+"_POCAPV_Mag","POCAPV_Mag",100,0,10./100,"POCAPV_Mag","Events");

  Phi_SVPV=HConfig.GetTH1D(Name+"_Phi_SVPV","Phi_SVPV",32,-3.14159265359,3.14159265359,"Phi_SVPV","Events");
  Phi_genTauh=HConfig.GetTH1D(Name+"_Phi_genTauh","Phi_genTauh",32,-3.14159265359,3.14159265359,"Phi_genTauh","Events");
  Theta_SVPV=HConfig.GetTH1D(Name+"_Theta_SVPV","Theta_SVPV",32,0.,3.14159265359,"Theta_SVPV","Events");
  Theta_genTauh=HConfig.GetTH1D(Name+"_Theta_genTauh","Theta_genTauh",32,0.,3.14159265359,"Theta_genTauh","Events");
  dPhi_SVPV_genTauh=HConfig.GetTH1D(Name+"_dPhi_SVPV_genTauh","dPhi_SVPV_genTauh",128,-3.14159265359/32,3.14159265359/32,"dPhi_SVPV_genTauh","Events");
  dTheta_SVPV_genTauh=HConfig.GetTH1D(Name+"_dTheta_SVPV_genTauh","dTheta_SVPV_genTauh",64,-3.14159265359/32,3.14159265359/32,"dTheta_SVPV_genTauh","Events");
  Angle_SVPV_genTauh=HConfig.GetTH1D(Name+"_Angle_SVPV_genTauh","Angle_SVPV_genTauh",64,0,3.14159265359/32,"Angle_SVPV_genTauh","Events");

  Phi_POCAPV=HConfig.GetTH1D(Name+"_Phi_POCAPV","Phi_POCAPV",32,-3.14159265359,3.14159265359,"Phi_POCAPV","Events");
  Phi_genTaumu=HConfig.GetTH1D(Name+"_Phi_genTaumu","Phi_genTaumu",32,-3.14159265359,3.14159265359,"Phi_genTaumu","Events");
  Theta_POCAPV=HConfig.GetTH1D(Name+"_Theta_POCAPV","Theta_POCAPV",32,0.,3.14159265359,"Theta_POCAPV","Events");
  Theta_genTaumu=HConfig.GetTH1D(Name+"_Theta_genTaumu","Theta_genTaumu",32,0.,3.14159265359,"Theta_genTaumu","Events");
  dPhi_POCAPV_genTaumu=HConfig.GetTH1D(Name+"_dPhi_POCAPV_genTaumu","dPhi_POCAPV_genTaumu",128,-3.14159265359,3.14159265359,"dPhi_POCAPV_genTaumu","Events");
  dTheta_POCAPV_genTaumu=HConfig.GetTH1D(Name+"_dTheta_POCAPV_genTaumu","dTheta_POCAPV_genTaumu",128,-3.14159265359,3.14159265359,"dTheta_POCAPV_genTaumu","Events");

  dPhi_MinusSVPV_genTaumu=HConfig.GetTH1D(Name+"_dPhi_MinusSVPV_genTaumu","dPhi_MinusSVPV_genTaumu",128,-3.14159265359/2,3.14159265359/2,"dPhi_MinusSVPV_genTaumu","Events");
  dTheta_MinusSVPV_genTaumu=HConfig.GetTH1D(Name+"_dTheta_MinusSVPV_genTaumu","dTheta_MinusSVPV_genTaumu",128,-3.14159265359/2,3.14159265359/2,"dTheta_MinusSVPV_genTaumu","Events");
  Angle_MinusSVPV_genTaumu=HConfig.GetTH1D(Name+"_Angle_MinusSVPV_genTaumu","Angle_MinusSVPV_genTaumu",128,0,3.14159265359,"Angle_MinusSVPV_genTaumu","Events");

  GJ_Tauh=HConfig.GetTH1D(Name+"_GJ_Tauh","GJ_Tauh",100,0,.05,"GJ_Tauh","Events");
  GJ_Taumu=HConfig.GetTH1D(Name+"_GJ_Taumu","GJ_Taumu",100,0,.05,"GJ_Taumu","Events");

  dPhi_DiTauGen=HConfig.GetTH1D(Name+"_dPhi_DiTauGen","dPhi_DiTauGen",256,0,2*3.14159265359,"dPhi_DiTauGen","Events");
  Pt_DiTauGen=HConfig.GetTH1D(Name+"_Pt_DiTauGen","Pt_DiTauGen",30,0,30,"Pt_DiTauGen","Events");
  Pt_ZGen=HConfig.GetTH1D(Name+"_Pt_ZGen","Pt_ZGen",30,0,30,"Pt_ZGen","Events");
  M_ZGen=HConfig.GetTH1D(Name+"_M_ZGen","M_ZGen",180,60,150,"M_ZGen","Events");
  M_DiTauPtBalance=HConfig.GetTH1D(Name+"_M_DiTauPtBalance","M_DiTauPtBalance",280,20,160,"M_DiTauPtBalance","Events");
  dM_DiTau=HConfig.GetTH1D(Name+"_dM_DiTau","dM_DiTau",100,-50,50,"dM_DiTau","Events");
  dPt_GenTaumuPtBalance=HConfig.GetTH1D(Name+"_dPt_GenTaumuPtBalance","dPt_GenTaumuPtBalance",100,-50,50,"dPt_GenTaumuPtBalance","Events");
  dP_GenTaumuPtBalance=HConfig.GetTH1D(Name+"_dP_GenTaumuPtBalance","dP_GenTaumuPtBalance",100,-50,50,"dP_GenTaumuPtBalance","Events");
  dP_GenTauh=HConfig.GetTH1D(Name+"_dP_GenTauh","dP_GenTauh",50,0,50,"dP_GenTauh","Events");

  dP_GenTauMuPtBalance_vs_dPTauh=HConfig.GetTH2D(Name+"_dP_GenTauMuPtBalance_vs_dPTauh","dP_GenTauMuPtBalance_vs_dPTauh",50,-50,50,50,0,50,"p^{'}_{#tau_{#mu}} - p_{#tau_{#mu}} ","2*D");
  Pt_vs_dPhi_DiTauGen=HConfig.GetTH2D(Name+"_Pt_vs_dPhi_DiTauGen","Pt_vs_dPhi_DiTauGen",10,0,30,256,0,2*3.14159265359,"Pt ditau","dPhi(tau,tau)");

  TauFLSigmaCut_vs_Res=HConfig.GetTH2D(Name+"_TauFLSigmaCut_vs_Res","TauFLSigmaCut_vs_Res",80,-10,30,2000,-0.1,0.1,"TauFLSigmaCut","dPhi(SVPV,genTauh)");
  TauFLSigma_vs_Res=HConfig.GetTH2D(Name+"_TauFLSigma_vs_Res","TauFLSigma_vs_Res",80,-10,30,200,-0.1,0.1,"TauFLSigma","dPhi(SVPV,genTauh)");

  //QCD Histos
  NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",4,0.5,4.5,"NQCD in ABCD","Events");

  QCD_MT_MuMET_A=HConfig.GetTH1D(Name+"_QCD_MT_MuMET_A","QCD_MT_MuMET_A",75,0,150.,"m_{T}(#mu,MET) in A","Events");
  QCD_MT_MuMET_B=HConfig.GetTH1D(Name+"_QCD_MT_MuMET_B","QCD_MT_MuMET_B",75,0,150.,"m_{T}(#mu,MET) in B","Events");
  QCD_MT_MuMET_C=HConfig.GetTH1D(Name+"_QCD_MT_MuMET_C","QCD_MT_MuMET_C",75,0,150.,"m_{T}(#mu,MET) in C","Events");
  QCD_MT_MuMET_D=HConfig.GetTH1D(Name+"_QCD_MT_MuMET_D","QCD_MT_MuMET_D",75,0,150.,"m_{T}(#mu,MET) in D","Events");

  QCD_MET_A=HConfig.GetTH1D(Name+"_QCD_MET_A","QCD_MET_A",75,0,150.,"MVA MET in A","Events");
  QCD_MET_B=HConfig.GetTH1D(Name+"_QCD_MET_B","QCD_MET_B",75,0,150.,"MVA MET in B","Events");
  QCD_MET_C=HConfig.GetTH1D(Name+"_QCD_MET_C","QCD_MET_C",75,0,150.,"MVA MET in C","Events");
  QCD_MET_D=HConfig.GetTH1D(Name+"_QCD_MET_D","QCD_MET_D",75,0,150.,"MVA MET in D","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}

void  ZToTaumuTauh::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist1d.push_back(&NMtTauMET);
 Extradist1d.push_back(&NMvis);

 Extradist1d.push_back(&Mvis_SignalOnly);
 Extradist1d.push_back(&Mvis_SignalOnly_genMu);
 Extradist1d.push_back(&Mvis_SignalOnly_genA1);
 Extradist1d.push_back(&Mvis_SignalOnly_genTaumu);
 Extradist1d.push_back(&Mvis_SignalOnly_genTauh);

 Extradist1d.push_back(&NSB);
 Extradist1d.push_back(&Mu_pt);
 Extradist1d.push_back(&Mu_phi);
 Extradist1d.push_back(&Mu_eta);
 Extradist1d.push_back(&Tau_pt);
 Extradist1d.push_back(&Tau_phi);
 Extradist1d.push_back(&Tau_eta);
 Extradist1d.push_back(&MET_phi);
 Extradist1d.push_back(&TauFL);
 Extradist1d.push_back(&TauFLSigned);
 Extradist1d.push_back(&TauFLSigmaSigned);
 Extradist1d.push_back(&TauFLSigmaUnsigned);
 Extradist1d.push_back(&A1mass);
 Extradist1d.push_back(&A1mass10GeV);

 Extradist1d.push_back(&Mvis3Prong);
 Extradist1d.push_back(&Mvis1Prong);
 Extradist1d.push_back(&MvisIncl);
 Extradist1d.push_back(&MTMuMET3Prong);
 Extradist1d.push_back(&MTMuMET1Prong);
 Extradist1d.push_back(&MTMuMETIncl);

 Extradist1d.push_back(&POCAPV_Mag);
 Extradist1d.push_back(&dR_selTauh_genTauh);
 Extradist1d.push_back(&dR_selMu_genMu);
 Extradist1d.push_back(&Phi_SVPV);
 Extradist1d.push_back(&Phi_genTauh);
 Extradist1d.push_back(&Theta_SVPV);
 Extradist1d.push_back(&Theta_genTauh);
 Extradist1d.push_back(&dPhi_SVPV_genTauh);
 Extradist1d.push_back(&dTheta_SVPV_genTauh);
 Extradist1d.push_back(&Angle_SVPV_genTauh);
 Extradist1d.push_back(&Phi_POCAPV);
 Extradist1d.push_back(&Phi_genTaumu);
 Extradist1d.push_back(&Theta_POCAPV);
 Extradist1d.push_back(&Theta_genTaumu);
 Extradist1d.push_back(&dPhi_POCAPV_genTaumu);
 Extradist1d.push_back(&dTheta_POCAPV_genTaumu);
 Extradist1d.push_back(&GJ_Tauh);
 Extradist1d.push_back(&GJ_Taumu);
 Extradist1d.push_back(&dPhi_DiTauGen);
 Extradist1d.push_back(&Pt_DiTauGen);
 Extradist1d.push_back(&Pt_ZGen);
 Extradist1d.push_back(&M_ZGen);
 Extradist1d.push_back(&M_DiTauPtBalance);
 Extradist1d.push_back(&dM_DiTau);
 Extradist1d.push_back(&dPt_GenTaumuPtBalance);
 Extradist1d.push_back(&dP_GenTaumuPtBalance);
 Extradist1d.push_back(&dP_GenTauh);

 Extradist1d.push_back(&dPhi_MinusSVPV_genTaumu);
 Extradist1d.push_back(&dTheta_MinusSVPV_genTaumu);
 Extradist1d.push_back(&Angle_MinusSVPV_genTaumu);

 Extradist2d.push_back(&dP_GenTauMuPtBalance_vs_dPTauh);
 Extradist2d.push_back(&Pt_vs_dPhi_DiTauGen);
 Extradist2d.push_back(&TauFLSigmaCut_vs_Res);
 Extradist2d.push_back(&TauFLSigma_vs_Res);

 Extradist1d.push_back(&NQCD);
 Extradist1d.push_back(&QCD_MT_MuMET_A);
 Extradist1d.push_back(&QCD_MT_MuMET_B);
 Extradist1d.push_back(&QCD_MT_MuMET_C);
 Extradist1d.push_back(&QCD_MT_MuMET_D);
 Extradist1d.push_back(&QCD_MET_A);
 Extradist1d.push_back(&QCD_MET_B);
 Extradist1d.push_back(&QCD_MET_C);
 Extradist1d.push_back(&QCD_MET_D);

 Extradist1d_OS.push_back(&NVtx);
 Extradist1d_OS.push_back(&NGoodVtx);
 Extradist1d_OS.push_back(&NTrackperVtx);
 Extradist1d_OS.push_back(&NMtTauMET);
 Extradist1d_OS.push_back(&NMvis);
 Extradist1d_OS.push_back(&Mu_pt);
 Extradist1d_OS.push_back(&Mu_phi);
 Extradist1d_OS.push_back(&Mu_eta);
 Extradist1d_OS.push_back(&Tau_pt);
 Extradist1d_OS.push_back(&Tau_phi);
 Extradist1d_OS.push_back(&Tau_eta);
 Extradist1d_OS.push_back(&MET_phi);
 Extradist1d_OS.push_back(&TauFL);
 Extradist1d_OS.push_back(&TauFLSigned);
 Extradist1d_OS.push_back(&TauFLSigmaSigned);
 Extradist1d_OS.push_back(&TauFLSigmaUnsigned);
 Extradist1d_OS.push_back(&A1mass);
 Extradist1d_OS.push_back(&A1mass10GeV);
 Extradist1d_OS.push_back(&Mvis3Prong);
 Extradist1d_OS.push_back(&Mvis1Prong);
 Extradist1d_OS.push_back(&MvisIncl);
 Extradist1d_OS.push_back(&MTMuMET3Prong);
 Extradist1d_OS.push_back(&MTMuMET1Prong);
 Extradist1d_OS.push_back(&MTMuMETIncl);
}

void  ZToTaumuTauh::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}

  int selVertex = selVertexDummy;
  int selMuon_Iso = selMuon_IsoDummy;
  int selMuon_AntiIso = selMuon_AntiIsoDummy;
  int selTau = selTauDummy;
  Charge = ChargeSumDummy;

  if(Ntp->GetMCID() == DataMCType::DY_tautau || (Ntp->GetMCID()>=10 && Ntp->GetMCID()<= 13)) tau_corr = "scalecorr";
  else tau_corr = "";

  // Apply Selection
  if(verbose) std::cout << "Cut on good vertex" << std::endl;
  unsigned int nGoodVtx=0;
  for(unsigned int i_vtx=0;i_vtx<Ntp->NVtx();i_vtx++){
    if(Ntp->isVtxGood(i_vtx)){
    	if(selVertex == selVertexDummy) selVertex = i_vtx; // selected vertex = first vertex (highest sum[pT^2]) to fulfill vertex requirements
    	nGoodVtx++;
    }
  }
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  if(verbose){
	  std::cout << "value at Primevtx: " <<value.at(PrimeVtx) << std::endl;
	  std::cout << "pass at Primevtx: " <<pass.at(PrimeVtx) << std::endl;
  }
  
  // Trigger
  if(verbose) std::cout << "Cut on Trigger" << std::endl;
  value.at(TriggerOk) = TriggerOkDummy;
  for (std::vector<TString>::iterator it_trig = cTriggerNames.begin(); it_trig != cTriggerNames.end(); ++it_trig){
	  if(Ntp->TriggerAccept(*it_trig)){
		  if ( value.at(TriggerOk) == TriggerOkDummy )
			  value.at(TriggerOk) = it_trig - cTriggerNames.begin();
		  else // more than 1 trigger fired, save this separately
			  value.at(TriggerOk) = cTriggerNames.size();
	  }
  }
  pass.at(TriggerOk) = (value.at(TriggerOk) >= cut.at(TriggerOk));
  if(id == DataMCType::DY_mutau_embedded) pass.at(TriggerOk) = true;
  if(verbose){
	  std::cout << "value at TriggerOk: " <<value.at(TriggerOk) << std::endl;
	  std::cout << "pass at TriggerOk: " <<pass.at(TriggerOk) << std::endl;
  }
  
  // Muon cuts
  if(verbose) std::cout << "Cut on MuonID" << std::endl;
  std::vector<int> selectedMuonsId;
  selectedMuonsId.clear();
  for(unsigned i_mu=0;i_mu<Ntp->NMuons();i_mu++){
	  if( selectMuon_Id(i_mu,selVertex) ) {
		  selectedMuonsId.push_back(i_mu);
	  }
  }
  value.at(NMuId)=selectedMuonsId.size();
  pass.at(NMuId)=(value.at(NMuId)>=cut.at(NMuId));
  if(verbose){
	  std::cout << "Number of Muons: " << Ntp->NMuons() << std::endl;
	  std::cout << "value at NMuId: " <<value.at(NMuId) << std::endl;
	  std::cout << "pass at NMuId: " <<pass.at(NMuId) << std::endl;
  }

  if(verbose) std::cout << "Cut on Muon Kinematics" << std::endl;
  std::vector<int> selectedMuonsKin;
  selectedMuonsKin.clear();
  for(std::vector<int>::iterator it_mu = selectedMuonsId.begin(); it_mu != selectedMuonsId.end(); ++it_mu){
	  if( selectMuon_Kinematics(*it_mu) ) {
		  selectedMuonsKin.push_back(*it_mu);
	  }
  }
  value.at(NMuKin)=selectedMuonsKin.size();
  pass.at(NMuKin)=(value.at(NMuKin)>=cut.at(NMuKin));
  if(verbose){
	  std::cout << "value at NMuKin: " <<value.at(NMuKin) << std::endl;
	  std::cout << "pass at NMuKin: " <<pass.at(NMuKin) << std::endl;
  }

  if(verbose) std::cout << "Cut on Muon Isolation (Iso)" << std::endl;
  std::vector<int> selectedMuonsIso;
  selectedMuonsIso.clear();
  for(std::vector<int>::iterator it_mu = selectedMuonsKin.begin(); it_mu != selectedMuonsKin.end(); ++it_mu){
	  if( selectMuon_Isolation(*it_mu) ) {
		  if(selMuon_Iso == selMuon_IsoDummy) selMuon_Iso = *it_mu;
		  selectedMuonsIso.push_back(*it_mu);
	  }
  }
  value.at(NMuIso)=selectedMuonsIso.size();
  pass.at(NMuIso)=(value.at(NMuIso)>=cut.at(NMuIso));
  if(verbose){
	  std::cout << "value at NMuIso: " <<value.at(NMuIso) << std::endl;
	  std::cout << "pass at NMuIso: " <<pass.at(NMuIso) << std::endl;
  }

  if(verbose) std::cout << "Cut on Muon Isolation (Anti Iso)" << std::endl;
  std::vector<int> selectedMuonsAntiIso;
  selectedMuonsAntiIso.clear();
  for(std::vector<int>::iterator it_mu = selectedMuonsKin.begin(); it_mu != selectedMuonsKin.end(); ++it_mu){
	  if( selectMuon_AntiIsolation(*it_mu) ) {
		  if(selMuon_AntiIso == selMuon_AntiIsoDummy && selMuon_Iso == selMuon_IsoDummy) selMuon_AntiIso = *it_mu;
		  selectedMuonsAntiIso.push_back(*it_mu);
	  }
  }

  // Tau cuts
  if(verbose) std::cout << "Cut on TauID" << std::endl;
  std::vector<int> selectedTausId;
  selectedTausId.clear();
  for(unsigned i_tau=0; i_tau < Ntp->NPFTaus(); i_tau++){
	  if ( selectPFTau_Id(i_tau,selectedMuonsId) ){
		  selectedTausId.push_back(i_tau);
	  }
  }
  value.at(NTauId)=selectedTausId.size();
  pass.at(NTauId)=(value.at(NTauId)>=cut.at(NTauId));
  if(verbose){
	  std::cout << "value at NTauId: " <<value.at(NTauId) << std::endl;
	  std::cout << "pass at NTauId: " <<pass.at(NTauId) << std::endl;
  }

  if(verbose) std::cout << "Cut on Tau Kinematics" << std::endl;
  std::vector<int> selectedTausKin;
  selectedTausKin.clear();
  for(std::vector<int>::iterator it_tau = selectedTausId.begin(); it_tau != selectedTausId.end(); ++it_tau){
	  if ( selectPFTau_Kinematics(*it_tau) ){
		  selectedTausKin.push_back(*it_tau);
	  }
  }
  value.at(NTauKin)=selectedTausKin.size();
  pass.at(NTauKin)=(value.at(NTauKin)>=cut.at(NTauKin));
  if(verbose){
	  std::cout << "value at NTauKin: " <<value.at(NTauKin) << std::endl;
	  std::cout << "pass at NTauKin: " <<pass.at(NTauKin) << std::endl;
  }

  if(verbose) std::cout << "Cut on Tau Isolation" << std::endl;
  std::vector<int> selectedTausIso;
  selectedTausIso.clear();
  for(std::vector<int>::iterator it_tau = selectedTausKin.begin(); it_tau != selectedTausKin.end(); ++it_tau){
	  if ( selectPFTau_Isolation(*it_tau) ){
		  if(selTau == selTauDummy) selTau = *it_tau;
		  selectedTausIso.push_back(*it_tau);
	  }
  }
  value.at(NTauIso)=selectedTausIso.size();
  pass.at(NTauIso)=(value.at(NTauIso)>=cut.at(NTauIso));
  if(verbose){
	  std::cout << "value at NTauIso: " <<value.at(NTauIso) << std::endl;
	  std::cout << "pass at NTauIso: " <<pass.at(NTauIso) << std::endl;
  }

  // Charge of MuTau
  if(verbose) std::cout << "Cut on Charge of MuTau System" << std::endl;
  value.at(ChargeSum) = ChargeSumDummy;
  if(selTau != -1){
	  if(selMuon_Iso != selMuon_IsoDummy && selMuon_AntiIso == selMuon_AntiIsoDummy){
		  Charge = Ntp->Muon_Charge(selMuon_Iso) + Ntp->PFTau_Charge(selTau);
		  value.at(ChargeSum) = Charge;
  	  }
  	  else if(selMuon_Iso == selMuon_IsoDummy && selMuon_AntiIso != selMuon_AntiIsoDummy){
  		  Charge = Ntp->Muon_Charge(selMuon_AntiIso) + Ntp->PFTau_Charge(selTau);
  		  value.at(ChargeSum) = Charge;
  	  }
  	  else{
  		  Charge = ChargeSumDummy;
  	  }
  }
  else{
	  Charge = ChargeSumDummy;
  }
  pass.at(ChargeSum)=(value.at(ChargeSum)==cut.at(ChargeSum));
  if(verbose){
	  std::cout << "value at ChargeSum: " <<value.at(ChargeSum) << std::endl;
	  std::cout << "pass at ChargeSum: " <<pass.at(ChargeSum) << std::endl;
  }

  // MT calculation
  if(verbose) std::cout << "Calculation and Cut on MT distribution" << std::endl;
  double pT,phi,eTmiss,eTmPhi;
  double MT_TauMET, MT_MuMET_AntiIso;

  if(selMuon_Iso == selMuon_IsoDummy && selMuon_AntiIso == selMuon_AntiIsoDummy){
	  value.at(MT_MuMET) = MTDummy;
	  if(verbose) std::cout << "No Muon selected: neither isolated or anti isolated" << std::endl;
  }
  else if(selMuon_Iso != selMuon_IsoDummy && selMuon_AntiIso != selMuon_AntiIsoDummy){
	  value.at(MT_MuMET) = MTDummy;
	  std::cout << "CRITICAL: SELECTED MUON PASSED ISOLATION AND ANTI-ISOLATION CUT --> FIX" << std::endl;
  }
  else if(selMuon_Iso != selMuon_IsoDummy && selMuon_AntiIso == selMuon_AntiIsoDummy){
	  //std::cout << "selMuon_Iso is " << selMuon_Iso << std::endl;
	  eTmiss				 	= Ntp->MET_CorrMVAMuTau_et();
	  eTmPhi 					= Ntp->MET_CorrMVAMuTau_phi();
	  pT 						= Ntp->Muon_p4(selMuon_Iso).Pt();
	  phi						= Ntp->Muon_p4(selMuon_Iso).Phi();
	  value.at(MT_MuMET)		= Ntp->transverseMass(pT,phi,eTmiss,eTmPhi);
	  //std::cout << "MT_MuMET is " << MT_MuMET << std::endl;
  }
  else if(selMuon_AntiIso != selMuon_AntiIsoDummy && selMuon_Iso == selMuon_IsoDummy){
	  eTmiss				 	= Ntp->MET_CorrMVAMuTau_et();
	  eTmPhi 					= Ntp->MET_CorrMVAMuTau_phi();
	  pT 						= Ntp->Muon_p4(selMuon_AntiIso).Pt();
	  phi						= Ntp->Muon_p4(selMuon_AntiIso).Phi();
	  MT_MuMET_AntiIso			= Ntp->transverseMass(pT,phi,eTmiss,eTmPhi);
	  value.at(MT_MuMET) 		= MT_MuMET_AntiIso;
  }
  if(value.at(MT_MuMET) != MTDummy) pass.at(MT_MuMET)=(value.at(MT_MuMET)<cut.at(MT_MuMET));
  if(verbose){
	  std::cout << "value at MT_MuMET: " <<value.at(MT_MuMET) << std::endl;
	  std::cout << "pass at MT_MuMET: " <<pass.at(MT_MuMET) << std::endl;
  }

  // Tau Decay Mode
  if(verbose) std::cout << "Cut on Tau Decay Mode" << std::endl;
  if(selTau != selTauDummy){
	  value.at(TauDecayMode) = Ntp->PFTau_hpsDecayMode(selTau);
  }
  pass.at(TauDecayMode)=(value.at(TauDecayMode)>=cut.at(TauDecayMode));
  if(verbose){
	  std::cout << "value at TauDecayMode: " <<value.at(TauDecayMode) << std::endl;
	  std::cout << "pass at TauDecayMode: " <<pass.at(TauDecayMode) << std::endl;
  }

  // Tau FlightLength Significance
  if(verbose) std::cout << "Cut on Tau Flight Length Significance" << std::endl;
  value.at(TauFLSigma) = TauFLSigmaDummy;
  if(pass.at(TauDecayMode) && selTau != selTauDummy){
	  //std::cout << "selTau" << selTau << std::endl;
	  //std::cout << "Ntp->PFTau_TIP_primaryVertex_pos(selTau).Mag() " << Ntp->PFTau_TIP_primaryVertex_pos(selTau).Mag() << std::endl;
	  //std::cout << "Ntp->PFTau_TIP_hassecondaryVertex(selTau) " << Ntp->PFTau_TIP_hassecondaryVertex(selTau) << std::endl;
	  //std::cout << "before hassecondaryVertex(selTau)" << std::endl;
	  if(Ntp->PFTau_TIP_hassecondaryVertex(selTau) && pass.at(PrimeVtx)){
		  //std::cout << "after hassecondaryVertex(selTau)" << std::endl;
		  //std::cout << "Ntp->PFTau_TIP_secondaryVertex_pos(selTau).Mag()" << Ntp->PFTau_TIP_secondaryVertex_pos(selTau).Mag() << std::endl;
		  //std::cout << "Ntp->PFTau_FlightLength(selTau) " << Ntp->PFTau_FlightLength(selTau) << std::endl;
		  //std::cout << "Ntp->PF_Tau_FlightLegth3d_TauFrame_cov(selTau) " << Ntp->PF_Tau_FlightLegth3d_TauFrame_cov(selTau)(LorentzVectorParticle::vz,LorentzVectorParticle::vz) << std::endl;
		  //Ntp->PFTau_TIP_secondaryVertex_cov(selTau).Print();
		  //Ntp->PFTau_TIP_primaryVertex_cov(selTau).Print();
		  //Ntp->PFTau_FlightLength3d_cov(selTau).Print();
		  //Ntp->PF_Tau_FlightLegth3d_TauFrame_cov(selTau).Print();
		  //std::cout << "Ntp->PFTau_FlightLength(selTau) " << Ntp->PFTau_FlightLength(selTau) << std::endl;
		  //std::cout << "Ntp->PFTau_FlightLength_error(selTau) " << Ntp->PFTau_FlightLength_error(selTau) << std::endl;
		  //std::cout << "Ntp->PFTau_FlightLength(selTau)/Ntp->PFTau_FlightLength_error(selTau) " << Ntp->PFTau_FlightLength(selTau)/Ntp->PFTau_FlightLength_error(selTau) << std::endl;
		  //std::cout << "Ntp->PFTau_FlightLength_error(selTau) " << Ntp->PFTau_FlightLength_error(selTau) << std::endl;
		  if(Ntp->PFTau_p4(selTau, tau_corr).Vect().Dot(Ntp->PFTau_FlightLength3d(selTau)) < 0){
			  value.at(TauFLSigma) = -Ntp->PFTau_FlightLength_significance(selTau);
		  }
		  else{
			  value.at(TauFLSigma) = Ntp->PFTau_FlightLength_significance(selTau);
		  }
	  }
  }
  pass.at(TauFLSigma) = (value.at(TauFLSigma)>=cut.at(TauFLSigma));
  if(verbose){
	  std::cout << "value at TauFLSigma: " <<value.at(TauFLSigma) << std::endl;
	  std::cout << "pass at TauFLSigma: " <<pass.at(TauFLSigma) << std::endl;
  }

  //////////////////////////////////////////////////////////////////////
  //*************************END OF SELECTION*************************//
  //////////////////////////////////////////////////////////////////////


  if(selTau == selTauDummy){
	  //std::cout << "selTau is " << selTau << std::endl;
	  MT_TauMET = MTDummy;
  }
  else{
	  //std::cout << "selTau is " << selTau << std::endl;
	  eTmiss 		= Ntp->MET_CorrMVAMuTau_et();
	  eTmPhi 		= Ntp->MET_CorrMVAMuTau_phi();
	  pT 			= Ntp->PFTau_p4(selTau, tau_corr).Pt();
	  phi			= Ntp->PFTau_p4(selTau, tau_corr).Phi();
	  MT_TauMET		= Ntp->transverseMass(pT,phi,eTmiss,eTmPhi);
	  //std::cout << "MT_TauMET is " << MT_TauMET << std::endl;
  }

  // Mvis
  if(verbose) std::cout << "Calculation of Mvis" << std::endl;
  double Mvis;

  if(selTau != selTauDummy){
	  if(selMuon_Iso != selMuon_IsoDummy && selMuon_AntiIso == selMuon_AntiIsoDummy){
		  Mvis = (Ntp->PFTau_p4(selTau, tau_corr) + Ntp->Muon_p4(selMuon_Iso)).M();
	  }
	  else if(selMuon_Iso == selMuon_IsoDummy && selMuon_AntiIso != selMuon_AntiIsoDummy){
		  Mvis = (Ntp->PFTau_p4(selTau, tau_corr) + Ntp->Muon_p4(selMuon_AntiIso)).M();
	  }
	  else{
		  Mvis = MvisDummy;
	  }
  }
  else{
	  Mvis = MvisDummy;
  }

  // Weights
  double wobs=1;
  double w=1;
  if(!Ntp->isData()){
	  if(id != DataMCType::DY_mutau_embedded){
		  w *= Ntp->PUWeightFineBins();
	  }
	  if(id == DataMCType::DY_mutau_embedded){
		  w *= Ntp->TauSpinnerWeight();
		  w *= Ntp->MinVisPtFilter();
	  }
	  if(selMuon_Iso != selMuon_IsoDummy){
		  w *= RSF->HiggsTauTau_MuTau_Id_Mu(Ntp->Muon_p4(selMuon_Iso));
		  w *= RSF->HiggsTauTau_MuTau_Iso_Mu(Ntp->Muon_p4(selMuon_Iso));
		  w *= RSF->HiggsTauTau_MuTau_Trigger_Mu(Ntp->Muon_p4(selMuon_Iso));
	  }
	  if(selTau != selTauDummy){
		  w *= RSF->HiggsTauTau_MuTau_Trigger_Tau(Ntp->PFTau_p4(selTau));
		  if(Ntp->PFTau_hpsDecayMode(selTau) == 0){
			  w *= OneProngNoPiWeight;
		  }
	  }
  }
  else{w=1;}

  // W+Jets BG Method
  if(verbose) std::cout << "W+Jets BG Method" << std::endl;
  std::vector<int> exclude_cuts;
  exclude_cuts.clear();
  exclude_cuts.push_back(ChargeSum);
  exclude_cuts.push_back(MT_MuMET);
  if(passAllBut(exclude_cuts)){
  	  if(pass.at(ChargeSum)){ //Opposite Sign WJets yield (bin #1 with value 1)
  		  if(value.at(MT_MuMET) > SB_lowerLimit && value.at(MT_MuMET) < SB_upperLimit){
  			  NSB.at(t).Fill(1.,w);
  		  }
  	  }
  	  if(!pass.at(ChargeSum) && Charge != ChargeSumDummy){ //Same Sign WJets yield (bin #2 with value 2)
  		  if(value.at(MT_MuMET) > SB_lowerLimit && value.at(MT_MuMET) < SB_upperLimit){
  			  NSB.at(t).Fill(2.,w);
  		  }
  	  }
  }

  // QCD ABCD BG Method
  /*******************
   *        |        *
   *    C   |    D   *  SS
   *        |        *       S
   * ----------------*------ i
   *        |        *       g
   *    A   |    B   *  OS   n
   *        |        *
   *******************
   *  Iso   | AntiIso
   *
   *     relIso(mu)
   */
  if(verbose) std::cout << "QCD ABCD BG Method" << std::endl;
  exclude_cuts.push_back(NMuIso);
  if(passAllBut(exclude_cuts)){
	  if(pass.at(NMuIso) && selMuon_Iso != selMuon_IsoDummy){
		  if(pass.at(ChargeSum)){ //A --> Signal-Selection (pass all cuts except MT_MuMET)
			  QCD_MT_MuMET_A.at(t).Fill(value.at(MT_MuMET),w);
	 	 	  QCD_MET_A.at(t).Fill(Ntp->MET_CorrMVAMuTau_et(), w);
	 	 	  if(pass.at(MT_MuMET)){
	 	 		  NQCD.at(t).Fill(1.,w);
	 	 	  }
	 	  }
	 	  if(!pass.at(ChargeSum) && value.at(ChargeSum) != ChargeSumDummy){ //C
	 	 	  QCD_MT_MuMET_C.at(t).Fill(value.at(MT_MuMET),w);
	 	      QCD_MET_C.at(t).Fill(Ntp->MET_CorrMVAMuTau_et(),w);
	 	 	  if(pass.at(MT_MuMET)){
	 	 		  NQCD.at(t).Fill(3.,w);
	 	 	  }
	 	  }
	  }
	  if(!pass.at(NMuIso) && selMuon_AntiIso != selMuon_AntiIsoDummy){
		  if(pass.at(ChargeSum)){ //B
	 		  QCD_MT_MuMET_B.at(t).Fill(MT_MuMET_AntiIso,w);
	 		  QCD_MET_B.at(t).Fill(Ntp->MET_CorrMVAMuTau_et(),w);
	 	 	  if(pass.at(MT_MuMET)){
	 	 		  NQCD.at(t).Fill(2.,w);
	 	 	  }
	 	  }
	 	  if(!pass.at(ChargeSum) && Charge != ChargeSumDummy){ //D
	 	 	  QCD_MT_MuMET_D.at(t).Fill(MT_MuMET_AntiIso,w);
	 	 	  QCD_MET_D.at(t).Fill(Ntp->MET_CorrMVAMuTau_et(),w);
	 	 	  if(pass.at(MT_MuMET)){
	 	 		  NQCD.at(t).Fill(4.,w);
	 	 	  }
			  if(id == DataMCType::Data){
				 QCD_MT_MuMET_D.at(HConfig.GetType(DataMCType::QCD)).Fill(MT_MuMET_AntiIso,w);
				 QCD_MET_D.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->MET_CorrMVAMuTau_et(),w);
				 QCD_MT_MuMET_A.at(HConfig.GetType(DataMCType::QCD)).Fill(MT_MuMET_AntiIso,w);
				 QCD_MET_A.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->MET_CorrMVAMuTau_et(),w);
			  }
	 	  }
	  }
  }

  if(!pass.at(NMuIso) && selMuon_AntiIso != selMuon_AntiIsoDummy){
	  if(!pass.at(ChargeSum) && Charge != ChargeSumDummy){
		  if(id == DataMCType::Data){
			 pass.at(ChargeSum) = true;
			 pass.at(NMuIso) = true;
			 bool QCD_status = AnalysisCuts(HConfig.GetType(DataMCType::QCD),w,wobs);
			 if(QCD_status){
				NVtx.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->NVtx(),w);
				unsigned int nGoodVtx=0;
				for(unsigned int i=0;i<Ntp->NVtx();i++){
					NTrackperVtx.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->Vtx_Track_idx(i).size(),w);
					if(Ntp->isVtxGood(i))nGoodVtx++;
				}
				NQCD.at(HConfig.GetType(DataMCType::QCD)).Fill(1.,w);
				NGoodVtx.at(HConfig.GetType(DataMCType::QCD)).Fill(nGoodVtx,w);
				NMtTauMET.at(HConfig.GetType(DataMCType::QCD)).Fill(MT_TauMET,w);
				NMvis.at(HConfig.GetType(DataMCType::QCD)).Fill(Mvis,w);
				Mu_pt.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->Muon_p4(selMuon_AntiIso).Pt(),w);
				Mu_phi.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->Muon_p4(selMuon_AntiIso).Phi(),w);
				Mu_eta.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->Muon_p4(selMuon_AntiIso).Eta(),w);
				Tau_pt.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_p4(selTau, tau_corr).Pt(),w);
				Tau_phi.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_p4(selTau, tau_corr).Phi(),w);
				Tau_eta.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_p4(selTau, tau_corr).Eta(),w);
				MET_phi.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->MET_CorrMVAMuTau_phi(),w);
				A1mass.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_p4(selTau, tau_corr).M(),w);
				A1mass10GeV.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_p4(selTau, tau_corr).M(),w);
				if(value.at(TauFLSigma) != TauFLSigmaDummy){
					TauFL.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_FlightLength(selTau),w);
					if(value.at(TauFLSigma) < 0){
						TauFLSigned.at(HConfig.GetType(DataMCType::QCD)).Fill(-Ntp->PFTau_FlightLength(selTau),w);
					}
					else if(value.at(TauFLSigma) >= 0){
						TauFLSigned.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_FlightLength(selTau),w);
					}
				}
			}
			pass.at(ChargeSum) = false;
			pass.at(NMuIso) = false;
		  }
	  }
  }

  bool status=AnalysisCuts(t,w,wobs);
  if(verbose){
	  std::cout << "------------------------" << std::endl;
	  if(status){
		  std::cout << "!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		  std::cout << "Event passed all cuts" << std::endl;
		  std::cout << "!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	  }
	  else std::cout << "Event failed selection" << std::endl;
	  std::cout << "------------------------" << std::endl;
  }

  //DiTau Reco
  if(status && value.at(TauFLSigma) != TauFLSigmaDummy){
	  for(unsigned Ambiguity=0; Ambiguity<3; Ambiguity++){
		  LorentzVectorParticle theTau;
		  std::vector<LorentzVectorParticle> daughter;
		  double LC_chi2;
		  double phisign;
		  LorentzVectorParticle theZ;
		  LorentzVectorParticle TauA1 = Ntp->PFTau_a1_lvp(selTau);
		  std::vector<LorentzVectorParticle> Daughters;
		  int Niterat;
		  double csum;

		  if(Ntp->ThreeProngTauFit(selTau, Ambiguity, theTau, daughter, LC_chi2, phisign)){
			  Ntp->EventFit(Ambiguity, selMuon_Iso, TauA1, theZ, Daughters, LC_chi2, Niterat, csum);
			  if(theZ.Mass() >= 0){
				  std::cout << theZ.Mass() << std::endl;
			  	  std::cout << Niterat << std::endl;
			  	  std::cout << csum << std::endl;
			  	  std::cout << LC_chi2 << std::endl;
			  }
		  }
	  }
  }

  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
	  NVtx.at(t).Fill(Ntp->NVtx(),w);
	  unsigned int nGoodVtx=0;
	  for(unsigned int i=0;i<Ntp->NVtx();i++){
		  NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
		  if(Ntp->isVtxGood(i))nGoodVtx++;
	  }
	  NGoodVtx.at(t).Fill(nGoodVtx,w);
	  NMtTauMET.at(t).Fill(MT_TauMET,w);
	  NMvis.at(t).Fill(Mvis,w);
	  Mu_pt.at(t).Fill(Ntp->Muon_p4(selMuon_Iso).Pt(),w);
	  Mu_phi.at(t).Fill(Ntp->Muon_p4(selMuon_Iso).Phi(),w);
	  Mu_eta.at(t).Fill(Ntp->Muon_p4(selMuon_Iso).Eta(),w);
	  Tau_pt.at(t).Fill(Ntp->PFTau_p4(selTau, tau_corr).Pt(),w);
	  Tau_phi.at(t).Fill(Ntp->PFTau_p4(selTau, tau_corr).Phi(),w);
	  Tau_eta.at(t).Fill(Ntp->PFTau_p4(selTau, tau_corr).Eta(),w);
	  MET_phi.at(t).Fill(Ntp->MET_CorrMVAMuTau_phi(),w);
	  A1mass.at(t).Fill(Ntp->PFTau_p4(selTau, tau_corr).M(),w);
	  A1mass10GeV.at(t).Fill(Ntp->PFTau_p4(selTau, tau_corr).M(),w);
  }

  //Gen Studies
  if(id != DataMCType::Data && (id == DataMCType::DY_tautau || (id >= DataMCType::H_tautau && id <= DataMCType::H_tautau_WHZHTTH))){
	  /*
	  if(id == DataMCType::DY_tautau) std::cout << "Ditau event" <<std::endl;;
	  int NTaus = 0;
	  for(unsigned i=0; i<Ntp->NMCParticles(); i++){
		  if(fabs(Ntp->MCParticle_pdgid(i))==15) NTaus++;
	  }
	  std::cout << "Number of Gen Taus" << NTaus <<std::endl;
	  std::cout << "Number of Gen Taus in MCTaus" << Ntp->NMCTaus() <<std::endl;
	  */
	  TLorentzVector GenTaumu, GenTauh, GenZH, GenMu, GenA1;
	  bool hasA1 = false;
	  bool hasMu = false;
	  bool hasZH = false;
	  if(Ntp->NMCTaus() == 2){
		  for(int i=0; i<Ntp->NMCTaus(); i++){
			  if(Ntp->MCTau_JAK(i) == 2){//Tau->Muon
				  GenTaumu = Ntp->MCTau_p4(i);
				  for(int j=0; j<Ntp->NMCTauDecayProducts(i); j++){
					  if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::mu_minus){
						  GenMu = Ntp->MCTauandProd_p4(i,j);
						  hasMu = true;
					  }
				  }
			  }
			  else if(Ntp->MCTau_JAK(i) != 2){
				  GenTauh = Ntp->MCTau_p4(i);
				  //std::cout << "--------" << std::endl;
				  for(int j=0; j<Ntp->NMCTauDecayProducts(i); j++){
					  //std::cout << "PDG ID of decay particle " << j << ": " << Ntp->MCTauandProd_pdgid(i,j) <<std::endl;
					  if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::a_1_plus){
						  GenA1 = Ntp->MCTauandProd_p4(i,j);
						  hasA1 = true;
					  }
				  }
				  if(status && value.at(TauFLSigma) != TauFLSigmaDummy){
					  double dPhi = Ntp->PFTau_FlightLength3d(selTau).Phi() - GenTauh.Phi();
					  TauFLSigma_vs_Res.at(t).Fill(value.at(TauFLSigma), dPhi);
					  //std::cout << "--------" << std::endl;
					  for(unsigned j=0; j<80; j++){
						  double TauFLSigmaCut = double(j)*0.5 -10;
						  if(value.at(TauFLSigma) >= TauFLSigmaCut){
							  TauFLSigmaCut_vs_Res.at(t).Fill(TauFLSigmaCut, dPhi);
						  }
					  }
				  }
			  }
		  }
		  for(unsigned i=0; i<Ntp->NMCParticles(); i++){
			  if(Ntp->MCParticle_pdgid(i) == PDGInfo::Z0){
				  GenZH = Ntp->MCParticle_p4(i);
				  hasZH = true;
			  }
			  else if(Ntp->MCParticle_pdgid(i) == PDGInfo::Higgs0){
				  GenZH = Ntp->MCParticle_p4(i);
				  hasZH = true;
			  }
		  }
		  if(hasMu && hasA1 && hasZH){
			  GJ_Taumu.at(t).Fill(GenMu.Vect().Angle(GenTaumu.Vect()),w);
			  GJ_Tauh.at(t).Fill(GenA1.Vect().Angle(GenTauh.Vect()),w);
			  Phi_genTaumu.at(t).Fill(GenTaumu.Phi(),w);
			  Theta_genTaumu.at(t).Fill(GenTaumu.Theta(),w);
			  Phi_genTauh.at(t).Fill(GenTauh.Phi(),w);
			  Theta_genTauh.at(t).Fill(GenTauh.Theta(),w);

			  double dPhi_genTaus = GenTauh.Phi() - GenTaumu.Phi();
			  //if(dPhi_genTaus > TMath::Pi()) dPhi_genTaus = 2*TMath::Pi() - dPhi_genTaus;
			  TLorentzVector DiTau = GenTauh + GenTaumu;
			  Pt_ZGen.at(t).Fill(GenZH.Pt(),w);
			  dPhi_DiTauGen.at(t).Fill(dPhi_genTaus,w);
			  Pt_DiTauGen.at(t).Fill(DiTau.Pt(),w);
			  M_ZGen.at(t).Fill(GenZH.M(),w);
			  Pt_vs_dPhi_DiTauGen.at(t).Fill(DiTau.Pt(), dPhi_genTaus);

			  TLorentzVector GenTaumu_PtBalance = GenTauh;
			  GenTaumu_PtBalance.SetX(GenTauh.X()*-1);
			  GenTaumu_PtBalance.SetY(GenTauh.Y()*-1);

			  double M_DiTau = (GenTaumu_PtBalance + GenTauh).M();
			  double Delta_M_DiTau = M_DiTau - GenZH.M();

			  M_DiTauPtBalance.at(t).Fill(M_DiTau,w);
			  dM_DiTau.at(t).Fill(Delta_M_DiTau, w);

			  double dPt = GenTaumu_PtBalance.Pt() - GenTaumu.Pt();
			  double dP_GenTaumu = GenTaumu_PtBalance.P() - GenTaumu.P();
			  double GJ_Angle_Tauh = GenA1.Vect().Angle(GenTauh.Vect());
			  double dP_Tauh = GenA1.E()*sqrt(pow(pow(GenTauh.M(),2) - pow(GenA1.M(),2),2) - pow(2*GenA1.P()*GenTauh.M()*sin(GJ_Angle_Tauh),2))/(pow(GenA1.M(),2) + pow(sin(GJ_Angle_Tauh)*GenA1.P(),2));

			  dPt_GenTaumuPtBalance.at(t).Fill(dPt,w);
			  dP_GenTaumuPtBalance.at(t).Fill(dP_GenTaumu,w);
			  dP_GenTauh.at(t).Fill(dP_Tauh,w);
			  dP_GenTauMuPtBalance_vs_dPTauh.at(t).Fill(dP_GenTaumu, dP_Tauh);

			  if(status && value.at(TauFLSigma) != TauFLSigmaDummy){
				  TVector3 POCAPV_dir = Ntp->Muon_Poca(selMuon_Iso) - Ntp->PFTau_TIP_primaryVertex_pos(selTau);
				  POCAPV_Mag.at(t).Fill(POCAPV_dir.Mag());
				  double Mvis_recoA1genMu = (Ntp->PFTau_p4(selTau, tau_corr) + GenMu).M();
				  double Mvis_genA1recoMu = (GenA1 + Ntp->Muon_p4(selMuon_Iso)).M();
				  double Mvis_recoA1genTaumu = (Ntp->PFTau_p4(selTau, tau_corr) + GenTaumu).M();
				  double Mvis_genTauhrecoMu = (GenTauh + Ntp->Muon_p4(selMuon_Iso)).M();

				  Mvis_SignalOnly.at(t).Fill(Mvis,w);
				  Mvis_SignalOnly_genMu.at(t).Fill(Mvis_recoA1genMu,w);
				  Mvis_SignalOnly_genA1.at(t).Fill(Mvis_genA1recoMu,w);
				  Mvis_SignalOnly_genTaumu.at(t).Fill(Mvis_recoA1genTaumu,w);
				  Mvis_SignalOnly_genTauh.at(t).Fill(Mvis_genTauhrecoMu,w);

				  Phi_POCAPV.at(t).Fill(POCAPV_dir.Phi(),w);
				  Theta_POCAPV.at(t).Fill(POCAPV_dir.Theta(),w);
				  double dPhi = POCAPV_dir.Phi() - GenTaumu.Phi();
				  dPhi_POCAPV_genTaumu.at(t).Fill(dPhi, w);
				  double dTheta = POCAPV_dir.Theta() - GenTaumu.Theta();
				  dTheta_POCAPV_genTaumu.at(t).Fill(dTheta,w);

				  Phi_SVPV.at(t).Fill(Ntp->PFTau_FlightLength3d(selTau).Phi(),w);
				  Theta_SVPV.at(t).Fill(Ntp->PFTau_FlightLength3d(selTau).Theta(),w);
				  dPhi = Ntp->PFTau_FlightLength3d(selTau).Phi() - GenTauh.Phi();
				  dPhi_SVPV_genTauh.at(t).Fill(dPhi, w);
				  dTheta = Ntp->PFTau_FlightLength3d(selTau).Theta() - GenTauh.Theta();
				  dTheta_SVPV_genTauh.at(t).Fill(dTheta,w);
				  Angle_SVPV_genTauh.at(t).Fill(Ntp->PFTau_FlightLength3d(selTau).Angle(GenTauh.Vect()));

				  TVector3 MinusPFTau_FlightLength3d = - Ntp->PFTau_FlightLength3d(selTau);
				  dPhi = MinusPFTau_FlightLength3d.Phi() - GenTaumu.Phi();
				  dPhi_MinusSVPV_genTaumu.at(t).Fill(dPhi, w);
				  dTheta = MinusPFTau_FlightLength3d.Theta() - GenTaumu.Theta();
				  dTheta_MinusSVPV_genTaumu.at(t).Fill(dTheta, w);
				  Angle_MinusSVPV_genTaumu.at(t).Fill(MinusPFTau_FlightLength3d.Angle(GenTaumu.Vect()));
			  }
		  }
	  }
  }
  if(passAllBut(TauFLSigma)){
	  TauFLSigmaSigned.at(t).Fill(value.at(TauFLSigma),w);
	  if(value.at(TauFLSigma)<0 && value.at(TauFLSigma) != TauFLSigmaDummy){
		  TauFLSigmaUnsigned.at(t).Fill(-value.at(TauFLSigma),w);
		  TauFL.at(t).Fill(Ntp->PFTau_FlightLength(selTau),w);
		  TauFLSigned.at(t).Fill(-Ntp->PFTau_FlightLength(selTau),w);
	  }
	  else if(value.at(TauFLSigma)>0){
		  TauFLSigmaUnsigned.at(t).Fill(value.at(TauFLSigma),w);
		  TauFL.at(t).Fill(Ntp->PFTau_FlightLength(selTau),w);
		  TauFLSigned.at(t).Fill(Ntp->PFTau_FlightLength(selTau),w);
	  }
  }
  std::vector<int> exclude_TauFLQCD;
  exclude_TauFLQCD.clear();
  exclude_TauFLQCD.push_back(NMuIso);
  exclude_TauFLQCD.push_back(ChargeSum);
  exclude_TauFLQCD.push_back(TauFLSigma);
  if(passAllBut(exclude_TauFLQCD)){
	  if(!pass.at(NMuIso) && selMuon_AntiIso != selMuon_AntiIsoDummy){
		  if(!pass.at(ChargeSum) && Charge != ChargeSumDummy){
			  if(id == DataMCType::Data){
				  TauFLSigmaSigned.at(HConfig.GetType(DataMCType::QCD)).Fill(value.at(TauFLSigma),w);
				  if(value.at(TauFLSigma)<0 && value.at(TauFLSigma) != TauFLSigmaDummy){
					  TauFLSigmaUnsigned.at(HConfig.GetType(DataMCType::QCD)).Fill(-value.at(TauFLSigma),w);
					  TauFL.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_FlightLength(selTau),w);
					  TauFLSigned.at(HConfig.GetType(DataMCType::QCD)).Fill(-Ntp->PFTau_FlightLength(selTau),w);
				  }
				  else if(value.at(TauFLSigma)>0){
					  TauFLSigmaUnsigned.at(HConfig.GetType(DataMCType::QCD)).Fill(value.at(TauFLSigma),w);
					  TauFL.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_FlightLength(selTau),w);
					  TauFLSigned.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_FlightLength(selTau),w);
				  }
			  }
		  }
	  }
  }
  std::vector<int> exclude_MTDecay;
  exclude_MTDecay.clear();
  exclude_MTDecay.push_back(NMuIso);
  exclude_MTDecay.push_back(ChargeSum);
  exclude_MTDecay.push_back(MT_MuMET);
  exclude_MTDecay.push_back(TauDecayMode);
  exclude_MTDecay.push_back(TauFLSigma);
  if(passAllBut(exclude_MTDecay)){
	  if(pass.at(NMuIso) && pass.at(ChargeSum)){
		  if(pass.at(MT_MuMET)) MvisIncl.at(t).Fill(Mvis,w);
		  MTMuMETIncl.at(t).Fill(value.at(MT_MuMET),w);
		  if(pass.at(TauDecayMode)){
			  if(pass.at(MT_MuMET))  Mvis3Prong.at(t).Fill(Mvis,w);
			  MTMuMET3Prong.at(t).Fill(value.at(MT_MuMET),w);
		  }
	  	  if(!pass.at(TauDecayMode)){
	  		  if(pass.at(MT_MuMET)) Mvis1Prong.at(t).Fill(Mvis,w);
			  MTMuMET1Prong.at(t).Fill(value.at(MT_MuMET),w);
	  	  }
	  }
	  else if(!pass.at(NMuIso) && selMuon_AntiIso != selMuon_AntiIsoDummy && !pass.at(ChargeSum) && Charge != ChargeSumDummy){
		  if(pass.at(MT_MuMET)) MvisIncl.at(HConfig.GetType(DataMCType::QCD)).Fill(Mvis,w);
		  MTMuMETIncl.at(HConfig.GetType(DataMCType::QCD)).Fill(value.at(MT_MuMET),w);
		  if(pass.at(TauDecayMode)){
			  if(pass.at(MT_MuMET)) Mvis3Prong.at(HConfig.GetType(DataMCType::QCD)).Fill(Mvis,w);
			  MTMuMET3Prong.at(HConfig.GetType(DataMCType::QCD)).Fill(value.at(MT_MuMET),w);
		  }
	  	  if(!pass.at(TauDecayMode)){
	  		  if(pass.at(MT_MuMET)) Mvis1Prong.at(HConfig.GetType(DataMCType::QCD)).Fill(Mvis,w);
			  MTMuMET1Prong.at(HConfig.GetType(DataMCType::QCD)).Fill(value.at(MT_MuMET),w);
	  	  }
	  }
  }

}//final bracket of DoEvent


void  ZToTaumuTauh::Finish(){
  std::cout << "Starting Finish" << std::endl;
  if(mode==RECONSTRUCT){
	  std::cout << "Enter mode==RECONSTRUCT" << std::endl;
	  SkimConfig SC;
	  SC.ApplySkimEfficiency(types,Npassed,Npassed_noweight);
		//for(unsigned i=0; i<Npassed.size();i++){
		//	nevents_noweight_default.push_back(Npassed_noweight.at(i).GetBinContent(1));
		//}

  	  // Scale DY embedded sample to the yield of DY MC and scale DY MC to 0 afterwards
	  double NDY_tautau(Npassed.at(HConfig.GetType(DataMCType::DY_tautau)).GetBinContent(NCuts)*CrossSectionandAcceptance.at(HConfig.GetType(DataMCType::DY_tautau))*Lumi/Npassed.at(HConfig.GetType(DataMCType::DY_tautau)).GetBinContent(0));
	  double NDY_tautau_Emb(Npassed.at(HConfig.GetType(DataMCType::DY_mutau_embedded)).GetBinContent(NCuts));

	  std::cout << "NDY_tautau: " << NDY_tautau << std::endl;
	  std::cout << "NDY_tautau_Emb: " << NDY_tautau_Emb << std::endl;
	  std::cout << "NDY_tautau/NDY_tautau_Emb: " << NDY_tautau/NDY_tautau_Emb << std::endl;

	  int Exclude_DY_ID;
	  if(Use_Embedded){
		  if(NDY_tautau>0 && NDY_tautau_Emb>0){
			  ScaleAllHistOfType(HConfig.GetType(DataMCType::DY_mutau_embedded),NDY_tautau/NDY_tautau_Emb);
	  	  	  ScaleAllHistOfType(HConfig.GetType(DataMCType::DY_tautau),0);
	  	  	  Exclude_DY_ID = DataMCType::DY_tautau;
		  }
	  }
	  else{
		  ScaleAllHistOfType(HConfig.GetType(DataMCType::DY_mutau_embedded),0);
  	  	  Exclude_DY_ID = DataMCType::DY_mutau_embedded;
	  }
  	  //CrossSectionandAcceptance.at(HConfig.GetType(DataMCType::DY_tautau)) = 0;
  	  //HConfig.SetCrossSection(DataMCType::DY_tautau,0);
  	  //HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);

  	  std::vector<double> SB_Integral;
  	  double SB_Integral_Data_minus_MC = 0;
  	  double SB_Integral_WJets = 0;

  	  std::vector<double> SB_Counting_OS;
  	  double SB_Counting_Data_minus_MC_OS = 0;
  	  double SB_Counting_WJets_OS = 0;

  	  std::vector<double> SB_Counting_SS;
  	  double SB_Counting_Data_minus_MC_SS = 0;
  	  double SB_Counting_WJets_SS = 0;

  	  std::vector<double> QCD_Integral_B;
  	  double QCD_Integral_B_Data_minus_MC = 0;
  	  std::vector<double> QCD_Integral_C;
  	  double QCD_Integral_C_Data_minus_MC = 0;
  	  std::vector<double> QCD_Integral_D;
  	  double QCD_Integral_D_Data_minus_MC = 0;

  	  // Get Yields in W+Jets-Sideband for OS/SS
  	  std::cout << "W+Jets Background Method " << std::endl;
  	  for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
  	  	  int Bin_SB_lowerLimit = Nminus1.at(MT_MuMET).at(i).FindFixBin(SB_lowerLimit);
  	  	  int Bin_SB_upperLimit = Nminus1.at(MT_MuMET).at(i).FindFixBin(SB_upperLimit) -1;
  		  SB_Integral.push_back(Nminus1.at(MT_MuMET).at(i).Integral(Bin_SB_lowerLimit, Bin_SB_upperLimit));
  		  SB_Counting_OS.push_back(NSB.at(i).GetBinContent(1));
  		  SB_Counting_SS.push_back(NSB.at(i).GetBinContent(2));
  		  if(CrossSectionandAcceptance.at(i)>0){
			  SB_Integral.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
			  SB_Counting_OS.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
			  SB_Counting_SS.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
		  }
		  std::cout << "SB integral 70-140GeV from Nminus1 plot for Sample " << i << " : " << SB_Integral.at(i) << std::endl;
		  std::cout << "SB integral 70-140GeV from Counting for OS for Sample " << i << " : " << SB_Counting_OS.at(i) << std::endl;
		  std::cout << "SB integral 70-140GeV from Counting for SS for Sample " << i << " : " << SB_Counting_SS.at(i) << std::endl;
  	  }
  	  for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
  		  if(HConfig.GetID(i) == DataMCType::W_lnu || HConfig.GetID(i) == DataMCType::W_taunu){
  			  SB_Integral_WJets				+= SB_Integral.at(i);
  			  SB_Counting_WJets_OS			+= SB_Counting_OS.at(i);
  			  SB_Counting_WJets_SS			+= SB_Counting_SS.at(i);
  		  }
  		  else if(HConfig.GetID(i) == DataMCType::Data){
  			  SB_Integral_Data_minus_MC 	+= SB_Integral.at(i);
  			  SB_Counting_Data_minus_MC_OS	+= SB_Counting_OS.at(i);
  			  SB_Counting_Data_minus_MC_SS	+= SB_Counting_SS.at(i);
  		  }
  		  else if(HConfig.GetID(i) != Exclude_DY_ID){
  			  SB_Integral_Data_minus_MC		-= SB_Integral.at(i);
  			  SB_Counting_Data_minus_MC_OS	-= SB_Counting_OS.at(i);
  			  SB_Counting_Data_minus_MC_SS	-= SB_Counting_SS.at(i);
  		  }
  	  }
  	  if(SB_Integral_Data_minus_MC > 0 && SB_Integral_WJets > 0){
  		  std::cout << "Scale Factor for W+Jets Sample with Histo.Integral Method: " << SB_Integral_Data_minus_MC/SB_Integral_WJets << std::endl;
  		  if(!Scaleby_Counting){
  	  		  std::cout << "Scaling by Histo.Integral Method" << std::endl;
  			  ScaleAllHistOfType(HConfig.GetType(DataMCType::W_lnu), SB_Integral_Data_minus_MC/SB_Integral_WJets);
  			  ScaleAllHistOfType(HConfig.GetType(DataMCType::W_taunu), SB_Integral_Data_minus_MC/SB_Integral_WJets);
  		  }
  	  }
  	  else{
		  std::cout << "SB_Integral_WJets is: " << SB_Integral_WJets << std::endl;
		  std::cout << "SB_Integral_Data_minus_MC is: " << SB_Integral_Data_minus_MC << std::endl;
  	  }
  	  if(SB_Counting_Data_minus_MC_OS > 0 && SB_Counting_WJets_OS > 0){
		  std::cout << "Scale Factor for W+Jets Sample with Counting Method for OS: " << SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS << std::endl;
		  std::cout << "Scaleby_Counting boolean is set to " << Scaleby_Counting << std::endl;
  		  if(Scaleby_Counting){
  	  		  std::cout << "Scaling by Counting Method" << std::endl;
  	  		  QCD_MT_MuMET_A.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  QCD_MT_MuMET_A.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  QCD_MT_MuMET_B.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  QCD_MT_MuMET_B.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  QCD_MET_A.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  QCD_MET_A.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  QCD_MET_B.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  QCD_MET_B.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  	  for(unsigned int i=0; i<Nminus1.size(); i++){
  		  		  Nminus0.at(i).at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  		  Nminus0.at(i).at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  		  Nminus1.at(i).at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  		  Nminus1.at(i).at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  	  }
  		  	  for(unsigned int i=0; i<Extradist1d_OS.size(); i++){
  		  		  Extradist1d_OS.at(i)->at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  		  Extradist1d_OS.at(i)->at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  	  }
  		  }
  	  }
  	  else{
		  std::cout << "SB_Counting_WJets_OS is: " << SB_Counting_WJets_OS << std::endl;
		  std::cout << "SB_Counting_Data_minus_MC_OS is: " << SB_Counting_Data_minus_MC_OS << std::endl;
  	  }
  	  if(SB_Counting_Data_minus_MC_SS > 0 && SB_Counting_WJets_SS > 0){
		  std::cout << "Scale Factor for W+Jets Sample with Counting Method for SS: " << SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS << std::endl;
		  std::cout << "Scaleby_Counting boolean is set to " << Scaleby_Counting << std::endl;
		  QCD_MT_MuMET_C.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
		  QCD_MT_MuMET_C.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  QCD_MT_MuMET_D.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  QCD_MT_MuMET_D.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  QCD_MET_C.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  QCD_MET_C.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  QCD_MET_D.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  QCD_MET_D.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
  	  }
  	  else{
		  std::cout << "SB_Counting_WJets_SS is: " << SB_Counting_WJets_SS << std::endl;
		  std::cout << "SB_Counting_Data_minus_MC_SS is: " << SB_Counting_Data_minus_MC_SS << std::endl;
  	  }
  	  // Get Yields in ABCD for QCD Scalefactor
  	  for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
  		  QCD_Integral_B.push_back(NQCD.at(i).GetBinContent(2));
  		  QCD_Integral_C.push_back(NQCD.at(i).GetBinContent(3));
  		  QCD_Integral_D.push_back(NQCD.at(i).GetBinContent(4));
  		  //QCD_Integral_B.push_back(QCD_MT_MuMET_B.at(i).Integral());
  		  //QCD_Integral_C.push_back(QCD_MT_MuMET_C.at(i).Integral());
  		  //QCD_Integral_D.push_back(QCD_MT_MuMET_D.at(i).Integral());
  		  if(HConfig.GetID(i) == DataMCType::W_lnu || HConfig.GetID(i) == DataMCType::W_taunu){
	  		  QCD_Integral_B.at(i) *= SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS;
	  		  QCD_Integral_C.at(i) *= SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS;
	  		  QCD_Integral_D.at(i) *= SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS;
  		  }
  		  if(CrossSectionandAcceptance.at(i)>0){
	  		  QCD_Integral_B.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	  		  QCD_Integral_C.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	  		  QCD_Integral_D.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
		  }
  	  }
  	  for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
  		  if(HConfig.GetID(i) == DataMCType::Data){
  			  QCD_Integral_B_Data_minus_MC	+= QCD_Integral_B.at(i);
  			  QCD_Integral_C_Data_minus_MC	+= QCD_Integral_C.at(i);
  			  QCD_Integral_D_Data_minus_MC	+= QCD_Integral_D.at(i);
  		  }
  		  else if(HConfig.GetID(i) != DataMCType::QCD && HConfig.GetID(i) != Exclude_DY_ID){
  			  QCD_Integral_B_Data_minus_MC	-= QCD_Integral_B.at(i);
  			  QCD_Integral_C_Data_minus_MC	-= QCD_Integral_C.at(i);
  			  QCD_Integral_D_Data_minus_MC	-= QCD_Integral_D.at(i);
  		  }
  	  }
  	  if(QCD_Integral_B_Data_minus_MC > 0 && QCD_Integral_C_Data_minus_MC > 0 && QCD_Integral_D_Data_minus_MC > 0){
		  std::cout << "Factor AntiIso OS/SS QCD Sample: " << QCD_Integral_B_Data_minus_MC/QCD_Integral_D_Data_minus_MC << std::endl;
		  std::cout << "Scale Factor for QCD Sample: " << QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC << std::endl;
		  double QCD_ScaleFactor = QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC;
		  QCD_MT_MuMET_A.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_ScaleFactor);
		  QCD_MET_A.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_ScaleFactor);
		  NQCD.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_ScaleFactor);
		  for(unsigned i=0; i<Nminus1.size(); i++){
			  Nminus0.at(i).at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_ScaleFactor);
		  	  Nminus1.at(i).at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_ScaleFactor);
		  }
		  for(unsigned i=0; i<Extradist1d_OS.size(); i++){
			  Extradist1d_OS.at(i)->at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_ScaleFactor);
		  }
  	  }
  	  else{
		  std::cout << "No QCD Scaling. Reason: " << std::endl;
		  std::cout << "QCD_Integral_B_Data_minus_MC is: " << QCD_Integral_B_Data_minus_MC << std::endl;
		  std::cout << "QCD_Integral_C_Data_minus_MC is: " << QCD_Integral_C_Data_minus_MC << std::endl;
		  std::cout << "QCD_Integral_D_Data_minus_MC is: " << QCD_Integral_D_Data_minus_MC << std::endl;
  	  }
  }
  // weight all Histograms
  std::cout << "Starting Selection::Finish" << std::endl;
  Selection::Finish();
}


/////////////////////////////////////////
// Definition of selection and helper functions
/////////////////////////////////////////

///////// Muons

bool ZToTaumuTauh::selectMuon_Id(unsigned i, unsigned vertex){
	if(	Ntp->isSelectedMuon(i,vertex,cMu_dxy,cMu_dz) &&
		(Ntp->GetMCID() == DataMCType::DY_mutau_embedded || Ntp->matchTrigger(Ntp->Muon_p4(i),cTriggerNames,"muon") < cMu_dRHltMatch)
			){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectMuon_Kinematics(unsigned i){
	if(	Ntp->Muon_p4(i).Pt() >= cMu_pt &&
		fabs(Ntp->Muon_p4(i).Eta()) <= cMu_eta
			){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectMuon_Isolation(unsigned i){
	if(Ntp->Muon_RelIso(i) < cMu_relIso){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectMuon_AntiIsolation(unsigned i){
	if(		Ntp->Muon_RelIso(i) < 0.5 &&
			Ntp->Muon_RelIso(i) > 0.2){
		return true;
	}
	return false;
}

///////// Taus
bool ZToTaumuTauh::selectPFTau_Id(unsigned i){
	if ( 	Ntp->PFTau_isHPSByDecayModeFinding(i) &&
			Ntp->PFTau_isHPSAgainstElectronsLoose(i) &&
			Ntp->PFTau_isHPSAgainstMuonTight(i)
			){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectPFTau_Id(unsigned i, std::vector<int> muonCollection){
	// check if tau is matched to a muon, if so this is not a good tau
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013#Sync_Issues
	for(std::vector<int>::iterator it_mu = muonCollection.begin(); it_mu != muonCollection.end(); ++it_mu){
	  if( Ntp->PFTau_p4(i).DeltaR(Ntp->Muon_p4(*it_mu)) < cMuTau_dR ) {
		  return false;
	  }
	}
	// trigger matching
	if(Ntp->GetMCID() != DataMCType::DY_mutau_embedded){
		if (Ntp->matchTrigger(Ntp->PFTau_p4(i),cTriggerNames,"tau") > cTau_dRHltMatch) {
			return false;
		}
	}

	if ( 	selectPFTau_Id(i) ){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectPFTau_Isolation(unsigned i){
	if(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(i) < cTau_IsoRaw){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectPFTau_Kinematics(unsigned i){
	if(	Ntp->PFTau_p4(i).Pt() >= cTau_pt &&
		fabs(Ntp->PFTau_p4(i).Eta()) <= cTau_eta
			){
		return true;
	}
	return false;
}
double ZToTaumuTauh::Reconstruct_hadronicTauEnergy(unsigned i){
	double TauEnergy,TauMomentumPlus, TauMomentumMinus, TauMomentum;
	TLorentzVector a1 = Ntp->PFTau_p4(i,tau_corr);
	double GJ_angle = a1.Angle(Ntp->PFTau_FlightLength3d(i));
	double val1 = (pow(PDG_Var::Tau_mass(),2.) + pow(a1.M(),2.))*a1.P()*cos(GJ_angle);
	double val2 = a1.Energy()*sqrt(pow(pow(PDG_Var::Tau_mass(),2.) - pow(a1.M(),2.),2.) - 4.*pow(a1.P()*PDG_Var::Tau_mass()*sin(GJ_angle),2.));
	TauMomentumPlus = (val1 + val2)/2./(pow(a1.M(),2) + pow(a1.P()*sin(GJ_angle),2.));
	TauMomentumMinus = (val1 - val2)/2./(pow(a1.M(),2) + pow(a1.P()*sin(GJ_angle),2.));
	TauEnergy = sqrt(pow(TauMomentum,2.) + pow(PDG_Var::Tau_mass(),2.));
	return TauEnergy;
}

