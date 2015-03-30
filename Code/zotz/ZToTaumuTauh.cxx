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
#include "SimpleFits/FitSoftware/interface/TrackTools.h"
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
	selMuonDummy = -1;
	selTauDummy = -1;
	ChargeSumDummy = -999;
	MTDummy = -999;
	MvisDummy = -999;
	TauFLSigmaDummy = -999;

	//WJets BG Method
	SB_lowerLimit = 70; //in GeV, Region > SB_lowerLimit dominated by W+Jets
	SB_upperLimit = 140; //in GeV, Region < SB_upperLimit ... not used currently (open end)
	Scaleby_Counting = true; // = false --> Scale by Integral

	//Set verbose boolean
	verbose = false;

	Use_Embedded = true;
}

ZToTaumuTauh::~ZToTaumuTauh(){
	for(unsigned int j=0; j<Npassed.size(); j++){
		std::cout << "ZToTaumuTauh::~ZToTaumuTauh Selection Summary before: "
		<< Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		<< Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts+1) << std::endl;
	}
	std::cout << "ZToTaumuTauh::~ZToTaumuTauh()" << std::endl;
}

void  ZToTaumuTauh::Configure(){
	// Setup Cut Values
	for(int i=0; i<NCuts;i++){
		cut.push_back(0);
		value.push_back(0);
		pass.push_back(false);
		if(i==TriggerOk)		cut.at(TriggerOk)=0;
		if(i==PrimeVtx)			cut.at(PrimeVtx)=1;
		if(i==NMuId)			cut.at(NMuId)=1;
		if(i==NMuKin)			cut.at(NMuKin)=1;
		if(i==NMuIso)			cut.at(NMuIso)=1;
		if(i==NTauId)			cut.at(NTauId)=1;
		if(i==NTauKin)			cut.at(NTauKin)=1;
		if(i==NTauIso)			cut.at(NTauIso)=1;
		if(i==ChargeSum)		cut.at(ChargeSum)=0;
		if(i==TauDecayMode)		cut.at(TauDecayMode)=10;//10
		if(i==TauFLSigma)		cut.at(TauFLSigma)=3;//3
		if(i==MT_MuMET)			cut.at(MT_MuMET)=30;
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
			title.at(i_cut)="TauFLSigma >= ";
			title.at(i_cut)+=cut.at(TauFLSigma);
			htitle=title.at(i_cut);
			htitle.ReplaceAll("$","");
			htitle.ReplaceAll("\\","#");
			hlabel="TauFLSigma";
			Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauFLSigma_",htitle,80,-10,30,hlabel,"Events"));
			Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauFLSigma_",htitle,80,-10,30,hlabel,"Events"));
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

	Tau_Mass_Inclusive=HConfig.GetTH1D(Name+"_Tau_Mass_Inclusive","Tau_Mass_Inclusive",170,-0.2,1.5,"PFTau mass","Events");
	Tau_Mass_sq_Inclusive=HConfig.GetTH1D(Name+"_Tau_Mass_sq_Inclusive","Tau_Mass_sq_Inclusive",170,-0.05,0.05,"PFTau mass squared","Events");
	Tau_Mass_Inclusive_NoTLV=HConfig.GetTH1D(Name+"_Tau_Mass_Inclusive_NoTLV","Tau_Mass_Inclusive_NoTLV",170,-0.2,1.5,"PFTau mass No TLV","Events");
	Tau_Mass_Inclusive_UnFitTracks=HConfig.GetTH1D(Name+"_Tau_Mass_Inclusive_UnFitTracks","Tau_Mass_Inclusive_UnFitTracks",170,-0.2,1.5,"Tau_Mass_Inclusive_UnFitTracks","Events");
	Tau_Mass_Inclusive_ReFitTracks=HConfig.GetTH1D(Name+"_Tau_Mass_Inclusive_ReFitTracks","Tau_Mass_Inclusive_ReFitTracks",170,-0.2,1.5,"Tau_Mass_Inclusive_ReFitTracks","Events");

	Tau_Mass_Difference_PFTau_UnFitTracks_3PS=HConfig.GetTH1D(Name+"_Tau_Mass_Difference_PFTau_UnFitTracks_3PS","Tau_Mass_Difference_PFTau_UnFitTracks_3PS",50,-0.1,0.1,"Tau_Mass_Difference_PFTau_UnFitTracks_3PS","Events");

	MET_phi=HConfig.GetTH1D(Name+"_MET_phi","MET_phi",30,-3.14159265359,3.14159265359,"MET #phi","Events");

	TauFL_NoTauFLSigmaCut=HConfig.GetTH1D(Name+"_TauFL_NoTauFLSigmaCut","TauFL_NoTauFLSigmaCut",100,-2,4,"tau flight length without tauflsigma cut","Events");
	TauFLSigned_NoTauFLSigmaCut=HConfig.GetTH1D(Name+"_TauFLSigned_NoTauFLSigmaCut","TauFLSigned_NoTauFLSigmaCut",100,-2,4,"tau flight length Signed without tauflsigma cut","Events");

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
	dPhi_SVPV_genTauh_vs_TauFL=HConfig.GetTH2D(Name+"_dPhi_SVPV_genTauh_vs_TauFL","dPhi_SVPV_genTauh_vs_TauFL",64,-3.14159265359/32,3.14159265359/32,50,0.,2,"dPhi_SVPV_genTauh","TauFL");
	dPhi_SVPV_genTauhPlus_vs_TauFL=HConfig.GetTH2D(Name+"_dPhi_SVPV_genTauhPlus_vs_TauFL","dPhi_SVPV_genTauhPlus_vs_TauFL",64,-3.14159265359/32,3.14159265359/32,50,0,2,"dPhi_SVPV_genTauhPlus","TauFL");
	dPhi_SVPV_genTauhMinus_vs_TauFL=HConfig.GetTH2D(Name+"_dPhi_SVPV_genTauhMinus_vs_TauFL","dPhi_SVPV_genTauhMinus_vs_TauFL",64,-3.14159265359/32,3.14159265359/32,50,0,2,"dPhi_SVPV_genTauhMinus","TauFL");

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

	Gen_TauA1_GJ=HConfig.GetTH1D(Name+"_Gen_TauA1_GJ","Gen_TauA1_GJ",100,0,.05,"Gen_TauA1_GJ","Events");
	Gen_TauMu_GJ=HConfig.GetTH1D(Name+"_Gen_TauMu_GJ","Gen_TauMu_GJ",100,0,.05,"Gen_TauMu_GJ","Events");

	Gen_DiTau_dPhi=HConfig.GetTH1D(Name+"_Gen_DiTau_dPhi","Gen_DiTau_dPhi",256,0,2*3.14159265359,"Gen_DiTau_dPhi","Events");
	Gen_DiTau_Pt=HConfig.GetTH1D(Name+"_Gen_DiTau_Pt","Gen_DiTau_Pt",30,0,30,"Gen_DiTau_Pt","Events");
	Gen_Z_Pt=HConfig.GetTH1D(Name+"_Gen_Z_Pt","Gen_Z_Pt",30,0,30,"Gen_Z_Pt","Events");
	Gen_Z_M=HConfig.GetTH1D(Name+"_Gen_Z_M","Gen_Z_M",180,60,150,"Gen_Z_M","Events");
	Gen_DiTau_PtBalance_M=HConfig.GetTH1D(Name+"_Gen_DiTau_PtBalance_M","Gen_DiTau_PtBalance_M",280,20,160,"Gen_DiTau_PtBalance_M","Events");
	Gen_DiTau_PtBalance_dM=HConfig.GetTH1D(Name+"_Gen_DiTau_PtBalance_dM","Gen_DiTau_PtBalance_dM",100,-50,50,"Gen_DiTau_PtBalance_dM","Events");
	Gen_TauMu_PtBalance_Pt=HConfig.GetTH1D(Name+"_Gen_TauMu_PtBalance_Pt","Gen_TauMu_PtBalance_Pt",100,-50,50,"Gen_TauMu_PtBalance_Pt","Events");
	Gen_TauMu_PtBalance_dP=HConfig.GetTH1D(Name+"_Gen_TauMu_PtBalance_dP","Gen_TauMu_PtBalance_dP",100,-50,50,"Gen_TauMu_PtBalance_dP","Events");
	Gen_TauA1_dP=HConfig.GetTH1D(Name+"_Gen_TauA1_dP","Gen_TauA1_dP",50,0,50,"Gen_TauA1_dP","Events");

	Gen_TPTF_TauA1_Solution_NoSelection=HConfig.GetTH1D(Name+"_Gen_TPTF_TauA1_RightSolution_NoSelection","Gen_TPTF_TauA1_RightSolution_NoSelection",3,-0.5,2.5,"Gen_TPTF_TauA1_RightSolution_NoSelection","Events");
	Gen_TPTF_TauA1_Solution_WithSelection=HConfig.GetTH1D(Name+"_Gen_TPTF_TauA1_RightSolution_WithSelection","Gen_TPTF_TauA1_RightSolution_WithSelection",3,-0.5,2.5,"Gen_TPTF_TauA1_RightSolution_WithSelection","Events");

	dP_GenTauMuPtBalance_vs_dPTauh=HConfig.GetTH2D(Name+"_dP_GenTauMuPtBalance_vs_dPTauh","dP_GenTauMuPtBalance_vs_dPTauh",50,-50,50,50,0,50,"p^{'}_{#tau_{#mu}} - p_{#tau_{#mu}} ","2*D");
	Pt_vs_dPhi_DiTauGen=HConfig.GetTH2D(Name+"_Pt_vs_dPhi_DiTauGen","Pt_vs_dPhi_DiTauGen",10,0,30,256,0,2*3.14159265359,"Pt ditau","dPhi(tau,tau)");

	TauFL_WithTauFLSigmaCut=HConfig.GetTH1D(Name+"_TauFL_WithTauFLSigmaCut","TauFL_WithTauFLSigmaCut",100,-2,4,"tau flight length with tauflsigma cut","Events");
	TauFLSigmaCut_vs_Res=HConfig.GetTH2D(Name+"_TauFLSigmaCut_vs_Res","TauFLSigmaCut_vs_Res",80,-10,30,2000,-0.1,0.1,"TauFLSigmaCut","dPhi(SVPV,genTauh)");
	TauFLSigma_vs_Res=HConfig.GetTH2D(Name+"_TauFLSigma_vs_Res","TauFLSigma_vs_Res",80,-10,30,200,-0.1,0.1,"TauFLSigma","dPhi(SVPV,genTauh)");
	TauFLSigma_vs_UnphysicalAll=HConfig.GetTH2D(Name+"_TauFLSigma_vs_UnphysicalAll","TauFLSigma_vs_UnphysicalAll",40,-10,70,2,-0.5,1.5,"TauFLSigma","UnphysicalAll");
	TauFLSigma_vs_UnphysicalProb=HConfig.GetTH1D(Name+"_TauFLSigma_vs_UnphysicalProb","TauFLSigma_vs_UnphysicalProb",40,-10,70,"TauFLSigma","Unphysical Prob");
	TauFL_vs_UnphysicalAll=HConfig.GetTH2D(Name+"_TauFL_vs_UnphysicalAll","TauFL_vs_UnphysicalAll",40,0,2,2,-0.5,1.5,"TauFL_NoTauFLSigmaCut","UnphysicalAll");
	TauFL_vs_UnphysicalProb=HConfig.GetTH1D(Name+"_TauFL_vs_UnphysicalProb","TauFL_vs_UnphysicalProb",40,0,2,"TauFL_NoTauFLSigmaCut","Unphysical Prob");

	//ThreeProngTauFit Studies
	TPTF_TauA1_pRes_BestSolution=HConfig.GetTH1D(Name+"_TPTF_TauA1_pRes_BestSolution","TPTF_TauA1_pRes_BestSolution",100,-50,50,"TPTF_TauA1_pRes_BestSolution","Events");
	TPTF_TauA1_pRes_FitSolution=HConfig.GetTH1D(Name+"_TPTF_TauA1_pRes_FitSolution","TPTF_TauA1_pRes_FitSolution",100,-50,50,"TPTF_TauA1_pRes_FitSolution","Events");
	TPTF_TauA1_BestSolution_vs_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_BestSolution_vs_FitSolution","TPTF_TauA1_BestSolution_vs_FitSolution",3,-0.5,2.5,3,-0.5,2.5,"TPTF_TauA1_BestSolution","FitSolution");
	TPTF_TauA1_RightSolution_vs_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_RightSolution_vs_FitSolution","TPTF_TauA1_RightSolution_vs_FitSolution",3,-0.5,2.5,3,-0.5,2.5,"TPTF_TauA1_RightSolution","FitSolution");
	TPTF_TauA1_pRes_vs_GenA1Mass_BestSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pRes_vs_GenA1Mass_BestSolution","TPTF_TauA1_pRes_vs_GenA1Mass_BestSolution",51,-51,51,7,0.8,1.5,"TPTF_TauA1_pRes_BestSolution","GenA1_Mass");
	TPTF_TauA1_pRes_vs_GenA1Mass_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pRes_vs_GenA1Mass_FitSolution","TPTF_TauA1_pRes_vs_GenA1Mass_FitSolution",51,-51,51,7,0.8,1.5,"TPTF_TauA1_pRes_FitSolution","GenA1_Mass");
	TPTF_TauA1_pRes_vs_GenGJAngle_BestSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pRes_vs_GenGJAngle_BestSolution","TPTF_TauA1_pRes_vs_GenGJAngle_BestSolution",51,-51,51,50,0.0,0.035,"TPTF_TauA1_pRes_BestSolution","GenGJAngle");
	TPTF_TauA1_pRes_vs_GenGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pRes_vs_GenGJAngle_FitSolution","TPTF_TauA1_pRes_vs_GenGJAngle_FitSolution",51,-51,51,50,0.0,0.035,"TPTF_TauA1_pRes_FitSolution","GenGJAngle");
	TPTF_A1_pRes_vs_GenGJAngle=HConfig.GetTH2D(Name+"_TPTF_A1_pRes_vs_GenGJAngle","TPTF_A1_pRes_vs_GenGJAngle",51,-11,11,50,0.0,0.035,"TPTF_A1_pRes","GenGJAngle");

	TPTF_TauA1_pRes_vs_RecoA1Mass_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pRes_vs_RecoA1Mass_FitSolution","TPTF_TauA1_pRes_vs_RecoA1Mass_FitSolution",51,-51,51,7,0.8,1.5,"TPTF_TauA1_pRes_FitSolution","RecoA1_Mass");
	TPTF_TauA1_p_iRes_sq_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_p_iRes_sq_vs_RecoGJAngle_FitSolution","TPTF_TauA1_p_iRes_sq_vs_RecoGJAngle_FitSolution",51,-51,201,50,0.0,0.035,"TPTF_TauA1_p_iRes_sq_FitSolution","RecoGJAngle");
	TPTF_TauA1_pRes_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pRes_vs_RecoGJAngle_FitSolution","TPTF_TauA1_pRes_vs_RecoGJAngle_FitSolution",51,-51,51,50,0.0,0.035,"TPTF_TauA1_pRes_FitSolution","RecoGJAngle");
	TPTF_TauA1_pxRes_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pxRes_vs_RecoGJAngle_FitSolution","TPTF_TauA1_pxRes_vs_RecoGJAngle_FitSolution",51,-51,51,50,0.0,0.035,"TPTF_TauA1_pxRes_FitSolution","RecoGJAngle");
	TPTF_TauA1_pyRes_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pyRes_vs_RecoGJAngle_FitSolution","TPTF_TauA1_pyRes_vs_RecoGJAngle_FitSolution",51,-51,51,50,0.0,0.035,"TPTF_TauA1_pyRes_FitSolution","RecoGJAngle");
	TPTF_TauA1_pzRes_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pzRes_vs_RecoGJAngle_FitSolution","TPTF_TauA1_pzRes_vs_RecoGJAngle_FitSolution",51,-51,51,50,0.0,0.035,"TPTF_TauA1_pzRes_FitSolution","RecoGJAngle");
	TPTF_TauA1_ptRes_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_ptRes_vs_RecoGJAngle_FitSolution","TPTF_TauA1_ptRes_vs_RecoGJAngle_FitSolution",51,-51,51,50,0.0,0.035,"TPTF_TauA1_ptRes_FitSolution","RecoGJAngle");
	TPTF_TauA1_pRes_vs_RecoChi2_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pRes_vs_RecoChi2_FitSolution","TPTF_TauA1_pRes_vs_RecoChi2_FitSolution",51,-51,51,50,0,20,"TPTF_TauA1_pRes_FitSolution","Reco_Chi2");
	TPTF_TauA1_pxsqRes_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pxsqRes_vs_RecoGJAngle_FitSolution","TPTF_TauA1_pxsqRes_vs_RecoGJAngle_FitSolution",51,-201,201,50,0.0,0.035,"TPTF_TauA1_pxsqRes_FitSolution","RecoGJAngle");
	TPTF_TauA1_pysqRes_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pysqRes_vs_RecoGJAngle_FitSolution","TPTF_TauA1_pysqRes_vs_RecoGJAngle_FitSolution",51,-201,201,50,0.0,0.035,"TPTF_TauA1_pysqRes_FitSolution","RecoGJAngle");
	TPTF_TauA1_pzsqRes_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_pzsqRes_vs_RecoGJAngle_FitSolution","TPTF_TauA1_pzsqRes_vs_RecoGJAngle_FitSolution",51,-201,201,50,0.0,0.035,"TPTF_TauA1_pzsqRes_FitSolution","RecoGJAngle");
	TPTF_TauA1_p_orthoRes_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_p_orthoRes_vs_RecoGJAngle_FitSolution","TPTF_TauA1_p_orthoRes_vs_RecoGJAngle_FitSolution",51,-51,51,50,0.0,0.035,"TPTF_TauA1_p_orthoRes_FitSolution","RecoGJAngle");
	TPTF_TauA1_p_paralRes_vs_RecoGJAngle_FitSolution=HConfig.GetTH2D(Name+"_TPTF_TauA1_p_paralRes_vs_RecoGJAngle_FitSolution","TPTF_TauA1_p_paralRes_vs_RecoGJAngle_FitSolution",51, -1, 1,50,0.0,0.035,"TPTF_TauA1_p_paralRes_FitSolution","RecoGJAngle");
	TPTF_A1_pRes_vs_RecoGJAngle=HConfig.GetTH2D(Name+"_TPTF_A1_pRes_vs_RecoGJAngle","TPTF_A1_pRes_vs_RecoGJAngle",51,-11,11,50,0.0,0.035,"TPTF_A1_pRes","RecoGJAngle");

	TPTF_TauA1_p_Reco=HConfig.GetTH1D(Name+"_TPTF_TauA1_p_Reco","TPTF_TauA1_p_Reco",75, 0,150,"TPTF_TauA1_p_Reco","Events");
	TPTF_TauA1_pt_Reco=HConfig.GetTH1D(Name+"_TPTF_TauA1_pt_Reco","TPTF_TPTF_TauA1_pt_Reco",75, 0,150,"TPTF_TauA1_pt_Reco","Events");
	TPTF_TauA1_px_Reco=HConfig.GetTH1D(Name+"_TPTF_TauA1_px_Reco","TPTF_TPTF_TauA1_px_Reco",100,-100,100,"TPTF_TauA1_px_Reco","Events");
	TPTF_TauA1_py_Reco=HConfig.GetTH1D(Name+"_TPTF_TauA1_py_Reco","TPTF_TPTF_TauA1_py_Reco",100,-100,100,"TPTF_TauA1_py_Reco","Events");
	TPTF_TauA1_pz_Reco=HConfig.GetTH1D(Name+"_TPTF_TauA1_pz_Reco","TPTF_TPTF_TauA1_pz_Reco",100,-100,100,"TPTF_TauA1_pz_Reco","Events");
	TPTF_TauA1_p_Gen=HConfig.GetTH1D(Name+"_TPTF_TauA1_p_Gen","TPTF_TauA1_p_Gen",75, 0,150,"TPTF_TauA1_p_Gen","Events");
	TPTF_TauA1_pt_Gen=HConfig.GetTH1D(Name+"_TPTF_TauA1_pt_Gen","TPTF_TauA1_pt_Gen",75, 0,150,"TPTF_TauA1_pt_Gen","Events");
	TPTF_TauA1_px_Gen=HConfig.GetTH1D(Name+"_TPTF_TauA1_px_Gen","TPTF_TauA1_px_Gen",100,-100,100,"TPTF_TauA1_px_Gen","Events");
	TPTF_TauA1_py_Gen=HConfig.GetTH1D(Name+"_TPTF_TauA1_py_Gen","TPTF_TauA1_py_Gen",100,-100,100,"TPTF_TauA1_py_Gen","Events");
	TPTF_TauA1_pz_Gen=HConfig.GetTH1D(Name+"_TPTF_TauA1_pz_Gen","TPTF_TauA1_pz_Gen",100,-100,100,"TPTF_TauA1_pz_Gen","Events");

	TPTF_TauA1_pxsq_Reco=HConfig.GetTH1D(Name+"_TPTF_TauA1_pxsq_Reco","TPTF_TauA1_pxsq_Reco",100,0,4000,"TPTF_TauA1_pxsq_Reco","Events");
	TPTF_TauA1_pysq_Reco=HConfig.GetTH1D(Name+"_TPTF_TauA1_pysq_Reco","TPTF_TauA1_pysq_Reco",100,0,4000,"TPTF_TauA1_pysq_Reco","Events");
	TPTF_TauA1_pzsq_Reco=HConfig.GetTH1D(Name+"_TPTF_TauA1_pzsq_Reco","TPTF_TauA1_pzsq_Reco",100,0,4000,"TPTF_TauA1_pzsq_Reco","Events");
	TPTF_TauA1_pxsq_Gen=HConfig.GetTH1D(Name+"_TPTF_TauA1_pxsq_Gen","TPTF_TauA1_pxsq_Gen",100,0,4000,"TPTF_TauA1_pxsq_Gen","Events");
	TPTF_TauA1_pysq_Gen=HConfig.GetTH1D(Name+"_TPTF_TauA1_pysq_Gen","TPTF_TauA1_pysq_Gen",100,0,4000,"TPTF_TauA1_pysq_Gen","Events");
	TPTF_TauA1_pzsq_Gen=HConfig.GetTH1D(Name+"_TPTF_TauA1_pzsq_Gen","TPTF_TauA1_pzsq_Gen",100,0,4000,"TPTF_TauA1_pzsq_Gen","Events");

	TPTF_TauA1_ptRes_vs_ptGen=HConfig.GetTH2D(Name+"_TPTF_TauA1_ptRes_vs_ptGen","TPTF_TauA1_ptRes_vs_ptGen",100,-50,50,50,0.0,50,"TPTF_TauA1_ptRes","ptGen");
	TPTF_TauA1_ptRes_vs_ptReco=HConfig.GetTH2D(Name+"_TPTF_TauA1_ptRes_vs_ptReco","TPTF_TauA1_ptRes_vs_ptReco",100,-50,50,50,0.0,50,"TPTF_TauA1_ptRes","ptReco");

	TPTF_Neutrino_UnFitTracks_Mass=HConfig.GetTH1D(Name+"_TPTF_Neutrino_UnFitTracks_Mass","TPTF_Neutrino_UnFitTracks_Mass",50,-1,1,"TPTF_Neutrino_UnFitTracks_Mass","Events");
	TPTF_Neutrino_UnFitTracks_Mass_vs_TauFL=HConfig.GetTH2D(Name+"_TPTF_Neutrino_UnFitTracks_Mass_vs_TauFL","TPTF_Neutrino_UnFitTracks_Mass_vs_TauFL",50,-1,1,50,0.0,2,"TPTF_Neutrino_UnFitTracks_Mass","TauFL");
	TPTF_Neutrino_ReFitTracks_Mass=HConfig.GetTH1D(Name+"_TPTF_Neutrino_ReFitTracks_Mass","TPTF_Neutrino_ReFitTracks_Mass",50,-1,1,"TPTF_Neutrino_ReFitTracks_Mass","Events");
	TPTF_Neutrino_ReFitTracks_Mass_vs_TauFL=HConfig.GetTH2D(Name+"_TPTF_Neutrino_ReFitTracks_Mass_vs_TauFL","TPTF_Neutrino_ReFitTracks_Mass_vs_TauFL",50,-1,1,50,0.0,2,"TPTF_Neutrino_ReFitTracks_Mass","TauFL");
	TPTF_Neutrino_PFTau_Mass=HConfig.GetTH1D(Name+"_TPTF_Neutrino_PFTau_Mass","TPTF_Neutrino_PFTau_Mass",50,-1,1,"TPTF_Neutrino_PFTau_Mass","Events");
	TPTF_Neutrino_PFTau_Mass_vs_TauFL=HConfig.GetTH2D(Name+"_TPTF_Neutrino_PFTau_Mass_vs_TauFL","TPTF_Neutrino_PFTau_Mass_vs_TauFL",50,-1,1,50,0.0,2,"TPTF_Neutrino_PFTau_Mass","TauFL");

	TransTrk_Failure_withSelection=HConfig.GetTH1D(Name+"_TransTrk_Failure_withSelection","TransTrk_Failure_withSelection",2,-0.5,1.5,"TransTrk_Failure_withSelection","Events");
	TransTrk_Failure_noSelection=HConfig.GetTH1D(Name+"_TransTrk_Failure_noSelection","TransTrk_Failure_noSelection",2,-0.5,1.5,"TransTrk_Failure_noSelection","Events");

	//DiTau Reco
	Reco_ZMass=HConfig.GetTH1D(Name+"_Reco_ZMass","Reco_ZMass",180,60,150,"Reco_ZMass","Events");
	Reco_ZMass_UnboostedGenZ=HConfig.GetTH1D(Name+"_Reco_ZMass_UnboostedGenZ","Reco_ZMass_UnboostedGenZ",180,60,150,"Reco_ZMass_UnboostedGenZ","Events");
	Reco_EventFit_Solution=HConfig.GetTH1D(Name+"_Reco_EventFit_Solution","Reco_EventFit_Solution",5,-2.5,2.5,"Reco_EventFit_Solution","Events");
	Reco_A1Fit_Solution=HConfig.GetTH1D(Name+"_Reco_A1Fit_Solution","Reco_A1Fit_Solution",5,-2.5,2.5,"Reco_A1Fit_Solution","Events");
	Reco_Chi2=HConfig.GetTH1D(Name+"_Reco_Chi2","Reco_Chi2",30,0,30,"Reco_Chi2","Events");
	Reco_Chi2_FitSolutionOnly=HConfig.GetTH1D(Name+"_Reco_Chi2_FitSolutionOnly","Reco_Chi2_FitSolutionOnly",30,0,30,"Reco_Chi2_FitSolutionOnly","Events");
	Reco_Chi2_FitSolutionOnlyLargeScale=HConfig.GetTH1D(Name+"_Reco_Chi2_FitSolutionOnlyLargeScale","Reco_Chi2_FitSolutionOnlyLargeScale",100,0,500000,"Reco_Chi2_FitSolutionOnlyLargeScale","Events");
	Reco_ConstrainedDeltaSum=HConfig.GetTH1D(Name+"_Reco_ConstrainedDeltaSum","Reco_ConstrainedDeltaSum",50,0,0.1,"Reco_ConstrainedDeltaSum","Events");
	Reco_ConstrainedDeltaMass=HConfig.GetTH1D(Name+"_Reco_ConstrainedDeltaMass","Reco_ConstrainedDeltaMass",100,0,40,"Reco_ConstrainedDeltaMass","Events");
	Reco_ConstrainedDeltaPt=HConfig.GetTH1D(Name+"_Reco_ConstrainedDeltaPt","Reco_ConstrainedDeltaPt",100,0,1,"Reco_ConstrainedDeltaPt","Events");
	Reco_NIter=HConfig.GetTH1D(Name+"_Reco_NIter","Reco_NIter",100,0,100,"Reco_NIter","Events");
	Reco_TauMu_DeltaPX_FitImpact=HConfig.GetTH1D(Name+"_Reco_TauMu_DeltaPX_FitImpact","Reco_TauMu_DeltaPX_FitImpact",100,-30,30,"Reco_TauMu_DeltaPX_FitImpact","Events");
	Reco_TauMu_DeltaPY_FitImpact=HConfig.GetTH1D(Name+"_Reco_TauMu_DeltaPY_FitImpact","Reco_TauMu_DeltaPY_FitImpact",100,-30,30,"Reco_TauMu_DeltaPY_FitImpact","Events");
	Reco_TauMu_DeltaPZ_FitImpact=HConfig.GetTH1D(Name+"_Reco_TauMu_DeltaPZ_FitImpact","Reco_TauMu_DeltaPZ_FitImpact",100,-30,30,"Reco_TauMu_DeltaPZ_FitImpact","Events");
	Reco_TauA1_DeltaPX_FitImpact=HConfig.GetTH1D(Name+"_Reco_TauA1_DeltaPX_FitImpact","Reco_TauA1_DeltaPX_FitImpact",100,-30,30,"Reco_TauA1_DeltaPX_FitImpact","Events");
	Reco_TauA1_DeltaPY_FitImpact=HConfig.GetTH1D(Name+"_Reco_TauA1_DeltaPY_FitImpact","Reco_TauA1_DeltaPY_FitImpact",100,-30,30,"Reco_TauA1_DeltaPY_FitImpact","Events");
	Reco_TauA1_DeltaPZ_FitImpact=HConfig.GetTH1D(Name+"_Reco_TauA1_DeltaPZ_FitImpact","Reco_TauA1_DeltaPZ_FitImpact",100,-30,30,"Reco_TauA1_DeltaPZ_FitImpact","Events");
	Reco_Z_Energy_Res=HConfig.GetTH1D(Name+"_Reco_Z_Energy_Res","Reco_Z_Energy_Res",100,-100,100,"Reco_Z_Energy_Res","Events");

	Reco_PtRes_TauA1=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1","Reco_PtRes_TauA1",100,-50,50,"Reco_PtRes_TauA1","Events");
	Reco_PtRes_TauA1_AmbPoint0=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_AmbPoint0","Reco_PtRes_TauA1_AmbPoint0",100,-50,50,"Reco_PtRes_TauA1_AmbPoint0","Events");
	Reco_PtRes_TauA1_AmbPoint12=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_AmbPoint12","Reco_PtRes_TauA1_AmbPoint12",100,-50,50,"Reco_PtRes_TauA1_AmbPoint12","Events");
	Reco_PtRes_TauA1_AmbPoint1=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_AmbPoint1","Reco_PtRes_TauA1_AmbPoint1",100,-50,50,"Reco_PtRes_TauA1_AmbPoint1","Events");

	Reco_PtRes_TauMu=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu","Reco_PtRes_TauMu",100,-50,50,"Reco_PtRes_TauMu","Events");
	Reco_PtRes_TauMu_AmbPoint0=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_AmbPoint0","Reco_PtRes_TauMu_AmbPoint0",100,-50,50,"Reco_PtRes_TauMu_AmbPoint0","Events");
	Reco_PtRes_TauMu_AmbPoint12=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_AmbPoint12","Reco_PtRes_TauMu_AmbPoint12",100,-50,50,"Reco_PtRes_TauMu_AmbPoint12","Events");
	Reco_PtRes_TauMu_AmbPoint1=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_AmbPoint1","Reco_PtRes_TauMu_AmbPoint1",100,-50,50,"Reco_PtRes_TauMu_AmbPoint1","Events");

	Reco_dPhi_TauMuTauA1_AfterFit=HConfig.GetTH1D(Name+"_Reco_dPhi_TauMuTauA1_AfterFit","Reco_dPhi_TauMuTauA1_AfterFit",100,-7,7,"Reco_dPhi_TauMuTauA1_AfterFit","Events");
	Reco_dPhi_TauMuTauA1_BeforeFit=HConfig.GetTH1D(Name+"_Reco_dPhi_TauMuTauA1_BeforeFit","Reco_dPhi_TauMuTauA1_BeforeFit",100,-7,7,"Reco_dPhi_TauMuTauA1_BeforeFit","Events");

	Reco_PtRes_TauA1_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_NoFit","Reco_PtRes_TauA1_NoFit",100,-50,50,"Reco_PtRes_TauA1_NoFit","Events");
	Reco_PtRes_TauA1_AmbPoint0_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_AmbPoint0_NoFit","Reco_PtRes_TauA1_AmbPoint0_NoFit",100,-50,50,"Reco_PtRes_TauA1_AmbPoint0_NoFit","Events");
	Reco_PtRes_TauA1_AmbPoint12_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_AmbPoint12_NoFit","Reco_PtRes_TauA1_AmbPoint12_NoFit",100,-50,50,"Reco_PtRes_TauA1_AmbPoint12_NoFit","Events");

	Reco_PtRes_TauMu_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_NoFit","Reco_PtRes_TauMu_NoFit",100,-50,50,"Reco_PtRes_TauMu_NoFit","Events");
	Reco_PtRes_TauMu_AmbPoint0_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_AmbPoint0_NoFit","Reco_PtRes_TauMu_AmbPoint0_NoFit",100,-50,50,"Reco_PtRes_TauMu_AmbPoint0_NoFit","Events");
	Reco_PtRes_TauMu_AmbPoint12_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_AmbPoint12_NoFit","Reco_PtRes_TauMu_AmbPoint12_NoFit",100,-50,50,"Reco_PtRes_TauMu_AmbPoint12_NoFit","Events");

	Reco_ZMass_MassScan=HConfig.GetTH1D(Name+"_Reco_ZMass_MassScan","Reco_ZMass_MassScan",50,0,250,"Reco_ZMass_MassScan","Events");
	Reco_ZMass_MassScanUnboosted=HConfig.GetTH1D(Name+"_Reco_ZMass_MassScanUnboosted","Reco_ZMass_MassScanUnboosted",50,0,250,"Reco_ZMass_MassScanUnboosted","Events");
	Reco_ZMasswithProbWeight_MassScan=HConfig.GetTH1D(Name+"_Reco_ZMasswithProbWeight_MassScan","Reco_ZMasswithProbWeight_MassScan",50,0,250,"Reco_ZMasswithProbWeight_MassScan","Events");
	Reco_ProbStack_MassScan=HConfig.GetTH1D(Name+"_Reco_ProbStack_MassScan","Reco_ProbStack_MassScan",50,0,1,"Reco_ProbStack_MassScan","Events");
	Reco_ZMass_PDF=HConfig.GetTH1D(Name+"_Reco_ZMass_PDF","Reco_ZMass_PDF",50,0,250,"Reco_ZMass_PDF","Events");
	GenZ_Pt_Unboosted=HConfig.GetTH1D(Name+"_GenZ_Pt_Unboosted","GenZ_Pt_Unboosted",40,0,40,"GenZ_Pt_Unboosted","Events");
	RecoZ_Pt=HConfig.GetTH1D(Name+"_RecoZ_Pt","RecoZ_Pt",40,0,40,"RecoZ_Pt","Events");
	RecoZ_Pt_Unboosted=HConfig.GetTH1D(Name+"_RecoZ_Pt_Unboosted","RecoZ_Pt_Unboosted",40,0,40,"RecoZ_Pt_Unboosted","Events");

	Reco_TauMu_ResCosTheta=HConfig.GetTH1D(Name+"_RecoTaumu_ResCosTheta","RecoTaumu_ResCosTheta",31,-2,2,"RecoTaumu_ResCosTheta","Events");
	Reco_TauMu_ResPhi=HConfig.GetTH1D(Name+"_Reco_Taumu_ResPhi","Reco_Taumu_ResPhi",31,-3.14159265359*2,3.14159265359*2,"Reco_Taumu_ResPhi","Events");

	//DiTau Reco with generator information/particles
	GenReco_ZMass=HConfig.GetTH1D(Name+"_GenReco_ZMass","GenReco_ZMass",180,60,150,"GenReco_ZMass","Events");
	GenReco_EventFit_Solution=HConfig.GetTH1D(Name+"_GenReco_EventFit_Solution","GenReco_EventFit_Solution",5,-2.5,2.5,"GenReco_EventFit_Solution","Events");
	GenReco_A1Fit_Solution=HConfig.GetTH1D(Name+"_GenReco_A1Fit_Solution","GenReco_A1Fit_Solution",5,-2.5,2.5,"GenReco_A1Fit_Solution","Events");
	GenReco_Chi2=HConfig.GetTH1D(Name+"_GenReco_Chi2","GenReco_Chi2",30,0,30,"GenReco_Chi2","Events");
	GenReco_Chi2_FitSolutionOnly=HConfig.GetTH1D(Name+"_GenReco_Chi2_FitSolutionOnly","GenReco_Chi2_FitSolutionOnly",30,0,30,"GenReco_Chi2_FitSolutionOnly","Events");
	GenReco_ConstrainedDeltaSum=HConfig.GetTH1D(Name+"_GenReco_ConstrainedDeltaSum","GenReco_ConstrainedDeltaSum",50,0,0.1,"GenReco_ConstrainedDeltaSum","Events");
	GenReco_NIter=HConfig.GetTH1D(Name+"_GenReco_NIter","GenReco_NIter",100,0,100,"GenReco_NIter","Events");

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
	Extradist1d.push_back(&Tau_Mass_Inclusive);
	Extradist1d.push_back(&Tau_Mass_sq_Inclusive);
	Extradist1d.push_back(&Tau_Mass_Inclusive_NoTLV);
	Extradist1d.push_back(&Tau_Mass_Inclusive_UnFitTracks);
	Extradist1d.push_back(&Tau_Mass_Inclusive_ReFitTracks);
	Extradist1d.push_back(&Tau_Mass_Difference_PFTau_UnFitTracks_3PS);
	Extradist1d.push_back(&MET_phi);
	Extradist1d.push_back(&TauFL_NoTauFLSigmaCut);
	Extradist1d.push_back(&TauFLSigned_NoTauFLSigmaCut);
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
	Extradist2d.push_back(&dPhi_SVPV_genTauh_vs_TauFL);
	Extradist2d.push_back(&dPhi_SVPV_genTauhPlus_vs_TauFL);
	Extradist2d.push_back(&dPhi_SVPV_genTauhMinus_vs_TauFL);
	Extradist1d.push_back(&dTheta_SVPV_genTauh);
	Extradist1d.push_back(&Angle_SVPV_genTauh);
	Extradist1d.push_back(&Phi_POCAPV);
	Extradist1d.push_back(&Phi_genTaumu);
	Extradist1d.push_back(&Theta_POCAPV);
	Extradist1d.push_back(&Theta_genTaumu);
	Extradist1d.push_back(&dPhi_POCAPV_genTaumu);
	Extradist1d.push_back(&dTheta_POCAPV_genTaumu);
	Extradist1d.push_back(&Gen_TauA1_GJ);
	Extradist1d.push_back(&Gen_TauMu_GJ);
	Extradist1d.push_back(&Gen_DiTau_dPhi);
	Extradist1d.push_back(&Gen_DiTau_Pt);
	Extradist1d.push_back(&Gen_Z_Pt);
	Extradist1d.push_back(&Gen_Z_M);
	Extradist1d.push_back(&Gen_DiTau_PtBalance_M);
	Extradist1d.push_back(&Gen_DiTau_PtBalance_dM);
	Extradist1d.push_back(&Gen_TauMu_PtBalance_Pt);
	Extradist1d.push_back(&Gen_TauMu_PtBalance_dP);
	Extradist1d.push_back(&Gen_TauA1_dP);
	Extradist1d.push_back(&Gen_TPTF_TauA1_Solution_NoSelection);
	Extradist1d.push_back(&Gen_TPTF_TauA1_Solution_WithSelection);
	Extradist1d.push_back(&dPhi_MinusSVPV_genTaumu);
	Extradist1d.push_back(&dTheta_MinusSVPV_genTaumu);
	Extradist1d.push_back(&Angle_MinusSVPV_genTaumu);

	Extradist1d.push_back(&TPTF_TauA1_pRes_BestSolution);
	Extradist1d.push_back(&TPTF_TauA1_pRes_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_BestSolution_vs_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_RightSolution_vs_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pRes_vs_GenA1Mass_BestSolution);
	Extradist2d.push_back(&TPTF_TauA1_pRes_vs_GenA1Mass_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pRes_vs_GenGJAngle_BestSolution);
	Extradist2d.push_back(&TPTF_TauA1_pRes_vs_GenGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_A1_pRes_vs_GenGJAngle);
	Extradist2d.push_back(&TPTF_TauA1_pRes_vs_RecoChi2_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pRes_vs_RecoA1Mass_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pRes_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_p_iRes_sq_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pxRes_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pyRes_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pzRes_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_ptRes_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pRes_vs_RecoChi2_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pxsqRes_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pysqRes_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_pzsqRes_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_p_orthoRes_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_TauA1_p_paralRes_vs_RecoGJAngle_FitSolution);
	Extradist2d.push_back(&TPTF_A1_pRes_vs_RecoGJAngle);
	Extradist2d.push_back(&TPTF_TauA1_ptRes_vs_ptGen);
	Extradist2d.push_back(&TPTF_TauA1_ptRes_vs_ptReco);

	Extradist1d.push_back(&TPTF_TauA1_p_Reco);
	Extradist1d.push_back(&TPTF_TauA1_pt_Reco);
	Extradist1d.push_back(&TPTF_TauA1_px_Reco);
	Extradist1d.push_back(&TPTF_TauA1_py_Reco);
	Extradist1d.push_back(&TPTF_TauA1_pz_Reco);
	Extradist1d.push_back(&TPTF_TauA1_p_Gen);
	Extradist1d.push_back(&TPTF_TauA1_pt_Gen);
	Extradist1d.push_back(&TPTF_TauA1_px_Gen);
	Extradist1d.push_back(&TPTF_TauA1_py_Gen);
	Extradist1d.push_back(&TPTF_TauA1_pz_Gen);

	Extradist1d.push_back(&TPTF_TauA1_pxsq_Reco);
	Extradist1d.push_back(&TPTF_TauA1_pysq_Reco);
	Extradist1d.push_back(&TPTF_TauA1_pzsq_Reco);
	Extradist1d.push_back(&TPTF_TauA1_pxsq_Gen);
	Extradist1d.push_back(&TPTF_TauA1_pysq_Gen);
	Extradist1d.push_back(&TPTF_TauA1_pzsq_Gen);

	Extradist1d.push_back(&TPTF_Neutrino_UnFitTracks_Mass);
	Extradist2d.push_back(&TPTF_Neutrino_UnFitTracks_Mass_vs_TauFL);
	Extradist1d.push_back(&TPTF_Neutrino_ReFitTracks_Mass);
	Extradist2d.push_back(&TPTF_Neutrino_ReFitTracks_Mass_vs_TauFL);
	Extradist1d.push_back(&TPTF_Neutrino_PFTau_Mass);
	Extradist2d.push_back(&TPTF_Neutrino_PFTau_Mass_vs_TauFL);

	Extradist1d.push_back(&TransTrk_Failure_withSelection);
	Extradist1d.push_back(&TransTrk_Failure_noSelection);

	Extradist1d.push_back(&Reco_ZMass);
	Extradist1d.push_back(&Reco_ZMass_UnboostedGenZ);
	Extradist1d.push_back(&Reco_ZMass_MassScan);
	Extradist1d.push_back(&Reco_ZMass_MassScanUnboosted);
	Extradist1d.push_back(&Reco_ZMasswithProbWeight_MassScan);
	Extradist1d.push_back(&Reco_ProbStack_MassScan);
	Extradist1d.push_back(&Reco_EventFit_Solution);
	Extradist1d.push_back(&Reco_A1Fit_Solution);
	Extradist1d.push_back(&Reco_Chi2);
	Extradist1d.push_back(&Reco_Chi2_FitSolutionOnly);
	Extradist1d.push_back(&Reco_Chi2_FitSolutionOnlyLargeScale);
	Extradist1d.push_back(&Reco_ConstrainedDeltaSum);
	Extradist1d.push_back(&Reco_ConstrainedDeltaMass);
	Extradist1d.push_back(&Reco_ConstrainedDeltaPt);
	Extradist1d.push_back(&Reco_NIter);

	Extradist1d.push_back(&Reco_PtRes_TauA1);
	Extradist1d.push_back(&Reco_PtRes_TauA1_AmbPoint0);
	Extradist1d.push_back(&Reco_PtRes_TauA1_AmbPoint12);
	Extradist1d.push_back(&Reco_PtRes_TauA1_AmbPoint1);
	Extradist1d.push_back(&Reco_PtRes_TauMu);
	Extradist1d.push_back(&Reco_PtRes_TauMu_AmbPoint0);
	Extradist1d.push_back(&Reco_PtRes_TauMu_AmbPoint12);
	Extradist1d.push_back(&Reco_PtRes_TauMu_AmbPoint1);

	Extradist1d.push_back(&Reco_dPhi_TauMuTauA1_AfterFit);
	Extradist1d.push_back(&Reco_dPhi_TauMuTauA1_BeforeFit);

	Extradist1d.push_back(&Reco_PtRes_TauA1_NoFit);
	Extradist1d.push_back(&Reco_PtRes_TauA1_AmbPoint0_NoFit);
	Extradist1d.push_back(&Reco_PtRes_TauA1_AmbPoint12_NoFit);

	Extradist1d.push_back(&Reco_PtRes_TauMu_NoFit);
	Extradist1d.push_back(&Reco_PtRes_TauMu_AmbPoint0_NoFit);
	Extradist1d.push_back(&Reco_PtRes_TauMu_AmbPoint12_NoFit);

	Extradist1d.push_back(&Reco_TauMu_DeltaPX_FitImpact);
	Extradist1d.push_back(&Reco_TauMu_DeltaPY_FitImpact);
	Extradist1d.push_back(&Reco_TauMu_DeltaPZ_FitImpact);
	Extradist1d.push_back(&Reco_TauA1_DeltaPX_FitImpact);
	Extradist1d.push_back(&Reco_TauA1_DeltaPY_FitImpact);
	Extradist1d.push_back(&Reco_TauA1_DeltaPZ_FitImpact);
	Extradist1d.push_back(&GenZ_Pt_Unboosted);
	Extradist1d.push_back(&RecoZ_Pt_Unboosted);
	Extradist1d.push_back(&Reco_TauMu_ResCosTheta);
	Extradist1d.push_back(&Reco_TauMu_ResPhi);

	Extradist1d.push_back(&Reco_ZMass_PDF);
	Extradist1d.push_back(&Reco_Z_Energy_Res);
	Extradist1d.push_back(&RecoZ_Pt);

	Extradist1d.push_back(&GenReco_ZMass);
	Extradist1d.push_back(&GenReco_EventFit_Solution);
	Extradist1d.push_back(&GenReco_A1Fit_Solution);
	Extradist1d.push_back(&GenReco_Chi2);
	Extradist1d.push_back(&GenReco_Chi2_FitSolutionOnly);
	Extradist1d.push_back(&GenReco_ConstrainedDeltaSum);
	Extradist1d.push_back(&GenReco_NIter);

	Extradist1d.push_back(&NQCD);
	Extradist1d.push_back(&QCD_MT_MuMET_A);
	Extradist1d.push_back(&QCD_MT_MuMET_B);
	Extradist1d.push_back(&QCD_MT_MuMET_C);
	Extradist1d.push_back(&QCD_MT_MuMET_D);
	Extradist1d.push_back(&QCD_MET_A);
	Extradist1d.push_back(&QCD_MET_B);
	Extradist1d.push_back(&QCD_MET_C);
	Extradist1d.push_back(&QCD_MET_D);

	Extradist2d.push_back(&dP_GenTauMuPtBalance_vs_dPTauh);
	Extradist2d.push_back(&Pt_vs_dPhi_DiTauGen);
	Extradist1d.push_back(&TauFL_WithTauFLSigmaCut);
	Extradist2d.push_back(&TauFLSigmaCut_vs_Res);
	Extradist2d.push_back(&TauFLSigma_vs_Res);
	Extradist2d.push_back(&TauFLSigma_vs_UnphysicalAll);
	Extradist1d.push_back(&TauFLSigma_vs_UnphysicalProb);
	Extradist2d.push_back(&TauFL_vs_UnphysicalAll);
	Extradist1d.push_back(&TauFL_vs_UnphysicalProb);

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
	Extradist1d_OS.push_back(&Tau_Mass_Inclusive);
	Extradist1d_OS.push_back(&Tau_Mass_sq_Inclusive);
	Extradist1d_OS.push_back(&Tau_Mass_Inclusive_NoTLV);
	Extradist1d_OS.push_back(&Tau_Mass_Inclusive_UnFitTracks);
	Extradist1d_OS.push_back(&Tau_Mass_Inclusive_ReFitTracks);
	Extradist1d_OS.push_back(&MET_phi);
	Extradist1d_OS.push_back(&TauFL_NoTauFLSigmaCut);
	Extradist1d_OS.push_back(&TauFLSigned_NoTauFLSigmaCut);
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
	Extradist1d_OS.push_back(&TPTF_TauA1_pRes_BestSolution);
	Extradist1d_OS.push_back(&TPTF_TauA1_pRes_FitSolution);
	Extradist1d_OS.push_back(&Reco_ZMass);
	Extradist1d_OS.push_back(&Reco_ZMass_UnboostedGenZ);
	Extradist1d_OS.push_back(&Reco_ZMass_MassScan);
	Extradist1d_OS.push_back(&Reco_ZMasswithProbWeight_MassScan);
	Extradist1d_OS.push_back(&Reco_ProbStack_MassScan);
	Extradist1d_OS.push_back(&Reco_EventFit_Solution);
	Extradist1d_OS.push_back(&Reco_A1Fit_Solution);
	Extradist1d_OS.push_back(&Reco_Chi2);
	Extradist1d_OS.push_back(&Reco_Chi2_FitSolutionOnly);
	Extradist1d_OS.push_back(&Reco_Chi2_FitSolutionOnlyLargeScale);
	Extradist1d_OS.push_back(&Reco_ConstrainedDeltaSum);
	Extradist1d_OS.push_back(&Reco_ConstrainedDeltaMass);
	Extradist1d_OS.push_back(&Reco_ConstrainedDeltaPt);
	Extradist1d_OS.push_back(&Reco_NIter);
	Extradist1d_OS.push_back(&Reco_PtRes_TauA1);
	Extradist1d_OS.push_back(&Reco_PtRes_TauA1_AmbPoint0);
	Extradist1d_OS.push_back(&Reco_PtRes_TauA1_AmbPoint12);
	Extradist1d_OS.push_back(&Reco_PtRes_TauMu);
	Extradist1d_OS.push_back(&Reco_PtRes_TauMu_AmbPoint0);
	Extradist1d_OS.push_back(&Reco_PtRes_TauMu_AmbPoint12);
	Extradist1d_OS.push_back(&Reco_TauMu_DeltaPX_FitImpact);
	Extradist1d_OS.push_back(&Reco_TauMu_DeltaPY_FitImpact);
	Extradist1d_OS.push_back(&Reco_TauMu_DeltaPZ_FitImpact);
	Extradist1d_OS.push_back(&GenZ_Pt_Unboosted);
	Extradist1d_OS.push_back(&RecoZ_Pt_Unboosted);
	Extradist1d_OS.push_back(&Reco_ZMass_PDF);
	Extradist1d_OS.push_back(&TauFL_WithTauFLSigmaCut);
}

void  ZToTaumuTauh::doEvent(){
	unsigned int t;
	int id(Ntp->GetMCID());
	if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}

	static const double TPTF_TauA1_pBiasArr[] = {
		-7.47058823529 ,
		-7.79166666667 ,
		-4.16058394161 ,
		-3.8875 ,
		-1.75294117647 ,
		-1.5 ,
		-1.81707317073 ,
		-1.65217391304 ,
		-0.333333333333 ,
		-0.542857142857 ,
		-1.2972972973 ,
		-1.48275862069 ,
		-1.16 ,
		-3.21739130435 ,
		-2.48979591837 ,
		-1.74468085106 ,
		-2.22222222222 ,
		-2.66666666667 ,
		-4.85714285714 ,
		-0.153846153846 ,
		-10.0 ,
		3.0 ,
		0.0 ,
		-6.0,
		0.0
	};
	std::vector<double> TPTF_TauA1_pBias (TPTF_TauA1_pBiasArr, TPTF_TauA1_pBiasArr + sizeof(TPTF_TauA1_pBiasArr) / sizeof(TPTF_TauA1_pBiasArr[0]) );

	int selVertex(selVertexDummy);
	int selMuon_Iso(selMuonDummy);
	int selMuon_AntiIso(selMuonDummy);
	int selMuon(selMuonDummy);
	int selTau(selTauDummy);
	Charge = ChargeSumDummy;

	if(Ntp->GetMCID() == DataMCType::DY_tautau || (Ntp->GetMCID()>=10 && Ntp->GetMCID()<= 13)) Ntp->SetTauCorrections("scalecorr");

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
		if(selectMuon_Isolation(*it_mu)){
			if(selMuon_Iso == selMuonDummy) selMuon_Iso = *it_mu;
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
		if(selectMuon_AntiIsolation(*it_mu)){
			if(selMuon_AntiIso == selMuonDummy && selMuon_Iso == selMuonDummy) selMuon_AntiIso = *it_mu;
			selectedMuonsAntiIso.push_back(*it_mu);
		}
	}
	if(selMuon_Iso != selMuonDummy && selMuon_AntiIso != selMuonDummy){
		std::cout << "CRITICAL: SELECTED MUON PASSED ISOLATION AND ANTI-ISOLATION CUT --> FIX" << std::endl;
		return;
	}
	if(id == DataMCType::Data && selMuon_AntiIso != selMuonDummy) selMuon = selMuon_AntiIso;
	else selMuon = selMuon_Iso;

	// Tau cuts
	if(verbose) std::cout << "Cut on TauID" << std::endl;
	std::vector<int> selectedTausId;
	selectedTausId.clear();
	for(unsigned i_tau=0; i_tau < Ntp->NPFTaus(); i_tau++){
		if (selectPFTau_Id(i_tau,selectedMuonsId)){
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
	if(selTau != selTauDummy && selMuon != selMuonDummy){
		Charge = Ntp->Muon_Charge(selMuon) + Ntp->PFTau_Charge(selTau);
		value.at(ChargeSum) = Charge;
	}
	else{
		Charge = ChargeSumDummy;
	}
	pass.at(ChargeSum)=(value.at(ChargeSum)==cut.at(ChargeSum));
	if(verbose){
		std::cout << "value at ChargeSum: " <<value.at(ChargeSum) << std::endl;
		std::cout << "pass at ChargeSum: " <<pass.at(ChargeSum) << std::endl;
	}

	// Tau Decay Mode
	if(verbose) std::cout << "Cut on Tau Decay Mode" << std::endl;
	if(selTau != selTauDummy){
		value.at(TauDecayMode) = Ntp->PFTau_hpsDecayMode(selTau);
	}
	pass.at(TauDecayMode)= (value.at(TauDecayMode)>=cut.at(TauDecayMode));
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
			if(Ntp->PFTau_p4(selTau).Vect().Dot(Ntp->PFTau_FlightLength3d(selTau)) < 0){
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

	// MT calculation
	if(verbose) std::cout << "Calculation and Cut on MT distribution" << std::endl;
	double pT,phi,eTmiss,eTmPhi;
	double MT_TauMET;

	if(selMuon == selMuonDummy){
		value.at(MT_MuMET) = MTDummy;
		if(verbose) std::cout << "No Muon selected: neither isolated or anti isolated" << std::endl;
	}
	else if(selMuon_Iso != selMuonDummy && selMuon_AntiIso != selMuonDummy){
		value.at(MT_MuMET) = MTDummy;
		std::cout << "CRITICAL: SELECTED MUON PASSED ISOLATION AND ANTI-ISOLATION CUT --> FIX" << std::endl;
	}
	else if(selMuon != selMuonDummy){
		eTmiss					= Ntp->MET_CorrMVAMuTau_et();
		eTmPhi					= Ntp->MET_CorrMVAMuTau_phi();
		pT						= Ntp->Muon_p4(selMuon).Pt();
		phi						= Ntp->Muon_p4(selMuon).Phi();
		value.at(MT_MuMET)		= Ntp->transverseMass(pT,phi,eTmiss,eTmPhi);
	}
	if(value.at(MT_MuMET) != MTDummy) pass.at(MT_MuMET)=(value.at(MT_MuMET)<cut.at(MT_MuMET));
	if(verbose){
		std::cout << "value at MT_MuMET: " <<value.at(MT_MuMET) << std::endl;
		std::cout << "pass at MT_MuMET: " <<pass.at(MT_MuMET) << std::endl;
	}

	//////////////////////////////////////////////////////////////////////
	//*************************END OF SELECTION*************************//
	//////////////////////////////////////////////////////////////////////


	if(selTau == selTauDummy) MT_TauMET = MTDummy;
	else{
		eTmiss			= Ntp->MET_CorrMVAMuTau_et();
		eTmPhi			= Ntp->MET_CorrMVAMuTau_phi();
		pT				= Ntp->PFTau_p4(selTau).Pt();
		phi				= Ntp->PFTau_p4(selTau).Phi();
		MT_TauMET		= Ntp->transverseMass(pT,phi,eTmiss,eTmPhi);
	}

	// Mvis
	if(verbose) std::cout << "Calculation of Mvis" << std::endl;
	double Mvis;

	if(selTau != selTauDummy){
		if(selMuon_Iso != selMuonDummy && selMuon_AntiIso == selMuonDummy){
			Mvis = (Ntp->PFTau_p4(selTau) + Ntp->Muon_p4(selMuon_Iso)).M();
		}
		else if(selMuon_Iso == selMuonDummy && selMuon_AntiIso != selMuonDummy){
			Mvis = (Ntp->PFTau_p4(selTau) + Ntp->Muon_p4(selMuon_AntiIso)).M();
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
			if(selMuon_Iso != selMuonDummy){
				w *= RSF->HiggsTauTau_MuTau_Trigger_Tau_ScaleMCtoData(Ntp->Muon_p4(selMuon_Iso));
				w *= RSF->HiggsTauTau_MuTau_Iso_Mu(Ntp->Muon_p4(selMuon_Iso));
			}
			if(selTau != selTauDummy){
				w *= RSF->HiggsTauTau_MuTau_Trigger_Mu_ScaleMCtoData(Ntp->PFTau_p4(selTau));
				if(Ntp->PFTau_hpsDecayMode(selTau) == 0){
					w *= OneProngNoPiWeight;
				}
			}
		}
		if(id == DataMCType::DY_mutau_embedded){
			w *= Ntp->Embedding_TauSpinnerWeight();
			w *= Ntp->Embedding_MinVisPtFilter();
			//w *= Ntp->Embedding_SelEffWeight(); not used by Bastian
			if(selTau != selTauDummy) w *= RSF->HiggsTauTau_MuTau_Trigger_Mu_Eff_Data(Ntp->PFTau_p4(selTau));
			if(selMuon_Iso != selMuonDummy) w *= RSF->HiggsTauTau_MuTau_Trigger_Tau_Eff_Data(Ntp->Muon_p4(selMuon_Iso));
		}
		if(selMuon_Iso != selMuonDummy){
			w *= RSF->HiggsTauTau_MuTau_Id_Mu(Ntp->Muon_p4(selMuon_Iso));
		}
	}
	else{w=1;}

	// W+Jets BG Method
	if(verbose) std::cout << "W+Jets BG Method" << std::endl;
	std::vector<unsigned int> exclude_cuts;
	exclude_cuts.clear();
	exclude_cuts.push_back(ChargeSum);
	exclude_cuts.push_back(MT_MuMET);
	//exclude_cuts.push_back(TauDecayMode);
	//exclude_cuts.push_back(TauFLSigma);

	if(passAllBut(exclude_cuts)){
		if(pass.at(ChargeSum)){ //Opposite Sign WJets yield (bin #1 with value 1)
			if(value.at(MT_MuMET) > SB_lowerLimit){// && value.at(MT_MuMET) < SB_upperLimit){
				NSB.at(t).Fill(1., w);
			}
		}
		if(!pass.at(ChargeSum) && Charge != ChargeSumDummy){ //Same Sign WJets yield (bin #2 with value 2)
			if(value.at(MT_MuMET) > SB_lowerLimit){// && value.at(MT_MuMET) < SB_upperLimit){
				NSB.at(t).Fill(2., w);
			}
		}
	}

	// QCD ABCD BG Method
	/*******************
	 *        |        *
	 *    A   |    B   *  OS
	 *        |        *       S
	 * ----------------*------ i
	 *        |        *       g
	 *    C   |    D   *  SS   n
	 *        |        *
	 *******************
	 *  Iso   | AntiIso
	 *
	 *     relIso(mu)
	 */
	bool IsQCD_Event(false);
	if(!pass.at(NMuIso) && selMuon_AntiIso != selMuonDummy){
		if(!pass.at(ChargeSum) && value.at(ChargeSum) != ChargeSumDummy){
			if(id == DataMCType::Data){
				t = HConfig.GetType(DataMCType::QCD);
				IsQCD_Event = true;
			}
		}
	}
	if(verbose) std::cout << "QCD ABCD BG Method" << std::endl;
	exclude_cuts.push_back(NMuIso);
	if(passAllBut(exclude_cuts)){
		if(pass.at(NMuIso) && selMuon_Iso != selMuonDummy){
			if(pass.at(ChargeSum)){ //A --> Signal-Selection (pass all cuts except MT_MuMET)
				QCD_MT_MuMET_A.at(t).Fill(value.at(MT_MuMET), w);
				QCD_MET_A.at(t).Fill(Ntp->MET_CorrMVAMuTau_et(), w);
				if(pass.at(MT_MuMET)) NQCD.at(t).Fill(1., w);
			}
			if(!pass.at(ChargeSum) && value.at(ChargeSum) != ChargeSumDummy){ //C
				QCD_MT_MuMET_C.at(t).Fill(value.at(MT_MuMET), w);
				QCD_MET_C.at(t).Fill(Ntp->MET_CorrMVAMuTau_et(), w);
				if(pass.at(MT_MuMET)) NQCD.at(t).Fill(3., w);
			}
		}
		if(!pass.at(NMuIso) && selMuon_AntiIso != selMuonDummy){
			if(pass.at(ChargeSum)){ //B
				QCD_MT_MuMET_B.at(t).Fill(value.at(MT_MuMET), w);
				QCD_MET_B.at(t).Fill(Ntp->MET_CorrMVAMuTau_et(), w);
				if(pass.at(MT_MuMET)) NQCD.at(t).Fill(2., w);
			}
			if(!pass.at(ChargeSum) && value.at(ChargeSum) != ChargeSumDummy){ //D
				QCD_MT_MuMET_D.at(t).Fill(value.at(MT_MuMET), w);
				QCD_MET_D.at(t).Fill(Ntp->MET_CorrMVAMuTau_et(), w);
				if(pass.at(MT_MuMET)) NQCD.at(t).Fill(4., w);
				if(id == DataMCType::Data){
					if(pass.at(MT_MuMET)){
						NQCD.at(t).Fill(1., w);
						NQCD.at(HConfig.GetType(DataMCType::Data)).Fill(4., w);
					}
					QCD_MT_MuMET_A.at(t).Fill(value.at(MT_MuMET), w);
					QCD_MET_A.at(t).Fill(Ntp->MET_CorrMVAMuTau_et(), w);
				}
			}
		}
	}
	if(IsQCD_Event){
		pass.at(ChargeSum) = true;
		pass.at(NMuIso) = true;
	}

	bool status = AnalysisCuts(t,w,wobs);

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

	//single fit at default mass = 91.5
	std::vector<bool> A1Fit, EventFit; A1Fit.clear(); EventFit.clear();
	std::vector<double> Probs, Chi2s, Csums; Probs.clear(); Chi2s.clear(); Csums.clear();
	std::vector<LorentzVectorParticle> ZFits, daughter; ZFits.clear(); daughter.clear();
	std::vector<std::vector <LorentzVectorParticle> > RefitDaughters, InitDaughters; RefitDaughters.clear(), InitDaughters.clear();
	int IndexToReturn(-1), Niterat(-1);
	bool AmbiguityPoint(false);
	bool AmbiguitySolvable(false);
	std::vector<LorentzVectorParticle> TPTF_TausA1; TPTF_TausA1.clear();
	//std::vector<TVectorD> par_0, par; par_0.clear(), par.clear();

	if(status && value.at(TauFLSigma) != TauFLSigmaDummy){
		for(unsigned Ambiguity=0; Ambiguity<3; Ambiguity++){
			double LC_chi2(-1), phisign(0), csum(-1);
			LorentzVectorParticle Reco_TauA1, CorrectedReco_TauA1, Reco_Z;
			TVectorD tmp_par(3), tmp_par_0(3);
			std::vector <LorentzVectorParticle> tmp_Daughters, tmp_Daughters0; tmp_Daughters.clear(), tmp_Daughters0.clear();
			A1Fit.push_back(Ntp->ThreeProngTauFit(selTau, Ambiguity, Reco_TauA1, daughter, LC_chi2, phisign));
			//CorrectedReco_TauA1 = CorrectRecoTauMomentumBias(Reco_TauA1, Ntp->PFTau_p4(selTau), TPTF_TauA1_pBias);
			//Reco_TauA1 = CorrectedReco_TauA1;
			TPTF_TausA1.push_back(Reco_TauA1);
			if(A1Fit.at(Ambiguity)){
				Reco_A1Fit_Solution.at(t).Fill(Ambiguity, w);
				std::cout << "Chi2 for ambiguity " << Ambiguity << " : " << LC_chi2 << std::endl;
				bool EventFit_bool = Ntp->EventFit(selTau, selMuon, TPTF_TausA1.at(Ambiguity), Reco_Z, tmp_Daughters, tmp_Daughters0, LC_chi2, Niterat, csum, 91.5);//, tmp_par_0, tmp_par);
				//if(LC_chi2>= 0) EventFit.push_back(EventFit_bool);
				//else EventFit.push_back(false);
				EventFit.push_back(EventFit_bool);
				//std::cout << "Chi2 for ambiguity " << Ambiguity << " : " << LC_chi2 << std::endl;
				Reco_Chi2.at(t).Fill(LC_chi2, w);
				Reco_NIter.at(t).Fill(Niterat, w);
				Probs.push_back(TMath::Prob(LC_chi2, 1));
			}
			else{
				EventFit.push_back(false);
				Probs.push_back(0);
			}
			ZFits.push_back(Reco_Z);
			Chi2s.push_back(LC_chi2);
			Csums.push_back(csum);
			RefitDaughters.push_back(tmp_Daughters);
			InitDaughters.push_back(tmp_Daughters0);
			//par_0.push_back(tmp_par_0);
			//par.push_back(tmp_par);
			//std::cout << "Chi2 for ambiguity " << Ambiguity << " : " << Chi2s.at(Ambiguity) << std::endl;
			//std::cout << "TPTF status: " << A1Fit.at(Ambiguity) << " and DiTauFit status: " << EventFit.at(Ambiguity) << std::endl;
			//std::cout << "Chi2 probability for ambiguity " << Ambiguity << " : " << Probs.at(Ambiguity) << std::endl;
			//std::cout << "p(tau) for  " << Ambiguity << " : " << TPTF_TausA1.at(Ambiguity).LV().P() << std::endl;
		}
		for(unsigned int i=0;i<A1Fit.size();i++){
			if(A1Fit.at(i)){
				Reco_A1Fit_Solution.at(t).Fill(-1, w);
				break;
			}
			if(i==2){
				Reco_A1Fit_Solution.at(t).Fill(-2, w);
			}
		}
		std::cout << "Ambiguitysolver: " << std::endl;
		if(Ntp->AmbiguitySolverByChi2(A1Fit, EventFit, Chi2s, IndexToReturn, AmbiguityPoint)){
			int IndexToReturnTEST(-1); bool AmbiguityPointTEST(false);
			Ntp->AmbiguitySolver(A1Fit, EventFit, Probs, IndexToReturnTEST, AmbiguityPointTEST);
			AmbiguitySolvable = true;
			if(ZFits.at(IndexToReturn).Mass() >=0) Reco_ZMass.at(t).Fill(ZFits.at(IndexToReturn).LV().M(), w);
			Reco_EventFit_Solution.at(t).Fill(IndexToReturn, w);
			Reco_EventFit_Solution.at(t).Fill(-1, w); //all solutions
			Reco_ConstrainedDeltaSum.at(t).Fill(Csums.at(IndexToReturn), w);
			Reco_Chi2_FitSolutionOnly.at(t).Fill(Chi2s.at(IndexToReturn), w);
			std::cout << "Single Mass Fit; Fit Mass:  "<< ZFits.at(IndexToReturn).LV().M() << std::endl;
			std::cout << "Picked Solution with Chi2 for ambiguity " << IndexToReturn << " : " << Chi2s.at(IndexToReturn) << std::endl;
			std::cout << "Picked Solution with Ambiguity value (by Chi2s): " << IndexToReturn << std::endl;
			std::cout << "Picked Solution with Ambiguity value (by probs): " << IndexToReturnTEST << std::endl;
			Reco_TauMu_DeltaPX_FitImpact.at(t).Fill(RefitDaughters.at(IndexToReturn).at(1).LV().Px() - InitDaughters.at(IndexToReturn).at(1).LV().Px(), w);
			Reco_TauMu_DeltaPY_FitImpact.at(t).Fill(RefitDaughters.at(IndexToReturn).at(1).LV().Py() - InitDaughters.at(IndexToReturn).at(1).LV().Py(), w);
			Reco_TauMu_DeltaPZ_FitImpact.at(t).Fill(RefitDaughters.at(IndexToReturn).at(1).LV().Pz() - InitDaughters.at(IndexToReturn).at(1).LV().Pz(), w);
			Reco_TauA1_DeltaPX_FitImpact.at(t).Fill(RefitDaughters.at(IndexToReturn).at(0).LV().Px() - InitDaughters.at(IndexToReturn).at(0).LV().Px(), w);
			Reco_TauA1_DeltaPY_FitImpact.at(t).Fill(RefitDaughters.at(IndexToReturn).at(0).LV().Py() - InitDaughters.at(IndexToReturn).at(0).LV().Py(), w);
			Reco_TauA1_DeltaPZ_FitImpact.at(t).Fill(RefitDaughters.at(IndexToReturn).at(0).LV().Pz() - InitDaughters.at(IndexToReturn).at(0).LV().Pz(), w);
			//TVector3 par_Vec(par.at(IndexToReturn)(0),par.at(IndexToReturn)(1),par.at(IndexToReturn)(2));
			//TVector3 par0_Vec(par_0.at(IndexToReturn)(0),par_0.at(IndexToReturn)(1),par_0.at(IndexToReturn)(2));
			Reco_TauMu_ResPhi.at(t).Fill(RefitDaughters.at(IndexToReturn).at(1).LV().Phi() - InitDaughters.at(IndexToReturn).at(1).LV().Phi(), w);
			Reco_TauMu_ResCosTheta.at(t).Fill(RefitDaughters.at(IndexToReturn).at(1).LV().CosTheta() - InitDaughters.at(IndexToReturn).at(1).LV().CosTheta(), w);
			RecoZ_Pt.at(t).Fill(ZFits.at(IndexToReturn).LV().Pt(), w);
			Reco_dPhi_TauMuTauA1_BeforeFit.at(t).Fill(RefitDaughters.at(IndexToReturn).at(1).LV().Phi() - RefitDaughters.at(IndexToReturn).at(0).LV().Phi(), w);
			Reco_dPhi_TauMuTauA1_AfterFit.at(t).Fill(InitDaughters.at(IndexToReturn).at(1).LV().Phi() - InitDaughters.at(IndexToReturn).at(0).LV().Phi(), w);
		}
		else{
			std::cout << "Failed" << std::endl;
			Reco_EventFit_Solution.at(t).Fill(-2, w); //not able to solve ambiguity/no solution
		}
		std::cout << "-----------------End of Fit-----------------" << std::endl;
		if(id == DataMCType::DY_mutau_embedded){
			if(IndexToReturn == 0){
				TauFLSigma_vs_UnphysicalAll.at(t).Fill(value.at(TauFLSigma), 0);
				TauFL_vs_UnphysicalAll.at(t).Fill(Ntp->PFTau_FlightLength(selTau), 0);
			}
			TauFLSigma_vs_UnphysicalAll.at(t).Fill(value.at(TauFLSigma), 1);
			TauFL_vs_UnphysicalAll.at(t).Fill(Ntp->PFTau_FlightLength(selTau), 1);
		}
	}

	//multiple fits with different masses
	std::vector<bool> A1Fit_MassScan; A1Fit_MassScan.clear();
	std::vector<std::vector<bool> > EventFit_MassScan; EventFit_MassScan.clear();
	std::vector<std::vector<double> > Probs_MassScan; Probs_MassScan.clear();
	std::vector<std::vector<double> > Chi2s_MassScan; Chi2s_MassScan.clear();
	std::vector<std::vector<LorentzVectorParticle> > ZFits_MassScan, daughter_MassScan, Daughters_MassScan;
	ZFits_MassScan.clear(); daughter_MassScan.clear(); Daughters_MassScan.clear();
	std::vector<int> IndicesToReturn; IndicesToReturn.clear();
	std::vector<bool> AmbiguitySolver_MassScan; AmbiguitySolver_MassScan.clear();
	int FinalIndex(-1);

	unsigned N_mass = 50;
	std::vector<double> Masses;
	Masses.push_back(91.5);Masses.push_back(125.0);
	if(status && value.at(TauFLSigma) != TauFLSigmaDummy){
		//for(unsigned i_mass=0; i_mass<N_mass; i_mass++){
			//double MassConstraint = (double)i_mass*4.;
		for(unsigned i_mass=0; i_mass<Masses.size(); i_mass++){
			double MassConstraint = Masses.at(i_mass);
			std::vector<LorentzVectorParticle> ZFits; ZFits.clear();
			std::vector<double> Probs; Probs.clear();
			std::vector<double> Chi2s; Chi2s.clear();
			int IndexToReturn;
			LorentzVectorParticle Reco_Z;
			std::vector<bool> A1Fit, EventFit;
			for(unsigned Ambiguity=0; Ambiguity<3; Ambiguity++){
				LorentzVectorParticle Reco_TauA1; std::vector<LorentzVectorParticle> daughter;
				double LC_chi2(-1), phisign(0), csum(-1);
				A1Fit.push_back(Ntp->ThreeProngTauFit(selTau, Ambiguity, Reco_TauA1, daughter, LC_chi2, phisign));
				if(A1Fit.at(Ambiguity)){
					//Reco_A1Fit_Solution.at(t).Fill(Ambiguity, w);
					std::vector<LorentzVectorParticle> RefitDaughters, InitDaughters;
					int Niterat(-1);
					//TVectorD tmp1(3),tmp2(3);
					bool SingleEventFit = Ntp->EventFit(selTau, selMuon, Reco_TauA1, Reco_Z, RefitDaughters, InitDaughters, LC_chi2, Niterat, csum, MassConstraint);//), tmp2, tmp2);
					if(LC_chi2>= 0) EventFit.push_back(SingleEventFit);
					else EventFit.push_back(false);
					//Reco_Chi2.at(t).Fill(LC_chi2, w);
					//Reco_NIter.at(t).Fill(Niterat, w);
					Probs.push_back(TMath::Prob(LC_chi2, 1));
					//std::cout << "Chi2 probability for ambiguity " << Ambiguity << " : " << Probs.at(Ambiguity) << std::endl;
				}
				else{
					EventFit.push_back(false);
					Probs.push_back(0);
				}
				Chi2s.push_back(LC_chi2);
				ZFits.push_back(Reco_Z);
			}
			//std::cout << "Ambiguitysolver: " << std::endl;
			bool AmbiguityPoint;
			AmbiguitySolver_MassScan.push_back(Ntp->AmbiguitySolverByChi2(A1Fit, EventFit, Chi2s, IndexToReturn, AmbiguityPoint));
			if(AmbiguitySolver_MassScan.at(i_mass)){
				//if(ZFits.at(IndexToReturn).Mass() >=0) Reco_ZMass.at(t).Fill(ZFits.at(IndexToReturn).Mass(), w);
				//Reco_EventFit_Solution.at(t).Fill(IndexToReturn, w);
				//Reco_EventFit_Solution.at(t).Fill(-1, w); //all solutions
				//Reco_ConstrainedDeltaSum.at(t).Fill(csum, w);
				//std::cout << "Picked Solution with Ambiguity value: " << IndexToReturn << std::endl;
				Reco_ZMass_PDF.at(t).Fill(MassConstraint,w*Probs.at(IndexToReturn));
				IndicesToReturn.push_back(IndexToReturn);
			}
			else{
				//std::cout << "Failed" << std::endl;
				//Reco_EventFit_Solution.at(t).Fill(-2, w); //not able to solve ambiguity/no solution
				IndicesToReturn.push_back(-1);
			}
			ZFits_MassScan.push_back(ZFits);
			Probs_MassScan.push_back(Probs);
			Chi2s_MassScan.push_back(Chi2s);
			//std::cout << "Probs size" << Probs.size() << std::endl;
		}
		double maxProb(0);
		//std::cout << "Probs_MassScan size" << Probs_MassScan.size() << std::endl;
		//std::cout << "IndicesToReturn size" << IndicesToReturn.size() << std::endl;
		for(unsigned i=0; i<IndicesToReturn.size(); i++){
			if(AmbiguitySolver_MassScan.at(i)){
				Reco_ProbStack_MassScan.at(t).Fill(Probs_MassScan.at(i).at(IndicesToReturn.at(i)), w);
				if(Probs_MassScan.at(i).at(IndicesToReturn.at(i)) > maxProb){
					maxProb = Probs_MassScan.at(i).at(IndicesToReturn.at(i));
					FinalIndex = i;
				}
			}
		}
		//std::cout << "picked solution number: " << FinalIndex << std::endl;
		//std::cout << "highest prob: " << maxProb << std::endl;
		if(FinalIndex != -1){
			if(ZFits_MassScan.at(FinalIndex).at(IndicesToReturn.at(FinalIndex)).Mass() >=0){
				Reco_ZMass_MassScan.at(t).Fill(ZFits_MassScan.at(FinalIndex).at(IndicesToReturn.at(FinalIndex)).Mass(), w);
				Reco_ZMasswithProbWeight_MassScan.at(t).Fill(ZFits_MassScan.at(FinalIndex).at(IndicesToReturn.at(FinalIndex)).Mass(),w*Probs_MassScan.at(FinalIndex).at(IndicesToReturn.at(FinalIndex)));
				std::cout << "Multiple Masses Fit; Fit Mass:  " << ZFits_MassScan.at(FinalIndex).at(IndicesToReturn.at(FinalIndex)).Mass() << std::endl;
				std::cout << "Picked Solution with Chi2 for ambiguity " << IndicesToReturn.at(FinalIndex) << " : " << Chi2s_MassScan.at(FinalIndex).at(IndicesToReturn.at(FinalIndex)) << std::endl;
				std::cout << "Picked Solution with Ambiguity value: " << IndicesToReturn.at(FinalIndex) << std::endl;
			}
		}
	}

	///////////////////////////////////////////////////////////
	// Add plots
	if(status){
		NVtx.at(t).Fill(Ntp->NVtx(), w);
		unsigned int nGoodVtx=0;
		for(unsigned int i=0;i<Ntp->NVtx();i++){
			NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(), w);
			if(Ntp->isVtxGood(i))nGoodVtx++;
		}
		NGoodVtx.at(t).Fill(nGoodVtx, w);
		NMtTauMET.at(t).Fill(MT_TauMET, w);
		NMvis.at(t).Fill(Mvis, w);
		Mu_pt.at(t).Fill(Ntp->Muon_p4(selMuon).Pt(), w);
		Mu_phi.at(t).Fill(Ntp->Muon_p4(selMuon).Phi(), w);
		Mu_eta.at(t).Fill(Ntp->Muon_p4(selMuon).Eta(), w);
		Tau_pt.at(t).Fill(Ntp->PFTau_p4(selTau).Pt(), w);
		Tau_phi.at(t).Fill(Ntp->PFTau_p4(selTau).Phi(), w);
		Tau_eta.at(t).Fill(Ntp->PFTau_p4(selTau).Eta(), w);
		MET_phi.at(t).Fill(Ntp->MET_CorrMVAMuTau_phi(), w);
		A1mass.at(t).Fill(Ntp->PFTau_p4(selTau).M(), w);
		A1mass10GeV.at(t).Fill(Ntp->PFTau_p4(selTau).M(), w);
		TauFL_WithTauFLSigmaCut.at(t).Fill(Ntp->PFTau_FlightLength(selTau) , w);
		//std::cout << "PV Coord: " << std::endl;
		  //  Ntp->PFTau_TIP_primaryVertex_pos(selTau).Print();
		//std::cout << "SV Coord: " << std::endl;
		  //  Ntp->PFTau_TIP_secondaryVertex_pos(selTau).Print();

		TLorentzVector PFTau_UnFitTracks;
		(Ntp->PFTau_NdaughterTracks(selTau) == 3) ? TransTrk_Failure_withSelection.at(t).Fill(0.,w) : TransTrk_Failure_withSelection.at(t).Fill(1.,w);
		if(Ntp->PFTau_NdaughterTracks(selTau) == 3){
			for(unsigned i = 0; i<Ntp->PFTau_daughterTracks(selTau).size(); i++){
				//std::cout << "Ntp->PFTau_daughterTracks(selTau).size()" << Ntp->PFTau_daughterTracks(selTau).size() << std::endl;
				TrackParticle tmpTP = Ntp->PFTau_daughterTracks(selTau).at(i);
				TVector3 SV = Ntp->PFTau_TIP_secondaryVertex_pos(selTau);
				TLorentzVector tmpLV = (TrackTools::LorentzParticleAtPosition(tmpTP, SV)).LV();
				PFTau_UnFitTracks += tmpLV;
			}
		}
		Tau_Mass_Difference_PFTau_UnFitTracks_3PS.at(t).Fill(PFTau_UnFitTracks.M() - Ntp->PFTau_p4(selTau).M(), w);
	}

	//single fit at default mass = 91.5 with generator info
	std::vector<bool> GenA1Fit, GenEventFit; GenA1Fit.clear(); GenEventFit.clear();
	std::vector<double> GenProbs, GenChi2s, GenCsums; GenProbs.clear(); GenChi2s.clear(); GenCsums.clear();
	std::vector<LorentzVectorParticle> GenZFits, Gendaughter, GenDaughters; GenZFits.clear(); Gendaughter.clear(); GenDaughters.clear();
	int GenIndexToReturn(-1), GenNiterat(-1);
	bool GenAmbiguityPoint(false);
	bool GenAmbiguitySolvable(false);
	std::vector<LorentzVectorParticle> GenTPTF_TausA1; GenTPTF_TausA1.clear();

	//Gen Studies
	if(id != DataMCType::Data && (id == DataMCType::DY_tautau  || id== DataMCType::DY_mutau_embedded || (id >= DataMCType::H_tautau && id <= DataMCType::H_tautau_WHZHTTH))){
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
		int TauMu_charge(0);
		int Tauh_charge(0);
		if(Ntp->NMCTaus() == 2){
			for(int i=0; i<Ntp->NMCTaus(); i++){
				if(Ntp->MCTau_JAK(i) == 2){//Tau->Muon
					GenTaumu = Ntp->MCTau_p4(i);
					if(Ntp->MCTau_pdgid(i) == PDGInfo::tau_minus) TauMu_charge = -1;
					else if(Ntp->MCTau_pdgid(i) == PDGInfo::tau_plus) TauMu_charge = 1;
					for(int j=0; j<Ntp->NMCTauDecayProducts(i); j++){
						if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::mu_minus){
							GenMu = Ntp->MCTauandProd_p4(i,j);
							hasMu = true;
						}
					}
				}
				else if(Ntp->MCTau_JAK(i) != 2){
					GenTauh = Ntp->MCTau_p4(i);
					if(Ntp->MCTau_pdgid(i) == PDGInfo::tau_minus) Tauh_charge = -1;
					else if(Ntp->MCTau_pdgid(i) == PDGInfo::tau_plus) Tauh_charge = 1;
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
				Gen_TauMu_GJ.at(t).Fill(GenMu.Vect().Angle(GenTaumu.Vect()), w);
				Gen_TauA1_GJ.at(t).Fill(GenA1.Vect().Angle(GenTauh.Vect()), w);
				Phi_genTaumu.at(t).Fill(GenTaumu.Phi(), w);
				Theta_genTaumu.at(t).Fill(GenTaumu.Theta(), w);
				Phi_genTauh.at(t).Fill(GenTauh.Phi(), w);
				Theta_genTauh.at(t).Fill(GenTauh.Theta(), w);

				double dPhi_genTaus = GenTauh.Phi() - GenTaumu.Phi();
				//if(dPhi_genTaus > TMath::Pi()) dPhi_genTaus = 2*TMath::Pi() - dPhi_genTaus;
				TLorentzVector DiTau = GenTauh + GenTaumu;
				Gen_Z_Pt.at(t).Fill(GenZH.Pt(), w);
				Gen_DiTau_dPhi.at(t).Fill(dPhi_genTaus, w);
				Gen_DiTau_Pt.at(t).Fill(DiTau.Pt(), w);
				Gen_Z_M.at(t).Fill(GenZH.M(), w);
				Pt_vs_dPhi_DiTauGen.at(t).Fill(DiTau.Pt(), dPhi_genTaus);

				TLorentzVector GenTaumu_PtBalance = GenTauh;
				GenTaumu_PtBalance.SetX(GenTauh.X()*-1);
				GenTaumu_PtBalance.SetY(GenTauh.Y()*-1);

				double M_DiTau = (GenTaumu_PtBalance + GenTauh).M();
				double Delta_M_DiTau = M_DiTau - GenZH.M();

				Gen_DiTau_PtBalance_M.at(t).Fill(M_DiTau, w);
				Gen_DiTau_PtBalance_dM.at(t).Fill(Delta_M_DiTau, w);

				double dPt = GenTaumu_PtBalance.Pt() - GenTaumu.Pt();
				double dP_GenTaumu = GenTaumu_PtBalance.P() - GenTaumu.P();
				double GJ_Angle_Tauh = GenA1.Vect().Angle(GenTauh.Vect());
				double dP_Tauh = GenA1.E()*sqrt(pow(pow(GenTauh.M(),2) - pow(GenA1.M(),2),2) - pow(2*GenA1.P()*GenTauh.M()*sin(GJ_Angle_Tauh),2))/(pow(GenA1.M(),2) + pow(sin(GJ_Angle_Tauh)*GenA1.P(),2));

				Gen_TauMu_PtBalance_Pt.at(t).Fill(dPt, w);
				Gen_TauMu_PtBalance_dP.at(t).Fill(dP_GenTaumu, w);
				Gen_TauA1_dP.at(t).Fill(dP_Tauh, w);
				dP_GenTauMuPtBalance_vs_dPTauh.at(t).Fill(dP_GenTaumu, dP_Tauh);

				TLorentzVector GenA1_boosted(BoostToRestFrame(GenTauh, GenA1));
				if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) < 0){
					Gen_TPTF_TauA1_Solution_NoSelection.at(t).Fill(2., w);
				}
				else if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) > 0){
					Gen_TPTF_TauA1_Solution_NoSelection.at(t).Fill(1., w);
				}
				else if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) == 0){
					Gen_TPTF_TauA1_Solution_NoSelection.at(t).Fill(0., w);
				}

				if(status && value.at(TauFLSigma) != TauFLSigmaDummy){
					TVector3 POCAPV_dir = Ntp->Muon_Poca(selMuon_Iso) - Ntp->PFTau_TIP_primaryVertex_pos(selTau);
					POCAPV_Mag.at(t).Fill(POCAPV_dir.Mag());
					double Mvis_recoA1genMu = (Ntp->PFTau_p4(selTau) + GenMu).M();
					double Mvis_genA1recoMu = (GenA1 + Ntp->Muon_p4(selMuon_Iso)).M();
					double Mvis_recoA1genTaumu = (Ntp->PFTau_p4(selTau) + GenTaumu).M();
					double Mvis_genTauhrecoMu = (GenTauh + Ntp->Muon_p4(selMuon_Iso)).M();

					Mvis_SignalOnly.at(t).Fill(Mvis, w);
					Mvis_SignalOnly_genMu.at(t).Fill(Mvis_recoA1genMu, w);
					Mvis_SignalOnly_genA1.at(t).Fill(Mvis_genA1recoMu, w);
					Mvis_SignalOnly_genTaumu.at(t).Fill(Mvis_recoA1genTaumu, w);
					Mvis_SignalOnly_genTauh.at(t).Fill(Mvis_genTauhrecoMu, w);

					Phi_POCAPV.at(t).Fill(POCAPV_dir.Phi(), w);
					Theta_POCAPV.at(t).Fill(POCAPV_dir.Theta(), w);
					double dPhi = POCAPV_dir.Phi() - GenTaumu.Phi();
					dPhi_POCAPV_genTaumu.at(t).Fill(dPhi, w);
					double dTheta = POCAPV_dir.Theta() - GenTaumu.Theta();
					dTheta_POCAPV_genTaumu.at(t).Fill(dTheta, w);

					Phi_SVPV.at(t).Fill(Ntp->PFTau_FlightLength3d(selTau).Phi(), w);
					Theta_SVPV.at(t).Fill(Ntp->PFTau_FlightLength3d(selTau).Theta(), w);
					dPhi = Ntp->PFTau_FlightLength3d(selTau).Phi() - GenTauh.Phi();
					dPhi_SVPV_genTauh.at(t).Fill(dPhi, w);
					dPhi_SVPV_genTauh_vs_TauFL.at(t).Fill(dPhi, Ntp->PFTau_FlightLength(selTau));
					if(Tauh_charge == 1) dPhi_SVPV_genTauhPlus_vs_TauFL.at(t).Fill(dPhi, Ntp->PFTau_FlightLength(selTau));
					else if(Tauh_charge == -1) dPhi_SVPV_genTauhMinus_vs_TauFL.at(t).Fill(dPhi, Ntp->PFTau_FlightLength(selTau));
					dTheta = Ntp->PFTau_FlightLength3d(selTau).Theta() - GenTauh.Theta();
					dTheta_SVPV_genTauh.at(t).Fill(dTheta, w);
					Angle_SVPV_genTauh.at(t).Fill(Ntp->PFTau_FlightLength3d(selTau).Angle(GenTauh.Vect()));

					TVector3 MinusPFTau_FlightLength3d = - Ntp->PFTau_FlightLength3d(selTau);
					dPhi = MinusPFTau_FlightLength3d.Phi() - GenTaumu.Phi();
					dPhi_MinusSVPV_genTaumu.at(t).Fill(dPhi, w);
					dTheta = MinusPFTau_FlightLength3d.Theta() - GenTaumu.Theta();
					dTheta_MinusSVPV_genTaumu.at(t).Fill(dTheta, w);
					Angle_MinusSVPV_genTaumu.at(t).Fill(MinusPFTau_FlightLength3d.Angle(GenTaumu.Vect()));

					/*
					if(status && value.at(TauFLSigma) != TauFLSigmaDummy){
						for(unsigned Ambiguity=0; Ambiguity<3; Ambiguity++){
							double GenLC_chi2(-1), Genphisign(0), Gencsum(-1);
							LorentzVectorParticle GenReco_TauA1, GenReco_Z;
							GenA1Fit.push_back(Ntp->ThreeProngTauFit(selTau, Ambiguity, GenReco_TauA1, Gendaughter, GenLC_chi2, Genphisign));
							GenTPTF_TausA1.push_back(GenReco_TauA1);
							if(GenA1Fit.at(Ambiguity)){
								GenReco_A1Fit_Solution.at(t).Fill(Ambiguity, w);
								GenEventFit.push_back(Ntp->EventFit(selTau, selMuon, GenTPTF_TausA1.at(Ambiguity), GenReco_Z, GenDaughters, GenLC_chi2, GenNiterat, Gencsum, 91.5));
								GenReco_Chi2.at(t).Fill(GenLC_chi2, w);
								GenReco_NIter.at(t).Fill(GenNiterat, w);
								GenProbs.push_back(TMath::Prob(GenLC_chi2, 1));
								//std::cout << "Chi2 for ambiguity " << Ambiguity << " : " << Chi2s.at(Ambiguity) << std::endl;
								//std::cout << "TPTF status: " << A1Fit.at(Ambiguity) << " and DiTauFit status: " << EventFit.at(Ambiguity) << std::endl;
								//std::cout << "Chi2 probability for ambiguity " << Ambiguity << " : " << Probs.at(Ambiguity) << std::endl;
							}
							else{
								GenEventFit.push_back(false);
								GenProbs.push_back(0);
							}
							GenZFits.push_back(GenReco_Z);
							GenChi2s.push_back(GenLC_chi2);
							GenCsums.push_back(Gencsum);
						}
						for(unsigned int i=0;i<GenA1Fit.size();i++){
							if(GenA1Fit.at(i)){
							  GenReco_A1Fit_Solution.at(t).Fill(-1, w);
								break;
							}
							if(i==2){
								Reco_A1Fit_Solution.at(t).Fill(-2, w);
							}
						}
						//std::cout << "Ambiguitysolver: " << std::endl;
						if(Ntp->AmbiguitySolver(GenA1Fit, GenEventFit, GenProbs, GenIndexToReturn, GenAmbiguityPoint)){
							GenAmbiguitySolvable = true;
							if(GenZFits.at(GenIndexToReturn).Mass() >=0) Reco_ZMass.at(t).Fill(GenZFits.at(GenIndexToReturn).Mass(), w);
							GenReco_EventFit_Solution.at(t).Fill(GenIndexToReturn, w);
							GenReco_EventFit_Solution.at(t).Fill(-1, w); //all solutions
							GenReco_ConstrainedDeltaSum.at(t).Fill(GenCsums.at(GenIndexToReturn), w);
							GenReco_Chi2_FitSolutionOnly.at(t).Fill(GenChi2s.at(GenIndexToReturn), w);
							//std::cout << "Picked Solution with Chi2 for ambiguity " << IndexToReturn << " : " << Chi2s.at(IndexToReturn) << std::endl;
							//std::cout << "Picked Solution with Ambiguity value: " << IndexToReturn << std::endl;
						}
						else{
							//std::cout << "Failed" << std::endl;
							GenReco_EventFit_Solution.at(t).Fill(-2, w); //not able to solve ambiguity/no solution
						}
						//std::cout << "-----------------End of Fit-----------------" << std::endl;
					}
					*/

					TLorentzVector Neutrino_UnFitTracks(GenTauh);
					(Ntp->PFTau_NdaughterTracks(selTau) == 3) ? TransTrk_Failure_withSelection.at(t).Fill(0.,w) : TransTrk_Failure_withSelection.at(t).Fill(1.,w);
					if(Ntp->PFTau_NdaughterTracks(selTau) == 3){
						for(unsigned i = 0; i<Ntp->PFTau_daughterTracks(selTau).size(); i++){
							//std::cout << "Ntp->PFTau_daughterTracks(selTau).size()" << Ntp->PFTau_daughterTracks(selTau).size() << std::endl;
							TrackParticle tmpTP = Ntp->PFTau_daughterTracks(selTau).at(i);
							TVector3 SV = Ntp->PFTau_TIP_secondaryVertex_pos(selTau);
							TLorentzVector tmpLV = (TrackTools::LorentzParticleAtPosition(tmpTP, SV)).LV();
							Neutrino_UnFitTracks-=tmpLV;
						}
					}
					TPTF_Neutrino_UnFitTracks_Mass.at(t).Fill(Neutrino_UnFitTracks.M(), w);
					TPTF_Neutrino_UnFitTracks_Mass_vs_TauFL.at(t).Fill(Neutrino_UnFitTracks.M(), Ntp->PFTau_FlightLength(selTau));

					TLorentzVector Neutrino_ReFitTracks(GenTauh);
					if(Ntp->PFTau_NdaughtersReFitTracks_p4(selTau) == 3){
						for(unsigned i = 0; i<Ntp->PFTau_daughterReFitTracks_p4(selTau).size(); i++){
							//std::cout << "Ntp->PFTau_daughterTracks(selTau).size()" << Ntp->PFTau_daughterTracks(selTau).size() << std::endl;
							TLorentzVector tmpTLV = Ntp->PFTau_daughterReFitTracks_p4(selTau).at(i);
							Neutrino_ReFitTracks-=tmpTLV;
						}
					}
					TPTF_Neutrino_ReFitTracks_Mass.at(t).Fill(Neutrino_ReFitTracks.M(), w);
					TPTF_Neutrino_ReFitTracks_Mass_vs_TauFL.at(t).Fill(Neutrino_ReFitTracks.M(), Ntp->PFTau_FlightLength(selTau));

					TLorentzVector Neutrino_PFTau(GenTauh - Ntp->PFTau_p4(selTau));
					TPTF_Neutrino_PFTau_Mass.at(t).Fill(Neutrino_PFTau.M(), w);
					TPTF_Neutrino_PFTau_Mass_vs_TauFL.at(t).Fill(Neutrino_PFTau.M(), Ntp->PFTau_FlightLength(selTau));

					if(AmbiguitySolvable){
						if(GenZH.Pt() < 5.){
							if(ZFits.at(IndexToReturn).Mass() >=0){
								Reco_ZMass_UnboostedGenZ.at(t).Fill(ZFits.at(IndexToReturn).Mass(), w);
								RecoZ_Pt_Unboosted.at(t).Fill(ZFits.at(IndexToReturn).Parameter(LorentzVectorParticle::pt), w);
								GenZ_Pt_Unboosted.at(t).Fill(GenZH.Pt(), w);
							}
							//if(FinalIndex != -1) Reco_ZMass_MassScanUnboosted.at(t).Fill(ZFits_MassScan.at(FinalIndex).at(IndicesToReturn.at(FinalIndex)).Mass(),w);
						}
						int IndexNoFit(IndexToReturn); //always solution 1
						//if(IndexToReturn == 2) IndexNoFit = 1;
						//else IndexNoFit = IndexToReturn;
						TLorentzVector TauA1_TLV_FitSolution = RefitDaughters.at(IndexToReturn).at(0).LV();
						TLorentzVector Z_TLV = ZFits.at(IndexToReturn).LV();
						TLorentzVector TauMu_TLV = RefitDaughters.at(IndexToReturn).at(1).LV();
						//TVector3 par_Vec(par.at(IndexToReturn)(0),par.at(IndexToReturn)(1),par.at(IndexToReturn)(2));
						//TVector3 par0_Vec(par_0.at(IndexToReturn)(0),par_0.at(IndexToReturn)(1),par_0.at(IndexToReturn)(2));
						double dPt_TauA1 = TauA1_TLV_FitSolution.Pt() - GenTauh.Pt();
						double dPt_TauMu = TauMu_TLV.Pt() - GenTaumu.Pt();
						double dPt_TauA1_NoFit = RefitDaughters.at(IndexNoFit).at(0).LV().Pt() - GenTauh.Pt();
						double dPt_TauMu_NoFit = InitDaughters.at(IndexNoFit).at(1).LV().Pt() - GenTaumu.Pt();

						double dP_TauA1_FitSolution = TauA1_TLV_FitSolution.P() - GenTauh.P();
						double dE_Z = Z_TLV.E() - GenZH.E();
						double dP_TauA1_BestSolution(999); int i_BestSolution(-1);
						for(unsigned i=0; i<A1Fit.size(); i++){
							TLorentzVector TauA1_TLV_tmp = TPTF_TausA1.at(i).LV();
							if(fabs(dP_TauA1_BestSolution) > fabs(TauA1_TLV_tmp.P() - GenTauh.P())){
								dP_TauA1_BestSolution = TauA1_TLV_tmp.P() - GenTauh.P();
								i_BestSolution = i;
							}
						}

						if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) < 0){
							TPTF_TauA1_RightSolution_vs_FitSolution.at(t).Fill(2., IndexToReturn);
							Gen_TPTF_TauA1_Solution_WithSelection.at(t).Fill(2., w);
						}
						else if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) > 0){
							TPTF_TauA1_RightSolution_vs_FitSolution.at(t).Fill(1., IndexToReturn);
							Gen_TPTF_TauA1_Solution_WithSelection.at(t).Fill(1., w);
						}
						else if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) == 0){
							TPTF_TauA1_RightSolution_vs_FitSolution.at(t).Fill(0., IndexToReturn);
							Gen_TPTF_TauA1_Solution_WithSelection.at(t).Fill(0., w);
						}
						Reco_Z_Energy_Res.at(t).Fill(dE_Z, w);
						Reco_PtRes_TauMu.at(t).Fill(dPt_TauMu, w);
						Reco_PtRes_TauA1.at(t).Fill(dPt_TauA1, w);
						Reco_PtRes_TauMu_NoFit.at(t).Fill(dPt_TauMu_NoFit, w);
						Reco_PtRes_TauA1_NoFit.at(t).Fill(dPt_TauA1_NoFit, w);

						if(IndexToReturn == 0){
							Reco_PtRes_TauMu_AmbPoint0.at(t).Fill(dPt_TauMu, w);
							Reco_PtRes_TauA1_AmbPoint0.at(t).Fill(dPt_TauA1, w);
							Reco_PtRes_TauMu_AmbPoint0_NoFit.at(t).Fill(dPt_TauMu_NoFit, w);
							Reco_PtRes_TauA1_AmbPoint0_NoFit.at(t).Fill(dPt_TauA1_NoFit, w);
						}

						if(IndexToReturn == 1 || IndexToReturn == 2){
							Reco_PtRes_TauMu_AmbPoint12.at(t).Fill(dPt_TauMu, w);
							Reco_PtRes_TauA1_AmbPoint12.at(t).Fill(dPt_TauA1, w);
							Reco_PtRes_TauMu_AmbPoint12_NoFit.at(t).Fill(dPt_TauMu_NoFit, w);
							Reco_PtRes_TauA1_AmbPoint12_NoFit.at(t).Fill(dPt_TauA1_NoFit, w);
						}

						if(IndexToReturn == 1){
							Reco_PtRes_TauMu_AmbPoint1.at(t).Fill(dPt_TauMu, w);
							Reco_PtRes_TauA1_AmbPoint1.at(t).Fill(dPt_TauA1, w);
						}

						double dP_A1 = Ntp->PFTau_p4(selTau).P() - GenA1.P();
						double dPx_TauA1 = TauA1_TLV_FitSolution.Px() - GenTauh.Px();
						double dPy_TauA1 = TauA1_TLV_FitSolution.Py() - GenTauh.Py();
						double dPz_TauA1 = TauA1_TLV_FitSolution.Pz() - GenTauh.Pz();
						double dPxsq_TauA1 = pow(TauA1_TLV_FitSolution.Px(),2.) - pow(GenTauh.Px(),2.);
						double dPysq_TauA1 = pow(TauA1_TLV_FitSolution.Py(),2.) - pow(GenTauh.Py(),2.);
						double dPzsq_TauA1 = pow(TauA1_TLV_FitSolution.Pz(),2.) - pow(GenTauh.Pz(),2.);
						double dP_isq_TauA1 = pow(dPx_TauA1,2) + pow(dPy_TauA1,2) + pow(dPz_TauA1,2);
						//double dPt_TauA1 = TauA1_TLV_FitSolution.Pt() - GenTauh.Pt();
						double GenGJAngle = GenA1.Vect().Angle(GenTauh.Vect());
						double RecoGJAngle = Ntp->PFTau_p4(selTau).Vect().Angle(TauA1_TLV_FitSolution.Vect());
						double GJ_GenTauh_RecoA1 = Ntp->PFTau_p4(selTau).Vect().Angle(GenTauh.Vect());
						double dP_ortho_TauA1_RecoA1 = TauA1_TLV_FitSolution.P()*cos(RecoGJAngle) - GenTauh.P()*cos(GJ_GenTauh_RecoA1);
						double dP_paral_TauA1_RecoA1 = TauA1_TLV_FitSolution.P()*sin(RecoGJAngle) - GenTauh.P()*sin(GJ_GenTauh_RecoA1);

						TPTF_TauA1_pRes_FitSolution.at(t).Fill(dP_TauA1_FitSolution, w);
						TPTF_TauA1_pRes_BestSolution.at(t).Fill(dP_TauA1_BestSolution, w);
						TPTF_TauA1_BestSolution_vs_FitSolution.at(t).Fill(i_BestSolution, IndexToReturn);
						TPTF_TauA1_pRes_vs_GenA1Mass_FitSolution.at(t).Fill(dP_TauA1_FitSolution, GenA1.M());
						TPTF_TauA1_pRes_vs_GenA1Mass_BestSolution.at(t).Fill(dP_TauA1_BestSolution, GenA1.M());
						TPTF_TauA1_pRes_vs_GenGJAngle_FitSolution.at(t).Fill(dP_TauA1_FitSolution, GenGJAngle);
						TPTF_TauA1_pRes_vs_GenGJAngle_BestSolution.at(t).Fill(dP_TauA1_BestSolution, GenGJAngle);
						TPTF_A1_pRes_vs_GenGJAngle.at(t).Fill(dP_A1, GenGJAngle);

						TPTF_TauA1_pRes_vs_RecoChi2_FitSolution.at(t).Fill(dP_TauA1_FitSolution, Chi2s.at(IndexToReturn));
						TPTF_TauA1_pRes_vs_RecoA1Mass_FitSolution.at(t).Fill(dP_TauA1_FitSolution, Ntp->PFTau_p4(selTau).M());
						TPTF_TauA1_pRes_vs_RecoGJAngle_FitSolution.at(t).Fill(dP_TauA1_FitSolution, RecoGJAngle);
						TPTF_TauA1_p_iRes_sq_vs_RecoGJAngle_FitSolution.at(t).Fill(dP_isq_TauA1, RecoGJAngle);
						TPTF_TauA1_pxRes_vs_RecoGJAngle_FitSolution.at(t).Fill(dPx_TauA1, RecoGJAngle);
						TPTF_TauA1_pyRes_vs_RecoGJAngle_FitSolution.at(t).Fill(dPy_TauA1, RecoGJAngle);
						TPTF_TauA1_pzRes_vs_RecoGJAngle_FitSolution.at(t).Fill(dPz_TauA1, RecoGJAngle);
						TPTF_TauA1_ptRes_vs_RecoGJAngle_FitSolution.at(t).Fill(dPt_TauA1, RecoGJAngle);
						TPTF_TauA1_pxsqRes_vs_RecoGJAngle_FitSolution.at(t).Fill(dPxsq_TauA1, RecoGJAngle);
						TPTF_TauA1_pysqRes_vs_RecoGJAngle_FitSolution.at(t).Fill(dPysq_TauA1, RecoGJAngle);
						TPTF_TauA1_pzsqRes_vs_RecoGJAngle_FitSolution.at(t).Fill(dPzsq_TauA1, RecoGJAngle);
						TPTF_TauA1_p_orthoRes_vs_RecoGJAngle_FitSolution.at(t).Fill(dP_ortho_TauA1_RecoA1, RecoGJAngle);
						TPTF_TauA1_p_paralRes_vs_RecoGJAngle_FitSolution.at(t).Fill(dP_paral_TauA1_RecoA1, RecoGJAngle);
						TPTF_A1_pRes_vs_RecoGJAngle.at(t).Fill(dP_A1, RecoGJAngle);

						TPTF_TauA1_p_Reco.at(t).Fill(TauA1_TLV_FitSolution.P(), w);
						TPTF_TauA1_pt_Reco.at(t).Fill(TauA1_TLV_FitSolution.Pt(), w);
						TPTF_TauA1_px_Reco.at(t).Fill(TauA1_TLV_FitSolution.Px(), w);
						TPTF_TauA1_py_Reco.at(t).Fill(TauA1_TLV_FitSolution.Py(), w);
						TPTF_TauA1_pz_Reco.at(t).Fill(TauA1_TLV_FitSolution.Pz(), w);
						TPTF_TauA1_p_Gen.at(t).Fill(GenTauh.P(), w);
						TPTF_TauA1_pt_Gen.at(t).Fill(GenTauh.Pt(), w);
						TPTF_TauA1_px_Gen.at(t).Fill(GenTauh.Px(), w);
						TPTF_TauA1_py_Gen.at(t).Fill(GenTauh.Py(), w);
						TPTF_TauA1_pz_Gen.at(t).Fill(GenTauh.Pz(), w);

						TPTF_TauA1_pxsq_Reco.at(t).Fill(pow(TauA1_TLV_FitSolution.Px(),2.), w);
						TPTF_TauA1_pysq_Reco.at(t).Fill(pow(TauA1_TLV_FitSolution.Py(),2.), w);
						TPTF_TauA1_pzsq_Reco.at(t).Fill(pow(TauA1_TLV_FitSolution.Pz(),2.), w);
						TPTF_TauA1_pxsq_Gen.at(t).Fill(pow(GenTauh.Px(),2.), w);
						TPTF_TauA1_pysq_Gen.at(t).Fill(pow(GenTauh.Py(),2.), w);
						TPTF_TauA1_pzsq_Gen.at(t).Fill(pow(GenTauh.Pz(),2.), w);

						TPTF_TauA1_ptRes_vs_ptGen.at(t).Fill(dPt_TauA1, GenTauh.Pt());
						TPTF_TauA1_ptRes_vs_ptReco.at(t).Fill(dPt_TauA1, TauA1_TLV_FitSolution.Pt());
					}
				}
			}
		}
	}

	if(passAllBut(TauFLSigma)){
		TauFLSigmaSigned.at(t).Fill(value.at(TauFLSigma), w);
		if(value.at(TauFLSigma)<0 && value.at(TauFLSigma) != TauFLSigmaDummy){
			TauFLSigmaUnsigned.at(t).Fill(-value.at(TauFLSigma), w);
			TauFL_NoTauFLSigmaCut.at(t).Fill(Ntp->PFTau_FlightLength(selTau), w);
			TauFLSigned_NoTauFLSigmaCut.at(t).Fill(-Ntp->PFTau_FlightLength(selTau), w);
		}
		else if(value.at(TauFLSigma)>0){
			TauFLSigmaUnsigned.at(t).Fill(value.at(TauFLSigma), w);
			TauFL_NoTauFLSigmaCut.at(t).Fill(Ntp->PFTau_FlightLength(selTau), w);
			TauFLSigned_NoTauFLSigmaCut.at(t).Fill(Ntp->PFTau_FlightLength(selTau), w);
		}
	}

	std::vector<unsigned int> exclude_MTDecay;
	exclude_MTDecay.clear();
	exclude_MTDecay.push_back(MT_MuMET);
	exclude_MTDecay.push_back(TauDecayMode);
	exclude_MTDecay.push_back(TauFLSigma);
	if(passAllBut(exclude_MTDecay)){
		if(pass.at(MT_MuMET)){
			Tau_Mass_Inclusive.at(t).Fill(Ntp->PFTau_p4(selTau).M(), w);
			Tau_Mass_sq_Inclusive.at(t).Fill(Ntp->PFTau_p4(selTau).M2(), w);
			Tau_Mass_Inclusive_NoTLV.at(t).Fill(Ntp->PFTau_Mass(selTau), w);
			std::cout << "Ntp->PFTau_NdaughtersReFitTracks_p4(selTau) " << Ntp->PFTau_NdaughtersReFitTracks_p4(selTau) << std::endl;
			for(unsigned int i=0; i<Ntp->PFTau_NdaughtersReFitTracks_p4(selTau);i++){
				std::cout << "Ntp->PFTau_daughterReFitTracks_p4(selTau).at(i).M() " << Ntp->PFTau_daughterReFitTracks_p4(selTau).at(i).M() << std::endl;
				Tau_Mass_Inclusive_ReFitTracks.at(t).Fill(Ntp->PFTau_daughterReFitTracks_p4(selTau).at(i).M(), w);

			TLorentzVector PFTau_UnFitTracks;
			for(unsigned int i=0; i<Ntp->PFTau_NdaughterTracks(selTau);i++){
				//std::cout << "Ntp->PFTau_daughterTracks(selTau).size()" << Ntp->PFTau_daughterTracks(selTau).size() << std::endl;
				TrackParticle tmpTP = Ntp->PFTau_daughterTracks(selTau).at(i);
				TVector3 SV = Ntp->PFTau_TIP_secondaryVertex_pos(selTau);
				TLorentzVector tmpLV = (TrackTools::LorentzParticleAtPosition(tmpTP, SV)).LV();
				Tau_Mass_Inclusive_UnFitTracks.at(t).Fill(tmpLV.M(), w);
			}
			}
			MvisIncl.at(t).Fill(Mvis, w);
		}
		MTMuMETIncl.at(t).Fill(value.at(MT_MuMET), w);
		if(pass.at(TauDecayMode)){
			if(pass.at(MT_MuMET))  Mvis3Prong.at(t).Fill(Mvis, w);
			MTMuMET3Prong.at(t).Fill(value.at(MT_MuMET), w);
		}
		if(!pass.at(TauDecayMode)){
			if(pass.at(MT_MuMET)) Mvis1Prong.at(t).Fill(Mvis, w);
			MTMuMET1Prong.at(t).Fill(value.at(MT_MuMET), w);
		}
	}

	for(unsigned i=0; i<Ntp->NPFTaus(); i++){
		if(Ntp->PFTau_hpsDecayMode(i) == 10){
			if(Ntp->PFTau_NdaughterTracks(i) == 3){
				TransTrk_Failure_noSelection.at(t).Fill(0., w);
			}
			else TransTrk_Failure_noSelection.at(t).Fill(1., w);
		}
	}

}//final bracket of DoEvent


void  ZToTaumuTauh::Finish(){
	std::cout << "Starting Finish" << std::endl;
	if(mode==RECONSTRUCT){
		std::cout << "Enter mode==RECONSTRUCT" << std::endl;
		SkimConfig SC;
		SC.ApplySkimEfficiency(types,Npassed,Npassed_noweight);


		// Scale DY embedded sample to the yield of DY MC and scale DY MC to 0 afterwards
		double NDY_tautau(Npassed.at(HConfig.GetType(DataMCType::DY_tautau)).GetBinContent(MT_MuMET+1)*CrossSectionandAcceptance.at(HConfig.GetType(DataMCType::DY_tautau))*Lumi/Npassed.at(HConfig.GetType(DataMCType::DY_tautau)).GetBinContent(0));
		double NDY_tautau_Emb(Npassed.at(HConfig.GetType(DataMCType::DY_mutau_embedded)).GetBinContent(MT_MuMET+1));

		std::cout << "NDY_tautau: " << NDY_tautau << std::endl;
		std::cout << "NDY_tautau_Emb: " << NDY_tautau_Emb << std::endl;
		std::cout << "NDY_tautau/NDY_tautau_Emb: " << NDY_tautau/NDY_tautau_Emb << std::endl;

		int Exclude_DY_ID;
		if(Use_Embedded){
			if(NDY_tautau>0 && NDY_tautau_Emb>0){
				ScaleAllHistOfType(HConfig.GetType(DataMCType::DY_mutau_embedded),NDY_tautau/NDY_tautau_Emb);
				//ScaleAllHistOfType(HConfig.GetType(DataMCType::DY_tautau),0);
				suppressDrawingHistOfType(HConfig.GetType(DataMCType::DY_tautau));
				Exclude_DY_ID = DataMCType::DY_tautau;
			}
		}
		else{
			//ScaleAllHistOfType(HConfig.GetType(DataMCType::DY_mutau_embedded),0);
			suppressDrawingHistOfType(HConfig.GetType(DataMCType::DY_mutau_embedded));
			Exclude_DY_ID = DataMCType::DY_mutau_embedded;
		}

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
				QCD_Integral_B_Data_minus_MC += QCD_Integral_B.at(i);
				QCD_Integral_C_Data_minus_MC += QCD_Integral_C.at(i);
				QCD_Integral_D_Data_minus_MC += QCD_Integral_D.at(i);
		  }
		  else if(HConfig.GetID(i) != DataMCType::QCD && HConfig.GetID(i) != Exclude_DY_ID){
				QCD_Integral_B_Data_minus_MC -= QCD_Integral_B.at(i);
				QCD_Integral_C_Data_minus_MC -= QCD_Integral_C.at(i);
				QCD_Integral_D_Data_minus_MC -= QCD_Integral_D.at(i);
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
		for(unsigned i=1; i<41; i++){
			double ratioTFLSigma = TauFLSigma_vs_UnphysicalAll.at(HConfig.GetType(DataMCType::DY_mutau_embedded)).GetBinContent(i, 1)/TauFLSigma_vs_UnphysicalAll.at(HConfig.GetType(DataMCType::DY_mutau_embedded)).GetBinContent(i, 2);
			if(ratioTFLSigma >= 0){
				TauFLSigma_vs_UnphysicalProb.at(HConfig.GetType(DataMCType::DY_mutau_embedded)).SetBinContent(i,ratioTFLSigma);
			}
			else TauFLSigma_vs_UnphysicalProb.at(HConfig.GetType(DataMCType::DY_mutau_embedded)).SetBinContent(i,0.);

			double ratioTFL = TauFL_vs_UnphysicalAll.at(HConfig.GetType(DataMCType::DY_mutau_embedded)).GetBinContent(i, 1)/TauFL_vs_UnphysicalAll.at(HConfig.GetType(DataMCType::DY_mutau_embedded)).GetBinContent(i, 2);
			if(ratioTFL >= 0) TauFL_vs_UnphysicalProb.at(HConfig.GetType(DataMCType::DY_mutau_embedded)).SetBinContent(i,ratioTFL);
			else TauFL_vs_UnphysicalProb.at(HConfig.GetType(DataMCType::DY_mutau_embedded)).SetBinContent(i,0.);
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
	if (	Ntp->PFTau_isHPSByDecayModeFinding(i) &&
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
		if(Ntp->PFTau_p4(i).DeltaR(Ntp->Muon_p4(*it_mu)) < cMuTau_dR){
			return false;
		}
	}
	// trigger matching
	if(Ntp->GetMCID() != DataMCType::DY_mutau_embedded){
		if(Ntp->matchTrigger(Ntp->PFTau_p4(i),cTriggerNames,"tau") > cTau_dRHltMatch){
		  return false;
		}
	}
	if(selectPFTau_Id(i)){
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
	if(		Ntp->PFTau_p4(i).Pt() >= cTau_pt &&
			fabs(Ntp->PFTau_p4(i).Eta()) <= cTau_eta
			){
		return true;
	}
	return false;
}
double ZToTaumuTauh::Reconstruct_hadronicTauEnergy(unsigned i){
	double TauEnergy,TauMomentumPlus, TauMomentumMinus;
	TLorentzVector a1 = Ntp->PFTau_p4(i);
	double GJ_angle = a1.Angle(Ntp->PFTau_FlightLength3d(i));
	double val1 = (pow(PDG_Var::Tau_mass(),2.) + pow(a1.M(),2.))*a1.P()*cos(GJ_angle);
	double val2 = a1.Energy()*sqrt(pow(pow(PDG_Var::Tau_mass(),2.) - pow(a1.M(),2.),2.) - 4.*pow(a1.P()*PDG_Var::Tau_mass()*sin(GJ_angle),2.));
	TauMomentumPlus = (val1 + val2)/2./(pow(a1.M(),2) + pow(a1.P()*sin(GJ_angle),2.));
	TauMomentumMinus = (val1 - val2)/2./(pow(a1.M(),2) + pow(a1.P()*sin(GJ_angle),2.));
	return TauEnergy;
}
LorentzVectorParticle ZToTaumuTauh::CorrectRecoTauMomentumBias(LorentzVectorParticle RecoTau, TLorentzVector RecoA1, std::vector<double> BiasInGJAngleBins){
	double Angle = RecoA1.Vect().Angle(RecoTau.LV().Vect());
	double BinWidth = 0.035/BiasInGJAngleBins.size();
	unsigned int i=0;
	while((i*BinWidth < Angle) && (i < BiasInGJAngleBins.size())) i++;
	double alpha = (RecoTau.Parameter(LorentzVectorParticle::p) + BiasInGJAngleBins.at(i))/RecoTau.Parameter(LorentzVectorParticle::p);
	TMatrixT<double> correctedPar = RecoTau.getParMatrix();
	TMatrixTSym<double> correctedCov = RecoTau.getCovMatrix();
	correctedPar(LorentzVectorParticle::px,0) = RecoTau.Parameter(LorentzVectorParticle::px)*alpha;
	correctedPar(LorentzVectorParticle::py,0) = RecoTau.Parameter(LorentzVectorParticle::py)*alpha;
	correctedPar(LorentzVectorParticle::pz,0) = RecoTau.Parameter(LorentzVectorParticle::pz)*alpha;
	LorentzVectorParticle CorrectedRecoTau(correctedPar, correctedCov, RecoTau.PDGID(), RecoTau.Charge(), RecoTau.BField());
	return CorrectedRecoTau;
}

//gets two TLorentzVectors (1. defines RF/Boost, 2. is boosted), creates a copy of the second and boosts it
TLorentzVector ZToTaumuTauh::BoostToRestFrame(TLorentzVector TLV1, TLorentzVector TLV2){
	TVector3 boostvector = TLV1.BoostVector();
	TLorentzVector boosted_TLV2(TLV2);
	boosted_TLV2.Boost(-boostvector);
	return boosted_TLV2;
}
TLorentzVector ZToTaumuTauh::TauHelixP4AtSV(unsigned int selTau, TLorentzVector Tau, TVector3 PV, TVector3 SV){
	TLorentzVector Helix_TLV(Tau);
	double pTau = Tau.P();
	TVector3 PVSV = TVector3(SV - PV);
	TVector3 p0 = TVector3(PVSV); double p0Mag = p0.Mag(); double p0Scale = pTau/p0Mag; p0 = p0*p0Scale;
	double Tau_SinLambda = p0.Pt()/p0.Mag();
	double Tau_Lambda = asin(Tau_SinLambda);
	double alpha = -Ntp->PFTau_Charge(selTau)*3*10**(-3)*3.8;
	double Tau_R = pTau/alpha;
	double kappa = 1/Tau_R/2.;
	double Tau_Phi0 = atan(p0.Y()/p0.X()) - (TMath::PiOver2() - acos(PVSV.Mag()*kappa));
	double Tau_Phi02 = atan(p0.Y()/p0.X()) + (TMath::PiOver2() - acos(PVSV.Mag()*kappa));
	TLorentzVector Helix_TLV = TLorentzVector(p0.Pt()*cos(Tau_Phi02), p0.Pt()*sin(Tau_Phi02), p0.Z(), Tau.E());
	return Helix_TLV;
}
