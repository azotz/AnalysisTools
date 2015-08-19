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
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "Objects.h"

ZToTaumuTauh::ZToTaumuTauh(TString Name_, TString id_):
	Selection(Name_,id_),
	cMu_dxy(0.045),// 100 is dummy value
	cMu_dz(0.2),
	cMu_relIso(0.1),
	cMu_pt(20.),
	cMu_eta(2.1),
	cMu_dRHltMatch(0.5),
	cDiMuVeto_pt(15.),
	cDiMuVeto_eta(2.4),
	cDiMuVeto_dz(0.2),
	cDiMuVeto_dBetaCombIso(0.3),
	cDiMuVeto_dR(0.15),
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
	selection_verbose = false;

	Use_Embedded = true;

	Logger::Instance()->SetLevel(Logger::Info);
}

ZToTaumuTauh::~ZToTaumuTauh(){
	for(unsigned int j=0; j<Npassed.size(); j++){
		Logger(Logger::Info) << "Selection Summary before: "
		<< Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		<< Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts+1) << std::endl;
	}
	Logger(Logger::Info) << "" << std::endl;
	Logger::Instance()->SetLevel(Logger::Debug);
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
		if(i==DiMuonVeto)		cut.at(DiMuonVeto)=0.15;
		if(i==NTauId)			cut.at(NTauId)=1;
		if(i==NTauKin)			cut.at(NTauKin)=1;
		if(i==NTauIso)			cut.at(NTauIso)=1;
		if(i==ChargeSum)		cut.at(ChargeSum)=0;
		if(i==TauDecayMode)		cut.at(TauDecayMode)=10;//10
		if(i==TauFLSigma)		cut.at(TauFLSigma)=3.;//3
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
		else if(i_cut==DiMuonVeto){
			title.at(i_cut)="Di muon veto $\\Delta R <$";
	        char buffer[50];
	        sprintf(buffer,"%5.2f",cut.at(DiMuonVeto));
	    	title.at(i_cut)+=buffer;
			htitle=title.at(i_cut);
			htitle.ReplaceAll("$","");
			htitle.ReplaceAll("\\","#");
			hlabel="#Delta R(#mu_{veto}^{+},#mu_{veto}^{-})";
			Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DiMuonVeto_",htitle,100,-1.,4.,hlabel,"Events"));
			Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DiMuonVeto_",htitle,100,-1.,4,hlabel,"Events"));
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
			title.at(i_cut)="TauFLSigma $>=$";
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

	Tau_pt=HConfig.GetTH1D(Name+"_Tau_pt","Tau_pt",50,0,100,"vis. Tau p_{t}","Events");
	Tau_phi=HConfig.GetTH1D(Name+"_Tau_phi","Tau_phi",30,-3.14159265359,3.14159265359,"vis. Tau #phi","Events");
	Tau_eta=HConfig.GetTH1D(Name+"_Tau_eta","Tau_eta",50,-2.5,2.5,"vis. Tau #eta","Events");

	Tau_pt_wo_FLSigmaCut=HConfig.GetTH1D(Name+"_Tau_pt_wo_FLSigmaCut","Tau_pt_wo_FLSigmaCut",50,0,100,"vis. Tau p_{t} w/o #sigma_{FL} cut","Events");
	Tau_phi_wo_FLSigmaCut=HConfig.GetTH1D(Name+"_Tau_phi_wo_FLSigmaCut","Tau_phi_wo_FLSigmaCut",30,-3.14159265359,3.14159265359,"vis. Tau #phi w/o #sigma_{FL} cut","Events");
	Tau_eta_wo_FLSigmaCut=HConfig.GetTH1D(Name+"_Tau_eta_wo_FLSigmaCut","Tau_eta_wo_FLSigmaCut",50,-2.5,2.5,"vis. Tau #eta w/o #sigma_{FL} cut","Events");

	Tau_Mass_Inclusive=HConfig.GetTH1D(Name+"_Tau_Mass_Inclusive","Tau_Mass_Inclusive",170,-0.2,1.5,"PFTau mass","Events");
	Tau_Mass_sq_Inclusive=HConfig.GetTH1D(Name+"_Tau_Mass_sq_Inclusive","Tau_Mass_sq_Inclusive",170,-0.05,0.05,"PFTau mass squared","Events");
	Tau_Mass_Inclusive_UnFitTracks=HConfig.GetTH1D(Name+"_Tau_Mass_Inclusive_UnFitTracks","Tau_Mass_Inclusive_UnFitTracks",170,-0.2,1.5,"Tau_Mass_Inclusive_UnFitTracks","Events");
	Tau_Mass_Inclusive_ReFitTracks=HConfig.GetTH1D(Name+"_Tau_Mass_Inclusive_ReFitTracks","Tau_Mass_Inclusive_ReFitTracks",170,-0.2,1.5,"Tau_Mass_Inclusive_ReFitTracks","Events");

	Tau_Mass_Difference_PFTau_UnFitTracks_3PS=HConfig.GetTH1D(Name+"_Tau_Mass_Difference_PFTau_UnFitTracks_3PS","Tau_Mass_Difference_PFTau_UnFitTracks_3PS",50,-0.1,0.1,"Tau_Mass_Difference_PFTau_UnFitTracks_3PS","Events");
	Tau_Mass_Difference_PFTau_ReFitTracks_3PS=HConfig.GetTH1D(Name+"_Tau_Mass_Difference_PFTau_ReFitTracks_3PS","Tau_Mass_Difference_PFTau_ReFitTracks_3PS",50,-0.1,0.1,"Tau_Mass_Difference_PFTau_ReFitTracks_3PS","Events");

	MET_phi=HConfig.GetTH1D(Name+"_MET_phi","MET_phi",30,-3.14159265359,3.14159265359,"MET #phi","Events");

	TauFL_NoTauFLSigmaCut=HConfig.GetTH1D(Name+"_TauFL_NoTauFLSigmaCut","TauFL_NoTauFLSigmaCut",100,-2,4,"tau flight length without tauflsigma cut","Events");
	TauFLSigned_NoTauFLSigmaCut=HConfig.GetTH1D(Name+"_TauFLSigned_NoTauFLSigmaCut","TauFLSigned_NoTauFLSigmaCut",100,-2,4,"tau flight length Signed without tauflsigma cut","Events");

	TauFLSigmaSigned=HConfig.GetTH1D(Name+"_TauFLSigmaSigned","TauFLSigma Signed",80,-10,30,"TauFLSigma Signed","Events");
	TauFLSigmaUnsigned=HConfig.GetTH1D(Name+"_TauFLSigmaUnsigned","TauFLSigma Unsigned",80,-10,30,"TauFLSigma Unsigned","Events");

	A1mass=HConfig.GetTH1D(Name+"_A1mass","A1mass",100,0,2,"A1mass","Events");
	A1mass10GeV=HConfig.GetTH1D(Name+"_A1massGeV","A1massGeV",100,0,10,"A1mass","Events");
	A1massRefit=HConfig.GetTH1D(Name+"_A1massRefit","A1massRefit",100,0,2,"A1massRefit","Events");
	dA1mass_PFTau_Refit=HConfig.GetTH1D(Name+"_dA1mass_PFTau_Refit","dA1mass_PFTau_Refit",100,-0.05,0.05,"dA1mass_PFTau_Refit","Events");
	A1_Phi_Res=HConfig.GetTH1D(Name+"_A1_Phi_Res","A1_Phi_Res",100,-0.001,0.001,"A1_Phi_Res","Events");
	A1_Theta_Res=HConfig.GetTH1D(Name+"_A1_Theta_Res","A1_Theta_Res",100,-0.001,0.001,"A1_Theta_Res","Events");

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
	dPhi_GenTauMu_GenMu=HConfig.GetTH1D(Name+"_dPhi_GenTauMu_GenMu","dPhi_GenTauMu_GenMu",128,-3.14159265359/32,3.14159265359/32,"dPhi_GenTauMu_GenMu","Events");
	dPhi_GenTauMu_RecoMu=HConfig.GetTH1D(Name+"_dPhi_GenTauMu_RecoMu","dPhi_GenTauMu_RecoMu",128,-3.14159265359/32,3.14159265359/32,"dPhi_GenTauMu_RecoMu","Events");
	dTheta_GenTauMu_GenMu=HConfig.GetTH1D(Name+"_dTheta_GenTauMu_GenMu","dTheta_GenTauMu_GenMu",128,-3.14159265359/32,3.14159265359/32,"dTheta_GenTauMu_GenMu","Events");
	dTheta_GenTauMu_RecoMu=HConfig.GetTH1D(Name+"_dTheta_GenTauMu_RecoMu","dTheta_GenTauMu_RecoMu",128,-3.14159265359/32,3.14159265359/32,"dTheta_GenTauMu_RecoMu","Events");

	GJAngle_Over_GJAngleMax_StraightTau=HConfig.GetTH1D(Name+"_GJAngle_Over_GJAngleMax_StraightTau","GJAngle_Over_GJAngleMax_StraightTau",50,0.,5.,"GJAngle_Over_GJAngleMax_StraightTau","Events");
	GJAngle_Over_GJAngleMax_HelixTau=HConfig.GetTH1D(Name+"_GJAngle_Over_GJAngleMax_HelixTau","GJAngle_Over_GJAngleMax_HelixTau",50,0.,5.,"GJAngle_Over_GJAngleMax_HelixTau","Events");
	dGJAngle_GJAngleMAX_StraightTau=HConfig.GetTH1D(Name+"_dGJAngle_GJAngleMAX_StraightTau","dGJAngle_GJAngleMAX_StraightTau",50,-0.01,0.01,"dGJAngle_GJAngleMAX_StraightTau","Events");
	dGJAngle_GJAngleMAX_HelixTau=HConfig.GetTH1D(Name+"_dGJAngle_GJAngleMAX_HelixTau","dGJAngle_GJAngleMAX_HelixTau",50,-0.01,0.01,"dGJAngle_GJAngleMAX_HelixTau","Events");
	dGJAngle_HelixTau_StraightTau=HConfig.GetTH1D(Name+"_dGJAngle_HelixTau_StraightTau","dGJAngle_HelixTau_StraightTau",51,-0.0001,0.0001,"dGJAngle_HelixTau_StraightTau","Events");
	dGJAngle_HelixTau_StraightTauOverGJAngle=HConfig.GetTH1D(Name+"_dGJAngle_HelixTau_StraightTauOverGJAngle","dGJAngle_HelixTau_StraightTauOverGJAngle",51,-0.04,0.04,"dGJAngle_HelixTau_StraightTauOverGJAngle","Events");
	Angle_HelixTau_StraightTau=HConfig.GetTH1D(Name+"_Angle_HelixTau_StraightTau","Angle_HelixTau_StraightTau",100,-0.00005,0.0002,"Angle_HelixTau_StraightTau","Events");

	NUnphysical_StraightTau_HelixTau=HConfig.GetTH1D(Name+"_NUnphysical_StraightTau_HelixTau","NUnphysical_StraightTau_HelixTau",7,-2.5,4.5,"NUnphysical_StraightTau_HelixTau","Events");

	//dGJAngle_StraightTau_Gen
	//dGJAngle_HelixTau_Gen

	TauA1_Reco_Solution_StraightTau=HConfig.GetTH1D(Name+"_TauA1_Reco_Solution_StraightTau","TauA1_Reco_Solution_StraightTau",5,-2.5,2.5,"TauA1_Reco_Solution_StraightTau","Events");
	TauA1_Reco_Solution_HelixTau=HConfig.GetTH1D(Name+"_TauA1_Reco_Solution_HelixTau","TauA1_Reco_Solution_HelixTau",5,-2.5,2.5,"TauA1_Reco_Solution_HelixTau","Events");

	Gen_TauA1_GJ=HConfig.GetTH1D(Name+"_Gen_TauA1_GJ","Gen_TauA1_GJ",100,0,.05,"Gen_TauA1_GJ","Events");
	Gen_TauMu_GJ=HConfig.GetTH1D(Name+"_Gen_TauMu_GJ","Gen_TauMu_GJ",100,0,.05,"Gen_TauMu_GJ","Events");

	Gen_DiTau_dPhi=HConfig.GetTH1D(Name+"_Gen_DiTau_dPhi","Gen_DiTau_dPhi",128,-3.14159265359,3.14159265359,"Gen_DiTau_dPhi","Events");
	Gen_DiTau_Pt=HConfig.GetTH1D(Name+"_Gen_DiTau_Pt","Gen_DiTau_Pt",30,0,30,"Gen_DiTau_Pt","Events");

	Gen_DiTau_PtBalance_M=HConfig.GetTH1D(Name+"_Gen_DiTau_PtBalance_M","Gen_DiTau_PtBalance_M",280,20,160,"Gen_DiTau_PtBalance_M","Events");
	Gen_DiTau_PtBalance_dM=HConfig.GetTH1D(Name+"_Gen_DiTau_PtBalance_dM","Gen_DiTau_PtBalance_dM",100,-50,50,"Gen_DiTau_PtBalance_dM","Events");
	Gen_TauMu_PtBalance_Pt=HConfig.GetTH1D(Name+"_Gen_TauMu_PtBalance_Pt","Gen_TauMu_PtBalance_Pt",100,-50,50,"Gen_TauMu_PtBalance_Pt","Events");
	Gen_TauMu_PtBalance_dP=HConfig.GetTH1D(Name+"_Gen_TauMu_PtBalance_dP","Gen_TauMu_PtBalance_dP",100,-50,50,"Gen_TauMu_PtBalance_dP","Events");
	Gen_TauA1_dP=HConfig.GetTH1D(Name+"_Gen_TauA1_dP","Gen_TauA1_dP",50,0,50,"Gen_TauA1_dP","Events");

	Gen_TauA1_P=HConfig.GetTH1D(Name+"_Gen_TauA1_P","Gen_TauA1_P",75, 0,150,"Gen_TauA1_P","Events");
	Gen_TauA1_Pt=HConfig.GetTH1D(Name+"_Gen_TauA1_Pt","Gen_TauA1_Pt",75, 0,150,"Gen_TauA1_Pt","Events");
	Gen_TauA1_Px=HConfig.GetTH1D(Name+"_Gen_TauA1_Px","Gen_TauA1_Px",100,-100,100,"Gen_TauA1_Px","Events");
	Gen_TauA1_Py=HConfig.GetTH1D(Name+"_Gen_TauA1_Py","Gen_TauA1_Py",100,-100,100,"Gen_TauA1_Py","Events");
	Gen_TauA1_Pz=HConfig.GetTH1D(Name+"_Gen_TauA1_Pz","Gen_TauA1_Pz",100,-100,100,"Gen_TauA1_Pz","Events");
	Gen_TauA1_Phi=HConfig.GetTH1D(Name+"_Gen_TauA1_Phi","Gen_TauA1_Phi",64,-3.14159265359,3.14159265359,"Gen_TauA1_Phi","Events");
	Gen_TauA1_Eta=HConfig.GetTH1D(Name+"_Gen_TauA1_Eta","Gen_TauA1_Eta",50,-2.5, 2.5,"Gen_TauA1_Eta","Events");

	Gen_TauMu_P=HConfig.GetTH1D(Name+"_Gen_TauMu_P","Gen_TauMu_P",75, 0,150,"Gen_TauMu_P","Events");
	Gen_TauMu_Pt=HConfig.GetTH1D(Name+"_Gen_TauMu_Pt","Gen_TauMu_Pt",75, 0,150,"Gen_TauMu_Pt","Events");
	Gen_TauMu_Px=HConfig.GetTH1D(Name+"_Gen_TauMu_Px","Gen_TauMu_Px",100,-100,100,"Gen_TauMu_Px","Events");
	Gen_TauMu_Py=HConfig.GetTH1D(Name+"_Gen_TauMu_Py","Gen_TauMu_Py",100,-100,100,"Gen_TauMu_Py","Events");
	Gen_TauMu_Pz=HConfig.GetTH1D(Name+"_Gen_TauMu_Pz","Gen_TauMu_Pz",100,-100,100,"Gen_TauMu_Pz","Events");
	Gen_TauMu_Phi=HConfig.GetTH1D(Name+"_Gen_TauMu_Phi","Gen_TauMu_Phi",64,-3.14159265359,3.14159265359,"Gen_TauMu_Phi","Events");
	Gen_TauMu_Eta=HConfig.GetTH1D(Name+"_Gen_TauMu_Eta","Gen_TauMu_Eta",50,-2.5, 2.5,"Gen_TauMu_Eta","Events");

	Gen_Z_M=HConfig.GetTH1D(Name+"_Gen_Z_M","Gen_Z_M",180,60,150,"Gen_Z_M","Events");
	Gen_Z_Pt=HConfig.GetTH1D(Name+"_Gen_Z_Pt","Gen_Z_Pt",100,0,100,"Gen_Z_Pt","Events");
	Gen_Z_Phi=HConfig.GetTH1D(Name+"_Gen_Z_Phi","Gen_Z_Phi",64,-3.14159265359,3.14159265359,"Gen_Z_Phi","Events");
	Gen_Z_Eta=HConfig.GetTH1D(Name+"_Gen_Z_Eta","Gen_Z_Eta",50,-10,10,"Gen_Z_Eta","Events");

	Gen_TauA1_P_noSel=HConfig.GetTH1D(Name+"_Gen_TauA1_P_noSel","Gen_TauA1_P_noSel",75, 0,150,"Gen_TauA1_P_noSel","Events");
	Gen_TauA1_Pt_noSel=HConfig.GetTH1D(Name+"_Gen_TauA1_Pt_noSel","Gen_TauA1_Pt_noSel",75, 0,150,"Gen_TauA1_Pt_noSel","Events");
	Gen_TauA1_Px_noSel=HConfig.GetTH1D(Name+"_Gen_TauA1_Px_noSel","Gen_TauA1_Px_noSel",100,-100,100,"Gen_TauA1_Px_noSel","Events");
	Gen_TauA1_Py_noSel=HConfig.GetTH1D(Name+"_Gen_TauA1_Py_noSel","Gen_TauA1_Py_noSel",100,-100,100,"Gen_TauA1_Py_noSel","Events");
	Gen_TauA1_Pz_noSel=HConfig.GetTH1D(Name+"_Gen_TauA1_Pz_noSel","Gen_TauA1_Pz_noSel",100,-100,100,"Gen_TauA1_Pz_noSel","Events");
	Gen_TauA1_Phi_noSel=HConfig.GetTH1D(Name+"_Gen_TauA1_Phi_noSel","Gen_TauA1_Phi_noSel",64,-3.14159265359,3.14159265359,"Gen_TauA1_Phi_noSel","Events");
	Gen_TauA1_Eta_noSel=HConfig.GetTH1D(Name+"_Gen_TauA1_Eta_noSel","Gen_TauA1_Eta_noSel",50,-2.5, 2.5,"Gen_TauA1_Eta_noSel","Events");

	Gen_TauMu_P_noSel=HConfig.GetTH1D(Name+"_Gen_TauMu_P_noSel","Gen_TauMu_P_noSel",75, 0,150,"Gen_TauMu_P_noSel","Events");
	Gen_TauMu_Pt_noSel=HConfig.GetTH1D(Name+"_Gen_TauMu_Pt_noSel","Gen_TauMu_Pt_noSel",75, 0,150,"Gen_TauMu_Pt_noSel","Events");
	Gen_TauMu_Px_noSel=HConfig.GetTH1D(Name+"_Gen_TauMu_Px_noSel","Gen_TauMu_Px_noSel",100,-100,100,"Gen_TauMu_Px_noSel","Events");
	Gen_TauMu_Py_noSel=HConfig.GetTH1D(Name+"_Gen_TauMu_Py_noSel","Gen_TauMu_Py_noSel",100,-100,100,"Gen_TauMu_Py_noSel","Events");
	Gen_TauMu_Pz_noSel=HConfig.GetTH1D(Name+"_Gen_TauMu_Pz_noSel","Gen_TauMu_Pz_noSel",100,-100,100,"Gen_TauMu_Pz_noSel","Events");
	Gen_TauMu_Phi_noSel=HConfig.GetTH1D(Name+"_Gen_TauMu_Phi_noSel","Gen_TauMu_Phi_noSel",64,-3.14159265359,3.14159265359,"Gen_TauMu_Phi_noSel","Events");
	Gen_TauMu_Eta_noSel=HConfig.GetTH1D(Name+"_Gen_TauMu_Eta_noSel","Gen_TauMu_Eta_noSel",50,-2.5, 2.5,"Gen_TauMu_Eta_noSel","Events");

	Gen_Z_Pt_noSel=HConfig.GetTH1D(Name+"_Gen_Z_Pt_noSel","Gen_Z_Pt_noSel",100,0,100,"Gen_Z_Pt_noSel","Events");
	Gen_Z_Phi_noSel=HConfig.GetTH1D(Name+"_Gen_Z_Phi_noSel","Gen_Z_Phi_noSel",64,-3.14159265359,3.14159265359,"Gen_Z_Phi_noSel","Events");
	Gen_Z_Eta_noSel=HConfig.GetTH1D(Name+"_Gen_Z_Eta_noSel","Gen_Z_Eta_noSel",50,-10,10,"Gen_Z_Eta_noSel","Events");


	Gen_TPTF_TauA1_Solution_NoSelection=HConfig.GetTH1D(Name+"_Gen_TPTF_TauA1_RightSolution_NoSelection","Gen_TPTF_TauA1_RightSolution_NoSelection",3,-0.5,2.5,"Gen_TPTF_TauA1_RightSolution_NoSelection","Events");
	Gen_TPTF_TauA1_Solution_WithSelection=HConfig.GetTH1D(Name+"_Gen_TPTF_TauA1_RightSolution_WithSelection","Gen_TPTF_TauA1_RightSolution_WithSelection",3,-0.5,2.5,"Gen_TPTF_TauA1_RightSolution_WithSelection","Events");

	Gen_Z_Pt_vs_MET=HConfig.GetTH2D(Name+"_Gen_Z_Pt_vs_MET","Gen_Z_Pt_vs_MET",120,0,30,120,0,30,"Gen_Z_Pt","MET");
	Gen_Z_Pt_vs_VtxTracksPt=HConfig.GetTH2D(Name+"_Gen_Z_Pt_vs_VtxTracksPt","Gen_Z_Pt_vs_VtxTracksPt",30,0,30,30,0,30,"Gen_Z_Pt","VtxTracksPt");
	Gen_Z_Phi_vs_VtxTracksPhi=HConfig.GetTH2D(Name+"_Gen_Z_Phi_vs_VtxTracksPhi","Gen_Z_Phi_vs_VtxTracksPhi",32,-3.14159265359,3.14159265359,32,-3.14159265359,3.14159265359,"Gen_Z_Phi","VtxTracksPhi");

	VtxTracksPtRes=HConfig.GetTH1D(Name+"_VtxTracksPtRes","VtxTracksPtRes",120,-30,30,"VtxTracksPtRes","Events");
	VtxTracksPhiCorrectedRes=HConfig.GetTH1D(Name+"_VtxTracksPhiCorrectedRes","VtxTracksPhiCorrectedRes",128,2*-3.14159265359,2*3.14159265359,"VtxTracksPhiCorrectedRes","Events");

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
	TPTF_Neutrino_RefitPFTau_HelixGenTau_Mass=HConfig.GetTH1D(Name+"_TPTF_Neutrino_RefitPFTau_HelixGenTau_Mass","TPTF_Neutrino_RefitPFTau_HelixGenTau_Mass",50,-1,1,"TPTF_Neutrino_RefitPFTau_HelixGenTau_Mass","Events");
	TPTF_Neutrino_GenA1_StraightGenTau_Mass=HConfig.GetTH1D(Name+"_TPTF_Neutrino_GenA1_StraightGenTau_Mass","TPTF_Neutrino_GenA1_StraightGenTau_Mass",50,-0.04,0.04,"TPTF_Neutrino_GenA1_StraightGenTau_Mass","Events");
	TPTF_Neutrino_GenA1_HelixGenTau_Mass=HConfig.GetTH1D(Name+"_TPTF_Neutrino_GenA1_HelixGenTau_Mass","TPTF_Neutrino_GenA1_HelixGenTau_Mass",50,-0.2,0.2,"TPTF_Neutrino_GenA1_HelixGenTau_Mass","Events");

	TransTrk_Failure_withSelection=HConfig.GetTH1D(Name+"_TransTrk_Failure_withSelection","TransTrk_Failure_withSelection",2,-0.5,1.5,"TransTrk_Failure_withSelection","Events");
	TransTrk_Failure_noSelection=HConfig.GetTH1D(Name+"_TransTrk_Failure_noSelection","TransTrk_Failure_noSelection",2,-0.5,1.5,"TransTrk_Failure_noSelection","Events");

	Est_Z_Pt_wTruth=HConfig.GetTH1D(Name+"_Est_Z_Pt_wTruth","Est_Z_Pt_wTruth",40,0,40,"Est_Z_Pt_wTruth","Events");
	Est_Z_PtRes_wTruth=HConfig.GetTH1D(Name+"_Est_Z_PtRes_wTruth","Est_Z_PtRes_wTruth",80,-40,40,"Est_Z_PtRes_wTruth","Events");
	Est_Z_Pt_wTruth_vs_GenZ_Pt=HConfig.GetTH2D(Name+"_Est_Z_Pt_wTruth_vs_GenZ_Pt","Est_Z_Pt_wTruth_vs_GenZ_Pt",40,0,40,40,0,40,"Est_Z_Pt_wTruth","GenZ_Pt");
	Est_TauMu_PtRes_wTruth=HConfig.GetTH1D(Name+"_Est_TauMu_PtRes_wTruth","Est_TauMu_PtRes_wTruth",100,-50,50,"Est_TauMu_PtRes_wTruth","Events");
	Est_Z_Energy_wTruth=HConfig.GetTH1D(Name+"_Est_Z_Energy_wTruth","Est_Z_Energy_wTruth",100,0,1000,"Est_Z_Energy_wTruth","Events");
	Est_Z_EnergyRes_wTruth=HConfig.GetTH1D(Name+"_Est_Z_EnergyRes_wTruth","Est_Z_EnergyRes_wTruth",100,-100,100,"Est_Z_EnergyRes_wTruth","Events");
	Est_Z_Pt_alwaysMinus=HConfig.GetTH1D(Name+"_Est_Z_Pt_alwaysMinus","Est_Z_Pt_alwaysMinus",40,0,40,"Est_Z_Pt_alwaysMinus","Events");
	Est_Z_PtRes_alwaysMinus=HConfig.GetTH1D(Name+"_Est_Z_PtRes_alwaysMinus","Est_Z_PtRes_alwaysMinus",80,-40,40,"Est_Z_PtRes_alwaysMinus","Events");
	Est_Z_Pt_alwaysMinus_vs_GenZ_Pt=HConfig.GetTH2D(Name+"_Est_Z_Pt_alwaysMinus_vs_GenZ_Pt","Est_Z_Pt_alwaysMinus_vs_GenZ_Pt",40,0,40,40,0,40,"Est_Z_Pt_alwaysMinus","GenZ_Pt");
	Est_Z_Energy_alwaysMinus=HConfig.GetTH1D(Name+"_Est_Z_Energy_alwaysMinus","Est_Z_Energy_alwaysMinus",100,0,1000,"Est_Z_Energy_alwaysMinus","Events");
	Est_Z_EnergyRes_alwaysMinus=HConfig.GetTH1D(Name+"_Est_Z_EnergyRes_alwaysMinus","Est_Z_EnergyRes_alwaysMinus",100,-100,100,"Est_Z_EnergyRes_alwaysMinus","Events");

	Est_TauMu_PtRes_wTruth2=HConfig.GetTH1D(Name+"_Est_TauMu_PtRes_wTruth2","Est_TauMu_PtRes_wTruth2",100,-50,50,"Est_TauMu_PtRes_wTruth2","Events");
	Est_Z_EnergyRes_wTruth2=HConfig.GetTH1D(Name+"_Est_Z_EnergyRes_wTruth2","Est_Z_EnergyRes_wTruth2",100,-100,100,"Est_Z_EnergyRes_wTruth2","Events");
	Est_Z_EnergyRes_alwaysMinus2=HConfig.GetTH1D(Name+"_Est_Z_EnergyRes_alwaysMinus2","Est_Z_EnergyRes_alwaysMinus2",100,-100,100,"Est_Z_EnergyRes_alwaysMinus2","Events");
	Est_Z_M_wTruth2=HConfig.GetTH1D(Name+"_Est_Z_M_wTruth2","Est_Z_M_wTruth2",100,0,300,"Est_Z_M_wTruth2","Events");

	Est_TauMu_PtRes_wTruth_NoZMass=HConfig.GetTH1D(Name+"_Est_TauMu_PtRes_wTruth_NoZMass","Est_TauMu_PtRes_wTruth_NoZMass",100,-100,100,"Est_TauMu_PtRes_wTruth_NoZMass","Events");
	Est_Z_M_wTruth_NoZMass=HConfig.GetTH1D(Name+"_Est_Z_M_wTruth_NoZMass","Est_Z_M_wTruth_NoZMass",100,0,300,"Est_Z_M_wTruth_NoZMass","Events");
	Est_Z_EnergyRes_wTruth_NoZMass=HConfig.GetTH1D(Name+"_Est_Z_EnergyRes_wTruth_NoZMass","Est_Z_EnergyRes_wTruth_NoZMass",100,-100,100,"Est_Z_EnergyRes_wTruth_NoZMass","Events");

	Est_TauMu_wMET_PtRes=HConfig.GetTH1D(Name+"_Est_TauMu_wMET_PtRes","Est_TauMu_wMET_PtRes",100,-40,40,"Est_TauMu_wMET_PtRes","Events");
	Est_TauMu_wMET_PhiRes=HConfig.GetTH1D(Name+"_Est_TauMu_wMET_PhiRes","Est_TauMu_wMET_PhiRes",100,-7,7,"Est_TauMu_wMET_PhiRes","Events");
	Est_TauMu_wMET_EtaRes=HConfig.GetTH1D(Name+"_Est_TauMu_wMET_EtaRes","Est_TauMu_wMET_EtaRes",100,-5,5,"Est_TauMu_wMET_EtaRes","Events");

	Est_Z_wMET_PtRes=HConfig.GetTH1D(Name+"_Est_Z_wMET_PtRes","Est_Z_wMET_PtRes",100,-40,40,"Est_Z_wMET_PtRes","Events");
	Est_Z_wMET_PhiRes=HConfig.GetTH1D(Name+"_Est_Z_wMET_PhiRes","Est_Z_wMET_PhiRes",64,-3.14159265359,3.14159265359,"Est_Z_wMET_PhiRes","Events");

	//DiTau Reco

	Reco_ZMass=HConfig.GetTH1D(Name+"_Reco_ZMass","Reco_ZMass",180,60,150,"Reco_ZMass","Events");
	Reco_ZMass_UnboostedGenZ=HConfig.GetTH1D(Name+"_Reco_ZMass_UnboostedGenZ","Reco_ZMass_UnboostedGenZ",180,60,150,"Reco_ZMass_UnboostedGenZ","Events");
	Reco_EventFit_Solution=HConfig.GetTH1D(Name+"_Reco_EventFit_Solution","Reco_EventFit_Solution",5,-2.5,2.5,"Reco_EventFit_Solution","Events");
	Reco_A1Fit_Solution=HConfig.GetTH1D(Name+"_Reco_A1Fit_Solution","Reco_A1Fit_Solution",5,-2.5,2.5,"Reco_A1Fit_Solution","Events");
	Reco_Chi2_FitSolutionOnlyLargeScale=HConfig.GetTH1D(Name+"_Reco_Chi2_FitSolutionOnlyLargeScale","Reco_Chi2_FitSolutionOnlyLargeScale",100,0,500000,"Reco_Chi2_FitSolutionOnlyLargeScale","Events");
	Reco_ConstrainedDeltaSum=HConfig.GetTH1D(Name+"_Reco_ConstrainedDeltaSum","Reco_ConstrainedDeltaSum",150,0,300,"Reco_ConstrainedDeltaSum","Events");
	Reco_ConstrainedDeltaMass=HConfig.GetTH1D(Name+"_Reco_ConstrainedDeltaMass","Reco_ConstrainedDeltaMass",100,0,40,"Reco_ConstrainedDeltaMass","Events");
	Reco_ConstrainedDeltaPt=HConfig.GetTH1D(Name+"_Reco_ConstrainedDeltaPt","Reco_ConstrainedDeltaPt",100,0,1,"Reco_ConstrainedDeltaPt","Events");
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

	Reco_dPhi_TauMuTauA1_AfterFit=HConfig.GetTH1D(Name+"_Reco_dPhi_TauMuTauA1_AfterFit","Reco_dPhi_TauMuTauA1_AfterFit",64,-3.14159265359,3.14159265359,"Reco_dPhi_TauMuTauA1_AfterFit","Events");
	Reco_dPhi_TauMuTauA1_BeforeFit=HConfig.GetTH1D(Name+"_Reco_dPhi_TauMuTauA1_BeforeFit","Reco_dPhi_TauMuTauA1_BeforeFit",64,-3.14159265359,3.14159265359,"Reco_dPhi_TauMuTauA1_BeforeFit","Events");

	Reco_dPhi_TauMuTauA1_AfterFit_lowBoost=HConfig.GetTH1D(Name+"_Reco_dPhi_TauMuTauA1_AfterFit_lowBoost","Reco_dPhi_TauMuTauA1_AfterFit_lowBoost",64,-3.14159265359,3.14159265359,"Reco_dPhi_TauMuTauA1_AfterFit_lowBoost","Events");
	Reco_dPhi_TauMuTauA1_BeforeFit_lowBoost=HConfig.GetTH1D(Name+"_Reco_dPhi_TauMuTauA1_BeforeFit_lowBoost","Reco_dPhi_TauMuTauA1_BeforeFit_lowBoost",64,-3.14159265359,3.14159265359,"Reco_dPhi_TauMuTauA1_BeforeFit_lowBoost","Events");

	Reco_PtRes_TauA1_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_NoFit","Reco_PtRes_TauA1_NoFit",100,-50,50,"Reco_PtRes_TauA1_NoFit","Events");
	Reco_PtRes_TauA1_AmbPoint0_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_AmbPoint0_NoFit","Reco_PtRes_TauA1_AmbPoint0_NoFit",100,-50,50,"Reco_PtRes_TauA1_AmbPoint0_NoFit","Events");
	Reco_PtRes_TauA1_AmbPoint12_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_AmbPoint12_NoFit","Reco_PtRes_TauA1_AmbPoint12_NoFit",100,-50,50,"Reco_PtRes_TauA1_AmbPoint12_NoFit","Events");

	Reco_PtRes_TauMu_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_NoFit","Reco_PtRes_TauMu_NoFit",100,-50,50,"Reco_PtRes_TauMu_NoFit","Events");
	Reco_PtRes_TauMu_AmbPoint0_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_AmbPoint0_NoFit","Reco_PtRes_TauMu_AmbPoint0_NoFit",100,-50,50,"Reco_PtRes_TauMu_AmbPoint0_NoFit","Events");
	Reco_PtRes_TauMu_AmbPoint12_NoFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_AmbPoint12_NoFit","Reco_PtRes_TauMu_AmbPoint12_NoFit",100,-50,50,"Reco_PtRes_TauMu_AmbPoint12_NoFit","Events");
	Reco_PhiRes_TauMu_PreFit=HConfig.GetTH1D(Name+"_Reco_PhiRes_TauMu_PreFit","Reco_PhiRes_TauMu_PreFit",64,-2*3.14159265359,2*3.14159265359,"Reco_PhiRes_TauMu_PreFit","Events");
	Reco_ThetaRes_TauMu_PreFit=HConfig.GetTH1D(Name+"_Reco_ThetaRes_TauMu_PreFit","Reco_ThetaRes_TauMu_PreFit",64,-2*3.14159265359,2*3.14159265359,"Reco_ThetaRes_TauMu_PreFit","Events");

	Reco_TauA1_P=HConfig.GetTH1D(Name+"_Reco_TauA1_P","Reco_TauA1_P",75, 0,150,"Reco_TauA1_P","Events");
	Reco_TauA1_Pt=HConfig.GetTH1D(Name+"_Reco_TauA1_Pt","Reco_TauA1_Pt",75, 0,150,"Reco_TauA1_Pt","Events");
	Reco_TauA1_Px=HConfig.GetTH1D(Name+"_Reco_TauA1_Px","Reco_TauA1_Px",100,-100,100,"Reco_TauA1_Px","Events");
	Reco_TauA1_Py=HConfig.GetTH1D(Name+"_Reco_TauA1_Py","Reco_TauA1_Py",100,-100,100,"Reco_TauA1_Py","Events");
	Reco_TauA1_Pz=HConfig.GetTH1D(Name+"_Reco_TauA1_Pz","Reco_TauA1_Pz",100,-100,100,"Reco_TauA1_Pz","Events");
	Reco_TauA1_Phi=HConfig.GetTH1D(Name+"_Reco_TauA1_Phi","Reco_TauA1_Phi",64,-3.14159265359,3.14159265359,"Reco_TauA1_Phi","Events");
	Reco_TauA1_Eta=HConfig.GetTH1D(Name+"_Reco_TauA1_eta","Reco_TauA1_Eta",50,-2.5,2.5,"Reco_TauA1_Eta","Events");

	Reco_TauMu_P=HConfig.GetTH1D(Name+"_Reco_TauMu_P","Reco_TauMu_P",75, 0,150,"Reco_TauMu_P","Events");
	Reco_TauMu_Pt=HConfig.GetTH1D(Name+"_Reco_TauMu_Pt","Reco_TauMu_Pt",75, 0,150,"Reco_TauMu_Pt","Events");
	Reco_TauMu_Px=HConfig.GetTH1D(Name+"_Reco_TauMu_Px","Reco_TauMu_Px",100,-100,100,"Reco_TauMu_Px","Events");
	Reco_TauMu_Py=HConfig.GetTH1D(Name+"_Reco_TauMu_Py","Reco_TauMu_Py",100,-100,100,"Reco_TauMu_Py","Events");
	Reco_TauMu_Pz=HConfig.GetTH1D(Name+"_Reco_TauMu_Pz","Reco_TauMu_Pz",100,-100,100,"Reco_TauMu_Pz","Events");
	Reco_TauMu_Phi=HConfig.GetTH1D(Name+"_Reco_TauMu_Phi","Reco_TauMu_Phi",64,-3.14159265359,3.14159265359,"Reco_TauMu_Phi","Events");
	Reco_TauMu_Eta=HConfig.GetTH1D(Name+"_Reco_TauMu_eta","Reco_TauMu_Eta",50,-2.5,2.5,"Reco_TauMu_Eta","Events");

	Reco_Z_Px=HConfig.GetTH1D(Name+"_Reco_Z_Px","Reco_Z_Px",100,-100,100,"Reco_Z_Px","Events");
	Reco_Z_Py=HConfig.GetTH1D(Name+"_Reco_Z_Py","Reco_Z_Py",100,-100,100,"Reco_Z_Py","Events");
	Reco_Z_Pz=HConfig.GetTH1D(Name+"_Reco_Z_Pz","Reco_Z_Pz",100,-100,100,"Reco_Z_Pz","Events");
	Reco_Z_Pt=HConfig.GetTH1D(Name+"_Reco_Z_Pt","Reco_Z_Pt",100,0,100,"Reco_Z_Pt","Events");
	Reco_Z_PtRes=HConfig.GetTH1D(Name+"_Reco_Z_PtRes","Reco_Z_PtRes",100,-100,100,"Reco_Z_PtRes","Events");
	Reco_Z_Phi=HConfig.GetTH1D(Name+"_Reco_Z_Phi","Reco_Z_Phi",64,-3.14159265359,3.14159265359,"Reco_Z_Phi","Events");
	Reco_Z_Eta=HConfig.GetTH1D(Name+"_Reco_Z_Eta","Reco_Z_Eta",100,-10,10,"Reco_Z_Eta","Events");

	Reco_ZMass_MassScan=HConfig.GetTH1D(Name+"_Reco_ZMass_MassScan","Reco_ZMass_MassScan",125,0,250,"Reco_ZMass_MassScan","Events");
	Reco_ZMass_MassScanUnboosted=HConfig.GetTH1D(Name+"_Reco_ZMass_MassScanUnboosted","Reco_ZMass_MassScanUnboosted",50,0,250,"Reco_ZMass_MassScanUnboosted","Events");
	Reco_ZMasswithProbWeight_MassScan=HConfig.GetTH1D(Name+"_Reco_ZMasswithProbWeight_MassScan","Reco_ZMasswithProbWeight_MassScan",50,0,250,"Reco_ZMasswithProbWeight_MassScan","Events");
	Reco_ProbStack_MassScan=HConfig.GetTH1D(Name+"_Reco_ProbStack_MassScan","Reco_ProbStack_MassScan",50,0,1,"Reco_ProbStack_MassScan","Events");
	Reco_ZMass_PDF=HConfig.GetTH1D(Name+"_Reco_ZMass_PDF","Reco_ZMass_PDF",50,0,250,"Reco_ZMass_PDF","Events");
	GenZ_Pt_Unboosted=HConfig.GetTH1D(Name+"_GenZ_Pt_Unboosted","GenZ_Pt_Unboosted",40,0,40,"GenZ_Pt_Unboosted","Events");
	RecoZ_Pt=HConfig.GetTH1D(Name+"_RecoZ_Pt","RecoZ_Pt",40,0,40,"RecoZ_Pt","Events");
	RecoZ_Pt_Unboosted=HConfig.GetTH1D(Name+"_RecoZ_Pt_Unboosted","RecoZ_Pt_Unboosted",40,0,40,"RecoZ_Pt_Unboosted","Events");

	Reco_TauMu_ResCosTheta=HConfig.GetTH1D(Name+"_RecoTaumu_ResCosTheta","RecoTaumu_ResCosTheta",31,-2,2,"RecoTaumu_ResCosTheta","Events");
	Reco_TauMu_DeltaPhi_FitImpact=HConfig.GetTH1D(Name+"_Reco_TauMu_DeltaPhi_FitImpact","Reco_TauMu_DeltaPhi_FitImpact",31,-3.14159265359*2,3.14159265359*2,"Reco_TauMu_DeltaPhi_FitImpact","Events");

	Reco_Z_PhiRes=HConfig.GetTH1D(Name+"_Reco_Z_PhiRes","Reco_Z_PhiRes",64,-3.14159265359,3.14159265359,"Reco_Z_PhiRes","Events");
	Reco_Z_PhiRes_noAmb=HConfig.GetTH1D(Name+"_Reco_Z_PhiRes_noAmb","Reco_Z_PhiRes_noAmb",64,-3.14159265359,3.14159265359,"Reco_Z_PhiRes_noAmb","Events");
	Reco_Z_PhiRes_wAmb=HConfig.GetTH1D(Name+"_Reco_Z_PhiRes_wAmb","Reco_Z_PhiRes_wAmb",64,-3.14159265359,3.14159265359,"Reco_Z_PhiRes_wAmb","Events");
	Reco_Z_PhiRes_pickedrightAmb=HConfig.GetTH1D(Name+"_Reco_Z_PhiRes_pickedrightAmb","Reco_Z_PhiRes_pickedrightAmb",64,-3.14159265359,3.14159265359,"Reco_Z_PhiRes_pickedrightAmb","Events");
	Reco_Z_PhiRes_pickedwrongAmb=HConfig.GetTH1D(Name+"_Reco_Z_PhiRes_pickedwrongAmb","Reco_Z_PhiRes_pickedwrongAmb",64,-3.14159265359,3.14159265359,"Reco_Z_PhiRes_pickedwrongAmb","Events");
	Reco_Z_EtaRes=HConfig.GetTH1D(Name+"_Reco_Z_EtaRes","Reco_Z_EtaRes",100,-10,10,"Reco_Z_EtaRes","Events");
	Reco_Z_EtaRes_noAmb=HConfig.GetTH1D(Name+"_Reco_Z_EtaRes_noAmb","Reco_Z_EtaRes_noAmb",100,-10,10,"Reco_Z_EtaRes_noAmb","Events");
	Reco_Z_EtaRes_wAmb=HConfig.GetTH1D(Name+"_Reco_Z_EtaRes_wAmb","Reco_Z_EtaRes_wAmb",100,-10,10,"Reco_Z_EtaRes_wAmb","Events");
	Reco_Z_EtaRes_pickedrightAmb=HConfig.GetTH1D(Name+"_Reco_Z_EtaRes_pickedrightAmb","Reco_Z_EtaRes_pickedrightAmb",100,-10,10,"Reco_Z_EtaRes_pickedrightAmb","Events");
	Reco_Z_EtaRes_pickedwrongAmb=HConfig.GetTH1D(Name+"_Reco_Z_EtaRes_pickedwrongAmb","Reco_Z_EtaRes_pickedwrongAmb",100,-10,10,"Reco_Z_EtaRes_pickedwrongAmb","Events");
	Reco_Z_PRes=HConfig.GetTH1D(Name+"_Reco_Z_PRes","Reco_Z_PRes",100,-100,100,"Reco_Z_PRes","Events");
	Reco_Z_PRes_noAmb=HConfig.GetTH1D(Name+"_Reco_Z_PRes_noAmb","Reco_Z_PRes_noAmb",100,-100,100,"Reco_Z_PRes_noAmb","Events");
	Reco_Z_PRes_wAmb=HConfig.GetTH1D(Name+"_Reco_Z_PRes_wAmb","Reco_Z_PRes_wAmb",100,-100,100,"Reco_Z_PRes_wAmb","Events");
	Reco_Z_PRes_pickedrightAmb=HConfig.GetTH1D(Name+"_Reco_Z_PRes_pickedrightAmb","Reco_Z_PRes_pickedrightAmb",100,-100,100,"Reco_Z_PRes_pickedrightAmb","Events");
	Reco_Z_PRes_pickedwrongAmb=HConfig.GetTH1D(Name+"_Reco_Z_PRes_pickedwrongAmb","Reco_Z_PRes_pickedwrongAmb",100,-100,100,"Reco_Z_PRes_pickedwrongAmb","Events");
	Reco_NIter=HConfig.GetTH1D(Name+"_Reco_NIter","Reco_NIter",100,0,100,"Reco_NIter all ambiguities","Events");
	Reco_NIter_noAmb=HConfig.GetTH1D(Name+"_Reco_NIter_noAmb","Reco_NIter_noAmb",100,0,100,"Reco_NIter_noAmb","Events");
	Reco_NIter_wAmb=HConfig.GetTH1D(Name+"_Reco_NIter_wAmb","Reco_NIter_wAmb",100,0,100,"Reco_NIter_wAmb","Events");
	Reco_NIter_pickedrightAmb=HConfig.GetTH1D(Name+"_Reco_NIter_pickedrightAmb","Reco_NIter_pickedrightAmb",100,0,100,"Reco_NIter_pickedrightAmb","Events");
	Reco_NIter_pickedwrongAmb=HConfig.GetTH1D(Name+"_Reco_NIter_pickedwrongAmb","Reco_NIter_pickedwrongAmb",100,0,100,"Reco_NIter_pickedwrongAmb","Events");
	Reco_Chi2=HConfig.GetTH1D(Name+"_Reco_Chi2","Reco_Chi2",100,0,1,"Reco_Chi2","Events");
	Reco_Chi2_noAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_noAmb","Reco_Chi2_noAmb",100,0,1,"Reco_Chi2_noAmb","Events");
	Reco_Chi2_wAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_wAmb","Reco_Chi2_wAmb",100,0,1,"Reco_Chi2_wAmb","Events");
	Reco_Chi2_rightAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_rightAmb","Reco_Chi2_rightAmb",100,0,1,"Reco_Chi2_rightAmb","Events");
	Reco_Chi2_wrongAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_wrongAmb","Reco_Chi2_wrongAmb",100,0,1,"Reco_Chi2_wrongAmb","Events");
	Reco_Chi2_pickedrightAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_pickedrightAmb","Reco_Chi2_pickedrightAmb",100,0,1,"Reco_Chi2_pickedrightAmb","Events");
	Reco_Chi2_pickedwrongAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_pickedwrongAmb","Reco_Chi2_pickedwrongAmb",100,0,1,"Reco_Chi2_pickedwrongAmb","Events");
	Reco_Chi2_orig=HConfig.GetTH1D(Name+"_Reco_Chi2_orig","Reco_Chi2_orig",100,0,1,"Reco_Chi2_orig","Events");
	Reco_Chi2_orig_noAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_orig_noAmb","Reco_Chi2_orig_noAmb",100,0,1,"Reco_Chi2_orig_noAmb","Events");
	Reco_Chi2_orig_wAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_orig_wAmb","Reco_Chi2_orig_wAmb",100,0,1,"Reco_Chi2_orig_wAmb","Events");
	Reco_Chi2_orig_rightAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_orig_rightAmb","Reco_Chi2_orig_rightAmb",100,0,1,"Reco_Chi2_orig_rightAmb","Events");
	Reco_Chi2_orig_wrongAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_orig_wrongAmb","Reco_Chi2_orig_wrongAmb",100,0,1,"Reco_Chi2_orig_wrongAmb","Events");
	Reco_Chi2_orig_pickedrightAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_orig_pickedrightAmb","Reco_Chi2_orig_pickedrightAmb",100,0,1,"Reco_Chi2_orig_pickedrightAmb","Events");
	Reco_Chi2_orig_pickedwrongAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_orig_pickedwrongAmb","Reco_Chi2_orig_pickedwrongAmb",100,0,1,"Reco_Chi2_orig_pickedwrongAmb","Events");
	Reco_Chi2_SC=HConfig.GetTH1D(Name+"_Reco_Chi2_SC","Reco_Chi2_SC",100,0,1,"Reco_Chi2_SC","Events");
	Reco_Chi2_SC_noAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_SC_noAmb","Reco_Chi2_SC_noAmb",100,0,1,"Reco_Chi2_SC_noAmb","Events");
	Reco_Chi2_SC_wAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_SC_wAmb","Reco_Chi2_SC_wAmb",100,0,1,"Reco_Chi2_SC_wAmb","Events");
	Reco_Chi2_SC_rightAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_SC_rightAmb","Reco_Chi2_SC_rightAmb",100,0,1,"Reco_Chi2_SC_rightAmb","Events");
	Reco_Chi2_SC_wrongAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_SC_wrongAmb","Reco_Chi2_SC_wrongAmb",100,0,1,"Reco_Chi2_SC_wrongAmb","Events");
	Reco_Chi2_SC_pickedrightAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_SC_pickedrightAmb","Reco_Chi2_SC_pickedrightAmb",100,0,1,"Reco_Chi2_SC_pickedrightAmb","Events");
	Reco_Chi2_SC_pickedwrongAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_SC_pickedwrongAmb","Reco_Chi2_SC_pickedwrongAmb",100,0,1,"Reco_Chi2_SC_pickedwrongAmb","Events");
	Reco_Chi2_HC=HConfig.GetTH1D(Name+"_Reco_Chi2_HC","Reco_Chi2_HC",100,0,1,"Reco_Chi2_HC","Events");
	Reco_Chi2_HC_noAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_HC_noAmb","Reco_Chi2_HC_noAmb",100,0,1,"Reco_Chi2_HC_noAmb","Events");
	Reco_Chi2_HC_wAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_HC_wAmb","Reco_Chi2_HC_wAmb",100,0,1,"Reco_Chi2_HC_wAmb","Events");
	Reco_Chi2_HC_rightAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_HC_rightAmb","Reco_Chi2_HC_rightAmb",100,0,1,"Reco_Chi2_HC_rightAmb","Events");
	Reco_Chi2_HC_wrongAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_HC_wrongAmb","Reco_Chi2_HC_wrongAmb",100,0,1,"Reco_Chi2_HC_wrongAmb","Events");
	Reco_Chi2_HC_pickedrightAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_HC_pickedrightAmb","Reco_Chi2_HC_pickedrightAmb",100,0,1,"Reco_Chi2_HC_pickedrightAmb","Events");
	Reco_Chi2_HC_pickedwrongAmb=HConfig.GetTH1D(Name+"_Reco_Chi2_HC_pickedwrongAmb","Reco_Chi2_HC_pickedwrongAmb",100,0,1,"Reco_Chi2_HC_pickedwrongAmb","Events");

	Reco_Chi2_diff=HConfig.GetTH1D(Name+"_Reco_Chi2_diff","Reco_Chi2_diff",100,-2,2,"Reco_Chi2_diff","Events");
	Reco_Chi2_orig_diff=HConfig.GetTH1D(Name+"_Reco_Chi2_orig_diff","Reco_Chi2_orig_diff",100,-2,2,"Reco_Chi2_orig_diff","Events");
	Reco_Chi2_HC_diff=HConfig.GetTH1D(Name+"_Reco_Chi2_HC_diff","Reco_Chi2_HC_diff",100,-1,1,"Reco_Chi2_HC_diff","Events");
	Reco_Chi2_SC_diff=HConfig.GetTH1D(Name+"_Reco_Chi2_SC_diff","Reco_Chi2_SC_diff",100,-2,2,"Reco_Chi2_SC_diff","Events");
	Reco_Chi2_origplusSC_diff=HConfig.GetTH1D(Name+"_Reco_Chi2_origplusSC_diff","Reco_Chi2_origplusSC_diff",100,-2,2,"Reco_Chi2_origplusSC_diff","Events");

	Reco_Chi2_diff_vs_correctAssignment=HConfig.GetTH2D(Name+"_Reco_Chi2_diff_vs_correctAssignment","Reco_Chi2_diff_vs_correctAssignment",40,-2,2,2,-0.5,1.5,"Reco_Chi2_diff","correctAssignment");

	Reco_FitSolution_byChi2_Full_vs_RightSolution=HConfig.GetTH2D(Name+"_Reco_FitSolution_byChi2_Full_vs_RightSolution","Reco_FitSolution_byChi2_Full_vs_RightSolution",3,-0.5,2.5,3,-0.5,2.5,"RightSolution","Reco_FitSolution_byChi2_Full_vs_RightSolution");
	Reco_FitSolution_byChi2_orig_vs_RightSolution=HConfig.GetTH2D(Name+"_Reco_FitSolution_byChi2_orig_vs_RightSolution","Reco_FitSolution_byChi2_orig_vs_RightSolution",3,-0.5,2.5,3,-0.5,2.5,"RightSolution","Reco_FitSolution_byChi2_orig_vs_RightSolution");
	Reco_FitSolution_byChi2_SC_vs_RightSolution=HConfig.GetTH2D(Name+"_Reco_FitSolution_byChi2_SC_vs_RightSolution","Reco_FitSolution_byChi2_SC_vs_RightSolution",3,-0.5,2.5,3,-0.5,2.5,"RightSolution","Reco_FitSolution_byChi2_SC_vs_RightSolution");
	Reco_FitSolution_byChi2_HC_vs_RightSolution=HConfig.GetTH2D(Name+"_Reco_FitSolution_byChi2_HC_vs_RightSolution","Reco_FitSolution_byChi2_HC_vs_RightSolution",3,-0.5,2.5,3,-0.5,2.5,"RightSolution","Reco_FitSolution_byChi2_HC_vs_RightSolution");
	Reco_FitSolution_byChi2_origplusSC_vs_RightSolution=HConfig.GetTH2D(Name+"_Reco_FitSolution_byChi2_origplusSC_vs_RightSolution","Reco_FitSolution_byChi2_origplusSC_vs_RightSolution",3,-0.5,2.5,3,-0.5,2.5,"RightSolution","Reco_FitSolution_byChi2_origplusSC");

	//DiTau Reco with Pt recoil estimate

	Reco_PtRes_TauA1_wRecoil=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_wRecoil","Reco_PtRes_TauA1_wRecoil",100,-50,50,"Reco_PtRes_TauA1_wRecoil","Events");
	Reco_PtRes_TauA1_wRecoil_AmbZero=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_wRecoil_AmbZero","Reco_PtRes_TauA1_wRecoil_AmbZero",100,-50,50,"Reco_PtRes_TauA1_wRecoil_AmbZero","Events");
	Reco_PtRes_TauA1_wRecoil_wAmb=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_wRecoil_wAmb","Reco_PtRes_TauA1_wRecoil_wAmb",100,-50,50,"Reco_PtRes_TauA1_wRecoil_wAmb","Events");
	Reco_PtRes_TauMu_wRecoil=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_wRecoil","Reco_PtRes_TauMu_wRecoil",100,-50,50,"Reco_PtRes_TauMu_wRecoil","Events");
	Reco_PtRes_TauMu_wRecoil_AmbZero=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_wRecoil_AmbZero","Reco_PtRes_TauMu_wRecoil_AmbZero",100,-50,50,"Reco_PtRes_TauMu_wRecoil_AmbZero","Events");
	Reco_PtRes_TauMu_wRecoil_wAmb=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_wRecoil_wAmb","Reco_PtRes_TauMu_wRecoil_wAmb",100,-50,50,"Reco_PtRes_TauMu_wRecoil_wAmb","Events");

	Reco_PtRes_TauA1_wRecoil_PreFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauA1_wRecoil_PreFit","Reco_PtRes_TauA1_wRecoil_PreFit",100,-50,50,"Reco_PtRes_TauA1_wRecoil_PreFit","Events");
	Reco_PtRes_TauMu_wRecoil_PreFit=HConfig.GetTH1D(Name+"_Reco_PtRes_TauMu_wRecoil_PreFit","Reco_PtRes_TauMu_wRecoil_PreFit",100,-50,50,"Reco_PtRes_TauMu_wRecoil_PreFit","Events");

	Reco_dPhi_TauMuTauA1_AfterFit_wRecoil=HConfig.GetTH1D(Name+"_Reco_dPhi_TauMuTauA1_AfterFit_wRecoil","Reco_dPhi_TauMuTauA1_AfterFit_wRecoil",64,-3.14159265359,3.14159265359,"Reco_dPhi_TauMuTauA1_AfterFit_wRecoil","Events");
	Reco_dPhi_TauMuTauA1_BeforeFit_wRecoil=HConfig.GetTH1D(Name+"_Reco_dPhi_TauMuTauA1_BeforeFit_wRecoil","Reco_dPhi_TauMuTauA1_BeforeFit_wRecoil",64,-3.14159265359,3.14159265359,"Reco_dPhi_TauMuTauA1_BeforeFit_wRecoil","Events");

	Reco_EventFit_Solution_wRecoil=HConfig.GetTH1D(Name+"_Reco_EventFit_Solution_wRecoil","Reco_EventFit_Solution_wRecoil",5,-2.5,2.5,"Reco_EventFit_Solution_wRecoil","Events");
	Reco_Chi2_FitSolutionOnly_wRecoil=HConfig.GetTH1D(Name+"_Reco_Chi2_FitSolutionOnly_wRecoil","Reco_Chi2_FitSolutionOnly_wRecoil",30,0,30,"Reco_Chi2_FitSolutionOnly_wRecoil","Events");
	Reco_Chi2_Full_wRecoil=HConfig.GetTH1D(Name+"_Reco_Chi2_Full_wRecoil","Reco_Chi2_Full_wRecoil",50,0,1,"Reco_Chi2_Full_wRecoil","Events");
	Reco_Chi2_Orig_wRecoil=HConfig.GetTH1D(Name+"_Reco_Chi2_Orig_wRecoil","Reco_Chi2_Orig_wRecoil",50,0,1,"Reco_Chi2_Orig_wRecoil","Events");
	Reco_Chi2_SC_wRecoil=HConfig.GetTH1D(Name+"_Reco_Chi2_SC_wRecoil","Reco_Chi2_SC_wRecoil",50,0,1,"Reco_Chi2_SC_wRecoil","Events");
	Reco_Chi2_HC_wRecoil=HConfig.GetTH1D(Name+"_Reco_Chi2_HC_wRecoil","Reco_Chi2_HC_wRecoil",50,0,1,"Reco_Chi2_HC_wRecoil","Events");
	Reco_Chi2_OrigProb_wRecoil=HConfig.GetTH1D(Name+"_Reco_Chi2_OrigProb_wRecoil","Reco_Chi2_OrigProb_wRecoil",50,0,1,"Reco_Chi2_OrigProb_wRecoil","Events");
	Reco_Chi2_diff_wRecoil=HConfig.GetTH1D(Name+"_Reco_Chi2_diff_wRecoil","Reco_Chi2_diff_wRecoil",100,-2,2,"Reco_Chi2_diff_wRecoil","Events");
	Reco_Chi2_orig_diff_wRecoil=HConfig.GetTH1D(Name+"_Reco_Chi2_orig_diff_wRecoil","Reco_Chi2_orig_diff_wRecoil",100,-2,2,"Reco_Chi2_orig_diff_wRecoil","Events");
	Reco_ZMass_wRecoil=HConfig.GetTH1D(Name+"_Reco_ZMass_wRecoil","Reco_ZMass_wRecoil",180,60,150,"Reco_ZMass_wRecoil","Events");
	Reco_NIter_wRecoil=HConfig.GetTH1D(Name+"_Reco_NIter_wRecoil","Reco_NIter_wRecoil",100,0,100,"Reco_NIter_wRecoil","Events");
	Reco_Z_Energy_Res_wRecoil=HConfig.GetTH1D(Name+"_Reco_Z_Energy_Res_wRecoil","Reco_Z_Energy_Res_wRecoil",100,-100,100,"Reco_Z_Energy_Res_wRecoil","Events");

	Reco_TauA1_ptRes_vs_ptGen_wRecoil=HConfig.GetTH2D(Name+"_Reco_TauA1_ptRes_vs_ptGen_wRecoil","Reco_TauA1_ptRes_vs_ptGen_wRecoil",100,-50,50,50,0.0,50,"Reco_TauA1_ptRes_wRecoil","ptGen");
	Reco_TauA1_ptRes_vs_ptReco_wRecoil=HConfig.GetTH2D(Name+"_Reco_TauA1_ptRes_vs_ptReco_wRecoil","Reco_TauA1_ptRes_vs_ptReco_wRecoil",100,-50,50,50,0.0,100,"Reco_TauA1_ptRes_wRecoil","ptReco");

	Reco_TauMu_ptRes_vs_ptGen_wRecoil=HConfig.GetTH2D(Name+"_Reco_TauMu_ptRes_vs_ptGen_wRecoil","Reco_TauMu_ptRes_vs_ptGen_wRecoil",100,-50,50,50,0.0,50,"Reco_TauMu_ptRes_wRecoil","ptGen");
	Reco_TauMu_ptRes_vs_ptReco_wRecoil=HConfig.GetTH2D(Name+"_Reco_TauMu_ptRes_vs_ptReco_wRecoil","Reco_TauMu_ptRes_vs_ptReco_wRecoil",100,-50,50,50,0.0,100,"Reco_TauMu_ptRes_wRecoil","ptReco");

	Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr=HConfig.GetTH2D(Name+"_Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr","Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr",100,-50,50,50,0.0,100,"Reco_TauA1_ptRes_wRecoilCorr","ptReco");
	Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr=HConfig.GetTH2D(Name+"_Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr","Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr",100,-50,50,50,0.0,100,"Reco_TauMu_ptRes_wRecoilCorr","ptReco");

	Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_wAmb=HConfig.GetTH2D(Name+"_Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_wAmb","Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_wAmb",100,-50,50,50,0.0,100,"Reco_TauA1_ptRes_wRecoilCorr_wAmb","ptReco");
	Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_AmbZero=HConfig.GetTH2D(Name+"_Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_AmbZero","Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_AmbZero",100,-50,50,50,0.0,100,"Reco_TauA1_ptRes_wRecoilCorr_AmbZero","ptReco");

	Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_wAmb=HConfig.GetTH2D(Name+"_Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_wAmb","Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_wAmb",100,-50,50,50,0.0,100,"Reco_TauMu_ptRes_wRecoilCorr_wAmb","ptReco");
	Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_AmbZero=HConfig.GetTH2D(Name+"_Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_AmbZero","Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_AmbZero",100,-50,50,50,0.0,100,"Reco_TauMu_ptRes_wRecoilCorr_AmbZero","ptReco");

	Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit=HConfig.GetTH2D(Name+"_Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit","Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit",100,-50,50,50,0.0,100,"Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit","ptReco");
	Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit=HConfig.GetTH2D(Name+"_Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit","Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit",100,-50,50,50,0.0,100,"Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit","ptReco");
	Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit=HConfig.GetTH2D(Name+"_Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit","Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit",100,-50,50,50,0.0,100,"Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit","ptReco");
	Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit=HConfig.GetTH2D(Name+"_Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit","Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit",100,-50,50,50,0.0,100,"Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit","ptReco");
	Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit=HConfig.GetTH2D(Name+"_Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit","Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit",100,-50,50,50,0.0,100,"Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit","ptReco");
	Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit=HConfig.GetTH2D(Name+"_Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit","Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit",100,-50,50,50,0.0,100,"Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit","ptReco");

	Reco_TauA1_P_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauA1_P_wRecoil","Reco_TauA1_P_wRecoil",75, 0,150,"Reco_TauA1_P_wRecoil","Events");
	Reco_TauA1_Pt_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauA1_Pt_wRecoil","Reco_TauA1_Pt_wRecoil",75, 0,150,"Reco_TauA1_Pt_wRecoil","Events");
	Reco_TauA1_Px_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauA1_Px_wRecoil","Reco_TauA1_Px_wRecoil",100,-100,100,"Reco_TauA1_Px_wRecoil","Events");
	Reco_TauA1_Py_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauA1_Py_wRecoil","Reco_TauA1_Py_wRecoil",100,-100,100,"Reco_TauA1_Py_wRecoil","Events");
	Reco_TauA1_Pz_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauA1_Pz_wRecoil","Reco_TauA1_Pz_wRecoil",100,-100,100,"Reco_TauA1_Pz_wRecoil","Events");
	Reco_TauA1_Phi_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauA1_Phi_wRecoil","Reco_TauA1_Phi_wRecoil",64,-3.14159265359,3.14159265359,"Reco_TauA1_Phi_wRecoil","Events");
	Reco_TauA1_Eta_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauA1_Eta_wRecoil","Reco_TauA1_Eta_wRecoil",50,-2.5,2.5,"Reco_TauA1_Eta_wRecoil","Events");

	Reco_TauMu_P_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauMu_P_wRecoil","Reco_TauMu_P_wRecoil",75, 0,150,"Reco_TauMu_P_wRecoil","Events");
	Reco_TauMu_Pt_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauMu_Pt_wRecoil","Reco_TauMu_Pt_wRecoil",75, 0,150,"Reco_TauMu_Pt_wRecoil","Events");
	Reco_TauMu_Px_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauMu_Px_wRecoil","Reco_TauMu_Px_wRecoil",100,-100,100,"Reco_TauMu_Px_wRecoil","Events");
	Reco_TauMu_Py_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauMu_Py_wRecoil","Reco_TauMu_Py_wRecoil",100,-100,100,"Reco_TauMu_Py_wRecoil","Events");
	Reco_TauMu_Pz_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauMu_Pz_wRecoil","Reco_TauMu_Pz_wRecoil",100,-100,100,"Reco_TauMu_Pz_wRecoil","Events");
	Reco_TauMu_Phi_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauMu_Phi_wRecoil","Reco_TauMu_Phi_wRecoil",64,-3.14159265359,3.14159265359,"Reco_TauMu_Phi_wRecoil","Events");
	Reco_TauMu_Eta_wRecoil=HConfig.GetTH1D(Name+"_Reco_TauMu_Eta_wRecoil","Reco_TauMu_Eta_wRecoil",50,-2.5,2.5,"Reco_TauMu_Eta_wRecoil","Events");

	Reco_Z_Px_wRecoil=HConfig.GetTH1D(Name+"_Reco_Z_Px_wRecoil","Reco_Z_Px_wRecoil",100,-100,100,"Reco_Z_Px_wRecoil","Events");
	Reco_Z_Py_wRecoil=HConfig.GetTH1D(Name+"_Reco_Z_Py_wRecoil","Reco_Z_Py_wRecoil",100,-100,100,"Reco_Z_Py_wRecoil","Events");
	Reco_Z_Pz_wRecoil=HConfig.GetTH1D(Name+"_Reco_Z_Pz_wRecoil","Reco_Z_Pz_wRecoil",100,-100,100,"Reco_Z_Pz_wRecoil","Events");
	Reco_Z_Pt_wRecoil=HConfig.GetTH1D(Name+"_Reco_Z_Pt_wRecoil","Reco_Z_Pt_wRecoil",100,0,100,"Reco_Z_Pt_wRecoil","Events");
	Reco_Z_PtRes_wRecoil=HConfig.GetTH1D(Name+"_Reco_Z_PtRes_wRecoil","Reco_Z_PtRes_wRecoil",100,-100,100,"Reco_Z_PtRes_wRecoil","Events");
	Reco_Z_Phi_wRecoil=HConfig.GetTH1D(Name+"_Reco_Z_Phi_wRecoil","Reco_Z_Phi_wRecoil",64,-3.14159265359,3.14159265359,"Reco_Z_Phi_wRecoil","Events");
	Reco_Z_PhiRes_wRecoil=HConfig.GetTH1D(Name+"_Reco_Z_PhiRes_wRecoil","Reco_Z_PhiRes_wRecoil",64,-3.14159265359,3.14159265359,"Reco_Z_PhiRes_wRecoil","Events");
	Reco_Z_Eta_wRecoil=HConfig.GetTH1D(Name+"_Reco_Z_Eta_wRecoil","Reco_Z_Eta_wRecoil",100,-10,10,"Reco_Z_Eta_wRecoil","Events");
	Reco_Z_EtaRes_wRecoil=HConfig.GetTH1D(Name+"_Reco_Z_EtaRes_wRecoil","Reco_Z_EtaRes_wRecoil",100,-10,10,"Reco_Z_EtaRes_wRecoil","Events");

	Reco_FitSolution_byChi2_Full_vs_RightSolution_wRecoil=HConfig.GetTH2D(Name+"_Reco_FitSolution_byChi2_Full_vs_RightSolution_wRecoil","Reco_FitSolution_byChi2_Full_vs_RightSolution_wRecoil",3,-0.5,2.5,3,-0.5,2.5,"RightSolution","Reco_FitSolution_byChi2_Full_vs_RightSolution_wRecoil");
	Reco_FitSolution_byChi2_orig_vs_RightSolution_wRecoil=HConfig.GetTH2D(Name+"_Reco_FitSolution_byChi2_orig_vs_RightSolution_wRecoil","Reco_FitSolution_byChi2_orig_vs_RightSolution_wRecoil",3,-0.5,2.5,3,-0.5,2.5,"RightSolution","Reco_FitSolution_byChi2_orig_vs_RightSolution_wRecoil");
	Reco_FitSolution_byChi2_SC_vs_RightSolution_wRecoil=HConfig.GetTH2D(Name+"_Reco_FitSolution_byChi2_SC_vs_RightSolution_wRecoil","Reco_FitSolution_byChi2_SC_vs_RightSolution_wRecoil",3,-0.5,2.5,3,-0.5,2.5,"RightSolution","Reco_FitSolution_byChi2_SC_vs_RightSolution_wRecoil");
	Reco_FitSolution_byChi2_HC_vs_RightSolution_wRecoil=HConfig.GetTH2D(Name+"_Reco_FitSolution_byChi2_HC_vs_RightSolution_wRecoil","Reco_FitSolution_byChi2_HC_vs_RightSolution_wRecoil",3,-0.5,2.5,3,-0.5,2.5,"RightSolution","Reco_FitSolution_byChi2_HC_vs_RightSolution_wRecoil");
	Reco_FitSolution_byChi2_origplusSC_vs_RightSolution_wRecoil=HConfig.GetTH2D(Name+"_Reco_FitSolution_byChi2_origplusSC_vs_RightSolution_wRecoil","Reco_FitSolution_byChi2_origplusSC_vs_RightSolution_wRecoil",3,-0.5,2.5,3,-0.5,2.5,"RightSolution","Reco_FitSolution_byChi2_origplusSC_wRecoil");

	//DiTau Reco with generator information/particles
	GenReco_ZMass=HConfig.GetTH1D(Name+"_GenReco_ZMass","GenReco_ZMass",180,60,150,"GenReco_ZMass","Events");
	GenReco_EventFit_Solution=HConfig.GetTH1D(Name+"_GenReco_EventFit_Solution","GenReco_EventFit_Solution",5,-2.5,2.5,"GenReco_EventFit_Solution","Events");
	GenReco_A1Fit_Solution=HConfig.GetTH1D(Name+"_GenReco_A1Fit_Solution","GenReco_A1Fit_Solution",5,-2.5,2.5,"GenReco_A1Fit_Solution","Events");
	GenReco_Chi2=HConfig.GetTH1D(Name+"_GenReco_Chi2","GenReco_Chi2",30,0,30,"GenReco_Chi2","Events");
	GenReco_Chi2_FitSolutionOnly=HConfig.GetTH1D(Name+"_GenReco_Chi2_FitSolutionOnly","GenReco_Chi2_FitSolutionOnly",30,0,30,"GenReco_Chi2_FitSolutionOnly","Events");
	GenReco_ConstrainedDeltaSum=HConfig.GetTH1D(Name+"_GenReco_ConstrainedDeltaSum","GenReco_ConstrainedDeltaSum",50,0,0.1,"GenReco_ConstrainedDeltaSum","Events");
	GenReco_NIter=HConfig.GetTH1D(Name+"_GenReco_NIter","GenReco_NIter",100,0,100,"GenReco_NIter","Events");

	//TauMu starting point
	TauMu_Start_Collinear_PtRes=HConfig.GetTH1D(Name+"_TauMu_Start_Collinear_PtRes","TauMu_Start_Collinear_PtRes",100,-40,40,"TauMu_Start_Collinear_PtRes","Events");
	TauMu_Start_Collinear_PtReco_vs_PtGen=HConfig.GetTH2D(Name+"_TauMu_Start_Collinear_PtReco_vs_PtGen","TauMu_Start_Collinear_PtReco_vs_PtGen",50,10.,50,50,0.0,50,"TauMu_Start_Collinear_PtReco","PtGen");

	TauMu_Start_MET_PtRes=HConfig.GetTH1D(Name+"_TauMu_Start_MET_PtRes","TauMu_Start_MET_PtRes",100,-40,40,"TauMu_Start_MET_PtRes","Events");
	TauMu_Start_MET_PtReco_vs_PtGen=HConfig.GetTH2D(Name+"_TauMu_Start_MET_PtReco_vs_PtGen","TauMu_Start_MET_PtReco_vs_PtGen",50,0.,50,50,0.0,50,"TauMu_Start_MET_PtReco","PtGen");

	TauMu_Start_PtBalance_PtRes=HConfig.GetTH1D(Name+"_TauMu_Start_PtBalance_PtRes","TauMu_Start_PtBalance_PtRes",100,-40,40,"TauMu_Start_PtBalance_PtRes","Events");
	TauMu_Start_PtBalance_PtReco_vs_PtGen=HConfig.GetTH2D(Name+"_TauMu_Start_PtBalance_PtReco_vs_PtGen","TauMu_Start_PtBalance_PtReco_vs_PtGen",50,0.,50,50,0.0,50,"TauMu_Start_PtBalance_PtReco","PtGen");

	TauMu_Start_EventRecoil_PtRes=HConfig.GetTH1D(Name+"_TauMu_Start_EventRecoil_PtRes","TauMu_Start_EventRecoil_PtRes",100,-40,40,"TauMu_Start_EventRecoil_PtRes","Events");
	TauMu_Start_EventRecoil_PtReco_vs_PtGen=HConfig.GetTH2D(Name+"_TauMu_Start_EventRecoil_PtReco_vs_PtGen","TauMu_Start_EventRecoil_PtReco_vs_PtGen",50,0.,50,50,0.0,50,"TauMu_Start_EventRecoil_PtReco","PtGen");

	TauMu_Start_MET_PtRes_AfterMC=HConfig.GetTH1D(Name+"_TauMu_Start_MET_PtRes_AfterMC","TauMu_Start_MET_PtRes_AfterMC",100,-40,40,"TauMu_Start_MET_PtRes_AfterMC","Events");
	TauMu_Start_PtBalance_PtRes_AfterMC=HConfig.GetTH1D(Name+"_TauMu_Start_PtBalance_PtRes_AfterMC","TauMu_Start_PtBalance_PtRes_AfterMC",100,-40,40,"TauMu_Start_PtBalance_PtRes_AfterMC","Events");
	TauMu_Start_EventRecoil_PtRes_AfterMC=HConfig.GetTH1D(Name+"_TauMu_Start_EventRecoil_PtRes_AfterMC","TauMu_Start_EventRecoil_PtRes_AfterMC",100,-40,40,"TauMu_Start_EventRecoil_PtRes_AfterMC","Events");

	TauMu_Start_dPhi_TauMuTauH_MET=HConfig.GetTH1D(Name+"_TauMu_Start_dPhi_TauMuTauH_MET","TauMu_Start_dPhi_TauMuTauH_MET",64,-3.14159265359,3.14159265359,"TauMu_Start_dPhi_TauMuTauH_MET","Events");
	TauMu_Start_dPhi_TauMuTauH_EventRecoil=HConfig.GetTH1D(Name+"_TauMu_Start_dPhi_TauMuTauH_EventRecoil","TauMu_Start_dPhi_TauMuTauH_EventRecoil",64,-3.14159265359,3.14159265359,"TauMu_Start_dPhi_TauMuTauH_EventRecoil","Events");
	TauMu_Start_dPhi_TauMuTauH_PtBalance=HConfig.GetTH1D(Name+"_TauMu_Start_dPhi_TauMuTauH_PtBalance","TauMu_Start_dPhi_TauMuTauH_PtBalance",64,-3.14159265359,3.14159265359,"TauMu_Start_dPhi_TauMuTauH_PtBalance","Events");

	Z_Start_MET_PtRes=HConfig.GetTH1D(Name+"_Z_Start_MET_PtRes","Z_Start_MET_PtRes",100,-40,40,"Z_Start_MET_PtRes","Events");
	Z_Start_MET_PhiRes=HConfig.GetTH1D(Name+"_Z_Start_MET_PhiRes","Z_Start_MET_PhiRes",64,-3.14159265359,3.14159265359,"Z_Start_MET_PhiRes","Events");

	Z_Start_PtBalance_PtRes=HConfig.GetTH1D(Name+"_Z_Start_PtBalance_PtRes","Z_Start_PtBalance_PtRes",100,-40,40,"Z_Start_PtBalance_PtRes","Events");
	Z_Start_PtBalance_PhiRes=HConfig.GetTH1D(Name+"_Z_Start_PtBalance_PhiRes","Z_Start_PtBalance_PhiRes",64,-3.14159265359,3.14159265359,"Z_Start_PtBalance_PhiRes","Events");

	Z_Start_EventRecoil_PtRes=HConfig.GetTH1D(Name+"_Z_Start_EventRecoil_PtRes","Z_Start_EventRecoil_PtRes",100,-40,40,"Z_Start_EventRecoil_PtRes","Events");
	Z_Start_EventRecoil_PhiRes=HConfig.GetTH1D(Name+"_Z_Start_EventRecoil_PhiRes","Z_Start_EventRecoil_PhiRes",64,-3.14159265359,3.14159265359,"Z_Start_EventRecoil_PhiRes","Events");

	Z_Start_MET_noNu_PtRes=HConfig.GetTH1D(Name+"_Z_Start_MET_noNu_PtRes","Z_Start_MET_noNu_PtRes",100,-40,40,"Z_Start_MET_noNu_PtRes","Events");
	Z_Start_MET_noNu_PhiRes=HConfig.GetTH1D(Name+"_Z_Start_MET_noNu_PhiRes","Z_Start_MET_noNu_PhiRes",64,-3.14159265359,3.14159265359,"Z_Start_MET_noNu_PhiRes","Events");
	Z_Start_MET_noNu_M=HConfig.GetTH1D(Name+"_Z_Start_MET_noNu_M","Z_Start_MET_noNu_M",100,0,250,"Z_Start_MET_noNu_M","Events");
	Z_Start_MET_noNu_dE=HConfig.GetTH1D(Name+"_Z_Start_MET_noNu_dE","Z_Start_MET_noNu_dE",100,-100,100,"Z_Start_MET_noNu_dE","Events");

	Z_Start_MET_noNu_PtRes_withRot=HConfig.GetTH1D(Name+"_Z_Start_MET_noNu_PtRes_withRot","Z_Start_MET_noNu_PtRes_withRot",100,-40,40,"Z_Start_MET_noNu_PtRes_withRot","Events");
	Z_Start_MET_noNu_PhiRes_withRot=HConfig.GetTH1D(Name+"_Z_Start_MET_noNu_PhiRes_withRot","Z_Start_MET_noNu_PhiRes_withRot",64,-3.14159265359,3.14159265359,"Z_Start_MET_noNu_PhiRes_withRot","Events");
	Z_Start_MET_noNu_M_withRot=HConfig.GetTH1D(Name+"_Z_Start_MET_noNu_M_withRot","Z_Start_MET_noNu_M_withRot",100,0,250,"Z_Start_MET_noNu_M_withRot","Events");

	Z_Start_B2B_M=HConfig.GetTH1D(Name+"_Z_Start_B2B_M","Z_Start_B2B_M",100,0,250,"Z_Start_B2B_M","Events");
	Z_Start_B2B_dE=HConfig.GetTH1D(Name+"_Z_Start_B2B_dE","Z_Start_B2B_dE",100,-100,100,"Z_Start_B2B_dE","Events");

	//Muon Track Parameters

	Mu_TP_phi0=HConfig.GetTH1D(Name+"_Mu_TP_phi0","Mu_TP_phi0",100,-7,7,"Mu_TP_phi0","Events");
	Mu_TP_lambda=HConfig.GetTH1D(Name+"_Mu_TP_lambda","Mu_TP_lambda",100,-7,7,"Mu_TP_lambda","Events");
	Mu_TP_dxy=HConfig.GetTH1D(Name+"_Mu_TP_dxy","Mu_TP_dxy",100,-1,1,"Mu_TP_dxy","Events");
	Mu_TP_dz=HConfig.GetTH1D(Name+"_Mu_TP_dz","Mu_TP_dz",100,-20,20,"Mu_TP_dz","Events");
	Mu_TP_kappa=HConfig.GetTH1D(Name+"_Mu_TP_kappa","Mu_TP_kappa",100,-0.001,0.001,"Mu_TP_kappa","Events");
	Mu_TP_POCA_quadrant=HConfig.GetTH1D(Name+"_Mu_TP_POCA_quadrant","Mu_TP_POCA_quadrant",4,0.5,4.5,"Mu_TP_POCA_quadrant","Events");
	Mu_TP_POCA_quadrantVlad=HConfig.GetTH1D(Name+"_Mu_TP_POCA_quadrantVlad","Mu_TP_POCA_quadrantVlad",4,0.5,4.5,"Mu_TP_POCA_quadrantVlad","Events");
	Mu_TP_POCA_quadrantby_dxyphi0=HConfig.GetTH1D(Name+"_Mu_TP_POCA_quadrantby_dxyphi0","Mu_TP_POCA_quadrantby_dxyphi0",4,0.5,4.5,"Mu_TP_POCA_quadrantby_dxyphi0","Events");

	Mu_TP_Poca_quadrantData=HConfig.GetTH1D(Name+"_Mu_TP_Poca_quadrantData","Mu_TP_Poca_quadrantData",4,0.5,4.5,"Mu_TP_Poca_quadrantData","Events");
	Mu_TP_Poca_quadrantMCDY=HConfig.GetTH1D(Name+"_Mu_TP_Poca_quadrantMCDY","Mu_TP_Poca_quadrantMCDY",4,0.5,4.5,"Mu_TP_Poca_quadrantMCDY","Events");

	Mu_TP_Poca_xy=HConfig.GetTH2D(Name+"_Mu_TP_Poca_xy","Mu_TP_Poca_xy",100,-0.5,0.5,100,-0.5,0.5,"Mu_TP_Poca_x","Mu_TP_Poca_y");
	Mu_TP_Vertex_xy=HConfig.GetTH2D(Name+"_Mu_TP_Vertex_xy","Vertex_xy",100,-0.5,0.5,100,-0.5,0.5,"Mu_TP_Vertex_x","Mu_TP_Vertex_y");
	Mu_TP_RefitVertex_xy=HConfig.GetTH2D(Name+"_Mu_TP_RefitVertex_xy","Mu_TP_RefitVertex_xy",100,-0.5,0.5,100,-0.5,0.5,"Mu_TP_RefitVertex_x","Mu_TP_RefitVertex_y");
	Mu_TP_BeamSpot_xy=HConfig.GetTH2D(Name+"_Mu_TP_BeamSpot_xy","Vertex_xy",100,-0.5,0.5,100,-0.5,0.5,"Mu_TP_BeamSpot_x","Mu_TP_BeamSpot_y");
	Mu_TP_NTP_Poca_xy=HConfig.GetTH2D(Name+"_Mu_TP_NTP_Poca_xy","Mu_TP_NTP_Poca_xy",100,-0.5,0.5,100,-0.5,0.5,"Mu_TP_NTP_Poca_x","Mu_TP_NTP_Poca_y");

	//MET
	MVAMET_metobject_XX=HConfig.GetTH1D(Name+"_MVAMET_metobject_XX","MVAMET_metobject_XX",100,-10,90,"MVAMET_metobject_XX","Events");
	MVAMET_ptobject_XX=HConfig.GetTH1D(Name+"_MVAMET_ptobject_XX","MVAMET_ptobject_XX",100,-10,90,"MVAMET_ptobject_XX","Events");

	//A1
	A1_Par_Px=HConfig.GetTH1D(Name+"_A1_Par_Px","A1_Par_Px",100,-100,100,"A1_Par_Px","Events");
	A1_Par_Py=HConfig.GetTH1D(Name+"_A1_Par_Py","A1_Par_Py",100,-100,100,"A1_Par_Py","Events");
	A1_Par_Pz=HConfig.GetTH1D(Name+"_A1_Par_Pz","A1_Par_Pz",100,-100,100,"A1_Par_Pz","Events");
	A1_Par_M=HConfig.GetTH1D(Name+"_A1_Par_M","A1_Par_M",100,-5,5,"A1_Par_M","Events");
	A1_Cov_Pxx=HConfig.GetTH1D(Name+"_A1_Cov_Pxx","A1_Cov_Pxx",100,0,0.05,"A1_Cov_Pxx","Events");
	A1_Cov_Pyy=HConfig.GetTH1D(Name+"_A1_Cov_Pyy","A1_Cov_Pyy",100,0,0.05,"A1_Cov_Pyy","Events");
	A1_Cov_Pzz=HConfig.GetTH1D(Name+"_A1_Cov_Pzz","A1_Cov_Pzz",100,0,0.05,"A1_Cov_Pzz","Events");
	A1_Cov_Pxy=HConfig.GetTH1D(Name+"_A1_Cov_Pxy","A1_Cov_Pxy",100,-0.05,0.05,"A1_Cov_Pxy","Events");
	A1_Cov_Pxz=HConfig.GetTH1D(Name+"_A1_Cov_Pxz","A1_Cov_Pxz",100,-0.05,0.05,"A1_Cov_Pxz","Events");
	A1_Cov_Pyz=HConfig.GetTH1D(Name+"_A1_Cov_Pyz","A1_Cov_Pyz",100,-0.05,0.05,"A1_Cov_Pyz","Events");
	A1_Pull_Px=HConfig.GetTH1D(Name+"_A1_Pull_Px","A1_Pull_Px",100,-5,5,"A1_Pull_Px","Events");
	A1_Pull_Py=HConfig.GetTH1D(Name+"_A1_Pull_Py","A1_Pull_Py",100,-5,5,"A1_Pull_Py","Events");
	A1_Pull_Pz=HConfig.GetTH1D(Name+"_A1_Pull_Pz","A1_Pull_Pz",100,-5,5,"A1_Pull_Pz","Events");
	A1_Pull_M=HConfig.GetTH1D(Name+"_A1_Pull_M","A1_Pull_M",100,-5,5,"A1_Pull_M","Events");
	A1_Pull_Px_AmbZero=HConfig.GetTH1D(Name+"_A1_Pull_Px_AmbZero","A1_Pull_Px_AmbZero",100,-5,5,"A1_Pull_Px_AmbZero","Events");
	A1_Pull_Py_AmbZero=HConfig.GetTH1D(Name+"_A1_Pull_Py_AmbZero","A1_Pull_Py_AmbZero",100,-5,5,"A1_Pull_Py_AmbZero","Events");
	A1_Pull_Pz_AmbZero=HConfig.GetTH1D(Name+"_A1_Pull_Pz_AmbZero","A1_Pull_Pz_AmbZero",100,-5,5,"A1_Pull_Pz_AmbZero","Events");
	A1_Pull_M_AmbZero=HConfig.GetTH1D(Name+"_A1_Pull_M_AmbZero","A1_Pull_M_AmbZero",100,-5,5,"A1_Pull_M_AmbZero","Events");
	A1_Pull_Px_wAmb=HConfig.GetTH1D(Name+"_A1_Pull_Px_wAmb","A1_Pull_Px_wAmb",100,-5,5,"A1_Pull_Px_wAmb","Events");
	A1_Pull_Py_wAmb=HConfig.GetTH1D(Name+"_A1_Pull_Py_wAmb","A1_Pull_Py_wAmb",100,-5,5,"A1_Pull_Py_wAmb","Events");
	A1_Pull_Pz_wAmb=HConfig.GetTH1D(Name+"_A1_Pull_Pz_wAmb","A1_Pull_Pz_wAmb",100,-5,5,"A1_Pull_Pz_wAmb","Events");
	A1_Pull_M_wAmb=HConfig.GetTH1D(Name+"_A1_Pull_M_wAmb","A1_Pull_M_wAmb",100,-5,5,"A1_Pull_M_wAmb","Events");

	//SV
	SV_Par_x=HConfig.GetTH1D(Name+"_SV_Par_x","SV_Par_x",100,-2,2,"SV_Par_x","Events");
	SV_Par_y=HConfig.GetTH1D(Name+"_SV_Par_y","SV_Par_y",100,-2,2,"SV_Par_y","Events");
	SV_Par_z=HConfig.GetTH1D(Name+"_SV_Par_z","SV_Par_z",100,-20,20,"SV_Par_z","Events");
	SV_Cov_xx=HConfig.GetTH1D(Name+"_SV_Cov_xx","SV_Cov_xx",100,0,0.01,"SV_Cov_xx","Events");
	SV_Cov_yy=HConfig.GetTH1D(Name+"_SV_Cov_yy","SV_Cov_yy",100,0,0.01,"SV_Cov_yy","Events");
	SV_Cov_zz=HConfig.GetTH1D(Name+"_SV_Cov_zz","SV_Cov_zz",100,0,0.2,"SV_Cov_zz","Events");
	SV_Cov_xy=HConfig.GetTH1D(Name+"_SV_Cov_xy","SV_Cov_xy",100,-0.01,0.01,"SV_Cov_xy","Events");
	SV_Cov_xz=HConfig.GetTH1D(Name+"_SV_Cov_xz","SV_Cov_xz",100,-0.01,0.01,"SV_Cov_xz","Events");
	SV_Cov_yz=HConfig.GetTH1D(Name+"_SV_Cov_yz","SV_Cov_yz",100,-0.01,0.01,"SV_Cov_yz","Events");
	SV_Pull_Px=HConfig.GetTH1D(Name+"_SV_Pull_Px","SV_Pull_Px",100,-5,5,"SV_Pull_Px","Events");
	SV_Pull_Py=HConfig.GetTH1D(Name+"_SV_Pull_Py","SV_Pull_Py",100,-5,5,"SV_Pull_Py","Events");
	SV_Pull_Pz=HConfig.GetTH1D(Name+"_SV_Pull_Pz","SV_Pull_Pz",100,-5,5,"SV_Pull_Pz","Events");
	SV_Pull_Px_AmbZero=HConfig.GetTH1D(Name+"_SV_Pull_Px_AmbZero","SV_Pull_Px_AmbZero",100,-5,5,"SV_Pull_Px_AmbZero","Events");
	SV_Pull_Py_AmbZero=HConfig.GetTH1D(Name+"_SV_Pull_Py_AmbZero","SV_Pull_Py_AmbZero",100,-5,5,"SV_Pull_Py_AmbZero","Events");
	SV_Pull_Pz_AmbZero=HConfig.GetTH1D(Name+"_SV_Pull_Pz_AmbZero","SV_Pull_Pz_AmbZero",100,-5,5,"SV_Pull_Pz_AmbZero","Events");
	SV_Pull_Px_wAmb=HConfig.GetTH1D(Name+"_SV_Pull_Px_wAmb","SV_Pull_Px_wAmb",100,-5,5,"SV_Pull_Px_wAmb","Events");
	SV_Pull_Py_wAmb=HConfig.GetTH1D(Name+"_SV_Pull_Py_wAmb","SV_Pull_Py_wAmb",100,-5,5,"SV_Pull_Py_wAmb","Events");
	SV_Pull_Pz_wAmb=HConfig.GetTH1D(Name+"_SV_Pull_Pz_wAmb","SV_Pull_Pz_wAmb",100,-5,5,"SV_Pull_Pz_wAmb","Events");

	//TauA1
	TauA1_Par_Px_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Par_Px_AmbZero","TauA1_Par_Px_AmbZero",100,-100,100,"TauA1_Par_Px_AmbZero","Events");
	TauA1_Par_Py_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Par_Py_AmbZero","TauA1_Par_Py_AmbZero",100,-100,100,"TauA1_Par_Py_AmbZero","Events");
	TauA1_Par_Pz_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Par_Pz_AmbZero","TauA1_Par_Pz_AmbZero",100,-100,100,"TauA1_Par_Pz_AmbZero","Events");
	TauA1_Cov_Pxx_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Cov_Pxx_AmbZero","TauA1_Cov_Pxx_AmbZero",100,-10,490,"TauA1_Cov_Pxx_AmbZero","Events");
	TauA1_Cov_Pyy_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Cov_Pyy_AmbZero","TauA1_Cov_Pyy_AmbZero",100,-10,490,"TauA1_Cov_Pyy_AmbZero","Events");
	TauA1_Cov_Pzz_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Cov_Pzz_AmbZero","TauA1_Cov_Pzz_AmbZero",100,-10,490,"TauA1_Cov_Pzz_AmbZero","Events");
	TauA1_Cov_Pxy_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Cov_Pxy_AmbZero","TauA1_Cov_Pxy_AmbZero",100,-100,100,"TauA1_Cov_Pxy_AmbZero","Events");
	TauA1_Cov_Pxz_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Cov_Pxz_AmbZero","TauA1_Cov_Pxz_AmbZero",100,-100,100,"TauA1_Cov_Pxz_AmbZero","Events");
	TauA1_Cov_Pyz_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Cov_Pyz_AmbZero","TauA1_Cov_Pyz_AmbZero",100,-100,100,"TauA1_Cov_Pyz_AmbZero","Events");
	TauA1_Pull_Px_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Pull_Px_AmbZero","TauA1_Pull_Px_AmbZero",100,-5,5,"TauA1_Pull_Px_AmbZero","Events");
	TauA1_Pull_Py_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Pull_Py_AmbZero","TauA1_Pull_Py_AmbZero",100,-5,5,"TauA1_Pull_Py_AmbZero","Events");
	TauA1_Pull_Pz_AmbZero=HConfig.GetTH1D(Name+"_TauA1_Pull_Pz_AmbZero","TauA1_Pull_Pz_AmbZero",100,-5,5,"TauA1_Pull_Pz_AmbZero","Events");

	TauA1_Par_Px_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Par_Px_CorrectAmb","TauA1_Par_Px_CorrectAmb",100,-100,100,"TauA1_Par_Px_CorrectAmb","Events");
	TauA1_Par_Py_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Par_Py_CorrectAmb","TauA1_Par_Py_CorrectAmb",100,-100,100,"TauA1_Par_Py_CorrectAmb","Events");
	TauA1_Par_Pz_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Par_Pz_CorrectAmb","TauA1_Par_Pz_CorrectAmb",100,-100,100,"TauA1_Par_Pz_CorrectAmb","Events");
	TauA1_Cov_Pxx_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pxx_CorrectAmb","TauA1_Cov_Pxx_CorrectAmb",100,-10,490,"TauA1_Cov_Pxx_CorrectAmb","Events");
	TauA1_Cov_Pyy_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pyy_CorrectAmb","TauA1_Cov_Pyy_CorrectAmb",100,-10,490,"TauA1_Cov_Pyy_CorrectAmb","Events");
	TauA1_Cov_Pzz_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pzz_CorrectAmb","TauA1_Cov_Pzz_CorrectAmb",100,-10,490,"TauA1_Cov_Pzz_CorrectAmb","Events");
	TauA1_Cov_Pxy_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pxy_CorrectAmb","TauA1_Cov_Pxy_CorrectAmb",100,-100,100,"TauA1_Cov_Pxy_CorrectAmb","Events");
	TauA1_Cov_Pxz_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pxz_CorrectAmb","TauA1_Cov_Pxz_CorrectAmb",100,-100,100,"TauA1_Cov_Pxz_CorrectAmb","Events");
	TauA1_Cov_Pyz_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pyz_CorrectAmb","TauA1_Cov_Pyz_CorrectAmb",100,-100,100,"TauA1_Cov_Pyz_CorrectAmb","Events");
	TauA1_Pull_Px_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Pull_Px_CorrectAmb","TauA1_Pull_Px_CorrectAmb",100,-5,5,"TauA1_Pull_Px_CorrectAmb","Events");
	TauA1_Pull_Py_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Pull_Py_CorrectAmb","TauA1_Pull_Py_CorrectAmb",100,-5,5,"TauA1_Pull_Py_CorrectAmb","Events");
	TauA1_Pull_Pz_CorrectAmb=HConfig.GetTH1D(Name+"_TauA1_Pull_Pz_CorrectAmb","TauA1_Pull_Pz_CorrectAmb",100,-5,5,"TauA1_Pull_Pz_CorrectAmb","Events");

	TauA1_Par_Px_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Par_Px_WrongAmb","TauA1_Par_Px_WrongAmb",100,-100,100,"TauA1_Par_Px_WrongAmb","Events");
	TauA1_Par_Py_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Par_Py_WrongAmb","TauA1_Par_Py_WrongAmb",100,-100,100,"TauA1_Par_Py_WrongAmb","Events");
	TauA1_Par_Pz_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Par_Pz_WrongAmb","TauA1_Par_Pz_WrongAmb",100,-100,100,"TauA1_Par_Pz_WrongAmb","Events");
	TauA1_Cov_Pxx_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pxx_WrongAmb","TauA1_Cov_Pxx_WrongAmb",100,-10,490,"TauA1_Cov_Pxx_WrongAmb","Events");
	TauA1_Cov_Pyy_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pyy_WrongAmb","TauA1_Cov_Pyy_WrongAmb",100,-10,490,"TauA1_Cov_Pyy_WrongAmb","Events");
	TauA1_Cov_Pzz_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pzz_WrongAmb","TauA1_Cov_Pzz_WrongAmb",100,-10,490,"TauA1_Cov_Pzz_WrongAmb","Events");
	TauA1_Cov_Pxy_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pxy_WrongAmb","TauA1_Cov_Pxy_WrongAmb",100,-100,100,"TauA1_Cov_Pxy_WrongAmb","Events");
	TauA1_Cov_Pxz_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pxz_WrongAmb","TauA1_Cov_Pxz_WrongAmb",100,-100,100,"TauA1_Cov_Pxz_WrongAmb","Events");
	TauA1_Cov_Pyz_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Cov_Pyz_WrongAmb","TauA1_Cov_Pyz_WrongAmb",100,-100,100,"TauA1_Cov_Pyz_WrongAmb","Events");
	TauA1_Pull_Px_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Pull_Px_WrongAmb","TauA1_Pull_Px_WrongAmb",100,-5,5,"TauA1_Pull_Px_WrongAmb","Events");
	TauA1_Pull_Py_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Pull_Py_WrongAmb","TauA1_Pull_Py_WrongAmb",100,-5,5,"TauA1_Pull_Py_WrongAmb","Events");
	TauA1_Pull_Pz_WrongAmb=HConfig.GetTH1D(Name+"_TauA1_Pull_Pz_WrongAmb","TauA1_Pull_Pz_WrongAmb",100,-5,5,"TauA1_Pull_Pz_WrongAmb","Events");

	//Tau flight length significance studies

	TauFLSigmaVlad=HConfig.GetTH1D(Name+"_TauFLSigmaVlad","TauFLSigmaVlad",80,-10,30,"TauFLSigmaVlad","Events");
	TauFLSigmaVlad_PhiA1=HConfig.GetTH1D(Name+"_TauFLSigmaVlad_PhiA1","TauFLSigmaVlad_PhiA1",64,-3.14159265359,3.14159265359,"TauFLSigmaVlad_PhiA1","Events");
	TauFLSigmaVlad_PhiTau=HConfig.GetTH1D(Name+"_TauFLSigmaVlad_PhiTau","TauFLSigmaVlad_PhiTau",64,-3.14159265359,3.14159265359,"TauFLSigmaVlad_PhiTau","Events");
	TauFLSigmaVlad_PhiTauNoCorr=HConfig.GetTH1D(Name+"_TauFLSigmaVlad_PhiTauNoCorr","TauFLSigmaVlad_PhiTauNoCorr",64,-3.14159265359,3.14159265359,"TauFLSigmaVlad_PhiTauNoCorr","Events");
	TauFLSigmaVlad_PhiZnoCorr=HConfig.GetTH1D(Name+"_TauFLSigmaVlad_PhiZnoCorr","TauFLSigmaVlad_PhiZnoCorr",64,-3.14159265359,3.14159265359,"TauFLSigmaVlad_PhiZnoCorr","Events");
	TauFLSigmaVlad_PhiZwCorr=HConfig.GetTH1D(Name+"_TauFLSigmaVlad_PhiZwCorr","TauFLSigmaVlad_PhiZwCorr",64,-3.14159265359,3.14159265359,"TauFLSigmaVlad_PhiZwCorr","Events");

	TauFLSigmaAlex=HConfig.GetTH1D(Name+"_TauFLSigmaAlex","TauFLSigmaAlex",80,-10,30,"TauFLSigmaAlex","Events");
	TauFLSigmaAlex_PhiA1=HConfig.GetTH1D(Name+"_TauFLSigmaAlex_PhiA1","TauFLSigmaAlex_PhiA1",64,-3.14159265359,3.14159265359,"TauFLSigmaAlex_PhiA1","Events");
	TauFLSigmaAlex_PhiTau=HConfig.GetTH1D(Name+"_TauFLSigmaAlex_PhiTau","TauFLSigmaAlex_PhiTau",64,-3.14159265359,3.14159265359,"TauFLSigmaAlex_PhiTau","Events");
	TauFLSigmaAlex_PhiTauNoCorr=HConfig.GetTH1D(Name+"_TauFLSigmaAlex_PhiTauNoCorr","TauFLSigmaAlex_PhiTauNoCorr",64,-3.14159265359,3.14159265359,"TauFLSigmaAlex_PhiTauNoCorr","Events");
	TauFLSigmaAlex_PhiZnoCorr=HConfig.GetTH1D(Name+"_TauFLSigmaAlex_PhiZnoCorr","TauFLSigmaAlex_PhiZnoCorr",64,-3.14159265359,3.14159265359,"TauFLSigmaAlex_PhiZnoCorr","Events");
	TauFLSigmaAlex_PhiZwCorr=HConfig.GetTH1D(Name+"_TauFLSigmaAlex_PhiZwCorr","TauFLSigmaAlex_PhiZwCorr",64,-3.14159265359,3.14159265359,"TauFLSigmaAlex_PhiZwCorr","Events");

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
	Extradist1d.push_back(&Tau_pt_wo_FLSigmaCut);
	Extradist1d.push_back(&Tau_phi_wo_FLSigmaCut);
	Extradist1d.push_back(&Tau_eta_wo_FLSigmaCut);
	Extradist1d.push_back(&Tau_Mass_Inclusive);
	Extradist1d.push_back(&Tau_Mass_sq_Inclusive);
	Extradist1d.push_back(&Tau_Mass_Inclusive_UnFitTracks);
	Extradist1d.push_back(&Tau_Mass_Inclusive_ReFitTracks);
	Extradist1d.push_back(&Tau_Mass_Difference_PFTau_UnFitTracks_3PS);
	Extradist1d.push_back(&Tau_Mass_Difference_PFTau_ReFitTracks_3PS);
	Extradist1d.push_back(&MET_phi);
	Extradist1d.push_back(&TauFL_NoTauFLSigmaCut);
	Extradist1d.push_back(&TauFLSigned_NoTauFLSigmaCut);
	Extradist1d.push_back(&TauFLSigmaSigned);
	Extradist1d.push_back(&TauFLSigmaUnsigned);
	Extradist1d.push_back(&A1mass);
	Extradist1d.push_back(&A1mass10GeV);
	Extradist1d.push_back(&A1massRefit);
	Extradist1d.push_back(&dA1mass_PFTau_Refit);
	Extradist1d.push_back(&A1_Phi_Res);
	Extradist1d.push_back(&A1_Theta_Res);

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
	Extradist1d.push_back(&dPhi_GenTauMu_GenMu);
	Extradist1d.push_back(&dPhi_GenTauMu_RecoMu);
	Extradist1d.push_back(&dTheta_GenTauMu_GenMu);
	Extradist1d.push_back(&dTheta_GenTauMu_RecoMu);

	Extradist1d.push_back(&GJAngle_Over_GJAngleMax_StraightTau);
	Extradist1d.push_back(&GJAngle_Over_GJAngleMax_HelixTau);
	Extradist1d.push_back(&dGJAngle_GJAngleMAX_StraightTau);
	Extradist1d.push_back(&dGJAngle_GJAngleMAX_HelixTau);
	Extradist1d.push_back(&dGJAngle_HelixTau_StraightTau);
	Extradist1d.push_back(&dGJAngle_HelixTau_StraightTauOverGJAngle);
	Extradist1d.push_back(&Angle_HelixTau_StraightTau);
	Extradist1d.push_back(&NUnphysical_StraightTau_HelixTau);
	Extradist1d.push_back(&TauA1_Reco_Solution_StraightTau);
	Extradist1d.push_back(&TauA1_Reco_Solution_HelixTau);

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
	Extradist1d.push_back(&Gen_DiTau_PtBalance_M);
	Extradist1d.push_back(&Gen_DiTau_PtBalance_dM);
	Extradist1d.push_back(&Gen_TauMu_PtBalance_Pt);
	Extradist1d.push_back(&Gen_TauMu_PtBalance_dP);
	Extradist1d.push_back(&Gen_TauA1_dP);
	Extradist1d.push_back(&Gen_TPTF_TauA1_Solution_NoSelection);
	Extradist1d.push_back(&Gen_TPTF_TauA1_Solution_WithSelection);
	Extradist2d.push_back(&Gen_Z_Pt_vs_MET);
	Extradist2d.push_back(&Gen_Z_Pt_vs_VtxTracksPt);
	Extradist2d.push_back(&Gen_Z_Phi_vs_VtxTracksPhi);
	Extradist1d.push_back(&Gen_TauA1_P);
	Extradist1d.push_back(&Gen_TauA1_Pt);
	Extradist1d.push_back(&Gen_TauA1_Px);
	Extradist1d.push_back(&Gen_TauA1_Py);
	Extradist1d.push_back(&Gen_TauA1_Pz);
	Extradist1d.push_back(&Gen_TauA1_Phi);
	Extradist1d.push_back(&Gen_TauA1_Eta);
	Extradist1d.push_back(&Gen_TauMu_P);
	Extradist1d.push_back(&Gen_TauMu_Pt);
	Extradist1d.push_back(&Gen_TauMu_Px);
	Extradist1d.push_back(&Gen_TauMu_Py);
	Extradist1d.push_back(&Gen_TauMu_Pz);
	Extradist1d.push_back(&Gen_TauMu_Phi);
	Extradist1d.push_back(&Gen_TauMu_Eta);
	Extradist1d.push_back(&Gen_TauA1_P_noSel);
	Extradist1d.push_back(&Gen_TauA1_Pt_noSel);
	Extradist1d.push_back(&Gen_TauA1_Px_noSel);
	Extradist1d.push_back(&Gen_TauA1_Py_noSel);
	Extradist1d.push_back(&Gen_TauA1_Pz_noSel);
	Extradist1d.push_back(&Gen_TauA1_Phi_noSel);
	Extradist1d.push_back(&Gen_TauA1_Eta_noSel);
	Extradist1d.push_back(&Gen_TauMu_P_noSel);
	Extradist1d.push_back(&Gen_TauMu_Pt_noSel);
	Extradist1d.push_back(&Gen_TauMu_Px_noSel);
	Extradist1d.push_back(&Gen_TauMu_Py_noSel);
	Extradist1d.push_back(&Gen_TauMu_Pz_noSel);
	Extradist1d.push_back(&Gen_TauMu_Phi_noSel);
	Extradist1d.push_back(&Gen_TauMu_Eta_noSel);
	Extradist1d.push_back(&Gen_Z_Pt);
	Extradist1d.push_back(&Gen_Z_M);
	Extradist1d.push_back(&Gen_Z_Phi);
	Extradist1d.push_back(&Gen_Z_Eta);
	Extradist1d.push_back(&Gen_Z_Pt_noSel);
	Extradist1d.push_back(&Gen_Z_Phi_noSel);
	Extradist1d.push_back(&Gen_Z_Eta_noSel);

	Extradist1d.push_back(&VtxTracksPtRes);
	Extradist1d.push_back(&VtxTracksPhiCorrectedRes);
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
	Extradist1d.push_back(&TPTF_Neutrino_RefitPFTau_HelixGenTau_Mass);
	Extradist1d.push_back(&TPTF_Neutrino_GenA1_StraightGenTau_Mass);
	Extradist1d.push_back(&TPTF_Neutrino_GenA1_HelixGenTau_Mass);

	Extradist1d.push_back(&TransTrk_Failure_withSelection);
	Extradist1d.push_back(&TransTrk_Failure_noSelection);

	Extradist1d.push_back(&Est_Z_Pt_wTruth);
	Extradist1d.push_back(&Est_Z_PtRes_wTruth);
	Extradist2d.push_back(&Est_Z_Pt_wTruth_vs_GenZ_Pt);
	Extradist1d.push_back(&Est_Z_Pt_alwaysMinus);
	Extradist1d.push_back(&Est_Z_PtRes_alwaysMinus);
	Extradist2d.push_back(&Est_Z_Pt_alwaysMinus_vs_GenZ_Pt);
	Extradist1d.push_back(&Est_Z_Energy_wTruth);
	Extradist1d.push_back(&Est_Z_Energy_alwaysMinus);
	Extradist1d.push_back(&Est_Z_EnergyRes_wTruth);
	Extradist1d.push_back(&Est_Z_EnergyRes_alwaysMinus);
	Extradist1d.push_back(&Est_Z_EnergyRes_wTruth2);
	Extradist1d.push_back(&Est_Z_EnergyRes_alwaysMinus2);
	Extradist1d.push_back(&Est_TauMu_PtRes_wTruth);
	Extradist1d.push_back(&Est_TauMu_PtRes_wTruth2);
	Extradist1d.push_back(&Est_Z_M_wTruth2);
	Extradist1d.push_back(&Est_TauMu_PtRes_wTruth_NoZMass);
	Extradist1d.push_back(&Est_Z_M_wTruth_NoZMass);
	Extradist1d.push_back(&Est_Z_EnergyRes_wTruth_NoZMass);
	Extradist1d.push_back(&Est_TauMu_wMET_PtRes);
	Extradist1d.push_back(&Est_TauMu_wMET_PhiRes);
	Extradist1d.push_back(&Est_TauMu_wMET_EtaRes);
	Extradist1d.push_back(&Est_Z_wMET_PtRes);
	Extradist1d.push_back(&Est_Z_wMET_PhiRes);

	Extradist1d.push_back(&Reco_ZMass);
	Extradist1d.push_back(&Reco_ZMass_UnboostedGenZ);
	Extradist1d.push_back(&Reco_ZMass_MassScan);
	Extradist1d.push_back(&Reco_ZMass_MassScanUnboosted);
	Extradist1d.push_back(&Reco_ZMasswithProbWeight_MassScan);
	Extradist1d.push_back(&Reco_ProbStack_MassScan);
	Extradist1d.push_back(&Reco_EventFit_Solution);
	Extradist1d.push_back(&Reco_A1Fit_Solution);
	Extradist1d.push_back(&Reco_Chi2_FitSolutionOnlyLargeScale);
	Extradist1d.push_back(&Reco_ConstrainedDeltaSum);
	Extradist1d.push_back(&Reco_ConstrainedDeltaMass);
	Extradist1d.push_back(&Reco_ConstrainedDeltaPt);

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

	Extradist1d.push_back(&Reco_dPhi_TauMuTauA1_AfterFit_lowBoost);
	Extradist1d.push_back(&Reco_dPhi_TauMuTauA1_BeforeFit_lowBoost);

	Extradist1d.push_back(&Reco_PtRes_TauA1_NoFit);
	Extradist1d.push_back(&Reco_PtRes_TauA1_AmbPoint0_NoFit);
	Extradist1d.push_back(&Reco_PtRes_TauA1_AmbPoint12_NoFit);

	Extradist1d.push_back(&Reco_PtRes_TauMu_NoFit);
	Extradist1d.push_back(&Reco_PtRes_TauMu_AmbPoint0_NoFit);
	Extradist1d.push_back(&Reco_PtRes_TauMu_AmbPoint12_NoFit);

	Extradist1d.push_back(&Reco_PhiRes_TauMu_PreFit);
	Extradist1d.push_back(&Reco_ThetaRes_TauMu_PreFit);

	Extradist1d.push_back(&Reco_TauMu_DeltaPX_FitImpact);
	Extradist1d.push_back(&Reco_TauMu_DeltaPY_FitImpact);
	Extradist1d.push_back(&Reco_TauMu_DeltaPZ_FitImpact);
	Extradist1d.push_back(&Reco_TauA1_DeltaPX_FitImpact);
	Extradist1d.push_back(&Reco_TauA1_DeltaPY_FitImpact);
	Extradist1d.push_back(&Reco_TauA1_DeltaPZ_FitImpact);
	Extradist1d.push_back(&GenZ_Pt_Unboosted);
	Extradist1d.push_back(&RecoZ_Pt_Unboosted);
	Extradist1d.push_back(&Reco_TauMu_ResCosTheta);
	Extradist1d.push_back(&Reco_TauMu_DeltaPhi_FitImpact);

	Extradist1d.push_back(&Reco_ZMass_PDF);
	Extradist1d.push_back(&Reco_Z_Energy_Res);
	Extradist1d.push_back(&RecoZ_Pt);

	Extradist1d.push_back(&Reco_TauA1_P);
	Extradist1d.push_back(&Reco_TauA1_Pt);
	Extradist1d.push_back(&Reco_TauA1_Px);
	Extradist1d.push_back(&Reco_TauA1_Py);
	Extradist1d.push_back(&Reco_TauA1_Pz);
	Extradist1d.push_back(&Reco_TauA1_Phi);
	Extradist1d.push_back(&Reco_TauA1_Eta);
	Extradist1d.push_back(&Reco_TauMu_P);
	Extradist1d.push_back(&Reco_TauMu_Pt);
	Extradist1d.push_back(&Reco_TauMu_Px);
	Extradist1d.push_back(&Reco_TauMu_Py);
	Extradist1d.push_back(&Reco_TauMu_Pz);
	Extradist1d.push_back(&Reco_TauMu_Phi);
	Extradist1d.push_back(&Reco_TauMu_Eta);
	Extradist1d.push_back(&Reco_Z_Px);
	Extradist1d.push_back(&Reco_Z_Py);
	Extradist1d.push_back(&Reco_Z_Pz);
	Extradist1d.push_back(&Reco_Z_Pt);
	Extradist1d.push_back(&Reco_Z_PtRes);
	Extradist1d.push_back(&Reco_Z_Phi);
	Extradist1d.push_back(&Reco_Z_Eta);

	Extradist1d.push_back(&Reco_Z_PhiRes);
	Extradist1d.push_back(&Reco_Z_PhiRes_noAmb);
	Extradist1d.push_back(&Reco_Z_PhiRes_wAmb);
	Extradist1d.push_back(&Reco_Z_PhiRes_pickedrightAmb);
	Extradist1d.push_back(&Reco_Z_PhiRes_pickedwrongAmb);
	Extradist1d.push_back(&Reco_Z_EtaRes);
	Extradist1d.push_back(&Reco_Z_EtaRes_noAmb);
	Extradist1d.push_back(&Reco_Z_EtaRes_wAmb);
	Extradist1d.push_back(&Reco_Z_EtaRes_pickedrightAmb);
	Extradist1d.push_back(&Reco_Z_EtaRes_pickedwrongAmb);
	Extradist1d.push_back(&Reco_Z_PRes);
	Extradist1d.push_back(&Reco_Z_PRes_noAmb);
	Extradist1d.push_back(&Reco_Z_PRes_wAmb);
	Extradist1d.push_back(&Reco_Z_PRes_pickedrightAmb);
	Extradist1d.push_back(&Reco_Z_PRes_pickedwrongAmb);
	Extradist1d.push_back(&Reco_NIter);
	Extradist1d.push_back(&Reco_NIter_noAmb);
	Extradist1d.push_back(&Reco_NIter_wAmb);
	Extradist1d.push_back(&Reco_NIter_pickedrightAmb);
	Extradist1d.push_back(&Reco_NIter_pickedwrongAmb);
	Extradist1d.push_back(&Reco_Chi2);
	Extradist1d.push_back(&Reco_Chi2_noAmb);
	Extradist1d.push_back(&Reco_Chi2_wAmb);
	Extradist1d.push_back(&Reco_Chi2_rightAmb);
	Extradist1d.push_back(&Reco_Chi2_wrongAmb);
	Extradist1d.push_back(&Reco_Chi2_pickedrightAmb);
	Extradist1d.push_back(&Reco_Chi2_pickedwrongAmb);
	Extradist1d.push_back(&Reco_Chi2_orig);
	Extradist1d.push_back(&Reco_Chi2_orig_noAmb);
	Extradist1d.push_back(&Reco_Chi2_orig_wAmb);
	Extradist1d.push_back(&Reco_Chi2_orig_rightAmb);
	Extradist1d.push_back(&Reco_Chi2_orig_wrongAmb);
	Extradist1d.push_back(&Reco_Chi2_orig_pickedrightAmb);
	Extradist1d.push_back(&Reco_Chi2_orig_pickedwrongAmb);
	Extradist1d.push_back(&Reco_Chi2_SC);
	Extradist1d.push_back(&Reco_Chi2_SC_noAmb);
	Extradist1d.push_back(&Reco_Chi2_SC_wAmb);
	Extradist1d.push_back(&Reco_Chi2_SC_rightAmb);
	Extradist1d.push_back(&Reco_Chi2_SC_wrongAmb);
	Extradist1d.push_back(&Reco_Chi2_SC_pickedrightAmb);
	Extradist1d.push_back(&Reco_Chi2_SC_pickedwrongAmb);
	Extradist1d.push_back(&Reco_Chi2_HC);
	Extradist1d.push_back(&Reco_Chi2_HC_noAmb);
	Extradist1d.push_back(&Reco_Chi2_HC_wAmb);
	Extradist1d.push_back(&Reco_Chi2_HC_rightAmb);
	Extradist1d.push_back(&Reco_Chi2_HC_wrongAmb);
	Extradist1d.push_back(&Reco_Chi2_HC_pickedrightAmb);
	Extradist1d.push_back(&Reco_Chi2_HC_pickedwrongAmb);

	Extradist1d.push_back(&Reco_Chi2_diff);
	Extradist1d.push_back(&Reco_Chi2_orig_diff);
	Extradist1d.push_back(&Reco_Chi2_HC_diff);
	Extradist1d.push_back(&Reco_Chi2_SC_diff);
	Extradist1d.push_back(&Reco_Chi2_origplusSC_diff);

	Extradist2d.push_back(&Reco_Chi2_diff_vs_correctAssignment);

	Extradist2d.push_back(&Reco_FitSolution_byChi2_Full_vs_RightSolution);
	Extradist2d.push_back(&Reco_FitSolution_byChi2_orig_vs_RightSolution);
	Extradist2d.push_back(&Reco_FitSolution_byChi2_SC_vs_RightSolution);
	Extradist2d.push_back(&Reco_FitSolution_byChi2_HC_vs_RightSolution);
	Extradist2d.push_back(&Reco_FitSolution_byChi2_origplusSC_vs_RightSolution);

	Extradist1d.push_back(&Reco_PtRes_TauA1_wRecoil);
	Extradist1d.push_back(&Reco_PtRes_TauA1_wRecoil_AmbZero);
	Extradist1d.push_back(&Reco_PtRes_TauA1_wRecoil_wAmb);
	Extradist1d.push_back(&Reco_PtRes_TauMu_wRecoil);
	Extradist1d.push_back(&Reco_PtRes_TauMu_wRecoil_AmbZero);
	Extradist1d.push_back(&Reco_PtRes_TauMu_wRecoil_wAmb);
	Extradist1d.push_back(&Reco_PtRes_TauA1_wRecoil_PreFit);
	Extradist1d.push_back(&Reco_PtRes_TauMu_wRecoil_PreFit);
	Extradist1d.push_back(&Reco_dPhi_TauMuTauA1_AfterFit_wRecoil);
	Extradist1d.push_back(&Reco_dPhi_TauMuTauA1_BeforeFit_wRecoil);
	Extradist1d.push_back(&Reco_Chi2_FitSolutionOnly_wRecoil);
	Extradist1d.push_back(&Reco_Chi2_Full_wRecoil);
	Extradist1d.push_back(&Reco_Chi2_Orig_wRecoil);
	Extradist1d.push_back(&Reco_Chi2_SC_wRecoil);
	Extradist1d.push_back(&Reco_Chi2_HC_wRecoil);
	Extradist1d.push_back(&Reco_Chi2_OrigProb_wRecoil);
	Extradist1d.push_back(&Reco_Chi2_diff_wRecoil);
	Extradist1d.push_back(&Reco_Chi2_orig_diff_wRecoil);
	Extradist1d.push_back(&Reco_ZMass_wRecoil);
	Extradist1d.push_back(&Reco_NIter_wRecoil);
	Extradist1d.push_back(&Reco_EventFit_Solution_wRecoil);
	Extradist1d.push_back(&Reco_Z_Energy_Res_wRecoil);
	Extradist1d.push_back(&Reco_TauA1_P_wRecoil);
	Extradist1d.push_back(&Reco_TauA1_Pt_wRecoil);
	Extradist1d.push_back(&Reco_TauA1_Px_wRecoil);
	Extradist1d.push_back(&Reco_TauA1_Py_wRecoil);
	Extradist1d.push_back(&Reco_TauA1_Pz_wRecoil);
	Extradist1d.push_back(&Reco_TauA1_Phi_wRecoil);
	Extradist1d.push_back(&Reco_TauA1_Eta_wRecoil);
	Extradist1d.push_back(&Reco_TauMu_P_wRecoil);
	Extradist1d.push_back(&Reco_TauMu_Pt_wRecoil);
	Extradist1d.push_back(&Reco_TauMu_Px_wRecoil);
	Extradist1d.push_back(&Reco_TauMu_Py_wRecoil);
	Extradist1d.push_back(&Reco_TauMu_Pz_wRecoil);
	Extradist1d.push_back(&Reco_TauMu_Phi_wRecoil);
	Extradist1d.push_back(&Reco_TauMu_Eta_wRecoil);
	Extradist1d.push_back(&Reco_Z_Px_wRecoil);
	Extradist1d.push_back(&Reco_Z_Py_wRecoil);
	Extradist1d.push_back(&Reco_Z_Pz_wRecoil);
	Extradist1d.push_back(&Reco_Z_Pt_wRecoil);
	Extradist1d.push_back(&Reco_Z_PtRes_wRecoil);
	Extradist1d.push_back(&Reco_Z_Phi_wRecoil);
	Extradist1d.push_back(&Reco_Z_PhiRes_wRecoil);
	Extradist1d.push_back(&Reco_Z_Eta_wRecoil);
	Extradist1d.push_back(&Reco_Z_EtaRes_wRecoil);

	Extradist2d.push_back(&Reco_TauA1_ptRes_vs_ptGen_wRecoil);
	Extradist2d.push_back(&Reco_TauA1_ptRes_vs_ptReco_wRecoil);
	Extradist2d.push_back(&Reco_TauMu_ptRes_vs_ptGen_wRecoil);
	Extradist2d.push_back(&Reco_TauMu_ptRes_vs_ptReco_wRecoil);

	Extradist2d.push_back(&Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr);
	Extradist2d.push_back(&Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr);

	Extradist2d.push_back(&Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_wAmb);
	Extradist2d.push_back(&Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_AmbZero);
	Extradist2d.push_back(&Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_wAmb);
	Extradist2d.push_back(&Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_AmbZero);

	Extradist2d.push_back(&Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit);
	Extradist2d.push_back(&Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit);
	Extradist2d.push_back(&Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit);
	Extradist2d.push_back(&Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit);
	Extradist2d.push_back(&Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit);
	Extradist2d.push_back(&Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit);

	Extradist2d.push_back(&Reco_FitSolution_byChi2_Full_vs_RightSolution_wRecoil);
	Extradist2d.push_back(&Reco_FitSolution_byChi2_orig_vs_RightSolution_wRecoil);
	Extradist2d.push_back(&Reco_FitSolution_byChi2_SC_vs_RightSolution_wRecoil);
	Extradist2d.push_back(&Reco_FitSolution_byChi2_HC_vs_RightSolution_wRecoil);
	Extradist2d.push_back(&Reco_FitSolution_byChi2_origplusSC_vs_RightSolution_wRecoil);

	Extradist1d.push_back(&GenReco_ZMass);
	Extradist1d.push_back(&GenReco_EventFit_Solution);
	Extradist1d.push_back(&GenReco_A1Fit_Solution);
	Extradist1d.push_back(&GenReco_Chi2);
	Extradist1d.push_back(&GenReco_Chi2_FitSolutionOnly);
	Extradist1d.push_back(&GenReco_ConstrainedDeltaSum);
	Extradist1d.push_back(&GenReco_NIter);

	Extradist1d.push_back(&TauMu_Start_Collinear_PtRes);
	Extradist1d.push_back(&TauMu_Start_MET_PtRes);
	Extradist1d.push_back(&TauMu_Start_PtBalance_PtRes);
	Extradist1d.push_back(&TauMu_Start_EventRecoil_PtRes);
	Extradist2d.push_back(&TauMu_Start_Collinear_PtReco_vs_PtGen);
	Extradist2d.push_back(&TauMu_Start_MET_PtReco_vs_PtGen);
	Extradist2d.push_back(&TauMu_Start_PtBalance_PtReco_vs_PtGen);
	Extradist2d.push_back(&TauMu_Start_EventRecoil_PtReco_vs_PtGen);

	Extradist1d.push_back(&Z_Start_MET_PtRes);
	Extradist1d.push_back(&Z_Start_MET_PhiRes);
	Extradist1d.push_back(&Z_Start_PtBalance_PtRes);
	Extradist1d.push_back(&Z_Start_PtBalance_PhiRes);
	Extradist1d.push_back(&Z_Start_EventRecoil_PtRes);
	Extradist1d.push_back(&Z_Start_EventRecoil_PhiRes);

	Extradist1d.push_back(&Z_Start_MET_noNu_PtRes);
	Extradist1d.push_back(&Z_Start_MET_noNu_PhiRes);
	Extradist1d.push_back(&Z_Start_MET_noNu_M);
	Extradist1d.push_back(&Z_Start_MET_noNu_dE);

	Extradist1d.push_back(&Z_Start_MET_noNu_PtRes_withRot);
	Extradist1d.push_back(&Z_Start_MET_noNu_PhiRes_withRot);
	Extradist1d.push_back(&Z_Start_MET_noNu_M_withRot);

	Extradist1d.push_back(&Z_Start_B2B_M);
	Extradist1d.push_back(&Z_Start_B2B_dE);

	Extradist1d.push_back(&TauMu_Start_MET_PtRes_AfterMC);
	Extradist1d.push_back(&TauMu_Start_PtBalance_PtRes_AfterMC);
	Extradist1d.push_back(&TauMu_Start_EventRecoil_PtRes_AfterMC);

	Extradist1d.push_back(&TauMu_Start_dPhi_TauMuTauH_MET);
	Extradist1d.push_back(&TauMu_Start_dPhi_TauMuTauH_EventRecoil);
	Extradist1d.push_back(&TauMu_Start_dPhi_TauMuTauH_PtBalance);

	Extradist1d.push_back(&Mu_TP_phi0);
	Extradist1d.push_back(&Mu_TP_lambda);
	Extradist1d.push_back(&Mu_TP_dxy);
	Extradist1d.push_back(&Mu_TP_dz);
	Extradist1d.push_back(&Mu_TP_kappa);
	Extradist1d.push_back(&Mu_TP_POCA_quadrant);
	Extradist1d.push_back(&Mu_TP_POCA_quadrantVlad);
	Extradist1d.push_back(&Mu_TP_POCA_quadrantby_dxyphi0);

	Extradist2d.push_back(&Mu_TP_Poca_xy);
	Extradist2d.push_back(&Mu_TP_Vertex_xy);
	Extradist2d.push_back(&Mu_TP_RefitVertex_xy);
	Extradist2d.push_back(&Mu_TP_BeamSpot_xy);
	Extradist2d.push_back(&Mu_TP_NTP_Poca_xy);

	Extradist1d.push_back(&MVAMET_metobject_XX);
	Extradist1d.push_back(&MVAMET_ptobject_XX);

	Extradist1d.push_back(&A1_Par_Px);
	Extradist1d.push_back(&A1_Par_Py);
	Extradist1d.push_back(&A1_Par_Pz);
	Extradist1d.push_back(&A1_Par_M);
	Extradist1d.push_back(&A1_Cov_Pxx);
	Extradist1d.push_back(&A1_Cov_Pyy);
	Extradist1d.push_back(&A1_Cov_Pzz);
	Extradist1d.push_back(&A1_Cov_Pxy);
	Extradist1d.push_back(&A1_Cov_Pxz);
	Extradist1d.push_back(&A1_Cov_Pyz);
	Extradist1d.push_back(&A1_Pull_Px);
	Extradist1d.push_back(&A1_Pull_Py);
	Extradist1d.push_back(&A1_Pull_Pz);
	Extradist1d.push_back(&A1_Pull_M);
	Extradist1d.push_back(&A1_Pull_Px_AmbZero);
	Extradist1d.push_back(&A1_Pull_Py_AmbZero);
	Extradist1d.push_back(&A1_Pull_Pz_AmbZero);
	Extradist1d.push_back(&A1_Pull_M_AmbZero);
	Extradist1d.push_back(&A1_Pull_Px_wAmb);
	Extradist1d.push_back(&A1_Pull_Py_wAmb);
	Extradist1d.push_back(&A1_Pull_Pz_wAmb);
	Extradist1d.push_back(&A1_Pull_M_wAmb);

	Extradist1d.push_back(&SV_Par_x);
	Extradist1d.push_back(&SV_Par_y);
	Extradist1d.push_back(&SV_Par_z);
	Extradist1d.push_back(&SV_Cov_xx);
	Extradist1d.push_back(&SV_Cov_yy);
	Extradist1d.push_back(&SV_Cov_zz);
	Extradist1d.push_back(&SV_Cov_xy);
	Extradist1d.push_back(&SV_Cov_xz);
	Extradist1d.push_back(&SV_Cov_yz);
	Extradist1d.push_back(&SV_Pull_Px);
	Extradist1d.push_back(&SV_Pull_Py);
	Extradist1d.push_back(&SV_Pull_Pz);
	Extradist1d.push_back(&SV_Pull_Px_AmbZero);
	Extradist1d.push_back(&SV_Pull_Py_AmbZero);
	Extradist1d.push_back(&SV_Pull_Pz_AmbZero);
	Extradist1d.push_back(&SV_Pull_Px_wAmb);
	Extradist1d.push_back(&SV_Pull_Py_wAmb);
	Extradist1d.push_back(&SV_Pull_Pz_wAmb);

	Extradist1d.push_back(&TauA1_Par_Px_AmbZero);
	Extradist1d.push_back(&TauA1_Par_Py_AmbZero);
	Extradist1d.push_back(&TauA1_Par_Pz_AmbZero);
	Extradist1d.push_back(&TauA1_Cov_Pxx_AmbZero);
	Extradist1d.push_back(&TauA1_Cov_Pyy_AmbZero);
	Extradist1d.push_back(&TauA1_Cov_Pzz_AmbZero);
	Extradist1d.push_back(&TauA1_Cov_Pxy_AmbZero);
	Extradist1d.push_back(&TauA1_Cov_Pxz_AmbZero);
	Extradist1d.push_back(&TauA1_Cov_Pyz_AmbZero);
	Extradist1d.push_back(&TauA1_Pull_Px_AmbZero);
	Extradist1d.push_back(&TauA1_Pull_Py_AmbZero);
	Extradist1d.push_back(&TauA1_Pull_Pz_AmbZero);

	Extradist1d.push_back(&TauA1_Par_Px_CorrectAmb);
	Extradist1d.push_back(&TauA1_Par_Py_CorrectAmb);
	Extradist1d.push_back(&TauA1_Par_Pz_CorrectAmb);
	Extradist1d.push_back(&TauA1_Cov_Pxx_CorrectAmb);
	Extradist1d.push_back(&TauA1_Cov_Pyy_CorrectAmb);
	Extradist1d.push_back(&TauA1_Cov_Pzz_CorrectAmb);
	Extradist1d.push_back(&TauA1_Cov_Pxy_CorrectAmb);
	Extradist1d.push_back(&TauA1_Cov_Pxz_CorrectAmb);
	Extradist1d.push_back(&TauA1_Cov_Pyz_CorrectAmb);
	Extradist1d.push_back(&TauA1_Pull_Px_CorrectAmb);
	Extradist1d.push_back(&TauA1_Pull_Py_CorrectAmb);
	Extradist1d.push_back(&TauA1_Pull_Pz_CorrectAmb);

	Extradist1d.push_back(&TauA1_Par_Px_WrongAmb);
	Extradist1d.push_back(&TauA1_Par_Py_WrongAmb);
	Extradist1d.push_back(&TauA1_Par_Pz_WrongAmb);
	Extradist1d.push_back(&TauA1_Cov_Pxx_WrongAmb);
	Extradist1d.push_back(&TauA1_Cov_Pyy_WrongAmb);
	Extradist1d.push_back(&TauA1_Cov_Pzz_WrongAmb);
	Extradist1d.push_back(&TauA1_Cov_Pxy_WrongAmb);
	Extradist1d.push_back(&TauA1_Cov_Pxz_WrongAmb);
	Extradist1d.push_back(&TauA1_Cov_Pyz_WrongAmb);
	Extradist1d.push_back(&TauA1_Pull_Px_WrongAmb);
	Extradist1d.push_back(&TauA1_Pull_Py_WrongAmb);
	Extradist1d.push_back(&TauA1_Pull_Pz_WrongAmb);

	Extradist1d.push_back(&TauFLSigmaVlad);
	Extradist1d.push_back(&TauFLSigmaVlad_PhiA1);
	Extradist1d.push_back(&TauFLSigmaVlad_PhiTau);
	Extradist1d.push_back(&TauFLSigmaVlad_PhiTauNoCorr);
	Extradist1d.push_back(&TauFLSigmaVlad_PhiZnoCorr);
	Extradist1d.push_back(&TauFLSigmaVlad_PhiZwCorr);

	Extradist1d.push_back(&TauFLSigmaAlex);
	Extradist1d.push_back(&TauFLSigmaAlex_PhiA1);
	Extradist1d.push_back(&TauFLSigmaAlex_PhiTau);
	Extradist1d.push_back(&TauFLSigmaAlex_PhiTauNoCorr);
	Extradist1d.push_back(&TauFLSigmaAlex_PhiZnoCorr);
	Extradist1d.push_back(&TauFLSigmaAlex_PhiZwCorr);

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
	Extradist1d_OS.push_back(&Tau_Mass_Inclusive_UnFitTracks);
	Extradist1d_OS.push_back(&Tau_Mass_Inclusive_ReFitTracks);
	Extradist1d_OS.push_back(&MET_phi);
	Extradist1d_OS.push_back(&TauFL_NoTauFLSigmaCut);
	Extradist1d_OS.push_back(&TauFLSigned_NoTauFLSigmaCut);
	Extradist1d_OS.push_back(&TauFLSigmaSigned);
	Extradist1d_OS.push_back(&TauFLSigmaUnsigned);
	Extradist1d_OS.push_back(&A1mass);
	Extradist1d_OS.push_back(&A1mass10GeV);
	Extradist1d_OS.push_back(&A1massRefit);
	Extradist1d_OS.push_back(&dA1mass_PFTau_Refit);
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
	Extradist1d_OS.push_back(&Reco_Chi2_FitSolutionOnlyLargeScale);
	Extradist1d_OS.push_back(&Reco_ConstrainedDeltaSum);
	Extradist1d_OS.push_back(&Reco_ConstrainedDeltaMass);
	Extradist1d_OS.push_back(&Reco_ConstrainedDeltaPt);
	Extradist1d_OS.push_back(&Reco_NIter);
	Extradist1d_OS.push_back(&Reco_Chi2);
	Extradist1d_OS.push_back(&Reco_Chi2_orig);
	Extradist1d_OS.push_back(&Reco_Chi2_SC);
	Extradist1d_OS.push_back(&Reco_Chi2_HC);
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

	Extradist1d_OS.push_back(&Reco_TauA1_P);
	Extradist1d_OS.push_back(&Reco_TauA1_Pt);
	Extradist1d_OS.push_back(&Reco_TauA1_Px);
	Extradist1d_OS.push_back(&Reco_TauA1_Py);
	Extradist1d_OS.push_back(&Reco_TauA1_Pz);
	Extradist1d_OS.push_back(&Reco_TauA1_Phi);
	Extradist1d_OS.push_back(&Reco_TauA1_Eta);
	Extradist1d_OS.push_back(&Reco_TauMu_P);
	Extradist1d_OS.push_back(&Reco_TauMu_Pt);
	Extradist1d_OS.push_back(&Reco_TauMu_Px);
	Extradist1d_OS.push_back(&Reco_TauMu_Py);
	Extradist1d_OS.push_back(&Reco_TauMu_Pz);
	Extradist1d_OS.push_back(&Reco_TauMu_Phi);
	Extradist1d_OS.push_back(&Reco_TauMu_Eta);
	Extradist1d_OS.push_back(&Reco_Z_Px);
	Extradist1d_OS.push_back(&Reco_Z_Py);
	Extradist1d_OS.push_back(&Reco_Z_Pz);
	Extradist1d_OS.push_back(&Reco_Z_Pt);
	Extradist1d_OS.push_back(&Reco_Z_Phi);
	Extradist1d_OS.push_back(&Reco_Z_Eta);

	Extradist1d_OS.push_back(&Mu_TP_phi0);
	Extradist1d_OS.push_back(&Mu_TP_lambda);
	Extradist1d_OS.push_back(&Mu_TP_dxy);
	Extradist1d_OS.push_back(&Mu_TP_dz);
	Extradist1d_OS.push_back(&Mu_TP_kappa);
	Extradist1d_OS.push_back(&Mu_TP_POCA_quadrant);
	Extradist1d_OS.push_back(&Mu_TP_POCA_quadrantVlad);
	Extradist1d_OS.push_back(&Mu_TP_POCA_quadrantby_dxyphi0);

	Extradist1d_OS.push_back(&MVAMET_metobject_XX);
	Extradist1d_OS.push_back(&MVAMET_ptobject_XX);
	Extradist1d_OS.push_back(&A1_Par_Px);
	Extradist1d_OS.push_back(&A1_Par_Py);
	Extradist1d_OS.push_back(&A1_Par_Pz);
	Extradist1d_OS.push_back(&A1_Par_M);
	Extradist1d_OS.push_back(&A1_Cov_Pxx);
	Extradist1d_OS.push_back(&A1_Cov_Pyy);
	Extradist1d_OS.push_back(&A1_Cov_Pzz);
	Extradist1d_OS.push_back(&A1_Cov_Pxy);
	Extradist1d_OS.push_back(&A1_Cov_Pxz);
	Extradist1d_OS.push_back(&A1_Cov_Pyz);
	Extradist1d_OS.push_back(&SV_Par_x);
	Extradist1d_OS.push_back(&SV_Par_y);
	Extradist1d_OS.push_back(&SV_Par_z);
	Extradist1d_OS.push_back(&SV_Cov_xx);
	Extradist1d_OS.push_back(&SV_Cov_yy);
	Extradist1d_OS.push_back(&SV_Cov_zz);
	Extradist1d_OS.push_back(&SV_Cov_xy);
	Extradist1d_OS.push_back(&SV_Cov_xz);
	Extradist1d_OS.push_back(&SV_Cov_yz);

	Extradist1d_OS.push_back(&TauA1_Par_Px_AmbZero);
	Extradist1d_OS.push_back(&TauA1_Par_Py_AmbZero);
	Extradist1d_OS.push_back(&TauA1_Par_Pz_AmbZero);
	Extradist1d_OS.push_back(&TauA1_Cov_Pxx_AmbZero);
	Extradist1d_OS.push_back(&TauA1_Cov_Pyy_AmbZero);
	Extradist1d_OS.push_back(&TauA1_Cov_Pzz_AmbZero);
	Extradist1d_OS.push_back(&TauA1_Cov_Pxy_AmbZero);
	Extradist1d_OS.push_back(&TauA1_Cov_Pxz_AmbZero);
	Extradist1d_OS.push_back(&TauA1_Cov_Pyz_AmbZero);

	Extradist1d_OS.push_back(&TauA1_Par_Px_WrongAmb);
	Extradist1d_OS.push_back(&TauA1_Par_Py_WrongAmb);
	Extradist1d_OS.push_back(&TauA1_Par_Pz_WrongAmb);
	Extradist1d_OS.push_back(&TauA1_Cov_Pxx_WrongAmb);
	Extradist1d_OS.push_back(&TauA1_Cov_Pyy_WrongAmb);
	Extradist1d_OS.push_back(&TauA1_Cov_Pzz_WrongAmb);
	Extradist1d_OS.push_back(&TauA1_Cov_Pxy_WrongAmb);
	Extradist1d_OS.push_back(&TauA1_Cov_Pxz_WrongAmb);
	Extradist1d_OS.push_back(&TauA1_Cov_Pyz_WrongAmb);

	Extradist1d_OS.push_back(&Reco_TauA1_P_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauA1_Pt_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauA1_Px_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauA1_Py_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauA1_Pz_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauA1_Phi_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauA1_Eta_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauMu_P_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauMu_Pt_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauMu_Px_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauMu_Py_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauMu_Pz_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauMu_Phi_wRecoil);
	Extradist1d_OS.push_back(&Reco_TauMu_Eta_wRecoil);
	Extradist1d_OS.push_back(&Reco_Z_Px_wRecoil);
	Extradist1d_OS.push_back(&Reco_Z_Py_wRecoil);
	Extradist1d_OS.push_back(&Reco_Z_Pz_wRecoil);
	Extradist1d_OS.push_back(&Reco_Z_Pt_wRecoil);
	Extradist1d_OS.push_back(&Reco_Z_Phi_wRecoil);
	Extradist1d_OS.push_back(&Reco_Z_Eta_wRecoil);
	Extradist1d_OS.push_back(&Reco_dPhi_TauMuTauA1_AfterFit_wRecoil);
	Extradist1d_OS.push_back(&Reco_dPhi_TauMuTauA1_BeforeFit_wRecoil);
	Extradist1d_OS.push_back(&Reco_Chi2_FitSolutionOnly_wRecoil);
	Extradist1d_OS.push_back(&Reco_Chi2_Full_wRecoil);
	Extradist1d_OS.push_back(&Reco_Chi2_Orig_wRecoil);
	Extradist1d_OS.push_back(&Reco_Chi2_SC_wRecoil);
	Extradist1d_OS.push_back(&Reco_Chi2_HC_wRecoil);
	Extradist1d_OS.push_back(&Reco_Chi2_OrigProb_wRecoil);
	Extradist1d_OS.push_back(&Reco_Chi2_diff_wRecoil);
	Extradist1d_OS.push_back(&Reco_Chi2_orig_diff_wRecoil);
	Extradist1d_OS.push_back(&Reco_ZMass_wRecoil);
	Extradist1d_OS.push_back(&Reco_NIter_wRecoil);
	Extradist1d_OS.push_back(&Reco_EventFit_Solution_wRecoil);

	Extradist1d_OS.push_back(&TauFLSigmaVlad);
	Extradist1d_OS.push_back(&TauFLSigmaVlad_PhiA1);
	Extradist1d_OS.push_back(&TauFLSigmaVlad_PhiTau);
	Extradist1d_OS.push_back(&TauFLSigmaVlad_PhiTauNoCorr);
	Extradist1d_OS.push_back(&TauFLSigmaVlad_PhiZnoCorr);
	Extradist1d_OS.push_back(&TauFLSigmaVlad_PhiZwCorr);

	Extradist1d_OS.push_back(&TauFLSigmaAlex);
	Extradist1d_OS.push_back(&TauFLSigmaAlex_PhiA1);
	Extradist1d_OS.push_back(&TauFLSigmaAlex_PhiTau);
	Extradist1d_OS.push_back(&TauFLSigmaAlex_PhiTauNoCorr);
	Extradist1d_OS.push_back(&TauFLSigmaAlex_PhiZnoCorr);
	Extradist1d_OS.push_back(&TauFLSigmaAlex_PhiZwCorr);

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
	if(selection_verbose) Logger(Logger::Verbose) << "Cut on good vertex" << std::endl;
	unsigned int nGoodVtx=0;
	for(unsigned int i_vtx=0;i_vtx<Ntp->NVtx();i_vtx++){
		if(Ntp->isGoodVtx(i_vtx)){
			if(selVertex == selVertexDummy) selVertex = i_vtx; // selected vertex = first vertex (highest sum[pT^2]) to fulfill vertex requirements
			nGoodVtx++;
		}
	}
	value.at(PrimeVtx)=nGoodVtx;
	pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at Primevtx: " <<value.at(PrimeVtx) << std::endl;
		Logger(Logger::Verbose) << "pass at Primevtx: " <<pass.at(PrimeVtx) << std::endl;
	}

	// Trigger
	if(selection_verbose) Logger(Logger::Verbose) << "Cut on Trigger" << std::endl;
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
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at TriggerOk: " <<value.at(TriggerOk) << std::endl;
		Logger(Logger::Verbose) << "pass at TriggerOk: " <<pass.at(TriggerOk) << std::endl;
	}

	// Muon cuts
	if(selection_verbose) Logger(Logger::Verbose) << "Cut on MuonID" << std::endl;
	std::vector<int> selectedMuonsId;
	selectedMuonsId.clear();
	for(unsigned i_mu=0;i_mu<Ntp->NMuons();i_mu++){
		if( selectMuon_Id(i_mu,selVertex) ) {
			selectedMuonsId.push_back(i_mu);
		}
	}
	value.at(NMuId)=selectedMuonsId.size();
	pass.at(NMuId)=(value.at(NMuId)>=cut.at(NMuId));
	if(selection_verbose){
		Logger(Logger::Verbose) << "Number of Muons: " << Ntp->NMuons() << std::endl;
		Logger(Logger::Verbose) << "value at NMuId: " <<value.at(NMuId) << std::endl;
		Logger(Logger::Verbose) << "pass at NMuId: " <<pass.at(NMuId) << std::endl;
	}

	if(selection_verbose) Logger(Logger::Verbose) << "Cut on Muon Kinematics" << std::endl;
	std::vector<int> selectedMuonsKin;
	selectedMuonsKin.clear();
	for(std::vector<int>::iterator it_mu = selectedMuonsId.begin(); it_mu != selectedMuonsId.end(); ++it_mu){
		if( selectMuon_Kinematics(*it_mu) ) {
			selectedMuonsKin.push_back(*it_mu);
		}
	}
	value.at(NMuKin)=selectedMuonsKin.size();
	pass.at(NMuKin)=(value.at(NMuKin)>=cut.at(NMuKin));
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at NMuKin: " <<value.at(NMuKin) << std::endl;
		Logger(Logger::Verbose) << "pass at NMuKin: " <<pass.at(NMuKin) << std::endl;
	}

	if(selection_verbose) Logger(Logger::Verbose) << "Cut on Muon Isolation (Iso)" << std::endl;
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
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at NMuIso: " <<value.at(NMuIso) << std::endl;
		Logger(Logger::Verbose) << "pass at NMuIso: " <<pass.at(NMuIso) << std::endl;
	}

	if(selection_verbose) Logger(Logger::Verbose) << "Cut on Muon Isolation (Anti Iso)" << std::endl;
	std::vector<int> selectedMuonsAntiIso;
	selectedMuonsAntiIso.clear();
	for(std::vector<int>::iterator it_mu = selectedMuonsKin.begin(); it_mu != selectedMuonsKin.end(); ++it_mu){
		if(selectMuon_AntiIsolation(*it_mu)){
			if(selMuon_AntiIso == selMuonDummy && selMuon_Iso == selMuonDummy) selMuon_AntiIso = *it_mu;
			selectedMuonsAntiIso.push_back(*it_mu);
		}
	}
	if(selMuon_Iso != selMuonDummy && selMuon_AntiIso != selMuonDummy){
		Logger(Logger::Error) << "CRITICAL: SELECTED MUON PASSED ISOLATION AND ANTI-ISOLATION CUT --> FIX" << std::endl;
		return;
	}
	if(id == DataMCType::Data && selMuon_AntiIso != selMuonDummy) selMuon = selMuon_AntiIso;
	else selMuon = selMuon_Iso;

	//Di muon veto
	if(selection_verbose) Logger(Logger::Verbose) << "Di Muon Veto" << std::endl;
	std::vector<int> PosMu_DiMuVeto;
	PosMu_DiMuVeto.clear();
	std::vector<int> NegMu_DiMuVeto;
	NegMu_DiMuVeto.clear();
	for(unsigned i_mu=0; i_mu<Ntp->NMuons(); i_mu++){
		if(selectMuon_DiMuonVeto(i_mu, selVertex)){
			if(Ntp->Muon_Charge(i_mu) == -1) NegMu_DiMuVeto.push_back(i_mu);
			else if(Ntp->Muon_Charge(i_mu) == 1) PosMu_DiMuVeto.push_back(i_mu);
		}
	}
	double dRMax(-1.);
	for(unsigned i_posmu = 0; i_posmu < PosMu_DiMuVeto.size(); i_posmu++){
		for(unsigned i_negmu = 0; i_negmu < NegMu_DiMuVeto.size(); i_negmu++){
			double dR_tmp = Ntp->Muon_p4(PosMu_DiMuVeto.at(i_posmu)).DeltaR(Ntp->Muon_p4(NegMu_DiMuVeto.at(i_negmu)));
			if(dR_tmp > dRMax) dRMax = dR_tmp;
		}
	}
	value.at(DiMuonVeto) = dRMax;
	pass.at(DiMuonVeto) = (value.at(DiMuonVeto)<cut.at(DiMuonVeto));
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at DiMuonVeto: " <<value.at(DiMuonVeto) << std::endl;
		Logger(Logger::Verbose) << "pass at DiMuonVeto: " <<pass.at(DiMuonVeto) << std::endl;
	}

	// Tau cuts
	if(selection_verbose) Logger(Logger::Verbose) << "Cut on TauID" << std::endl;
	std::vector<int> selectedTausId;
	selectedTausId.clear();
	for(unsigned i_tau=0; i_tau < Ntp->NPFTaus(); i_tau++){
		if (selectPFTau_Id(i_tau,selectedMuonsId)){
			selectedTausId.push_back(i_tau);
		}
	}
	value.at(NTauId)=selectedTausId.size();
	pass.at(NTauId)=(value.at(NTauId)>=cut.at(NTauId));
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at NTauId: " <<value.at(NTauId) << std::endl;
		Logger(Logger::Verbose) << "pass at NTauId: " <<pass.at(NTauId) << std::endl;
	}

	if(selection_verbose) Logger(Logger::Verbose) << "Cut on Tau Kinematics" << std::endl;
	std::vector<int> selectedTausKin;
	selectedTausKin.clear();
	for(std::vector<int>::iterator it_tau = selectedTausId.begin(); it_tau != selectedTausId.end(); ++it_tau){
		if ( selectPFTau_Kinematics(*it_tau) ){
			selectedTausKin.push_back(*it_tau);
		}
	}
	value.at(NTauKin)=selectedTausKin.size();
	pass.at(NTauKin)=(value.at(NTauKin)>=cut.at(NTauKin));
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at NTauKin: " <<value.at(NTauKin) << std::endl;
		Logger(Logger::Verbose) << "pass at NTauKin: " <<pass.at(NTauKin) << std::endl;
	}

	if(selection_verbose) Logger(Logger::Verbose) << "Cut on Tau Isolation" << std::endl;
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
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at NTauIso: " <<value.at(NTauIso) << std::endl;
		Logger(Logger::Verbose) << "pass at NTauIso: " <<pass.at(NTauIso) << std::endl;
	}

	// Charge of MuTau
	if(selection_verbose) Logger(Logger::Verbose) << "Cut on Charge of MuTau System" << std::endl;
	value.at(ChargeSum) = ChargeSumDummy;
	if(selTau != selTauDummy && selMuon != selMuonDummy){
		Charge = Ntp->Muon_Charge(selMuon) + Ntp->PFTau_Charge(selTau);
		value.at(ChargeSum) = Charge;
	}
	else{
		Charge = ChargeSumDummy;
	}
	pass.at(ChargeSum)=(value.at(ChargeSum)==cut.at(ChargeSum));
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at ChargeSum: " <<value.at(ChargeSum) << std::endl;
		Logger(Logger::Verbose) << "pass at ChargeSum: " <<pass.at(ChargeSum) << std::endl;
	}

	// Tau Decay Mode
	if(selection_verbose) Logger(Logger::Verbose) << "Cut on Tau Decay Mode" << std::endl;
	if(selTau != selTauDummy){
		value.at(TauDecayMode) = Ntp->PFTau_hpsDecayMode(selTau);
	}
	pass.at(TauDecayMode)= (value.at(TauDecayMode)>=cut.at(TauDecayMode));
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at TauDecayMode: " <<value.at(TauDecayMode) << std::endl;
		Logger(Logger::Verbose) << "pass at TauDecayMode: " <<pass.at(TauDecayMode) << std::endl;
	}

	// Tau FlightLength Significance
	if(selection_verbose) Logger(Logger::Verbose) << "Cut on Tau Flight Length Significance" << std::endl;
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
			double FLSigma = Ntp->PFTau_FlightLength_significance(Ntp->PFTau_TIP_primaryVertex_pos(selTau), Ntp->PFTau_TIP_primaryVertex_cov(selTau), Ntp->PFTau_TIP_secondaryVertex_pos(selTau), Ntp->PFTau_TIP_secondaryVertex_cov(selTau));
			//Logger(Logger::Debug) << "FLSigma Ich: " << Ntp->PFTau_FlightLength_significance(selTau) << std::endl;
			Logger(Logger::Debug) << "FLSigma Vladimir: " << FLSigma << std::endl;
			if(Ntp->PFTau_3PS_A1_LV(selTau).Vect().Dot(Ntp->PFTau_FlightLength3d(selTau)) < 0){
				value.at(TauFLSigma) = -FLSigma;
			}
			else{
				value.at(TauFLSigma) = FLSigma;
			}
		}
	}

	pass.at(TauFLSigma) = (value.at(TauFLSigma)>=cut.at(TauFLSigma));
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at TauFLSigma: " <<value.at(TauFLSigma) << std::endl;
		Logger(Logger::Verbose) << "pass at TauFLSigma: " <<pass.at(TauFLSigma) << std::endl;
	}

	// MT calculation
	if(selection_verbose) Logger(Logger::Verbose) << "Calculation and Cut on MT distribution" << std::endl;
	double pT,phi,eTmiss,eTmPhi;
	double MT_TauMET;

	if(selMuon == selMuonDummy){
		value.at(MT_MuMET) = MTDummy;
		if(selection_verbose) Logger(Logger::Verbose) << "No Muon selected: neither isolated nor anti isolated" << std::endl;
	}
	else if(selMuon_Iso != selMuonDummy && selMuon_AntiIso != selMuonDummy){
		value.at(MT_MuMET) = MTDummy;
		Logger(Logger::Error) << "CRITICAL: SELECTED MUON PASSED ISOLATION AND ANTI-ISOLATION CUT --> FIX" << std::endl;
	}
	else if(selMuon != selMuonDummy){
		eTmiss					= Ntp->MET_CorrMVAMuTau_et();
		eTmPhi					= Ntp->MET_CorrMVAMuTau_phi();
		pT						= Ntp->Muon_p4(selMuon).Pt();
		phi						= Ntp->Muon_p4(selMuon).Phi();
		value.at(MT_MuMET)		= Ntp->transverseMass(pT,phi,eTmiss,eTmPhi);
	}
	if(value.at(MT_MuMET) != MTDummy) pass.at(MT_MuMET)=(value.at(MT_MuMET)<cut.at(MT_MuMET));
	if(selection_verbose){
		Logger(Logger::Verbose) << "value at MT_MuMET: " <<value.at(MT_MuMET) << std::endl;
		Logger(Logger::Verbose) << "pass at MT_MuMET: " <<pass.at(MT_MuMET) << std::endl;
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
	if(selection_verbose) Logger(Logger::Verbose) << "Calculation of Mvis" << std::endl;
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
	if(selection_verbose) Logger(Logger::Verbose) << "W+Jets BG Method" << std::endl;
	std::vector<unsigned int> exclude_cuts;
	exclude_cuts.clear();
	exclude_cuts.push_back(ChargeSum);
	exclude_cuts.push_back(MT_MuMET);
	exclude_cuts.push_back(TauDecayMode);
	exclude_cuts.push_back(TauFLSigma);

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
	if(selection_verbose) Logger(Logger::Verbose) << "QCD ABCD BG Method" << std::endl;
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

	if(selection_verbose){
		Logger(Logger::Verbose) << "------------------------" << std::endl;
		if(status){
			Logger(Logger::Verbose) << "!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			Logger(Logger::Verbose) << "Event passed all cuts" << std::endl;
			Logger(Logger::Verbose) << "!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		}
		else Logger(Logger::Verbose) << "Event failed selection" << std::endl;
		Logger(Logger::Verbose) << "------------------------" << std::endl;
	}

	if(passAllBut(TauFLSigma) && value.at(TauFLSigma) != TauFLSigmaDummy){
		Tau_pt_wo_FLSigmaCut.at(t).Fill(Ntp->PFTau_p4(selTau).Pt(), w);
		Tau_phi_wo_FLSigmaCut.at(t).Fill(Ntp->PFTau_p4(selTau).Phi(), w);
		Tau_eta_wo_FLSigmaCut.at(t).Fill(Ntp->PFTau_p4(selTau).Eta(), w);
	}

	GEFObject FitResults, FitResultsNoCorr, FitResultswithFullRecoil;
	std::vector<GEFObject> FitResults_wRecoil_MassScan;
	TPTRObject TPResults;
	TVector2 Pt_Z;
	objects::MET MET(Ntp, "CorrMVAMuTau");
	TMatrixT<double> METpar(2,1); METpar(0,0) = MET.ex(); METpar(1,0) = MET.ey();
	TMatrixTSym<double> METCov; METCov.ResizeTo(2,2); METCov = MET.significanceMatrix< TMatrixTSym<double> >();
	PTObject MET2(METpar, METCov);
	TLorentzVector Recoil;
	std::vector<bool> fitstatuses, fitstatusesNoCorr, fitstatuses_wRecoil;
	TLorentzVector Reco_Z_Corr, Reco_Z_Corr_wRecoil;

	if(status && value.at(TauFLSigma) != TauFLSigmaDummy){

		MVAMET_metobject_XX.at(t).Fill(MET.significanceXX(), w);
		MVAMET_ptobject_XX.at(t).Fill(MET2.Cov()(0,0), w);

		LorentzVectorParticle A1 = Ntp->PFTau_a1_lvp(selTau);
		TrackParticle MuonTP = Ntp->Muon_TrackParticle(selMuon);
		TVector3 PV = Ntp->PFTau_TIP_primaryVertex_pos(selTau);
		TMatrixTSym<double> PVCov = Ntp->PFTau_TIP_primaryVertex_cov(selTau);

		for(unsigned i; i<Ntp->Vtx_nTrk(selVertex); i++){
			Recoil += Ntp->Vtx_TracksP4(selVertex, i);
		}
		Recoil -= Ntp->Muon_p4(selMuon);
		Recoil -= Ntp->PFTau_3PS_A1_LV(selTau);
		double Phi_Res = (Recoil.Phi() > 0) ? Recoil.Phi() - TMath::Pi() : Recoil.Phi() + TMath::Pi();
		Pt_Z = Recoil.Vect().XYvector();
		Pt_Z.Set(-Pt_Z.X(), -Pt_Z.Y()); //rotation by pi

		GlobalEventFit GEF(MuonTP, A1, Phi_Res, PV, PVCov);
		GEF.SetCorrectPt(false);
		//GEF.setMaxDelta(30.);
		//GEF.setMassConstraint(90.);
		TPResults = GEF.getTPTRObject();
		FitResults = GEF.Fit();
		fitstatuses = GEF.getFitStatuses();

		GlobalEventFit GEFNoCorr(MuonTP, A1, Phi_Res, PV, PVCov);
		GEFNoCorr.SetCorrectPt(false);
		//GEF.setMaxDelta(30.);
		//GEF.setMassConstraint(90.);
		FitResultsNoCorr = GEF.Fit();
		fitstatusesNoCorr = GEF.getFitStatuses();

		GlobalEventFit GEFwithFullRecoil(MuonTP, A1, MET2, PV, PVCov);
		GEFwithFullRecoil.SetCorrectPt(false);
		FitResultswithFullRecoil = GEFwithFullRecoil.Fit();
		fitstatuses_wRecoil = GEFwithFullRecoil.getFitStatuses();

		/*
		for(unsigned i_mass=0;i_mass<50;i_mass++){
		  GlobalEventFit GEFwithFullRecoilMassScan(MuonTP, A1, MET2, PV, PVCov);
		  GEFwithFullRecoilMassScan.setMassConstraint(28.5 + 2.5*i_mass);
		  GEFObject tmp = GEFwithFullRecoilMassScan.Fit();
		  FitResults_wRecoil_MassScan.push_back(tmp);
		}
		*/
		unsigned i_minchi2=999;
		double minchi2=999;
		bool MassScanConverged(false);
		/*
		for(unsigned i_mass=0;i_mass<50;i_mass++){
		  if(FitResults_wRecoil_MassScan.at(i_mass).Fitconverged()){
			if(FitResults_wRecoil_MassScan.at(i_mass).getChi2()<minchi2){
			  i_minchi2 = i_mass;
			  minchi2 = FitResults_wRecoil_MassScan.at(i_mass).getChi2();
			}
		  }
		}
		if(i_minchi2 != 999) MassScanConverged = true;
		*/

		Logger(Logger::Debug) << "Results.getTauH().LV().Pt(): " << FitResults.getTauH().LV().Pt() << std::endl;
		Logger(Logger::Debug) << "ResultswithFullRecoil.getTauH().LV().Pt(): " << FitResultswithFullRecoil.getTauH().LV().Pt() << std::endl;

		Logger(Logger::Debug) << "Results.getTauMu().LV().Pt(): " << FitResults.getTauMu().LV().Pt() << std::endl;
		Logger(Logger::Debug) << "ResultswithFullRecoil.getTauMu().LV().Pt(): " << FitResultswithFullRecoil.getTauMu().LV().Pt() << std::endl;

		if(FitResults.Fitconverged()){
			Reco_Z_Corr = FitResults.getTauH().LV() + FitResults.getTauMu().LV();

			Reco_EventFit_Solution.at(t).Fill(FitResults.getIndex(), w);
			Reco_EventFit_Solution.at(t).Fill(-1, w); //all solutions
			Reco_ConstrainedDeltaSum.at(t).Fill(FitResults.getCsum(), w);
			Reco_Chi2.at(t).Fill(FitResults.getChi2(), w);
			Reco_TauMu_DeltaPX_FitImpact.at(t).Fill(FitResults.getInitTauMu().LV().Px() - FitResults.getTauMu().LV().Px(), w);
			Reco_TauMu_DeltaPY_FitImpact.at(t).Fill(FitResults.getInitTauMu().LV().Py() - FitResults.getTauMu().LV().Py(), w);
			Reco_TauMu_DeltaPZ_FitImpact.at(t).Fill(FitResults.getInitTauMu().LV().Pz() - FitResults.getTauMu().LV().Pz(), w);
			Reco_TauA1_DeltaPX_FitImpact.at(t).Fill(FitResults.getInitTauH().LV().Px() - FitResults.getTauH().LV().Px(), w);
			Reco_TauA1_DeltaPY_FitImpact.at(t).Fill(FitResults.getInitTauH().LV().Py() - FitResults.getTauH().LV().Py(), w);
			Reco_TauA1_DeltaPZ_FitImpact.at(t).Fill(FitResults.getInitTauH().LV().Pz() - FitResults.getTauH().LV().Pz(), w);

			Reco_TauA1_P.at(t).Fill(FitResults.getTauH().LV().P(), w);
			Reco_TauA1_Pt.at(t).Fill(FitResults.getTauH().LV().Pt(), w);
			Reco_TauA1_Px.at(t).Fill(FitResults.getTauH().LV().Px(), w);
			Reco_TauA1_Py.at(t).Fill(FitResults.getTauH().LV().Py(), w);
			Reco_TauA1_Pz.at(t).Fill(FitResults.getTauH().LV().Pz(), w);
			Reco_TauA1_Phi.at(t).Fill(FitResults.getTauH().LV().Phi(), w);
			Reco_TauA1_Eta.at(t).Fill(FitResults.getTauH().LV().Eta(), w);

			Reco_TauMu_P.at(t).Fill(FitResults.getTauMu().LV().P(), w);
			Reco_TauMu_Pt.at(t).Fill(FitResults.getTauMu().LV().Pt(), w);
			Reco_TauMu_Px.at(t).Fill(FitResults.getTauMu().LV().Px(), w);
			Reco_TauMu_Py.at(t).Fill(FitResults.getTauMu().LV().Py(), w);
			Reco_TauMu_Pz.at(t).Fill(FitResults.getTauMu().LV().Pz(), w);
			Reco_TauMu_Phi.at(t).Fill(FitResults.getTauMu().LV().Phi(), w);
			Reco_TauMu_Eta.at(t).Fill(FitResults.getTauMu().LV().Eta(), w);

			Reco_Z_Px.at(t).Fill(Reco_Z_Corr.Px(), w);
			Reco_Z_Py.at(t).Fill(Reco_Z_Corr.Py(), w);
			Reco_Z_Pz.at(t).Fill(Reco_Z_Corr.Pz(), w);
			Reco_Z_Pt.at(t).Fill(Reco_Z_Corr.Pt(), w);
			Reco_Z_Phi.at(t).Fill(Reco_Z_Corr.Phi(), w);
			Reco_Z_Eta.at(t).Fill(Reco_Z_Corr.Eta(), w);

			Reco_TauMu_DeltaPhi_FitImpact.at(t).Fill(FitResults.getTauMu().LV().DeltaPhi(FitResults.getInitTauMu().LV()), w);
			Reco_TauMu_ResCosTheta.at(t).Fill(FitResults.getTauMu().LV().CosTheta() - FitResults.getInitTauMu().LV().CosTheta(), w);
			RecoZ_Pt.at(t).Fill(Reco_Z_Corr.Pt(), w);
			Reco_dPhi_TauMuTauA1_AfterFit.at(t).Fill(FitResults.getTauMu().LV().DeltaPhi(FitResults.getTauH().LV()), w);
			Reco_dPhi_TauMuTauA1_BeforeFit.at(t).Fill(FitResults.getInitTauMu().LV().DeltaPhi(FitResults.getInitTauH().LV()), w);
			Reco_NIter.at(t).Fill(FitResults.getNiterations(), w);
			Reco_Chi2.at(t).Fill(FitResults.getChi2Vector().Sum(), w);
			Reco_Chi2_orig.at(t).Fill(FitResults.getChi2Vector()(0), w);
			Reco_Chi2_SC.at(t).Fill(FitResults.getChi2Vector()(1), w);
			Reco_Chi2_HC.at(t).Fill(FitResults.getChi2Vector()(2), w);

			if(FitResults.getIndex() == 0){
				Reco_NIter_noAmb.at(t).Fill(FitResults.getNiterations(), w);
				Reco_Chi2_noAmb.at(t).Fill(FitResults.getChi2Vector().Sum(), w);
				Reco_Chi2_orig_noAmb.at(t).Fill(FitResults.getChi2Vector()(0), w);
				Reco_Chi2_SC_noAmb.at(t).Fill(FitResults.getChi2Vector()(1), w);
				Reco_Chi2_HC_noAmb.at(t).Fill(FitResults.getChi2Vector()(2), w);
			}
			else{
				Reco_NIter_wAmb.at(t).Fill(FitResults.getNiterations(), w);
				Reco_Chi2_wAmb.at(t).Fill(FitResults.getChi2Vector().Sum(), w);
				Reco_Chi2_orig_wAmb.at(t).Fill(FitResults.getChi2Vector()(0), w);
				Reco_Chi2_SC_wAmb.at(t).Fill(FitResults.getChi2Vector()(1), w);
				Reco_Chi2_HC_wAmb.at(t).Fill(FitResults.getChi2Vector()(2), w);
			}
		}
		else Reco_EventFit_Solution.at(t).Fill(-2., w);

		if(FitResultswithFullRecoil.Fitconverged()){
			Reco_dPhi_TauMuTauA1_AfterFit_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().DeltaPhi(FitResultswithFullRecoil.getTauH().LV()), w);
			Reco_dPhi_TauMuTauA1_BeforeFit_wRecoil.at(t).Fill(FitResultswithFullRecoil.getInitTauMu().LV().DeltaPhi(FitResultswithFullRecoil.getInitTauH().LV()), w);
			Reco_Chi2_FitSolutionOnly_wRecoil.at(t).Fill(FitResultswithFullRecoil.getChi2(), w);
			Reco_NIter_wRecoil.at(t).Fill(FitResultswithFullRecoil.getNiterations(), w);
			Reco_TauA1_P_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().P(), w);
			Reco_TauA1_Pt_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Pt(), w);
			Reco_TauA1_Px_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Px(), w);
			Reco_TauA1_Py_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Py(), w);
			Reco_TauA1_Pz_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Pz(), w);
			Reco_TauA1_Phi_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Phi(), w);
			Reco_TauA1_Eta_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Eta(), w);
			Reco_TauMu_P_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().P(), w);
			Reco_TauMu_Pt_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Pt(), w);
			Reco_TauMu_Px_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Px(), w);
			Reco_TauMu_Py_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Py(), w);
			Reco_TauMu_Pz_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Pz(), w);
			Reco_TauMu_Phi_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Phi(), w);
			Reco_TauMu_Eta_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Eta(), w);

			Reco_Z_Corr_wRecoil = FitResultswithFullRecoil.getTauMu().LV() + FitResultswithFullRecoil.getTauH().LV();
			Reco_Z_Px_wRecoil.at(t).Fill(Reco_Z_Corr_wRecoil.Px(), w);
			Reco_Z_Py_wRecoil.at(t).Fill(Reco_Z_Corr_wRecoil.Py(), w);
			Reco_Z_Pz_wRecoil.at(t).Fill(Reco_Z_Corr_wRecoil.Pz(), w);
			Reco_Z_Pt_wRecoil.at(t).Fill(Reco_Z_Corr_wRecoil.Pt(), w);
			Reco_Z_Phi_wRecoil.at(t).Fill(Reco_Z_Corr_wRecoil.Phi(), w);
			Reco_Z_Eta_wRecoil.at(t).Fill(Reco_Z_Corr_wRecoil.Eta(), w);

			Reco_EventFit_Solution_wRecoil.at(t).Fill(-1., w);
			Reco_EventFit_Solution_wRecoil.at(t).Fill(FitResultswithFullRecoil.getIndex(), w);

			Reco_Chi2_Full_wRecoil.at(t).Fill(FitResultswithFullRecoil.getChi2(), w);
			Reco_Chi2_Orig_wRecoil.at(t).Fill(FitResultswithFullRecoil.getChi2Vector()(0), w);
			Reco_Chi2_SC_wRecoil.at(t).Fill(FitResultswithFullRecoil.getChi2Vector()(1), w);
			Reco_Chi2_HC_wRecoil.at(t).Fill(FitResultswithFullRecoil.getChi2Vector()(2), w);
			Reco_Chi2_OrigProb_wRecoil.at(t).Fill(TMath::Prob(FitResultswithFullRecoil.getChi2Vector()(0), 1), w);

			if(fitstatuses_wRecoil.at(1) && fitstatuses_wRecoil.at(2)){
				Reco_Chi2_diff_wRecoil.at(t).Fill(FitResultswithFullRecoil.getChi2Vectors().at(2).Sum() - FitResultswithFullRecoil.getChi2Vectors().at(1).Sum(), w);
				Reco_Chi2_orig_diff_wRecoil.at(t).Fill(FitResultswithFullRecoil.getChi2Vectors().at(2)(0) - FitResultswithFullRecoil.getChi2Vectors().at(1)(0), w);
			}
		}
		else Reco_EventFit_Solution_wRecoil.at(t).Fill(-2., w);


		if(MassScanConverged){
			LorentzVectorParticle Resonance_MassScan = FitResults_wRecoil_MassScan.at(i_minchi2).getResonance();
			Reco_ZMass_MassScan.at(t).Fill(Resonance_MassScan.LV().M(), w);
		}

		A1_Par_Px.at(t).Fill(TPResults.getA1().LV().Px(), w);
		A1_Par_Py.at(t).Fill(TPResults.getA1().LV().Py(), w);
		A1_Par_Pz.at(t).Fill(TPResults.getA1().LV().Pz(), w);
		A1_Par_M.at(t).Fill(TPResults.getA1().LV().M(), w);
		A1_Cov_Pxx.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px), w);
		A1_Cov_Pyy.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py), w);
		A1_Cov_Pzz.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::pz,LorentzVectorParticle::pz), w);
		A1_Cov_Pxy.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py), w);
		A1_Cov_Pxz.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::pz), w);
		A1_Cov_Pyz.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::pz), w);

		SV_Par_x.at(t).Fill(TPResults.getA1().Vertex().X(), w);
		SV_Par_y.at(t).Fill(TPResults.getA1().Vertex().Y(), w);
		SV_Par_z.at(t).Fill(TPResults.getA1().Vertex().Z(), w);
		SV_Cov_xx.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::vx,LorentzVectorParticle::vx), w);
		SV_Cov_yy.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::vy,LorentzVectorParticle::vy), w);
		SV_Cov_zz.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::vz,LorentzVectorParticle::vz), w);
		SV_Cov_xy.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::vx,LorentzVectorParticle::vy), w);
		SV_Cov_xz.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::vx,LorentzVectorParticle::vz), w);
		SV_Cov_yz.at(t).Fill(TPResults.getA1().Covariance(LorentzVectorParticle::vy,LorentzVectorParticle::vz), w);

	}
	/*
	//DiTau Fit

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
		Logger(Logger::Verbose) << "" << std::endl;
		Logger(Logger::Verbose) << "Starting Constrained DiTau Fit \n" << std::endl;
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
				Logger(Logger::Debug) << "Chi2 for ambiguity " << Ambiguity << " : " << LC_chi2 << std::endl;
				bool EventFit_bool = Ntp->EventFit(selTau, selMuon, selVertex, TPTF_TausA1.at(Ambiguity), Reco_Z, tmp_Daughters, tmp_Daughters0, LC_chi2, Niterat, csum, 91.5);//, tmp_par_0, tmp_par);
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
		Logger(Logger::Verbose) << "Ambiguitysolver: " << std::endl;
		if(Ntp->AmbiguitySolverByChi2(A1Fit, EventFit, Chi2s, IndexToReturn, AmbiguityPoint)){
			int IndexToReturnTEST(-1); bool AmbiguityPointTEST(false);
			Ntp->AmbiguitySolver(A1Fit, EventFit, Probs, IndexToReturnTEST, AmbiguityPointTEST);
			AmbiguitySolvable = true;
			if(ZFits.at(IndexToReturn).Mass() >=0) Reco_ZMass.at(t).Fill(ZFits.at(IndexToReturn).LV().M(), w);
			Logger(Logger::Verbose) << "Single Mass Fit; Fit Mass:  "<< ZFits.at(IndexToReturn).LV().M() << std::endl;
			Logger(Logger::Verbose) << "Picked Solution with Chi2 for ambiguity " << IndexToReturn << " : " << Chi2s.at(IndexToReturn) << std::endl;
			Logger(Logger::Verbose) << "Picked Solution with Ambiguity value (by Chi2s): " << IndexToReturn << std::endl;
			Logger(Logger::Verbose) << "Picked Solution with Ambiguity value (by probs): " << IndexToReturnTEST << std::endl;
			Reco_EventFit_Solution.at(t).Fill(IndexToReturn, w);
			Reco_EventFit_Solution.at(t).Fill(-1, w); //all solutions
			Reco_ConstrainedDeltaSum.at(t).Fill(Csums.at(IndexToReturn), w);
			Reco_Chi2_FitSolutionOnly.at(t).Fill(Chi2s.at(IndexToReturn), w);
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
			Reco_dPhi_TauMuTauA1_AfterFit.at(t).Fill(RefitDaughters.at(IndexToReturn).at(1).LV().Phi() - RefitDaughters.at(IndexToReturn).at(0).LV().Phi(), w);
			Reco_dPhi_TauMuTauA1_BeforeFit.at(t).Fill(InitDaughters.at(IndexToReturn).at(1).LV().Phi() - InitDaughters.at(IndexToReturn).at(0).LV().Phi(), w);
		}
		else{
			Logger(Logger::Verbose) << "Failed" << std::endl;
			Reco_EventFit_Solution.at(t).Fill(-2, w); //not able to solve ambiguity/no solution
		}
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
					bool SingleEventFit = Ntp->EventFit(selTau, selMuon, selVertex, Reco_TauA1, Reco_Z, RefitDaughters, InitDaughters, LC_chi2, Niterat, csum, MassConstraint);//), tmp2, tmp2);
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
				Logger(Logger::Verbose) << "Multiple Masses Fit; Fit Mass:  " << ZFits_MassScan.at(FinalIndex).at(IndicesToReturn.at(FinalIndex)).Mass() << std::endl;
				Logger(Logger::Verbose) << "Picked Solution with Chi2 for ambiguity " << IndicesToReturn.at(FinalIndex) << " : " << Chi2s_MassScan.at(FinalIndex).at(IndicesToReturn.at(FinalIndex)) << std::endl;
				Logger(Logger::Verbose) << "Picked Solution with Ambiguity value: " << IndicesToReturn.at(FinalIndex) << std::endl;
			}
		}
		Logger(Logger::Verbose) << "-----------------End of Fit-----------------" << std::endl;
	}
	 */

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
		A1massRefit.at(t).Fill(Ntp->PFTau_3PS_A1_LV(selTau).M(), w);
		dA1mass_PFTau_Refit.at(t).Fill(Ntp->PFTau_p4(selTau).M() - Ntp->PFTau_3PS_A1_LV(selTau).M(), w);
		TauFL_WithTauFLSigmaCut.at(t).Fill(Ntp->PFTau_FlightLength(selTau) , w);
		//std::cout << "PV Coord: " << std::endl;
		  //  Ntp->PFTau_TIP_primaryVertex_pos(selTau).Print();
		//std::cout << "SV Coord: " << std::endl;
		  //  Ntp->PFTau_TIP_secondaryVertex_pos(selTau).Print();

		TLorentzVector PFTau_UnFitTracks;
		(Ntp->PFTau_NdaughterTracks(selTau) == 3) ? TransTrk_Failure_withSelection.at(t).Fill(0.,w) : TransTrk_Failure_withSelection.at(t).Fill(1.,w);
		if(Ntp->PFTau_NdaughterTracks(selTau) == 3){
			for(unsigned i = 0; i<Ntp->PFTau_NdaughterTracks(selTau); i++){
				//std::cout << "Ntp->PFTau_daughterTracks(selTau).size()" << Ntp->PFTau_daughterTracks(selTau).size() << std::endl;
				TrackParticle tmpTP = Ntp->PFTau_daughterTracks(selTau).at(i);
				TVector3 SV = Ntp->PFTau_TIP_secondaryVertex_pos(selTau);
				TLorentzVector tmpLV = (TrackTools::LorentzParticleAtPosition(tmpTP, SV)).LV();
				PFTau_UnFitTracks += tmpLV;
			}
		}
		Tau_Mass_Difference_PFTau_UnFitTracks_3PS.at(t).Fill(PFTau_UnFitTracks.M() - Ntp->PFTau_p4(selTau).M(), w);

		TLorentzVector PFTau_ReFitTracks;
		if(Ntp->PFTau_NdaughterTracks(selTau) == 3){
			for(unsigned i = 0; i<Ntp->PFTau_NdaughtersReFitTracks_p4(selTau); i++){
				//std::cout << "Ntp->PFTau_daughterTracks(selTau).size()" << Ntp->PFTau_daughterTracks(selTau).size() << std::endl;
				TLorentzVector tmpTLV = Ntp->PFTau_daughterReFitTracks_p4(selTau).at(i);
				PFTau_UnFitTracks += tmpTLV;
			}
		}
		Tau_Mass_Difference_PFTau_UnFitTracks_3PS.at(t).Fill(PFTau_ReFitTracks.M() - Ntp->PFTau_3PS_A1_LV(selTau).M(), w);

		TrackParticle MuonTP = Ntp->Muon_TrackParticle(selMuon);

		double dxy   =fabs(Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::dxy));
		double phi0  = Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::phi);
		double lam   = Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::lambda);
		double dz   = Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::dz);
		double signed_dxy = Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::dxy);

		Mu_TP_phi0.at(t).Fill(MuonTP.Parameter(TrackParticle::phi), w);
		Mu_TP_lambda.at(t).Fill(MuonTP.Parameter(TrackParticle::lambda), w);
		Mu_TP_dxy.at(t).Fill(MuonTP.Parameter(TrackParticle::dxy), w);
		Mu_TP_dz.at(t).Fill(MuonTP.Parameter(TrackParticle::dz), w);
		Mu_TP_kappa.at(t).Fill(MuonTP.Parameter(TrackParticle::kappa), w);

		double xPoca = -signed_dxy*sin(phi0);
		double yPoca = signed_dxy*cos(phi0);
		double xPocaVlad = dxy*sin(phi0);
		double yPocaVlad = dxy*cos(phi0);

		Mu_TP_Poca_xy.at(t).Fill(xPoca, yPoca);
		Mu_TP_Vertex_xy.at(t).Fill(Ntp->Vtx(selVertex).X(), Ntp->Vtx(selVertex).Y());
		Mu_TP_RefitVertex_xy.at(t).Fill(Ntp->PFTau_TIP_primaryVertex_pos(selTau).X(), Ntp->PFTau_TIP_primaryVertex_pos(selTau).Y());
		Mu_TP_BeamSpot_xy.at(t).Fill(Ntp->beamspot_par()(0), Ntp->beamspot_par()(1));
		Mu_TP_NTP_Poca_xy.at(t).Fill(Ntp->Muon_Poca(selMuon).X(), Ntp->Muon_Poca(selMuon).Y());

		if(xPoca >= 0){
			if(yPoca >= 0) Mu_TP_POCA_quadrant.at(t).Fill(1, w);
			else Mu_TP_POCA_quadrant.at(t).Fill(4, w);
		}
		else {
			if(yPoca >= 0) Mu_TP_POCA_quadrant.at(t).Fill(2, w);
			else Mu_TP_POCA_quadrant.at(t).Fill(3, w);
		}
		if(xPocaVlad >= 0){
			if(yPocaVlad >= 0) Mu_TP_POCA_quadrantVlad.at(t).Fill(1, w);
			else Mu_TP_POCA_quadrantVlad.at(t).Fill(4, w);
		}
		else {
			if(yPocaVlad >= 0) Mu_TP_POCA_quadrantVlad.at(t).Fill(2, w);
			else Mu_TP_POCA_quadrantVlad.at(t).Fill(3, w);
		}

		double pi = TMath::Pi();
		if(signed_dxy > 0){
		  if(phi0 > -pi && phi0 < -pi/2) Mu_TP_POCA_quadrantby_dxyphi0.at(t).Fill(4, w);
		  if(phi0 > -pi/2 && phi0 < 0) Mu_TP_POCA_quadrantby_dxyphi0.at(t).Fill(1, w);
		  if(phi0 > 0 && phi0 < pi/2) Mu_TP_POCA_quadrantby_dxyphi0.at(t).Fill(2, w);
		  if(phi0 > pi/2 && phi0 < pi) Mu_TP_POCA_quadrantby_dxyphi0.at(t).Fill(3, w);
		}
		else{
		  if(phi0 > -pi && phi0 < -pi/2) Mu_TP_POCA_quadrantby_dxyphi0.at(t).Fill(2, w);
		  if(phi0 > -pi/2 && phi0 < 0) Mu_TP_POCA_quadrantby_dxyphi0.at(t).Fill(3, w);
		  if(phi0 > 0 && phi0 < pi/2) Mu_TP_POCA_quadrantby_dxyphi0.at(t).Fill(4, w);
		  if(phi0 > pi/2 && phi0 < pi) Mu_TP_POCA_quadrantby_dxyphi0.at(t).Fill(1, w);
		}
		//Mu_TP_POCA_quadrantby_dxyphi0

		TVector3 MuonDir(cos(phi0)*cos(lam), sin(phi0)*cos(lam), sin(lam));
		TVector3 Muon0(-sin(phi0)*signed_dxy,cos(phi0)*signed_dxy,dz);

		TVector3 MuonPocaNTP = Ntp->Muon_Poca(selMuon);
		TVector3 Muon0PocaNTP = MuonPocaNTP - Muon0;
		TVector3 MinusMuon0PocaNTP = Muon0 - MuonPocaNTP;
		TVector3 Muon0PV = Muon0 - Ntp->Vtx(selVertex);
		/*
		Logger(Logger::Debug) << "------MUON POCA:------" << std::endl;
		MuonDir.Print();
		Muon0PocaNTP.Print();
		MinusMuon0PocaNTP.Print();
		*/
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
		TLorentzVector GenTaumu, GenTauh, GenZH, GenMu, GenA1, GenNuTau;
		bool hasA1 = false;
		bool hasMu = false;
		bool hasZH = false;
		int TauMu_charge(0);
		int Tauh_charge(0);
		int GenTauhIndex(-1);
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
					GenTauhIndex = i;
					if(Ntp->MCTau_pdgid(i) == PDGInfo::tau_minus) Tauh_charge = -1;
					else if(Ntp->MCTau_pdgid(i) == PDGInfo::tau_plus) Tauh_charge = 1;
					//std::cout << "--------" << std::endl;
					for(int j=0; j<Ntp->NMCTauDecayProducts(i); j++){
						//std::cout << "PDG ID of decay particle " << j << ": " << Ntp->MCTauandProd_pdgid(i,j) <<std::endl;
						if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::a_1_plus){
							GenA1 = Ntp->MCTauandProd_p4(i,j);
							hasA1 = true;
						}
						if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::nu_tau){
							GenNuTau = Ntp->MCTauandProd_p4(i,j);
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

				Gen_TauA1_P_noSel.at(t).Fill(GenTauh.P(), w);
				Gen_TauA1_Pt_noSel.at(t).Fill(GenTauh.Pt(), w);
				Gen_TauA1_Px_noSel.at(t).Fill(GenTauh.Px(), w);
				Gen_TauA1_Py_noSel.at(t).Fill(GenTauh.Py(), w);
				Gen_TauA1_Pz_noSel.at(t).Fill(GenTauh.Pz(), w);
				Gen_TauA1_Phi_noSel.at(t).Fill(GenTauh.Phi(), w);
				Gen_TauA1_Eta_noSel.at(t).Fill(GenTauh.Eta(), w);

				Gen_TauMu_P_noSel.at(t).Fill(GenTaumu.P(), w);
				Gen_TauMu_Pt_noSel.at(t).Fill(GenTaumu.Pt(), w);
				Gen_TauMu_Px_noSel.at(t).Fill(GenTaumu.Px(), w);
				Gen_TauMu_Py_noSel.at(t).Fill(GenTaumu.Py(), w);
				Gen_TauMu_Pz_noSel.at(t).Fill(GenTaumu.Pz(), w);
				Gen_TauMu_Phi_noSel.at(t).Fill(GenTaumu.Phi(), w);
				Gen_TauMu_Eta_noSel.at(t).Fill(GenTaumu.Eta(), w);

				Gen_TauMu_GJ.at(t).Fill(GenMu.Vect().Angle(GenTaumu.Vect()), w);
				Gen_TauA1_GJ.at(t).Fill(GenA1.Vect().Angle(GenTauh.Vect()), w);
				Phi_genTaumu.at(t).Fill(GenTaumu.Phi(), w);
				Theta_genTaumu.at(t).Fill(GenTaumu.Theta(), w);
				Phi_genTauh.at(t).Fill(GenTauh.Phi(), w);
				Theta_genTauh.at(t).Fill(GenTauh.Theta(), w);

				double dPhi_genTaus = GenTauh.DeltaPhi(GenTaumu);
				//if(dPhi_genTaus > TMath::Pi()) dPhi_genTaus = 2*TMath::Pi() - dPhi_genTaus;
				TLorentzVector DiTau = GenTauh + GenTaumu;
				Gen_Z_Pt_noSel.at(t).Fill(GenZH.Pt(), w);
				Gen_DiTau_dPhi.at(t).Fill(dPhi_genTaus, w);
				Gen_DiTau_Pt.at(t).Fill(DiTau.Pt(), w);
				Gen_Z_M.at(t).Fill(GenZH.M(), w);
				Gen_Z_Eta_noSel.at(t).Fill(GenZH.Eta(), w);
				Gen_Z_Phi_noSel.at(t).Fill(GenZH.Phi(), w);
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

				Gen_Z_Pt_vs_MET.at(t).Fill(GenZH.Pt(), Ntp->MET_CorrMVAMuTau_et());

				dPhi_GenTauMu_GenMu.at(t).Fill(GenTaumu.DeltaPhi(GenMu), w);
				dTheta_GenTauMu_GenMu.at(t).Fill(GenTaumu.Theta() - GenMu.Theta(), w);

				if(status && value.at(TauFLSigma) != TauFLSigmaDummy){

					int GenIndex, GenWrongIndex;
					if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) < 0){
						GenIndex = 2;
						GenWrongIndex = 1;
					}
					else if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) > 0){
						GenIndex = 1;
						GenWrongIndex = 2;
					}
					else if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) == 0){
						GenIndex = 0;
						GenWrongIndex = 0;
					}

					A1_Pull_Px.at(t).Fill((TPResults.getA1().LV().Px() - GenA1.Px())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px)), w);
					A1_Pull_Py.at(t).Fill((TPResults.getA1().LV().Py() - GenA1.Py())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py)), w);
					A1_Pull_Pz.at(t).Fill((TPResults.getA1().LV().Pz() - GenA1.Pz())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::pz,LorentzVectorParticle::pz)), w);
					A1_Pull_M.at(t).Fill((TPResults.getA1().LV().M() - GenA1.M())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::m,LorentzVectorParticle::m)), w);

					SV_Pull_Px.at(t).Fill((TPResults.getA1().Vertex().X() - Ntp->MCTauandProd_Vertex(GenTauhIndex,1).X())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::vx,LorentzVectorParticle::vx)), w);
					SV_Pull_Py.at(t).Fill((TPResults.getA1().Vertex().Y() - Ntp->MCTauandProd_Vertex(GenTauhIndex,1).Y())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::vy,LorentzVectorParticle::vy)), w);
					SV_Pull_Pz.at(t).Fill((TPResults.getA1().Vertex().Z() - Ntp->MCTauandProd_Vertex(GenTauhIndex,1).Z())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::vz,LorentzVectorParticle::vz)), w);


					if(!TPResults.isAmbiguous()){
						TauA1_Par_Px_AmbZero.at(t).Fill(TPResults.getTauZero().LV().Px(), w);
						TauA1_Par_Py_AmbZero.at(t).Fill(TPResults.getTauZero().LV().Py(), w);
						TauA1_Par_Pz_AmbZero.at(t).Fill(TPResults.getTauZero().LV().Pz(), w);
						TauA1_Cov_Pxx_AmbZero.at(t).Fill(TPResults.getTauZero().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px), w);
						TauA1_Cov_Pyy_AmbZero.at(t).Fill(TPResults.getTauZero().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py), w);
						TauA1_Cov_Pzz_AmbZero.at(t).Fill(TPResults.getTauZero().Covariance(LorentzVectorParticle::pz,LorentzVectorParticle::pz), w);
						TauA1_Cov_Pxy_AmbZero.at(t).Fill(TPResults.getTauZero().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py), w);
						TauA1_Cov_Pxz_AmbZero.at(t).Fill(TPResults.getTauZero().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::pz), w);
						TauA1_Cov_Pyz_AmbZero.at(t).Fill(TPResults.getTauZero().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::pz), w);
						TauA1_Pull_Px_AmbZero.at(t).Fill((TPResults.getTauZero().LV().Px() - GenTauh.Px())/sqrt(TPResults.getTauZero().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px)), w);
						TauA1_Pull_Py_AmbZero.at(t).Fill((TPResults.getTauZero().LV().Py() - GenTauh.Py())/sqrt(TPResults.getTauZero().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py)), w);
						TauA1_Pull_Pz_AmbZero.at(t).Fill((TPResults.getTauZero().LV().Pz() - GenTauh.Pz())/sqrt(TPResults.getTauZero().Covariance(LorentzVectorParticle::pz,LorentzVectorParticle::pz)), w);

						A1_Pull_Px_AmbZero.at(t).Fill((TPResults.getA1().LV().Px() - GenA1.Px())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px)), w);
						A1_Pull_Py_AmbZero.at(t).Fill((TPResults.getA1().LV().Py() - GenA1.Py())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py)), w);
						A1_Pull_Pz_AmbZero.at(t).Fill((TPResults.getA1().LV().Pz() - GenA1.Pz())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::pz,LorentzVectorParticle::pz)), w);
						A1_Pull_M_AmbZero.at(t).Fill((TPResults.getA1().LV().M() - GenA1.M())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::m,LorentzVectorParticle::m)), w);

						SV_Pull_Px_AmbZero.at(t).Fill((TPResults.getA1().Vertex().X() - Ntp->MCTauandProd_Vertex(GenTauhIndex,1).X())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::vx,LorentzVectorParticle::vx)), w);
						SV_Pull_Py_AmbZero.at(t).Fill((TPResults.getA1().Vertex().Y() - Ntp->MCTauandProd_Vertex(GenTauhIndex,1).Y())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::vy,LorentzVectorParticle::vy)), w);
						SV_Pull_Pz_AmbZero.at(t).Fill((TPResults.getA1().Vertex().Z() - Ntp->MCTauandProd_Vertex(GenTauhIndex,1).Z())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::vz,LorentzVectorParticle::vz)), w);
					}

					if(TPResults.isAmbiguous()){
						TauA1_Par_Px_CorrectAmb.at(t).Fill(TPResults.getTaus().at(GenIndex).LV().Px(), w);
						TauA1_Par_Py_CorrectAmb.at(t).Fill(TPResults.getTaus().at(GenIndex).LV().Py(), w);
						TauA1_Par_Pz_CorrectAmb.at(t).Fill(TPResults.getTaus().at(GenIndex).LV().Pz(), w);
						TauA1_Cov_Pxx_CorrectAmb.at(t).Fill(TPResults.getTaus().at(GenIndex).Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px), w);
						TauA1_Cov_Pyy_CorrectAmb.at(t).Fill(TPResults.getTaus().at(GenIndex).Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py), w);
						TauA1_Cov_Pzz_CorrectAmb.at(t).Fill(TPResults.getTaus().at(GenIndex).Covariance(LorentzVectorParticle::pz,LorentzVectorParticle::pz), w);
						TauA1_Cov_Pxy_CorrectAmb.at(t).Fill(TPResults.getTaus().at(GenIndex).Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py), w);
						TauA1_Cov_Pxz_CorrectAmb.at(t).Fill(TPResults.getTaus().at(GenIndex).Covariance(LorentzVectorParticle::px,LorentzVectorParticle::pz), w);
						TauA1_Cov_Pyz_CorrectAmb.at(t).Fill(TPResults.getTaus().at(GenIndex).Covariance(LorentzVectorParticle::py,LorentzVectorParticle::pz), w);
						TauA1_Pull_Px_CorrectAmb.at(t).Fill((TPResults.getTaus().at(GenIndex).LV().Px()- GenTauh.Px())/sqrt(TPResults.getTaus().at(GenIndex).Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px)), w);
						TauA1_Pull_Py_CorrectAmb.at(t).Fill((TPResults.getTaus().at(GenIndex).LV().Py()- GenTauh.Py())/sqrt(TPResults.getTaus().at(GenIndex).Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py)), w);
						TauA1_Pull_Pz_CorrectAmb.at(t).Fill((TPResults.getTaus().at(GenIndex).LV().Pz()- GenTauh.Pz())/sqrt(TPResults.getTaus().at(GenIndex).Covariance(LorentzVectorParticle::pz,LorentzVectorParticle::pz)), w);

						TauA1_Par_Px_WrongAmb.at(t).Fill(TPResults.getTaus().at(GenWrongIndex).LV().Px(), w);
						TauA1_Par_Py_WrongAmb.at(t).Fill(TPResults.getTaus().at(GenWrongIndex).LV().Py(), w);
						TauA1_Par_Pz_WrongAmb.at(t).Fill(TPResults.getTaus().at(GenWrongIndex).LV().Pz(), w);
						TauA1_Cov_Pxx_WrongAmb.at(t).Fill(TPResults.getTaus().at(GenWrongIndex).Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px), w);
						TauA1_Cov_Pyy_WrongAmb.at(t).Fill(TPResults.getTaus().at(GenWrongIndex).Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py), w);
						TauA1_Cov_Pzz_WrongAmb.at(t).Fill(TPResults.getTaus().at(GenWrongIndex).Covariance(LorentzVectorParticle::pz,LorentzVectorParticle::pz), w);
						TauA1_Cov_Pxy_WrongAmb.at(t).Fill(TPResults.getTaus().at(GenWrongIndex).Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py), w);
						TauA1_Cov_Pxz_WrongAmb.at(t).Fill(TPResults.getTaus().at(GenWrongIndex).Covariance(LorentzVectorParticle::px,LorentzVectorParticle::pz), w);
						TauA1_Cov_Pyz_WrongAmb.at(t).Fill(TPResults.getTaus().at(GenWrongIndex).Covariance(LorentzVectorParticle::py,LorentzVectorParticle::pz), w);
						TauA1_Pull_Px_WrongAmb.at(t).Fill((TPResults.getTaus().at(GenWrongIndex).LV().Px()- GenTauh.Px())/sqrt(TPResults.getTaus().at(GenWrongIndex).Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px)), w);
						TauA1_Pull_Py_WrongAmb.at(t).Fill((TPResults.getTaus().at(GenWrongIndex).LV().Py()- GenTauh.Py())/sqrt(TPResults.getTaus().at(GenWrongIndex).Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py)), w);
						TauA1_Pull_Pz_WrongAmb.at(t).Fill((TPResults.getTaus().at(GenWrongIndex).LV().Pz()- GenTauh.Pz())/sqrt(TPResults.getTaus().at(GenWrongIndex).Covariance(LorentzVectorParticle::pz,LorentzVectorParticle::pz)), w);

						A1_Pull_Px_wAmb.at(t).Fill((TPResults.getA1().LV().Px() - GenA1.Px())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px)), w);
						A1_Pull_Py_wAmb.at(t).Fill((TPResults.getA1().LV().Py() - GenA1.Py())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py)), w);
						A1_Pull_Pz_wAmb.at(t).Fill((TPResults.getA1().LV().Pz() - GenA1.Pz())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::pz,LorentzVectorParticle::pz)), w);
						A1_Pull_M_wAmb.at(t).Fill((TPResults.getA1().LV().M() - GenA1.M())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::m,LorentzVectorParticle::m)), w);

						SV_Pull_Px_wAmb.at(t).Fill((TPResults.getA1().Vertex().X() - Ntp->MCTauandProd_Vertex(GenTauhIndex,1).X())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::vx,LorentzVectorParticle::vx)), w);
						SV_Pull_Py_wAmb.at(t).Fill((TPResults.getA1().Vertex().Y() - Ntp->MCTauandProd_Vertex(GenTauhIndex,1).Y())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::vy,LorentzVectorParticle::vy)), w);
						SV_Pull_Pz_wAmb.at(t).Fill((TPResults.getA1().Vertex().Z() - Ntp->MCTauandProd_Vertex(GenTauhIndex,1).Z())/sqrt(TPResults.getA1().Covariance(LorentzVectorParticle::vz,LorentzVectorParticle::vz)), w);
					}

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
					double dPhi = POCAPV_dir.DeltaPhi(GenTaumu.Vect());
					dPhi_POCAPV_genTaumu.at(t).Fill(dPhi, w);
					double dTheta = POCAPV_dir.Theta() - GenTaumu.Theta();
					dTheta_POCAPV_genTaumu.at(t).Fill(dTheta, w);

					Phi_SVPV.at(t).Fill(Ntp->PFTau_FlightLength3d(selTau).Phi(), w);
					Theta_SVPV.at(t).Fill(Ntp->PFTau_FlightLength3d(selTau).Theta(), w);
					dPhi = Ntp->PFTau_FlightLength3d(selTau).DeltaPhi(GenTauh.Vect());
					dPhi_SVPV_genTauh.at(t).Fill(dPhi, w);
					dPhi_SVPV_genTauh_vs_TauFL.at(t).Fill(dPhi, Ntp->PFTau_FlightLength(selTau));
					if(Tauh_charge == 1) dPhi_SVPV_genTauhPlus_vs_TauFL.at(t).Fill(dPhi, Ntp->PFTau_FlightLength(selTau));
					else if(Tauh_charge == -1) dPhi_SVPV_genTauhMinus_vs_TauFL.at(t).Fill(dPhi, Ntp->PFTau_FlightLength(selTau));
					dTheta = Ntp->PFTau_FlightLength3d(selTau).Theta() - GenTauh.Theta();
					dTheta_SVPV_genTauh.at(t).Fill(dTheta, w);
					Angle_SVPV_genTauh.at(t).Fill(Ntp->PFTau_FlightLength3d(selTau).Angle(GenTauh.Vect()));

					TVector3 MinusPFTau_FlightLength3d = - Ntp->PFTau_FlightLength3d(selTau);
					dPhi = MinusPFTau_FlightLength3d.DeltaPhi(GenTaumu.Vect());
					dPhi_MinusSVPV_genTaumu.at(t).Fill(dPhi, w);
					dTheta = MinusPFTau_FlightLength3d.Theta() - GenTaumu.Theta();
					dTheta_MinusSVPV_genTaumu.at(t).Fill(dTheta, w);
					Angle_MinusSVPV_genTaumu.at(t).Fill(MinusPFTau_FlightLength3d.Angle(GenTaumu.Vect()));

					double A1_dPhi = Ntp->PFTau_3PS_A1_LV(selTau).DeltaPhi(GenA1);
					double A1_dTheta = Ntp->PFTau_3PS_A1_LV(selTau).Theta() - GenA1.Theta();
					A1_Phi_Res.at(t).Fill(A1_dPhi , w);
					A1_Theta_Res.at(t).Fill(A1_dTheta , w);

					/*
					std::vector< bool > TauA1Reco_StraightTau;
					std::vector< bool > TauA1Reco_HelixTau;
					for(unsigned Ambiguity=0; Ambiguity<3; Ambiguity++){
						double tmp_LC_chi2(-1), tmp_phisign(0);
						LorentzVectorParticle Reco_TauA1, Reco_Z;
						std::vector<LorentzVectorParticle> tmp_daughter;
						TVectorD tmp_par(3), tmp_par_0(3);
						TauA1Reco_StraightTau.push_back(Ntp->ThreeProngTauFit(selTau, Ambiguity, Reco_TauA1, tmp_daughter, tmp_LC_chi2, tmp_phisign));
						if(TauA1Reco_StraightTau.at(Ambiguity)){
							TLorentzVector Reco_TauA1HelixAtSV = TauHelixP4AtSV(selTau, Reco_TauA1.LV());
							TVector3 SVPV = TVector3(Ntp->PFTau_FlightLength3d(selTau));
							double GJAngleStraight = Ntp->PFTau_3PS_A1_LV(selTau).Angle(SVPV);
							double GJAngleHelix = Ntp->PFTau_3PS_A1_LV(selTau).Angle(Reco_TauA1HelixAtSV.Vect());
							double MaxGJAngle = GJAngleMax(Ntp->PFTau_3PS_A1_LV(selTau));
							double GJAngleStraightOverGJAngleMax = GJAngleStraight/MaxGJAngle;
							double GJAngleHelixOverGJAngleMax = GJAngleHelix/MaxGJAngle;

							Logger(Logger::Debug) << "GJAngleStraight " << GJAngleStraight << std::endl;
							Logger(Logger::Debug) << "GJAngleHelix " << GJAngleHelix << std::endl;
							Logger(Logger::Debug) << "GJAngleMax " << MaxGJAngle << std::endl;
							Logger(Logger::Debug) << "TauCharge " << Ntp->PFTau_Charge(selTau) << std::endl;

							if(GenTauhIndex != -1){
								Logger(Logger::Debug) << "GenTau1 Vertex :" << Ntp->MCTauandProd_Vertex(0,0).X() << ", " << Ntp->MCTauandProd_Vertex(0,0).Y() << ", " << Ntp->MCTauandProd_Vertex(0,0).Z() << std::endl;
								Logger(Logger::Debug) << "GenTau2 Vertex :" << Ntp->MCTauandProd_Vertex(1,0).X() << ", " << Ntp->MCTauandProd_Vertex(1,0).Y() << ", " << Ntp->MCTauandProd_Vertex(1,0).Z() << std::endl;
								Logger(Logger::Debug) << "GenTauh Sec Vertex :" << Ntp->MCTauandProd_Vertex(GenTauhIndex,1).X() << ", " << Ntp->MCTauandProd_Vertex(GenTauhIndex,1).Y() << ", " << Ntp->MCTauandProd_Vertex(GenTauhIndex,1).Z() << std::endl;
								Logger(Logger::Debug) << "GenTauh Sec Vertex :" << Ntp->MCTauandProd_Vertex(GenTauhIndex,2).X() << ", " << Ntp->MCTauandProd_Vertex(GenTauhIndex,2).Y() << ", " << Ntp->MCTauandProd_Vertex(GenTauhIndex,2).Z() << std::endl;
								Logger(Logger::Debug) << "GenTauh Sec Vertex :" << Ntp->MCTauandProd_Vertex(GenTauhIndex,3).X() << ", " << Ntp->MCTauandProd_Vertex(GenTauhIndex,3).Y() << ", " << Ntp->MCTauandProd_Vertex(GenTauhIndex,3).Z() << std::endl;
							}

							double dGJAngleStraight = GJAngleStraight - MaxGJAngle;
							double dGJAngleHelix = GJAngleHelix - MaxGJAngle;
							dGJAngle_GJAngleMAX_StraightTau.at(t).Fill(dGJAngleStraight, w);
							dGJAngle_GJAngleMAX_HelixTau.at(t).Fill(dGJAngleHelix, w);
							Angle_HelixTau_StraightTau.at(t).Fill(Reco_TauA1.LV().Angle(Reco_TauA1HelixAtSV.Vect()), w);
							dGJAngle_HelixTau_StraightTau.at(t).Fill(GJAngleHelix - GJAngleStraight, w);
							dGJAngle_HelixTau_StraightTauOverGJAngle.at(t).Fill((GJAngleHelix - GJAngleStraight)/GJAngleStraight, w);
							GJAngle_Over_GJAngleMax_StraightTau.at(t).Fill(GJAngleStraightOverGJAngleMax, w);
							GJAngle_Over_GJAngleMax_HelixTau.at(t).Fill(GJAngleHelixOverGJAngleMax, w);
							if(Ambiguity != 0) NUnphysical_StraightTau_HelixTau.at(t).Fill(Ambiguity, w);
							if(GJAngleStraight > MaxGJAngle){
								NUnphysical_StraightTau_HelixTau.at(t).Fill(3., w);
							}
							if(GJAngleHelix > MaxGJAngle){
								NUnphysical_StraightTau_HelixTau.at(t).Fill(4., w);
							}
						}
					}
					for(unsigned int i=0;i<TauA1Reco_StraightTau.size();i++){
						if(TauA1Reco_StraightTau.at(i)){
							TauA1_Reco_Solution_StraightTau.at(t).Fill(-1, w);
							NUnphysical_StraightTau_HelixTau.at(t).Fill(-1., w);
							break;
						}
						if(i==2){
							TauA1_Reco_Solution_StraightTau.at(t).Fill(-2, w);
							NUnphysical_StraightTau_HelixTau.at(t).Fill(-2., w);
						}
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

					Gen_Z_Pt_vs_VtxTracksPt.at(t).Fill(GenZH.Pt(),Pt_Z.Mod());
					Gen_Z_Phi_vs_VtxTracksPhi.at(t).Fill(GenZH.Phi(),Pt_Z.Phi());
					VtxTracksPtRes.at(t).Fill(Pt_Z.Mod() - GenZH.Pt(), w);
					//double phinew = (Recoil.Phi() > 0) ? Recoil.Phi() - TMath::Pi() : Recoil.Phi() + TMath::Pi();
					VtxTracksPhiCorrectedRes.at(t).Fill(Pt_Z.Phi_mpi_pi(Pt_Z.Phi()) - GenZH.Phi(), w);
					dPhi_GenTauMu_RecoMu.at(t).Fill(Ntp->Muon_p4(selMuon).DeltaPhi(GenTaumu));
					dTheta_GenTauMu_RecoMu.at(t).Fill(GenTaumu.Theta() - Ntp->Muon_p4(selMuon).Theta());

					TLorentzVector GenTauhHelix = GenTauHelixP4AtSV(GenTauhIndex, GenTauh);
					TLorentzVector GenNeutrinoHelixTau = GenTauhHelix - GenA1;
					TLorentzVector GenNeutrinoStraightTau = GenTauh - GenA1;
					TLorentzVector Neutrino_RefitPFTau_GenTauHelix = GenTauhHelix - Ntp->PFTau_3PS_A1_LV(selTau);

					TPTF_Neutrino_RefitPFTau_HelixGenTau_Mass.at(t).Fill(Neutrino_RefitPFTau_GenTauHelix.M(), w);
					TPTF_Neutrino_GenA1_HelixGenTau_Mass.at(t).Fill(GenNeutrinoHelixTau.M(), w);
					TPTF_Neutrino_GenA1_StraightGenTau_Mass.at(t).Fill(GenNeutrinoStraightTau.M(), w);
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


					double dxy   =fabs(Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::dxy));
					double phi0  = Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::phi);
					double lam   = Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::lambda);
					double dz   = Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::dz);
					double signed_dxy = Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::dxy);
					//  double phiAnot = inpar(4,0);

					LorentzVectorParticle A1 = Ntp->PFTau_a1_lvp(selTau);
					TVector3 MuonDir(cos(phi0)*cos(lam), sin(phi0)*cos(lam), sin(lam));
					TVector3 Muon0(-sin(phi0)*signed_dxy,cos(phi0)*signed_dxy,dz);
					TVector3 PV = Ntp->PFTau_TIP_primaryVertex_pos(selTau);
					TVector3 Muon0PV = Muon0 - Ntp->Vtx(selVertex);

					//Logger(Logger::Debug) << "Lambda: " << lam << std::endl;
					//Logger(Logger::Debug) << "phi0: " << phi0 << std::endl;
					//MuonDir.Print();
					TVector3 PFMuonDir = Ntp->Muon_p4(selMuon).Vect();
					PFMuonDir.SetMag(1.);
					//PFMuonDir.Print();
					//Muon0.Print();
					//Muon0PV.Print();
					//Ntp->PFTau_Poca(selTau).Print();
					//Ntp->Muon_Poca(selMuon).Print();
					//Ntp->Vtx(selVertex).Print();
					/*
					std::cout << "--------------" << std::endl;
					for(unsigned i_vtx = 0; i_vtx<Ntp->NVtx(); i_vtx++){
					  Ntp->Vtx(i_vtx).Print();
					}
					std::cout << "--------------" << std::endl;
					*/
					//Ntp->PFTau_TIP_primaryVertex_pos(selTau).Print();


					double dxy_ntp = -Ntp->Muon_p4(selMuon).Py()/Ntp->Muon_p4(selMuon).Pt()*Ntp->Muon_Poca(selMuon).X() + Ntp->Muon_p4(selMuon).Px()/Ntp->Muon_p4(selMuon).Pt()*Ntp->Muon_Poca(selMuon).Y();

					double alphaTP = Ntp->Muon_TrackParticle(selMuon).BField();

					//Logger(Logger::Debug) << "dxy_ntp: " << dxy_ntp << std::endl;
					//Logger(Logger::Debug) << "dxy trackparticle: " << signed_dxy << std::endl;
					//Logger(Logger::Debug) << "phi0_ntp with x: " << acos(Ntp->Muon_p4(selMuon).Px()/Ntp->Muon_p4(selMuon).Pt())*180/TMath::Pi() << std::endl;
					//Logger(Logger::Debug) << "phi0_ntp with y: " << asin(Ntp->Muon_p4(selMuon).Py()/Ntp->Muon_p4(selMuon).Pt())*180/TMath::Pi() << std::endl;
					//Logger(Logger::Debug) << "phi0 trackparticle: " << phi0*180/TMath::Pi() << std::endl;
					//Logger(Logger::Debug) << "kappa: " << Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::kappa) << "; 1/kappa: " << 1/Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::kappa) << std::endl;
					Logger(Logger::Debug) << "Muon Pt: " << Ntp->Muon_p4(selMuon).Pt() << "; Muon P: " << Ntp->Muon_p4(selMuon).P() << std::endl;
					Logger(Logger::Debug) << "alpha/kappa: " << alphaTP/Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::kappa) << std::endl;
					Logger(Logger::Debug) << "muonTP x, y: " << alphaTP/Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::kappa)*cos(Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::phi)) << ", " << alphaTP/Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::kappa)*sin(Ntp->Muon_TrackParticle(selMuon).Parameter(TrackParticle::phi)) << std::endl;
					Logger(Logger::Debug) << "Muon x, y: " << Ntp->Muon_p4(selMuon).X() << ", " << Ntp->Muon_p4(selMuon).Y() << std::endl;

					//Logger(Logger::Debug) << "alpha: " << alphaTP << std::endl;

					TVector3 SV_Ntp = Ntp->PFTau_TIP_secondaryVertex_pos(selTau);
					TVector3 SV = TPResults.getA1().Vertex();

					TVector3 TauDir = SV - PV;
					double TauA1OppositePhi = atan2(-TauDir.Y(), -TauDir.X());


					// reconstruct line of 1st tau
					double aPVSV = (PV.Y() - SV.Y())/(PV.X() - SV.X());
					double bPVSV = (PV.X()*SV.Y() - SV.X()*PV.Y())/(PV.X() - SV.X());

					double tmu  = (dxy*(cos(phi0) - aPVSV*sin(phi0)) - bPVSV)  / (aPVSV*cos(phi0) - sin(phi0)); //projection onto XY plane

					double xdoc = dxy*sin(phi0) + tmu*cos(phi0);
					double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
					double zdoc = dz + tmu*tan(lam);

					TVector3 PointGuess(xdoc, ydoc, zdoc);
					TVector3 TauMuDir = PV - PointGuess;

					double tNEW2  = (signed_dxy*(cos(phi0) + aPVSV*sin(phi0)) - bPVSV)  / (aPVSV*cos(phi0) - sin(phi0)); //projection onto XY plane

					double xdocNEW2 = -signed_dxy*sin(phi0) + tNEW2*cos(phi0);
					double ydocNEW2 = signed_dxy*cos(phi0) + tNEW2*sin(phi0);
					double zdocNEW2 = dz + tNEW2*tan(lam);

					TVector3 PointGuessNEW2(xdocNEW2, ydocNEW2, zdocNEW2);
					TVector3 TauMuDirNEW2 = PointGuessNEW2 - PV;

					LorentzVectorParticle RecoTauh_wTruth;
					RecoTauh_wTruth = (TPResults.isAmbiguous()) ? TPResults.getTaus().at(GenIndex) : TPResults.getTauZero();
					LorentzVectorParticle RecoNeutrino_wTruth;
					RecoNeutrino_wTruth = (TPResults.isAmbiguous()) ? TPResults.getNeutrinos().at(GenIndex) : TPResults.getNeutrinoZero();

					TrackParticle Muon = Ntp->Muon_TrackParticle(selMuon);

					TLorentzVector MET(Ntp->MET_CorrMVAMuTau_ex(), Ntp->MET_CorrMVAMuTau_ey(), 0, Ntp->MET_CorrMVAMuTau_et());
					TLorentzVector TauMu_wMET = MET - RecoNeutrino_wTruth.LV() + Ntp->Muon_p4(selMuon);
					TVector2 TauMu_wMETPtVec = TauMu_wMET.Vect().XYvector();
					TVector3 IS_wMET;
					TLorentzVector TauMu_wMET_afterMC = TauMuFullEstimate(PV, Muon, RecoTauh_wTruth, TauMu_wMETPtVec, IS_wMET);

					Est_TauMu_wMET_PtRes.at(t).Fill(TauMu_wMET.Pt() - GenTaumu.Pt());
					Est_TauMu_wMET_PhiRes.at(t).Fill(TauMu_wMET.DeltaPhi(GenTaumu));

					TLorentzVector Z_wMET = RecoTauh_wTruth.LV() + TauMu_wMET;
					Est_Z_wMET_PtRes.at(t).Fill(Z_wMET.Pt() - GenZH.Pt());
					Est_Z_wMET_PhiRes.at(t).Fill(Z_wMET.DeltaPhi(GenZH));

					TLorentzVector TauMuPtfromEventRecoil = Recoil - RecoTauh_wTruth.LV();
					TVector2 TauMu_EventRecoil_PtVec = TauMuPtfromEventRecoil.Vect().XYvector();
					TVector3 IS_EventRecoil;
					TLorentzVector TauMu_EventRecoil_afterMC = TauMuFullEstimate(PV, Muon, RecoTauh_wTruth, TauMu_EventRecoil_PtVec, IS_EventRecoil);

					TVector2 TauMu_PtBalance_PtVec(-RecoTauh_wTruth.LV().X(), -RecoTauh_wTruth.LV().Y());
					TVector3 IS_PtBalance;
					TLorentzVector TauMu_PtBalance_afterMC = TauMuFullEstimate(PV, Muon, RecoTauh_wTruth, TauMu_PtBalance_PtVec, IS_PtBalance);

					//TauMu_Start_Collinear_PtRes.at(t).Fill
					TauMu_Start_MET_PtRes.at(t).Fill(TauMu_wMET.Pt() - GenTaumu.Pt());
					TauMu_Start_PtBalance_PtRes.at(t).Fill(RecoTauh_wTruth.LV().Pt() - GenTaumu.Pt());
					TauMu_Start_EventRecoil_PtRes.at(t).Fill(TauMuPtfromEventRecoil.Pt() - GenTaumu.Pt());
					//TauMu_Start_Collinear_PtReco_vs_PtGen.at(t).Fill
					TauMu_Start_MET_PtReco_vs_PtGen.at(t).Fill(TauMu_wMET.Pt(), GenTaumu.Pt());
					TauMu_Start_PtBalance_PtReco_vs_PtGen.at(t).Fill(RecoTauh_wTruth.LV().Pt(), GenTaumu.Pt());
					TauMu_Start_EventRecoil_PtReco_vs_PtGen.at(t).Fill(TauMuPtfromEventRecoil.Pt(), GenTaumu.Pt());

					TauMu_Start_MET_PtRes_AfterMC.at(t).Fill(TauMu_wMET_afterMC.Pt() - GenTaumu.Pt());
					TauMu_Start_PtBalance_PtRes_AfterMC.at(t).Fill(TauMu_PtBalance_afterMC.Pt() - GenTaumu.Pt());
					TauMu_Start_EventRecoil_PtRes_AfterMC.at(t).Fill(TauMu_EventRecoil_afterMC.Pt() - GenTaumu.Pt());

					Z_Start_MET_PtRes.at(t).Fill(Z_wMET.Pt() - GenZH.Pt());
					Z_Start_MET_PhiRes.at(t).Fill(Z_wMET.DeltaPhi(GenZH));
					Z_Start_PtBalance_PtRes.at(t).Fill(- GenZH.Pt());
					//Z_Start_PtBalance_PhiRes.at(t).Fill(
					Z_Start_EventRecoil_PtRes.at(t).Fill(Pt_Z.Mod() - GenZH.Pt());
					Z_Start_EventRecoil_PhiRes.at(t).Fill(TVector2::Phi_mpi_pi(Pt_Z.Phi_mpi_pi(Pt_Z.Phi()) - GenZH.Phi()));


					TVector2 METVec(MET2.X(), MET2.Y());
					TLorentzVector Z_noNu_noRot = ZEstimatorWithMET(false, PV, Muon, A1, RecoTauh_wTruth, METVec);
					TLorentzVector Z_noNu_wRot = ZEstimatorWithMET(true, PV, Muon, A1, RecoTauh_wTruth, METVec);

					Z_Start_MET_noNu_PtRes.at(t).Fill(Z_noNu_noRot.Pt() - GenZH.Pt());
					Z_Start_MET_noNu_PhiRes.at(t).Fill(Z_noNu_noRot.DeltaPhi(GenZH));
					Z_Start_MET_noNu_M.at(t).Fill(Z_noNu_noRot.M());
					Z_Start_MET_noNu_dE.at(t).Fill(Z_noNu_noRot.E() - GenZH.E());

					Z_Start_MET_noNu_PtRes_withRot.at(t).Fill(Z_noNu_wRot.Pt() - GenZH.Pt());
					Z_Start_MET_noNu_PhiRes_withRot.at(t).Fill(Z_noNu_wRot.DeltaPhi(GenZH));
					Z_Start_MET_noNu_M_withRot.at(t).Fill(Z_noNu_wRot.M());

					TLorentzVector Z_B2B = ZEstimatorB2B(false, PV, Muon, A1, RecoTauh_wTruth);
					Z_Start_B2B_M.at(t).Fill(Z_B2B.M());
					Z_Start_B2B_dE.at(t).Fill(Z_B2B.E() - GenZH.E());


					if(FitResults.Fitconverged()){
						TauMu_Start_dPhi_TauMuTauH_MET.at(t).Fill(TauMu_wMET_afterMC.DeltaPhi(RecoTauh_wTruth.LV()));
						TauMu_Start_dPhi_TauMuTauH_EventRecoil.at(t).Fill(TauMu_EventRecoil_afterMC.DeltaPhi(RecoTauh_wTruth.LV()));
						TauMu_Start_dPhi_TauMuTauH_PtBalance.at(t).Fill(TauMu_PtBalance_afterMC.DeltaPhi(RecoTauh_wTruth.LV()));
					}

					TVector3 PV_RecoTau(RecoTauh_wTruth.Vertex());
					Logger(Logger::Debug) << "PV: " << PV.X() << ", " << PV.Y() << ", " << PV.Z() << std::endl;
					Logger(Logger::Debug) << "PV_RecoTau: " << PV_RecoTau.X() << ", " << PV_RecoTau.Y() << ", " << PV_RecoTau.Z() << std::endl;
					Logger(Logger::Debug) << "dxy, phi0, lam, dz: " << dxy << ", " << phi0 << ", " << lam << ", " << dz << std::endl;
					//Logger(Logger::Debug) << "acos(cosTheta2): " << acos(cosTheta2) << " TauMuDir.Theta(): " << TauMuDir.Theta() << std::endl;
					Logger(Logger::Debug) << "GEN Z.Pt(): " << GenZH.Pt() << " GenZH.Phi(): " << GenZH.Phi() << std::endl;

					Logger(Logger::Debug) << "TauDir.Phi()(SVPV): " << TauDir.Phi() << " TauMuDirNEW2.Phi(): " << TauMuDirNEW2.Phi() << std::endl;
					Logger(Logger::Debug) << "TauDir.Phi()(GEF): " << FitResults.getInitTauH().LV().Phi() << " Results.getInitTauMu().LV().Phi(): " << FitResults.getInitTauMu().LV().Phi() << std::endl;
					Logger(Logger::Debug) << "TauDir.Phi()(TPTR wtruth): " << RecoTauh_wTruth.LV().Phi() << std::endl;
					Logger(Logger::Debug) << "aPVSV: " << aPVSV << " bPVSV: " << bPVSV << std::endl;
					Logger(Logger::Debug) << "PV: " << PV.X() << " " << PV.Y() << " " << PV.Z() << std::endl;
					Logger(Logger::Debug) << "rot. SV: " << SV_Ntp.X() << " " << SV_Ntp.Y() << " " << SV_Ntp.Z() << std::endl;
					Logger(Logger::Debug) << "unrot. SV: " << SV.X() << " " << SV.Y() << " " << SV.Z() << std::endl;
					Logger(Logger::Debug) << "Mu POCA: " << Muon0.X() << " " << Muon0.Y() << " " << Muon0.Z() << std::endl;
					Logger(Logger::Debug) << "Mu POCA wrt PV: " << Ntp->Muon_Poca(selMuon).X() << " " << Ntp->Muon_Poca(selMuon).Y() << " " << Ntp->Muon_Poca(selMuon).Z() << std::endl;
					Logger(Logger::Debug) << "Mu Phi, Lambda, tNEW2: " << phi0 << " " << lam << " " << tNEW2 << std::endl;
					Logger(Logger::Debug) << "docNEW2: " << xdocNEW2 << " " << ydocNEW2 << " " << zdocNEW2 << std::endl;
					Logger(Logger::Debug) << "IS_PtBalance: " << IS_PtBalance.X() << " " << IS_PtBalance.Y() << " " << IS_PtBalance.Z() << std::endl;
					Logger(Logger::Debug) << "IS_EventRecoil: " << IS_EventRecoil.X() << " " << IS_EventRecoil.Y() << " " << IS_EventRecoil.Z() << std::endl;
					Logger(Logger::Debug) << "IS_wMET: " << IS_wMET.X() << " " << IS_wMET.Y() << " " << IS_wMET.Z() << std::endl;
					Logger(Logger::Debug) << "RecoTauh_wTruth.LV(): " << RecoTauh_wTruth.LV().X() << ", " << RecoTauh_wTruth.LV().Y() << std::endl;
					Logger(Logger::Debug) << "TauMu_PtBalance_PtVec: " << TauMu_PtBalance_PtVec.X() << ", " << TauMu_PtBalance_PtVec.Y() << std::endl;

					if(FitResults.Fitconverged()){
						if(GenZH.Pt() < 5.){
							if(FitResults.getResonance().Mass() >=0){
								Reco_dPhi_TauMuTauA1_AfterFit_lowBoost.at(t).Fill(FitResults.getTauMu().LV().DeltaPhi(FitResults.getTauH().LV()), w);
								Reco_dPhi_TauMuTauA1_BeforeFit_lowBoost.at(t).Fill(FitResults.getInitTauMu().LV().DeltaPhi(FitResults.getInitTauH().LV()), w);
								Reco_ZMass_UnboostedGenZ.at(t).Fill(FitResults.getResonance().Mass(), w);
								RecoZ_Pt_Unboosted.at(t).Fill(FitResults.getResonance().LV().Pt(), w);
								GenZ_Pt_Unboosted.at(t).Fill(GenZH.Pt(), w);
							}
							//if(FinalIndex != -1) Reco_ZMass_MassScanUnboosted.at(t).Fill(ZFits_MassScan.at(FinalIndex).at(IndicesToReturn.at(FinalIndex)).Mass(),w);
						}
						int IndexNoFit(FitResults.getIndex()); //always solution 1
						//if(IndexToReturn == 2) IndexNoFit = 1;
						//else IndexNoFit = IndexToReturn;
						TLorentzVector TauA1_TLV_FitSolution = FitResults.getTauH().LV();
						TLorentzVector Z_TLV = FitResults.getResonance().LV();
						TLorentzVector TauMu_TLV = FitResults.getTauMu().LV();
						//TVector3 par_Vec(par.at(IndexToReturn)(0),par.at(IndexToReturn)(1),par.at(IndexToReturn)(2));
						//TVector3 par0_Vec(par_0.at(IndexToReturn)(0),par_0.at(IndexToReturn)(1),par_0.at(IndexToReturn)(2));
						double dPt_TauA1 = TauA1_TLV_FitSolution.Pt() - GenTauh.Pt();
						double dPt_TauMu = TauMu_TLV.Pt() - GenTaumu.Pt();
						double dPt_TauA1_NoFit = FitResults.getInitTauH().LV().Pt() - GenTauh.Pt();
						double dPt_TauMu_NoFit = FitResults.getInitTauMu().LV().Pt() - GenTaumu.Pt();

						double dP_TauA1_FitSolution = TauA1_TLV_FitSolution.P() - GenTauh.P();
						double dE_Z = Z_TLV.E() - GenZH.E();
						/*
						double dP_TauA1_BestSolution(999); int i_BestSolution(-1);
						for(unsigned i=0; i<A1Fit.size(); i++){
							TLorentzVector TauA1_TLV_tmp = TPTF_TausA1.at(i).LV();
							if(fabs(dP_TauA1_BestSolution) > fabs(TauA1_TLV_tmp.P() - GenTauh.P())){
								dP_TauA1_BestSolution = TauA1_TLV_tmp.P() - GenTauh.P();
								i_BestSolution = i;
							}
						}
						*/
						TPTF_TauA1_RightSolution_vs_FitSolution.at(t).Fill(GenIndex, FitResults.getIndex());
						Gen_TPTF_TauA1_Solution_WithSelection.at(t).Fill(GenIndex, w);

						TrackParticle MuonTP = Ntp->Muon_TrackParticle(selMuon);
						/*
						if(A1Fit.at(0)){
							TVector2 ZEst = ZPtCollinearTauMuEstimator(MuonTP, TPTF_TausA1.at(0).LV(), Pt_Z.Phi_mpi_pi(Pt_Z.Phi()));
							Est_Z_Pt_wTruth.at(t).Fill(ZEst.Mod(), w);
							Est_Z_PtRes_wTruth.at(t).Fill(ZEst.Mod() - GenZH.Pt(), w);
							Est_Z_Pt_wTruth_vs_GenZ_Pt.at(t).Fill(ZEst.Mod(), GenZH.Pt());
							Est_Z_Pt_alwaysMinus.at(t).Fill(ZEst.Mod(), w);
							Est_Z_PtRes_alwaysMinus.at(t).Fill(ZEst.Mod() - GenZH.Pt(), w);
							Est_Z_Pt_alwaysMinus_vs_GenZ_Pt.at(t).Fill(ZEst.Mod(), GenZH.Pt());

							TLorentzVector TauMuEst = TauMuEstimator(TPTF_TausA1.at(0).LV(), Ntp->Muon_p4(selMuon));
							TLorentzVector FullZEst = TPTF_TausA1.at(0).LV() + TauMuEst;
							Est_Z_Energy_wTruth.at(t).Fill(FullZEst.E(), w);
							Est_Z_Energy_alwaysMinus.at(t).Fill(FullZEst.E(), w);
							Est_Z_EnergyRes_wTruth.at(t).Fill(FullZEst.E() - GenZH.E(), w);
							Est_Z_EnergyRes_alwaysMinus.at(t).Fill(FullZEst.E() - GenZH.E(), w);
							Est_TauMu_PtRes_wTruth.at(t).Fill(TauMuEst.Pt() - GenTaumu.Pt(), w);

							TLorentzVector TauMuEst2 = TauMuEstimator2(MuonTP, TPTF_TausA1.at(0).LV(), Pt_Z.Phi_mpi_pi(Pt_Z.Phi()));
							TLorentzVector FullZEst2 = TPTF_TausA1.at(0).LV() + TauMuEst2;
							Est_Z_EnergyRes_wTruth2.at(t).Fill(FullZEst2.E() - GenZH.E(), w);
							Est_Z_EnergyRes_alwaysMinus2.at(t).Fill(FullZEst2.E() - GenZH.E(), w);
							Est_TauMu_PtRes_wTruth2.at(t).Fill(TauMuEst2.Pt() - GenTaumu.Pt(), w);
							Est_Z_M_wTruth2.at(t).Fill(FullZEst2.M(), w);

							TLorentzVector TauMuEst_NoZMass = TauMuEstimatorNoZMass(MuonTP, TPTF_TausA1.at(0).LV(), Pt_Z.Phi_mpi_pi(Pt_Z.Phi()));
							TLorentzVector FullZEst_NoZMass = TPTF_TausA1.at(0).LV() + TauMuEst_NoZMass;
							Est_TauMu_PtRes_wTruth_NoZMass.at(t).Fill(TauMuEst_NoZMass.Pt() - GenTaumu.Pt(), w);
							Est_Z_EnergyRes_wTruth_NoZMass.at(t).Fill(FullZEst_NoZMass.E() - GenZH.E(), w);
							Est_Z_M_wTruth_NoZMass.at(t).Fill(FullZEst_NoZMass.M(), w);
						}
						else if(A1Fit.at(1) && A1Fit.at(2)){
							TVector2 ZEstMinus = ZPtCollinearTauMuEstimator(MuonTP, TPTF_TausA1.at(1).LV(), Pt_Z.Phi_mpi_pi(Pt_Z.Phi()));
							Est_Z_Pt_alwaysMinus.at(t).Fill(ZEstMinus.Mod(), w);
							Est_Z_PtRes_alwaysMinus.at(t).Fill(ZEstMinus.Mod() - GenZH.Pt(), w);
							Est_Z_Pt_alwaysMinus_vs_GenZ_Pt.at(t).Fill(ZEstMinus.Mod(), GenZH.Pt());

							TLorentzVector TauMuEstMinus = TauMuEstimator(TPTF_TausA1.at(1).LV(), Ntp->Muon_p4(selMuon));
							TLorentzVector FullZEstMinus = TPTF_TausA1.at(1).LV() + TauMuEstMinus;
							Est_Z_Energy_alwaysMinus.at(t).Fill(FullZEstMinus.E(), w);
							Est_Z_EnergyRes_alwaysMinus.at(t).Fill(FullZEstMinus.E() - GenZH.E(), w);

							TLorentzVector TauMuEst2Minus = TauMuEstimator2(MuonTP, TPTF_TausA1.at(1).LV(), Pt_Z.Phi_mpi_pi(Pt_Z.Phi()));
							TLorentzVector FullZEst2Minus = TPTF_TausA1.at(1).LV() + TauMuEst2Minus;
							Est_Z_EnergyRes_alwaysMinus2.at(t).Fill(FullZEst2Minus.E() - GenZH.E(), w);

							TLorentzVector TauMuEst_NoZMassMinus = TauMuEstimatorNoZMass(MuonTP, TPTF_TausA1.at(1).LV(), Pt_Z.Phi_mpi_pi(Pt_Z.Phi()));
							TLorentzVector FullZEst_NoZMassMinus = TPTF_TausA1.at(1).LV() + TauMuEst_NoZMassMinus;

							if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) < 0){
								TVector2 ZEstPlus = ZPtCollinearTauMuEstimator(MuonTP, TPTF_TausA1.at(2).LV(), Pt_Z.Phi_mpi_pi(Pt_Z.Phi()));
								Est_Z_Pt_wTruth.at(t).Fill(ZEstPlus.Mod(), w);
								Est_Z_PtRes_wTruth.at(t).Fill(ZEstPlus.Mod() - GenZH.Pt(), w);
								Est_Z_Pt_wTruth_vs_GenZ_Pt.at(t).Fill(ZEstPlus.Mod(), GenZH.Pt());

								TLorentzVector TauMuEstPlus = TauMuEstimator(TPTF_TausA1.at(2).LV(), Ntp->Muon_p4(selMuon));
								TLorentzVector FullZEstPlus = TPTF_TausA1.at(2).LV() + TauMuEstPlus;
								Est_Z_Energy_wTruth.at(t).Fill(FullZEstPlus.E(), w);
								Est_Z_EnergyRes_wTruth.at(t).Fill(FullZEstPlus.E() - GenZH.E(), w);
								Est_TauMu_PtRes_wTruth.at(t).Fill(TauMuEstPlus.Pt() - GenTaumu.Pt(), w);

								TLorentzVector TauMuEst2Plus = TauMuEstimator2(MuonTP, TPTF_TausA1.at(2).LV(), Pt_Z.Phi_mpi_pi(Pt_Z.Phi()));
								TLorentzVector FullZEst2Plus = TPTF_TausA1.at(2).LV() + TauMuEst2Plus;
								Est_Z_EnergyRes_wTruth2.at(t).Fill(FullZEst2Plus.E() - GenZH.E(), w);
								Est_TauMu_PtRes_wTruth2.at(t).Fill(TauMuEst2Plus.Pt() - GenTaumu.Pt(), w);
								Est_Z_M_wTruth2.at(t).Fill(FullZEst2Plus.M(), w);

								TLorentzVector TauMuEst_NoZMassPlus = TauMuEstimatorNoZMass(MuonTP, TPTF_TausA1.at(2).LV(), Pt_Z.Phi_mpi_pi(Pt_Z.Phi()));
								TLorentzVector FullZEst_NoZMassPlus = TPTF_TausA1.at(2).LV() + TauMuEst_NoZMassPlus;
								Est_TauMu_PtRes_wTruth_NoZMass.at(t).Fill(TauMuEst_NoZMassPlus.Pt() - GenTaumu.Pt(), w);
								Est_Z_EnergyRes_wTruth_NoZMass.at(t).Fill(FullZEst_NoZMassPlus.E() - GenZH.E(), w);
								Est_Z_M_wTruth_NoZMass.at(t).Fill(FullZEst_NoZMassPlus.M(), w);

							}
							else if(GenA1_boosted.Vect().Dot(GenTauh.Vect()) > 0){
								Est_Z_Pt_wTruth.at(t).Fill(ZEstMinus.Mod(), w);
								Est_Z_PtRes_wTruth.at(t).Fill(ZEstMinus.Mod() - GenZH.Pt(), w);
								Est_Z_Pt_wTruth_vs_GenZ_Pt.at(t).Fill(ZEstMinus.Mod(), GenZH.Pt());

								Est_Z_Energy_wTruth.at(t).Fill(FullZEstMinus.E(), w);
								Est_Z_EnergyRes_wTruth.at(t).Fill(FullZEstMinus.E() - GenZH.E(), w);
								Est_TauMu_PtRes_wTruth.at(t).Fill(TauMuEstMinus.Pt() - GenTaumu.Pt(), w);

								Est_Z_EnergyRes_wTruth2.at(t).Fill(FullZEst2Minus.E() - GenZH.E(), w);
								Est_TauMu_PtRes_wTruth2.at(t).Fill(TauMuEst2Minus.Pt() - GenTaumu.Pt(), w);
								Est_Z_M_wTruth2.at(t).Fill(FullZEst2Minus.M(), w);

								Est_TauMu_PtRes_wTruth_NoZMass.at(t).Fill(TauMuEst_NoZMassMinus.Pt() - GenTaumu.Pt(), w);
								Est_Z_EnergyRes_wTruth_NoZMass.at(t).Fill(FullZEst_NoZMassMinus.E() - GenZH.E(), w);
								Est_Z_M_wTruth_NoZMass.at(t).Fill(FullZEst_NoZMassMinus.M(), w);

							}
						}
						*/

						Reco_ZMass.at(t).Fill(Reco_Z_Corr.M(), w);
						Reco_Z_Energy_Res.at(t).Fill(Reco_Z_Corr.E() - GenZH.E(), w);
						Reco_PtRes_TauMu.at(t).Fill(dPt_TauMu, w);
						Reco_PtRes_TauA1.at(t).Fill(dPt_TauA1, w);
						Reco_PtRes_TauMu_NoFit.at(t).Fill(dPt_TauMu_NoFit, w);
						Reco_PtRes_TauA1_NoFit.at(t).Fill(dPt_TauA1_NoFit, w);

						double dPhi_TauMu_Prefit = FitResults.getInitTauMu().LV().DeltaPhi(GenTaumu);
						double dTheta_TauMu_Prefit = FitResults.getInitTauMu().LV().Theta() - GenTaumu.Theta();
						Reco_PhiRes_TauMu_PreFit.at(t).Fill(dPhi_TauMu_Prefit);
						Reco_ThetaRes_TauMu_PreFit.at(t).Fill(dTheta_TauMu_Prefit);

						if(FitResults.getIndex() == 0){
							Reco_PtRes_TauMu_AmbPoint0.at(t).Fill(dPt_TauMu, w);
							Reco_PtRes_TauA1_AmbPoint0.at(t).Fill(dPt_TauA1, w);
							Reco_PtRes_TauMu_AmbPoint0_NoFit.at(t).Fill(dPt_TauMu_NoFit, w);
							Reco_PtRes_TauA1_AmbPoint0_NoFit.at(t).Fill(dPt_TauA1_NoFit, w);
						}

						if(FitResults.getIndex() == 1 || FitResults.getIndex() == 2){
							Reco_PtRes_TauMu_AmbPoint12.at(t).Fill(dPt_TauMu, w);
							Reco_PtRes_TauA1_AmbPoint12.at(t).Fill(dPt_TauA1, w);
							Reco_PtRes_TauMu_AmbPoint12_NoFit.at(t).Fill(dPt_TauMu_NoFit, w);
							Reco_PtRes_TauA1_AmbPoint12_NoFit.at(t).Fill(dPt_TauA1_NoFit, w);
						}

						if(FitResults.getIndex() == 1){
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
						//TPTF_TauA1_pRes_BestSolution.at(t).Fill(dP_TauA1_BestSolution, w);
						//TPTF_TauA1_BestSolution_vs_FitSolution.at(t).Fill(i_BestSolution, IndexToReturn);
						TPTF_TauA1_pRes_vs_GenA1Mass_FitSolution.at(t).Fill(dP_TauA1_FitSolution, GenA1.M());
						//TPTF_TauA1_pRes_vs_GenA1Mass_BestSolution.at(t).Fill(dP_TauA1_BestSolution, GenA1.M());
						TPTF_TauA1_pRes_vs_GenGJAngle_FitSolution.at(t).Fill(dP_TauA1_FitSolution, GenGJAngle);
						//TPTF_TauA1_pRes_vs_GenGJAngle_BestSolution.at(t).Fill(dP_TauA1_BestSolution, GenGJAngle);
						TPTF_A1_pRes_vs_GenGJAngle.at(t).Fill(dP_A1, GenGJAngle);

						TPTF_TauA1_pRes_vs_RecoChi2_FitSolution.at(t).Fill(dP_TauA1_FitSolution, FitResults.getChi2());
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

						Gen_TauA1_P.at(t).Fill(GenTauh.P(), w);
						Gen_TauA1_Pt.at(t).Fill(GenTauh.Pt(), w);
						Gen_TauA1_Px.at(t).Fill(GenTauh.Px(), w);
						Gen_TauA1_Py.at(t).Fill(GenTauh.Py(), w);
						Gen_TauA1_Pz.at(t).Fill(GenTauh.Pz(), w);
						Gen_TauA1_Phi.at(t).Fill(GenTauh.Phi(), w);
						Gen_TauA1_Eta.at(t).Fill(GenTauh.Eta(), w);

						Gen_TauMu_P.at(t).Fill(GenTaumu.P(), w);
						Gen_TauMu_Pt.at(t).Fill(GenTaumu.Pt(), w);
						Gen_TauMu_Px.at(t).Fill(GenTaumu.Px(), w);
						Gen_TauMu_Py.at(t).Fill(GenTaumu.Py(), w);
						Gen_TauMu_Pz.at(t).Fill(GenTaumu.Pz(), w);
						Gen_TauMu_Phi.at(t).Fill(GenTaumu.Phi(), w);
						Gen_TauMu_Eta.at(t).Fill(GenTaumu.Eta(), w);

						Gen_Z_Pt.at(t).Fill(GenZH.Pt(), w);
						Gen_Z_Phi.at(t).Fill(GenZH.Phi(), w);
						Gen_Z_Eta.at(t).Fill(GenZH.Eta(), w);

						TPTF_TauA1_pxsq_Reco.at(t).Fill(pow(TauA1_TLV_FitSolution.Px(),2.), w);
						TPTF_TauA1_pysq_Reco.at(t).Fill(pow(TauA1_TLV_FitSolution.Py(),2.), w);
						TPTF_TauA1_pzsq_Reco.at(t).Fill(pow(TauA1_TLV_FitSolution.Pz(),2.), w);
						TPTF_TauA1_pxsq_Gen.at(t).Fill(pow(GenTauh.Px(),2.), w);
						TPTF_TauA1_pysq_Gen.at(t).Fill(pow(GenTauh.Py(),2.), w);
						TPTF_TauA1_pzsq_Gen.at(t).Fill(pow(GenTauh.Pz(),2.), w);

						TPTF_TauA1_ptRes_vs_ptGen.at(t).Fill(dPt_TauA1, GenTauh.Pt());
						TPTF_TauA1_ptRes_vs_ptReco.at(t).Fill(dPt_TauA1, TauA1_TLV_FitSolution.Pt());

						Reco_Z_PhiRes.at(t).Fill(Reco_Z_Corr.DeltaPhi(GenZH), w);
						Reco_Z_EtaRes.at(t).Fill(Reco_Z_Corr.Eta() - GenZH.Eta(), w);
						Reco_Z_PRes.at(t).Fill(Reco_Z_Corr.P() - GenZH.P(), w);
						Reco_Z_PtRes.at(t).Fill(Reco_Z_Corr.Pt() - GenZH.Pt(), w);

						if(FitResults.getIndex() == 0){
							Reco_Z_PhiRes_noAmb.at(t).Fill(Reco_Z_Corr.DeltaPhi(GenZH), w);
							Reco_Z_EtaRes_noAmb.at(t).Fill(Reco_Z_Corr.Eta() - GenZH.Eta(), w);
							Reco_Z_PRes_noAmb.at(t).Fill(Reco_Z_Corr.P() - GenZH.P(), w);
						}
						else{
							if(fitstatuses.at(1) && fitstatuses.at(2)){
								Reco_Z_PhiRes_wAmb.at(t).Fill(Reco_Z_Corr.DeltaPhi(GenZH), w);
								Reco_Z_EtaRes_wAmb.at(t).Fill(Reco_Z_Corr.Eta() - GenZH.Eta(), w);
								Reco_Z_PRes_wAmb.at(t).Fill(Reco_Z_Corr.P() - GenZH.P(), w);
								Reco_Chi2_rightAmb.at(t).Fill(FitResults.getChi2Vectors().at(GenIndex).Sum(), w);
								Reco_Chi2_orig_rightAmb.at(t).Fill(FitResults.getChi2Vectors().at(GenIndex)(0), w);
								Reco_Chi2_SC_rightAmb.at(t).Fill(FitResults.getChi2Vectors().at(GenIndex)(1), w);
								Reco_Chi2_HC_rightAmb.at(t).Fill(FitResults.getChi2Vectors().at(GenIndex)(2), w);
								Reco_Chi2_wrongAmb.at(t).Fill(FitResults.getChi2Vectors().at(GenWrongIndex).Sum(), w);
								Reco_Chi2_orig_wrongAmb.at(t).Fill(FitResults.getChi2Vectors().at(GenWrongIndex)(0), w);
								Reco_Chi2_SC_wrongAmb.at(t).Fill(FitResults.getChi2Vectors().at(GenWrongIndex)(1), w);
								Reco_Chi2_HC_wrongAmb.at(t).Fill(FitResults.getChi2Vectors().at(GenWrongIndex)(2), w);
							}
						}
						if(FitResults.getIndex() == GenIndex){
							Reco_Z_PhiRes_pickedrightAmb.at(t).Fill(Reco_Z_Corr.DeltaPhi(GenZH), w);
							Reco_Z_EtaRes_pickedrightAmb.at(t).Fill(Reco_Z_Corr.Eta() - GenZH.Eta(), w);
							Reco_Z_PRes_pickedrightAmb.at(t).Fill(Reco_Z_Corr.P() - GenZH.P(), w);
							Reco_NIter_pickedrightAmb.at(t).Fill(FitResults.getNiterations(), w);
							Reco_Chi2_pickedrightAmb.at(t).Fill(FitResults.getChi2(), w);
							Reco_Chi2_orig_pickedrightAmb.at(t).Fill(FitResults.getChi2Vector()(0), w);
							Reco_Chi2_SC_pickedrightAmb.at(t).Fill(FitResults.getChi2Vector()(1), w);
							Reco_Chi2_HC_pickedrightAmb.at(t).Fill(FitResults.getChi2Vector()(2), w);
						}

						if(FitResults.getIndex() != GenIndex){
							Reco_Z_PhiRes_pickedwrongAmb.at(t).Fill(Reco_Z_Corr.DeltaPhi(GenZH), w);
							Reco_Z_EtaRes_pickedwrongAmb.at(t).Fill(Reco_Z_Corr.Eta() - GenZH.Eta(), w);
							Reco_Z_PRes_pickedwrongAmb.at(t).Fill(Reco_Z_Corr.P() - GenZH.P(), w);
							Reco_NIter_pickedwrongAmb.at(t).Fill(FitResults.getNiterations(), w);
							Reco_Chi2_pickedwrongAmb.at(t).Fill(FitResults.getChi2(), w);
							Reco_Chi2_orig_pickedwrongAmb.at(t).Fill(FitResults.getChi2Vector()(0), w);
							Reco_Chi2_SC_pickedwrongAmb.at(t).Fill(FitResults.getChi2Vector()(1), w);
							Reco_Chi2_HC_pickedwrongAmb.at(t).Fill(FitResults.getChi2Vector()(2), w);
						}
						int index_bychi2_full, index_bychi2_orig, index_bySC, index_byHC, index_byorigplusSC;
						if(fitstatuses.at(1) && fitstatuses.at(2)){
							if(FitResults.getChi2Vectors().at(1).Sum() <= FitResults.getChi2Vectors().at(2).Sum()) index_bychi2_full = 1;
							if(FitResults.getChi2Vectors().at(1).Sum() > FitResults.getChi2Vectors().at(2).Sum()) index_bychi2_full = 2;
							if(FitResults.getChi2Vectors().at(1)(0) <= FitResults.getChi2Vectors().at(2)(0)) index_bychi2_orig = 1;
							if(FitResults.getChi2Vectors().at(1)(0) > FitResults.getChi2Vectors().at(2)(0)) index_bychi2_orig = 2;
							if(FitResults.getChi2Vectors().at(1)(1) <= FitResults.getChi2Vectors().at(2)(1)) index_bySC = 1;
							if(FitResults.getChi2Vectors().at(1)(1) > FitResults.getChi2Vectors().at(2)(1)) index_bySC = 2;
							if(FitResults.getChi2Vectors().at(1)(2) <= FitResults.getChi2Vectors().at(2)(2)) index_byHC = 1;
							if(FitResults.getChi2Vectors().at(1)(2) > FitResults.getChi2Vectors().at(2)(2)) index_byHC = 2;
							if( 	(FitResults.getChi2Vectors().at(1)(0) + FitResults.getChi2Vectors().at(1)(1))
								<=	(FitResults.getChi2Vectors().at(2)(0) + FitResults.getChi2Vectors().at(2)(1))
								) index_byorigplusSC = 1;
							if( 	(FitResults.getChi2Vectors().at(1)(0) + FitResults.getChi2Vectors().at(1)(1))
								>	(FitResults.getChi2Vectors().at(2)(0) + FitResults.getChi2Vectors().at(2)(1))
								) index_byorigplusSC = 2;
							Reco_FitSolution_byChi2_Full_vs_RightSolution.at(t).Fill(GenIndex, index_bychi2_full);
							Reco_FitSolution_byChi2_orig_vs_RightSolution.at(t).Fill(GenIndex, index_bychi2_orig);
							Reco_FitSolution_byChi2_SC_vs_RightSolution.at(t).Fill(GenIndex, index_bySC);
							Reco_FitSolution_byChi2_HC_vs_RightSolution.at(t).Fill(GenIndex, index_byHC);
							Reco_FitSolution_byChi2_origplusSC_vs_RightSolution.at(t).Fill(GenIndex, index_byorigplusSC);

							double chi2_diff = FitResults.getChi2Vectors().at(2).Sum() - FitResults.getChi2Vectors().at(1).Sum();
							Reco_Chi2_diff.at(t).Fill(chi2_diff, w);
							Reco_Chi2_orig_diff.at(t).Fill(FitResults.getChi2Vectors().at(2)(0) - FitResults.getChi2Vectors().at(1)(0), w);
							Reco_Chi2_SC_diff.at(t).Fill(FitResults.getChi2Vectors().at(2)(1) - FitResults.getChi2Vectors().at(1)(1), w);
							Reco_Chi2_HC_diff.at(t).Fill(FitResults.getChi2Vectors().at(2)(2) - FitResults.getChi2Vectors().at(1)(2), w);
							Reco_Chi2_origplusSC_diff.at(t).Fill(FitResults.getChi2Vectors().at(2)(0) + FitResults.getChi2Vectors().at(2)(1) - FitResults.getChi2Vectors().at(1)(0) - FitResults.getChi2Vectors().at(1)(1), w);

							if(chi2_diff >= 0){
								if(GenIndex == 1) Reco_Chi2_diff_vs_correctAssignment.at(t).Fill(chi2_diff, 1.);
								else if(GenIndex == 2) Reco_Chi2_diff_vs_correctAssignment.at(t).Fill(chi2_diff, 0);
							}
							else{
								if(GenIndex == 2) Reco_Chi2_diff_vs_correctAssignment.at(t).Fill(chi2_diff, 1.);
								else if(GenIndex == 1) Reco_Chi2_diff_vs_correctAssignment.at(t).Fill(chi2_diff, 0);
							}
						}


					}
					if(FitResultswithFullRecoil.Fitconverged()){
						Reco_PtRes_TauA1_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Pt() - GenTauh.Pt(), w);
						Reco_PtRes_TauMu_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Pt() - GenTaumu.Pt(), w);
						Reco_PtRes_TauA1_wRecoil_PreFit.at(t).Fill(FitResultswithFullRecoil.getInitTauH().LV().Pt() - GenTauh.Pt(), w);
						Reco_PtRes_TauMu_wRecoil_PreFit.at(t).Fill(FitResultswithFullRecoil.getInitTauMu().LV().Pt() - GenTaumu.Pt(), w);
						Reco_Z_Energy_Res_wRecoil.at(t).Fill(FitResultswithFullRecoil.getResonance().LV().E() - GenZH.E(), w);
						Reco_ZMass_wRecoil.at(t).Fill(Reco_Z_Corr_wRecoil.M(), w);
						Reco_TauA1_ptRes_vs_ptGen_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Pt() - GenTauh.Pt(), GenTauh.Pt());
						Reco_TauA1_ptRes_vs_ptReco_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Pt() - GenTauh.Pt(), FitResultswithFullRecoil.getTauH().LV().Pt());
						Reco_TauMu_ptRes_vs_ptGen_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Pt() - GenTaumu.Pt(), GenTaumu.Pt());
						Reco_TauMu_ptRes_vs_ptReco_wRecoil.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Pt() - GenTaumu.Pt(), FitResultswithFullRecoil.getInitTauMu().LV().Pt());

						Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr.at(t).Fill(GenTauh.Pt() - FitResultswithFullRecoil.getTauH().LV().Pt(), FitResultswithFullRecoil.getTauH().LV().Pt());
						Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr.at(t).Fill(GenTaumu.Pt() - FitResultswithFullRecoil.getTauMu().LV().Pt(), FitResultswithFullRecoil.getInitTauMu().LV().Pt());

						Reco_Z_PtRes_wRecoil.at(t).Fill(Reco_Z_Corr_wRecoil.Pt() - GenZH.Pt(), w);
						Reco_Z_EtaRes_wRecoil.at(t).Fill(Reco_Z_Corr_wRecoil.Eta() - GenZH.Eta(), w);
						Reco_Z_PhiRes_wRecoil.at(t).Fill(Reco_Z_Corr_wRecoil.DeltaPhi(GenZH), w);

						if(FitResultswithFullRecoil.getIndex() == 0){
							Reco_PtRes_TauA1_wRecoil_AmbZero.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Pt() - GenTauh.Pt(), w);
							Reco_PtRes_TauMu_wRecoil_AmbZero.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Pt() - GenTaumu.Pt(), w);
							Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_AmbZero.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Pt() - GenTaumu.Pt(), FitResultswithFullRecoil.getTauMu().LV().Pt());
							Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_AmbZero.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Pt() - GenTauh.Pt(), FitResultswithFullRecoil.getTauH().LV().Pt());
							Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit.at(t).Fill(FitResultswithFullRecoil.getInitTauHs().at(0).LV().Pt() - GenTauh.Pt(), FitResultswithFullRecoil.getInitTauHs().at(0).LV().Pt());
							Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb0_Prefit.at(t).Fill(FitResultswithFullRecoil.getInitTauMus().at(0).LV().Pt() - GenTaumu.Pt(), FitResultswithFullRecoil.getInitTauMus().at(0).LV().Pt());
						}
						if(FitResultswithFullRecoil.getIndex() != 0){
							Reco_PtRes_TauA1_wRecoil_wAmb.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Pt() - GenTauh.Pt(), w);
							Reco_PtRes_TauMu_wRecoil_wAmb.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Pt() - GenTaumu.Pt(), w);
							Reco_TauMu_ptRes_vs_ptReco_wRecoilCorr_wAmb.at(t).Fill(FitResultswithFullRecoil.getTauMu().LV().Pt() - GenTaumu.Pt(), FitResultswithFullRecoil.getTauMu().LV().Pt());
							Reco_TauA1_ptRes_vs_ptReco_wRecoilCorr_wAmb.at(t).Fill(FitResultswithFullRecoil.getTauH().LV().Pt() - GenTauh.Pt(), FitResultswithFullRecoil.getTauH().LV().Pt());

							if(fitstatuses_wRecoil.at(1) && fitstatuses_wRecoil.at(2)){
								if(GenIndex == 1){
									Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit.at(t).Fill(FitResultswithFullRecoil.getInitTauHs().at(GenIndex).LV().Pt() - GenTauh.Pt(), FitResultswithFullRecoil.getInitTauHs().at(GenIndex).LV().Pt());
									Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb1_Prefit.at(t).Fill(FitResultswithFullRecoil.getInitTauMus().at(GenIndex).LV().Pt() - GenTaumu.Pt(), FitResultswithFullRecoil.getInitTauMus().at(GenIndex).LV().Pt());
								}
								if(GenIndex == 2){
									Reco_TauA1_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit.at(t).Fill(FitResultswithFullRecoil.getInitTauHs().at(GenIndex).LV().Pt() - GenTauh.Pt(), FitResultswithFullRecoil.getInitTauHs().at(GenIndex).LV().Pt());
									Reco_TauMu_ptRes_vs_ptReco_wRecoilwTruth_Amb2_Prefit.at(t).Fill(FitResultswithFullRecoil.getInitTauMus().at(GenIndex).LV().Pt() - GenTaumu.Pt(), FitResultswithFullRecoil.getInitTauMus().at(GenIndex).LV().Pt());
								}
							}
						}
						int index_bychi2_full, index_bychi2_orig, index_bySC, index_byHC, index_byorigplusSC;
						if(fitstatuses_wRecoil.at(1) && fitstatuses_wRecoil.at(2)){
							if(FitResultswithFullRecoil.getChi2Vectors().at(1).Sum() <= FitResultswithFullRecoil.getChi2Vectors().at(2).Sum()) index_bychi2_full = 1;
							if(FitResultswithFullRecoil.getChi2Vectors().at(1).Sum() > FitResultswithFullRecoil.getChi2Vectors().at(2).Sum()) index_bychi2_full = 2;
							if(FitResultswithFullRecoil.getChi2Vectors().at(1)(0) <= FitResultswithFullRecoil.getChi2Vectors().at(2)(0)) index_bychi2_orig = 1;
							if(FitResultswithFullRecoil.getChi2Vectors().at(1)(0) > FitResultswithFullRecoil.getChi2Vectors().at(2)(0)) index_bychi2_orig = 2;
							if(FitResultswithFullRecoil.getChi2Vectors().at(1)(1) <= FitResultswithFullRecoil.getChi2Vectors().at(2)(1)) index_bySC = 1;
							if(FitResultswithFullRecoil.getChi2Vectors().at(1)(1) > FitResultswithFullRecoil.getChi2Vectors().at(2)(1)) index_bySC = 2;
							if(FitResultswithFullRecoil.getChi2Vectors().at(1)(2) <= FitResultswithFullRecoil.getChi2Vectors().at(2)(2)) index_byHC = 1;
							if(FitResultswithFullRecoil.getChi2Vectors().at(1)(2) > FitResultswithFullRecoil.getChi2Vectors().at(2)(2)) index_byHC = 2;
							if( 	(FitResultswithFullRecoil.getChi2Vectors().at(1)(0) + FitResultswithFullRecoil.getChi2Vectors().at(1)(1))
								<=	(FitResultswithFullRecoil.getChi2Vectors().at(2)(0) + FitResultswithFullRecoil.getChi2Vectors().at(2)(1))
								) index_byorigplusSC = 1;
							if( 	(FitResultswithFullRecoil.getChi2Vectors().at(1)(0) + FitResultswithFullRecoil.getChi2Vectors().at(1)(1))
								>	(FitResultswithFullRecoil.getChi2Vectors().at(2)(0) + FitResultswithFullRecoil.getChi2Vectors().at(2)(1))
								) index_byorigplusSC = 2;
							Reco_FitSolution_byChi2_Full_vs_RightSolution_wRecoil.at(t).Fill(GenIndex, index_bychi2_full);
							Reco_FitSolution_byChi2_orig_vs_RightSolution_wRecoil.at(t).Fill(GenIndex, index_bychi2_orig);
							Reco_FitSolution_byChi2_SC_vs_RightSolution_wRecoil.at(t).Fill(GenIndex, index_bySC);
							Reco_FitSolution_byChi2_HC_vs_RightSolution_wRecoil.at(t).Fill(GenIndex, index_byHC);
							Reco_FitSolution_byChi2_origplusSC_vs_RightSolution_wRecoil.at(t).Fill(GenIndex, index_byorigplusSC);
						}
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
			//std::cout << "Ntp->PFTau_NdaughtersReFitTracks_p4(selTau) " << Ntp->PFTau_NdaughtersReFitTracks_p4(selTau) << std::endl;
			for(unsigned int i=0; i<Ntp->PFTau_NdaughtersReFitTracks_p4(selTau);i++){
				//std::cout << "Ntp->PFTau_daughterReFitTracks_p4(selTau).at(i).M() " << Ntp->PFTau_daughterReFitTracks_p4(selTau).at(i).M() << std::endl;
				Tau_Mass_Inclusive_ReFitTracks.at(t).Fill(Ntp->PFTau_daughterReFitTracks_p4(selTau).at(i).M(), w);
			}
			if(pass.at(TauFLSigma)){
				TLorentzVector PFTau_UnFitTracks;
				if(Ntp->PFTau_NdaughterTracks(selTau) == 3){
					for(unsigned int i=0; i<Ntp->PFTau_NdaughterTracks(selTau);i++){
						//std::cout << "Ntp->PFTau_daughterTracks(selTau).size()" << Ntp->PFTau_daughterTracks(selTau).size() << std::endl;
						TrackParticle tmpTP = Ntp->PFTau_daughterTracks(selTau).at(i);
						TVector3 SV = Ntp->PFTau_TIP_secondaryVertex_pos(selTau);
						TLorentzVector tmpLV = (TrackTools::LorentzParticleAtPosition(tmpTP, SV)).LV();
						Tau_Mass_Inclusive_UnFitTracks.at(t).Fill(tmpLV.M(), w);
					}
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
	if(passAllBut(TauFLSigma) && Ntp->PFTau_TIP_hassecondaryVertex(selTau) && pass.at(PrimeVtx)){
	  double FLSigmaVlad = Ntp->PFTau_FlightLength_significance(Ntp->PFTau_TIP_primaryVertex_pos(selTau), Ntp->PFTau_TIP_primaryVertex_cov(selTau), Ntp->PFTau_TIP_secondaryVertex_pos(selTau), Ntp->PFTau_TIP_secondaryVertex_cov(selTau));
	  double sign(0);
	  if(Ntp->PFTau_3PS_A1_LV(selTau).Vect().Dot(Ntp->PFTau_FlightLength3d(selTau)) < 0){
			sign = -1;
		}
		else{
			sign = 1;
		}
	  TauFLSigmaVlad.at(t).Fill(sign*FLSigmaVlad, w);
	  //TauFLSigmaAlex.at(t).Fill(sign*Ntp->PFTau_FlightLength_significance(selTau), w);

	  if(FLSigmaVlad >= cut.at(TauFLSigma)){
		TauFLSigmaVlad_PhiA1.at(t).Fill(Ntp->PFTau_p4(selTau).Phi(), w);
		if(FitResults.Fitconverged()){
			TauFLSigmaVlad_PhiTau.at(t).Fill(FitResults.getTauH().LV().Phi(), w);
			TauFLSigmaVlad_PhiTauNoCorr.at(t).Fill(FitResultsNoCorr.getTauH().LV().Phi(), w);
			TauFLSigmaVlad_PhiZnoCorr.at(t).Fill(FitResults.getResonance().LV().Phi(), w);
			TauFLSigmaVlad_PhiZwCorr.at(t).Fill(Reco_Z_Corr.Phi(), w);
		}
	  }
	  /*
	  if(Ntp->PFTau_FlightLength_significance(selTau) >= cut.at(TauFLSigma)){
		TauFLSigmaAlex_PhiA1.at(t).Fill(Ntp->PFTau_p4(selTau).Phi(), w);
		if(FitResults.Fitconverged()){
			TauFLSigmaAlex_PhiTau.at(t).Fill(FitResults.getTauH().LV().Phi(), w);
			TauFLSigmaAlex_PhiTauNoCorr.at(t).Fill(FitResultsNoCorr.getTauH().LV().Phi(), w);
			TauFLSigmaAlex_PhiZnoCorr.at(t).Fill(FitResults.getResonance().LV().Phi(), w);
			TauFLSigmaAlex_PhiZwCorr.at(t).Fill(Reco_Z_Corr.Phi(), w);
		}
	  }
	  */
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
bool ZToTaumuTauh::selectMuon_DiMuonVeto(unsigned i, unsigned i_vtx){
	if(	Ntp->Muon_isGlobalMuon(i) &&
		Ntp->Muon_isPFMuon(i) &&
		Ntp->Muon_p4(i).Pt() >= cDiMuVeto_pt &&
		fabs(Ntp->Muon_p4(i).Eta()) <= cDiMuVeto_eta &&
		Ntp->Muon_isTrackerMuon(i) &&
		Ntp->Muon_RelIso(i) < cDiMuVeto_dBetaCombIso &&
		Ntp->dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx)) < cDiMuVeto_dz){
		return true;
	}
	return false;
}
/*
LorentzVectorParticle ZToTaumuTauh::Reconstruct_hadronicTau(unsigned i, unsigned int Ambiguity){

}
TMatrixT<double> ZToTaumuTauh::ComputeTauDirection(TMatrixT<double> inpar){
  TMatrixT<double> outpar(3,1);

}
double ZToTaumuTauh::Reconstruct_hadronicTauEnergy(unsigned i, unsigned int Ambiguity, bool UseHelix){
	double TauEnergy,TauMomentumPlus, TauMomentumMinus;
	TLorentzVector A1 = Ntp->PFTau_p4(i);
	double GJ_angle = A1.Angle(Ntp->PFTau_FlightLength3d(i));
	double GJ_angleMax = GJAngleMax(A1);
	if(GJ_angle > GJ_angleMax){

	}
	else{
		double val1 = (pow(PDG_Var::Tau_mass(),2.) + pow(A1.M(),2.))*A1.P()*cos(GJ_angle);
		double val2 = A1.Energy()*sqrt(pow(pow(PDG_Var::Tau_mass(),2.) - pow(A1.M(),2.),2.) - 4.*pow(A1.P()*PDG_Var::Tau_mass()*sin(GJ_angle),2.));
		TauMomentumPlus = (val1 + val2)/2./(pow(A1.M(),2) + pow(A1.P()*sin(GJ_angle),2.));
		TauMomentumMinus = (val1 - val2)/2./(pow(A1.M(),2) + pow(A1.P()*sin(GJ_angle),2.));
	}
	return TauEnergy;
}
*/
double ZToTaumuTauh::GJAngleMax(TLorentzVector A1){
	double arg = (pow(PDGInfo::tau_mass(),2.) - A1.M2())/2./A1.P()/PDGInfo::tau_mass();
	return asin(arg);
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
TLorentzVector ZToTaumuTauh::TauHelixP4AtSV(unsigned int selTau, TLorentzVector Tau){
	double pTau = Tau.P();
	TVector3 PV = Ntp->PFTau_TIP_primaryVertex_pos(selTau);
	TVector3 SV = Ntp->PFTau_TIP_secondaryVertex_pos(selTau);
	TVector3 PVSV = TVector3(SV - PV);
	TVector3 p0 = TVector3(PVSV); double p0Mag = p0.Mag(); double p0Scale = pTau/p0Mag; p0 = p0*p0Scale;
	//TVector3 p0 = TVector3(Tau.Vect());
	double charge = ((double)Ntp->PFTau_Charge(selTau));
	double alpha = -charge*3*pow(10.,-3.)*3.8;
	double Tau_R = Tau.Pt()/alpha;
	double kappa = 1/Tau_R/2.;
	double Tau_Phi0 = atan2(p0.Y(),p0.X()) + charge*(TMath::PiOver2() - acos(PVSV.Mag()*kappa));
	double Tau_Phi02 = atan2(p0.Y(),p0.X()) - charge*(TMath::PiOver2() - acos(PVSV.Mag()*kappa));
	TLorentzVector Helix_TLV = TLorentzVector(p0.Pt()*cos(Tau_Phi02), p0.Pt()*sin(Tau_Phi02), p0.Z(), Tau.E());
	//Logger(Logger::Debug) << "Helix_TLV.M(): " << Helix_TLV.M() << std::endl;
	return Helix_TLV;
}
TLorentzVector ZToTaumuTauh::GenTauHelixP4AtSV(unsigned int selTau, TLorentzVector Tau){
	double pTau = Tau.P();
	TVector3 PV = Ntp->MCTauandProd_Vertex(selTau,0);
	TVector3 SV = Ntp->MCTauandProd_Vertex(selTau,1);
	TVector3 PVSV = TVector3(SV - PV);
	TVector3 p0 = TVector3(PVSV); double p0Mag = p0.Mag(); double p0Scale = pTau/p0Mag; p0 = p0*p0Scale;
	//TVector3 p0 = TVector3(Tau.Vect());
	double charge = ((double)Ntp->MCTau_charge(selTau));
	double alpha = -charge*3*pow(10.,-3.)*3.8;
	double Tau_R = Tau.Pt()/alpha;
	double kappa = 1/Tau_R/2.;
	double Tau_Phi0 = atan2(p0.Y(),p0.X()) + charge*(TMath::PiOver2() - acos(PVSV.Mag()*kappa));
	double Tau_Phi02 = atan2(p0.Y(),p0.X()) - charge*(TMath::PiOver2() - acos(PVSV.Mag()*kappa));
	TLorentzVector Helix_TLV = TLorentzVector(p0.Pt()*cos(Tau_Phi02), p0.Pt()*sin(Tau_Phi02), p0.Z(), Tau.E());
	//Logger(Logger::Debug) << "Helix_TLV.M(): " << Helix_TLV.M() << std::endl;
	return Helix_TLV;
}
TVector2 ZToTaumuTauh::ZPtCollinearTauMuEstimator(TrackParticle Muon, TLorentzVector Tauh, double PhiRecoil){
	double dxy 			= Muon.Parameter(TrackParticle::dxy);
	double phi0 		= Muon.Parameter(TrackParticle::phi);
	double lambda 		= Muon.Parameter(TrackParticle::lambda);
	double dz 			= Muon.Parameter(TrackParticle::dz);
	double cosphi0 		= cos(phi0); 				double sinphi0 		= sin(phi0);

	double PhiR = (PhiRecoil > 0) ? PhiRecoil - TMath::Pi() : PhiRecoil + TMath::Pi();
	double cosPhiR = cos(PhiR); double sinPhiR = sin(PhiR);
	double Pt_Z = (Tauh.X()/cosphi0 - Tauh.Y()/sinphi0)/(cosPhiR/cosphi0 - sinPhiR/sinphi0);
	TVector2 Pt_ZVec = TVector2(Pt_Z*cosPhiR,Pt_Z*sinPhiR);

	return Pt_ZVec;
}
TLorentzVector ZToTaumuTauh::TauMuEstimator(TLorentzVector Tauh, TLorentzVector Muon){
	double Zmass= 91.2;
	double pTauMu = pow(Zmass, 2.)/2./Tauh.P()/(1 - cos(Tauh.Angle(Muon.Vect())));
	TLorentzVector TauMu = TLorentzVector(pTauMu*sin(Muon.Theta())*cos(Muon.Phi()), pTauMu*sin(Muon.Theta())*sin(Muon.Phi()), pTauMu*cos(Muon.Theta()), sqrt(pow(pTauMu, 2.) + pow(PDGInfo::tau_mass(), 2.)));
	Logger(Logger::Debug) << "TauMu1.Pt(): " << TauMu.Pt() << std::endl;
	return TauMu;
}
TLorentzVector ZToTaumuTauh::TauMuEstimator2(TrackParticle Muon, TLorentzVector Tauh, double PhiRecoil){
	double dxy 			= Muon.Parameter(TrackParticle::dxy);
	double phi0 		= Muon.Parameter(TrackParticle::phi);
	double lambda 		= Muon.Parameter(TrackParticle::lambda);
	double dz 			= Muon.Parameter(TrackParticle::dz);
	double cosphi0 		= cos(phi0); 				double sinphi0 		= sin(phi0);
	double Zmass= 91.2;
	TVector3 MuonDir = TVector3(cosphi0*cos(lambda), sinphi0*cos(lambda),sin(lambda));

	double PhiR = (PhiRecoil > 0) ? PhiRecoil - TMath::Pi() : PhiRecoil + TMath::Pi();
	double cosPhiR = cos(PhiR); double sinPhiR = sin(PhiR);
	double Pt_TauMu = -(Tauh.X()/cosPhiR - Tauh.Y()/sinPhiR)/(cosphi0/cosPhiR - sinphi0/sinPhiR);
	double P_TauMu = pow(Zmass, 2.)/2./Tauh.P()/(1 - cos(Tauh.Angle(MuonDir)));
	double Pz_TauMu = sqrt(pow(P_TauMu, 2.) - pow(Pt_TauMu, 2.));
	Logger(Logger::Debug) << "TauMu2.Pt(): " << Pt_TauMu << std::endl;
	TLorentzVector TauMu = TLorentzVector(Pt_TauMu*cosphi0, Pt_TauMu*sinphi0, Pz_TauMu, sqrt(pow(P_TauMu, 2.) + pow(PDGInfo::tau_mass(), 2.)));
	return TauMu;
}
TLorentzVector ZToTaumuTauh::TauMuEstimatorNoZMass(TrackParticle Muon, TLorentzVector Tauh, double PhiRecoil){
	double dxy 			= Muon.Parameter(TrackParticle::dxy);
	double phi0 		= Muon.Parameter(TrackParticle::phi);
	double lambda 		= Muon.Parameter(TrackParticle::lambda);
	double dz 			= Muon.Parameter(TrackParticle::dz);
	double cosphi0 		= cos(phi0); 				double sinphi0 		= sin(phi0);

	TVector3 MuonDir = TVector3(cosphi0, sinphi0, tan(lambda)); // with mag = 1/cos(lambda) for ez scaling

	double PhiR = (PhiRecoil > 0) ? PhiRecoil - TMath::Pi() : PhiRecoil + TMath::Pi();
	double cosPhiR = cos(PhiR); double sinPhiR = sin(PhiR);
	double Pt_TauMu = -(Tauh.X()/cosPhiR - Tauh.Y()/sinPhiR)/(cosphi0/cosPhiR - sinphi0/sinPhiR);

	TVector3 TauMuPVec = MuonDir*Pt_TauMu;
	double TauMuE = sqrt(TauMuPVec.Mag2() + pow(PDGInfo::tau_mass(), 2.));
	Logger(Logger::Debug) << "TauMu2.Pt(): " << Pt_TauMu << std::endl;
	TLorentzVector TauMu = TLorentzVector(TauMuPVec, TauMuE);
	return TauMu;
}
TLorentzVector ZToTaumuTauh::TauMuFullEstimate(TVector3 PV, TrackParticle Muon, LorentzVectorParticle Tauh, TVector2 TauMuPt, TVector3 &Intersection){
	double dxy 			= Muon.Parameter(TrackParticle::dxy);
	double phi0 		= Muon.Parameter(TrackParticle::phi);
	double lam	 		= Muon.Parameter(TrackParticle::lambda);
	double dz			= Muon.Parameter(TrackParticle::dz);

	double tanphi_tau = TauMuPt.Y()/TauMuPt.X();
	double a = tanphi_tau;
	double b = PV.Y() - a*PV.X();

	double deno = a*cos(phi0) - sin(phi0);
	double num = dxy*(cos(phi0) + a*sin(phi0)) - b;

	double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b)  / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

	double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
	double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
	double zdoc = dz + tmu*tan(lam);

	TVector3 PointGuess(xdoc, ydoc, zdoc);
	TVector3 TauMuDir = PointGuess - PV;
	TauMuDir.SetMag(1.);

	Intersection = PointGuess;

	double Zmass = 91.2;
	double P_TauMu = pow(Zmass, 2.)/2./Tauh.LV().P()/(1 - cos(Tauh.LV().Angle(TauMuDir)));

	TLorentzVector Taumu(P_TauMu*TauMuDir, sqrt(pow(1.777,2.) + pow(P_TauMu, 2.)));


	Logger(Logger::Debug) << "acos(TauMuDir.Z()/TauMuDir.Mag()): " << acos(TauMuDir.Z()/TauMuDir.Mag()) << " TauMuDir.Theta(): " << TauMuDir.Theta() << std::endl;
	Logger(Logger::Debug) << "a: " << a << " b: " << b << std::endl;
	Logger(Logger::Debug) << "tmu: " << tmu << std::endl;
	Logger(Logger::Debug) << "numerator of tmu: " << num << " denominator of tmu: " << deno << std::endl;

	Logger(Logger::Debug) << "Muon Charge: " << Muon.Charge() << " Kappa: "<< Muon.Parameter(TrackParticle::kappa) <<  std::endl;
	Logger(Logger::Debug) << "Muon dxy: " << dxy << " dz: "<< dz <<  std::endl;
	Logger(Logger::Debug) << "Muon Phi: " << phi0 << std::endl;
	Logger(Logger::Debug) << "Tauh.LV().Phi(): " << Tauh.LV().Phi() << " Tauh.LV().Pt(): " << Tauh.LV().Pt() << std::endl;
	Logger(Logger::Debug) << "TauMuPt.X(): "<< TauMuPt.X() << " TauMuPt.Y(): " << TauMuPt.Y() << std::endl;
	Logger(Logger::Debug) << "TauMuPt.Phi(): "<< TauMuPt.Phi_mpi_pi(TauMuPt.Phi()) << " TauMuPt.Mod(): " << TauMuPt.Mod() << std::endl;
	Logger(Logger::Debug) << "Taumu.Phi(): "<< Taumu.Phi() << " Taumu.Pt(): " << Taumu.Pt() << std::endl;
	Logger(Logger::Debug) << "P_TauMu: "<< P_TauMu << " Taumu.P(): " << Taumu.P() << std::endl;
	Logger(Logger::Debug) << "Taumu.M(): "<< Taumu.M() << std::endl;
	TLorentzVector Z = Taumu + Tauh.LV();
	Logger(Logger::Debug) << "Z.M(): " << Z.M() << " Z.Pt(): " << Z.Pt() <<  " Z.Phi(): " << Z.Phi() << std::endl;

	return Taumu;
}

TLorentzVector ZToTaumuTauh::ZEstimatorWithMET(bool rotate, TVector3 PV, TrackParticle Muon, LorentzVectorParticle A1, LorentzVectorParticle Tauh, TVector2 MET){
	double dxy 			= Muon.Parameter(TrackParticle::dxy);
	double phi0 		= Muon.Parameter(TrackParticle::phi);
	double lam	 		= Muon.Parameter(TrackParticle::lambda);
	double dz 			= Muon.Parameter(TrackParticle::dz);

	TVector3 MuonDir(cos(phi0), sin(phi0), tan(lam));
	TVector3 Muon0(-dxy*sin(phi0), dxy*cos(phi0), dz);
	TVector3 P_Tauh = Tauh.LV().Vect();

	double MuonPt = fabs(Muon.BField()/Muon.Parameter(TrackParticle::kappa));

	TVector2 MuonPtVec(MuonPt*MuonDir.XYvector());
	TVector2 A1PtVec(A1.LV().Vect().XYvector());
	TVector2 ZPtVec = MET + MuonPtVec + A1PtVec;
	TVector2 TauMuPtVec = ZPtVec - Tauh.LV().Vect().XYvector();

	double tanphi_tau = TauMuPtVec.Y()/TauMuPtVec.X();
	double a = tanphi_tau;
	double b = PV.Y() - a*PV.X();

	double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b) / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

	double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
	double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
	double zdoc = dz + tmu*tan(lam);

	TVector3 PointGuess(xdoc, ydoc, zdoc);
	TVector3 TauMuDir = PointGuess - PV;


	if(rotate){
	  if(fabs(TauMuDir.Phi()-P_Tauh.Phi())<1.5){
		TauMuPtVec.Set(-TauMuPtVec.X(),-TauMuPtVec.Y());
	  }
	}
	TLorentzVector TauMu; TauMu.SetXYZM(TauMuPtVec.X(),TauMuPtVec.Y(),TauMuPtVec.Mod()/tan(TauMuDir.Theta()), PDGInfo::tau_mass());
	TLorentzVector Z(TauMu + Tauh.LV());
	Logger(Logger::Debug) << "Reco Z Px with MET without neutrino only Pt: " << ZPtVec.X() << std::endl;
	Logger(Logger::Debug) << "Reco Z Py with MET without neutrino only Pt: " << ZPtVec.Y() << std::endl;
	Logger(Logger::Debug) << "Reco Z Px with MET without neutrino Full Z: " << Z.X() << std::endl;
	Logger(Logger::Debug) << "Reco Z Py with MET without neutrino Full Z: " << Z.Y() << std::endl;

	return Z;
}
TLorentzVector ZToTaumuTauh::ZEstimatorB2B(bool rotate, TVector3 PV, TrackParticle Muon, LorentzVectorParticle A1, LorentzVectorParticle Tauh){
	double dxy 			= Muon.Parameter(TrackParticle::dxy);
	double phi0 		= Muon.Parameter(TrackParticle::phi);
	double lam	 		= Muon.Parameter(TrackParticle::lambda);
	double dz 			= Muon.Parameter(TrackParticle::dz);

	TVector3 MuonDir(cos(phi0), sin(phi0), tan(lam));
	TVector3 Muon0(-dxy*sin(phi0), dxy*cos(phi0), dz);
	TVector3 P_Tauh = Tauh.LV().Vect();

	double MuonPt = fabs(Muon.BField()/Muon.Parameter(TrackParticle::kappa));

	TVector2 MuonPtVec(MuonPt*MuonDir.XYvector());
	TVector2 A1PtVec(A1.LV().Vect().XYvector());
	TVector2 TauMuPtVec(-P_Tauh.X(), -P_Tauh.Y());

	double tanphi_tau = TauMuPtVec.Y()/TauMuPtVec.X();
	double a = tanphi_tau;
	double b = PV.Y() - a*PV.X();

	double tmu  = (dxy*(cos(phi0) + a*sin(phi0)) - b) / (a*cos(phi0) - sin(phi0)); //projection onto XY plane

	double xdoc = -dxy*sin(phi0) + tmu*cos(phi0);
	double ydoc = dxy*cos(phi0) + tmu*sin(phi0);
	double zdoc = dz + tmu*tan(lam);

	TVector3 PointGuess(xdoc, ydoc, zdoc);
	TVector3 TauMuDir = PointGuess - PV;


	if(rotate){
	  if(fabs(TauMuDir.Phi()-P_Tauh.Phi())<1.5){
		TauMuPtVec.Set(-TauMuPtVec.X(),-TauMuPtVec.Y());
	  }
	}
	TLorentzVector TauMu; TauMu.SetXYZM(TauMuPtVec.X(),TauMuPtVec.Y(),TauMuPtVec.Mod()/tan(TauMuDir.Theta()), PDGInfo::tau_mass());
	TLorentzVector Z(TauMu + Tauh.LV());
	Logger(Logger::Debug) << "Reco Z Px B2B Full Z: " << Z.X() << std::endl;
	Logger(Logger::Debug) << "Reco Z Py B2B Full Z: " << Z.Y() << std::endl;

	return Z;
}
/*
double ZToTaumuTauh::GJErrorCalc(TVector3 SVPV, LorentzVectorParticle A1){
  double error;
  TVector3 A1.
  return error;
}
*/
