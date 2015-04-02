//Ntuple_Controller.cxx IMPLEMENTATION FILE


#include "Ntuple_Controller.h"
#include "Tools.h"
#include "PDG_Var.h"
#include "TF1.h"
#include "Parameters.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"


// External code
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"

///////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::Ntuple_Controller(std::vector<TString> RootFiles):
  copyTree(false)
  ,cannotObtainHiggsMass(false)
  ,ObjEvent(-1)
  ,isInit(false)
{
  // TChains the ROOTuple file
  TChain *chain = new TChain("t");
  Logger(Logger::Verbose) << "Loading " << RootFiles.size() << " files" << std::endl;
  for(unsigned int i=0; i<RootFiles.size(); i++){
    chain->Add(RootFiles[i]);
  }
  TTree *tree = (TTree*)chain;  
  if(chain==0){
	Logger(Logger::Error) << "chain points to NULL" << std::endl;
  }
  Logger(Logger::Info) << "Number of Events in Ntuple: " << chain->GetEntries() << std::endl;
  Ntp=new NtupleReader(tree);
  nbytes=0; 
  nb=0;
  Logger(Logger::Info) << "Ntuple Configured" << std::endl;

  // Fit setup 

  // Resolution uncertainty setup

  gRandom->SetSeed(1234);

  // Rochester muon momentum corrections

  rmcor = new rochcor2012(); // For systematics use rmcor = new rochcor2012(seed!=1234);

  // Set object correction flags to default values
  tauCorrection = "";
  muonCorrection = "";
  elecCorrection = "";
  jetCorrection = "";
}

///////////////////////////////////////////////////////////////////////
//
// Function: void InitEvent()
//
// Purpose: Initialize variables etc on event base
//
///////////////////////////////////////////////////////////////////////

void Ntuple_Controller::InitEvent(){
	Muon_corrected_p4.clear();
	Muon_corrected_p4.resize(NMuons());
	Muon_isCorrected = false;

	// after everything is initialized
	isInit = true;
}

///////////////////////////////////////////////////////////////////////
//
// Function: Int_t Get_Entries()
//
// Purpose: To get the number of events in the Ntuple
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_Entries(){
  return Int_t(Ntp->fChain->GetEntries());
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the event _jentry
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Get_Event(int _jentry){
  jentry=_jentry;
  Ntp->LoadTree(jentry);
  nb = Ntp->fChain->GetEntry(jentry);   nbytes += nb;
  isInit = false;
  InitEvent();
}


///////////////////////////////////////////////////////////////////////
//
// Function: void Get_EventIndex()
//
// Purpose: To get the event index (jentry)
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_EventIndex(){
  return jentry;
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the file name of the root file currently being 
//          accesses
//
///////////////////////////////////////////////////////////////////////
TString Ntuple_Controller::Get_File_Name(){
  return Ntp->fChain->GetCurrentFile()->GetName();
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Branch_Setup(TString B_Name, int type)
//
// Purpose: To setup a branch
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Branch_Setup(TString B_Name, int type){   
  Ntp->fChain->SetBranchStatus(B_Name,type);
}

///////////////////////////////////////////////////////////////////////
//
// destructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::~Ntuple_Controller() {
  Logger(Logger::Verbose) << "Cleaning up" << std::endl;
  delete Ntp;
  delete rmcor;
  Logger(Logger::Verbose) << "Complete." << std::endl;
}


void Ntuple_Controller::CloneTree(TString n){
  if(!copyTree){
	Logger(Logger::Info) << "Starting D3PD cloning" << std::endl;
    newfile = new TFile(n+".root","recreate");
    SkimmedTree=Ntp->fChain->CloneTree(0);
    copyTree=true;
  }
}

void Ntuple_Controller::SaveCloneTree(){
  if(copyTree){
    SkimmedTree->AutoSave();
    newfile->Close();
  }
  Logger(Logger::Info) << "Done"<< std::endl;
}

void Ntuple_Controller::ThinTree(){
  Logger(Logger::Warning) << "ThinTree not implemented." << std::endl;
}

int Ntuple_Controller::SetupSystematics(TString sys){
  return Default;
}


void Ntuple_Controller::ConfigureObjects(){
  if(ObjEvent!=EventNumber()){
    ObjEvent=EventNumber();
    doElectrons();
    doPhotons();
    doJets();
    doMuons();
    doTaus();
    doMET();
  }
}

void Ntuple_Controller::doElectrons(){
  electrons.clear();
  electrons_default.clear();
}

void Ntuple_Controller::doPhotons(){
  photons.clear();
  photons_default.clear();

}

void Ntuple_Controller::doJets(){
  jets.clear();
  jets_default.clear();
}

void Ntuple_Controller::doMuons(){
  muons.clear();
  muons_default.clear();
}

void Ntuple_Controller::doTaus(){
  taus.clear();
  taus_default.clear();
}

void Ntuple_Controller::doMET(){
}


//Physics get Functions
int64_t Ntuple_Controller::GetMCID(){
	int64_t DataMCTypeFromTupel = Ntp->DataMC_Type;

	// move JAK Id information 3 digits to the left
	int64_t jakid = DataMCTypeFromTupel - (DataMCTypeFromTupel%100);
	jakid *= 1000;
	DataMCTypeFromTupel = jakid + (DataMCTypeFromTupel%100);

	// specific to Vladimir's analysis
	if (DataMCTypeFromTupel == DataMCType::DY_ll_Signal && HistoC.hasID(DataMCType::DY_ll_Signal)) {
		for (unsigned int i = 0; i < NMCSignalParticles(); i++) {
			if (abs(MCSignalParticle_pdgid(i)) == PDGInfo::Z0) {
				if (fabs(MCSignalParticle_p4(i).M() - PDG_Var::Z_mass()) < 3 * PDG_Var::Z_width()) {
					return DataMCType::Signal;
				}
			}
		}
		return DataMCTypeFromTupel;
	}

	int dmcType = -999;

	// hack for Higgs mass splitting
	// Higgs mass is added to the MCId, such that the structure is JJJJJJAAABB (with JJJJJJ = JakID, AAA = mass, BB = DataMCType)
	if( (DataMCTypeFromTupel % 100) == DataMCType::H_tautau_ggF ||
		(DataMCTypeFromTupel % 100) == DataMCType::H_tautau_VBF ||
		(DataMCTypeFromTupel % 100) == DataMCType::H_tautau_WHZHTTH){
	  int mass = getSampleHiggsMass();
	  if (mass > 999)	Logger(Logger::Error) << "Read mass with more than 3 digits from sample: m = " << mass << std::endl;
	  if (mass > 0)		DataMCTypeFromTupel += mass*100;
	  // strip off JAK-Id from DataMCType
	  if (HistoC.hasID(DataMCTypeFromTupel % 100000)) {
	  		dmcType = DataMCTypeFromTupel % 100000;
	  	}
	}
	else {
		// strip off JAK-Id from DataMCType
		if (HistoC.hasID(DataMCTypeFromTupel % 100)) {
			dmcType = DataMCTypeFromTupel % 100;
		}
	}

	return dmcType;
}

// return DataMCType without mass information
int Ntuple_Controller::GetStrippedMCID(){
	return GetMCID() % 100;
}

// return path of input dataset (as given in the line InputNtuples in Input.txt)
TString Ntuple_Controller::GetInputNtuplePath(){
	Parameters Par; // assumes configured in Analysis.cxx
	TString dsPath;
	Par.GetString("InputNtuples:",dsPath);
	return dsPath;
}

// return name of input dataset
TString Ntuple_Controller::GetInputDatasetName(){
	TString dsPath = GetInputNtuplePath();
	return gSystem->BaseName( gSystem->DirName( gSystem->DirName(dsPath) ) );
}

//
TString Ntuple_Controller::GetInputPublishDataName(){
	TString dsPath = GetInputNtuplePath();
	return gSystem->BaseName( gSystem->DirName(dsPath) );
}

// determine Higgs mass (from Dataset name or fallback options)
int Ntuple_Controller::getSampleHiggsMass(){
	int mass = -999;

	// default method: analyze dataset name
	mass = readHiggsMassFromString( GetInputNtuplePath() );
	if (mass >= 0) return mass;

	// first fallback: analyze filename (only working when running on GRID)
	mass = readHiggsMassFromString( Get_File_Name() );
	if (mass >= 0) return mass;

	// second fallback: get Higgs mass from MC info
	Logger(Logger::Warning) << "Not able to obtain Higgs mass neither from dataset nor from file name."
	<< "\n\tMake sure the line InputNtuples is set correctly in your Input.txt. "
	<< "\n\tFor now, we will fall back to obtaining the Higgs mass from the generator information."
	<< "\n\tBe aware that SOME EVENTS WILL END UP IN THE WRONG HISTOGRAM!!!" << std::endl;
	return getHiggsSampleMassFromGenInfo();
}

int Ntuple_Controller::readHiggsMassFromString(TString input){
	// loop over possible masses
	for (int m = 100; m < 200; m = m+5){
		if ( input.Contains("M-" + TString::Itoa(m, 10) + "_") ) return m;
	}

	// mass not found
	return -999;
}

double Ntuple_Controller::getResonanceMassFromGenInfo(bool useZ0 /* = true*/, bool useHiggs0 /* = true*/, bool useW /* = true*/){
	for (unsigned int i = 0; i < NMCSignalParticles(); i++) {
		unsigned pdgid = abs(MCSignalParticle_pdgid(i));
		bool Z0 = useZ0 && (pdgid == PDGInfo::Z0);
		bool H0 = useHiggs0 && (pdgid == PDGInfo::Higgs0);
		bool W = useW && (pdgid == PDGInfo::W_plus);
		if (Z0 || H0 || W) {
			return MCSignalParticle_p4(i).M();
		}
	}
	return -999;
}

int Ntuple_Controller::getHiggsSampleMassFromGenInfo(){
	double resMass = getResonanceMassFromGenInfo(false, true, false);
	for (int m = 100; m < 200; m = m+5){
		if (fabs(resMass - m) < 2.5 ) {
			return m;
		}
	}
	return -999;
}

TMatrixF Ntuple_Controller::Vtx_Cov(unsigned int i){
  unsigned int dim=3;
  TMatrixF M(dim,dim);
  for(unsigned int j=0;j<dim;j++){
    for(unsigned int k=0;k<=j;k++){
      M[j][k]=Ntp->Vtx_Cov->at(i).at(j).at(k);
      M[k][j]=Ntp->Vtx_Cov->at(i).at(j).at(k);
    }
  }
  return M;
}

bool Ntuple_Controller::isVtxGood(unsigned int i){
  if(0<=i && i<NVtx()){
    if(Vtx_Track_idx(i).size()>4)return true;
  }
  return false;
}

bool Ntuple_Controller::isGoodVtx(unsigned int i){
	if(fabs(Vtx(i).z())>=24) return false;
	if(Vtx(i).Perp()>=2) return false;
	if(Vtx_ndof(i)<=4) return false;
	if(Vtx_isFake(i)) return false;
	return true;
}

bool Ntuple_Controller::isGoodMuon(unsigned int i){
  //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger 
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel       
  //  isGoodMuon_nooverlapremoval(i) with
  //  ΔR(μ,jet)>0.3 where jet is any jet passing the jet requirements not applied applied       
  if(isGoodMuon_nooverlapremoval(i)){
    unsigned int jet_idx=0;
    return !muonhasJetOverlap(i,jet_idx);
  }
  return false;
}

bool Ntuple_Controller::muonhasJetOverlap(unsigned int muon_idx,unsigned int &jet_idx){
  for(unsigned int j=0;j<NPFJets();j++){
    if(isGoodJet_nooverlapremoval(j)){
      if(Tools::dr(Muon_p4(muon_idx),PFJet_p4(j))>0.2 && Tools::dr(Muon_p4(muon_idx),PFJet_p4(j))<0.4){ jet_idx=j;return true;}
    }
  }
  return false;
}

bool Ntuple_Controller::muonhasJetMatch(unsigned int muon_idx,unsigned int &jet_idx){
  for(unsigned int j=0;j<NPFJets();j++){
    if(isGoodJet_nooverlapremoval(j)){
      if(Tools::dr(Muon_p4(muon_idx),PFJet_p4(j))<0.2){ jet_idx=j;return true;}
    }
  }
  return false;
}



bool Ntuple_Controller::isGoodMuon_nooverlapremoval(unsigned int i){
  //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger and no overlpar removal applied
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel  
  //  GlobalMuon && TrackerMuon applied applied                                         
  //  pT > 20 GeV                                                                       
  //  fasb(eta) < 2.4                                                                   
  //  normChi2 < 10                                                                     
  //  n Tracker Hits > 10                                                               
  //  n Muon Hit > 0                                                                    
  //  Transverse IP of the muon wrt the beam spot (cm) < 0.02                           
  //  relIso < 0.20                                                                      
  //  fabs( muon.vertex().z() - PV.z()) not applied                           
  //  muon.innerTrack()->hitPattern().pixelLayersWithMeasurement() not applied
  //  numberOfMatchedStations() not applied                              
  if(Muon_isGlobalMuon(i) && Muon_isStandAloneMuon(i)){
    if(Muon_p4(i).Pt()>15.0){
      if(fabs(Muon_p4(i).Eta())<2.4){
	if(Muon_normChi2(i)<10.0){
	  if(Muon_innerTrack_numberofValidHits(i)>10){
	    if(Muon_hitPattern_numberOfValidMuonHits(i)>0){
	      return true;
	    }
	  }
	}
      }
    }
  }
  return false;
}

void Ntuple_Controller::CorrectMuonP4(){
	if(isInit){
		for(unsigned int i=0;i<NMuons();i++){
			TLorentzVector mup4 = Muon_p4(i,"");
			int runopt = 0; // 0: no run-dependece
			float qter = 1.0; // 1.0: don't care about muon momentum uncertainty
			if(!isData() && GetStrippedMCID()!=DataMCType::DY_emu_embedded && GetStrippedMCID()!=DataMCType::DY_mutau_embedded){
				rmcor->momcor_mc(mup4,Muon_Charge(i),runopt,qter);
			}else{
				rmcor->momcor_data(mup4,Muon_Charge(i),runopt,qter);
			}
			Muon_corrected_p4.at(i) = mup4;
		}
		Muon_isCorrected = true;
	}else{
		Muon_isCorrected = false;
		Logger(Logger::Warning) << "No muon corrections applied" << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//
// Get muon four-vector
//
// Options:
//  - "roch": will correct the momentum in data and MC using Rochester corrections. For systematics look in constructor of this class.
//  - "scale": if you don't use momentum corrections, use this to estimate the systematics on the scale (only MC).
//             Add "down" in order to estimate it for downward variations.
//  - "res": if you don't use momentum corrections, use this to estimate systematics caused by momentum resolution (only MC)
//

TLorentzVector Ntuple_Controller::Muon_p4(unsigned int i, TString corr){
	TLorentzVector vec = TLorentzVector(Ntp->Muon_p4->at(i).at(1),Ntp->Muon_p4->at(i).at(2),Ntp->Muon_p4->at(i).at(3),Ntp->Muon_p4->at(i).at(0));
	if (corr == "default") corr = muonCorrection;
	if(corr.Contains("roch")){
		if(!Muon_isCorrected){
			CorrectMuonP4();
		}
		if(Muon_isCorrected){
			vec = Muon_corrected_p4.at(i);
		}
		else{
			Logger(Logger::Warning) << "No muon corrections applied" << std::endl;
		}
	}
	if(!isData() && GetStrippedMCID()!=DataMCType::DY_emu_embedded && GetStrippedMCID()!=DataMCType::DY_mutau_embedded){
		if(corr.Contains("scale")){
			if(!corr.Contains("down")) vec.SetPerp(vec.Perp()*1.002);
			else vec.SetPerp(vec.Perp()*0.998);
		}else if(corr.Contains("res")){
			vec.SetPerp(gRandom->Gaus(vec.Perp(),vec.Perp()*0.006));
		}
		if(corr.Contains("met")){
			if(!corr.Contains("down")) vec.SetPerp(vec.Perp() * 1.002);
			else vec.SetPerp(vec.Perp() * 0.998);
		}
	}
	return vec;
}

/////////////////////////////////////////////////////////////////////
//
// Official muon id code
//

bool Ntuple_Controller::isTightMuon(unsigned int i){
	if(!Muon_isGlobalMuon(i)) return false;
	if(!Muon_isPFMuon(i)) return false;
	if(Muon_normChi2(i)>=10.) return false;
	if(Muon_hitPattern_numberOfValidMuonHits(i)<=0) return false;
	if(Muon_numberOfMatchedStations(i)<=1) return false;
	if(Muon_numberofValidPixelHits(i)<=0) return false;
	if(Muon_trackerLayersWithMeasurement(i)<=5) return false;
	return true;
}

bool Ntuple_Controller::isTightMuon(unsigned int i, unsigned int j, TString corr){
	if(j<0 || j>=NVtx()) return false;
	if(!isTightMuon(i)) return false;
	if(dxy(Muon_p4(i,corr),Muon_Poca(i),Vtx(j))>=0.2) return false;
	if(dz(Muon_p4(i,corr),Muon_Poca(i),Vtx(j))>=0.5) return false;
	return true;
}

bool Ntuple_Controller::isLooseMuon(unsigned int i){
	if(!Muon_isPFMuon(i)) return false;
	if(!(Muon_isGlobalMuon(i) || Muon_isTrackerMuon(i))) return false;
	return true;
}

float Ntuple_Controller::Muon_RelIso(unsigned int i, TString corr){
	return (Muon_sumChargedHadronPt04(i)+std::max(0.,Muon_sumNeutralHadronEt04(i)+Muon_sumPhotonEt04(i)-0.5*Muon_sumPUPt04(i)))/Muon_p4(i,corr).Pt();
}

/////////////////////////////////////////////////////////////////////

bool Ntuple_Controller::isSelectedMuon(unsigned int i, unsigned int j, double impact_xy, double impact_z, TString corr){
	if(j<0 || j>=NVtx()) return false;
	if(!isTightMuon(i)) return false;
	if(dxy(Muon_p4(i,corr),Muon_Poca(i),Vtx(j))>=impact_xy) return false;
	if(dz(Muon_p4(i,corr),Muon_Poca(i),Vtx(j))>=impact_z) return false;
	return true;
}

/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//
// Get electron four-vector
//
// Options:
//  - "res": use this to estimate the impact of the electron energy resolution on your result (only MC).
//           resolution in barrel (endcaps) for 2012: 1.6% (4.1%). uncertainty on resolution: 10%
//           RegEnergy corresponds to the electron energy after regression but before scale correction (data)
//           or smearing (MC)
//

TLorentzVector Ntuple_Controller::Electron_p4(unsigned int i, TString corr){
	TLorentzVector vec = TLorentzVector(Ntp->Electron_p4->at(i).at(1),Ntp->Electron_p4->at(i).at(2),Ntp->Electron_p4->at(i).at(3),Ntp->Electron_p4->at(i).at(0));
	if(!isData() && GetStrippedMCID()!=DataMCType::DY_emu_embedded && GetStrippedMCID()!=DataMCType::DY_mutau_embedded){
		if (corr == "default") corr = elecCorrection;
		if(corr.Contains("scale") && Electron_RegEnergy(i)!=0){
			if(!corr.Contains("down")) vec.SetPerp(vec.Perp() * (1+Electron_RegEnergyError(i)/Electron_RegEnergy(i)));
			else vec.SetPerp(vec.Perp() * (1-Electron_RegEnergyError(i)/Electron_RegEnergy(i)));
		}
		if(corr.Contains("res")){
			if(Electron_RegEnergy(i)>0){
				if(fabs(Electron_supercluster_eta(i))<1.479){
					if(corr.Contains("down")) vec.SetPerp(vec.Perp() * gRandom->Gaus(Electron_RegEnergy(i),Electron_RegEnergy(i)*0.0144) / Electron_RegEnergy(i));
					else vec.SetPerp(gRandom->Gaus(vec.Perp() * Electron_RegEnergy(i),Electron_RegEnergy(i)*0.0176) / Electron_RegEnergy(i));
				}
				else if(fabs(Electron_supercluster_eta(i))<2.5){
					if(corr.Contains("down")) vec.SetPerp(vec.Perp() * gRandom->Gaus(Electron_RegEnergy(i),Electron_RegEnergy(i)*0.0369) / Electron_RegEnergy(i));
					else vec.SetPerp(gRandom->Gaus(vec.Perp() * Electron_RegEnergy(i),Electron_RegEnergy(i)*0.0451) / Electron_RegEnergy(i));
				}
				else{
					Logger(Logger::Warning) << "Eta out of range: " << Electron_supercluster_eta(i) << ". Returning fourvector w/o smearing for resolution uncertainties." << std::endl;
				}
			}
			else{
				Logger(Logger::Warning) << "Energy <= 0: " << Electron_RegEnergy(i) << ". Returning fourvector w/o smearing for resolution uncertainties." << std::endl;
			}
		}
		if(corr.Contains("met")){
			if(fabs(Electron_supercluster_eta(i))<1.479){
				if(!corr.Contains("down")) vec.SetPerp(vec.Perp() * 1.006);
				else vec.SetPerp(vec.Perp() * 0.994);
			}
			else if(fabs(Electron_supercluster_eta(i))<2.5){
				if(!corr.Contains("down")) vec.SetPerp(vec.Perp() * 1.015);
				else vec.SetPerp(vec.Perp() * 0.985);
			}
			else{
				Logger(Logger::Warning) << "Eta out of range: " << Electron_supercluster_eta(i) << ". Returning fourvector w/o scale corrections for MET uncertainties." << std::endl;
			}
		}
	}
	return vec;
}

/////////////////////////////////////////////////////////////////////
//
// Official electron id code
//

bool Ntuple_Controller::isTrigPreselElectron(unsigned int i){
	if(fabs(Electron_supercluster_eta(i))>2.5) return false;
	if(Electron_numberOfMissedHits(i)>0) return false;
	if(Electron_Gsf_dr03TkSumPt(i)/Electron_p4(i).Pt()>0.2) return false;
	if(Electron_Gsf_dr03EcalRecHitSumE(i)/Electron_p4(i).Pt()>0.2) return false;
	if(Electron_Gsf_dr03HcalTowerSumEt(i)/Electron_p4(i).Pt()>0.2) return false;
	if(fabs(Electron_supercluster_eta(i))<1.479){
		if(Electron_sigmaIetaIeta(i)>0.014) return false;
		if(Electron_hadronicOverEm(i)>0.15) return false;
	}else{
		if(Electron_sigmaIetaIeta(i)>0.035) return false;
		if(Electron_hadronicOverEm(i)>0.1) return false;
	}
	return true;
}

bool Ntuple_Controller::isMVATrigElectron(unsigned int i, TString corr){
	// !!! make sure to also apply Electron_RelIso<0.15 in your analysis !!!
	double mvapt = Electron_p4(i,corr).Pt();
	double mvaeta = fabs(Electron_supercluster_eta(i));
	if(mvapt<10.) return false;
	if(mvaeta>2.5) return false;
	if(Electron_numberOfMissedHits(i)>0) return false;
	if(Electron_HasMatchedConversions(i)) return false;
	if(!isTrigPreselElectron(i)) return false;
	if(mvapt>10. && mvapt<20.){
		if(mvaeta<0.8 && Electron_MVA_Trig_discriminator(i)<=0.00) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_Trig_discriminator(i)<=0.10) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_Trig_discriminator(i)<=0.62) return false;
	}
	if(mvapt>=20.){
		if(mvaeta<0.8 && Electron_MVA_Trig_discriminator(i)<=0.94) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_Trig_discriminator(i)<=0.85) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_Trig_discriminator(i)<=0.92) return false;
	}
	return true;
}

bool Ntuple_Controller::isTrigNoIPPreselElectron(unsigned int i){
	if(fabs(Electron_supercluster_eta(i))>2.5) return false;
	if(Electron_numberOfMissedHits(i)>0) return false;
	if(Electron_Gsf_dr03TkSumPt(i)/Electron_p4(i).Pt()>0.2) return false;
	if(Electron_Gsf_dr03HcalTowerSumEt(i)/Electron_p4(i).Pt()>0.2) return false;
	if(fabs(Electron_supercluster_eta(i))<1.479){
		if(Electron_sigmaIetaIeta(i)>0.01) return false;
		if(Electron_hadronicOverEm(i)>0.12) return false;
		if(fabs(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i))>0.007) return false;
		if(fabs(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i))>0.15) return false;
	}else{
		if(Electron_sigmaIetaIeta(i)>0.03) return false;
		if(Electron_hadronicOverEm(i)>0.1) return false;
		if(fabs(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i))>0.009) return false;
		if(fabs(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i))>0.1) return false;
	}
	return true;
}

bool Ntuple_Controller::isMVATrigNoIPElectron(unsigned int i, TString corr){
	// at present there are no recommendations on the isolation
	double mvapt = Electron_p4(i, corr).Pt();
	double mvaeta = fabs(Electron_supercluster_eta(i));
	if(mvaeta>2.5) return false;
	if(Electron_HasMatchedConversions(i)) return false;
	if(Electron_numberOfMissedHits(i)>0) return false;
	if(!isTrigNoIPPreselElectron(i)) return false;
	if(mvapt<20){
		if(mvaeta<0.8 && Electron_MVA_TrigNoIP_discriminator(i)<=-0.5375) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_TrigNoIP_discriminator(i)<=-0.375) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_TrigNoIP_discriminator(i)<=-0.025) return false;
	}
	if(mvapt>=20){
		if(mvaeta<0.8 && Electron_MVA_TrigNoIP_discriminator(i)<=0.325) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_TrigNoIP_discriminator(i)<=0.775) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_TrigNoIP_discriminator(i)<=0.775) return false;
	}
	return true;
}

bool Ntuple_Controller::isMVANonTrigElectron(unsigned int i, unsigned int j, TString corr){
	// !!! make sure to also apply Electron_RelIso<0.4 in your analysis !!!
	double mvapt = Electron_p4(i,corr).Pt();
	double mvaeta = fabs(Electron_supercluster_eta(i));
	if(mvapt<7.) return false;
	if(mvaeta>2.5) return false;
	if(Electron_numberOfMissedHits(i)>1) return false;
	if(vertexSignificance(Electron_Poca(i),j)>=4) return false;
	if(mvapt>7. && mvapt<10.){
		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=0.47) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=0.004) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.295) return false;
	}
	if(mvapt>=10.){
		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=-0.34) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=-0.65) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.6) return false;
	}
	return true;
}

bool Ntuple_Controller::isTightElectron(unsigned int i, TString corr){
	if(Electron_HasMatchedConversions(i)) return false;
	if(Electron_numberOfMissedHits(i)>0) return false;
	if(Electron_RelIso04(i,corr)>=0.1) return false;
	if(fabs(1/Electron_ecalEnergy(i)-1/Electron_trackMomentumAtVtx(i))>=0.05) return false;
	if(fabs(Electron_supercluster_eta(i))<=1.479){
		if(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.004) return false;
		if(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.03) return false;
		if(Electron_sigmaIetaIeta(i)>=0.01) return false;
		if(Electron_hadronicOverEm(i)>=0.12) return false;
	}
	if(fabs(Electron_supercluster_eta(i))>1.479 && fabs(Electron_supercluster_eta(i))<2.5){
		if(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.005) return false;
		if(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.02) return false;
		if(Electron_sigmaIetaIeta(i)>=0.03) return false;
		if(Electron_hadronicOverEm(i)>=0.10) return false;
		if(Electron_p4(i).Pt()<20 && Electron_RelIso04(i)>=0.07) return false;
	}
	return true;
}

bool Ntuple_Controller::isTightElectron(unsigned int i, unsigned int j, TString corr){
	if(j<0 || j>=NVtx()) return false;
	if(!isTightElectron(i,corr)) return false;
	if(dxy(Electron_p4(i,corr),Electron_Poca(i),Vtx(j))>=0.02) return false;
	if(dz(Electron_p4(i,corr),Electron_Poca(i),Vtx(j))>=0.1) return false;
	return true;
}

float Ntuple_Controller::Electron_RelIso03(unsigned int i, TString corr){
	return (Electron_chargedHadronIso(i)+std::max((double)0.,Electron_neutralHadronIso(i)+Electron_photonIso(i)-RhoIsolationAllInputTags()*Electron_Aeff_R03(Electron_supercluster_eta(i))))/Electron_p4(i,corr).Pt();
}

double Ntuple_Controller::Electron_RelIsoDep03(unsigned int i, TString corr){
	return (Electron_isoDeposits_chargedHadronIso03(i)+std::max((double)0.,Electron_isoDeposits_neutralHadronIso03(i)+Electron_isoDeposits_photonIso03(i)-RhoIsolationAllInputTags()*Electron_Aeff_R03(Electron_supercluster_eta(i))))/Electron_p4(i,corr).Pt();
}

float Ntuple_Controller::Electron_RelIso04(unsigned int i, TString corr){
	return (Electron_chargedHadronIso(i)+std::max((double)0.,Electron_neutralHadronIso(i)+Electron_photonIso(i)-RhoIsolationAllInputTags()*Electron_Aeff_R04(Electron_supercluster_eta(i))))/Electron_p4(i,corr).Pt();
}

double Ntuple_Controller::Electron_RelIsoDep04(unsigned int i, TString corr){
	return (Electron_isoDeposits_chargedHadronIso04(i)+std::max((double)0.,Electron_isoDeposits_neutralHadronIso04(i)+Electron_isoDeposits_photonIso04(i)-RhoIsolationAllInputTags()*Electron_Aeff_R04(Electron_supercluster_eta(i))))/Electron_p4(i,corr).Pt();
}

double Ntuple_Controller::Electron_Aeff_R04(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.208;
	else if(eta>=1. && eta<1.479) return 0.209;
	else if(eta>=1.479 && eta<2.) return 0.115;
	else if(eta>=2. && eta<2.2) return 0.143;
	else if(eta>=2.2 && eta<2.3) return 0.183;
	else if(eta>=2.3 && eta<2.4) return 0.194;
	else if(eta>=2.4) return 0.261;
	else {Logger(Logger::Error) << "Electron eta out of range: " << Eta << std::endl; return -1;}
}

double Ntuple_Controller::Electron_Aeff_R03(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.13;
	else if(eta>=1. && eta<1.479) return 0.14;
	else if(eta>=1.479 && eta<2.) return 0.07;
	else if(eta>=2. && eta<2.2) return 0.09;
	else if(eta>=2.2 && eta<2.3) return 0.11;
	else if(eta>=2.3 && eta<2.4) return 0.11;
	else if(eta>=2.4) return 0.14;
	else {Logger(Logger::Error) << "Electron eta out of range: " << Eta << std::endl; return -1;}
}

/////////////////////////////////////////////////////////////////////

bool Ntuple_Controller::isSelectedElectron(unsigned int i, unsigned int j, double impact_xy, double impact_z, TString corr){
	double mvapt = Electron_p4(i,corr).Pt();
	double mvaeta = fabs(Electron_supercluster_eta(i));
	if(Electron_numberOfMissedHits(i)>0) return false;
	if(Electron_HasMatchedConversions(i)) return false;
	if(dxy(Electron_p4(i,corr),Electron_Poca(i),Vtx(j))>=impact_xy) return false;
	if(dz(Electron_p4(i,corr),Electron_Poca(i),Vtx(j))>=impact_z) return false;
	if(mvapt<20.){
		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=0.925) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=0.915) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.965) return false;
	}
	if(mvapt>=20.){
		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=0.905) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=0.955) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.975) return false;
	}
	return true;
}

bool Ntuple_Controller::isGoodJet(unsigned int i){
  //  Top Dilepton Jet selection with pt 15GeV
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel 
  //  isGoodJet_nooverlapremoval(i) with:
  //  deltaR jet-electron cleaning < 0.4 (2 selected lepton only) < 0.3 (e+jets only) 0.3
  //  deltaR jet-muon cleaning < 0.4(2 selected lepton only) < 0.3 (mu+jets only), 0.1 for PF and JET 0.3
  if(isGoodJet_nooverlapremoval(i)){
    unsigned int muon_idx;
    return !jethasMuonOverlap(i,muon_idx);
  }
  return false;
}

bool Ntuple_Controller::jethasMuonOverlap(unsigned int jet_idx,unsigned int &muon_idx){
  for(unsigned int j=0;j<NMuons();j++){
    if(isGoodMuon_nooverlapremoval(j) && Muon_RelIso(j)<0.2){
      if(Tools::dr(Muon_p4(j),PFJet_p4(jet_idx))<0.4){ muon_idx=j;return true;}
    }
  }
  return false;
}


bool Ntuple_Controller::isGoodJet_nooverlapremoval(unsigned int i){
  //  Top Dilepton Jet selection with pt 15GeV
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel
  //  Jet ID defined in isJetID(i)               
  //  cut         dilepton    l+jet
  //  corrected pT 30 GeV 30 GeV    
  //  residual correction (data) applied applied
  //  abs(eta) < 2.5 < 2.4                      
  //  jet ID applied applied     
  if(isJetID(i)){
    if(PFJet_p4(i).Pt()>15.0){
      if(fabs(PFJet_p4(i).Eta())<2.4){
	return true;
      }
    }
  }
  return false;
}

bool Ntuple_Controller::isJetID(unsigned int i, TString corr){
  //  Top Dilepton Jet selection with pt and iso matching the muon and tau.
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel  
  //  Jet ID :
  //  corrected jet pt > 10 GeV, jet eta < 5.2
  //  number of constituents>1 (patJet->numberOfDaughters())
  //  NHF<0.99 (( patJet->neutralHadronEnergy() + patJet->HFHadronEnergy() ) / patJet->energy())
  //  NEF<0.99 (patJet->neutralEmEnergyFraction())
  //  if |η|<2.4, CEF<0.99 (patJet->chargedEmEnergyFraction())
  //  if |η|<2.4, CHF>0 (patJet->chargedHadronEnergyFraction())
  //  if |η|<2.4, NCH>0 (patJet->chargedMultiplicity()) 
  /////////////////////////////////////////////////////////////////////////
  // apply jet ID
  if(PFJet_p4(i,corr).Pt()<=10.) return false;
  if(fabs(PFJet_p4(i,corr).Eta())>=5.2) return false;
  if(PFJet_numberOfDaughters(i)<=1) return false;
  if((PFJet_neutralHadronEnergy(i)+PFJet_HFHadronEnergy(i))/PFJet_p4(i,corr).E()>=0.99) return false;
  if(PFJet_neutralEmEnergyFraction(i)>=0.99) return false;
  if(fabs(PFJet_p4(i,corr).Eta())<2.4){
	  if(PFJet_chargedEmEnergyFraction(i)>=0.99) return false;
	  if(PFJet_chargedHadronEnergyFraction(i)<=0.) return false;
	  if(PFJet_chargedMultiplicity(i)<=0) return false;
  }
  return true;
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECL2ResidualTimeStability#2012Rereco
double Ntuple_Controller::rundependentJetPtCorrection(double jeteta, int runnumber){
	if(!isData() && GetStrippedMCID()!=DataMCType::DY_emu_embedded && GetStrippedMCID()!=DataMCType::DY_mutau_embedded){
		return 1.;
	}
	const double corrs[5] = {0.0, -0.454e-6, -0.952e-6, 1.378e-6, 0.0};
	const int run0 = 201000;
	double eta = fabs(jeteta);
	double corr = 0.;
	if(eta<1.3) corr = corrs[0];
	else if(eta<2.0) corr = corrs[1];
	else if(eta<2.5) corr = corrs[2];
	else if(eta<3.0) corr = corrs[3];
	else if(eta<5.0) corr = corrs[4];
	return (1.+corr*(runnumber-run0));
}

double Ntuple_Controller::JERCorrection(TLorentzVector jet, double dr, TString corr){
	double sf = jet.Pt();
	if (corr == "default") corr = jetCorrection;
	if(isData() || GetStrippedMCID()==DataMCType::DY_emu_embedded || GetStrippedMCID()==DataMCType::DY_mutau_embedded
			|| jet.Pt()<=10
			|| PFJet_matchGenJet(jet,dr)==TLorentzVector(0.,0.,0.,0.)
			){
		return sf;
	}else{
		double c = JetEnergyResolutionCorr(jet.Eta());
		if(corr.Contains("up")) c += JetEnergyResolutionCorrErr(jet.Eta());
		if(corr.Contains("down")) c -= JetEnergyResolutionCorrErr(jet.Eta());
		sf = std::max(0.,c*jet.Pt()+(1.-c)*PFJet_matchGenJet(jet,dr).Pt());
	}
	return sf;
}

TLorentzVector Ntuple_Controller::PFJet_matchGenJet(TLorentzVector jet, double dr){
	TLorentzVector genjet(0.,0.,0.,0.);
	for(unsigned int i=0;i<PFJet_NGenJetsNoNu();i++){
		if(PFJet_GenJetNoNu_p4(i).Vect()!=TVector3(0.,0.,0.)
				&& jet.DeltaR(PFJet_GenJetNoNu_p4(i))<dr
				){
			genjet = PFJet_GenJetNoNu_p4(i);
		}
	}
	return genjet;
}

double Ntuple_Controller::JetEnergyResolutionCorr(double jeteta){
	double eta = fabs(jeteta);
	double corr = 1.;
	if(eta<0.5) corr = 1.079;
	else if(eta<1.1) corr = 1.099;
	else if(eta<1.7) corr = 1.121;
	else if(eta<2.3) corr = 1.208;
	else if(eta<2.8) corr = 1.254;
	else if(eta<3.2) corr = 1.395;
	else if(eta<5.0) corr = 1.056;
	return corr;
}

double Ntuple_Controller::JetEnergyResolutionCorrErr(double jeteta){
	double eta = fabs(jeteta);
	double err = 0.;
	if(eta<0.5) err = 0.026;
	else if(eta<1.1) err = 0.027;
	else if(eta<1.7) err = 0.029;
	else if(eta<2.3) err = 0.046;
	else if(eta<2.8) err = 0.062;
	else if(eta<3.2) err = 0.063;
	else if(eta<5.0) err = 0.191;
	return err;
}

/////////////////////////////////////////////////////////////////////
//
// Get jet four-vector
//
// Options:
//  - "run": corrects the jet pt to account for calorimeter degredation during data taking (only data).
//  - "JER": smears the jet pt in MC to match the resolution in data. Additional use of "up" or "down"
//           varies the correction by its uncertainty -> systematics
//  - "JEC": use this to estimate the impact of scale correction uncertainties on your result.
//           use "plus" for an upward variation. if you use nothing, the variation will be downward.
//

TLorentzVector Ntuple_Controller::PFJet_p4(unsigned int i, TString corr){
	TLorentzVector vec = TLorentzVector(Ntp->PFJet_p4->at(i).at(1),Ntp->PFJet_p4->at(i).at(2),Ntp->PFJet_p4->at(i).at(3),Ntp->PFJet_p4->at(i).at(0));
	if (corr == "default") corr = jetCorrection;
	// apply run-dependent pT corrections
	if (corr.Contains("run")){
		vec.SetPerp(vec.Pt() * rundependentJetPtCorrection(vec.Eta(), RunNumber()));
	}
	if(corr.Contains("JER")){
		vec.SetPerp(JERCorrection(vec,0.25,corr));
	}
	if(corr.Contains("JEC")){
		if(corr.Contains("plus")) vec.SetPerp(vec.Pt() * (1 + PFJet_JECuncertainty(i)));
		else vec.SetPerp(vec.Pt() * (1 - PFJet_JECuncertainty(i)));
	}
	return vec;
}

double Ntuple_Controller::TauSpinerGet(int SpinType){
#ifdef USE_TauSpinner
  if(!isData()){
    TauDecay taudecay;
    std::vector<SimpleParticle> tau_daughters, tau_daughters2;
    SimpleParticle tau, tau2;
    for(int i=0; i<NMCSignalParticles();i++){
      // check for signal Boson
      if(MCSignalParticle_pdgid(i)==25 || 
	 MCSignalParticle_pdgid(i)==36 || 
	 MCSignalParticle_pdgid(i)==22 || 
	 MCSignalParticle_pdgid(i)==23){
	if(MCSignalParticle_Tauidx(i).size()==2){
	  SimpleParticle X(MCSignalParticle_p4(i).Px(),MCSignalParticle_p4(i).Py(),MCSignalParticle_p4(i).Pz(),
			   MCSignalParticle_p4(i).E(),MCSignalParticle_pdgid(i));
	  bool tau1good(false),tau2good(false);
	  //first tau
	  unsigned int tauidx=MCSignalParticle_Tauidx(i).at(0);
	  Logger(Logger::Verbose) << "tau 1 indx " << tauidx  << " Number of Tau and Products " << NMCTauDecayProducts(tauidx) << std::endl;
	  if(tauidx<NMCTaus()){
	    for(int t=0;t<NMCTauDecayProducts(tauidx);t++){
	      int mypdgid=abs((int)MCTauandProd_pdgid(tauidx,t));
	      if(abs( mypdgid)==abs(PDGInfo::tau_plus) ){
		tau=SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),MCTauandProd_p4(tauidx,t).Py(),MCTauandProd_p4(tauidx,t).Pz(),
				   MCTauandProd_p4(tauidx,t).E(),MCTauandProd_pdgid(tauidx,t));
	      }
	      else if((taudecay.isTauFinalStateParticle(mypdgid) && mypdgid!=PDGInfo::gamma)){
		tau_daughters.push_back(SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),MCTauandProd_p4(tauidx,t).Py(),
						       MCTauandProd_p4(tauidx,t).Pz(),MCTauandProd_p4(tauidx,t).E(),
						       MCTauandProd_pdgid(tauidx,t)));
		tau1good=true;
	      }
	    }
	  }
	  // second tau
	  tauidx=MCSignalParticle_Tauidx(i).at(1);
	  Logger(Logger::Verbose) << "tau 2 indx " << tauidx  << " Number of Tau and Products " << NMCTauDecayProducts(tauidx) << std::endl;
	  if(tauidx<NMCTaus()){
	    for(int t=0;t<NMCTauDecayProducts(tauidx);t++){
	      int mypdgid=abs((int)MCTauandProd_pdgid(tauidx,t));
	      if(abs( mypdgid)==abs(PDGInfo::tau_plus)){
		    tau2=SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),
					MCTauandProd_p4(tauidx,t).Py(),
					MCTauandProd_p4(tauidx,t).Pz(),
					MCTauandProd_p4(tauidx,t).E(),
					MCTauandProd_pdgid(tauidx,t));
		    tau2good=true;
	      }
	      if((taudecay.isTauFinalStateParticle(mypdgid) && mypdgid!=PDGInfo::gamma)){
	    Logger(Logger::Verbose) << "isDaughter" << std::endl;
		if(tau_daughters.size()>0)tau2good=true;
		tau_daughters2.push_back(SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),
							MCTauandProd_p4(tauidx,t).Py(),
							MCTauandProd_p4(tauidx,t).Pz(),
							MCTauandProd_p4(tauidx,t).E(),
							MCTauandProd_pdgid(tauidx,t)));
	      }
	    }
	  }
	  if(tau1good && tau2good){
		Logger(Logger::Verbose)  << "Two Taus found: " << tau_daughters.size() << " " << tau_daughters2.size() << std::endl;
	    return TauSpinerInt.Get(SpinType,X,tau,tau_daughters,tau2,tau_daughters2);
	  }
	}
      }
    }
  }
#endif
  return 1.0;
}

/////////////////////////////////////////////////////////////////////
//
// Get tau four-vector
//
// Options:
//  - "scalecorr": corrects the tau energy scale depending on the decay mode (only MC and embedding).
//

TLorentzVector Ntuple_Controller::PFTau_p4(unsigned int i, TString corr){
	TLorentzVector vec = TLorentzVector(Ntp->PFTau_p4->at(i).at(1),Ntp->PFTau_p4->at(i).at(2),Ntp->PFTau_p4->at(i).at(3),Ntp->PFTau_p4->at(i).at(0));
	if (corr == "default") corr = tauCorrection;
	if(!isData() || GetStrippedMCID() == DataMCType::DY_mutau_embedded){
		if(corr.Contains("scalecorr")){
			if(PFTau_hpsDecayMode(i)>0 && PFTau_hpsDecayMode(i)<5){
				vec *= 1.025+0.001*min(max(vec.Pt()-45.,0.),10.);
			}
			else if(PFTau_hpsDecayMode(i)>=10){
				vec *= 1.012+0.001*min(max(vec.Pt()-32.,0.),18.);
			}
		}
		if(corr.Contains("met")){
			if(!corr.Contains("down")) vec.SetPerp(vec.Perp() * 1.03);
			else vec.SetPerp(vec.Perp() * 0.97);
		}
	}
	return vec;
}


bool Ntuple_Controller::hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &tau_idx){
  for(unsigned int i=0; i<NMCSignalParticles();i++){
    if(MCSignalParticle_pdgid(i)==(int)parent_pdgid){
      for(unsigned int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
	if((int)MCSignalParticle_Tauidx(i).at(j)>=NMCTaus()){
	  Logger(Logger::Warning) << "INVALID Tau index... Skipping event! MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << std::endl;
	  return false;
	}
      }
      for(unsigned int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
	unsigned int tauidx=MCSignalParticle_Tauidx(i).at(j);
	Logger(Logger::Verbose) << "MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << " " << Ntp->MCTau_JAK->size() << " " << Ntp->MCTauandProd_pdgid->size() << std::endl;
	if((int)MCTau_JAK(tauidx)==tau_jak){ tau_idx=tauidx;Boson_idx=i;return true;}
      }
    }
  }
  return false;
}

bool Ntuple_Controller::hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,unsigned int &tau1_idx, unsigned int &tau2_idx){
  for(unsigned int i=0; i<NMCSignalParticles();i++){
    if(MCSignalParticle_pdgid(i)==parent_pdgid){
      for(unsigned int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
        if((int)MCSignalParticle_Tauidx(i).at(j)>=NMCTaus()){
          Logger(Logger::Warning) << "INVALID Tau index... Skipping event! MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << std::endl;
          return false;
        }
      }
      if(MCSignalParticle_Tauidx(i).size()==2){
	tau1_idx=MCSignalParticle_Tauidx(i).at(0);
        tau2_idx=MCSignalParticle_Tauidx(i).at(1);
	Boson_idx=i;
	return true;
      }
    }
  }
  return false;
}

// calculate flight length significance from primary and secondary vertex info
double Ntuple_Controller::PFTau_FlightLenght_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov ){
  TVector3 SVPV = sv - pv;
  TVectorF FD;
  FD.ResizeTo(3);
  FD(0) = SVPV.X();
  FD(1) = SVPV.Y();
  FD(2) = SVPV.Z();

  TMatrixT<double> PVcv;
  PVcv.ResizeTo(3,3);
  for(int nr =0; nr<PVcov.GetNrows(); nr++){
    for(int nc =0; nc<PVcov.GetNcols(); nc++){
      PVcv(nr,nc) = PVcov(nr,nc);
    }
  }
  TMatrixT<double> SVcv;
  SVcv.ResizeTo(3,3);
  for(int nr =0; nr<SVcov.GetNrows(); nr++){
    for(int nc =0; nc<SVcov.GetNcols(); nc++){
      SVcv(nr,nc) = SVcov(nr,nc);
    }
  }

  TMatrixT<double> SVPVMatrix(3,1);
  for(int i=0; i<SVPVMatrix.GetNrows();i++){
    SVPVMatrix(i,0)=FD(i);
  }

  TMatrixT<double> SVPVMatrixT=SVPVMatrix;
  SVPVMatrixT.T();

  TMatrixT<double> lambda2 = SVPVMatrixT*(SVcv + PVcv)*SVPVMatrix;
  double sigmaabs = sqrt(lambda2(0,0))/SVPV.Mag();
  double sign = SVPV.Mag()/sigmaabs;

  return sign;
}

//// Generator Information
int Ntuple_Controller::matchTruth(TLorentzVector tvector){
	double testdr=0.3;
	int pdgid = 0;
	for(unsigned int i=0;i<NMCParticles();i++){
		if(MCParticle_p4(i).Pt()>0.){
			if(tvector.DeltaR(MCParticle_p4(i))<testdr){
				testdr = tvector.DeltaR(MCParticle_p4(i));
				pdgid = MCParticle_pdgid(i);
			}
		}
	}
	return pdgid;
}
bool Ntuple_Controller::matchTruth(TLorentzVector tvector, int pid, double dr){
	if (getMatchTruthIndex(tvector, pid, dr) >= 0) return true;
	return false;
}
int Ntuple_Controller::getMatchTruthIndex(TLorentzVector tvector, int pid, double dr){
	int index = -9;
	for(unsigned int i=0;i<NMCParticles();i++){
		if(MCParticle_p4(i).Pt()>0.){
			if(fabs(MCParticle_pdgid(i))==pid){
				if(tvector.DeltaR(MCParticle_p4(i))<dr) index = i;
			}
		}
	}
	return index;
}

int Ntuple_Controller::MCTau_true3prongAmbiguity(unsigned int i){
	int amb = -9;

	int j_a1 = MCTau_getDaughterOfType(i, PDGInfo::a_1_plus, true);
	if (j_a1 < 0){
		Logger(Logger::Warning) << "MCTau has no a1 decay product." << std::endl;
		return -9;
	}

	TLorentzVector genTauh	= MCTau_p4(i);
	TLorentzVector genA1	= MCTauandProd_p4(i, j_a1);

	TLorentzVector GenA1_boosted(BoostToRestFrame(genTauh, genA1));
	double dotProduct = GenA1_boosted.Vect().Dot( genTauh.Vect() );
	if(dotProduct < 0){
		amb = MultiProngTauSolver::plus;
	}
	else if(dotProduct > 0){
		amb = MultiProngTauSolver::minus;
	}
	else if(dotProduct == 0){
		amb = MultiProngTauSolver::zero;
	}

	return amb;
}
int Ntuple_Controller::MCTau_getDaughterOfType(unsigned int i_mcTau, int daughter_pdgid, bool ignoreCharge /*= true*/){
	int matchedIndex = -1;
	for(int i_dau=0; i_dau < NMCTauDecayProducts(i_mcTau); i_dau++){
		if( ignoreCharge ){
			if( abs(MCTauandProd_pdgid(i_mcTau, i_dau)) == abs(daughter_pdgid) )
				matchedIndex = i_dau;
		}
		else{
			if( MCTauandProd_pdgid(i_mcTau, i_dau) == daughter_pdgid )
				matchedIndex = i_dau;
		}
	}
	return matchedIndex;
}
int Ntuple_Controller::matchTauTruth(unsigned int i_hpsTau, bool onlyHadrDecays /*= false*/){
	int matchedIndex = -1;
	double minDr = 999;
	for(int i=0; i < NMCTaus(); i++){
		if( onlyHadrDecays && (MCTau_JAK(i) <= 2) ) continue; // exclude decays to electrons and muons
		double dr = MCTau_p4(i).DeltaR(PFTau_p4(i_hpsTau));
		if( dr < 0.5 && dr < minDr ){
			matchedIndex = i;
			minDr = dr;
		}
	}
	return matchedIndex;
}

//gets two TLorentzVectors (1. defines RF/Boost, 2. is boosted), creates a copy of the second and boosts it
TLorentzVector Ntuple_Controller::BoostToRestFrame(TLorentzVector TLV1, TLorentzVector TLV2){
	TVector3 boostvector = TLV1.BoostVector();
	TLorentzVector boosted_TLV2(TLV2);
	boosted_TLV2.Boost(-boostvector);
	return boosted_TLV2;
}

// get visible/invisible part of gen Tau 4-vector
TLorentzVector Ntuple_Controller::MCTau_invisiblePart(unsigned int i){
	TLorentzVector lv(0,0,0,0);
	for(int i_dau=1; i_dau < NMCTauDecayProducts(i); i_dau++){
		if( abs(MCTauandProd_pdgid(i, i_dau)) == PDGInfo::nu_e ||
			abs(MCTauandProd_pdgid(i, i_dau)) == PDGInfo::nu_mu ||
			abs(MCTauandProd_pdgid(i, i_dau)) == PDGInfo::nu_tau)
				lv += MCTauandProd_p4(i, i_dau);
	}
	if (lv.Pt() == 0.0)
		Logger(Logger::Warning) << "This tau decay has no neutrinos!?" << std::endl;

	return lv;
}
TLorentzVector Ntuple_Controller::MCTau_visiblePart(unsigned int i){
	TLorentzVector lv = MCTau_p4(i);
	lv -= MCTau_invisiblePart(i);
	return lv;
}

//// draw decay chain
void Ntuple_Controller::printMCDecayChainOfEvent(bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
	Logger(Logger::Info) << "=== Draw MC decay chain of event ===" << std::endl;
	for(unsigned int par = 0; par < NMCParticles(); par++){
		if ( !MCParticle_hasMother(par) )
			printMCDecayChainOfMother(par, printStatus, printPt, printEtaPhi, printQCD);
	}
}
void Ntuple_Controller::printMCDecayChainOfMother(unsigned int par, bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
	Logger(Logger::Info) << "Draw decay chain for mother particle at index " << par << " :" << std::endl;
	printMCDecayChain(par, 0, printStatus, printPt, printEtaPhi, printQCD);
}
void Ntuple_Controller::printMCDecayChain(unsigned int par, unsigned int level, bool printStatus, bool printPt, bool printEtaPhi, bool printQCD){
	std::ostringstream out;
	for(unsigned int l = 0; l < level; l++) out << "    ";
	out << MCParticleToString(par, printStatus, printPt, printEtaPhi);
	if (!printQCD && MCParticle_pdgid(par) == 92) {out << "    QCD stuff truncated ...\n"; std::cout << out.str(); return;}
	out << "\n";
	std::cout << out.str();
	for (unsigned int i_dau = 0; i_dau < MCParticle_childidx(par).size(); i_dau++){
		if (MCParticle_childidx(par).at(i_dau) != (int)par) // skip cases where particles are daughters of themselves
			printMCDecayChain(MCParticle_childidx(par).at(i_dau), level+1, printStatus, printPt, printEtaPhi, printQCD);
	}
}

std::string Ntuple_Controller::MCParticleToString(unsigned int par, bool printStatus, bool printPt, bool printEtaPhi){
	std::ostringstream out;
	out << "+-> ";
	out << PDGInfo::pdgIdToName( MCParticle_pdgid(par) );
	if (printStatus) out << " (status = " << MCParticle_status(par) << ", idx = " << par << ")";
	if (printPt || printEtaPhi) out <<  " [";
	if (printPt) out << "pT = " << MCParticle_p4(par).Pt() << "GeV";
	if (printEtaPhi) out << " eta = " << MCParticle_p4(par).Eta() << " phi = " << MCParticle_p4(par).Phi();
	if (printPt || printEtaPhi) out << "]";
	return out.str();
}


//// Trigger Information
bool Ntuple_Controller::TriggerAccept(TString n){
  unsigned int i=0;
  if(GetTriggerIndex(n,i))return TriggerAccept(i);
  return false;
}

unsigned int Ntuple_Controller::HLTPrescale(TString n){
  unsigned int i=0;
  if(GetTriggerIndex(n,i))return HLTPrescale(i);
  return 1;
}

unsigned int Ntuple_Controller::L1SEEDPrescale(TString n){
  unsigned int i=0;
  if(GetTriggerIndex(n,i))return L1SEEDPrescale(i);
  return 1;
}

bool Ntuple_Controller::GetTriggerIndex(TString n, unsigned int &i){
  for(i=0; i<Ntp->HTLTriggerName->size();i++){
    TString name=HTLTriggerName(i);
    if(name.Contains(n))return true;
  }
  return false;
}

double Ntuple_Controller::matchTrigger(TLorentzVector obj, std::vector<TString> trigger, std::string objectType){
  unsigned int id = 0;
  TLorentzVector triggerObj(0.,0.,0.,0.);
  if(objectType=="tau"){
    id = 84;
  }
  if(objectType=="muon"){
    id = 83;
  }
  if(objectType=="electron"){
    id = 82;
  }
  
  double minDR = 100.;
  for(unsigned int i_trig = 0; i_trig < trigger.size(); i_trig++){
    for(int i=0;i<NHLTTrigger_objs();i++){
      if(HLTTrigger_objs_trigger(i).find(trigger.at(i_trig)) != string::npos){
	for(int j=0;j<NHLTTrigger_objs(i);j++){
	  if(HLTTrigger_objs_Id(i,j)==(int)id){
	    triggerObj.SetPtEtaPhiE(HLTTrigger_objs_Pt(i,j),
				    HLTTrigger_objs_Eta(i,j),
				    HLTTrigger_objs_Phi(i,j),
				    HLTTrigger_objs_E(i,j));
	  }
	  if( triggerObj.Pt()>0. && obj.Pt()>0. ) {
	    double dr = obj.DeltaR(triggerObj);
	    if (dr < minDR) minDR = dr;
	  }
	}
      }
    }
  }
  return minDR;
}
bool Ntuple_Controller::matchTrigger(TLorentzVector obj, double dr_cut, std::vector<TString> trigger, std::string objectType){
	double dr = matchTrigger(obj, trigger, objectType);
	return dr < dr_cut;
}
bool Ntuple_Controller::matchTrigger(TLorentzVector obj, double dr_cut, TString trigger, std::string objectType){
	std::vector<TString> triggerVec;
	triggerVec.push_back(trigger);
	return matchTrigger(obj, dr_cut, triggerVec, objectType);
}


TMatrixTSym<double> Ntuple_Controller::PFTau_TIP_primaryVertex_cov(unsigned int i){
  TMatrixTSym<double> V_cov(LorentzVectorParticle::NVertex);
  int l=0;
  for(int j=0;j<LorentzVectorParticle::NVertex;j++){
    for(int k=j;k<LorentzVectorParticle::NVertex;k++){
      //if(j==k) V_cov(i,j)=pow(0.0001,2.0);
      V_cov(j,k)=Ntp->PFTau_TIP_primaryVertex_cov->at(i).at(l);
      V_cov(k,j)=Ntp->PFTau_TIP_primaryVertex_cov->at(i).at(l);
      l++;
    }
  }
  return  V_cov;
}

TMatrixTSym<double> Ntuple_Controller::PFTau_TIP_secondaryVertex_cov(unsigned int i){
  TMatrixTSym<double> V_cov(LorentzVectorParticle::NVertex);
  int l=0;
  for(int j=0;j<LorentzVectorParticle::NVertex;j++){
    for(int k=j;k<LorentzVectorParticle::NVertex;k++){
      V_cov(j,k)=Ntp->PFTau_TIP_secondaryVertex_cov->at(i).at(l);
      V_cov(k,j)=Ntp->PFTau_TIP_secondaryVertex_cov->at(i).at(l);
      l++;
    }
  }
  return  V_cov;
}

LorentzVectorParticle Ntuple_Controller::PFTau_a1_lvp(unsigned int i){
  TMatrixT<double>    a1_par(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixTSym<double> a1_cov(LorentzVectorParticle::NLorentzandVertexPar);
  int l=0;
  if(Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar){
    for(int k=0; k<LorentzVectorParticle::NLorentzandVertexPar; k++){
      a1_par(k,0)=Ntp->PFTau_a1_lvp->at(i).at(k);
      for(int j=k; j<LorentzVectorParticle::NLorentzandVertexPar; j++){
	a1_cov(k,j)=Ntp->PFTau_a1_cov->at(i).at(l);
	a1_cov(j,k)=Ntp->PFTau_a1_cov->at(i).at(l);
	l++;
      } 
    }
  }
  return LorentzVectorParticle(a1_par,a1_cov,Ntp->PFTau_a1_pdgid->at(i).at(0),Ntp->PFTau_a1_charge->at(i).at(0),Ntp->PFTau_a1_B->at(i).at(0));
}

std::vector<TrackParticle> Ntuple_Controller::PFTau_daughterTracks(unsigned int i){
  std::vector<TrackParticle> daughter;
  //std::cout << "selTau no. " << i << std::endl;
  //std::cout << "Ntp->PFTau_daughterTracks_poca->size() " << Ntp->PFTau_daughterTracks_poca->size() << std::endl;
  //std::cout << "Ntp->PFTau_daughterTracks->size() " << Ntp->PFTau_daughterTracks->size() << std::endl;
  for(unsigned int d=0;d<Ntp->PFTau_daughterTracks->at(i).size();d++){
    //std::cout << "Pion daughter no. " << d << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_poca->at(i).size() " << Ntp->PFTau_daughterTracks_poca->at(i).size() << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks->at(i).size() " << Ntp->PFTau_daughterTracks->at(i).size() << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_M->at(i).size() " << Ntp->PFTau_daughterTracks_M->at(i).size() << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_charge->at(i).size() " << Ntp->PFTau_daughterTracks_charge->at(i).size() << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_pdgid->at(i).size() " << Ntp->PFTau_daughterTracks_pdgid->at(i).size() << std::endl;
    TMatrixT<double>    a1_par(TrackParticle::NHelixPar,1);
    TMatrixTSym<double> a1_cov(TrackParticle::NHelixPar);
    int l=0;
    //std::cout << "Ntp->PFTau_daughterTracks_M->at(i).at(d)  Pion Mass: " << Ntp->PFTau_daughterTracks_M->at(i).at(d) << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_poca->at(i).at(d) " << Ntp->PFTau_daughterTracks_poca->at(i).at(d).at(0) << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_poca->at(i).at(d) " << Ntp->PFTau_daughterTracks_poca->at(i).at(d).at(1) << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_poca->at(i).at(d) " << Ntp->PFTau_daughterTracks_poca->at(i).at(d).at(2) << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_charge->at(i).at(d) " << Ntp->PFTau_daughterTracks_charge->at(i).at(d) << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_pdgid->at(i).at(d) " << Ntp->PFTau_daughterTracks_pdgid->at(i).at(d) << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_B->at(i).at(d) " << Ntp->PFTau_daughterTracks_B->at(i).at(d) << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks_M->at(i).at(d) " << Ntp->PFTau_daughterTracks_M->at(i).at(d) << std::endl;
    //std::cout << "Ntp->PFTau_daughterTracks->at(i).at(d).size() " << Ntp->PFTau_daughterTracks->at(i).at(d).size() << std::endl;
    if(Ntp->PFTau_daughterTracks->at(i).at(d).size() == TrackParticle::NHelixPar){
    for(int k=0; k<TrackParticle::NHelixPar; k++){
      a1_par(k,0)=Ntp->PFTau_daughterTracks->at(i).at(d).at(k);
      for(int j=k; j<TrackParticle::NHelixPar; j++){
	a1_cov(k,j)=Ntp->PFTau_daughterTracks_cov->at(i).at(d).at(l);
	a1_cov(j,k)=Ntp->PFTau_daughterTracks_cov->at(i).at(d).at(l);
	l++;
      }
    }
    }
    daughter.push_back(TrackParticle(a1_par,a1_cov,Ntp->PFTau_daughterTracks_pdgid->at(i).at(d),Ntp->PFTau_daughterTracks_M->at(i).at(d),Ntp->PFTau_daughterTracks_charge->at(i).at(d),Ntp->PFTau_daughterTracks_B->at(i).at(d)));
  }
  return daughter;
}

std::vector<TVector3> Ntuple_Controller::PFTau_daughterTracks_poca(unsigned int i){
  std::vector<TVector3> poca;
  for(unsigned int k=0;k<Ntp->PFTau_daughterTracks_poca->at(i).size();k++){
    poca.push_back(TVector3(Ntp->PFTau_daughterTracks_poca->at(i).at(k).at(0),Ntp->PFTau_daughterTracks_poca->at(i).at(k).at(1),Ntp->PFTau_daughterTracks_poca->at(i).at(k).at(2)));
  }
  return poca;
}

//All ReFitTracks are evaluated at the fitted secondary vertex TODO: save TrackParticle also for full helix parametrization
std::vector<TLorentzVector> Ntuple_Controller::PFTau_daughterReFitTracks_p4(unsigned int i){
  std::vector<TLorentzVector> refittracks_p4;
  for(unsigned int d=0; d<Ntp->PFTau_PionsP4->at(i).size(); d++){
    TLorentzVector TLV(Ntp->PFTau_PionsP4->at(i).at(d).at(1), Ntp->PFTau_PionsP4->at(i).at(d).at(2), Ntp->PFTau_PionsP4->at(i).at(d).at(3), Ntp->PFTau_PionsP4->at(i).at(d).at(0));
    refittracks_p4.push_back(TLV);
  }
  return refittracks_p4;
}

bool Ntuple_Controller::AmbiguitySolver(std::vector<bool> A1Fit, std::vector<bool> EventFit, std::vector<double> Probs, int &IndexToReturn, bool &AmbiguityPoint){

  if(EventFit.at(0) == true && EventFit.at(1) == false && EventFit.at(2) == false){IndexToReturn =0; AmbiguityPoint = true; return true;}
  if(EventFit.at(1) == true && EventFit.at(2) == false){ IndexToReturn = 1;AmbiguityPoint = false;return true;}
  if(EventFit.at(1) == false && EventFit.at(2) == true){ IndexToReturn = 2;AmbiguityPoint = false;return true;}

  if((A1Fit.at(1) == true && A1Fit.at(2) == true) && (EventFit.at(1) == true && EventFit.at(2) == true)){
    if(Probs.at(1)  >Probs.at(2) ){ IndexToReturn  =1;AmbiguityPoint = false;return true;}
    if(Probs.at(1)  <Probs.at(2) ){ IndexToReturn  =2;AmbiguityPoint = false;return true;}
    //  std::cout<< " probs 1,2  "<< Probs.at(1) <<"  "<< Probs.at(2)<<std::endl;
  }
  return false;
}
bool Ntuple_Controller::AmbiguitySolverByChi2(std::vector<bool> A1Fit, std::vector<bool> EventFit, std::vector<double> Chi2s, int &IndexToReturn, bool &AmbiguityPoint){

  if(EventFit.at(0) == true && EventFit.at(1) == false && EventFit.at(2) == false){IndexToReturn =0; AmbiguityPoint = true; return true;}
  if(EventFit.at(1) == true && EventFit.at(2) == false){ IndexToReturn = 1; AmbiguityPoint = false;return true;}
  if(EventFit.at(1) == false && EventFit.at(2) == true){ IndexToReturn = 2; AmbiguityPoint = false;return true;}

  if((A1Fit.at(1) == true && A1Fit.at(2) == true) && (EventFit.at(1) == true && EventFit.at(2) == true)){
    if(Chi2s.at(1) < Chi2s.at(2) && Chi2s.at(1) > 0){ IndexToReturn =1; AmbiguityPoint = false;return true;}
    if(Chi2s.at(1) > Chi2s.at(2) && Chi2s.at(2) > 0){ IndexToReturn =2; AmbiguityPoint = false;return true;}
    //  std::cout<< " probs 1,2  "<< Probs.at(1) <<"  "<< Probs.at(2)<<std::endl;
  }
  return false;
}

TMatrixTSym<double> Ntuple_Controller::PF_Tau_FlightLegth3d_TauFrame_cov(unsigned int i){
  TVector3 f=PFTau_FlightLength3d(i);
  TMatrixT<double> Res(5,1);
  Res(0,0)=f.X();
  Res(1,0)=f.Y();
  Res(2,0)=f.Z();
  Res(3,0)=f.Phi();
  Res(4,0)=f.Theta();
  TMatrixTSym<double> ResCov(5);
  TMatrixTSym<double> cov=PFTau_FlightLength3d_cov(i);
  for(int s=0;s<LorentzVectorParticle::NVertex;s++){
    for(int t=0;t<LorentzVectorParticle::NVertex;t++){
      ResCov(s,t)=cov(s,t);
    }
  }
  TMatrixT<double> Resp=MultiProngTauSolver::RotateToTauFrame(Res);
  TMatrixTSym<double> RespCov=ErrorMatrixPropagator::PropagateError(&MultiProngTauSolver::RotateToTauFrame,Res,ResCov);
  for(int s=0;s<LorentzVectorParticle::NVertex;s++){
    for(int t=0;t<LorentzVectorParticle::NVertex;t++){
      cov(s,t)=RespCov(s,t);
    }
  }
  return cov;
}

TVector3 Ntuple_Controller::PF_Tau_FlightLegth3d_TauFrame(unsigned int i){
  TVector3 f=PFTau_FlightLength3d(i);
  TMatrixT<double> Res(5,1);
  Res(0,0)=f.X();
  Res(1,0)=f.Y();
  Res(2,0)=f.Z();
  Res(3,0)=f.Phi();
  Res(4,0)=f.Theta();
  TMatrixT<double> Resp=MultiProngTauSolver::RotateToTauFrame(Res);
  return TVector3(Resp(0,0),Resp(1,0),Resp(2,0));
}

double Ntuple_Controller::PFTau_Mass(unsigned int i){
  return Ntp->PFTau_Mass(i);
}

double Ntuple_Controller::dxySigned(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return (-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt();
}
double Ntuple_Controller::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(dxySigned(fourvector, poca, vtx));
}


double Ntuple_Controller::dzSigned(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2);
}
double Ntuple_Controller::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(dzSigned(fourvector, poca, vtx));
}

double Ntuple_Controller::vertexSignificance(TVector3 vec, unsigned int vertex){
	if(vertex>=0 && vertex<NVtx()){
		const double elm[3] = {(vec.X()-Vtx(vertex).X()),(vec.Y()-Vtx(vertex).Y()),(vec.Z()-Vtx(vertex).Z())};
		TVectorD diff(3,elm);
		TMatrixF M(Vtx_Cov(vertex));
		if(M.IsValid()){
			double mag = diff.Norm2Sqr();
			double sim = M.Similarity(diff);
			return mag/sqrt(sim);
		}
	}
	return 999;
}

// check if given lepton was used for MVA-MET calculation
bool Ntuple_Controller::findCorrMVASrcMuon(unsigned int muon_idx, int &mvaSrcMuon_idx, float &dR ){
	float minDr = 1000;
	float dr = 1001;
	for(unsigned int i_mvaLep = 0; i_mvaLep < NMET_CorrMVA_srcMuons(); i_mvaLep++){
		 dr = Tools::dr(Muon_p4(muon_idx),MET_CorrMVA_srcMuon_p4(i_mvaLep));
		 if ((dr < 0.05) && (dr < minDr)) {
			 minDr = dr;
			 dR = dr;
			 mvaSrcMuon_idx = i_mvaLep;
		 }
	}
	if (minDr < 0.05) return true;
	else return false;
}
bool Ntuple_Controller::findCorrMVASrcElectron(unsigned int elec_idx, int &mvaSrcElectron_idx, float &dR ){
	float minDr = 1000;
	float dr = 1001;
	for (unsigned int i_mvaLep = 0; i_mvaLep < NMET_CorrMVA_srcElectrons(); i_mvaLep++){
		 dr = Tools::dr(Electron_p4(elec_idx),MET_CorrMVA_srcElectron_p4(i_mvaLep));
		 if ((dr < 0.05) && (dr < minDr)) {
			 minDr = dr;
			 dR = dr;
			 mvaSrcElectron_idx = i_mvaLep;
		 }
	}
	if (minDr < 0.05) return true;
	else return false;
}
bool Ntuple_Controller::findCorrMVASrcTau(unsigned int tau_idx, int &mvaSrcTau_idx, float &dR ){
	float minDr = 1000;
	float dr = 1001;
	for(unsigned int i_mvaLep = 0; i_mvaLep < NMET_CorrMVA_srcTaus(); i_mvaLep++){
		 dr = Tools::dr(PFTau_p4(tau_idx),MET_CorrMVA_srcTau_p4(i_mvaLep));
		 if ((dr < 0.05) && (dr < minDr)) {
			 minDr = dr;
			 dR = dr;
			 mvaSrcTau_idx = i_mvaLep;
		 }
	}
	if (minDr < 0.05) return true;
	else return false;
}
bool Ntuple_Controller::findCorrMVAMuTauSrcMuon(unsigned int muon_idx, int &mvaMuTauSrcMuon_idx, float &dR ){
	float minDr = 1000;
	float dr = 1001;
	for (unsigned int i_mvaLep = 0; i_mvaLep < NMET_CorrMVAMuTau_srcMuons(); i_mvaLep++){
		 dr = Tools::dr(Muon_p4(muon_idx),MET_CorrMVAMuTau_srcMuon_p4(i_mvaLep));
		 if ((dr < 0.05) && (dr < minDr)) {
			 minDr = dr;
			 dR = dr;
			 mvaMuTauSrcMuon_idx = i_mvaLep;
		 }
	}
	if (minDr < 0.05) return true;
	else return false;
}
bool Ntuple_Controller::findCorrMVAMuTauSrcTau(unsigned int tau_idx, int &mvaMuTauSrcTau_idx, float &dR ){
	float minDr = 1000;
	float dr = 1001;
	for (unsigned int i_mvaLep = 0; i_mvaLep < NMET_CorrMVAMuTau_srcTaus(); i_mvaLep++){
		 dr = Tools::dr(PFTau_p4(tau_idx),MET_CorrMVAMuTau_srcTau_p4(i_mvaLep));
		 if ((dr < 0.05) && (dr < minDr)) {
			 minDr = dr;
			 dR = dr;
			 mvaMuTauSrcTau_idx = i_mvaLep;
		 }
	}
	if (minDr < 0.05) return true;
	else return false;
}

// function to sort any objects by any value in descending order
std::vector<int> Ntuple_Controller::sortObjects(std::vector<int> indices, std::vector<double> values){
	if (indices.size() != values.size()){
		Logger(Logger::Warning) << "Please make sure indices and values have same size for sorting. Abort." << std::endl;
		return std::vector<int>();
	}
	// create vector of pairs to allow for sorting by value
	std::vector< std::pair<int, double> > pairs;
	for(unsigned int i = 0; i<values.size(); i++ ){
		pairs.push_back( std::make_pair(indices.at(i),values.at(i)) );
	}
	// sort vector of pairs
	std::sort(pairs.begin(), pairs.end(), sortIdxByValue());
	// create vector of indices in correct order
	std::vector<int> sortedIndices;
	for(unsigned int i = 0; i<pairs.size(); i++){
		sortedIndices.push_back(pairs.at(i).first);
	}
	return sortedIndices;
}

std::vector<int> Ntuple_Controller::sortDefaultObjectsByPt(TString objectType){
	std::vector<int> indices;
	std::vector<double> values;
	if (objectType == "Jets" || objectType == "PFJets"){
	  for(unsigned int i = 0; i<NPFJets(); i++ ){
	    indices.push_back(i);
	    values.push_back(PFJet_p4(i).Pt());
	  }
	}
	else if(objectType == "Taus" || objectType == "PFTaus"){
	  for(unsigned int i = 0; i<NPFTaus(); i++ ){
	    indices.push_back(i);
	    values.push_back(PFTau_p4(i).Pt());
	  }
	}
	else if(objectType == "Muons"){
	  for (unsigned int i = 0; i<NMuons(); i++ ){
	    indices.push_back(i);
	    values.push_back(Muon_p4(i).Pt());
	  }
	}
	else if(objectType == "Electrons"){
	  for (unsigned int i = 0; i<NElectrons(); i++ ){
	    indices.push_back(i);
	    values.push_back(Electron_p4(i).Pt());
	  }
	}
	else{
	  Logger(Logger::Warning) << "sortDefaultObjectsByPt is only implemented for Jets, Taus, Muons and Electrons. Abort." << std::endl;
	  return std::vector<int>();
	}
	return sortObjects(indices, values);
}



// obtain, or create and store, SVFit results from/on dCache
#ifdef USE_SVfit
SVFitObject* Ntuple_Controller::getSVFitResult_MuTauh(SVFitStorage& svFitStor, TString metType, unsigned muIdx, unsigned tauIdx, unsigned rerunEvery /* = 5000 */, TString suffix /* ="" */, double scaleMu /* =1 */, double scaleTau /* =1 */) {
	 // configure svfitstorage on first call
	if ( !svFitStor.isConfigured() ) svFitStor.Configure(GetInputDatasetName(), suffix);
	// get SVFit result from cache
	SVFitObject* svfObj = svFitStor.GetEvent(RunNumber(), LuminosityBlock(), EventNumber());
	// if obtained object is not valid, create and store it
	if (!svfObj->isValid()) {
		runAndSaveSVFit_MuTauh(svfObj, svFitStor, metType, muIdx, tauIdx, scaleMu, scaleTau);
	}
	else{
		// calculate every N'th event and compare with what is stored
		if( (EventNumber() % rerunEvery) == 123){
			SVFitObject* newSvfObj = new SVFitObject();
			runAndSaveSVFit_MuTauh(newSvfObj, svFitStor, metType, muIdx, tauIdx, scaleMu, scaleTau, false); // will not be saved in output files

			if (*svfObj == *newSvfObj){
				Logger(Logger::Info) << "Recalculation of SVFit object gave same result." << std::endl;
			}
			else {
				Logger(Logger::Warning) << "Recalculation of SVFit object gave DIFFERENT result!!" <<
				"\n\told: mass = " << svfObj->get_mass() << " +/- " << svfObj->get_massUncert() << ", pt = " << svfObj->get_pt() << " +/- " << svfObj->get_ptUncert() <<
				"\n\tnew: mass = " << newSvfObj->get_mass() << " +/- " << newSvfObj->get_massUncert() << ", pt = " << newSvfObj->get_pt() << " +/- " << newSvfObj->get_ptUncert() <<
				"\n\tSmall discrepancies could be caused by fitting details. It's up to you whether to ignore them." << std::endl;
			}
			delete newSvfObj;
		}
	}
	return svfObj;
}
SVFitObject* Ntuple_Controller::getSVFitResult_MuTau3p(SVFitStorage& svFitStor, TString metType, unsigned muIdx, TLorentzVector tauLV, LorentzVectorParticle neutrino, TString suffix /* ="" */, double scaleMu /* =1 */, double scaleTau /* =1 */) {
	 // configure svfitstorage on first call
	if ( !svFitStor.isConfigured() ) svFitStor.Configure(GetInputDatasetName(), suffix);
	// get SVFit result from cache
	SVFitObject* svfObj = svFitStor.GetEvent(RunNumber(), LuminosityBlock(), EventNumber());
	// if obtained object is not valid, create and store it
	if (!svfObj->isValid()) {
		runAndSaveSVFit_MuTau3p(svfObj, svFitStor, metType, muIdx, tauLV, neutrino, scaleMu, scaleTau);
	}
	return svfObj;
}

// create SVFitObject from standard muon and standard tau_h
void Ntuple_Controller::runAndSaveSVFit_MuTauh(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned muIdx, unsigned tauIdx, double scaleMu, double scaleTau, bool save /*= true*/) {
	objects::MET met(this, metType);
	SVfitProvider svfProv(this, met, "Mu", muIdx, "Tau", tauIdx, 1, scaleMu, scaleTau);
	*svfObj = svfProv.runAndMakeObject();
	if (svfObj->isValid()) {
		// store only if object is valid
		if (save) svFitStor.SaveEvent(RunNumber(), LuminosityBlock(), EventNumber(), svfObj);
	} else {
		Logger(Logger::Error) << "Unable to create a valid SVFit object." << std::endl;
	}
}
// create SVFitObject from standard muon and fully reconstructed 3prong tau
void Ntuple_Controller::runAndSaveSVFit_MuTau3p(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned muIdx, TLorentzVector tauLV, LorentzVectorParticle neutrino, double scaleMu, double scaleTau, bool save /*= true*/){
	objects::MET met(this, metType);
	met.subtractNeutrino(neutrino);
	SVfitProvider svfProv(this, met, "Mu", muIdx, tauLV, 1, scaleMu, scaleTau);
	*svfObj = svfProv.runAndMakeObject();
	if (svfObj->isValid()) {
		// store only if object is valid
		if (save) svFitStor.SaveEvent(RunNumber(), LuminosityBlock(), EventNumber(), svfObj);
	} else {
		Logger(Logger::Error) << "Unable to create a valid SVFit object." << std::endl;
	}
}

#endif // USE_SVfit
