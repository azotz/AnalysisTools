#include "Selection_Factory.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

#include "Example.h"
#include "TauSpinExample.h"
#ifdef USE_cherepanov
#include "cherepanov/Validation.h"
#include "cherepanov/Ztotautau_hadmu_ControlSample.h"
#include "cherepanov/Tau_momentum_calculation.h"
#endif
#ifdef USE_inugent
#include "inugent/Ztotautau_ControlSample.h"
#include "inugent/ChargedHiggs_dilepontic.h"
#include "inugent/ChargedHiggs_tauplusjet.h"
#include "inugent/ZDouble3prong.h"
#include "inugent/Ztomumu_ControlSample.h"
#include "inugent/TriggerStudy.h"
#include "inugent/TriggerStudyMC.h"
#include "inugent/TauLifeTime.h"
#endif
#ifdef USE_nehrkorn
#include "nehrkorn/ZtoEMu.h"
#include "nehrkorn/Tvariable_EE.h"
#include "nehrkorn/Tvariable_MuMu.h"
#include "nehrkorn/Tvariable_EMu.h"
#endif
#ifdef USE_kargoll
#include "kargoll/HToTaumuTauh.h"
//#include "kargoll/MuTauSync.h"
#include "kargoll/OneJetBoost.h"
#include "kargoll/OneJetHigh.h"
#include "kargoll/OneJetLow.h"
#include "kargoll/VBFTight.h"
#include "kargoll/VBFLoose.h"
#include "kargoll/ZeroJetHigh.h"
#include "kargoll/ZeroJetLow.h"
#include "kargoll/Inclusive.h"
//#include "kargoll/MCDecayChain.h"
#include "kargoll/ZeroJet3Prong.h"
#endif
#ifdef USE_pistone

#endif
#ifdef USE_zotz
#include "zotz/ZToTaumuTauh.h"
#endif
#ifdef USE_stahl

#endif

Selection_Factory::Selection_Factory(){
}

Selection_Factory::~Selection_Factory(){
}

Selection_Base* Selection_Factory::Factory(TString Analysis, TString UncertType,int mode, int runtype, double lumi){
  Selection_Base* s;
  Analysis.ToLower();

  // ensuring code will compile independently of user code
  // WARNING: be aware of the consequences of "Contains". Make sure that Class "foo" is put after "foobar".
  if(Analysis.Contains("example"))s=new Example(Analysis,UncertType);
  else if(Analysis.Contains("tauspin"))s=new TauSpinExample(Analysis,UncertType);
#ifdef USE_cherepanov
  else if(Analysis.Contains("validation"))s=new Validation(Analysis,UncertType);
  else if(Analysis.Contains("ztotautau_hadmu_controlsample"))s=new Ztotautau_hadmu_ControlSample(Analysis,UncertType);
  else if(Analysis.Contains("tau_momentum_calculation"))s=new Tau_momentum_calculation(Analysis,UncertType);
#endif
#ifdef USE_inugent
  else if(Analysis.Contains("ztotautau_controlsample"))s=new Ztotautau_ControlSample(Analysis,UncertType);
  else if(Analysis.Contains("chargedhiggs_dilepton"))s=new ChargedHiggs_dilepontic(Analysis,UncertType);
  else if(Analysis.Contains("chargedhiggs_tauplusjet"))s=new ChargedHiggs_tauplusjet(Analysis,UncertType);
  else if(Analysis.Contains("zdouble3prong"))s=new ZDouble3prong(Analysis,UncertType);
  else if(Analysis.Contains("ztomumu_controlsample"))s=new Ztomumu_ControlSample(Analysis,UncertType);
  else if(Analysis.Contains("triggerstudymc"))s=new TriggerStudyMC(Analysis,UncertType);
  else if(Analysis.Contains("triggerstudy"))s=new TriggerStudy(Analysis,UncertType);
  else if(Analysis.Contains("taulifetime"))s=new TauLifeTime(Analysis,UncertType);
#endif
#ifdef USE_nehrkorn
  else if(Analysis.Contains("ztoemu"))s=new ZtoEMu(Analysis,UncertType);
  else if(Analysis.Contains("tvariable_ee"))s=new Tvariable_EE(Analysis,UncertType);
  else if(Analysis.Contains("tvariable_mumu"))s=new Tvariable_MuMu(Analysis,UncertType);
  else if(Analysis.Contains("tvariable_emu"))s=new Tvariable_EMu(Analysis,UncertType);
#endif
#ifdef USE_kargoll
  else if(Analysis.Contains("zerojet3prong")) s=new ZeroJet3Prong(Analysis,UncertType);
  else if(Analysis.Contains("htotaumutauh")) s=new HToTaumuTauh(Analysis,UncertType);
  else if(Analysis.Contains("onejetboost")) s=new OneJetBoost(Analysis,UncertType);
  else if(Analysis.Contains("onejethigh")) s=new OneJetHigh(Analysis,UncertType);
  else if(Analysis.Contains("onejetlow")) s=new OneJetLow(Analysis,UncertType);
  else if(Analysis.Contains("vbftight")) s=new VBFTight(Analysis,UncertType);
  else if(Analysis.Contains("vbfloose")) s=new VBFLoose(Analysis,UncertType);
  else if(Analysis.Contains("zerojethigh")) s=new ZeroJetHigh(Analysis,UncertType);
  else if(Analysis.Contains("zerojetlow")) s=new ZeroJetLow(Analysis,UncertType);
  else if(Analysis.Contains("inclusive")) s=new Inclusive(Analysis,UncertType);
  //else if(Analysis.Contains("mutausync")) s=new MuTauSync(Analysis,UncertType);
  //else if(Analysis.Contains("mcdecaychain")) s=new MCDecayChain(Analysis,UncertType);
#endif
#ifdef USE_pistone

#endif
#ifdef USE_zotz
  else if(Analysis.Contains("ztotaumutauh")) s=new ZToTaumuTauh(Analysis,UncertType);
#endif
#ifdef USE_stahl

#endif
  else{
	Logger(Logger::Error)<< "Invalid Analysis type \"" << Analysis << "\". Using default <Example.h> " << std::endl;
    s=new Example(Analysis,UncertType);
  }
  s->SetMode(mode);
  s->SetRunType(runtype);
  s->SetLumi(lumi);
  s->Configure();
  return s;
}
