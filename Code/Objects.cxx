/*
 * Objects.cxx
 *
 *  Created on: Jan 19, 2015
 *      Author: kargoll
 */

#include "Objects.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "Ntuple_Controller.h"
#include "TVector2.h"

namespace objects {

// default constructor
MET::MET() {}

// constructor of MET object from Ntuple_Controller and type of MET
MET::MET(Ntuple_Controller* const Ntp, TString met_type) {
	metType_ = met_type;

	Init(Ntp);

}

void MET::Init(Ntuple_Controller* const Ntp){
	if(metType_ == "Uncorr"){
	   et_             = Ntp->MET_Uncorr_et();
	   phi_            = Ntp->MET_Uncorr_phi();
	   ex_             = Ntp->MET_Uncorr_ex();
	   ey_             = Ntp->MET_Uncorr_ey();
	   significance_   = Ntp->MET_Uncorr_significance();
	   significanceXX_ = Ntp->MET_Uncorr_significance_xx();
	   significanceXY_ = Ntp->MET_Uncorr_significance_xy();
	   significanceYY_ = Ntp->MET_Uncorr_significance_yy();
	   hasSignificance_= true;
	}
	if(metType_ == "CorrT0rt"){
	   et_             = Ntp->MET_CorrT0rt_et();
	   phi_            = Ntp->MET_CorrT0rt_phi();
	   ex_             = Ntp->MET_CorrT0rt_ex();
	   ey_             = Ntp->MET_CorrT0rt_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrT0rtT1"){
	   et_             = Ntp->MET_CorrT0rtT1_et();
	   phi_            = Ntp->MET_CorrT0rtT1_phi();
	   ex_             = Ntp->MET_CorrT0rtT1_ex();
	   ey_             = Ntp->MET_CorrT0rtT1_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrT0pc"){
	   et_             = Ntp->MET_CorrT0pc_et();
	   phi_            = Ntp->MET_CorrT0pc_phi();
	   ex_             = Ntp->MET_CorrT0pc_ex();
	   ey_             = Ntp->MET_CorrT0pc_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrT0pcT1"){
	   et_             = Ntp->MET_CorrT0pcT1_et();
	   phi_            = Ntp->MET_CorrT0pcT1_phi();
	   ex_             = Ntp->MET_CorrT0pcT1_ex();
	   ey_             = Ntp->MET_CorrT0pcT1_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrT0rtTxy"){
	   et_             = Ntp->MET_CorrT0rtTxy_et();
	   phi_            = Ntp->MET_CorrT0rtTxy_phi();
	   ex_             = Ntp->MET_CorrT0rtTxy_ex();
	   ey_             = Ntp->MET_CorrT0rtTxy_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrT0rtT1Txy"){
	   et_             = Ntp->MET_CorrT0rtT1Txy_et();
	   phi_            = Ntp->MET_CorrT0rtT1Txy_phi();
	   ex_             = Ntp->MET_CorrT0rtT1Txy_ex();
	   ey_             = Ntp->MET_CorrT0rtT1Txy_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrT0pcTxy"){
	   et_             = Ntp->MET_CorrT0pcTxy_et();
	   phi_            = Ntp->MET_CorrT0pcTxy_phi();
	   ex_             = Ntp->MET_CorrT0pcTxy_ex();
	   ey_             = Ntp->MET_CorrT0pcTxy_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrT0pcT1Txy"){
	   et_             = Ntp->MET_CorrT0pcT1Txy_et();
	   phi_            = Ntp->MET_CorrT0pcT1Txy_phi();
	   ex_             = Ntp->MET_CorrT0pcT1Txy_ex();
	   ey_             = Ntp->MET_CorrT0pcT1Txy_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrT1"){
	   et_             = Ntp->MET_CorrT1_et();
	   phi_            = Ntp->MET_CorrT1_phi();
	   ex_             = Ntp->MET_CorrT1_ex();
	   ey_             = Ntp->MET_CorrT1_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrT1Txy"){
	   et_             = Ntp->MET_CorrT1Txy_et();
	   phi_            = Ntp->MET_CorrT1Txy_phi();
	   ex_             = Ntp->MET_CorrT1Txy_ex();
	   ey_             = Ntp->MET_CorrT1Txy_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrCaloT1"){
	   et_             = Ntp->MET_CorrCaloT1_et();
	   phi_            = Ntp->MET_CorrCaloT1_phi();
	   ex_             = Ntp->MET_CorrCaloT1_ex();
	   ey_             = Ntp->MET_CorrCaloT1_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrCaloT1T2"){
	   et_             = Ntp->MET_CorrCaloT1T2_et();
	   phi_            = Ntp->MET_CorrCaloT1T2_phi();
	   ex_             = Ntp->MET_CorrCaloT1T2_ex();
	   ey_             = Ntp->MET_CorrCaloT1T2_ey();
	   hasSignificance_= false;
	}
	if(metType_ == "CorrMVA"){
	   et_             = Ntp->MET_CorrMVA_et();
	   phi_            = Ntp->MET_CorrMVA_phi();
	   ex_             = Ntp->MET_CorrMVA_ex();
	   ey_             = Ntp->MET_CorrMVA_ey();
	   significance_   = Ntp->MET_CorrMVA_significance();
	   significanceXX_ = Ntp->MET_CorrMVA_significance_xx();
	   significanceXY_ = Ntp->MET_CorrMVA_significance_xy();
	   significanceYY_ = Ntp->MET_CorrMVA_significance_yy();
	   hasSignificance_= true;
	}
	if(metType_ == "CorrMVAMuTau"){
	   et_             = Ntp->MET_CorrMVAMuTau_et();
	   phi_            = Ntp->MET_CorrMVAMuTau_phi();
	   ex_             = Ntp->MET_CorrMVAMuTau_ex();
	   ey_             = Ntp->MET_CorrMVAMuTau_ey();
	   significance_   = Ntp->MET_CorrMVAMuTau_significance();
	   significanceXX_ = Ntp->MET_CorrMVAMuTau_significance_xx();
	   significanceXY_ = Ntp->MET_CorrMVAMuTau_significance_xy();
	   significanceYY_ = Ntp->MET_CorrMVAMuTau_significance_yy();
	   hasSignificance_= true;
	}
	else
		Logger(Logger::Error) << "MET type " << metType_ << " unknown!" << std::endl;
}

MET::~MET() {
}

void MET::subtractNeutrino(LorentzVectorParticle neutrino){
	metType_ += "NeutrinoSubtracted";

	TVectorD metvec(2), nuvec(2);
	metvec[0] = ex_;
	metvec[1] = ey_;
	nuvec[0] = neutrino.LV().Px();
	nuvec[1] = neutrino.LV().Py();
	TVectorD diff = metvec - nuvec;
	ex_ = diff[0];
	ey_ = diff[1];
	et_ = sqrt(ex_*ex_ + ey_*ey_);
	TVector2 temp(ex_, ey_);
	phi_ = temp.Phi();

	significanceXX_ = significanceXX_ + neutrino.Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px);
	significanceXY_ = significanceXY_ + neutrino.Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py);
	significanceYY_ = significanceYY_ + neutrino.Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py);
	TMatrixD metmat = significanceMatrix<TMatrixD>();
	if (fabs(metmat.Determinant())>0.00001){
		metmat.Invert();
		significance_ = diff * (metmat * diff);
	}
	else{
		Logger(Logger::Warning) << "Cannot update significance of MET object." << std::endl;
		significance_ = -1;
		hasSignificance_ = false;
	}
}

Vector3D MET::met3D() const {
	return Vector3D(ex_, ey_, 0.);
}

} /* namespace objects */
