#include "muon_scare.hpp"

#include <boost/math/special_functions/erf.hpp>
#include <memory>
#include "TRandom3.h"


using namespace std;

namespace scarekit {

CrystalBall::CrystalBall():m(0),s(1),a(10),n(10){
    init();
}
CrystalBall::CrystalBall(double mean, double sigma, double alpha, double n_)
    :m(mean),s(sigma),a(alpha),n(n_){
    init();
}
void CrystalBall::init(){
    double fa = fabs(a);
    double ex = exp(-fa*fa/2);
    double A  = pow(n/fa, n) * ex;
    double C1 = n/fa/(n-1) * ex; 
    double D1 = 2 * sqrtPiOver2 * erf(fa/sqrt2);
    B = n/fa-fa;
    C = (D1+2*C1)/C1;   
    D = (D1+2*C1)/2;   
    N = 1.0/s/(D1+2*C1); 
    k = 1.0/(n-1);  
    NA = N*A;       
    Ns = N*s;       
    NC = Ns*C1;     
    F = 1-fa*fa/n; 
    G = s*n/fa;    
    cdfMa = cdf(m-a*s);
    cdfPa = cdf(m+a*s);
}
double CrystalBall::pdf(double x) const{ 
    double d=(x-m)/s;
    if(d<-a) return NA*pow(B-d, -n);
    if(d>a) return NA*pow(B+d, -n);
    return N*exp(-d*d/2);
}
double CrystalBall::pdf(double x, double ks, double dm) const{ 
    double d=(x-m-dm)/(s*ks);
    if(d<-a) return NA/ks*pow(B-d, -n);
    if(d>a) return NA/ks*pow(B+d, -n);
    return N/ks*exp(-d*d/2);

}
double CrystalBall::cdf(double x) const{
    double d = (x-m)/s;
    if(d<-a) return NC / pow(F-s*d/G, n-1);
    if(d>a) return NC * (C - pow(F+s*d/G, 1-n) );
    return Ns * (D - sqrtPiOver2 * erf(-d/sqrt2));
}
double CrystalBall::invcdf(double u) const{
    if(u<cdfMa) return m + G*(F - pow(NC/u, k));
    if(u>cdfPa) return m - G*(F - pow(C-u/NC, -k) );
    return m - sqrt2 * s * boost::math::erf_inv((D - u/Ns )/sqrtPiOver2);
}

double get_rndm(double eta, float nL, unique_ptr<correction::CorrectionSet>& cset) {

    // obtain parameters from correctionlib
    double mean = cset->at("cb_params")->evaluate({abs(eta), nL, 0});
    double sigma = cset->at("cb_params")->evaluate({abs(eta), nL, 1});
    double n = cset->at("cb_params")->evaluate({abs(eta), nL, 2});
    double alpha = cset->at("cb_params")->evaluate({abs(eta), nL, 3});
    
    // instantiate CB and get random number following the CB
    CrystalBall cb(mean, sigma, alpha, n);
    TRandom3 rnd(time(0));
    double rndm = gRandom->Rndm();
    return cb.invcdf(rndm);
}


double get_std(double pt, double eta, float nL, unique_ptr<correction::CorrectionSet>& cset) {

    // obtain paramters from correctionlib
    double param_0 = cset->at("poly_params")->evaluate({abs(eta), nL, 0});
    double param_1 = cset->at("poly_params")->evaluate({abs(eta), nL, 1});
    double param_2 = cset->at("poly_params")->evaluate({abs(eta), nL, 2});

    // calculate value and return max(0, val)
    double sigma = param_0 + param_1 * pt + param_2 * pt*pt;
    if (sigma < 0) sigma = 0;
    return sigma; 
}


double get_k(double eta, string var, unique_ptr<correction::CorrectionSet>& cset) {

    // obtain parameters from correctionlib
    double k_data = cset->at("k_data")->evaluate({abs(eta), var});
    double k_mc = cset->at("k_mc")->evaluate({abs(eta), var});

    // calculate residual smearing factor
    // return 0 if smearing in MC already larger than in data
    double k = 0;
    if (k_mc < k_data) k = sqrt(k_data*k_data - k_mc*k_mc);
    return k;
}


double pt_resol(double pt, double eta, float nL, unique_ptr<correction::CorrectionSet>& cset, double low_pt_threshold = 26) {

    // load correction values
    double rndm = static_cast<double>(get_rndm(eta, nL, cset));
    double std = static_cast<double>(get_std(pt, eta, nL, cset));
    double k = static_cast<double>(get_k(eta, "nom", cset));

    // calculate corrected value and return original value if a parameter is nan
    double ptc = pt * ( 1 + k * std * rndm);
    if (isnan(ptc)) ptc = pt;
    if(ptc / pt > 2 || ptc / pt < 0.1 || ptc < 0 || pt < low_pt_threshold || pt > 200){
	ptc = pt;
    }
    return ptc;
}

double pt_resol_var(double pt_woresol, double pt_wresol, double eta, string updn, unique_ptr<correction::CorrectionSet>& cset){
    
    double k = static_cast<double>(get_k(eta, "nom", cset));

    if (k==0) return pt_wresol;

    double k_unc = cset->at("k_mc")->evaluate({abs(eta), "stat"});

    double std_x_rndm = (pt_wresol / pt_woresol - 1) / k;

    double pt_var = pt_wresol;

    if (updn=="up"){
        pt_var = pt_woresol * (1 + (k+k_unc) * std_x_rndm);
    }
    else if (updn=="dn"){
        pt_var = pt_woresol * (1 + (k-k_unc) * std_x_rndm);
    }
    else {
        cout << "ERROR: updn must be 'up' or 'dn'" << endl;
    }
    if (pt_var / pt_woresol > 2 || pt_var / pt_woresol < 0.1 || pt_var < 0)
        pt_var = pt_woresol; 

    return pt_var;
}

double pt_scale(bool is_data, double pt, double eta, double phi, int charge, unique_ptr<correction::CorrectionSet>& cset, double low_pt_threshold = 26) {
        
    // use right correction
    string dtmc = "mc";
    if (is_data) dtmc = "data";

    double a = cset->at("a_"+dtmc)->evaluate({eta, phi, "nom"});
    double m = cset->at("m_"+dtmc)->evaluate({eta, phi, "nom"});
    if(pt < low_pt_threshold)
	   return pt;

    return 1. / (m/pt + charge * a);
}

double pt_scale_var(double pt, double eta, double phi, int charge, string updn, unique_ptr<correction::CorrectionSet>& cset) {
        
    double stat_a = cset->at("a_mc")->evaluate({eta, phi, "stat"});
    double stat_m = cset->at("m_mc")->evaluate({eta, phi, "stat"});
    double stat_rho = cset->at("m_mc")->evaluate({eta, phi, "rho_stat"});

    double unc = pt*pt*sqrt(stat_m*stat_m / (pt*pt) + stat_a*stat_a + 2*charge*stat_rho*stat_m/pt*stat_a);

    double pt_var = pt;
    
    if (updn=="up"){
        pt_var = pt + unc;
    }
    else if (updn=="dn"){
        pt_var = pt - unc;
    }

    return pt_var;
}

}

