#ifndef H_MUON_SCARE
#define H_MUON_SCARE

#include <cmath>
#include <memory>
#include <string>

#include "correction.hpp"

namespace scarekit {

struct CrystalBall{
    double pi=3.14159;
    double sqrtPiOver2=sqrt(pi/2.0);
    double sqrt2=sqrt(2.0);
    double m;
    double s;
    double a;
    double n;
    double B;
    double C;
    double D;
    double N;
    double NA;
    double Ns;
    double NC;
    double F;
    double G;
    double k;
    double cdfMa;
    double cdfPa;
    CrystalBall();
    CrystalBall(double mean, double sigma, double alpha, double n_);
    void init();
    double pdf(double x) const;
    double pdf(double x, double ks, double dm) const;
    double cdf(double x) const;
    double invcdf(double u) const;
};

double get_rndm(double eta, float nL, std::unique_ptr<correction::CorrectionSet>& cset);

double get_std(double pt, double eta, float nL, std::unique_ptr<correction::CorrectionSet>& cset);

double get_k(double eta, std::string var, std::unique_ptr<correction::CorrectionSet>& cset);

double pt_resol(double pt, double eta, float nL, std::unique_ptr<correction::CorrectionSet>& cset, double low_pt_threshold);

double pt_resol_var(double pt_woresol, double pt_wresol, double eta, std::string updn, std::unique_ptr<correction::CorrectionSet>& cset);

double pt_scale(bool is_data, double pt, double eta, double phi, int charge, std::unique_ptr<correction::CorrectionSet>& cset, double low_pt_threshold);

double pt_scale_var(double pt, double eta, double phi, int charge, std::string updn, std::unique_ptr<correction::CorrectionSet>& cset);

}

#endif
