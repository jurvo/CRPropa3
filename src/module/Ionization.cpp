#include "crpropa/module/Ionization.h"
#include "crpropa/Common.h"

using namespace crpropa;

Ionization::Ionization(ref_ptr<Density> density) : 
    density(density) {

}



double Ionization::sigma_p_ion(double E) {
    // fiting parameter
    double A = 0.71; 
    double B = 1.63;
    double C = 0.51;
    double D = 1.24;
    
    double x = mass_electron * E / mass_proton / IH;

    double factor = 4 * M_PI * pow_integer<2>(bohr_radius);
    double sigmaL = factor * C * pow(x, D);
    double sigmaH = factor * (A * log(1 + x) + B) / x;

    return 1 / (1 / sigmaL + 1 / sigmaH);
}

double Ionization::sigma_p_ec(double E) {
    double A = 1.044; 
    double B = 2.88;
    double C = 0.016; 
    double D = 0.136;
    double F = 5.86;
    int N = 2; // number of electrons
    
    double x = E / IH;

    double sigma = 4 * M_PI * pow_integer<2>(bohr_radius) * A * N * pow_integer<2>(IH / IH2 * x);
    sigma /= C + pow(x, B) + D * pow(x, F);
    return sigma;
}

double Ionization::sigma_p_diss(double E) {
    return sigma_e_diss(mass_electron * E / mass_proton);
}

double Ionization::sigma_p_doub(double E) {
    return sigma_e_doub(mass_electron * E / mass_proton);
}

double Ionization::sigma_e_ion(double E) {
    double t = E / IH2;
    int N = 2; // number of electrons

    double n = 2.4;
    double A1 = 0.74;
    double A2 = 0.87;
    double A3 = -0.60;

    double F = (1 - pow(t, 1 - n)) / (n - 1) * pow(2 / (1 + t), n/2.) * (1 - pow(t, 1 - n / 2)) / (n - 2);
    double G = (A1 * log(t) + A2 + A3 / t) / t;

    double sigma = 4 * M_PI * pow_integer<2>(bohr_radius) * N * pow_integer<2>(IH / IH2) * F * G;
    return sigma;
}

double Ionization::sigma_e_diss(double E) {
    if(E < 18.1 * eV) 
        return 0;
    
    double logE = log(E / eV);

    double logSigma = 0;
    for (int i = 0; i < 6; i++)
    {
        logSigma += anDiss[i] * pow_integer<i>(logE);
    }

    double sigma = exp(logSigma) * 1e-18 * cm * cm;
    return sigma;
}

double Ionization::sigma_e_doub(double E) {
    if(E < 51 * eV) 
        return 0;

    double logE = log(E / eV);

    double logSigma = 0;
    for (int i = 0; i < 6; i++)
    {
        logSigma += anDoub[i] * pow_integer<i>(logE);
    }

    double sigma = exp(logSigma) * 1e-18 * cm * cm;
    return sigma;
}