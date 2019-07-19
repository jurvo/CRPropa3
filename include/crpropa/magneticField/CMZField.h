#ifndef CRPROPA_CMZFIELD_H
#define CRPROPA_CMZFIELD_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"

namespace crpropa {



class CMZField: public MagneticField {


public:
	CMZField();
    double r1(double r, double z, double a) const;
    double Bz1(double r, double z, double B1, double a, double L) const;
    double Br(double r, double z, double B1, double a, double L) const;
    double Bz(double r, double z, double B1, double a, double L) const;
    double g_s(double r,double z, double a, double Lp, double Hp, double p) const;
    double Bz12(double r, double phi, double z, double B1, double m, double a, double L, double Lp, double Hp, double p, double phi_s) const;
    double Br2(double r, double phi, double z, double B1, double m, double a, double L, double Lp, double Hp, double p, double phi_s) const;
    double Bz2(double r, double phi, double z, double B1, double m, double a, double L, double Lp, double Hp, double p, double phi_s) const;
    double Bphi(double r, double phi, double z, double B1, double m, double a, double L, double Lp, double Hp, double p, double phi_s) const;
    Vector3d TotalField(const Vector3d& pos) const;


};

} // namespace crpropa

#endif // CRPROPA_CMZFIELD_H
