#ifndef CRPROPA_CMZFIELD_H
#define CRPROPA_CMZFIELD_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"

namespace crpropa {



class CMZField: public MagneticField {


public:
	CMZField();
    double scale(double x, double d) const;
    double Br(double r, double z, double B1, double a, double L) const;
    double Bz(double r, double z, double B1, double a, double L) const;
    double BrAz(double r, double phi, double z, double m, double B1, double eta, double R,double rr) const;
    double BphiAz(double r, double phi, double z, double m, double B1, double eta, double R, double rr) const;
    double PHI(double x, double y) const;
    double BxAz(double x, double y, double z, double m, double B1, double eta, double R, double rr) const;
    double ByAz(double x, double y, double z, double m, double B1, double eta, double R, double rr) const;
    double BzAz(double x,double y, double z, double m, double B1, double eta, double R, double rr) const;
    double ByPol(double x,double y, double z, double B1, double a1, double a2, double eta) const;
    double BxPol(double x,double y, double z, double B1, double a1, double a2, double eta) const;
    double BzPol(double x, double y, double z, double B1, double a1, double a2, double eta) const;
    
    Vector3d MCField(const Vector3d& pos) const;
    Vector3d ICField(const Vector3d& pos) const;
    Vector3d NTFField(const Vector3d& pos) const;
    Vector3d RadioArc(const Vector3d& pos) const;
    Vector3d getField(const Vector3d& pos) const;
};

} // namespace crpropa

#endif // CRPROPA_CMZFIELD_H
