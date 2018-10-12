#include "crpropa/magneticField/CMZField.h"
#include "crpropa/Units.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#include <iostream>

namespace crpropa {
//                                     Magnetic Field                                         
// Models are taken from Katia Ferriere and Philippe Terral 2014 "Analytical models of X-shape magnetic fields in galactic halos"
//                                                                                    
//   model C
CMZField::CMZField() {}
	
	
	
double CMZField::r1(double r, double z, double a) const{
    return r/(1+a*pow(z,2.));
}
double CMZField::Bz1(double r, double z, double B1, double a, double L) const{
    return B1*exp(-r1(r,z,a)/L);
}
double CMZField::Br(double r, double z, double B1, double a, double L) const{
    return 2*a*pow(r1(r,z,a),3.) *z/pow(r,2.) *Bz1(r,z,B1,a,L);
}
double CMZField::Bz(double r, double z, double B1, double a, double L) const{
    return pow(r1(r,z,a),2.) /pow(r,2.) *Bz1(r,z,B1,a,L);
}
///model C + azimuthal component
double CMZField::g_s(double r,double z, double a, double Lp, double Hp, double p) const{
    double f1= 1+ pow(fabs(z)/Hp,2.);
    double f2=r/Lp;
    return cos(p* M_PI/180)/sin(p* M_PI/180) * log(1-exp(-f2)) / f1;
}
double CMZField::Bz12(double r, double phi, double z, double B1, double m, double a, double L, double Lp, double Hp, double p, double phi_s) const{
    return Bz1(r,z,B1,a,L)* cos(M_PI/180 *m*(phi-g_s(r,z,a,Lp,Hp,p)-phi_s));// from the paper Philippe Terral and Katia Ferriere 2017 "Constraints from Faraday rotation on the magnetic field structure in the Galactic halo"
}
double CMZField::Br2(double r, double phi, double z, double B1, double m, double a, double L, double Lp, double Hp, double p, double phi_s) const{
    return 2*a* pow(r1(r,z,a),3.) * z/pow(r,2.) *Bz12(r,phi,z,B1,m,a,L,Lp,Hp,p,phi_s);
}
double CMZField::Bz2(double r, double phi, double z, double B1, double m, double a, double L, double Lp, double Hp, double p, double phi_s) const{
    return pow(r1(r,z,a),2.)/pow(r,2.) *Bz12(r,phi,z,B1,m,a,L,Lp,Hp,p,phi_s);
}
double CMZField::Bphi(double r, double phi, double z, double B1, double m, double a, double L, double Lp, double Hp, double p, double phi_s) const{//  azimuthal component   from the paper Philippe Terral and Katia Ferriere 2017 "Constraints from Faraday rotation on the magnetic field structure in the Galactic halo"
    double f1= 1.+ pow(fabs(z)/Hp,2.);
    double f2=r/L;
    return cos(p*M_PI / 180)/sin(p*M_PI / 180) *(f2*exp(-f2)/(f1*(1-exp(-f2))) * Br2(r,phi,z,B1,m,a,L,Lp,Hp,p,phi_s) - 2./pow(f1,2.) *z/pow(Hp,2.) *r*log(1.-exp(-f2))*Bz2(r,phi,z,B1,m,a,L,Lp,Hp,p,phi_s));//the first term is cot(p)=cos(p)/sin(p) and p in deg  cot is 0 at 90 deg, maximum at 0 deg and minimum at 180 deg
}
//~ ///end model C + azimuthal component
//  Parameter identification:
// r: radius in zylindrical coordinates
// z= z component in zylindircal coordintes
// phi: phi component in zylindircal coordinates
// B1: maximum magnetic field which can be found at z,r=0,0
// m: azimuthal wavenumber m=0 for axissymmetric, m=1 for bisymmetric and m=2 quadrosymmetric, ... field lines
// a: stricly positive free parameters governing the opening of field lines away from the z-axis
// L: exponential scale length length of the cloud radius Longitudinal extent
// Lp:scale length of the winding function set to a fraction of the cloud radius  Longitudes of peaks mainly at high |l|
// H: exponential scale high length of the cloud radius Latitudinal extent
// Hp: exponential scale high  set to a fraction of the cloud vertical size Longitude of peaks mainly at hight |b|
// p: pitch angle at the origin, i.e. the angle between the horizontal projection of a field line and the local azimuthal direction  (i.e., 10° or 20°) if you want the field to be tightly wound up, and hence predominantly azimuthal
// phi_s: orientation angle of the azimuthal pattern if g_s(r1(r,z,a),0,a,Lp,Hp,p), azimuthal angle at infinity of the crest surface if g_s(r,z,a,Lp,Hp,p)
//                                                                                     

Vector3d CMZField::TotalField(const Vector3d& pos) const {
	Vector3d b(0.);
    double pi= M_PI;
	double r = sqrt(pos.x * pos.x + pos.y * pos.y); //  transformation from zylindrical coordinates into cartesian
	double d = pos.getR(); // distance to galactic center
	double phi = pos.getPhi(); // azimuth
	double sinPhi = sin(phi);
	double cosPhi = cos(phi);
	// azimuthal component in dense clouds
    double m=0;//Parameter of X-Shape model m=0 for axissymmetric, m=1 for bisymmetric and m=2 quadrosymmetric, ... field lines
    double p=10.;
    double ka=3./40.;
    double kb=20.;
    double kc=15.;
    double B1=20./9.*1.e-3/(cos(p)/sin(p));//Parameter of X-Shape model, maximum B-field strenght
    double B2=5.e-5;
    double B3=1.e-3;//Parameter of X-Shape model, maximum B-field strenght
    //A=SgrC 
    double x_mid=0;//  x midpoint of A
    double y_mid=sin(359.45*pi/180)*8.5*kpc;//  y midpoint of A
    double z_mid=sin(-0.11*pi/180)*8.5*kpc;//  z midpoint of A
	double x=pos.x;
	double y=pos.y;
	double z=pos.z;
    double xx=x-x_mid;// shifted coordinates
    double yy=y-y_mid;// shifted coordinates
    double zz=z-z_mid;// shiftes coordinates
    double rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    double pphi=asin(yy/rr)*180/M_PI;// shifted coordinates
    double R = 1.7*pc;// Radius of A
    double L = R;// Parameter of X-Shape model
    double a=1./pow(L,2.) *ka;// Parameter of X-Shape model
    double L2=L/3.;
    double a2=a*4.;
    double Lp=R/kb;
    double Hp=Lp*kc;
    double phi_s=0;
    b.y +=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x +=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    // A=G0.253+0.016 Dust Ridge A
    x_mid=0;//  x midpoint of A
    y_mid=sin(0.253*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(-0.016*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180./pi;// shifted coordinates
    R=2.9*pc;// Radius of A
    L=R;// Parameter of X-Shape model
    a=1/pow(L,2.) *ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B1=5.5e-3;
    b.y+=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x+=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    //                                                             A=Sgr B1-off
    x_mid=0;//  x midpoint of A
    y_mid=sin(0.48*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(-0.0*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    R=3.6*pc;// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.) *ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B1=20./9.*1.e-3/(cos(p)/sin(p));
    b.y+=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x+=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    //                                                             A=Dust Ridge B
    x_mid=0;//  x midpoint of A
    y_mid=sin(0.34*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(0.055*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    R=1.9*pc;// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.) *ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    b.y+=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x+=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    //                                                             A=Dust Ridge C
    x_mid=0;//  x midpoint of A
    y_mid=sin(0.38*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(0.05*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    R=1.9*pc;// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.) *ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    b.y+=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x+=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    //                                                             A=Dust Ridge D
    x_mid=0;//  x midpoint of A
    y_mid=sin(0.41*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(0.05*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    R=3.2*pc;// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.) *ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    b.y+=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x+=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    //                                                             A=Dust Ridge E
    x_mid=0;//  x midpoint of A
    y_mid=sin(0.478*pi/180)*8.5*kpc;// y midpoint of A
    z_mid=sin(0.005*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    R=2.4*pc;// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.) *ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    b.y+=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x+=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    //                                                             A=Dust Ridge F
    x_mid=0;//  x midpoint of A
    y_mid=sin(0.496*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(0.02*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    R=2.0*pc;// Radius of A
    L=R;//Parameter of X-Shape model
    a=1./pow(L,2.) *ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    b.y+=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x+=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    //                                                             Sgr D
    x_mid=0;//  x midpoint of A
    y_mid=sin(1.12*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(-0.07*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy); //shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    //a=sin(7.5*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    //bbnsin(6.4*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    R=2*pc;// sqrt((a*b)/pi) Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.) *ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    b.y+=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x+=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    //                                                             Sgr B2
    x_mid=0;//  x midpoint of A
    y_mid=sin(0.66*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(-0.04*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid ;//shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    R=7.1*pc;// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.) *ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B1=20./9.*500.e-6/(cos(p)/sin(p));
    b.y+=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x+=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    //                                                                                    
    //                                                             A=Inner R=10pc + poloid
    x_mid=0;//  x midpoint of A
    y_mid=sin(0*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(0.0*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    R=10*pc;// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.)*ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B1=20./9.*3.e-3/(cos(p)/sin(p));
    b.y+=sin(pi/180*pphi)*Br2(rr,pphi,zz,B2,m,a2,L2,Lp,Hp,p,phi_s)+cos(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.x+=cos(pi/180*pphi)*Br2(rr,pphi,zz,B2,m,a2,L2,Lp,Hp,p,phi_s)-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B1,m,a,L,Lp,Hp,p,phi_s);
    b.z+=Bz2(rr,pphi,zz,B2,m,a2,L2,Lp,Hp,p,phi_s);
    std::cout <<"bx="<< b.x << ";by=" << b.y << ";bz" << b.z << ";x="<< x/pc << ";y="<< y/pc <<";z="<< z/pc << ";pc="<< pc << std::endl;
    //POLOIDAL COMPONENT
    ////                                                                                    
    //                                                             A=Intercloudmedium
    x_mid=0;//  x midpoint of A
    y_mid=sin(0*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(-0.0*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    R=200.*pc;// Radius of A
    L=R*0.5;// Parameter of X-Shape model
    a=1./(3000*pc*pc);// Parameter of X-Shape model
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0*pi/180;
    b.y+=sin(pi/180*pphi)*Br(rr,zz,B2,a,L);
    b.x+=cos(pi/180*pphi)*Br(rr,zz,B2,a,L);
    b.z+=Bz(rr,zz,B2,a,L);
    //                                                             
    //    Filaments->Polodial component up to 1 mG
    //     Inner10=SgrA is already considered
    //                           A=SgrC
    x_mid=0;//  x midpoint of A
    y_mid=sin(359.45*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(-0.01*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    double arcmin=1./60.;
    double aa=sin(11.1*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    double bb=sin(0.7*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    R=sqrt((aa*bb)/pi);// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.)*ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B3=1.e-4;
    b.y+=sin(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.x+=cos(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.z+=Bz2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    //                           A=G359.15-0.2 The Snake
    x_mid=0;//  x midpoint of A
    y_mid=sin(359.15*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(-0.17*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    aa=sin(5.2*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    bb=sin(0.9*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    R=sqrt((aa*bb)/pi);// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.)*ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B3=88.e-6;
    b.y+=sin(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.x+=cos(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.z+=Bz2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    //                           A=G359.54+0.18 Nonthermal Filament
    x_mid=0;//  x midpoint of A
    y_mid=sin(359.54*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(0.17*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    aa=sin(6.1*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    bb=sin(1.1*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    R=sqrt((aa*bb)/pi);// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.)*ka ;//Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B3=1.e-3;
    b.y+=sin(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.x+=cos(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.z+=Bz2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    //                           A=G359.79 +17 Nonthermal Filament
    x_mid=0;//  x midpoint of A
    y_mid=sin(359.79*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(0.16*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    aa=sin(6.5*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    bb=sin(1.4*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    R=sqrt((aa*bb)/pi);// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.)*ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B3=1.e-3;
    b.y+=sin(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.x+=cos(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.z+=Bz2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    //                           A=G359.85+0.47  Nonthermal Filament The Pelican  is not poloidal but azimuthal
    x_mid=0;//  x midpoint of A
    y_mid=sin(359.85*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(0.47*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    aa=sin(4.6*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    bb=sin(0.9*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    R=sqrt((aa*bb)/pi);// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.)*ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B3=70.e-6*20./9.*1./(cos(p)/sin(p));
    b.y+=cos(pi/180*pphi)*Bphi(rr,pphi,zz,B3,m,a,L,Lp,Hp,p,phi_s);
    b.x+=-sin(pi/180*pphi)*Bphi(rr,pphi,zz,B3,m,a,L,Lp,Hp,p,phi_s);
    //                           A=G359.91 -1.03 Possible Nonthermal Filament
    x_mid=0;//  x midpoint of A
    y_mid=sin(359.91*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(-1.03*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    aa=sin(2.3*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    bb=sin(0.6*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    R=sqrt((aa*bb)/pi);// Radius of A
    L=R ;//Parameter of X-Shape model
    a=1./pow(L,2.)*ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B3=1.e-3;
    b.y+=sin(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.x+=cos(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.z+=Bz2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    //                           A=G359.96 +0.09  Nonthermal Filament Southern Thread
    x_mid=0;//  x midpoint of A
    y_mid=sin(359.96*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(0.11*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    aa=sin(11.6*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    bb=sin(0.7*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    R=sqrt((aa*bb)/pi);// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.)*ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B3=1e-4;
    b.y+=sin(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.x+=cos(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.z+=Bz2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    //                           A=G0.09 +0.17  Nonthermal Filament Northern thread
    x_mid=0;//  x midpoint of A
    y_mid=sin(0.09*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(0.17*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    aa=sin(11.9*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    bb=sin(0.9*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    R=sqrt((aa*bb)/pi);// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.)*ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B3=140.e-6;
    b.y+=sin(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.x+=cos(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.z+=Bz2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    //                           A=GCRA  The Radio Arc
    x_mid=0;//  x midpoint of A
    y_mid=sin(0.18*pi/180)*8.5*kpc;//  y midpoint of A
    z_mid=sin(-0.07*pi/180)*8.5*kpc;//  z midpoint of A
    xx=x-x_mid;// shifted coordinates
    yy=y-y_mid;// shifted coordinates
    zz=z-z_mid;// shiftes coordinates
    rr=sqrt(xx*xx+yy*yy);// shifted coordinates
    pphi=asin(yy/rr)*180/pi;// shifted coordinates
    aa=sin(28.5*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    bb=sin(4.0*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    R=sqrt((aa*bb)/pi);// Radius of A
    L=R;// Parameter of X-Shape model
    a=1./pow(L,2.)*ka;// Parameter of X-Shape model
    L2=L/3.;
    a2=a*4.;
    Lp=R/kb;
    Hp=Lp*kc;
    phi_s=0;
    B3=1.e-3;
    b.y+=sin(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.x+=cos(pi/180*pphi)*Br2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    b.z+=Bz2(rr,pphi,zz,B3,m,a2,L2,Lp,Hp,p,phi_s);
    //                                                     
	return b;//1e.4 due to the conversion from Gauss to Tesla

}

} // namespace crpropa
