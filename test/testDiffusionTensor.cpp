#include "crpropa/DiffusionTensor.h"
#include "crpropa/magneticField/turbulentField/SimpleGridTurbulence.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/Grid.h"

#include "gtest/gtest.h"

using namespace crpropa;

TEST(testQLTDiffusion, TestGetSetFunction){
    QLTDiffusion* Tens = new QLTDiffusion();

    // check default values 
    EXPECT_DOUBLE_EQ(Tens-> getEpsilon(), 0.1);
    EXPECT_DOUBLE_EQ(Tens-> getKappa0(), 6.1e24);
    EXPECT_DOUBLE_EQ(Tens-> getAlpha(), 1./3.);

    // set new values and check them
    Tens -> setEpsilon(0.3);
    EXPECT_DOUBLE_EQ(Tens-> getEpsilon(), 0.3);
    Tens -> setKappa0(1e5);
    EXPECT_DOUBLE_EQ(Tens-> getKappa0(), 1e5);
    Tens -> setAlpha(0.2);
    EXPECT_DOUBLE_EQ(Tens-> getAlpha(), 0.2);
}

TEST(testQLTDiffusion, TestForEnergyScaling){
    QLTDiffusion* Tens = new QLTDiffusion();

    int id = 1000010010; // protons

    // checking at norm Value 4 GeV
    double energy = 4*GeV;
    Candidate cand(id, energy);
    double norm = Tens-> getKappa0();
    EXPECT_DOUBLE_EQ(Tens-> getKappaParallel(&cand), norm);
    EXPECT_DOUBLE_EQ(Tens-> getKappaPerpendicular(&cand), 0.1*norm); // with epsilon = 0.1 as default
    EXPECT_DOUBLE_EQ(Tens-> getKappaPerpendicular2(&cand), 0.1*norm);

    // checking for another energy
    energy = 15*GeV;
    cand.current.setEnergy(energy);
    EXPECT_NEAR(Tens-> getKappaParallel(&cand), 9.47706e24, 1e19);
    EXPECT_NEAR(Tens-> getKappaPerpendicular(&cand), 9.47706e23, 1e18);
    EXPECT_NEAR(Tens-> getKappaPerpendicular2(&cand), 9.47706e23, 1e18);

    // checking energy far away from norming
    energy = 15*PeV;
    cand.current.setEnergy(energy);
    EXPECT_NEAR(Tens-> getKappaParallel(&cand), 9.47706e26, 1e21);
    EXPECT_NEAR(Tens-> getKappaPerpendicular(&cand),  9.47706e25, 1e20);
    EXPECT_NEAR(Tens-> getKappaPerpendicular2(&cand),  9.47706e25, 1e20);
}

TEST(testQLTTurbulence, testQLTTurbulence){
    // background field
    double B = 10;  // field strenght
    UniformMagneticField* background = new UniformMagneticField(Vector3d(B, 0, 0));

    // turbulent field
    double b = 1;
    unsigned int seed = 42; // random number
    auto spectrum = SimpleTurbulenceSpectrum(b, 30*pc, 90*pc);
    auto gp = GridProperties(Vector3d(0.), 10, 10*pc);
    auto turbulentField = new SimpleGridTurbulence(spectrum, gp, seed);
    // diffusion tensor
    QLTTurbulent* Tens = new QLTTurbulent(background, turbulentField);
 
    // check position for turbulence scaling
    Vector3d pos1(-1*kpc, 0., -20*pc);
    double locTurb = turbulentField -> getField(pos1).getR() / B;
    double normTurb = turbulentField->getField(Vector3d(-8.5*kpc, 0., 0.)).getR()/B; // norm Value for the turbulence at earth

    // check default values
    EXPECT_DOUBLE_EQ(Tens ->getKappa0(), 6.1e24);
    EXPECT_DOUBLE_EQ(Tens ->getAlphaPara(), 1./3.);
    EXPECT_DOUBLE_EQ(Tens ->getAlphaPerp(), 1./3.);
    EXPECT_DOUBLE_EQ(Tens ->getNormTurbulence(), normTurb);
    
    normTurb = turbulentField->getField(Vector3d(0.)).getR()/B;
    Tens ->normToEarthPosition(Vector3d(0.));
    EXPECT_DOUBLE_EQ(Tens ->getNormTurbulence(), normTurb);

    // Candidate for checking
    int id = 1000010010; // protons
    double energy = 4*GeV; // checking at scaling energy -> no energy dependence 
    Candidate cand(id, energy);
    cand.current.setPosition(pos1);

    // check values
    double kappaPara = 6.1e24 * pow(locTurb/normTurb, -2);
    double kappaPerp = 6.1e24 * pow(locTurb*normTurb, 2);
    EXPECT_DOUBLE_EQ(Tens -> getKappaParallel(&cand), kappaPara);
    EXPECT_DOUBLE_EQ(Tens -> getKappaPerpendicular(&cand), kappaPerp);
    EXPECT_DOUBLE_EQ(Tens -> getKappaPerpendicular2(&cand), kappaPerp);

    // check energy scaling
    cand.current.setEnergy(4*PeV);  
    EXPECT_DOUBLE_EQ(Tens -> getKappaParallel(&cand), 100*kappaPara);
    EXPECT_DOUBLE_EQ(Tens -> getKappaPerpendicular(&cand), 100*kappaPerp);
    EXPECT_DOUBLE_EQ(Tens -> getKappaPerpendicular2(&cand), 100*kappaPerp);

    // check set function
    Tens ->setKappa0(1e15);
    EXPECT_DOUBLE_EQ(Tens ->getKappa0(), 1e15);
    Tens ->setAlphaPara(0.3);
    EXPECT_DOUBLE_EQ(Tens ->getAlphaPara(), 0.3);
    Tens ->setAlphaPerp(0.7);
    EXPECT_DOUBLE_EQ(Tens ->getAlphaPerp(), 0.7);
    Tens ->setAlpha(0.5);
    EXPECT_DOUBLE_EQ(Tens ->getAlphaPara(), 0.5);
    EXPECT_DOUBLE_EQ(Tens ->getAlphaPerp(), 0.5);
    Tens ->setNormTurbulence(0.01);
    EXPECT_DOUBLE_EQ(Tens ->getNormTurbulence(), 0.01);
}

TEST(testQLTRigidity, testQLTRigidity){
    // background field
    double B = 10;  // field strenght
    UniformMagneticField* background = new UniformMagneticField(Vector3d(B, 0, 0));

    // turbulent field
    double b = 1;
    unsigned int seed = 42; // random number
    auto spectrum = SimpleTurbulenceSpectrum(b, 30*pc, 90*pc);
    auto gp = GridProperties(Vector3d(0.), 10, 10*pc);
    auto turbulentField = new SimpleGridTurbulence(spectrum, gp, seed);

    QLTRigidity* Tens = new QLTRigidity(background, turbulentField);

    // candidate for testparticle
    Vector3d pos1(0.);
    Vector3d posEarth(-8.5*kpc, 0., 0.);
    int id = 1000010010; // protons
    Candidate cand(id, 4*GeV, pos1);

    // check default values
    EXPECT_DOUBLE_EQ(Tens -> getKappa0(), 6.1e24);
    EXPECT_DOUBLE_EQ(Tens -> getAlphaPara(), 1./3.);
    EXPECT_DOUBLE_EQ(Tens -> getAlphaPerp(), 1./3.);
    double normEta = turbulentField -> getField(posEarth).getR()/B;
    EXPECT_DOUBLE_EQ(Tens -> getNormEta(), normEta);
    double normRho = 4e9*volt / B / c_light / turbulentField-> getCorrelationLength();
    EXPECT_DOUBLE_EQ(Tens -> getNormRho(), normRho);

    //EXPECT_DOUBLE_EQ(Tens -> getKappaParallel(&cand), 6.1e24);
    //EXPECT_DOUBLE_EQ(Tens -> getKappaPerpendicular(&cand), pow(normEta,4)*6.1e24);
    //EXPECT_DOUBLE_EQ(Tens -> getKappaPerpendicular2(&cand), pow(normEta,4)*6.1e24);

    // energy scaling 
    double kPara = Tens -> getKappaParallel(&cand);
    double kPerp = Tens -> getKappaPerpendicular(&cand);
    cand.current.setEnergy(4*PeV); // 1e6 in energy -> 1e2 for kappa (scales with E^1/3)
    EXPECT_DOUBLE_EQ(Tens -> getKappaParallel(&cand), 100*kPara);
    EXPECT_DOUBLE_EQ(Tens -> getKappaPerpendicular(&cand), 100*kPerp);

    // test set function
    Tens -> setKappa0(100);
    EXPECT_DOUBLE_EQ(Tens -> getKappa0(), 100);
    Tens -> setAlphaPara(0.2);
    EXPECT_DOUBLE_EQ(Tens -> getAlphaPara(), 0.2);
    Tens -> setAlphaPerp(0.4);
    EXPECT_DOUBLE_EQ(Tens -> getAlphaPerp(), 0.4);
    Tens -> setAlpha(0.5);
    EXPECT_DOUBLE_EQ(Tens -> getAlphaPara(), 0.5);
    EXPECT_DOUBLE_EQ(Tens -> getAlphaPerp(), 0.5);
    Tens -> setNormEta(0.2);
    EXPECT_DOUBLE_EQ(Tens -> getNormEta(), 0.2);
    Tens -> setNormRho(0.8);
    EXPECT_DOUBLE_EQ(Tens -> getNormRho(), 0.8);
    Vector3d normPos(1,5,8);
    Tens -> normToPosition(normPos);
    Vector3d tensNormPos = Tens -> getNormPos();
    EXPECT_DOUBLE_EQ(tensNormPos.x, normPos.x);
    EXPECT_DOUBLE_EQ(tensNormPos.y, normPos.y);
    EXPECT_DOUBLE_EQ(tensNormPos.z, normPos.z);
}