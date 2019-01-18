//#include "crpropa/massDistribution/Massdistribution.h"
#include "crpropa/massDistribution/Cordes.h"
#include "crpropa/massDistribution/Ferriere.h"
#include "crpropa/massDistribution/Nakanishi.h"
#include "crpropa/Units.h"
#include "crpropa/massDistribution/Density.h"

#include "gtest/gtest.h"

#include <stdexcept>
#include <cmath>
#include <string>

namespace crpropa {

TEST(testConstantDensity, SimpleTest) {
	//test ConstantDensity in all types and in total density (output)
	ConstantDensity n(2/ccm,3/ccm, 2/ccm);
	Vector3d p(1*pc,2*pc,1*kpc); 	// random position for testing density
	EXPECT_DOUBLE_EQ(n.getHIDensity(p), 2e6);	// density output in m^-3
	EXPECT_DOUBLE_EQ(n.getHIIDensity(p), 3e6);
	EXPECT_DOUBLE_EQ( n.getH2Density(p), 2e6);
	EXPECT_DOUBLE_EQ(n.getDensity(p), 7e6);		// total density 2+3+2 = 7 (/ccm)
	EXPECT_DOUBLE_EQ(n.getNucleonDensity(p),9e6);	// nucleon density 2+3+2*2 = 9 (/ccm) factor 2 for molecular hydrogen

	//set density number to 500
	n.setHI(500.);
	n.setHII(500.);
	n.setH2(500.);

	//check if output is changed to 500 in types 
	EXPECT_DOUBLE_EQ(n.getHIDensity(p), 500.);
	EXPECT_DOUBLE_EQ(n.getHIIDensity(p), 500.);
	EXPECT_DOUBLE_EQ(n.getH2Density(p), 500.);  

}

/*
TEST(testDensityList, SimpleTest) {

	DensityList MS;
	MS.addDensity(new ConstantDensity(1,1,2));	//sum 4
	MS.addDensity(new ConstantDensity(2,3,1));	//sum 6

	Vector3d p(50*pc,10*pc,-30*pc);	//random position for testing density
	EXPECT_DOUBLE_EQ(MS.getHIDensity(p),3);
	EXPECT_DOUBLE_EQ(MS.getHIIDensity(p),4);
	EXPECT_DOUBLE_EQ(MS.getH2Density(p),3);
	EXPECT_DOUBLE_EQ(MS.getDensity(p),10);	//sum of sums
	EXPECT_DOUBLE_EQ(MS.getNucleonDensity(p),13); 	// 3+4+2*3 factor 2 for molecular hydrogen

}
*/

TEST(testCordes, checkValueAtCertainPoints) {

	Cordes n;

	Vector3d p(3.1*kpc,2.9*kpc,-30*pc);	//position for testing density

	EXPECT_NEAR(n.getHIIDensity(p), 184500.,1);	// output in m^-3 ; uncertainty of 1e-6 cm^-3
	EXPECT_NEAR(n.getDensity(p), 184500.,1);
	EXPECT_NEAR(n.getNucleonDensity(p), 184500,1);	// only HII component -> no differenz between density and nucleon density
	p.z=30*pc;			// invariant density for +/- z
	EXPECT_NEAR(n.getDensity(p),184500,1);
}

TEST(testNakanishi, checkValueAtCertainPoints) {

	Nakanishi n;
	//first position for testing density
	Vector3d p(4*kpc,-2.5*kpc,-0.85*kpc);

	EXPECT_NEAR(n.getHIDensity(p),914,1); 	//uncertainc 1e-6 cm^-3
	EXPECT_NEAR(n.getH2Density(p),0,1);

	//testing total Density
	EXPECT_NEAR(n.getDensity(p),914,2); //double uncertaincy for both type á 1cm^-3
	EXPECT_NEAR(n.getNucleonDensity(p),914,2);	// 914 + 0*2	factor 2 for molecular hydrogen

	//second position for testing density
	p = Vector3d(50*pc,100*pc,10*pc);

	EXPECT_NEAR(n.getHIDensity(p),540867,1);
	EXPECT_NEAR(n.getH2Density(p),10335137,1);

	//testing total Density
	EXPECT_NEAR(n.getDensity(p),10876004,2);
	EXPECT_NEAR(n.getNucleonDensity(p),21211141,2);	// factor 2 in molecular hydrogen
}

TEST(testFerriere, checkValueAtCertainPoints) {
	Ferriere n;

	//testing density in inner Ring (R <= 3*kpc)
	Vector3d p(60*pc,-60*pc,-20*pc);	//testing position in region of CMZ
	//testing density
	EXPECT_NEAR(n.getHIDensity(p),6237723,1); 	//uncertaincy 1e-6 cm^-3
	EXPECT_NEAR(n.getH2Density(p),35484825,1);
	EXPECT_NEAR(n.getHIIDensity(p),6243793,1);
	EXPECT_NEAR(n.getDensity(p),47966341,1);
	EXPECT_NEAR(n.getNucleonDensity(p),83451166,2);		//factor 2 in molecular hydrogen; double uncertaincy

	Vector3d p2(-500*pc,-900*pc,35*pc);	//testing position in region of the DISK
	EXPECT_NEAR(n.getHIIDensity(p2),48190,1);
	EXPECT_NEAR(n.getHIDensity(p2),5,1);
	EXPECT_NEAR(n.getH2Density(p2),0,1);
	EXPECT_NEAR(n.getDensity(p2),48195,1);
	EXPECT_NEAR(n.getNucleonDensity(p2),48195,1);	// no H2 component -> no difference between density and nucleon-density

	//testing the outer region R>3kpc
	Vector3d p3(5*kpc,4*kpc,-29*pc);	//testing position with 3kpc < R < R_sun
	EXPECT_NEAR(n.getHIDensity(p3),540607,1);
	EXPECT_NEAR(n.getHIIDensity(p3),66495 ,1);
	EXPECT_NEAR(n.getH2Density(p3),2492685,1);
	EXPECT_NEAR(n.getDensity(p3),3099787,1);
	EXPECT_NEAR(n.getNucleonDensity(p3), 5592472,1);

	Vector3d p4(10*kpc,2*kpc,50*pc);	//testing position with R > R_sun
	EXPECT_NEAR(n.getHIDensity(p4),431294,1);
	EXPECT_NEAR(n.getHIIDensity(p4),22109,1);
	EXPECT_NEAR(n.getH2Density(p4),54099,1);
	EXPECT_NEAR(n.getDensity(p4),507502,1);
	EXPECT_NEAR(n.getNucleonDensity(p4),561601,1);
}
	
} //namespace crpropa
