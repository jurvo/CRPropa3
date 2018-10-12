from crpropa import *
from pylab import *
import random as rd
import numpy as np
import math
HI=HadronicInteraction()
#Number of particles
particles = 100


#Candidate
c = Candidate()

#Radius of Sphere

r= 234*pc

#Steplength
steplength = 1

#Initialise Sourcelist
sl=SourceList()
origin = Vector3d(0)

#Create Sources
#Calculated by Integeral


for i in range (10):
    s=Source()
    s.add(SourcePosition(HI.Position(5*pc,50*pc)))
    s.add(SourceParticleType(nucleusId(1,1)))
    s.add(SourcePowerLawSpectrum(1*TeV, 1000*TeV, -2.0))
    s.add(SourceIsotropicEmission())
    sl.add(s)



#Create Observer Sphere
obs=Observer()
obs.add(ObserverLargeSphere(origin, r))


#Create output2
output2 = TextOutput("Test_Secondaries.txt")
output2.disableAll()
output2.enable(output2.SerialNumberColumn)
output2.enable(output2.CurrentIdColumn)
output2.enable(output2.CurrentDirectionColumn)
output2.enable(output2.CurrentEnergyColumn)
output2.enable(output2.CurrentPositionColumn)
output2.enable(output2.SourcePositionColumn)
output2.enable(output2.TrajectoryLengthColumn)
output2.setEnergyScale(TeV)
output2.setLengthScale(pc)

#Create Time Evolution Source
obs2=Observer()
obs2.add(ObserverTimeEvolution(0, 0* pc, 1))
obs2.add(ObserverNucleusVeto())
#obs2.setDeactivateOnDetection(False)
obs2.onDetection(output2)

#Create output2
output3 = TextOutput("Test_Traj.txt")
output3.disableAll()
output3.enable(output3.SerialNumberColumn)
output3.enable(output3.CurrentIdColumn)
output3.enable(output3.CurrentEnergyColumn)
output3.enable(output3.CurrentPositionColumn)
output3.enable(output3.SourcePositionColumn)
output3.enable(output3.SourceEnergyColumn)
output3.enable(output3.TrajectoryLengthColumn)
output3.setEnergyScale(TeV)
output3.setLengthScale(pc)

#Create Time Evolution Source
obs3=Observer()
obs3.add(ObserverTimeEvolution(0*pc, 10* pc, 23))
obs3.setDeactivateOnDetection(False)
obs3.onDetection(output3)

#Initiate turbulent magnetic field
MagField = MagneticFieldList()
randomSeed = 27
lMin=40. * pc/20 #(=lmax/20)
lMax=40.*pc #(1=lmax/(2*pi*R)) R_10Tev ca. 1.81 pc
l= turbulentCorrelationLength(lMin, lMax, -11./3.)
spacing=0.3*pc
vgrid = VectorGrid(origin, 512, spacing)

b=10*1e-6*gauss
initTurbulence(vgrid, b, lMin, lMax, -11./3., randomSeed)
TurbField = MagneticFieldGrid(vgrid)
MagField.addField(TurbField)

CM=CMZField()
JF=JF12Field()

MagField.addField(JF)

#ModuleList
m = ModuleList()
#Propagation
mod_Diffision=DiffusionSDE(MagField, 0.5, 0.1 * steplength * pc, 1* steplength * pc)
m.add(mod_Diffision)

#Interaction
m.add(HadronicInteraction(True, True, True))
m.add(PerformanceModule())
#m.add(EMInverseComptonScattering(Blackbody, True, 1))
#m.add(EMDoublePairProduction(Blackbody, True, 1))
#m.add(EMPairProduction(Blackbody, True, 1))
#m.add(EMTripletPairProduction(Blackbody, True, 1))
#m.add(ElasticScattering(Blackbody))
#m.add(ElectronPairProduction(Blackbody, True, 1))
#m.add(SynchrotronRadiation(TurbField, True, 1))
#m.add(PhotoDisintegration(Blackbody, True, 1))
#m.add(PhotoPionProduction(Blackbody, True, True, False, 1, False))

#Observer
m.add(obs)
m.add(obs2)
m.add(obs3)

m.add(MaximumTrajectoryLength(1*kpc))

m.setShowProgress(True)
m.run(sl, particles, True)

#Explicitly closing files.
#Important e.g. for the usage with jupyter notebooks

output2.close()
output3.close()







