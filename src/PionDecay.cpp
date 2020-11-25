#include "crpropa/module/PionDecay.h"
#include "crpropa/Units.h"

#include <cmath>
#include <string>
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <functional>

namespace crpropa {
    
    struct particle {
        int ID;
        double px;
        double py;
        double pz;
        double E;
        double m0;
    };

    static int LBARP[49] = { 1,   3,   2,   5,   4,   6,   8,   7,  10,   9,
                         11,  12, -13, -14,  16,  15,  18,  17,  13,  14,
                         22,  21,  23,  24,  26,  25,  27,  29,  28,  31,
                         30,  32,  33, -34, -35, -36, -37, -38, -39, -40,
                        -41, -42, -43, -44, -45, -46, -47, -48, -49 };

    static double AM[49] = { 0., 0.511e-3, 0.511e-3, 0.10566, 0.10566,
                            0.13497,  0.13957,  0.13957, 0.49365, 0.49365,
                            0.49767,  0.49767,  0.93827, 0.93957,      0.,
                                 0.,       0.,       0., 0.93827, 0.93957,
                            0.49767,  0.49767,  0.54880, 0.95750, 0.76830,
                            0.76830,  0.76860,  0.89183, 0.89183, 0.89610,
                            0.89610,  0.78195,  1.01941, 1.18937, 1.19255,
                            1.19743,  1.31490,  1.32132, 1.11563, 1.23100,
                            1.23500,  1.23400,  1.23300, 1.38280, 1.38370,
                            1.38720,  1.53180,  1.53500, 1.67243 };

    static int KDEC[612] = { 3,  1, 15,  2, 18,  0,  3,  1, 16,  3, 17,  0,
                            2,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,
                            2,  0,  4, 17,  0,  0,  2,  0,  5, 18,  0,  0,
                            2,  0,  4, 17,  0,  0,  2,  0,  7,  6,  0,  0,
                            3,  0,  7,  7,  8,  0,  3,  0,  7,  6,  6,  0,
                            3,  1, 17,  4,  6,  0,  3,  1, 15,  2,  6,  0,
                            2,  0,  5, 18,  0,  0,  2,  0,  8,  6,  0,  0,
                            3,  0,  8,  8,  7,  0,  3,  0,  8,  6,  6,  0,
                            3,  1, 18,  5,  6,  0,  3,  1, 16,  3,  6,  0,
                            3,  0,  6,  6,  6,  0,  3,  0,  7,  8,  6,  0,
                            3,  1, 18,  5,  7,  0,  3,  1, 17,  4,  8,  0,
                            3,  1, 16,  3,  7,  0,  3,  1, 15,  2,  8,  0,
                            2,  0,  7,  8,  0,  0,  2,  0,  6,  6,  0,  0,
                            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                            0,  0,  0,  0,  0,  0,  1,  0, 11,  0,  0,  0,
                            1,  0, 12,  0,  0,  0,  1,  0, 11,  0,  0,  0,
                            1,  0, 12,  0,  0,  0,  2,  0,  1,  1,  0,  0,
                            3,  0,  6,  6,  6,  0,  3,  0,  7,  8,  6,  0,
                            3,  0,  1,  7,  8,  0,  3,  0,  1,  3,  2,  0,
                            0,  0,  0,  0,  0,  0,  3,  0,  7,  8, 23,  0,
                            3,  0,  6,  6, 23,  0,  2,  0,  1, 27,  0,  0,
                            2,  0,  1, 32,  0,  0,  2,  0,  1,  1,  0,  0,
                            3,  0,  6,  6,  6,  0,  2,  0,  7,  6,  0,  0,
                            2,  0,  8,  6,  0,  0,  2,  0,  7,  8,  0,  0,
                            2,  0, 21,  7,  0,  0,  2,  0,  9,  6,  0,  0,
                            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                            0,  0,  0,  0,  0,  0,  2,  0, 22,  8,  0,  0,
                            2,  0, 10,  6,  0,  0,  2,  0,  9,  8,  0,  0,
                            2,  0, 21,  6,  0,  0,  2,  0, 10,  7,  0,  0,
                            2,  0, 22,  6,  0,  0,  3,  0,  7,  8,  6,  0,
                            2,  0,  1,  6,  0,  0,  2,  0,  7,  8,  0,  0,
                            2,  0,  9, 10,  0,  0,  2,  0, 11, 12,  0,  0,
                            3,  0,  7,  8,  6,  0,  2,  0,  1, 23,  0,  0,
                            2,  0, 13,  6,  0,  0,  2,  0, 14,  7,  0,  0,
                            2,  0, 39,  1,  0,  0,  2,  0, 14,  8,  0,  0,
                            2,  0, 39,  6,  0,  0,  2,  0, 39,  8,  0,  0,
                            2,  0, 13,  8,  0,  0,  2,  0, 14,  6,  0,  0,
                            2,  0, 13,  7,  0,  0,  2,  0, 13,  6,  0,  0,
                            2,  0, 14,  7,  0,  0,  2,  0, 13,  8,  0,  0,
                            2,  0, 14,  6,  0,  0,  2,  0, 14,  8,  0,  0,
                            2,  0, 39,  7,  0,  0,  2,  0, 34,  6,  0,  0,
                            2,  0, 35,  7,  0,  0,  2,  0, 39,  6,  0,  0,
                            2,  0, 34,  8,  0,  0,  2,  0, 36,  7,  0,  0,
                            2,  0, 39,  8,  0,  0,  2,  0, 35,  8,  0,  0,
                            2,  0, 36,  6,  0,  0,  2,  0, 37,  6,  0,  0,
                            2,  0, 38,  7,  0,  0,  2,  0, 37,  8,  0,  0,
                            2,  0, 38,  6,  0,  0,  2,  0, 39, 10,  0,  0,
                            2,  0, 37,  8,  0,  0,  2,  0, 38,  6,  0,  0 };

    static double CBR[102] = { 1.,     1.,     1.,     0.,     1.,     1., 0.6351, 0.8468, 0.9027, 0.9200,
                              0.9518,     1., 0.6351, 0.8468, 0.9027, 0.9200, 0.9518,     1., 0.2160, 0.3398,
                              0.4748, 0.6098, 0.8049,     1., 0.6861,     1.,     0.,     0.,     0.,    0.5,
                                  1.,    0.5,     1., 0.3890, 0.7080, 0.9440, 0.9930,     1.,     0., 0.4420,
                              0.6470, 0.9470, 0.9770, 0.9990,     1.,     1.,     1.,     1., 0.6670,     1.,
                                  0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0., 0.6670,
                                  1., 0.6670,     1., 0.6670,     1., 0.8880, 0.9730,     1., 0.4950, 0.8390,
                              0.9870,     1., 0.5160,     1.,     1.,     1.,     1.,     1., 0.6410,     1.,
                                  1.,   0.67,     1.,   0.33,     1.,     1.,   0.88,   0.94,     1.,   0.88,
                                0.94,     1.,   0.88,   0.94,     1.,   0.33,     1.,   0.67,     1.,  0.678,
                               0.914,    1. };

    static int IDB[49] = { 0,  0,  0,  1,  2,  3,  5,  6,   7, 13,
                          19, 25,  0,  0,  0,  0,  0,  0,   0,  0,
                          30, 32, 34, 40, 46, 47, 48, 49,  60, 62,
                          64, 66, 69, 73, 75, 76, 77, 78,  79, 81,
                          82, 84, 86, 87, 90, 93, 96, 98, 100 };

    struct DECPAR_nonZero_output {
        int ND;
        int LL[10];
        double P_out[5][10];
    };

    DECPAR_nonZero_output DECPAR_nonZero(particle part) {
        // ***********************************************************************
        // This subroutine generates the decay of a particle
        // with ID = LA, and 5-momentum P0(1:5)
        // into ND particles of 5-momenta P(j,1:5) (j=1:ND)
        // (taken from SIBYLL 1.7, muon decay corrected, R.E. 04/98)
        // ***********************************************************************
        int LL[10] = { 0 };
        double P_out[5][10] = { 0. };
        double PV[5][10] = { 0. };
        double RORD[10] = { 0. };
        double UE[3] = { 0. };
        double BE[3] = { 0. };
        const double FACN[8] = { 2., 5., 15., 60., 250., 1500., 12000., 120000. };  // was in FORTRAN counted from position 3 to 10.
        const double PI = 3.1415927;

        // c.m.s. Momentum in two particle decays
        std::function<double(double, double, double)> PAWT = [&](double A, double B, double C) {return std::sqrt((A * A - (B + C) * (B + C)) * (A * A - (B - C) * (B - C))) / (2. * A); };

        int LA = part.ID;
        double P0[5] = { part.px, part.py, part.pz, part.E, part.m0 };

        // Choose decay channel
        int L = std::abs(LA);
        int ND = 0;
        int IDC = IDB[L - 1] - 1;

        if (IDC + 1 <= 0) {
            DECPAR_nonZero_output dno;
            dno.ND = ND;
            for (int i = 0; i < ND; ++i) {
                dno.LL[i] = LL[i];
                for (int j = 0; j < 5; ++j) {
                    dno.P_out[j][i] = P_out[j][i];
                }
            }
            return dno;
        }
        Random &random = Random::instance();
        double RBR = random.rand();
        do {
            IDC++;
        } while (RBR > CBR[IDC - 1]);

        int KD = 6 * (IDC - 1) + 1;
        ND = KDEC[KD - 1];
        double MAT = KDEC[KD];

        int MBST = (MAT > 0 && P0[3] > 20. * P0[4]) ? 1 : 0;
        double PS = 0.;
        for (int J = 0; J < ND; ++J) {
            LL[J] = KDEC[KD + 1 + J];
            P_out[4][J] = AM[LL[J] - 1];
            PV[4][J] = AM[LL[J] - 1];
            PS += P_out[4][J];
        }
        for (int J = 0; J < 4; ++J) {
            PV[J][0] = 0.;
            if (MBST == 0) PV[J][0] = P0[J];
        }
        if (MBST == 1) PV[3][0] = P0[4];
        PV[4][0] = P0[4];

        double WT = 0.;
        double WWTMAX = 0.;
        bool goto280 = (ND == 2);
        if (goto280 == false) {
            if (ND == 1) {
                // return with one particle
                for (int J = 0; J < 4; ++J) {
                    P_out[J][0] = P0[J];
                }
                DECPAR_nonZero_output dno;
                dno.ND = ND;
                for (int i = 0; i < ND; ++i) {
                    dno.LL[i] = LL[i];
                    for (int j = 0; j < 5; ++j) {
                        dno.P_out[j][i] = P_out[j][i];
                    }
                }
                return dno;
            }

            // Calculate maximum weight for ND-particle decay
            WWTMAX = 1. / FACN[ND - 3];
            double PMAX = PV[4][0] - PS + P_out[4][ND - 1];
            double PMIN = 0.;
            for (int IL = ND - 1; IL > 0; --IL) {
                PMAX += P_out[4][IL - 1];
                PMIN += P_out[4][IL];
                WWTMAX *= PAWT(PMAX, PMIN, P_out[4][IL - 1]);
            }
        }

        // 280 continuing into loop
        // generation of the masses, compute weight, if rejected try again
        bool repeat240 = false;
        do {  // 240
            repeat240 = false;

            if (goto280 == false) {
                RORD[0] = 1.;
                for (int IL1 = 2; IL1 < ND; ++IL1) {
                    double RSAV = random.rand();
                    int IL2tmp = 0;
                    for (int IL2 = IL1 - 1; IL2 > 0; --IL2) {
                        IL2tmp = IL2;
                        if (RSAV <= RORD[IL2 - 1]) break;
                        RORD[IL2] = RORD[IL2 - 1];
                    }
                    RORD[IL2tmp] = RSAV;
                }

                RORD[ND - 1] = 0.;
                WT = 1.;
                for (int IL = ND - 1; IL > 0; --IL) {
                    PV[4][IL - 1] = PV[4][IL] + P_out[4][IL - 1] + (RORD[IL - 1] - RORD[IL]) * (PV[4][0] - PS);
                    WT *= PAWT(PV[4][IL - 1], PV[4][IL], P_out[4][IL - 1]);
                }
                if (WT < random.rand() * WWTMAX) { repeat240 = true; continue; }
            }
            goto280 = false;  // 280

            // Perform two particle decays in respective cm frame
            for (int IL = 1; IL < ND; ++IL) {
                double PA = PAWT(PV[4][IL - 1], PV[4][IL], P_out[4][IL - 1]);
                UE[2] = 2. *random.rand() - 1.;
                double PHI = 2. * PI * random.rand();
                double UT = std::sqrt(1. - UE[2] * UE[2]);
                UE[0] = UT * std::cos(PHI);
                UE[1] = UT * std::sin(PHI);
                for (int J = 0; J < 3; ++J) {
                    P_out[J][IL - 1] = PA * UE[J];
                    PV[J][IL] = -PA * UE[J];
                }
                P_out[3][IL - 1] = std::sqrt(PA * PA + P_out[4][IL - 1] * P_out[4][IL - 1]);
                PV[3][IL] = std::sqrt(PA * PA + PV[4][IL] * PV[4][IL]);
            }

            // Lorentz transform decay products to lab frame
            for (int J = 0; J < 4; ++J) {
                P_out[J][ND - 1] = PV[J][ND - 1];
            }

            for (int IL = ND - 1; IL > 0; --IL) {  // 340 (first)
                for (int J = 0; J < 3; ++J) {
                    BE[J] = PV[J][IL - 1] / PV[3][IL - 1];
                }
                double GA = PV[3][IL - 1] / PV[4][IL - 1];
                for (int I = IL; I < ND + 1; ++I) {  // 340 (second)
                    double BEP = BE[0] * P_out[0][I - 1] + BE[1] * P_out[1][I - 1] + BE[2] * P_out[2][I - 1];
                    for (int J = 0; J < 3; ++J) {
                        P_out[J][I - 1] += GA * (GA * BEP / (1. + GA) + P_out[3][I - 1]) * BE[J];
                    }
                    P_out[3][I - 1] = GA * (P_out[3][I - 1] + BEP);
                }
            }

            // Weak decays
            if (MAT == 1) {
                double F1 = P_out[3][1] * P_out[3][2] - P_out[0][1] * P_out[0][2] - P_out[1][1] * P_out[1][2] - P_out[2][1] * P_out[2][2];
                if (MBST == 1) WT = P0[4] * P_out[3][0] * F1;
                if (MBST == 0) WT = F1 * (P_out[3][0] * P0[3] - P_out[0][0] * P0[0] - P_out[1][0] * P0[1] - P_out[2][0] * P0[2]);
                double WTMAX = std::pow(P0[4], 4) / 16.;
                if (WT < random.rand() * WTMAX) { repeat240 = true; continue; }
            }
        } while (repeat240);

        // Boost back for rapidly moving particle
        if (MBST == 1) {
            for (int J = 0; J < 3; ++J) {
                BE[J] = P0[J] / P0[3];
            }
            double GA = P0[3] / P0[4];
            for (int I = 0; I < ND; ++I) {
                double BEP = BE[0] * P_out[0][I] + BE[1] * P_out[1][I] + BE[2] * P_out[2][I];
                for (int J = 0; J < 3; ++J) {
                    P_out[J][I] += GA * (GA * BEP / (1. + GA) + P_out[3][I]) * BE[J];
                }
                P_out[3][I] = GA * (P_out[3][I] + BEP);
            }
        }

        // labels for antiparticle decay
        if (LA < 0 && L > 18) {
            for (int J = 0; J < ND; ++J) {
                LL[J] = LBARP[LL[J] - 1];
            }
        }

        DECPAR_nonZero_output dno;
        dno.ND = ND;
        for (int i = 0; i < ND; ++i) {
            dno.LL[i] = LL[i];
            for (int j = 0; j < 5; ++j) {
                dno.P_out[j][i] = P_out[j][i];
            }
        }
        return dno;
    }


	void PionDecay::process(Candidate* candidate) const {

        int CandidateId = candidate.current.getID();
        if (fabs(CandidateID) != 211) {
            return; // decay only for pions
        }
        
        const double E = candidate.current.getEnergy(); // check Units!
        Vector3d P = candidate.current.getMomentum(); // check Units!
        for (int i = 0; i < 10; ++i) {

            std::vector<particle> particles;

            const int ID;
            if (CandidateID == 211) {   // pi+
                ID = 7;
            }
            else {  // pi-
                ID = 8;
            }

            //const double E = 2.e7;

            particle primary;
            primary.ID = ID;
            primary.E = E;
            primary.m0 = AM[ID - 1];
            primary.px = P.x;
            primary.py = P.y;
            primary.pz = P.z;
            particles.push_back(primary);

            int nParticlesIN = 0;
            int nParticlesOUT = 0;
            do {
                nParticlesIN = particles.size();
                for (int i = 0; i < nParticlesIN; ++i) {
                    DECPAR_nonZero_output dno = DECPAR_nonZero(particles[i]);
                    if (dno.ND > 1) {
                        for (int j = 0; j < dno.ND; ++j) {
                            particle secondary;
                            secondary.ID = dno.LL[j];
                            secondary.px = dno.P_out[0][j];
                            secondary.py = dno.P_out[1][j];
                            secondary.pz = dno.P_out[2][j];
                            secondary.E = dno.P_out[3][j];
                            secondary.m0 = dno.P_out[4][j];
                            particles.push_back(secondary);
                        }
                        particles.erase(particles.begin() + i);
                    }
                }
                nParticlesOUT = particles.size();
            } while (nParticlesIN != nParticlesOUT);

            // create Secondaries
            for (int i = 1; i < particles.size(); ++i) { // start with 1: do not create primary as secondary
                int SecID; 
                switch (particles[i].ID)
                {
                case 3: // electron
                    SecID = 11;
                    break;
                case 15: // electron-neutrino
                    SecID = 12;
                    break;
                case 16: // anti-elcectron-neutrino
                    SecID = -12;
                    break;
                case 17:  // muon-neutrino
                    SecID = 14;
                    break;
                case 18: // anit-muon-neutrino
                    SecID = -14;
                    break;
                }
                Vector3d P;
                P.x = particles[i].px;
                P.y = particles[i].py;
                P.z = particles[i].pz;
                Vector3d Pos = candidate.current.getPosition();
                double E = particles[i].E;

                candidate.addSecondary(new Candidate(SecID, E, Pos, P));

            }
        }  // repetitions

	}

} // namespace