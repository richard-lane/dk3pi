#include "cleo_ll.h"

void test_cleo_ll()
{
    double parameters[64];
    parameters[0]  = 20418.9;    // D0{K-,pi+,pi+,pi-}_D0{K+,pi-,pi-,pi+}_N
    parameters[1]  = 17920.6;    // D0{K-,pi+,pi+,pi-}_D0{K+,pi-}_N
    parameters[2]  = 4189;       // D0{K-,pi+}_D0{K+,pi-}_N
    parameters[3]  = 54259;      // D0{K-,pi+,pi0}_D0{K+,pi-,pi0}_N
    parameters[4]  = 63703.4;    // D0{K-,pi+,pi+,pi-}_D0{K+,pi-,pi0}_N
    parameters[5]  = 31044.6;    // D0{K-,pi+,pi0}_D0{K+,pi-}_N
    parameters[6]  = 0.0553431;  // r_k3pi
    parameters[7]  = 0.0440767;  // r_kpipi0
    parameters[8]  = 0.627289;   // R_k3pi_0
    parameters[9]  = 99.5899;    // d_k3pi_0
    parameters[10] = 0.817555;   // R_kpipi0
    parameters[11] = 202.984;    // d_kpipi0
    parameters[12] = 0.0590658;  // r_kpi
    parameters[13] = 1;          // R_kpi
    parameters[14] = 190.356;    // d_kpi
    parameters[15] = 0.999999;   // R_k3pi_1
    parameters[16] = 132.766;    // d_k3pi_1
    parameters[17] = 0.483193;   // R_k3pi_2
    parameters[18] = 165.398;    // d_k3pi_2
    parameters[19] = 0.260427;   // R_k3pi_3
    parameters[20] = 331.317;    // d_k3pi_3
    parameters[21] = 0.00143751; // D0{pi+,pi-}_BR
    parameters[22] = 0.00681448; // y
    parameters[23] = 0.0790972;  // D0{K-,pi+,pi+,pi-}_BR
    parameters[24] = 0.00406694; // D0{K+,K-}_BR
    parameters[25] = 3249.77;    // D0{K0S0,pi0}_N
    parameters[26] = 0.0119;     // D0{K0S0,pi0}_BR
    parameters[27] = 3467.6;     // D0{K0S0,omega(782)0{pi+,pi-,pi0}}_N
    parameters[28] = 0.0099;     // D0{K0S0,omega(782)0{pi+,pi-,pi0}}_BR
    parameters[29] = 3732.31;    // D0{K0S0,pi0,pi0}_N
    parameters[30] = 0.0091;     // D0{K0S0,pi0,pi0}_BR
    parameters[31] = 4203.02;    // D0{K0S0,phi(1020)0{K+,K-}}_N
    parameters[32] = 0.002;      // D0{K0S0,phi(1020)0{K+,K-}}_BR
    parameters[33] = 3523.35;    // D0{K0S0,eta0}_N
    parameters[34] = 0.00189;    // D0{K0S0,eta0}_BR
    parameters[35] = 2671.25;    // D0{K0S0,eta0{pi+,pi-,pi0}}_N
    parameters[36] = 0.0011;     // D0{K0S0,eta0{pi+,pi-,pi0}}_BR
    parameters[37] = 3082.67;    // D0{K0S0,eta'(958)0{eta0,pi+,pi-}}_N
    parameters[38] = 0.00159;    // D0{K0S0,eta'(958)0{eta0,pi+,pi-}}_BR
    parameters[39] = 3464.94;    // D0{K0L0,pi0}_N
    parameters[40] = 0.01;       // D0{K0L0,pi0}_BR
    parameters[41] = 3469.86;    // D0{K0L0,omega(782)0{pi+,pi-,pi0}}_N
    parameters[42] = 0.0099;     // D0{K0L0,omega(782)0{pi+,pi-,pi0}}_BR
    parameters[43] = 3500.22;    // D0{pi+,pi-,pi0}_N
    parameters[44] = 0.0147;     // D0{pi+,pi-,pi0}_BR
    parameters[45] = 0.973;      // f_pipipi0
    parameters[46] = 2889.28;    // D0{K-,pi+,pi0}_D0{K0S0,pi+,pi-}_N
    parameters[47] = 0.655;      // D0{K0S0,pi+,pi-}::c1
    parameters[48] = -0.025;     // D0{K0S0,pi+,pi-}::s1
    parameters[49] = 2174.13;    // D0{K0S0,pi+,pi-}_N
    parameters[50] = 0.511;      // D0{K0S0,pi+,pi-}::c2
    parameters[51] = 0.141;      // D0{K0S0,pi+,pi-}::s2
    parameters[52] = 0.024;      // D0{K0S0,pi+,pi-}::c3
    parameters[53] = 1.111;      // D0{K0S0,pi+,pi-}::s3
    parameters[54] = -0.569;     // D0{K0S0,pi+,pi-}::c4
    parameters[55] = 0.328;      // D0{K0S0,pi+,pi-}::s4
    parameters[56] = -0.903;     // D0{K0S0,pi+,pi-}::c5
    parameters[57] = -0.181;     // D0{K0S0,pi+,pi-}::s5
    parameters[58] = -0.616;     // D0{K0S0,pi+,pi-}::c6
    parameters[59] = -0.52;      // D0{K0S0,pi+,pi-}::s6
    parameters[60] = 0.1;        // D0{K0S0,pi+,pi-}::c7
    parameters[61] = -1.129;     // D0{K0S0,pi+,pi-}::s7
    parameters[62] = 0.422;      // D0{K0S0,pi+,pi-}::c8
    parameters[63] = -0.35;      // D0{K0S0,pi+,pi-}::s8

    std::cout << cleo_ll(parameters) << std::endl;
}
