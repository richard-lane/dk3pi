#include "BESIII_chi2.h"

void test_BESIII()
{

    Double_t parameters[78];
    parameters[0]  = 5.81976e-01;  // Rk3pi0
    parameters[1]  = 7.80911e-01;  // Rk3pi1
    parameters[2]  = 8.49155e-01;  // Rk3pi2
    parameters[3]  = 4.53797e-01;  // Rk3pi3
    parameters[4]  = 1.31418e+02;  // deltak3pi0
    parameters[5]  = 1.49528e+02;  // deltak3pi1
    parameters[6]  = 1.76478e+02;  // deltak3pi2
    parameters[7]  = 2.73767e+02;  // deltak3pi3
    parameters[8]  = 4.41481e-02;  // rkpipi0
    parameters[9]  = 7.95940e-01;  // Rkpipi0
    parameters[10] = 2.00422e+02;  // deltakpipi0
    parameters[11] = 9.12111e+03;  // Norm0kspipipi
    parameters[12] = 7.58751e+03;  // Norm1kspipipi
    parameters[13] = 7.86921e+03;  // Norm2kspipipi
    parameters[14] = 1.03944e+04;  // Norm3kspipipi
    parameters[15] = 1.97604e+02;  // deltakpi
    parameters[16] = 3.93451e-02;  // Bkpi_CF
    parameters[17] = 1.47674e-01;  // Bkpipi0_CF
    parameters[18] = 8.32820e-02;  // Bk3pi_CF
    parameters[19] = 2.14488e-03;  // Bkpipi0_DCS/Bkpipi0_CF
    parameters[20] = 3.22223e-03;  // Bk3pi_DCS/Bk3pi_CF
    parameters[21] = 3.44355e-03;  // rkpi^2
    parameters[22] = 5.43595e-02;  // rk3pi0
    parameters[23] = 5.81052e-02;  // rk3pi1
    parameters[24] = 5.75272e-02;  // rk3pi2
    parameters[25] = 5.09254e-02;  // rk3pi3
    parameters[26] = 3.65404e-03;  // x_mixing
    parameters[27] = 6.87437e-03;  // y_mixing
    parameters[28] = 6.99649e-01;  // ci_1
    parameters[29] = 6.45961e-01;  // ci_2
    parameters[30] = -2.58903e-03; // ci_3
    parameters[31] = -6.07998e-01; // ci_4
    parameters[32] = -9.54417e-01; // ci_5
    parameters[33] = -5.89476e-01; // ci_6
    parameters[34] = 5.91171e-02;  // ci_7
    parameters[35] = 4.11390e-01;  // ci_8
    parameters[36] = 8.54049e-02;  // si_1
    parameters[37] = 3.10419e-01;  // si_2
    parameters[38] = 1.00831e+00;  // si_3
    parameters[39] = 6.48371e-01;  // si_4
    parameters[40] = -3.37068e-02; // si_5
    parameters[41] = -5.68738e-01; // si_6
    parameters[42] = -8.32666e-01; // si_7
    parameters[43] = -4.31009e-01; // si_8
    parameters[44] = 1.73402e-01;  // ki_1
    parameters[45] = 8.76005e-02;  // ki_2
    parameters[46] = 6.91972e-02;  // ki_3
    parameters[47] = 2.55000e-02;  // ki_4
    parameters[48] = 8.50011e-02;  // ki_5
    parameters[49] = 5.91996e-02;  // ki_6
    parameters[50] = 1.26900e-01;  // ki_7
    parameters[51] = 1.33800e-01;  // ki_8
    parameters[52] = 7.93998e-02;  // ki_-1
    parameters[53] = 1.74000e-02;  // ki_-2
    parameters[54] = 2.02000e-02;  // ki_-3
    parameters[55] = 1.62001e-02;  // ki_-4
    parameters[56] = 5.12002e-02;  // ki_-5
    parameters[57] = 1.42999e-02;  // ki_-6
    parameters[58] = 1.32000e-02;  // ki_-7
    parameters[59] = 2.74999e-02;  // ki_-8
    parameters[60] = 9.72667e-01;  // Fpipipi0
    parameters[61] = 6.12508e+04;  // Norm_kspipivskpipi0
    parameters[62] = 1.75216e-01;  // kip_1
    parameters[63] = 8.81138e-02;  // kip_2
    parameters[64] = 6.97188e-02;  // kip_3
    parameters[65] = 2.54630e-02;  // kip_4
    parameters[66] = 8.51398e-02;  // kip_5
    parameters[67] = 5.86190e-02;  // kip_6
    parameters[68] = 1.26872e-01;  // kip_7
    parameters[69] = 1.33835e-01;  // kip_8
    parameters[70] = 8.01789e-02;  // kip_-1
    parameters[71] = 1.74489e-02;  // kip_-2
    parameters[72] = 2.01895e-02;  // kip_-3
    parameters[73] = 1.67866e-02;  // kip_-4
    parameters[74] = 5.19655e-02;  // kip_-5
    parameters[75] = 1.32375e-02;  // kip_-6
    parameters[76] = 1.32823e-02;  // kip_-7
    parameters[77] = 2.70838e-02;  // kip_-8

    std::cout << BESIII_chi2(parameters) << std::endl;
}
