#include "../charmFitter/include/fitter/CleoCombinationFitter.h"

/*
 *g++ tmp.cpp -L charmFitter/libcharmFitter.a -L common/libcommon.a -I ../charmFitter/include/ -I ../common/include/ -I ../charmFitter/CLEO/ `root-config --glibs --cflags` -o tmp.exe -std=c++17
 */
int main(void) {
    const FitterUtil::DecayParams_t params {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::cout << "likelihood: "  << CLEO::cleoLikelihood(CLEO::binNumber(-30.0), params) <<  std::endl;

    std::cout << "im trying my best" << std::endl;
    return 0;
}

