g++ test/manual/test_simulator_equations.cpp src/DecaySimulator.cpp \
src/util.cpp src/MCGenerator.cpp -o manual_tests.exe `root-config --cflags --glibs` \
-lboost_filesystem
