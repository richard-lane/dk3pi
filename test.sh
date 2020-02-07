# Ensure the repo is clean
set -e
echo "Cleaning repo"
git clean -xfd

# Build pull study and binning executables
echo "Building executables"
./build_bin_ampgen.sh
./build_pull_study.sh

# Build and run "unit" tests
echo "Running UT"
./test/ut/build_tests.sh
./ut.exe

# Build and run slower tests
echo "Running integration tests"
echo "Running integration tests"
./test/it/build_integration_tests.sh
./it.exe

# Manual tests
echo "Running manual tests"
./test/manual/build_manual_tests.sh
./manual_tests.exe
