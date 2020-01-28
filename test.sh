# Ensure the repo is clean
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

# Manual tests
echo "Running manual tests"
./test/manual/build_manual_tests.sh
./manual_tests.exe
