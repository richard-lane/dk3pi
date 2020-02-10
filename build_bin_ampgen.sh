#!/bin/bash
set -e

OUTNAME=$1
if [[ -z $OUTNAME ]] ; then
    OUTNAME=foo.exe
fi

# Build helper lib
if [ ! -f ./libD2K3pi.a ]; then
    ./lib/build_lib.sh
fi

# Build with debug symbols without optimisation
g++ src/bin_generated_decays.cpp libD2K3pi.a -o $OUTNAME `root-config --cflags --glibs` \
-lboost_filesystem -Wall -Wextra -Wformat-security -Werror -O0 -g
