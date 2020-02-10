#!/bin/bash
set -e

OUTNAME=$1
if [[ -z $OUTNAME ]] ; then
    OUTNAME=pull.exe
fi

# Build helper lib
if [ ! -f ./libD2K3pi.a ]; then
    ./lib/build_lib.sh
fi

# Build without debug symbols, including optimisation
g++ src/pull_study.cpp src/PullStudyHelpers.cpp libD2K3pi.a -o $OUTNAME `root-config --cflags --glibs` \
-lboost_filesystem -Wall -Wextra -Wformat-security -Werror -O3 -g
