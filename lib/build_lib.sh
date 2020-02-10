#!/bin/bash

for sourcefile in ./lib/*.cpp; do
    g++ -c $sourcefile -o lib/"$(basename $sourcefile .cpp)".o
done

ar rvs libD2K3pi.a ./lib/*.o
