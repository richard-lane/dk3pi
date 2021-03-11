#!/bin/bash

set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

g++ -fPIC -Ofast -rdynamic --std=c++14 -march=native -shared $DIR/cf_wrapper.cpp -I $DIR -o $DIR/cf_wrapper.so
g++ -fPIC -Ofast -rdynamic --std=c++14 -march=native -shared $DIR/dcs_wrapper.cpp -I $DIR -o $DIR/dcs_wrapper.so
