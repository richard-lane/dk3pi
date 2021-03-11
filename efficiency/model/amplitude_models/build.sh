#!/bin/bash

set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

g++ -fPIC -O3 -shared $DIR/cf_wrapper.cpp -I $DIR -o $DIR/cf_wrapper.so
g++ -fPIC -O3 -shared $DIR/dcs_wrapper.cpp -I $DIR -o $DIR/dcs_wrapper.so
