#!/bin/bash

set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
MODEL_DIR=$( dirname $DIR)/amplitude_models/

g++ -fPIC -shared $DIR/cf_wrapper.cpp -I $MODEL_DIR -o $DIR/cf_wrapper.so
g++ -fPIC -shared $DIR/dcs_wrapper.cpp -I $MODEL_DIR -o $DIR/dcs_wrapper.so
