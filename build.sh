#!/bin/bash

cmake -DCMAKE_BUILD_TYPE=DEBUG -S . -B build && \
cd build && \
make && \
cd ..
