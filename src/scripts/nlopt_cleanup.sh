#!/bin/sh

# Create directory if it doesnt exist
mkdir -p ../inst/include

cp nlopt${R_ARCH}/include/* ../inst/include/
rm -fr nlopt-src
