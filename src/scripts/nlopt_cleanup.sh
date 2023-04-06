#!/bin/sh

mkdir ../inst/include
cp nlopt${R_ARCH}/include/* ../inst/include/
rm -fr nlopt-src
