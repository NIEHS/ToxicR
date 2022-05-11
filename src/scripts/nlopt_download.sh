#! /bin/sh

RSCRIPT_BIN=$1
echo $1
echo $1
echo $1
echo $1
echo $1
echo "-----------------------------"
echo "-----------------------------"
echo "-----------------------------"
echo "-----------------------------"

# Uncompress NLOPT source
${RSCRIPT_BIN} -e "utils::untar(tarfile = 'nlopt-src.tar.gz')"
mv nlopt-2.7.1 nlopt-src
