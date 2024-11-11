#!/bin/sh

# Define the CMake version and source tarball
CMAKE_VERSION=3.31.0
CMAKE_TARBALL=cmake-${CMAKE_VERSION}.tar.gz
CMAKE_SOURCE_DIR=cmake-${CMAKE_VERSION}

# Extract the tarball
tar -xzf src/${CMAKE_TARBALL} -C src

# Change to the source directory
cd src/${CMAKE_SOURCE_DIR}

# Configure and install CMake
./bootstrap --prefix=$PWD/../cmake_install
make
make install

# Set the CMAKE_BIN variable to the installed CMake binary
CMAKE_BIN=$PWD/../cmake_install/bin/cmake

# Export the CMAKE_BIN variable
export CMAKE_BIN

# Return to the original directory
cd ../../
