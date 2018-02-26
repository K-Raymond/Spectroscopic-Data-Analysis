#!/bin/bash

# Download and build root

# ie. 6.10/00 -> 6100
VERSION="6.10.00"
URL="https://root.cern.ch/download/"
OS=$(lsb_release -si)

# Installing nessary packages
apt -v
if [ $? -eq 0 ]
then
    sudo apt-get install git dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev
fi
yum -v
if [ $? -eq 0 ]
then
    sudo yum install git cmake gcc-c++ gcc binutils libX11-devel libXpm-devel libXft-devel libXext-devel
fi

cd ../ # move to head directory

if [ ! -r Root ]
then
    mkdir Root
fi

cd Root

if [ -e root_v$VERSION.source.tar.gz ]
then
    echo "Using cached file"
else
    echo "Downloading ROOT"
    wget $URL"root_v"$VERSION".source.tar.gz"
fi

tar xf "root_v"$VERSION".source.tar.gz"
if [ ! $? -eq 0 ]
then
    exit
fi

if [ ! -e root_build ]
then
    mkdir root_build
fi
cd root_build

cmake ../root-$VERSION
cmake --build . -- -s -j $(nproc)

source ./bin/thisroot.sh

echo "Compiling compeleted sucessfully! Make sure to source ./bin/thisroot.sh"
