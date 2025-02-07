#!/bin/csh

echo "=== Downloading Source ==="
set SYSNAME = `uname -s`
set PROCTYPE = `uname -m`

echo $SYSNAME
if ( -f cspice.tar.Z || -f cspice.tar ) then
    wget "Found existing cspice file, skipping download."
else if ( $SYSNAME == "Linux" && $PROCTYPE == "x86_64" ) then
    wget https://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Linux_GCC_64bit/packages/cspice.tar.Z
else if ( $SYSNAME == "Linux" && $PROCTYPE == "i686" ) then
    wget https://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Linux_GCC_32bit/packages/cspice.tar.Z
else if ( $SYSNAME == "Darwin" && $PROCTYPE == "x86_64" ) then
    wget https://naif.jpl.nasa.gov/pub/naif/toolkit//C/MacIntel_OSX_AppleC_64bit/packages/cspice.tar.Z
else if ( $SYSNAME == "Darwin" && $PROCTYPE == "arm64" ) then
    wget https://naif.jpl.nasa.gov/pub/naif/toolkit//C/MacM1_OSX_clang_64bit/packages/cspice.tar.Z
else
    echo "Failed to identify machine to set which version of CSPICE to download.  Please get it yourself from"
    echo "https://naif.jpl.nasa.gov/naif/toolkit_C.html and place the file cspice.tar.Z in this directory,"
    echo "then run this script again.  Sorry about that."
    exit
endif

echo "=== Extracting Files ==="
if ( -f cspice.tar.Z) then
    gzip -d cspice.tar.Z
endif
tar xfv cspice.tar

echo "=== Compiling CSPICE ==="
cd cspice
csh makeall.csh
# Reformatting things like an actual shared library
cp ./lib/cspice.a ./lib/libcspice.a
mkdir ./include/cspice
mv ./include/*.h ./include/cspice/
