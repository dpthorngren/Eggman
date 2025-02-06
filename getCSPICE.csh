#!/bin/csh

echo "=== Downloading Source ==="
wget https://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Linux_GCC_64bit/packages/cspice.tar.Z
echo "=== Extracting Files ==="
gzip -d cspice.tar.Z
tar xfv cspice.tar
echo "=== Compiling CSPICE ==="
cd cspice
csh makeall.csh
