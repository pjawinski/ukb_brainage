#!/bin/bash

# =====================
# === prepare METAL ===
# =====================

# download file
cd /fast/software/
git clone https://github.com/statgen/METAL
cd METAL

# compile
cat CMakeLists.txt | awk 'NR==14 { print "include_directories(${ZLIB_INCLUDE_DIR})"; next} { print }' > tmp && \mv tmp CMakeLists.txt
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make test
make install

# create symbolic link in binary folder
cd /fast/software/bin
ln -s /fast/software/METAL/build/bin/metal metal