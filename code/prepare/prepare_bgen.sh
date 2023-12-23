#!/bin/bash

# ====================
# === Install BGEN ===
# ====================

# download BGEN
mkdir -p /fast/software/bgen
cd /fast/software/bgen
wget http://code.enkre.net/bgen/tarball/release/bgen.tgz
tar -xzf bgen.tgz; mv bgen.tgz bgen
cd bgen

# compile bgen
./waf configure
./waf

# test it
./build/test/unit/test_bgen
./build/apps/bgenix -g example/example.16bits.bgen -list

# make symbolic link for cat-bgen
ln -s /fast/software/bgen/build/apps/cat-bgen /fast/software/bin/cat-bgen
ln -s /fast/software/bgen/build/apps/bgenix /fast/software/bin/bgenix
