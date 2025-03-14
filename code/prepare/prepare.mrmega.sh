#!/bin/bash

# =======================
# === install mr-mega ===
# =======================

# set working directory
mkdir /fast/software/mrmega
cd /fast/software/mrmega

# download mrmega
wget https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip
wget https://tools.gi.ut.ee/tools/fixP.r
wget https://tools.gi.ut.ee/tools/manh.r
wget https://tools.gi.ut.ee/tools/qq.r

# install
unzip MR-MEGA_v*.zip 
make
ln -s /fast/software/mrmega/MR-MEGA /fast/software/bin/MR-MEGA