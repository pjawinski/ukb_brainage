# ===========================
# === prepare finemapping ===
# ===========================

# prepare CAVIAR
cd /fast/software
git clone https://github.com/fhormoz/caviar.git

# in order to compile: add /include path of conda environment with gsl to Makefile
cd CAVIAR-C++
addLDFLAGS="-I /slow/projects/ukb_brainage/envs/default/include"
awk -v addLDFLAGS="${addLDFLAGS}" '$1=="LDFLAGS=" { print $0, addLDFLAGS; next } { print }' Makefile > Makefile.tmp; \mv Makefile.tmp Makefile
make
ln -s /fast/software/caviar/CAVIAR-C++/CAVIAR /fast/software/bin/caviar

# prepare FINEMAP
cd /fast/software
mkdir finemap; cd finemap
wget http://christianbenner.com/finemap_v1.4.1_x86_64.tgz
wget http://christianbenner.com/finemap_v1.4_x86_64.tgz
tar -xvzf finemap_v1.4.1_x86_64.tgz
tar -xvzf finemap_v1.4_x86_64.tgz
ln -s /fast/software/finemap/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 /fast/software/bin/finemap
ln -s /fast/software/finemap/finemap_v1.4_x86_64/finemap_v1.4_x86_64 /fast/software/bin/finemap14

# prepare LDstore2
cd /fast/software
mkdir ldstore; cd ldstore
wget http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz
tar -xvzf ldstore_v2.0_x86_64.tgz
ln -s /fast/software/ldstore/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 /fast/software/bin/ldstore
