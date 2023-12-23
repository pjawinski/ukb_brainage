# ==================================
# === prepare enigma surfaceplot ===
# ==================================

# set working directory
cd /Users/philippe/software/
cd /fast/software

# create environment
mkdir -p ENIGMA/envs/
conda create -p ENIGMA/envs/surfplotPy -c anaconda -c conda-forge python=3.7 matplotlib scikit-learn numpy scipy vtk=8.1.2 xvfbwrapper nibabel nilearn pillow pandas=1.2.5 openssl=1.0.2 mscorefonts jupyter xvfbwrapper xorg-x11-server-xvfb-cos6-x86_64 x11-xkb-utils binutils 
conda activate ENIGMA/envs/surfplotPy

# get enigma toolbox
cd ENIGMA
git clone https://github.com/MICA-MNI/ENIGMA.git
mv ENIGMA/* .

# run setup
python setup.py install

# export python environment
conda env export -p envs/surfplotPy > envs/surfplotPy.yml
awk 'NR==1 { print "name: surfplotPy"; next } $1=="prefix:" { $2="surfplotPy"; print; next} { print }' envs/surfplotPy.yml > tmp && \mv tmp envs/surfplotPy.yml
ln -s $CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/bin/xvfb-run $CONDA_PREFIX/bin/xvfb-run
ln -s $CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/bin/Xvfb $CONDA_PREFIX/bin/Xvfb
ln -s $CONDA_PREFIX/aarch64-conda-cos7-linux-gnu/sysroot/usr/bin/xkbcomp /usr/bin/xkbcomp

ln -s $CONDA_PREFIX/lib/libcrypto.so.1.0.0 $CONDA_PREFIX/lib/libcrypto.so.10

# get libxfont
mkdir libxfont
cd libxfont
wget http://ftp.debian.org/debian/pool/main/libx/libxfont/libxfont1_1.5.1-1+deb8u1_amd64.deb
ar x libxfont1_1.5.1-1+deb8u1_amd64.deb
tar -xf data.tar.xz
ln -s /fast/software/ENIGMA/libxfont/usr/lib/x86_64-linux-gnu/libXfont.so.1.4.1 $CONDA_PREFIX/lib/libXfont.so.1

# get xkbcomp
wget http://ftp.de.debian.org/debian/pool/main/x/x11-xkb-utils/x11-xkb-utils_7.7+7_i386.deb
ar x x11-xkb-utils_7.7+7_i386.deb
tar -xf data.tar.gz
ln -s $(pwd)/usr/bin/xkbcomp $CONDA_PREFIX/bin/xkbcomp
scp /Users/philippe/Downloads/default.xkm jawinskp@cluster3.psychologie.hu-berlin.de:/fast/software/ENIGMA/usr/bin/xkbcomp

# get keyboard config
wget https://www.x.org/archive/individual/data/xkeyboard-config/xkeyboard-config-2.32.tar.bz2
tar -xf xkeyboard-config-2.32.tar.bz2
cd xkeyboard-config-2.32
aclocal --verbose -I $CONDA_PREFIX/share/aclocal
./configure $CONDA_PREFIX/bin

mkdir $HOME/tmp
Xvfb -xkbdir $CONDA_PREFIX/bin/
ls -l

ldd $CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/bin/Xvfbq


# get xorg macros
git clone https://github.com/freedesktop/xorg-macros
cd xorg-macros
./configure -prefix=$CONDA_PREFIX
make
make install
cd ..

pwd
# install xkbcomp
git clone https://gitlab.freedesktop.org/xorg/app/xkbcomp
cd xkbcomp
chmod 770 configure.ac 
aclocal --verbose -I $CONDA_PREFIX/share/aclocal
autoconf configure.ac[]
automake
nano configure.aclocal
./configure.ac --prefix=$CONDA_PREFIX


wget https://distrib-coffee.ipsl.jussieu.fr/pub/linux/altlinux/p10/branch/x86_64/RPMS.classic/libXfont-1.5.4-alt2.x86_64.rpm

ldd $CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/bin/Xvfb


git clone https://github.com/freedesktop/xorg-libXfont


ln -s $(pwd)/envs/surfplotPy/aarch64-conda-cos7-linux-gnu/sysroot/usr/lib64/libXfont2.so.2.0.0 $(pwd)/envs/surfplotPy/x86_64-conda-linux-gnu/sysroot/usr/bin/../../../../lib/libXfont.so.1

ln -s $(pwd)/envs/surfplotPy/aarch64-conda-cos7-linux-gnu/sysroot/usr/lib64/libXfont2.so.2.0.0 $(pwd)/envs/surfplotPy/lib/libXfont.so.1
ln -s $(pwd)/envs/surfplotPy/aarch64-conda-cos7-linux-gnu/sysroot/usr/lib64/libXfont2.so.2.0.0 $(pwd)/envs/surfplotPy/bin/libXfont.so.1
ln -s $(pwd)/envs/surfplotPy/aarch64-conda-cos7-linux-gnu/sysroot/usr/lib64/libXfont2.so.2.0.0 $(pwd)/envs/surfplotPy/x86_64-conda-linux-gnu/sysroot/usr/lib64/libXfont.so.1


${pwd}/envs/surfplotPy/x86_64-conda-linux-gnu/sysroot/usr/bin/../../../../lib/libcrypto.so.10
lib/libXfont.so.1

$(pwd)/envs/surfplotPy/lib/libcrypto.so.10
ln -s $(pwd)/envs/surfplotPy/x86_64-conda-linux-gnu/sysroot/usr/bin/Xvfb envs/surfplotPy/bin/Xvfb



wget https://www.x.org/archive/individual/lib/libXfont-1.5.4.tar.bz2 -P $(pwd)/envs/surfplotPy/lib/; tar -xf $(pwd)/envs/surfplotPy/lib/libXfont-1.5.4.tar.bz2



sudo ln -s /usr/lib/x86_64-linux-gnu/libssl.so.1.0.0 /usr/lib/libssl.so.10

sudo ln -s /usr/lib/x86_64-linux-gnu/libcrypto.so.1.0.0 /usr/lib/libcrypto.so.10


conda deactivate

# create R environment
conda create -p envs/surfplotR -c r r-base=4.1.3 r-ggplot2 r-ggpubr r-magick r-patchwork 
conda activate envs/surfplotR
conda env export -p envs/surfplotR > envs/surfplotR.yml
awk 'NR==1 { print "name: surfplotR"; next } $1=="prefix:" { $2="surfplotR"; print; next} { print }' envs/surfplotR.yml > tmp && \mv tmp envs/surfplotR.yml
conda deactivate

# craete symbolic links in project folder and copy .yml files
cd /slow/projects/ukb_brainage
ln -s /fast/software/ENIGMA/envs/surfplotPy envs/surfplotPy
ln -s /fast/software/ENIGMA/envs/surfplotR envs/surfplotR
cp /fast/software/ENIGMA/envs/surfplotPy.yml envs/surfplotPy.yml
cp /fast/software/ENIGMA/envs/surfplotR.yml envs/surfplotR.yml
