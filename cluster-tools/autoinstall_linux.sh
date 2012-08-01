#!/usr/bin/bash

# This script should install each of the packages in the current directory
# Could be buggy --TJL

# -- Set the installation directories --
WD=`pwd`                               
FFTW_PREFIX="$WD/fftw"
GMX_PREFIX="$WD/gromacs453"
CCTOOLS_PREFIX="$WD/cctools"
#---------------------------------------

# create appropriate directories
mkdir -p $FFTW_PREFIX
mkdir -p $GMX_PREFIX
mkdir -p $CCTOOLS_PREFIX

# (1) FFTW
echo Installing FFTW3...

wget http://www.fftw.org/fftw-3.3.1.tar.gz
tar -xzvf fftw-3.3.1.tar.gz
rm fftw-3.3.1.tar.gz
cd fftw-3.3.1

./configure --enable-float --enable-sse --enable-shared --enable-threads --prefix=$FFTW_PREFIX
make
make install
cd $WD

# (2) GROMACS
echo Installing GROMACS v4.5.3 ...

wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-4.5.3.tar.gz
tar -xzvf gromacs-4.5.3.tar.gz
rm gromacs-4.5.3.tar.gz
cd gromacs-4.5.3

export CPPFLAGS="-I$FFTW_PREFIX/include"
export LDFLAGS="-L$FFTW_PREFIX/lib"

# use this one to have mpi
#./configure --enable-mpi --enable-shared --enable-threads --without-x --with-fft=fftw3 --prefix=$GMX_PREFIX

# use this one for no mpi
./configure --enable-shared --enable-threads --without-x --with-fft=fftw3 --prefix=$GMX_PREFIX

make -j 12
make install
cd $WD

# (3) CCTools
echo Installing CCTools ...

wget http://www.cse.nd.edu/~ccl/software/files/cctools-3.4.2-src.tar.gz
tar -xzvf cctools-3.4.2-src.tar.gz
rm cctools-3.4.2-src.tar.gz
cd cctools-3.4.2-src

./configure --prefix $CCTOOLS_PREFIX

# IF NOT INSTALLING TRY THIS:
# the makefile wants us to use a flag that causes ld to trip up
# but we can just remove that flag (-lz) from the makefile with sed
#sed 's/-lz//' Makefile.config > Makefile.config.2
#mv Makefile.config.2 Makefile.config

make
make install
cd $WD

# (4) Echo out things to add to your bash profile

echo ""
echo "Add the following to your .bash_profile:"
echo ""
echo "export GMXBIN=$GMX_PREFIX/bin"
echo "export GMXLIB=$GMX_PREFIX/share/gromacs/top"
echo "export PATH=$CCTOOLS_PREFIX/bin:\$PATH"
echo "export PATH=\$GMXBIN:\$PATH"
echo "export PATH=`pwd`:\$PATH"
echo ""
echo "Done!"

