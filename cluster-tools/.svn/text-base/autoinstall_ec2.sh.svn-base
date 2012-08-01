#!/bin/bash
# Compile and install gromacs, work_queue & dependencies on ec2 amazon linux AMI
# This linux flavor is a very minimal envirnoment
#
# Simply executing this script should be enough to download, compile, build and
# install everything. Fingers crossed.


# install build tools that don't come on the amu
echo "y" | sudo yum install gcc-c++
echo "y" | sudo yum install gcc
echo "y" | sudo yum install make
echo "y" | sudo yum install perl-devel

# install fftw (fourrier transform library required by gromacs)
wget http://www.fftw.org/fftw-3.3.1.tar.gz
tar -xzvf fftw-3.3.1.tar.gz
rm fftw-3.3.1.tar.gz
cd fftw-3.3.1
./configure --enable-shared --enable-float --enable-threads
make
sudo make install
cd $HOME

# install gromacs
wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-4.5.2.tar.gz
tar -xzvf gromacs-4.5.2.tar.gz
rm gromacs-4.5.2.tar.gz
cd gromacs-4.5.2
./configure --enable-shared --enable-threads --without-x --with-fft=fftw3
make -j 8 # TJL added to use 8 threads
sudo make install
# gromacs was complaining that it couldn't find libfftw3f.so.3, so you need
# to just put it in the gromacs lib folder
# TJL - try this:
export CPPFLAGS="-I/usr/local/include"
export LDFLAGS="-L/usr/local/lib"
#sudo cp /usr/local/lib/libfftw3f.so.3 /usr/local/gromacs/lib/
cd $HOME

# install two perl modules required by cctools (actually, the installation
# of the one required by cctools requires this one)
wget http://search.cpan.org/CPAN/authors/id/M/MS/MSCHWERN/ExtUtils-MakeMaker-6.62.tar.gz
tar -xvzf ExtUtils-MakeMaker-6.62.tar.gz
cd ExtUtils-MakeMaker-6.62
perl Makefile.PL
make
sudo make install
cd $HOME

# install a perl module required by cctools
wget http://www.cpan.org/modules/by-module/ExtUtils/ExtUtils-Embed-1.14.tar.gz
tar -xvzf ExtUtils-Embed-1.14.tar.gz
rm ExtUtils-Embed-1.14.tar.gz
cd ExtUtils-Embed-1.14
perl Makefile.PL
make
sudo make install
cd $HOME

# install cctools (work_queue) itself
wget http://www.cse.nd.edu/~ccl/software/files/cctools-3.4.2-src.tar.gz
tar -xzvf cctools-3.4.2-src.tar.gz
rm cctools-3.4.2-src.tar.gz
cd cctools-3.4.2-src
./configure
# the makefile wants us to use a flag that causes ld to trip up
# but we can just remove that flag (-lz) from the makefile with sed
sed 's/-lz//' Makefile.config > Makefile.config.2
mv Makefile.config.2 Makefile.config
make
sudo make install
cd $HOME


# cctools installs by default to $HOME/cctools, so we just add its bin to
# our path
echo "export PATH=\$PATH:\$HOME/cctools/bin" >> .bashrc

# gromacs by default is installed to /usr/local/gromacs, so we
# need to add its stuff to the bashrc
echo "source /usr/local/gromacs/bin/GMXRC.bash" >> .bashrc
source .bashrc
