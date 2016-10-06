#!/bin/bash

########## LICENCE ##########
# Copyright (c) 2014-2016 Genome Research Ltd.
#
# Author: Cancer Genome Project <cgpit@sanger.ac.uk>
#
# This file is part of cgpBinCounts.
#
# cgpBinCounts is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’."
########## LICENCE ##########


# Bio::DB::HTS
SOURCE_BIOBDHTS="https://github.com/Ensembl/Bio-HTS/archive/2.3.tar.gz"
SOURCE_HTSLIB="https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2"

get_file () {
# output, source
  if hash curl 2>/dev/null; then
    curl --insecure -sS -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

get_distro () {
  EXT=""
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
  elif [[ $2 == *.zip* ]] ; then
    EXT="zip"
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi
  rm -f $1.$EXT
  if hash curl 2>/dev/null; then
    curl --retry 10 -sS -o $1.$EXT -L $2
  else
    wget --tries=10 -nv -O $1.$EXT $2
  fi
}

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/pancan"
  exit 0
fi

INST_PATH=$1

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

# get current directory
INIT_DIR=`pwd`

# re-initialise log file
echo > $INIT_DIR/setup.log

# log information about this system
(
    echo '============== System information ===='
    set -x
    lsb_release -a
    uname -a
    sw_vers
    system_profiler
    grep MemTotal /proc/meminfo
    set +x
    echo; echo
) >>$INIT_DIR/setup.log 2>&1

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

# make sure that build is self contained
unset PERL5LIB
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"

#add bin path for install tests
export PATH="$INST_PATH/bin:$PATH"

cd $INIT_DIR

## grab cpanm:
rm -f $SETUP_DIR/cpanm
get_file $SETUP_DIR/cpanm http://xrl.us/cpanm
chmod +x $SETUP_DIR/cpanm

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi

if [ -e $SETUP_DIR/basePerlDeps.success ]; then
  echo "Previously installed base perl deps..."
else
#"ExtUtils::CBuilder" "Module::Build~0.42" "File::ShareDir" "File::ShareDir::Install" "Const::Fast" "File::Which" "LWP::UserAgent"
  perlmods=("Bio::Root::Version~1.006009001")
  for i in "${perlmods[@]}" ; do
    $SETUP_DIR/cpanm -v --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH $i
  done
  touch $SETUP_DIR/basePerlDeps.success
fi

CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Bio::DB::HTS`
if [[ "x$CHK" == "x" ]] ; then
  echo -n "Get htslib ..."
  if [ -e $SETUP_DIR/htslibGet.success ]; then
    echo " already staged ...";
  else
    echo
    cd $SETUP_DIR
    get_distro "htslib" $SOURCE_HTSLIB
    echo " DONE";
    touch $SETUP_DIR/htslibGet.success
  fi
  echo -n "Building Bio::DB::HTS..."
  cd $SETUP_DIR
  rm -rf bioDbHts
  get_distro "bioDbHts" $SOURCE_BIOBDHTS

  mkdir -p bioDbHts/htslib
  tar --strip-components 1 -C bioDbHts -zxf bioDbHts.tar.gz
  tar --strip-components 1 -C bioDbHts/htslib -jxf $SETUP_DIR/htslib.tar.bz2
  cd bioDbHts/htslib
  perl -pi -e 'if($_ =~ m/^CFLAGS/ && $_ !~ m/\-fPIC/i){chomp; s/#.+//; $_ .= " -fPIC -Wno-unused -Wno-unused-result\n"};' Makefile
  make -j$CPU
  rm -f libhts.so*
  cd ../
  env HTSLIB_DIR=$SETUP_DIR/bioDbHts/htslib perl Build.PL --install_base=$INST_PATH
  ./Build test
  ./Build install
  cd $SETUP_DIR
  rm -f bioDbHts.tar.gz
  echo " DONE"
else
  echo "Bio::DB::HTS already installed"
fi

exit

set -x
$SETUP_DIR/cpanm -v --mirror http://cpan.metacpan.org --notest -l $INST_PATH/ --installdeps . < /dev/null
set +x

set -e
cd $INIT_DIR
perl Makefile.PL INSTALL_BASE=$INST_PATH
make
make test
make install

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo

exit 0
