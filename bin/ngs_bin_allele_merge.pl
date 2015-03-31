#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014,2015 Genome Research Ltd.
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

use FindBin;
use lib "$FindBin::Bin/../lib";

use strict;
use warnings FATAL => qw(all);
use autodie qw(:all);

use File::Path qw(make_path);
use File::Copy;
use Capture::Tiny qw(capture);

use Bio::DB::Sam;

use Sanger::CGP::BinAllele;

if(scalar @ARGV < 2) {
  warn 'VERSION: '.Sanger::CGP::BinAllele->VERSION."\n";
  die "USAGE: ./ngs_bin_allele_merge.pl control.bam tum_control_dir\n";
}

my $control_bam = shift;
my $indir = shift;

my $cntl_sample = Sanger::CGP::BinAllele::sample_name($control_bam);
my %tumour_samples;
my %file_sets;
opendir(my $dh, $indir);
while(my $item = readdir $dh) {
  next if($item =~ m/^[.]/);
  if($item =~ m/^(.*)\.(allelic_counts|bin_counts).*\.csv\.gz$/) {
    my $sample = $1;
    my $type = $2;
    push @{$file_sets{$type}{$sample}}, "$indir/$item";
    $tumour_samples{$sample} = 1 unless($sample eq $cntl_sample);
  }
}
closedir $dh;

my @types = qw(allelic_counts bin_counts);
for my $tum_samp(keys %tumour_samples) {
  my $path = "$indir/$tum_samp";
  make_path $path unless(-d $path);
  for my $type(@types) {
    merge_files($file_sets{$type}{$tum_samp}, "$path/$tum_samp.$type.csv.gz");
    merge_files($file_sets{$type}{$cntl_sample}, "$path/$cntl_sample.$type.csv.gz");
  }
  my $archive = "$indir/$tum_samp.binnedReadCounts.tar.gz";
  my $command = join ' ', 'tar -C', $indir, '-czf', $archive, $tum_samp;
  my ($stdout, $stderr, $exit) = capture { system($command); };
  die $stderr if($exit != 0);
  Sanger::CGP::BinAllele::md5file($archive);
}

sub merge_files {
  my ($files, $dest) = @_;
  if(scalar @{$files} == 1) {
    copy($files->[0], $dest);
  }
  else {
    my $command = join ' ', 'zcat', @{$files}, q{| sort -t ',' -k1,1 -k2,2n | gzip -c >}, $dest;
    my ($stdout, $stderr, $exit) = capture { system($command); };
    die $stderr if($exit != 0);
  }
}
