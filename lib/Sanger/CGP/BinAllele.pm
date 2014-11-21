package Sanger::CGP::BinAllele;

use strict;

use Const::Fast qw(const);
use base 'Exporter';
use Capture::Tiny qw(capture);

our $VERSION = '1.0.0';
our @EXPORT = qw($VERSION);

const my $LICENSE =>
'#################
# Copyright (c) 2014 Genome Research Ltd.
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
#################';

sub license {
  return sprintf $LICENSE, $VERSION;
}

sub sample_name {
  my $bam = shift;
  my @lines = split /\n/, Bio::DB::Sam->new(-bam => $bam)->header->text;
  my $sample;
  for(@lines) {
    if($_ =~ m/^\@RG.*\tSM:([^\t]+)/) {
      $sample = $1;
      last;
    }
  }
  die "Failed to determine sample from BAM header\n" unless(defined $sample);
  return $sample;
}

sub md5file {
  my $file = shift;
  my ($stdout, $stderr, $exit) = capture { system(qq{md5sum $file | awk '{print \$1}' > $file.md5}); };
  die $stderr if($exit != 0);
}

1;
