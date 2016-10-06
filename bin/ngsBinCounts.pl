#!/usr/bin/perl

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

use FindBin;
use lib "$FindBin::Bin/../lib";

use strict;
use warnings FATAL => qw(all);
use autodie qw(:all);

use File::Path qw(remove_tree make_path);
use Const::Fast qw(const);

use Bio::DB::HTS;

const my $MIN_MAPQ => 35;

use Sanger::CGP::BinAllele;

if(scalar @ARGV < 3) {
  warn 'VERSION: '.Sanger::CGP::BinAllele->VERSION."\n";
  die "USAGE: ./ngs_bin_allele.pl bins.bed outdir my.bam [this_chr_only]\n";
}

my $bin_file = shift;
my $outdir = shift;
my $bam_in = shift;
my $restrict_to;
$restrict_to = shift if(scalar @ARGV);

make_path $outdir unless(-d $outdir);

bin_counts($outdir, $bam_in, $bin_file, $restrict_to);

sub bin_counts {
  my ($out_l, $bam_l, $bins_l, $restrict_l) = @_;
  my $sam = Bio::DB::HTS->new(-bam => $bam_l);
  my $bam = $sam->hts_file;
  my $header = $sam->header;
  my $bai = Bio::DB::HTSfile->index_load($bam);

  my $sample = Sanger::CGP::BinAllele::sample_name($bam_l);
  my $outfile = "$out_l/$sample.bin_counts";
  $outfile .= ".$restrict_l" if(defined $restrict_l);
  $outfile .= '.bed';

  open my $OUT, '>', $outfile or die "$!: $outfile";
  open my $IN, '<', $bins_l or die "$!: $bins_l";
  while(my $line = <$IN>) {
    chomp $line;
    my ($chr, $start, $end) = split /\t/, $line;
print "$chr\r";
    next if(defined $restrict_l && $chr ne $restrict_l);
    my $count = 0;
    $bai->fetch($bam,
                $header->parse_region(sprintf "%s:%d-%d", $chr, $start+1, $end),
                sub {
                  my $a = shift;
                  return unless($a->flag & 2); # proper pair
                  return if($a->flag & 256); # secondary alignment
                  return if($a->flag & 512); # vendor fail
                  return if($a->flag & 1024); # duplicate
                  return if($a->flag & 2048); # supplimentary
                  return if($a->flag & 4); # unmapped
                  return if($a->flag & 8); # mate unmapped
                  return if($a->qual < $MIN_MAPQ);
                  return if($a->pos < $start || $a->pos > $end);
                  $count++;
                  1;
                }
    );
    print $OUT join(qq{\t}, $chr, $start, $end, $count),"\n" or die "Failed to write line to $outfile: $!\n";
  }
  close $IN;
  close $OUT;
}
