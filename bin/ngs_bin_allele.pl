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
use IO::Uncompress::Gunzip qw($GunzipError) ;
use IO::Compress::Gzip qw($GzipError) ;
use Const::Fast qw(const);

use Bio::DB::HTS;

const my $MIN_MAPQ => 35;

use Sanger::CGP::BinAllele;

if(scalar @ARGV < 4) {
  warn 'VERSION: '.Sanger::CGP::BinAllele->VERSION."\n";
  die "USAGE: ./ngs_bin_allele.pl cn_bins.csv[.gz] snp6.csv[.gz] outdir my.bam [this_chr_only]\n";
}

my $bin_file = shift;
my $allele_file = shift;
my $outdir = shift;
my $bam_in = shift;
my $restrict_to;
$restrict_to = shift if(scalar @ARGV);

make_path $outdir unless(-d $outdir);

bin_counts($outdir, $bam_in, $bin_file, $restrict_to);
allele_counts($outdir, $bam_in, $allele_file, $restrict_to);

sub allele_counts {
  my ($out_l, $bam_l, $allele_l, $restrict_l) = @_;
  my $sam = Bio::DB::HTS->new(-bam => $bam_l);
  my $header = $sam->header;

  my $sample = Sanger::CGP::BinAllele::sample_name($bam_l);
  my $outfile = "$out_l/$sample.allelic_counts";
  $outfile .= ".$restrict_l" if(defined $restrict_l);
  $outfile .= '.csv.gz';

  my $gzip = new IO::Compress::Gzip $outfile or die "gzip to $outfile failed: $GzipError\n";
  my $gunzip = new IO::Uncompress::Gunzip $allele_l or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
  while(my $line = <$gunzip>) {
    chomp $line;
    my ($chr, $start, undef, undef, $allA, $allB, $freq) = split /[,]/, $line;
    next if(defined $restrict_l && $chr ne $restrict_l);

    my($countA, $countB) = (0,0);
    $sam->fast_pileup("$chr:$start-$start",
                      sub {
                        my ($seqid, $pos, $pu) = @_;
                        return if($pos != $start);
                        foreach my $p (@{$pu}) {
                          next if($p->indel || $p->is_refskip);
                          my $a = $p->alignment;
                          next unless($a->flag & 2); # proper pair
                          next if($a->flag & 256); # secondary alignment
                          next if($a->flag & 512); # vendor fail
                          next if($a->flag & 1024); # duplicate
                          next if($a->flag & 2048); # supplimentary
                          next if($a->flag & 4); # unmapped
                          next if($a->flag & 8); # mate unmapped
                          next if($a->qual < $MIN_MAPQ);


                          # get the base at this pos
                          my $qbase  = substr($a->qseq, $p->qpos, 1);
                          if($qbase eq $allA) {
                            $countA++;
                          }
                          elsif($qbase eq $allB) {
                            $countB++;
                          }
                        }
                        1;
                      }
    );
    print $gzip join(q{,}, $chr, $start, $countA, $countB, $freq),"\n" or die "Failed to write line to $outfile: $!\n";
  }
  close $gzip;
  close $gunzip;
}

sub bin_counts {
  my ($out_l, $bam_l, $bins_l, $restrict_l) = @_;
  my $sam = Bio::DB::HTS->new(-bam => $bam_l);
  my $bam = $sam->hts_file;
  my $header = $sam->header;
  my $bai = Bio::DB::HTSfile->index_load($bam);

  my $sample = Sanger::CGP::BinAllele::sample_name($bam_l);
  my $outfile = "$out_l/$sample.bin_counts";
  $outfile .= ".$restrict_l" if(defined $restrict_l);
  $outfile .= '.csv.gz';

  my $gzip = new IO::Compress::Gzip $outfile or die "gzip to $outfile failed: $GzipError\n";
  my $gunzip = new IO::Uncompress::Gunzip $bins_l or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
  my $no_gc_warned = 0;
  while(my $line = <$gunzip>) {
    chomp $line;
    my ($chr, $start, $end, $gc) = split /[,\t]/, $line;
    next if(defined $restrict_l && $chr ne $restrict_l);
    my $count = 0;
    $bai->fetch($bam,
                $header->parse_region("$chr:$start-$end"),
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
    if(defined $gc) {
      print $gzip join(q{,}, $chr, $start, $end, $count, $gc),"\n" or die "Failed to write line to $outfile: $!\n";
    }
    else {
      warn "WARN: No GC column in file $bins_l, ommitting from output\n" unless($no_gc_warned++);
      print $gzip join(q{,}, $chr, $start, $end, $count),"\n" or die "Failed to write line to $outfile: $!\n";
    }
  }
  close $gzip;
  close $gunzip;
}
