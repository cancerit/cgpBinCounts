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
use warnings;

use Getopt::Long 'GetOptions';
use Pod::Usage;

use Bio::DB::HTS;
use Carp;

use Const::Fast qw(const);

use Sanger::CGP::BinAllele;

const my $DEFAULT_READS_PER_BIN => 200; # this gave the best representation when comparing against old human version.


{
  my $options = option_builder();
  process($options);
}

sub process {
  my ($options) = @_;
  open my $OUT, '>', $options->{'o'} or croak 'Failed to create '.$options->{'o'}."\n$!\n";

  my $sam = Bio::DB::HTS->new(-fasta => $options->{'g'},
                              -bam   => $options->{'i'},
                              -split_splices => 1
                              );

  my %feat_opts = ( -flags=>{MAP_PAIR=>1},
                    -iterator =>1,
                    -type => 'match');

  my $proper_pairs = $sam->features( %feat_opts );
  my ($last_chr, $bin_start);
  my $read_starts = 0;
  my $bin_counts = 0;
  while (my $align = $proper_pairs->next_seq) {
    $last_chr = $align->seq_id unless(defined $last_chr);
    $bin_start = $align->start unless(defined $bin_start);

    if($last_chr ne $align->seq_id) {
      $read_starts = 0;
      $last_chr = $align->seq_id;
      $bin_start = $align->start;
      $bin_counts = 0;
    }

    next if($bin_start > $align->start);

    $read_starts++;
    if($read_starts == $options->{'b'}) {
      my @entry = (  'bin'.++$bin_counts
                    , $align->seq_id
                    , $bin_start
                    , $align->start
                    , gc_frac($sam, $align->seq_id, $bin_start, $align->start));
      if($options->{'f'} eq 'csv') {
        print $OUT join ',',@entry;
      }
      else {
        print $OUT join "\t", ($entry[1],$entry[2]-1,$entry[3],$entry[1].':'.$entry[0],$entry[4]);
      }
      print $OUT "\n";
      $bin_start = $align->start + 1;
      $read_starts = 0;
#exit if($bin_counts > 100);
    }
  }
  close $OUT;
}

sub gc_frac {
  my ($sam, $seq_id, $start, $end) = @_;
  my $dna_str = $sam->fai->fetch($seq_id.':'.$start.'-'.$end);
  my $dna_len = length $dna_str;
  $dna_str =~ tr/atAT//d;
  return sprintf '%.4f', (length $dna_str) / $dna_len;
}

sub option_builder {
	my ($factory) = @_;

	my %opts = ();

	my $result = &GetOptions (
		'h|help' => \$opts{'h'},
		'd|debug' => \$opts{'d'},
		'i|input=s' => \$opts{'i'},
		'o|output=s' => \$opts{'o'},
		'g|genome=s' => \$opts{'g'},
		'f|format=s' => \$opts{'f'},
		'b|binsize=n' => \$opts{'b'},
		'r|restrict=s' => \$opts{'r'},
		'v|version' => \$opts{'v'},
	);
	pod2usage(1) unless ($result);
	pod2usage(1) if($opts{'h'});
	if($opts{'v'}) {
	  print 'VERSION: ', Sanger::CGP::BinAllele->VERSION,"\n";
	  exit;
	}
	validate_input(\%opts);
	return \%opts;
}

sub validate_input {
  my $options = shift;
  pod2usage(qq{\n'-i' must be defined.\n}) unless(defined $options->{'i'});
  pod2usage(qq{\n'-o' must be defined.\n}) unless(defined $options->{'o'});
  pod2usage(qq{\n'-g' must be defined.\n}) unless(defined $options->{'g'});

  pod2usage(qq{\n}.$options->{'i'}.qq{ is not a valid file.\n}) unless(-f $options->{'i'});
  pod2usage(qq{\n}.$options->{'g'}.qq{ is not a valid file.\n}) unless(-f $options->{'g'});

  if($options->{'f'}) {
    pod2usage(qq{\n}.$options->{'f'}.qq{ must be 'csv' or 'bed'.\n}) if($options->{'f'} ne 'csv' && $options->{'f'} ne 'bed');
  }
  else {
    $options->{'f'} = 'bed';
  }

  my $output = $options->{'o'};
  $output .= '.'.$options->{'r'} if(defined $options->{'r'});
  $output .= '.'.$options->{'f'};
  if(-f $output) {
    pod2usage(qq{\n}.$output.qq{ exists but cannot be written to.\n}) unless(-w $output);
  }
  $options->{'o'} = $output;

  $options->{'b'} = $DEFAULT_READS_PER_BIN unless($options->{'b'});

  return 1;
}


__END__

=head1 NAME

gcBinGeneration.pl - Generate GC bins

=head1 SYNOPSIS

gcBinGeneration.pl [-h] -i <BAM> -o <FILE> -f [csv|bed] -g <genome.fa> [-r <CHR>]

  General Options:

    --help      (-h)  Brief documentation

    --input     (-i)  InSilico Bam file

    --output    (-o)  Output (format extension will be appended)

    --format    (-f)  csv or bed output
                        - csv gives format used by NGS copynumber
                        - bed

    --genome    (-g)  Reference genome

    --binsize   (-b)  Bin size, defaults to 200 reads (based on 30x coverage)

    --restrict  (-r)  Restrict to this chromosome (*.bam.bai file required)
                        - restriction will be added between output and format.

    --version   (-v)  Print version and exit.

  Examples:

    gcBinGeneration.pl -i inSilico.bam -f tsv -o genome.gc

=cut
