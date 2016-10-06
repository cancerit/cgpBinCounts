#!/usr/bin/perl

use strict;
use warnings;

use IO::Uncompress::Gunzip qw($GunzipError);

my ($allele_or_bin, $bins_in, $bg_out) = @ARGV;

open my $BG, '>', $bg_out or die "$!: $bg_out\n";

my $gunzip = new IO::Uncompress::Gunzip $bins_in or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
while(my $line = <$gunzip>) {
  chomp $line;
  my @F = (split /,/, $line)[0..3];
  if($allele_or_bin eq 'A') {
    $F[3] += $F[2];
    $F[2] = $F[1];
  }
  $F[1]--;
  print $BG join("\t",@F),"\n";
}
close $gunzip;
close $BG;
