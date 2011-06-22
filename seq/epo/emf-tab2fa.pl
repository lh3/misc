#!/usr/bin/perl -w

use strict;
use warnings;

die("Usage: cat *.tab | sort -k1,1N -k2,2n | emf-tab2fa.pl ref.fa.fai\n") if (@ARGV == 0);

my %lenaux;
my $fnfai = shift(@ARGV);
my $fp;
open($fp, $fnfai) || die;
while (<$fp>) {
	my @t = split;
	$lenaux{$t[0]} = $t[1];
}
close($fp);

my $len = \%lenaux;
my $last = '';
my $seq;
while (<>) {
  my @t = split;
  next if (!defined($len->{$t[0]}) || $t[2] > $len->{$t[0]});
  if ($t[0] ne $last) {
	if ($last) {
	  print ">$last\n";
	  for (my $i = 0; $i < length($seq); $i += 60) {
		print substr($seq, $i, 60), "\n";
	  }
	}
	$last = $t[0]; $seq = 'N'x$len->{$t[0]};
  }
  substr($seq, $t[1]-1, $t[2]-$t[1]+1) = $t[3];
}
if ($last) {
  print ">$last\n";
  for (my $i = 0; $i < length($seq); $i += 60) {
	print substr($seq, $i, 60), "\n";
  }
}
