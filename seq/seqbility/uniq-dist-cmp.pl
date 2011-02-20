#!/usr/bin/perl -w

use strict;

die("Usage: ud-cmp.pl <in1.ud> <in2.ud>\n") if (@ARGV != 2);
# output: chr pos which score half-same?
my ($fh0, $fh1);
open($fh0, $ARGV[0]) || die;
open($fh1, $ARGV[1]) || die;

my (@x, $y, %z);
$y = 3;
while (1) {
  if ($y&1) {
	$_ = <$fh0>;
	last unless ($_);
	@{$x[0]} = split;
	$z{k50}[0] += $x[0][3] if ($x[0][3] >= 50000);
	$z{k80}[0] += $x[0][3] if ($x[0][3] >= 80000);
  }
  if ($y&2) {
	$_ = <$fh1>;
	last unless ($_);
	@{$x[1]} = split;
	$z{k50}[1] += $x[1][3] if ($x[1][3] >= 50000);
	$z{k80}[1] += $x[1][3] if ($x[1][3] >= 80000);
  }
  if ($x[0][0] lt $x[1][0] || ($x[0][0] eq $x[1][0] && $x[0][2] < $x[1][2])) {
	if ($x[0][0] eq $x[1][0]) { # first is wrong
	  print join("\t", $x[0][0], $x[0][2], 1, $x[1][3], ($x[0][1]==$x[1][1])? 1:0), "\n";
	  ++$z{k50}[2] if ($x[1][3] >= 50000);
	  ++$z{k80}[2] if ($x[1][3] >= 80000);
	}
	$y = 1;
  } elsif ($x[0][0] eq $x[1][0] && $x[0][2] == $x[1][2]) {
	$y = 3;
  } else {
	if ($x[0][0] eq $x[1][0]) { # second is wrong
	  print join("\t", $x[1][0], $x[1][2], 2, $x[0][3], ($x[0][1]==$x[1][1])? 1:0), "\n";
	  ++$z{k50}[3] if ($x[0][3] >= 50000);
	  ++$z{k80}[3] if ($x[0][3] >= 80000);
	}
	$y = 2;
  }
}
close($fh0); close($fh1);
print STDERR join("\t", "k50", $z{k50}[1]/$z{k50}[2], $z{k50}[0]/$z{k50}[3]), "\n";
print STDERR join("\t", "k80", $z{k80}[1]/$z{k80}[2], $z{k80}[0]/$z{k80}[3]), "\n";
