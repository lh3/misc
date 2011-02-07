#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

&main;
exit;

sub main {
  if (@ARGV == 0) {
	print STDERR qq(
Usage:    bedutils.pl <command> [arguments]

Commands: union   collapse overlaps in adjacent regions (sorted)
          cov2    compute regions covered by at least two regions (sorted)

);
	exit(1);
  }
  my $command = shift(@ARGV);
  my %func = (union=>\&f_union, cov2=>\&f_cov2);
  die("Unknown command \"$command\".\n") if (!defined($func{$command}));
  &{$func{$command}};
}

sub f_union {
  die("Usage: bedutils.pl union <in.sorted.bed>\n") if (@ARGV == 0 && -t STDIN);
  my ($c, $b, $e);
  $c = '';
  while (<>) {
	my @t = split;
	$t[1] = 0 if ($t[1] < 0);
	if ($t[0] ne $c || $t[1] > $e) {
	  print join("\t", $c, $b, $e), "\n" if ($c);
	  ($c, $b, $e) = @t[0..2];
	} else {
	  $e = $t[2] if ($e < $t[2]);
	}
  }
  print join("\t", $c, $b, $e), "\n";
}

sub f_cov2 {
  die("Usage: bedutils.pl cov2 <in.sorted.bed>\n") if (@ARGV == 0 && -t STDIN);
  my ($c, $e, $ob, $oe);
  # $e: rightmost pos; $ob/$oe: begin/end of the overlap to be ouputted
  # $oe <= $e; $ob < $oe
  $c = ''; $ob = $oe = -1;
  while (<>) {
	my @t = split;
	$t[1] = 0 if ($t[1] < 0);
	if ($t[0] ne $c || $t[1] >= $e) { # no overlap
	  print join("\t", $c, $ob, $oe), "\n" if ($ob != -1);
	  ($c, $e) = @t[0,2];
	  $ob = $oe = -1;
	} else { # overlap
	  if ($ob == -1 || $t[1] >= $oe) { # $oe <= $t[1] < $e here if $ob == -1
		print join("\t", $c, $ob, $oe), "\n" if ($ob != -1);
		$ob = $t[1]; $oe = $t[2] < $e? $t[2] : $e;
	  } else {
		my $new_oe = $t[2] < $e? $t[2] : $e;
		$oe = $new_oe if ($oe < $new_oe);
	  }
	  $e = $t[2] if ($e < $t[2]);
	}
  }
  print join("\t", $c, $ob, $oe), "\n" if ($ob != -1);
}
