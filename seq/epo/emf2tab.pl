#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = (0=>'homo_sapiens', 1=>'pan_troglodytes');
getopts('0:1:', \%opts);
die("Usage: emf2fa.pl [-0 $opts{0}] [-1 $opts{1}] <in.emf>\n") if (@ARGV == 0 && -t STDIN);

my @species = ($opts{0}, $opts{1});
my (@hdr, %dup, $flag, $id, @col, $seq, $strand);
$id = $flag = 0; $strand = 1;
@col = (-1, -1);
while (<>) {
  if (/^SEQ/) {
	my $h = \@{$hdr[$id]};
	@$h = split;
	shift(@$h);
	++$dup{$h->[0]};
	if ($species[0] eq $h->[0]) {
	  $strand = $h->[4];
	  $col[0] = $id;
	} elsif ($species[1] eq $h->[0]) {
	  $col[1] = $id;
	}
	++$id;
  } elsif (/^DATA/) {
	$flag = 1 if ($col[0] >= 0 && $col[1] >= 0 && $dup{$species[0]} == 1 && $dup{$species[1]} == 1);
	$seq = '';
  } elsif (/^\/\//) {
	my $h = \@{$hdr[$col[0]]};
	if ($strand < 0) {
	  $seq = reverse($seq);
	  $seq =~ tr/acgtrymkACGTRYMK/tgcayrkmTGCAYRKM/;
	}
	if ($seq && length($seq) != $h->[3] - $h->[2] + 1) { # bad region
	  warn("$h->[1]:$h->[2]-$h->[3]\n");
	} else {
	  print join("\t", @$h[1..3], $seq), "\n" if ($flag);
	}
	%dup = (); @hdr = (); @col = (-1, -1);
	$flag = $id = 0;
  } elsif ($flag) {
	my @c = (substr($_, $col[0], 1), substr($_, $col[1], 1));
	if ($c[0] ne '-') {
	  $seq .= ($c[1] eq '-')? 'n' : $c[1];
	}
  }
}
