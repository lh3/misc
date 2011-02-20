#!/usr/bin/perl -w

use strict;
use warnings;

&main;
exit;

sub main {
  die("Usage: gen_mask.pl <bwa.sam>\n") if (@ARGV == 0 && -t STDIN);
  my @conv = (0 .. 127);
  # calculate conv table
  for (1 .. 127) {
	$conv[$_] = &int_log2($_) + 1;
  }
  $conv[0] = 0;

  # core loop
  my ($last, $seq) = ('', '');
  my $lineno = 0;
  while (<>) {
	if (/^@/) {
	  ++$lineno;
	  next;
	} elsif (/^(\S+)_(\d+)/) {
	  ++$lineno;
	  my ($a0, $a1) = (0, 0);
	  if ($last ne $1) {
		&print_seq($last, \$seq) if ($last);
		$last = $1; $seq = '';
	  }
	  $a0 = $1 > 127? 127 : $1 if (/X0:i:(\d+)/);
	  $a1 = $1 > 127? 127 : $1 if (/X1:i:(\d+)/);
	  $seq .= chr(63 + ($conv[$a0]<<3 | $conv[$a1]));
	}
	warn("$lineno lines processed\n") if ($lineno % 1000000 == 0);
  }
  &print_seq($last, \$seq);
}

sub print_seq {
  my ($last, $seq) = @_;
  print ">$last\n";
  for (my $i = 0; $i < length($$seq); $i += 60) {
	print substr($$seq, $i, 60), "\n";
  }
}

sub int_log2 {
  my $v = shift;
  my $c = 0;
  if ($v & 0xffff0000) { $v >>= 16; $c |= 16; }
  if ($v & 0xff00) { $v >>= 8; $c |= 8; }
  if ($v & 0xf0) { $v >>= 4; $c |= 4; }
  if ($v & 0xc) { $v >>= 2; $c |= 2; }
  if ($v & 0x2) { $c |= 1; }
  return $c;
}
