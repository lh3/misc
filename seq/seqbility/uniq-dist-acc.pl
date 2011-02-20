#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = (n=>0);
getopts('n:', \%opts);
die("Usage: uniq-dist-acc.pl <uniq-dist-output>\n") if (@ARGV == 0 && -t STDIN);
my @a;
my $sum = 0;
while (<>) {
  my @t = split;
  push(@a, $t[3]); $sum += $t[3];
}
@a = reverse(sort {$a<=>$b} @a);
my $z = 1.0;
if ($opts{n} > 0) {
  my $s = 0;
  for my $x (@a) {
	if ($s/$sum > 0.5) {
	  $z = $opts{n} / $x;
	  last;
	}
	$s += $x;
  }
}
my $s = 0;
for (@a) {
  printf("%d\t%.4g\n", int($_ * $z + .499), $s/$sum);
  $s += $_;
  last if ($s/$sum + 1e-6 > 1);
}
print "1\t1\n";
