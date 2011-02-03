#!/usr/bin/env perl

# this script is translated from the Lua version

use strict;
use warnings;
use Getopt::Std;

my %opts;
getopts('s', \%opts);
my $site_only = defined($opts{s})? 1 : 0;

if (@ARGV == 0) {
	print("Usage: perl vcf2afs.pl [-s] <samples.txt> <in.vcf>\n");
	exit(1);
}

my (%pop, %sample, $fp);
my $fn = shift(@ARGV);
open($fp, $fn) || die;
while (<$fp>) {
	if (/^(\S+)\s+(\S+)/) {
		$sample{$1} = $2;
		push(@{$pop{$2}}, $1);
	}
}
close($fp);

my (%col, %cnt);
for my $k (keys %pop) {
	@{$col{$k}} = ();
	@{$cnt{$k}} = (0);
}
while (<>) {
	if (/^##/) {
	} elsif (/^#/) {
		my @t = split("\t");
		for my $i (9 .. $#t) {
			my $k = $sample{$t[$i]};
			if ($k) {
				push(@{$col{$k}}, $i);
				push(@{$cnt{$k}}, 0, 0);
			}
		}
	} else {
		my @t = split("\t");
		if ($t[4] ne '.' && $t[4] !~ /,/) {
			if ($site_only) {
				print join("\t", @t[0,1,3,4]);
			}
			for my $k (keys %col) {
				my ($ac, $an) = (0, 0);
				my $v = $col{$k};
				for my $i (0 .. @$v-1) {
					if ($t[$$v[$i]] =~ /^(\d).(\d)/) {
						$ac += $1 + $2;
						$an += 2;
					}
				}
				print "\t$k:$an:$ac" if ($site_only);
				++$cnt{$k}[$ac] if ($an == @{$cnt{$k}} - 1);
			}
			print "\n" if ($site_only);
		}
	}	
}

if (!$site_only) {
	for my $k (keys %cnt) {
		my $v = $cnt{$k};
		print "$k\t", (@$v-1), "\t", join("\t", @$v), "\n";
	}
}
