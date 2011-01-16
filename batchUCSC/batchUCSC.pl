#!/usr/bin/perl -w

# Author: lh3
# Version: 0.2.0

use strict;
use warnings;
use DBI;
use Getopt::Std;

&main;
exit;

sub main {
  my %preset =
	(snpAnc=>['snp129OrthoPt2Pa2Rm2', '', '',
			  'chrom,chromStart,chromEnd-chromStart,humanAllele,chimpAllele,orangAllele,macaqueAllele'],
	 snp=>['snp129', '', '', '*'],
	 refGeneBrief=>['refGene', 'txStart', 'txEnd',
					'name,strand,chrom,txStart,txEnd,cdsStart,cdsEnd,exonCount,name2']);

  my %opts = (d=>'hg18', s=>'genomep:password@genome-mysql.cse.ucsc.edu:3306', p=>'refGeneBrief', S=>'');
  getopts('d:s:p:Ch', \%opts);
  &help_message if (defined $opts{h});
  &usage(\%opts, \%preset) if (@ARGV == 0 && -t STDIN && !$opts{S});
  $opts{C} = (defined $opts{C})? 1 : 0;

  # parse the server string
  $_ = $opts{s};
  s/\s//g;
  die("(EE) fail to parse the server string.\n") unless (/^(\S+):(\S*)\@(\S+)$/);
  my ($user, $passwd, $server) = ($1, $2, $3);

  # retrieving data
  my $dbh = DBI->connect_cached("dbi:mysql:$opts{d}:$server", $user, $passwd);
  &config($dbh, \%opts, \%preset);
  &retrieve($dbh, \%opts);
  $dbh->disconnect;
}

sub config {
  my ($dbh, $opts, $preset) = @_;
  if ($opts->{p} =~ /^\S+(:\S*){3}$/) {
	($opts->{t}, $opts->{beg}, $opts->{end}, $opts->{c}) = split(':', $opts->{p});
  } else {
	die("(EE) undefined preset.\n") unless (defined $preset->{$opts->{p}});
	($opts->{t}, $opts->{beg}, $opts->{end}, $opts->{c}) = @{$preset->{$opts->{p}}};
  }
  $opts->{b} = 0;
  $opts->{c} = '*' unless ($opts->{c});
  my $sth = $dbh->prepare("DESCRIBE $opts->{t}");
  $sth->execute;
  while ((@_ = $sth->fetchrow_array)) {
	if ($_[0] eq 'bin' && $_[1] =~ /int/) {
	  $opts->{b} = 1;
	} elsif ($_[0] =~ /Start$/) {
	  $opts->{beg} = $_[0] unless ($opts->{beg});
	  warn("(WW) multiple start columns: $opts->{beg} and $_[0]; use $opts->{beg}\n") if ($opts->{beg} ne $_[0]);
	} elsif ($_[0] =~ /End$/) {
	  $opts->{end} = $_[0] unless ($opts->{end});
	  warn("(WW) multiple end columns: $opts->{end} and $_[0]; use $opts->{end}\n") if ($opts->{end} ne $_[0]);
	}
  }
  die("(EE) fail to detect the column name for the start or the end of a region.\n") unless ($opts->{beg} && $opts->{end});
  $sth->finish;
}

sub retrieve {
  my ($dbh, $opts) = @_;
  my (@sth_cache, $sth, @bin);
  my $sql_comm = qq/SELECT $opts->{c} FROM $opts->{t} WHERE chrom=? AND /;
  $sql_comm .= $opts->{C}? qq/$opts->{beg}>=? AND $opts->{end}<=?/ : qq/$opts->{end}>=? AND $opts->{beg}<=?/;
  $sql_comm .= ' AND (bin=?' if ($opts->{b});
  while (<>) {
	my @t = split;
	$t[0] = "chr$t[0]" if ($opts->{d} =~ /^hg/ && $t[0] !~ /^chr/);
	if ($t[1] > $t[2]) { my $tt = $t[1]; $t[1] = $t[2]; $t[2] = $tt; }
	@bin = ();
	if ($opts->{b}) {
	  @bin = &region2bin($t[1], $t[2]);
	  $sth_cache[@bin] ||= $dbh->prepare($sql_comm . (' OR bin=?' x (@bin-1)) . ')');
	} else {
	  $sth_cache[0] ||= $dbh->prepare($sql_comm);
	}
	$sth = $sth_cache[@bin];
	$sth->execute(@t[0..2], @bin);
	while ((@_ = $sth->fetchrow_array)) {
	  print join("\t", @_), "\n";
	}
  }
  for (@sth_cache) {
	$_->finish if (defined $_);
  }
}

sub region2bin {
  my ($beg, $end) = @_;
  my @bin = (1);
  push(@bin, (  1 + ($beg>>26) ..   1 + (($end-1)>>26)));
  push(@bin, (  9 + ($beg>>23) ..   9 + (($end-1)>>23)));
  push(@bin, ( 73 + ($beg>>20) ..  73 + (($end-1)>>20)));
  push(@bin, (585 + ($beg>>17) .. 585 + (($end-1)>>17)));
  return @bin;
}

sub usage {
  my ($opts, $preset) = @_;
  my $tmp = join(", ", keys(%$preset));
  print qq{
Program: batchUCSC.pl (Batch data downloader from UCSC)
Contact: Heng Li <lh3\@sanger.ac.uk>
Usage:   batchUCSC.pl [options] <in.region_file>

Options: -s STR     server [$opts->{s}]
         -d STR     database [$opts->{d}]
         -p STR     query: ${tmp} [$opts->{p}]
         -C         contained (by default overlap)
         -h         longer help

Examples:

  batchUCSC.pl -p snpAnc inp.txt
  batchUCSC.pl -p snp129::: inp.txt
  batchUCSC.pl -p 'refGene:cdsStart:cdsEnd:*' inp.txt
  batchUCSC.pl -p simpleRepeat::: inp.txt

IMPORTANT NOTES:

  DO NOT retrieve the data from a whole chromosome, not to speak from
  the whole genome. Doing in this way is extremely slow and will put a
  lot of burden on the UCSC MySQL server. If you want the whole genome
  data, please download from UCSC's FTP site. You can also import UCSC
  data to your local MySQL database, but you need to change '-s' and
  '-d' according to your local configuration of MySQL.

};
  exit;
}

sub help_message {
  print qq{
Composing queries:

  a) batchUCSC.pl -p snp inp.txt

     You can use presetting. The above command is equivalent to (for
     explanation, see below):

     batchUCSC.pl -p 'snp129:::' inp.txt

  b) batchUCSC.pl -p 'knownGene:txStart:txEnd:*' inp.txt

     For a customized query, the query consists of four fields,
     separated by colons. The four fields are: the table name
     (required), the column name for the start of a region, name for the
     end of a region and the list of fields to be displayed. When the
     second the third fields are missing, this script will infer from
     the table schema. The default value for the last column is '*',
     which shows all the columns. You can even do calculation between
     cloumns at the last field, as long as that fields follows the SQL
     syntax.

     The list of available UCSC tables can be found on its FTP site or
     here:

     http://genome.ucsc.edu/cgi-bin/hgTables?command=start

  c) batchUCSC.pl -p 'simpleRepeat:::' inp.txt

     The above command is equivalent to:

     batchUCSC.pl -p 'simpleRepeat:chromStart:chromEnd:*' inp.txt

     If the script finds multiple start or end columns (for example in
     table knownGene), it will use the first such column in the
     table and give a warning.

};
  exit;
}
