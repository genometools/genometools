#!/usr/bin/env perl

use strict;
use warnings;

unless(@ARGV)
{
  print STDERR "Usage: $0 <inputfile>\n";
  exit 1;
}

if((scalar @ARGV) eq 0)
{
  print STDERR "$0: missing argument\n";
  exit 1;
}

my($inputfile) = $ARGV[0];

unless ( -e $inputfile)
{
  print STDERR "$0: file \"$inputfile\" does not exist\n";
  exit 1;
}

unless ( open(INPUTFILEHANDLE, $inputfile) )
{
  print STDERR "Cannot open file \"$inputfile\"";
  exit 1;
}

while (my $in = <INPUTFILEHANDLE>)
{
  if($in =~ m/^\\EXECUTE\{([^\}]*)\}/)
  {
    runtheprogram($1);
  } else
  {
    print $in;
  }
}

sub determinemaxfilewidth
{
  my ($filename) = @_;
  my $maxwidth = 0;
  my $width;

  unless(open(FILEHANDLE,$filename))
  {
    print STDERR "Cannot open file \"$filename\"\n";
    exit 1;
  }
  while(my $line = <FILEHANDLE>)
  {
    $width = length $line;
    if($maxwidth < $width)
    {
      $maxwidth = $width;
    }
  }
  close(FILEHANDLE);
  return $maxwidth;
}

sub runtheprogram
{
  my($argcommand) = @_;
  my(@argv) = split(' ',$argcommand);

  my $program = $argv[0];
  my $outfileprefix = $argcommand;
  $outfileprefix =~ s/ \|.*//;   # delete possible pipe at end
  $outfileprefix =~ s/\s+/-/g;     # delete white spaces
  $outfileprefix =~ s/\*//g;     # delete white spaces
  $outfileprefix =~ s/\'.*\'//g;     # ' ' expressions
  my $outfilename = "output/${outfileprefix}" . ".out";

  if(-f ${outfilename})
  {
    print STDERR "# file ${outfilename} already exists\n";
  } else
  {
    my $fh;
    my @skippipe = split(/\|/,$argcommand);
    unless(open($fh, ">" . $outfilename))
    {
      print STDERR "Cannot open file $outfilename\n";
      exit 1;
    }
    print $fh "\$ $skippipe[0]\n";
    close($fh);
    push(@argv,'>>');
    push(@argv,$outfilename);
    my $argstring = "../../bin/" . join(' ',@argv);
    print STDERR "run $argstring\n";
    my($retcode) = system($argstring);
    $retcode = $? >> 8;
    if($retcode ne 0)
    {
      print STDERR "failure: $argstring\n";
      unlink($outfilename);
      exit 1;
    }
  }
  my $maxwidth = determinemaxfilewidth($outfilename);
  my $env;
  if($maxwidth < 78)
  {
    $env = "footnotesize";
  } else
  {
    $env = "scriptsize";
  }
  print "\\begin{$env}";
  print "\\verbatiminput{$outfilename}";
  print "\\end{$env}\n";
}
