#!/usr/local/bin/perl

use strict;

my $values=join('',<STDIN>);
$values=~s/[\s\n]//g;
my ($x0,$x1);
my $n=0;
print "int [][2]={\n";
for (split(/,/,$values)) {
  if(m/^\s*"\\u([0-9A-F]+)"\s*$/i ) {
    $x0=hex($1); $x1=$x0;
  } elsif(/^\s*"\\u([0-9A-F]+)"\s*-\s*"\\u([0-9A-F]+)"\s*$/i) {
    $x0=hex($1); $x1=hex($2);
  } else {
    next;
  }
  printf "  {0x%x,0x%x},\n",$x0,$x1;
  ++$n;
}
print "};\n";
