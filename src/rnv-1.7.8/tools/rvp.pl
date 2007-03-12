#!/usr/local/bin/perl
# $Id: rvp.pl,v 1.6 2004/01/14 16:42:38 dvd Exp $

# embedding sample for RVP, a part of RNV, http://davidashen.net/rnv.html
# code kept simple to show the technique, not to provide a general purpose
# module.
#
# details of the protocol are in a long comment near the start of rvp.c
#

use FileHandle; use IPC::Open2;
use XML::Parser::Expat;

use strict;

my @RVP=("rvp");

my ($parser,$errors,$prevline,$prevcol); # declared here for use in resp()

# rvp wrapper
$|=1;  # using pipes
$/="\0"; # the zero byte is separator

# write queries to RVPIN, get responses from RVPOUT
open2(\*RVPOUT,\*RVPIN,@RVP,@ARGV);

sub resp {
  $_=<RVPOUT>;
  chop;
  /^ok (\d+).*/ and return $1;
  /^error (\d+) (\d+) (.*)/ and do {
      my ($pat,$msg)=($1,$3);

     # if the message is empty, don't print the message,
     # the error has occured in already erroneous state
      if($msg) {
        my ($line,$col)=($parser->current_line(),$parser->current_column());
	if($line!=$prevline || $col!=$prevcol) {
          printf STDERR "%u,%u: %s\n",$line,$col,$msg;
	  $prevline=$line; $prevcol=$col;
	}
      }
      $errors=1;
      return $pat;
    };
  /^er (\d+).*/ and do {$errors=1; return $1;};
  die "protocol error ($_), stopped ";
}

sub start_tag_open {
  my($cur,$name)=@_;
  syswrite RVPIN,"start-tag-open $cur $name\0";
  return resp();
}

sub attribute {
  my($cur,$name,$val)=@_;
  syswrite RVPIN,"attribute $cur $name $val\0";
  return resp();
}

sub start_tag_close {
  my($cur,$name)=@_;
  syswrite RVPIN,"start-tag-close $cur $name\0";
  return resp();
}

sub end_tag {
  my($cur,$name)=@_;
  syswrite RVPIN,"end-tag $cur $name\0";
  return resp();
}

sub text {
  my($cur,$text)=@_;
  syswrite RVPIN,"text $cur $text\0";
  return resp();
}

# in mixed content, whitespace is simply discarded, and any
# non-whitespace is equal; but this optimization gives only only
# 5% increase in speed at most in practical cases
sub mixed {
  my($cur,$text)=@_;
  if($text=~m/[^\t\n ]/s) {
    syswrite RVPIN,"mixed $cur .\0";
    return resp();
  } else {
    return $cur;
  }
}

sub start {
  my $no=shift @_ or 0;
  syswrite RVPIN,"start $no\0";
  return resp();
}

sub quit {
  syswrite RVPIN,"quit\0";
  return resp();
}

my ($current,$text,$mixed)=(0,"",0);

# Expat handlers

# Expat does not merge cdata into one text node;
# application has to do it explicitely (see characters())
sub flush_text {
  $current=$mixed?mixed($current,$text):text($current,$text);
  $text="";
}

# last colon in a name separates the local name from the URI
sub qname {
   my ($p,$name)=@_;
   return join(':',$p->namespace($name),$name);
}

sub start_element {
  my ($p,$el,%attrs)=@_;

  $mixed=1;
  flush_text();
  $current=start_tag_open($current,qname($p,$el));
  $mixed=0;
  for my $atname (keys %attrs) {
    $current=attribute($current,qname($p,$atname),$attrs{$atname});
  }
  $current=start_tag_close($current,qname($p,$el));
}

sub end_element {
  my ($p,$el)=@_;
  flush_text();
  $current=end_tag($current,qname($p,$el));
  $mixed=1;
}

sub characters {
  my ($p,$chars)=@_;
  $text=$text.$chars;
}

# Main

$errors=0;
$parser=new XML::Parser::Expat(
  Namespaces=>1);

$parser->setHandlers(
  Start=>\&start_element,
  End=>\&end_element,
  Char=>\&characters);

$current=start(0);
$parser->parse(*STDIN);
quit();

exit($errors);
