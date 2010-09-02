#!/bin/sh

set -e -x

if test $# -ne 1
then
  echo "Usage: $0 <runfileprefix>"
  exit 1
fi

prefix=$1

echo "prefix=${prefix}"

${GTDIR}/scripts/run2tex.rb time+space ${prefix}.run > ${prefix}.tex
latex ${prefix}.tex
dvips -P www -o ${prefix}.ps ${prefix}.dvi
