#!/bin/sh
if test $# -ne 1
then
  echo "Usage: $0 <gzipped trembl file>"
  exit 1
fi
gunzip -c $1 | grep '^>' | sed -e 's/>tr\|\(.*\)\|.*/\1/g'
