#!/bin/sh
#
# Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

localpath=src/match

if test "XX${GT_ENV_OPTIONS}" != "XX"
then
  unset GT_ENV_OPTIONS
fi

TMPFILE=`mktemp /tmp/skproto-all.XXXXXX` || exit 1
cat << END_OF_TEXT > ${TMPFILE}
alphabet.c
apmeoveridx.c
bcktab.c
cutendpfx.c
cgr_spacedseed.c
dist-short.c
eis-blockcomp.c
eis-blockcomp-param.c
eis-blockenc_params.c
eis-bwtconstruct_params.c
eis-bwtseq.c
eis-bwtseqcreate.c
eis-bwtseq-extinfo.c
eis-bwtseqlocate.c
eis-bwtseq-param.c
eis-bwtseq-sass.c
eis-encidxseq.c
eis-encidxseq-param.c
eis-sa-common.c
eis-seqblocktranslate.c
eis-seqranges.c
eis-sequencemultiread.c
eis-suffixarray-interface.c
eis-suffixerator-interface.c
eis-suffixeratorreader.c
eis-voiditf.c
encodedseq.c
encseq-specialsrank.c
enum-patt.c
esa-dfs.c
esa-lcpval.c
esa-map.c
esa-mmsearch.c
esa-ppbuckwid.c
esa-seqread.c
esa-splititv.c
echoseq.c
fillsci.c
greedyfwdmat.c
idx-limdfs.c
idxlocali.c
idxlocalidp.c
initeqsvec.c
inl-encseq.c
iter-window.c
mapspec-gen.c
measure-time.c
merger-trie.c
mssufpat.c
myersapm.c
optionargmode.c
overallseq.c
opensfxfile.c
pckbucket.c
pos2seqnum.c
pssm.c
nullcols.c
revcompl.c
sfx-bentsedg.c
sfx-enumcodes.c
sfx-input.c
sfx-outlcp.c
sfx-partssuf.c
sfx-readint.c
sfx-run.c
sfx-suffixer.c
spaced-seeds.c
substriter.c
tagerator.c
test-pairwise.c
trie-ssort.c
turnwheels.c
tyr-mkindex.c
tyr-map.c
tyr-mersplit.c
tyr-occratio.c
tyr-search.c
uniquesub.c
verbose.c
END_OF_TEXT

# SKPROTO=./bin/skproto
SKPROTO='./bin/gt dev skproto'

if test ! -f ./bin/gt
then
  echo "$0: ${SKPROTO} does not exist"
  exit 1
fi

for filename in `ls ${localpath}/*.c | grep -v -f ${TMPFILE}`
do
prfile="${localpath}/`basename ${filename} .c`.pr"
  if test ! -f ${prfile} ||
     test ! -s ${prfile} ||
     test ${filename} -nt ${prfile}
  then
    echo "create ${prfile}"
    ${SKPROTO} ${filename} > ${prfile}
  fi
done

rm -f ${TMPFILE}
