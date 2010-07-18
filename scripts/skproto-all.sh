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
cgr_spacedseed.c
cutendpfx.c
chainofin.c
chain2dim.c
dist-short.c
echoseq.c
eis-blockcomp.c
eis-blockcomp-param.c
eis-bwtseq.c
eis-bwtseq-extinfo.c
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
eis-voiditf.c
encodedseq.c
encseq-specialsrank.c
enum-patt.c
esa-dfs.c
esa-maxpairs.c
esa-merge.c
esa-lcpval.c
esa-map.c
esa-mmsearch.c
esa-ppbuckwid.c
esa-seqread.c
esa-splititv.c
giextract.c
greedyfwdmat.c
idx-limdfs.c
idxlocali.c
idxlocalidp.c
idxlocalisw.c
initbasepower.c
initeqsvec.c
inl-encseq.c
iter-window.c
kmer2string.c
mapspec-gen.c
measure-time.c
merger-trie.c
mssufpat.c
myersapm.c
nullcols.c
opensfxfile.c
optionargmode.c
pckbucket.c
pos2seqnum.c
pssm.c
prsqualint.c
querymatch.c
revcompl.c
sfx-apfxlen.c
sfx-bentsedg.c
sfx-bltrie.c
sfx-diffcov.c
sfx-copysort.c
sfx-enumcodes.c
sfx-input.c
sfx-linlcp.c
sfx-mappedstr.c
sfx-partssuf.c
sfx-progress.c
sfx-readint.c
sfx-readmode.c
sfx-remainsort.c
sfx-run.c
sfx-suffixer.c
sfx-suftaborder.c
spaced-seeds.c
specialrank.c
substriter.c
tagerator.c
test-pairwise.c
turnwheels.c
tyr-map.c
tyr-mersplit.c
tyr-mkindex.c
tyr-occratio.c
tyr-search.c
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
