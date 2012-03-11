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
chain2dim.c
chainofin.c
cutendpfx.c
dist-short.c
echoseq.c
eis-blockcomp-param.c
eis-blockcomp.c
eis-bwtseq-extinfo.c
eis-bwtseq-param.c
eis-bwtseq-sass.c
eis-bwtseq.c
eis-encidxseq-param.c
eis-encidxseq.c
eis-sa-common.c
eis-seqblocktranslate.c
eis-seqranges.c
eis-sequencemultiread.c
eis-specialsrank.c
eis-suffixarray-interface.c
eis-suffixerator-interface.c
eis-voiditf.c
encodedseq.c
enum-patt.c
esa-bottomup.c
esa-dfs.c
esa-lcpintervals.c
esa-lcpval.c
esa-map.c
esa-maxpairs.c
esa-merge.c
esa-mmsearch.c
esa-ppbuckwid.c
esa-scanprj.c
esa-seqread.c
esa-shulen.c
esa-splititv.c
esa-spmitvs.c
esa-spmsk.c
esa_lcpintervals_visitor.c
esa_spmitvs_visitor.c
esa_visitor.c
firstcodes-accum.c
firstcodes-insert.c
firstcodes-psbuf.c
firstcodes-scan.c
firstcodes-spacelog.c
firstcodes-tab.c
firstcodes.c
giextract.c
greedyedist.c
greedyfwdmat.c
hashfirstcodes.c
idx-limdfs.c
idxlocali.c
idxlocalidp.c
idxlocalisw.c
index_options.c
initbasepower.c
initeqsvec.c
inl-encseq.c
iter-window.c
kmer2string.c
lua_tools.c
mapspec-gen.c
marksubstring.c
measure-time.c
merger-trie.c
mssufpat.c
myersapm.c
nullcols.c
opensfxfile.c
optionargmode.c
pck-count-nodes.c
pckbucket.c
pckdfs.c
pos2seqnum.c
prsqualint.c
pssm.c
querymatch.c
rdj-cntlist.c
rdj-contfind-bottomup.c
rdj-contfinder.c
rdj-contigpaths.c
rdj-contigs-writer.c
rdj-errfind.c
rdj-ovlfind-bf.c
rdj-ovlfind-dp.c
rdj-ovlfind-gusfield.c
rdj-ovlfind-kmp.c
rdj-pairwise.c
rdj-radixsort.c
rdj-spmfind.c
rdj-spmlist.c
rdj-spmproc.c
rdj-ssar.c
rdj-strgraph.c
rdj-strgraph.c
revcompl.c
seqnumrelpos.c
sfx-apfxlen.c
sfx-bentsedg.c
sfx-bltrie.c
sfx-copysort.c
sfx-diffcov.c
sfx-enumcodes.c
sfx-input.c
sfx-lcpvalues.c
sfx-linlcp.c
sfx-mappedstr.c
sfx-maprange.c
sfx-outprj.c
sfx-opt.c
sfx-partssuf.c
sfx-progress.c
sfx-readint.c
sfx-readmode.c
sfx-remainsort.c
sfx-run.c
sfx-shortreadsort.c
sfx-suffixer.c
sfx-suffixgetset.c
sfx-suftaborder.c
shu-dfs.c
shu-divergence.c
shu-encseq-gc.c
shu-genomediff-pck-simple.c
shu-genomediff.c
shu_unitfile.c
spaced-seeds.c
specialrank.c
spmsuftab.c
squarededist.c
substriter.c
tagerator.c
test-pairwise.c
turnwheels.c
twobits2kmers.c
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

Rm -f ${TMPFILE}
