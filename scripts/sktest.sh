#!/bin/sh
#
# Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

set -e -x

USAGE="Usage: $0 [-memcheck]"
program="./testsuite.rb -threads 2"

if test $# -eq 1
then
  if test "$1" = "-memcheck"
  then
    program="${program} -memcheck" 
  else
    echo ${USAGE}
    exit 1
  fi 
fi

# optional -memcheck   (run valgrind)
#          -select 253 (run testcase 253)

startdate="`date`"
cd testsuite
env -i GT_MEM_BOOKKEEPING=on ${program} -keywords 'gt_encseq'

env -i GT_MEM_BOOKKEEPING=on ${program} -keywords 'gt_suffixerator'

env -i GT_MEM_BOOKKEEPING=on ${program} -keywords 'gt_encseq2spm'

env -i GT_MEM_BOOKKEEPING=on \
       ${program} -keywords 'gt_suffixerator and gttestdata' \
       -gttestdata ${GTTESTDATA}

env -i GT_MEM_BOOKKEEPING=on ${program} -keywords 'gt_sortbench'

env -i GT_MEM_BOOKKEEPING=on ${program} -keywords 'gt_extractseq' \
       -gttestdata ${GTTESTDATA}

env -i GT_MEM_BOOKKEEPING=on ${program} -keywords 'gt_checkprjfiles' \
       -gttestdata ${GTTESTDATA}

env -i GT_MEM_BOOKKEEPING=on ${program} -keywords 'gt_trieins'

env -i GT_MEM_BOOKKEEPING=on ${program} -keywords 'gt_mergeesa'

env -i GT_MEM_BOOKKEEPING=on ${program} -keywords 'gt_packedindex' \
       -gttestdata ${GTTESTDATA}

cd ..

sktest-match.sh
echo "start at ${startdate}"
echo "end at `date`"
