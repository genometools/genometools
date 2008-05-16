/*
  Copyright (c) 2003-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef GTCORE_H
#define GTCORE_H

/* The core GenomeTools library (libgtcore) header */
#include "libgtcore/alpha.h"          /* alphabet class */
#include "libgtcore/array.h"          /* array class */
#include "libgtcore/array2dim.h"      /* 2-dimensional array class */
#include "libgtcore/bioseq.h"         /* biosequence class */
#include "libgtcore/bittab.h"         /* bittab class */
#include "libgtcore/bitpackarray.h"   /* packed integer array class */
#include "libgtcore/bitpackstring.h"  /* packed integer module */
#include "libgtcore/bsearch.h"        /* bsearch module */
#include "libgtcore/countingsort.h"   /* countingsort module */
#include "libgtcore/cstr.h"           /* C-string class */
#include "libgtcore/disc_distri.h"    /* discrete distribution class */
#include "libgtcore/dlist.h"          /* double-linked list class */
#include "libgtcore/ensure.h"         /* defines the ensure macro */
#include "libgtcore/error.h"          /* error class */
#include "libgtcore/fasta.h"          /* fasta module */
#include "libgtcore/fileutils.h"      /* file utilities module */
#include "libgtcore/gc_content.h"     /* GC-content module */
#include "libgtcore/grep.h"           /* grep module */
#include "libgtcore/hashtable.h"      /* hashtable class */
#include "libgtcore/mathsupport.h"    /* math support module */
#include "libgtcore/minmax.h"         /* min/max module */
#include "libgtcore/msort.h"          /* msort module */
#include "libgtcore/option.h"         /* option parser class */
#include "libgtcore/outputfile.h"     /* outputfile module */
#include "libgtcore/orf.h"            /* ORF module */
#include "libgtcore/parseutils.h"     /* parse utilities module */
#include "libgtcore/phase.h"          /* phase module */
#include "libgtcore/progressbar.h"    /* progressbar module */
#include "libgtcore/queue.h"          /* queue class */
#include "libgtcore/range.h"          /* range class */
#include "libgtcore/safearith.h"      /* safe arithmetics module */
#include "libgtcore/score_function.h" /* score function class */
#include "libgtcore/score_matrix.h"   /* score matrix class */
#include "libgtcore/splitter.h"       /* splitter class */
#include "libgtcore/str.h"            /* string class */
#include "libgtcore/strand.h"         /* strand module */
#include "libgtcore/tokenizer.h"      /* tokenizer class */
#include "libgtcore/translate.h"      /* translate module */
#include "libgtcore/undef.h"          /* undef module */
#include "libgtcore/unused.h"         /* unused module */
#include "libgtcore/versionfunc.h"    /* version module */
#include "libgtcore/warning.h"        /* warning module */
#include "libgtcore/xansi.h"          /* ANSI wrapper module */
#include "libgtcore/xposix.h"         /* POSIX wrapper module */

#endif
