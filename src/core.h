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

#ifndef CORE_H
#define CORE_H

/* The core GenomeTools library  header */
#include "core/array2dim.h"      /* 2-dimensional array class */
#include "core/bioseq.h"         /* biosequence class */
#include "core/bittab.h"         /* bittab class */
#include "core/bitpackarray.h"   /* packed integer array class */
#include "core/bitpackstring.h"  /* packed integer module */
#include "core/bsearch.h"        /* bsearch module */
#include "core/countingsort.h"   /* countingsort module */
#include "core/cstr.h"           /* C-string class */
#include "core/disc_distri.h"    /* discrete distribution class */
#include "core/dlist.h"          /* double-linked list class */
#include "core/ensure.h"         /* defines the ensure macro */
#include "core/error.h"          /* error class */
#include "core/fasta.h"          /* fasta module */
#include "core/fileutils.h"      /* file utilities module */
#include "core/gc_content.h"     /* GC-content module */
#include "core/grep.h"           /* grep module */
#include "core/hashtable.h"      /* hashtable class */
#include "core/mathsupport.h"    /* math support module */
#include "core/minmax.h"         /* min/max module */
#include "core/msort.h"          /* msort module */
#include "core/option.h"         /* option parser class */
#include "core/outputfile.h"     /* outputfile module */
#include "core/orf.h"            /* ORF module */
#include "core/parseutils.h"     /* parse utilities module */
#include "core/phase.h"          /* phase module */
#include "core/progressbar.h"    /* progressbar module */
#include "core/queue.h"          /* queue class */
#include "core/range.h"          /* range class */
#include "core/safearith.h"      /* safe arithmetics module */
#include "core/score_function.h" /* score function class */
#include "core/score_matrix.h"   /* score matrix class */
#include "core/splitter.h"       /* splitter class */
#include "core/str.h"            /* string class */
#include "core/strand.h"         /* strand module */
#include "core/tokenizer.h"      /* tokenizer class */
#include "core/translate.h"      /* translate module */
#include "core/undef.h"          /* undef module */
#include "core/unused.h"         /* unused module */
#include "core/versionfunc.h"    /* version module */
#include "core/warning.h"        /* warning module */
#include "core/xansi.h"          /* ANSI wrapper module */
#include "core/xposix.h"         /* POSIX wrapper module */

#endif
