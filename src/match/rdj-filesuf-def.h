/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_FILESUF_DEF_H
#define RDJ_FILESUF_DEF_H

/*
 * Filename suffices relevant to Readjoiner are defined here
 *
 */

/* -- main output -- */

#define GT_READJOINER_SUFFIX_SPMLIST           ".spm"
#define GT_READJOINER_SUFFIX_CONTIG_PATHS      ".paths"
#define GT_READJOINER_SUFFIX_CONTIGS           ".contigs.fas"

/* -- prefilter -- */

#define GT_READJOINER_SUFFIX_PREFILTERED_FAS   ".pf.fas"
#define GT_READJOINER_SUFFIX_READSCOPYNUM      ".rcn"

/* -- test / debug / development output -- */

#define GT_READJOINER_SUFFIX_CNTLIST           ".cnt"
#define GT_READJOINER_SUFFIX_SPMCOUNTS         ".nofspm"
#define GT_READJOINER_SUFFIX_SPMCOUNTS_DISTRI  ".nofspm.distri"
#define GT_READJOINER_SUFFIX_ELEN_DISTRI       ".elen.distri"
#define GT_READJOINER_SUFFIX_SEPPOS            ".sep"
#define GT_READJOINER_SUFFIX_TWOBIT            ".2bit"
#define GT_READJOINER_SUFFIX_SEQNUMS           ".seqnums"
#define GT_READJOINER_SUFFIX_WSIZE_DISTRI      ".wsize.distri"

/* string graph: */

/* before transitive reduction (u=unreduced) */
#define GT_READJOINER_SUFFIX_U_SG              ".u.sg"
#define GT_READJOINER_SUFFIX_U_SG_DOT          ".u.dot"
#define GT_READJOINER_SUFFIX_U_SG_SPMLIST      ".u.sg.spm"
#define GT_READJOINER_SUFFIX_U_SG_ADJLIST      ".u.adj"

/* without submax edges reduction (m=multigraph) */
#define GT_READJOINER_SUFFIX_M_SG              ".m.sg"
#define GT_READJOINER_SUFFIX_M_SG_DOT          ".m.dot"
#define GT_READJOINER_SUFFIX_M_SG_SPMLIST      ".m.sg.spm"
#define GT_READJOINER_SUFFIX_M_SG_ADJLIST      ".m.adj"

/* reduced graph */
#define GT_READJOINER_SUFFIX_SG                ".sg"
#define GT_READJOINER_SUFFIX_SG_DOT            ".dot"
#define GT_READJOINER_SUFFIX_SG_SPMLIST        ".sg.spm"
#define GT_READJOINER_SUFFIX_SG_ADJLIST        ".adj"

/* after error correction */
#define GT_READJOINER_SUFFIX_C_SG              ".c.sg"
#define GT_READJOINER_SUFFIX_C_SG_DOT          ".c.dot"
#define GT_READJOINER_SUFFIX_C_SG_SPMLIST      ".c.sg.spm"
#define GT_READJOINER_SUFFIX_C_SG_ADJLIST      ".c.adj"

#endif
