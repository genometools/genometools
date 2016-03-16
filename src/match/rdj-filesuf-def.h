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

/* Suffix of Readjoiner output and intermediate files */

#define GT_READJOINER_SUFFIX_SPMLIST            ".spm"
#define GT_READJOINER_SUFFIX_CONTIGS            ".contigs.fas"
#define GT_READJOINER_SUFFIX_CONTIG_PATHS       ".paths"
#define GT_READJOINER_SUFFIX_READSCOPYNUM       ".rcn"
#define GT_READJOINER_SUFFIX_READSLIBRARYTABLE  ".rlt"
#define GT_READJOINER_SUFFIX_SPMCOUNTS          ".nofspm"
#define GT_READJOINER_SUFFIX_SPMCOUNTS_DISTRI   ".sc.dsr"
#define GT_READJOINER_SUFFIX_WSIZE_DISTRI       ".ws.dsr"
#define GT_READJOINER_SUFFIX_PREFILTERED_FAS    ".p.fas"
#define GT_READJOINER_SUFFIX_SEQNUMS            ".sn"
#define GT_READJOINER_SUFFIX_CNTLIST            ".cnt"
#define GT_READJOINER_SUFFIX_CONNECTING_PATHS   ".connect.fas"

/* string graph */
#define GT_READJOINER_SUFFIX_SG                 ".sg"
#define GT_READJOINER_SUFFIX_SG_MONO_DOT        ".m.dot"
#define GT_READJOINER_SUFFIX_SG_BI_DOT          ".b.dot"
#define GT_READJOINER_SUFFIX_SG_SUB_DOT         ".sub.dot"
#define GT_READJOINER_SUFFIX_SG_SPMLIST         ".sg.spm"
#define GT_READJOINER_SUFFIX_SG_ADJLIST         ".adj"
#define GT_READJOINER_SUFFIX_SG_ELEN_DISTRI     ".el.dsr"
#define GT_READJOINER_SUFFIX_SG_ASQG            ".asqg"
#define GT_READJOINER_SUFFIX_SG_ASQG_GZ         ".asqg.gz"

/* contig graph */
#define GT_READJOINER_SUFFIX_JUNCTIONS          ".jnc"
#define GT_READJOINER_SUFFIX_CJ_I_LINKS         ".cji"
#define GT_READJOINER_SUFFIX_CJ_O_LINKS         ".cjo"
#define GT_READJOINER_SUFFIX_DEPTHINFO          ".dpt"
#define GT_READJOINER_SUFFIX_CG_DOT             ".cg.dot"
#define GT_READJOINER_SUFFIX_CG_SUB_DOT         ".cg.sub.dot"
#define GT_READJOINER_SUFFIX_CG_PATHS           ".cg.paths"

#endif
