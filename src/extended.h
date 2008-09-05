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

#ifndef EXTENDED_H
#define EXTENDED_H

/* The extended GenomeTools library header */
#include "extended/add_introns_stream.h"        /* add introns stream */
#include "extended/affinealign.h"               /* affine align module */
#include "extended/align.h"                     /* align module */
#include "extended/alignment.h"                 /* alignment class */
#include "extended/cds_stream.h"                /* CDS stream */
#include "extended/consensus_sa.h"              /* consensus spl. align. mod. */
#include "extended/csa_stream.h"                /* consensus spl. align. str. */
#include "extended/coin_hmm.h"                  /* the coin HMM class */
#include "extended/compare.h"                   /* compare module */
#include "extended/dice_hmm.h"                  /* the dice HMM class */
#include "extended/extract_feat_stream.h"       /* extract feat. stream class */
#include "extended/evaluator.h"                 /* evaluator class */
#include "extended/feature_type_factory_any.h"  /* feature type factory class */
#include "extended/filter_stream.h"             /* filter stream class */
#include "extended/genome_stream.h"             /* genome stream class */
#include "extended/gff3_in_stream.h"            /* GFF3 input stream class */
#include "extended/gff3_out_stream.h"           /* GFF3 output stream class */
#include "extended/gtf_in_stream.h"             /* GTF input stream class */
#include "extended/gtf_out_stream.h"            /* GTF output stream class */
#include "extended/hmm.h"                       /* HMM class */
#include "extended/linearalign.h"               /* linear alignment module */
#include "extended/linearedist.h"               /* linear edit distance mod. */
#include "extended/merge_stream.h"              /* merge stream class */
#include "extended/mergefeat_stream_sorted.h"   /* merge feat. stream class */
#include "extended/mergefeat_stream_unsorted.h" /* merge feat. stream class */
#include "extended/msa.h"                       /* multiple seq. align. class */
#include "extended/multiset_matching.h"         /* multiset matching module */
#include "extended/mutate.h"                    /* mutate module */
#include "extended/neighborjoining.h"           /* the Neighbor-Joining class */
#include "extended/qgramdist.h"                 /* q-gram distance module */
#include "extended/regioncov_visitor.h"         /* regioncov visitor class */
#include "extended/seqid2file.h"                /* seqid2file module */
#include "extended/sort_stream.h"               /* sort stream class */
#include "extended/splicedseq.h"                /* splicedseq class */
#include "extended/stat_visitor.h"              /* status visitor class */
#include "extended/stream_evaluator.h"          /* string class */
#include "extended/swalign.h"                   /* Smith-Waterman align. mod. */
#include "extended/toolbox.h"                   /* toolbox class */
#include "extended/upgma.h"                     /* UPGMA class */

#endif
