/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef CHAINING_H
#define CHAINING_H

#include "core/unused_api.h"
#include "gth/call_info.h"
#include "gth/chain_collection.h"
#include "gth/gthmatch.h"
#include "gth/input.h"
#include "gth/matcher.h"
#include "gth/plugins.h"

void gth_chaining(GthChainCollection *chain_collection,
                  unsigned long gen_file_num,
                  unsigned long ref_file_num,
                  GthCallInfo*,
                  GthInput*,
                  GthStat*,
                  bool directmatches,
                  const GthPlugins *plugins);

typedef struct {
  bool directmatches,
       refseqisindex, /* (inverse || !refseqisdna) */
       jtdebug;
  GthCallInfo *call_info;
  GthInput *input;
  GthStat *stat;
  unsigned long gen_file_num,
                ref_file_num,
                bucketnum,
                maxbucketnum;
} GthChainingInfo;

typedef struct {
  GthChainCollection *chain_collection;
  GtArray *matches;
  bool directmatches,
       refseqisdna,
       online,
       refseqisindex; /* (inverse || !refseqisdna) */
  GthStat *stat;
  GthChainingInfo *chaining_info;
  unsigned long *matchnumcounter,
                maxnumofmatches,
                rare,
                lastrefseqnum;
  double fragweightfactor;

  GthJumpTableNew jump_table_new;
  GthJumpTableNewReverse jump_table_new_reverse;
  GthJumpTableDelete jump_table_delete;

  /* can be filled and used by matcher */
  GthSeqCon *gen_seq_con,
            *ref_seq_con;
} GthMatchProcessorInfo;

/* the matcher has to call this */
int gth_match_processor(GthMatchProcessorInfo *info, GthSeqCon *gen_seq_con,
                        GthSeqCon *ref_seq_con, GthMatch *match);

#endif
