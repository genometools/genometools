/*
  Copyright (c) 2003-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef GTHCHAIN_H
#define GTHCHAIN_H

#include "core/array_api.h"
#include "core/file.h"
#include "core/range_api.h"
#include "extended/globalchaining.h"
#include "gth/jump_table.h"
#include "gth/stat.h"

/*
  A GthChain is an array of genomic ranges (denoting the potential exons)
  plus some additional information
*/

typedef struct {
  unsigned long gen_file_num,
                gen_seq_num,
                ref_file_num,
                ref_seq_num;
  GtArray *forwardranges,  /* refering forward strand */
          *reverseranges;  /* refering to reverse complement strand */

  double refseqcoverage; /* the reference sequence coverage of this DP range */
  GthJumpTableDelete jump_table_delete;
  GthJumpTable *forward_jump_table,
               *reverse_jump_table;
} GthChain;

#include "gth/chain_collection.h"

GthChain* gth_chain_new(void);
void      gth_chain_delete(GthChain*);
void      gth_chain_copy(GthChain*, const GthChain*);
void      gth_chain_extend_borders(GthChain *dprange,
                                   const GtRange *gen_seq_bounds,
                                   const GtRange *gen_seq_bounds_rc,
                                   unsigned long gen_total_length,
                                   unsigned long gen_offset);
/* shorten the potential introns contained in the chain */
void      gth_chain_shorten_introns(GthChain *dprange,
                                    unsigned long icdelta,
                                    unsigned long icminremintronlength,
                                    unsigned long gen_total_length,
                                    unsigned long gen_offset,
                                    bool comments,
                                    GtFile *outfp);
void      gth_chain_contract(GthChain*, const GthChain*);

typedef struct {
  GthChainCollection *chain_collection;
  unsigned long gen_file_num,
                gen_seq_num,
                ref_file_num,
                ref_seq_num,
                numofremainingbuckets,
                gen_total_length,
                gen_offset,
                ref_total_length,
                ref_offset,
                referencelength,   /* needed by globalchainislongenough() */
                gcmincoverage;     /* needed by globalchainislongenough() */
  bool comments,
       stopafterchaining,
       directmatches,     /* needed by stopafterchaining */
       paralogs,
       enrichchains,
       jump_table;
  GthJumpTableNew jump_table_new;
  GthJumpTableNewReverse jump_table_new_reverse;
  GthJumpTableDelete jump_table_delete;
  bool jtdebug;
  GthStat *stat;
  GtFile *outfp;
} GthSaveChainInfo;

void gth_save_chain(GtChain*, GtFragment*, unsigned long, unsigned long, void*);

#endif
