/*
  Copyright (c) 2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef UNIQUE_ENCSEQ_REP_H
#define UNIQUE_ENCSEQ_REP_H

#include "core/encseq_api.h"
#include "core/hashmap_api.h"
#include "core/logger_api.h"
#include "match/xdrop.h"
#include "match/sfx-mappedstr.h"
#include "extended/alignment.h"
#include "core/score_matrix.h"
#include "core/score_function.h"
#include "extended/editscript.h"
#include "core/str_api.h"

typedef struct
{
  GtStr       *indexname_option;
  GtUword      minalignlength_option,
               nkmers_option,
               udbsize_option;
  GtWord       xdrop_option;
  unsigned int kmersize_option,
               nhits_option,
               windowsize_option;
  int          arbitscores_del_option,
               arbitscores_ins_option,
               arbitscores_mat_option,
               arbitscores_mis_option;
  bool         debug_logger_option;
} GtUniqueEncseqArguments;

/* TODO do not use a pointer to this in the union, use the value itself */
typedef struct
{
  GtUword compressed_start;
} GtUniqueEncseqUniqueEntry;

typedef enum {
  Unique,
  Link
} GtUniqueEncseqEntrytype;

typedef struct {
  GtEditscript *editscript;
  GtUword       aligned_seq_idx,
                offset;
} GtUniqueEncseqLinkEntry;

typedef struct {
  union {
    GtUniqueEncseqUniqueEntry *unique;
    GtUniqueEncseqLinkEntry *link;
  } entry;
  GtUword orig_endpos,
          orig_startpos;
  GtUniqueEncseqEntrytype type;
/* TODO do not use a pointer to entry in the union, use the value itself */
} GtUniqueEncseqDBentry;

typedef struct {
  GtUniqueEncseqDBentry *fragmentdb;
  GtUword               *sde,
                        *ssp;
  char                  *desc;
  GtUword                cumulength,
                         maxelems,
                         nelems,
                         nseq;
} GtUniqueEncseqDB;

typedef struct
{
  GtAlphabet              *alphabet;
  GtHashmap               *hashmap;
  GtKmercodeiterator      *kmercodeitAdd,
                          *kmercodeitMain;
  GtLogger                *logger;
  GtSeqabstract           *useq,
                          *vseq;
  GtStr                   *indexname;
  GtUniqueEncseqDB        *uniqueencseqdb;
  GtXdropArbitraryscores  *arbitscores;
  GtXdropresources        *res;
  char                    *buffer;
  const GtEncseq          *encseq;
  const GtKmercode        *kmercodeAdd,
                          *kmercodeExt,
                          *kmercodeMain;
  const GtUchar           *characters;
  GtXdropscore             xdropbelowscore;
  unsigned int             kmersize,
                           maxkmerhits,
                           minkmerhits,
                           numofchars,
                           windowsize;
  GtUword                  alignmentcount,
                           currentposition,
                           initKmerCounter,
                           initUniqueDBsize,
                           kmercount,
                           kmerhitcount,
                           maxlen,
                           minalignlength,
                           nkmers,
                           nPosSinceInsert,
                           nSequences,
                           nextUniqueStartpos,
                           seqlen,
                           seqnum,
                           seqstartpos,
                           totallength,
                           unique_cumulen;
} GtUniqueEncseqInfo;

#endif
