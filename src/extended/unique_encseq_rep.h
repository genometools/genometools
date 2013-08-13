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
  unsigned int kmersize_option;
  bool debug_logger_option;
  unsigned int windowsize_option;
  unsigned int nhits_option;
  GtUword minalignlength_option;
  int arbitscores_mat_option;
  int arbitscores_mis_option;
  int arbitscores_ins_option;
  int arbitscores_del_option;
  GtWord xdrop_option;
  GtUword udbsize_option;
  GtStr *indexname_option;
  GtUword nkmers_option;
} GtUniqueEncseqArguments;

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
  GtUword offset;
  GtUword aligned_seq_idx;
} GtUniqueEncseqLinkEntry;

typedef struct {
  GtUword orig_startpos;
  GtUword orig_endpos;
  GtUniqueEncseqEntrytype type;
  union {
    GtUniqueEncseqUniqueEntry *unique;
    GtUniqueEncseqLinkEntry *link;
  } entry;
} GtUniqueEncseqDBentry;

typedef struct {
  GtUword nelems;
  GtUword maxelems;
  GtUword cumulength;
  GtUword nseq;
  GtUword *ssp;
  GtUword *sde;
  char *desc;
  GtUniqueEncseqDBentry *fragmentdb;
} GtUniqueEncseqDB;

typedef struct
{
  const GtEncseq * encseq;
  unsigned int kmersize;
  GtUword nSequences;
  GtUword seqnum;
  GtUword seqlen;
  GtUword seqstartpos;
  GtUword totallength;
  unsigned int numofchars;
  const GtUchar * characters;
  GtUword maxlen;
  GtHashmap * hashmap;
  GtLogger * debug_logger;
  GtUword kmercount;
  GtUword kmerhitcount;
  GtUword minkmerhits;
  GtUword maxkmerhits;
  GtUword alignmentcount;
  GtUniqueEncseqArguments * arguments;
  GtUword currentposition;
  GtXdropArbitraryscores * arbitscores;
  GtXdropscore xdropbelowscore;
  char *buffer;
  GtSeqabstract *useq;
  GtSeqabstract *vseq;
  GtXdropresources *res;
  GtUword nextUniqueStartpos;
  GtUword initUniqueDBsize;
  GtUword initKmerCounter;
  GtUniqueEncseqDB *uniqueencseqdb;
  GtKmercodeiterator *kmercodeitMain;
  GtKmercodeiterator *kmercodeitAdd;
  const GtKmercode *kmercodeMain;
  const GtKmercode *kmercodeExt;
  const GtKmercode *kmercodeAdd;
  GtAlphabet *alphabet;
  GtUword unique_cumulen;
  GtUword nPosSinceInsert;
} GtUniqueEncseqInfo;

#endif
