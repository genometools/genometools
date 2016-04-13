/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#ifndef CONDENSEQ_REP_H
#define CONDENSEQ_REP_H

#include "core/alphabet_api.h"
#include "core/arraydef.h"
#include "core/encseq_api.h"
#include "core/types_api.h"
#include "extended/editscript.h"
#include "extended/intset.h"

#include "extended/condenseq.h"

/* this switch is meant to test performance differences and memory consumption
   issues */
#ifdef GT_CONDENSEQ_64_BIT
typedef uint64_t ces_unsigned;
#define CES_UNSIGNED_MAX ((ces_unsigned) UINT64_MAX)
#else
typedef uint32_t ces_unsigned;
#define CES_UNSIGNED_MAX ((ces_unsigned) UINT32_MAX)
#endif

/*
  The contents of this file is to be considered private implementation detail.
*/

/* TODO DW: maybe uint32_t for len, unique_id and offset would be all right?
   */
typedef struct {
  GtEditscript *editscript;
  GtUword       orig_startpos;
  ces_unsigned  len,
                unique_id,
                unique_offset;
} GtCondenseqLink;

typedef struct {
  GtArrayuint32_t links;
  GtUword         orig_startpos;
  ces_unsigned    len;
} GtCondenseqUnique;

struct GtCondenseq {
  GtAlphabet        *alphabet;
  GtCondenseqLink   *links;
  GtCondenseqUnique *uniques;
  GtEncseq          *unique_es;
  GtIntset          *sdstab,
                    *ssptab;
  GtUchar           *ubuffer;
  char              *buffer,
                    *filename,
                    *orig_ids;
  GtUword buffsize,
          id_len, /* GT_UNDEF_UWORD if sdstab != NULL */
          ids_total_len,
          ldb_allocated,
          ldb_nelems,
          orig_len,
          orig_num_seq,
          ubuffsize,
          udb_allocated,
          udb_nelems;
};

/* Returns a new GtCondenseq object, which is empty and can be filled by
   GtCondenseqCreator */
GtCondenseq* gt_condenseq_new(const GtEncseq *orig_es, GtLogger *logger);

/* Returns index of the unique element with the biggest orig_startpos smaller
   than <position>. if smallest is larger: return first. */
GtUword      gt_condenseq_uniques_position_binsearch(
                                                   const GtCondenseq *condenseq,
                                                   GtUword position);

/* Add <link> to <condenseq>, fails if links are not added in a sorted (by
   position) manner. */
void         gt_condenseq_add_link_to_db(GtCondenseq *condenseq,
                                         GtCondenseqLink link);

/* Add unique substring to <condenseq>. */
void         gt_condenseq_add_unique_to_db(GtCondenseq *condenseq,
                                           GtUword orig_startpos,
                                           ces_unsigned len);
#endif
