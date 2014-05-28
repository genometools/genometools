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
#ifndef INTSET_REP_H
#define INTSET_REP_H

#include <stdbool.h>
#include <stdlib.h>

#include "core/error_api.h"
#include "core/types_api.h"
#include "extended/intset.h"

#define GT_BITS_FOR_TYPE(TYPE)     ((sizeof (TYPE)) * ((size_t) CHAR_BIT))
#define GT_ELEM2SECTION(X, LOGVAL) ((X) >> (LOGVAL))
#define GT_SECTIONMINELEM(S)       ((S) << members->logsectionsize)

typedef struct GtIntsetClass GtIntsetClass;
typedef struct GtIntsetMembers GtIntsetMembers;

typedef void      (*GtIntsetAddFunc)       (GtIntset*, GtUword);
typedef bool      (*GtIntsetFileIsTypeFunc)(GtUword);
typedef GtUword   (*GtIntsetGetFunc)       (GtIntset*, GtUword);
typedef GtIntset* (*GtIntsetIOFunc)        (GtIntset*, FILE*, GtError*);
typedef GtUword   (*GtIntsetIdxSmGeqFunc)  (GtIntset*, GtUword);
typedef bool      (*GtIntsetIsMemberFunc)  (GtIntset*, GtUword);
typedef size_t    (*GtIntsetRepSizeFunc)   (GtUword, GtUword);
typedef GtUword   (*GtIntsetSizeFunc)      (GtIntset*);
typedef size_t    (*GtIntsetStructSizeFunc)(void);
typedef GtIntset* (*GtIntsetWriteFunc)     (GtIntset*, FILE*, GtError*);
typedef void      (*GtIntsetDeleteFunc)    (GtIntset*);

struct GtIntset {
  const GtIntsetClass *c_class;
  GtIntsetMembers *members;
};

struct GtIntsetClass {
  size_t                 size;
  GtIntsetAddFunc        add_func;
  GtIntsetFileIsTypeFunc file_is_type_func;
  GtIntsetGetFunc        get_func;
  GtIntsetIOFunc         io_func;
  GtIntsetIdxSmGeqFunc   idx_sm_geq_func;
  GtIntsetIsMemberFunc   is_member_func;
  GtIntsetRepSizeFunc    rep_size_func;
  GtIntsetSizeFunc       size_func;
  GtIntsetStructSizeFunc struct_size_func;
  GtIntsetWriteFunc      write_func;
  GtIntsetDeleteFunc     delete_func;
};

struct GtIntsetMembers {
  GtUword *sectionstart;
  size_t   logsectionsize;
  GtUword  currentsectionnum,
           maxelement,
           nextfree,
           num_of_elems,
           numofsections,
           previouselem,
           refcount;
};

const GtIntsetClass* gt_intset_class_new(size_t size,
                                         GtIntsetAddFunc,
                                         GtIntsetFileIsTypeFunc,
                                         GtIntsetGetFunc,
                                         GtIntsetIOFunc,
                                         GtIntsetIdxSmGeqFunc,
                                         GtIntsetIsMemberFunc,
                                         GtIntsetRepSizeFunc,
                                         GtIntsetSizeFunc,
                                         GtIntsetStructSizeFunc,
                                         GtIntsetWriteFunc,
                                         GtIntsetDeleteFunc);

GtIntset*            gt_intset_create(const GtIntsetClass*);

void*                gt_intset_cast(const GtIntsetClass*, GtIntset*);

/* Function for unit tests within implementations of this class. Fails if
   <gt_intset_is_member()> called with any number between and including <start>
   and <end> returns true.
   */
int gt_intset_unit_test_notinset(GtIntset *intset, GtUword start,
                                 GtUword end, GtError *err);

/* Function for unit tests within implementations of this class. Fails if
   <gt_intset_get_idx_smaller_geq()> called with any number between and
   including <start> and <end> returns any number different than <num>. */
int gt_intset_unit_test_check_seqnum(GtIntset *intset, GtUword start,
                                     GtUword end, GtUword num, GtError *err);

#endif
