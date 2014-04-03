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

#ifndef INTSET_IMPL_H
#define INTSET_IMPL_H

#define BITS_FOR_SIZE(SIZE)     ((SIZE) * CHAR_BIT)
#define ELEM2SECTION(LOGVAL,X)  ((X) >> (LOGVAL))
#define SECTIONMINELEM(S)       ((S) << intset->logsectionsize)

#define GT_IS_PASTER(A,B) A ## B
#define GT_IS_EVALUATOR(A,B) GT_IS_PASTER(A,B)

#define GT_INTSET_APPEND_S(STR) GT_IS_EVALUATOR(STR,GT_INTSET_SIZE)

#define GT_IS_ELEM_T_H(PRE,MID,SUF) PRE##MID##SUF
#define GT_IS_ELEM_T_(PRE,MID,SUF) GT_IS_ELEM_T_H(PRE, MID, SUF)

#define GT_IS_ELEM_T(PRE,SUF) \
  GT_IS_EVALUATOR(GT_IS_EVALUATOR(PRE, GT_INTSET_SIZE),SUF)

#define GT_IS_E_TYPE GT_IS_ELEM_T(uint,_t)
#define GT_IS_TYPE GT_INTSET_APPEND_S(GtIntset)

#define GT_INTSET_SIZE 8
#include "extended/intset.gen"
#undef GT_INTSET_SIZE

#define GT_INTSET_SIZE 16
#include "extended/intset.gen"
#undef GT_INTSET_SIZE

#define GT_INTSET_SIZE 32
#include "extended/intset.gen"
#undef GT_INTSET_SIZE

#endif
