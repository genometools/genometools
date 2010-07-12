/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2000-2004 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef GTHALIGNMENT_H
#define GTHALIGNMENT_H

#include <stdio.h>
#include "core/arraydef.h"
#include "core/file.h"
#include "core/score_matrix.h"
#include "gth/editoperation.h"

#define CONCRETEINTRONSYMBOL    '.'
#define CONCRETEGAPSYMBOL       '-'
#define ABSTRACTINTRONSYMBOL    (SEPARATOR-4) /* the abstract intron symbol */
#define ABSTRACTGAPSYMBOL       (SEPARATOR-3) /* the abstract gap symbol */

/*
  For the display of long introns in short form we need the following
  structure.
*/

typedef struct
{
  unsigned long start,
                length;
} ShortIntronInfo;

GT_DECLAREARRAYSTRUCT(ShortIntronInfo);

unsigned long gthfillthetwoalignmentlines(GtUchar **firstline,
                                          GtUchar **secondline,
                                          const GtUchar *useq,
                                          unsigned long ulen,
                                          const GtUchar *vseq,
                                          unsigned long vlen,
                                          Editoperation *alignment,
                                          unsigned long lenalg,
                                          unsigned long linewidth,
                                          unsigned long showintronmaxlen,
                                          GtArrayShortIntronInfo
                                          *shortintroninfo,
                                          unsigned long indelcount);

unsigned long gthfillthethreealignmentlines(GtUchar **firstline,
                                            GtUchar **secondline,
                                            GtUchar **thirdline,
                                            Editoperation *alignment,
                                            unsigned long lenalg,
                                            unsigned long indelcount,
                                            const GtUchar *genseqorig,
                                            unsigned long genseqlen,
                                            const GtUchar *refseqorig,
                                            unsigned long refseqlen,
                                            unsigned long
                                            translationschemenumber);

void gthshowalignmentprotein(GtFile *outfp,
                             unsigned long linewidth,
                             Editoperation *alignment,
                             unsigned long lenalg,
                             unsigned long indelcount,
                             const GtUchar *genseqorig,
                             unsigned long genseqlen,
                             const GtUchar *refseqorig,
                             unsigned long refseqlen,
                             unsigned long startfirst,
                             unsigned long startsecond,
                             unsigned long totalulen,
                             unsigned long showintronmaxlen,
                             GtAlphabet *alphabet,
                             unsigned long translationschemenumber,
                             GtScoreMatrix *score_matrix,
                             GtAlphabet *score_matrix_alphabet,
                             bool reverse_subject_pos,
                             bool wildcardimplosion);

void gthshowalignmentdna(GtFile *outfp,
                         unsigned long linewidth,
                         Editoperation *alignment,
                         unsigned long lenalg,
                         unsigned long indelcount,
                         const GtUchar *useqorig,
                         unsigned long ulen,
                         const GtUchar *vseqorig,
                         unsigned long vlen,
                         unsigned long startfirst,
                         unsigned long startsecond,
                         unsigned long totalulen,
                         unsigned long showintronmaxlen,
                         GtAlphabet *alphabet,
                         bool reverse_subject_pos,
                         bool wildcardimplosion);

#endif
