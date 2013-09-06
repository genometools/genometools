/*
  Copyright (c) 2003-2009 Gordon Gremme <gordon@gremme.org>
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
  GtUword start,
                length;
} ShortIntronInfo;

GT_DECLAREARRAYSTRUCT(ShortIntronInfo);

GtUword gthfillthetwoalignmentlines(GtUchar **firstline,
                                          GtUchar **secondline,
                                          const GtUchar *useq,
                                          GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vlen,
                                          Editoperation *alignment,
                                          GtUword lenalg,
                                          GtUword linewidth,
                                          GtUword showintronmaxlen,
                                          GtArrayShortIntronInfo
                                          *shortintroninfo,
                                          GtUword indelcount);

GtUword gthfillthethreealignmentlines(GtUchar **firstline,
                                            GtUchar **secondline,
                                            GtUchar **thirdline,
                                            Editoperation *alignment,
                                            GtUword lenalg,
                                            GtUword indelcount,
                                            const GtUchar *genseqorig,
                                            GtUword genseqlen,
                                            const GtUchar *refseqorig,
                                            GtUword refseqlen,
                                            GtUword
                                            translationschemenumber);

void gthshowalignmentprotein(GtFile *outfp,
                             GtUword linewidth,
                             Editoperation *alignment,
                             GtUword lenalg,
                             GtUword indelcount,
                             const GtUchar *genseqorig,
                             GtUword genseqlen,
                             const GtUchar *refseqorig,
                             GtUword refseqlen,
                             GtUword startfirst,
                             GtUword startsecond,
                             GtUword totalulen,
                             GtUword showintronmaxlen,
                             GtAlphabet *alphabet,
                             GtUword translationschemenumber,
                             GtScoreMatrix *score_matrix,
                             GtAlphabet *score_matrix_alphabet,
                             bool reverse_subject_pos,
                             bool wildcardimplosion);

void gthshowalignmentdna(GtFile *outfp,
                         GtUword linewidth,
                         Editoperation *alignment,
                         GtUword lenalg,
                         GtUword indelcount,
                         const GtUchar *useqorig,
                         GtUword ulen,
                         const GtUchar *vseqorig,
                         GtUword vlen,
                         GtUword startfirst,
                         GtUword startsecond,
                         GtUword totalulen,
                         GtUword showintronmaxlen,
                         GtAlphabet *alphabet,
                         bool reverse_subject_pos,
                         bool wildcardimplosion);

#endif
