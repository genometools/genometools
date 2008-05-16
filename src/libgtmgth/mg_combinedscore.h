/*
  Copyright (c) 2007 David Schmitz-Huebsch <dschmitz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef MG_COMBINEDSCORE_H
#define MG_COMBINEDSCORE_H

#include "libgtcore/codon.h"
#include "metagenomethreader.h"
#include "mg_computepath.h"

#define LONG_VALUE(VALUE, INDEX)\
                *(long*)array_get((VALUE), (INDEX))

#define POSITION(QUERY_FROM, HIT_NUMBER, POSITION, K)\
           LONG_VALUE(QUERY_FROM, HIT_NUMBER)+(POSITION)+(K)-1

/* Funktion zum Eintragen des Scores in die Combined-Score Matrix
   Parameter: Combined-Score-Matrix, Hit-AS, Query-AS, aktueller
              Leserahmen der Query-DNA, Position in der Query-Seq.,
              Position in der Hit-Seq., Laenge des Blast-Hits, Laenge der
              Contig-Seq, Hit-Nr, Zeiger auf die ParseStruct-Struktur,
              Hilfezeile fuer die Combined-Scores, Hilfszeile der Counts,
              Query-DNA-Seq, Hit-DNA-Seq, Zeiger auf die
              HitInformation-Struktur
   Returnwert: void */
static void fill_matrix(CombinedScoreMatrixEntry **,
                        char *,
                        char *,
                        short,
                        unsigned long,
                        unsigned long,
                        unsigned long,
                        unsigned long,
                        unsigned long,
                        ParseStruct *,
                        double *,
                        unsigned long *,
                        char *,
                        char *,
                        HitInformation *);

/* Funktion zur Berechnung des Matrixscore an der entsprechenden Position
   der Combined-Score-Matrix
   Parameter: Zeiger auf ParseStruct, Hilfszeilen matrix_row u. count_row,
              aktuelle Matrix-Zeile, Hit-Number, Position in der Query-DNA,
              Position im Triplet, Score - abh. vom Fall (syn, nonsyn,
              stop-codon in Query- oder Hit-Sequence, Blast-Hit-Ende)
   Returnwert: void */
static void add_scores(ParseStruct *,
                       double *,
                       unsigned long *,
                       short,
                       unsigned long,
                       unsigned long,
                       unsigned short,
                       double);

#endif
