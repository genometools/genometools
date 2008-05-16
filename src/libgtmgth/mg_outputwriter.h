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

#ifndef MG_OUTPUTWRITER_H
#define MG_OUTPUTWRITER_H

#include "metagenomethreader.h"
#include "mg_codon.h"

/* Funktion zur Ausgabe der berechneten Ergebnisse
   Parameter:  Zeiger auf ParseStruct-Struktur, CombinedScore-Matrix,
               die HitInformation-Struktur, die RegionStruct-Struktur,
               das Char-Zeichen, um welchen Bereich es sich handelt
   Returnwert: void */
void mg_outputwriter(ParseStruct *,
                     CombinedScoreMatrixEntry **,
                     HitInformation *, RegionStruct **, char, Error *);
#endif
