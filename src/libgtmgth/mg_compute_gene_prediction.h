
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

#include <limits.h>
#include <float.h>
#include "libgtmgth/metagenomethreader.h"

#define GENEPREDSTRUCT(path)    parsestruct_ptr->geneprediction_static.path

/* Funktion zur Berechnung der Start-Stop Informationen der genkodierenden
   Bereiche sowie deren Leserahmen
   Parameter: row, column, max-Wert der letzten Spalte der Opt-Path-Matrix,
              Path-Matrix, CombinedScore-Matrix, Array des opt. Pfades
              (Leserahmen), Hit-Information, Zeiger auf ParseStruct, Zeiger
              auf RegionStruct, Zeiger auf den Speicherbereich zum
              Abspeichern der beteiligten Hits an einem Ergebnis
   Returnwert: void */
static void gene_prediction(unsigned short,
                            unsigned long,
                            double,
                            PathMatrixEntry **,
                            CombinedScoreMatrixEntry **,
                            Array *,
                            HitInformation *,
                            ParseStruct *,
                            RegionStruct **, unsigned long *, Error *);

/* Funktion zur Identifizierung von Frameshifts
   Parameter: Zeiger auf ParseStruct, Zeiger auf RegionStruct,
              der "reale" Leserahmen
   Returnwert: void */
int frameshiftprocessing(ParseStruct *, RegionStruct **, short, Error *);

/* Funktion zur Vereinigung von genkodierenden Bereichen innerhalb des
   selben Leserahmen
   Parameter: Zeiger auf ParseStruct, Zeiger auf RegionStruct,
              der "reale" Leserahmen
   Returnwert: void */
static void genemergeprocessing(ParseStruct *, RegionStruct **, Error *);

/* Funktion zur Ueberpruefung von Sequenzbereichen auf beinhaltende
   Stop-Codons
   Parameter: Zeiger auf ParseStruct, from-Wert, to-Wert,
              aktueller Leserahmen
   Returnwert: 0 - kein Stop-Codon; 1 - Stop-Codon */
static int check_coding(ParseStruct *,
                        unsigned long, unsigned long, short, Error *);

/* Funktion zum sortierten Einfuegen neuer Bereichsgrenzen kodierender
   Abschnitte in das real-Frame-Array
   Parameter: Zeiger auf RegionStruct-Struktur, Arrays mit den From- und
              To-Werten des real-Frames,Arrays mit den From- und To-Werten
              der einzufuegenden Abschnitte, Index des Real-Frames,Index
              des tmp-Frames, real-Frame, Zeilenindex;
   Returnwert: void */
static void merge_array(RegionStruct **,
                        Array *,
                        Array *,
                        Array *,
                        Array *,
                        unsigned long,
                        unsigned long, short, unsigned short);

/* Funktion zum sortierten der Arrays mit den neu einzufuegenden
   Bereichsgrenzen
   Parameter: Arrays der sortierten From- und To-Werte, Arrays mit den
              From- und To-Werten der zu sortierenden Arrays;
   Returnwert: void */
static void sort_realtmp(Array *, Array *, Array *, Array *);
