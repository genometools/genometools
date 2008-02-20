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

#include "libgtmgth/metagenomethreader.h"

/* Funktion, die nacheinander die erforderlichen Funktionen zur Erstellung
   des txt-Dokuments aufruft
   Parameter:  Zeiger auf ParseStruct-Struktur, CombinedScore-Matrix,
               die HitInformation-Struktur, die RegionStruct-Struktur,
               das Char-Zeichen, um welchen Bereich es sich handelt
   Returnwert: void */
static void outputwriter_txt(ParseStruct *,
                             CombinedScoreMatrixEntry **,
                             HitInformation *,
                             RegionStruct **, char, Error *);

/* Funktion, die nacheinander die erforderlichen Funktionen zur Erstellung
   des html-Dokuments aufruft
   Parameter: Zeiger auf ParseStruct-Struktur, CombinedScore-Matrix,
              die HitInformation-Struktur, die RegionStruct-Struktur,
              das Char-Zeichen, um welchen Bereich es sich handelt
   Returnwert: void */
static void outputwriter_html(ParseStruct *,
                              CombinedScoreMatrixEntry **,
                              HitInformation *,
                              RegionStruct **, char, Error *);

/* Funktion, die nacheinander die erforderlichen Funktionen zur Erstellung
   des xml-Dokuments aufruft
   Parameter: Zeiger auf ParseStruct-Struktur, CombinedScore-Matrix,
              die HitInformation-Struktur, die RegionStruct-Struktur,
              das Char-Zeichen, um welchen Bereich es sich handelt
   Returnwert: void */
static void outputwriter_xml(ParseStruct *,
                             CombinedScoreMatrixEntry **,
                             HitInformation *,
                             RegionStruct **, char, Error *);

/* Funktion zur Berechnung der AS-Sequenz
   Parameter: Zeiger auf die DNA-Sequnez, Str fuer die Proteinsequenz,
              from-Wert, to-Wert, aktueller Frame
   Returnwert: void */
static void as_coding(ParseStruct *,
                      char *,
                      Str *,
                      unsigned long,
                      unsigned long, unsigned short, Error *);

/* Funktion zum Schreiben des Output-Header der txt-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_header_txt(ParseStruct *);

/* Funktion zum Schreiben des Output-Header der html-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_header_html(ParseStruct *);

/* Funktion zum Schreiben des Output-Header der xml-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_header_xml(ParseStruct *);

/* Funktion zum Schreiben des Query-DNA Abschnittes der txt-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_querydna_txt(ParseStruct *);

/* Funktion zum Schreiben des Query-DNA Abschnittes der html-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_querydna_html(ParseStruct *);

/* Funktion zum Schreiben des Query-DNA Abschnittes der xml-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_querydna_xml(ParseStruct *);

/* Funktion zur Berechnung der Ausgabemenge fuer den Hit-Information
   Abschnittes;
   Parameter: Zeiger auf ParseStruct-Struktur, CombinedScore-Matrix,
              die HitInformation-Struktur, die RegionStruct-Struktur,
              das Char-Zeichen, um welchen Bereich es sich handelt
   Returnwert: void */
void output_hitdna(ParseStruct *,
                   CombinedScoreMatrixEntry **,
                   HitInformation *, RegionStruct **, Error *);

/* Funktion zum Schreiben des Coding-DNA-Abschnittes
   Parameter:  Zeiger auf ParseStruct, Zeiger auf die DNA-Seq,
               Str der Protein-Seq, Ausgabe-Format
   Returnwert: void */
static void print_codingheader(ParseStruct *, char *, Str *);

/* Funktion zum Schreiben des Hit-Information-Abschnittes
   Parameter: Zeiger auf die ParseStruct, Zeiger auf die
              HitInformation-Struktur, Zeiger auf die kodierende-Seq,
              Ausgabeformat, Sequenzposition
   Returnwert: void */
static void print_hitinformation(ParseStruct *,
                                 HitInformation *, unsigned long);

/* Funktion zum Schreiben des Output-Footer der HTML-Datei
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_footer_html(ParseStruct *);

/* Funktion zum Schreiben des Output-Footer der XML-Datei
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_footer_xml(ParseStruct *);

/* Funktion zum Schreiben der Metagenomethreader-Statistik der Text-Datei
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_statistics_txt(ParseStruct *);

/* Funktion zum Schreiben der Metagenomethreader-Statistik der HTML-Datei
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_statistics_html(ParseStruct *);

/* Funktion zum Schreiben der Metagenomethreader-Statistik der XML-Datei
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_statistics_xml(ParseStruct *);

/* Funktion zum uebersetzen eines Codons in eine Aminosaeure
   Parameter: die drei Zeichen des Codons
   Returnwert: die AS, die durch das uebergebene Codon kodiert wurde */
char mg_codon2amino(char, char, char);

/* Funktion zum Schreiben des Statistic-Headers
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_statistics_header(ParseStruct *);

/* Funktion zum Abschluss eines Iterations-Bereich im XML-File
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_close_iteration_xml(ParseStruct *);

#endif
