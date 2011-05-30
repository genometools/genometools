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

#ifndef METAGENOMETHREADER_H
#define METAGENOMETHREADER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "core/array.h"
#include "core/array2dim_api.h"
#include "core/bioseq.h"
#include "core/dlist_api.h"
#include "core/error.h"
#include "core/hashmap-generic.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/versionfunc.h"

/* jeweils die Anzahl der zu betrachtenden XML-Tags; Definitionen werden auch
 * in Schleifenkoepfen verwendet */
#define MGTH_QUERY_SIZE           1
#define MGTH_HIT_SIZE             3
#define MGTH_HIT_TO_QUERY_SIZE    9

#define MGTH_SET             1
#define MGTH_UNSET           0
#define MGTH_MEMORY_SIZE     250

/* Makros zum Zugriff auf die Parsestruct-Strukturen und Substrukturen */
#define ARGUMENTS(PATH)\
             parsestruct.metagenomethreader_arguments.PATH
#define ARGUMENTSSTRUCT(PATH)\
             parsestruct_ptr->metagenomethreader_arguments.PATH
#define MATRIXSTRUCT(PATH)      parsestruct_ptr->matrix_info.PATH
#define XMLPARSERSTRUCT(PATH)   parsestruct_ptr->xmlparser_static.PATH
#define HITSTRUCT(PATH)         parsestruct_ptr->hits_statistics.PATH
#define FILEPOINTEROUT          parsestruct_ptr->fp_outputfile
#define HITFILEOUT              parsestruct_ptr->fp_giexp_file
#define PARSESTRUCT(PATH)       parsestruct_ptr->PATH

/* Strukturdefinition fuer die Eintraege in der Combined-Score-Matrix */
typedef struct
{
  double matrix_score;
  unsigned long count;
  GtArray *hit_number;
} CombinedScoreMatrixEntry;

/* Struktur zur Speicherung der Bereichsgrenzen kodierender Abschnitte */
typedef struct
{
  GtArray *from;
  GtArray *to;
} RegionStruct;

/* Strukturdefinition fuer die Struktur zur Speicherung ausschliesslich
   der Hit-Informationen, die zu einem syn-Score beitragen */
typedef struct
{
  GtStrArray *hit_gi,
   *hit_def,
   *hit_hsp_nr,
   *hit_from,
   *hit_to;
} HitInformation;

/* Strukturdefinition fuer die Eintraege in die Path-Matrix */
typedef struct
{
  double score;
  unsigned short path_frame;
} PathMatrixEntry;

/* Strukturdefinition, die zur Aufnahme der Kommandozeilenargumente
   dient */
typedef struct
{
  double synonomic_value,
    nonsynonomic_value,
    blasthit_end_value,
    stopcodon_queryseq,
    stopcodon_hitseq,
    frameshift_span,
    prediction_span,
    leavegene_value,
    percent_value;
  GtStr *curl_fcgi_db,
   *outputtextfile_name,
   *giexpfile_name;
  int outputfile_format,
    codon_mode;
  bool hitfile_bool,
    homology_mode,
    extended_mode,
    testmodus_mode;
  unsigned long min_as;
} MetagenomeThreaderArguments;

/* Typdefinition der Struktur zur Aufnahme von Informationen aus dem
   XML-File und den DNA-Informationen von Query und Hit */
typedef struct
{
  GtArray *query_frame,
   *hit_frame,
   *query_from,
   *query_to;
  GtStr *query_dna,
   *query_def;
  GtStrArray *hit_gi_nr,
   *hit_num,
   *hit_gi_def,
   *hit_acc,
   *fasta_row,
   *hit_from,
   *hit_to,
   *hit_dna,
   *hsp_qseq,
   *hsp_hseq;
} MatrixMemory;

/* Struktur fuer den Statistikbereich */
typedef struct
{
  GtStrArray *hits_statistic;
  unsigned long *hitsnum,
   *memory,
    hitsnumber,
    stat_pos;
} HitsStatistic;

/* Struktur fuer Variablen in der mg_xmlparser-Datei */
typedef struct
{
  unsigned short query_array_index_start,
    hit_array_index_start,
    hit_hsp_array_index_start,
    query_array_index_end,
    hit_array_index_end,
    hit_hsp_array_index_end;
  unsigned long hit_counter;
} XMLparser_static;

/* Struktur fuer Variablen in der mg_compute_gene_prediction-Datei */
typedef struct
{
  double matrixscore,
    matrixscore_before;
  short current_frame,
    frame_before,
    function_stop;
  unsigned long noncodingcounter,
    codingcounter;
} GenePrediction_static;

/* Typdefinition der Struktur, die als Data an die Expat-Funktionen
   uebergeben werden kann Hauptstruktur, die auch an die wesentlichen
   anderen Funktionen uebergeben wird */

typedef struct
{
  GtStrArray *query_array,
   *hit_array,
   *hit_hsp_array,
   *query_frame_tmp,
   *hit_frame_tmp,
   *key_tmp;
  GtStr *xml_tag,
   *buf_ptr,
   *hit_fastafile,
   *xmlfile,
   *gi_def_tmp,
   *gi_acc_tmp,
   *hit_gi_nr_tmp,
   *fasta_row,
   *result_hits;
  GtArray *value_tmp;
  GtBioseq *queryseq,
   *hitseq;
  GtFile *fp_outputfile,
   *fp_blasthit_file,
   *fp_giexp_file;
  GtHashtable *queryhash,
   *hithash,
   *resulthits;
  GtDlist *outlist;
  GtError *err;
  int had_err;
  unsigned short def_flag,
    hit_flag,
    hit_hsp_flag,
    xml_tag_flag,
    giexp_flag,
    gi_flag;
  unsigned long xml_linenumber,
    hits_memory;
  double syn,
    non_syn;
  MatrixMemory matrix_info;
  MetagenomeThreaderArguments metagenomethreader_arguments;
  HitsStatistic hits_statistics;
  XMLparser_static xmlparser_static;
  GenePrediction_static geneprediction_static;
} ParseStruct;

/* specific access mode of queryhash */
DECLARE_HASHMAP(char *, gt_cstr_nofree, unsigned long *, ulp,,)
DECLARE_SAFE_DEREF(unsigned long *, ulp)

/* Funktion, mit der der Metagenomethreader gestartet wird
   Parameter: Anzahl der Kommandozeilenargumente, Zeiger auf die
              Kommandozeilenargumente, GtError-Variable
   Returnwert: Fehlercode */
int gt_metagenomethreader(int argc, const char **argv, GtError *);

/* Funktion zur Ueberpruefung auf ein Stop-Codon
   Parameter:  Zeiger auf ein Triplet von Zeichen
   Returnwert: 0 = kein Stop-Codon, 1 = Stop-Codon */
short gt_check_stopcodon(char *);

/* Funktion zum Auslesen der Query-DNA Sequenz zu einem GI-Def Eintrag des
   XML-Files
   Parameter:  Zeiger auf StringArray (Sequenzen der Treffer zu einer
               Hit-GI-Nr), String-Zeiger auf die Hit-GI-Nr, String-Zeiger
               auf Hit-From und String-Zeiger auf den Hit-To String,
               GtError-Variable
   Returnwert: void */
int mg_curl(ParseStruct *, unsigned long, GtError *);

/* Funktion zur Berechnung der Combined-Scores
   Parameter: Zeiger auf die Parsestruct; Anzahl der Hits zur betrachteten
              Hit-GI-Nr; GtError-Variable
   Returnwert: void */
int mg_combinedscore(ParseStruct *, unsigned long, GtError *);

/* Funktion zur Ausgabe der berechneten Ergebnisse
   Parameter:  Zeiger auf ParseStruct-Struktur, CombinedScore-Matrix,
               die HitInformation-Struktur, die RegionStruct-Struktur,
               das Char-Zeichen, um welchen Bereich es sich handelt
   Returnwert: void */
void mg_outputwriter(ParseStruct *,
                     CombinedScoreMatrixEntry **,
                     HitInformation *, RegionStruct **, char, GtError *);

/* Funktion zur Bestimmung der zu einem Lesereahmen gehoerigen Matrix-Zeile
   Parameter: aktueller Leserahmen
   Returnwert: Matrixzeile */
short get_matrix_row(long);

/* Funktion zur Bestimmung des Leserahmens zu einer Matrix-Zeile
   Parameter: Matrixzeile
   Returnwert: Leserahmen */
short get_current_frame(long);

/* Funktion zur Berechnung des reversen Komplements
   Parameter: Zeiger auf eine Seq., Seq-Laenge, GtError-Variable
   Returnwert: had_err */
int mg_reverse_complement(char *, unsigned long, GtError *);

#endif
