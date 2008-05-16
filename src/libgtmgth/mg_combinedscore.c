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

#include "mg_combinedscore.h"

int mg_combinedscore(ParseStruct *parsestruct_ptr,
                     unsigned long hit_counter, Error * err)
{
  int had_err = 0;

  unsigned short current_row,
    k;

  unsigned long contig_len = 0,
    contig_seq_diff = 0,
    hit_len = 0,
    hit_seq_diff = 0,
    i = 0,
    j = 0,
    contig_index = 0,
    hit_index = 0,
    gap_len = 0,
    ulhit_from = 0,
    ulhit_to = 0,
    mod = 0,
    *count_row;                        /* Hilfszaehler fuer die
                                          involvierten Hits fuer den jew.
                                          combined-Score */

  /* Matrix der Combined-Scores */
  CombinedScoreMatrixEntry **combinedscore_matrix;

  /* Hit-Informationen der relevanten Blast-Hits (relevant:
     Synonym/Nicht-Synonym Rate >=1) */
  HitInformation hit_information;

  /* Hilfszeile zur Abarbeitungs eines Hit */
  double *matrix_row;

  char contig_triplet[3],              /* Query-Triplet von dem die AS
                                          bestimmt werden soll */
    contig_as,                         /* Aminosaeure des uebergebenen
                                          Query-Triplets */
    hit_triplet[3],                    /* Hit-Triplet von dem die AS
                                          bestimmt werden soll */
    hit_as,                            /* Aminosaeure des uebergebenen
                                          Hit-Triplets */
   *contig_seq = NULL,                 /* Contig-DNA-Sequenz */
    *hit_seq = NULL,                   /* Hit-DNA-Sequenz */
    *contig_seq_tri = NULL,
    *hit_seq_tri = NULL,
    *contig_seq_ptr;

  const char *contig_as_ptr,
   *hit_as_ptr;

  /* Check Umgebungsvariablen */
  error_check(err);

  /* Zeiger auf den vollstaendigen Query-DNA Eintrag */
  contig_seq_ptr = str_get(MATRIXSTRUCT(query_dna));

  /* contig_len legt die Anzahl der Spalten der Matrix fest */
  contig_len = str_length(MATRIXSTRUCT(query_dna));

  /* Speicherplatzreservierung fuer die Combined-Score Matrix */
  array2dim_calloc(combinedscore_matrix, 7, contig_len);

  /* Erstellen der Arrays fuer die Hit-Nummer fuer jeden Matrix-Eintrag
     Ueber die Hit-Nummer erfolgt die Zuordnung zu den entsprechenden
     Blast-Hits in der Hit-Information-Struktur */
  for (k = 0; k < 7; k++)
  {
    for (j = 0; j < contig_len; j++)
    {
      combinedscore_matrix[k][j].hit_number =
        array_new(sizeof (unsigned long));
    }
  }

  /* String-Arrays der Hit-Information Struktur anlegen */
  hit_information.hit_gi = strarray_new();
  hit_information.hit_def = strarray_new();
  hit_information.hit_hsp_nr = strarray_new();
  hit_information.hit_from = strarray_new();
  hit_information.hit_to = strarray_new();

  /* fuer jeden Hit erfolgt eine Berechnung der Combined-Scores */
  for (i = 0; i < hit_counter; i++)
  {
    /* Speicherplatz fuer die Hilfszeilen wird reserviert */
    matrix_row = ma_calloc(contig_len, sizeof (double));
    count_row = ma_calloc(contig_len, sizeof (unsigned long));
    contig_seq_tri = ma_calloc(4, sizeof (char));
    hit_seq_tri = ma_calloc(4, sizeof (char));

    /* Zeiger auf die Proteinsequenzen von Hit und Query */
    contig_as_ptr = strarray_get(MATRIXSTRUCT(hsp_qseq), i);
    hit_as_ptr = strarray_get(MATRIXSTRUCT(hsp_hseq), i);

    /* Funktion zur Bestimmung der Zeile in der Matrix, die fuer den
       aktuellen Leserahmen steht */
    current_row = get_matrix_row(LONG_VALUE(MATRIXSTRUCT(query_frame), i));

    /* to-from+1 ist die anzahl der zu betrachtenden DNA-Basen,
       zusaetzlich +1 fuer das Stringendezeichen */
    contig_seq_diff =
      LONG_VALUE(MATRIXSTRUCT(query_to), i)
                 - LONG_VALUE(MATRIXSTRUCT(query_from), i) + 2;

    hit_len = strlen(strarray_get(MATRIXSTRUCT(hit_dna), i));

    ulhit_from = atol(strarray_get(MATRIXSTRUCT(hit_from), i));
    ulhit_to = atol(strarray_get(MATRIXSTRUCT(hit_to), i));
    hit_seq_diff = ulhit_to - ulhit_from + 2;

    /* Ueberpruefen der Vereinbarkeit von Query- und Hit-Sequenz
       unterschiedliche Laengen sind nur in 3er Schritten erlaubt -> Gaps
       in der AS-Sequenz */
    if (contig_seq_diff - hit_seq_diff > 0)
    {
      mod = (contig_seq_diff - hit_seq_diff) % 3;
    }
    else
      mod = (hit_seq_diff - contig_seq_diff) % 3;

    if ((LONG_VALUE(MATRIXSTRUCT(query_to), i) > contig_len)
        || (hit_seq_diff - 1 != hit_len))
      mod = 1;

    /* Fehlermeldung bei Unvereinbarkeit */
    if (mod != 0)
    {
      error_set(err,
                "sequences error: matching sequences do not fit in length.\
                 wrong FASTA-files or please delete entry %s!?",
                strarray_get(MATRIXSTRUCT(hit_gi_def), i));
      had_err = -1;
    }

    /* Speicherplatreservierung fuer die Query-DNA-Seq. */
    contig_seq = ma_calloc(contig_seq_diff, sizeof (char));
    /* kopieren von contig_seq_diff-1 Zeichen */
    (void) snprintf(contig_seq, contig_seq_diff, "%s",
                    contig_seq_ptr +
                    (LONG_VALUE(MATRIXSTRUCT(query_from), i) - 1));

    /* die laenge der Hit-Sequenz kann max. der Laenge der QueryDNA Seq
       entsprechen */
    hit_seq = ma_calloc(hit_seq_diff, sizeof (char));
    /* kopieren von hit_seq_diff-1 Zeichen */
    (void) snprintf(hit_seq, hit_seq_diff, "%s",
                    strarray_get(MATRIXSTRUCT(hit_dna), i));

    /* Bei einem negativen Leserahmen muss das Reverse-Komplement der
       Sequenz gebildet werden */
    if ((LONG_VALUE(MATRIXSTRUCT(query_frame), i) < 0) && !had_err)
    {
      /* bestimmen des Reverse-Komplement */
      had_err =
        mg_reverse_complement(contig_seq, contig_seq_diff - 1, err);
    }

    /* Gleiches Vorgehen bei negativen Leserahmen wie bei der Query-DNA */
    if ((LONG_VALUE(MATRIXSTRUCT(hit_frame), i) < 0) && !had_err)
    {
      had_err = mg_reverse_complement(hit_seq, hit_seq_diff - 1, err);
    }

    if (!had_err)
    {
      for (j = 0, contig_index = 0, hit_index = 0; j < hit_len - 2;
           j += 3, contig_index += 3, hit_index += 3)
      {
        if (contig_index < contig_len && hit_index < hit_len)
        {
          contig_as = contig_as_ptr[j / 3];
          hit_as = hit_as_ptr[j / 3];

          if (contig_as == '-')
          {
            gap_len = strspn(contig_as_ptr + j / 3, "-");
            hit_index += 3 * gap_len;
          }
          if (hit_as == '-')
          {
            gap_len = strspn(hit_as_ptr + j / 3, "-");
            contig_index += 3 * gap_len;
          }
          if (hit_as != '-' && contig_as != '-')
          {
            (void) snprintf(contig_seq_tri, 4, "%s",
                            contig_seq + contig_index);
            (void) snprintf(hit_seq_tri, 4, "%s", hit_seq + hit_index);

            if ((strspn(contig_seq_tri, "acgtuACGTU") == 3)
                && (strspn(hit_seq_tri, "acgtuACGTU") == 3))
            {
              /* Aktuelles Hit-Triplet */
              hit_triplet[0] = hit_seq[hit_index];
              hit_triplet[1] = hit_seq[hit_index + 1];
              hit_triplet[2] = hit_seq[hit_index + 2];

              /* Aktuelles Query-Triplet */
              contig_triplet[0] = contig_seq[contig_index];
              contig_triplet[1] = contig_seq[contig_index + 1];
              contig_triplet[2] = contig_seq[contig_index + 2];

              /* Bestimmen der AS der jeweiligen Triplets */
              contig_as =
                codon2amino(contig_triplet[0], contig_triplet[1],
                            contig_triplet[2]);
              hit_as =
                codon2amino(hit_triplet[0], hit_triplet[1],
                            hit_triplet[2]);
            }
          }

          /* Aufruf der Funktion fill_matrix; berechnet die
             Combined-Scores uebergeben werden: die Combined-Score-Matrix
             Hit-AS des aktuellen Triplets Query-AS des aktuelen Triplets
             der aktuelle Leserahmen j - hier: die Position in der
             Query-Sequenz hit-len und contig-len i - hier: die aktuell
             betrachtete Hit-Sequenz Zeiger auf die parsestruct die
             Hilfs-Matrix-Zeile Contig- und Hit-Sequenz */
          fill_matrix(combinedscore_matrix,
                      &hit_as,
                      &contig_as,
                      current_row,
                      contig_index,
                      hit_index,
                      hit_len,
                      contig_len,
                      i,
                      parsestruct_ptr,
                      matrix_row,
                      count_row, contig_seq, hit_seq, &hit_information);
        }
      }
      ma_free(contig_seq);
      ma_free(hit_seq);
      ma_free(contig_seq_tri);
      ma_free(hit_seq_tri);
      ma_free(matrix_row);
      ma_free(count_row);
    }
    else
    {
      ma_free(contig_seq);
      ma_free(hit_seq);
      ma_free(contig_seq_tri);
      ma_free(hit_seq_tri);
      ma_free(matrix_row);
      ma_free(count_row);

      break;
    }
  }

  if (!had_err)
  {
    /* nach dem erstellen der Combined-Score Matrix muessen die einzelnen
       Positionen noch mit der Anzahl der beteiligten Hits normalisiert
       werden */
    for (j = 0; j < contig_len; j++)
    {
      for (k = 0; k < 7; k++)
      {
        if (combinedscore_matrix[k][j].count != 0)
        {
          combinedscore_matrix[k][j].
            matrix_score /= combinedscore_matrix[k][j].count;
        }
      }
    }

    /* Aufruf der Funktion computepath - DP-Ansatz uebergeben werden: die
       Combined-Score-Matrix Anzahl Zeilen Anzahl Spalten */
    had_err = mg_computepath(combinedscore_matrix,
                             &hit_information,
                             7, contig_len, parsestruct_ptr, err);
  }

  for (i = 0; i < 7; i++)
  {
    for (j = 0; j < contig_len; j++)
    {
      array_delete(combinedscore_matrix[i][j].hit_number);
    }
  }

  array2dim_delete(combinedscore_matrix);

  strarray_delete(hit_information.hit_gi);
  strarray_delete(hit_information.hit_def);
  strarray_delete(hit_information.hit_hsp_nr);
  strarray_delete(hit_information.hit_from);
  strarray_delete(hit_information.hit_to);

  return had_err;
}

/* Funktion zur Bestimmung der dem Leserahmen entsprechenden Matrix-Zeile */
short get_matrix_row(long frame_fct)
{
  switch (frame_fct)
  {
    case -3:
      return 6;
    case -2:
      return 5;
    case -1:
      return 4;
    case 0:
      return 3;
    case 1:
      return 2;
    case 2:
      return 1;
    case 3:
      return 0;
    default:
      return 0;
  }
}

/* Umkehrfunktion zu get_matrix_row - aus der Matrix-Zeile wird der Leserahmen
   bestimmt */
short get_current_frame(long row_fct)
{
  switch (row_fct)
  {
    case 6:
      return -3;
    case 5:
      return -2;
    case 4:
      return -1;
    case 3:
      return 0;
    case 2:
      return 1;
    case 1:
      return 2;
    case 0:
      return 3;
    default:
      return 0;
  }
}

static void fill_matrix(CombinedScoreMatrixEntry **combinedscore_matrix,
                        char *hit_amino,
                        char *query_amino,
                        short current_row_fct,
                        unsigned long position_contig,
                        unsigned long position_hit,
                        unsigned long hit_len,
                        unsigned long contig_len,
                        unsigned long hit_number,
                        ParseStruct *parsestruct_ptr,
                        double *matrix_row,
                        unsigned long *count_row_fct,
                        char *contig_seq,
                        char *hit_seq, HitInformation *hit_information)
{
  unsigned long j = 0,
    nr_of_strings = 0;
  unsigned short k = 0;

  unsigned long query_from;
  unsigned long query_to;

  query_from = LONG_VALUE(MATRIXSTRUCT(query_from), hit_number) - 1;
  query_to = LONG_VALUE(MATRIXSTRUCT(query_to), hit_number) - 1;

  /* das Ende von Blast-Hits innerhalb der Query-Sequenz wird mit -10
     bewertet */
  if (position_hit == hit_len - 3
      && *(long *) MATRIXSTRUCT(query_to) != contig_len && k == 3)
  {
    /* Wenn der Leserahmen negativ ist, muessen die Combined-Scores von
       Rechts nach Links in die Combined-Score Matrix eingetragen werden -
       also ausgehend vom Query-to Wert */
    if (current_row_fct > 3)
    {
      matrix_row[POSITION
                 (MATRIXSTRUCT(query_to), hit_number, -position_contig,
                  -k + 1)] += ARGUMENTSSTRUCT(blasthit_end_value);
    }
    else
    {
      matrix_row[POSITION
                 (MATRIXSTRUCT(query_from), hit_number, position_contig,
                  k - 1)] += ARGUMENTSSTRUCT(blasthit_end_value);
    }
  }
  /* Stop-Codon in der Query-Sequenz, kein Stop-Codon in der Hit-Sequenz */
  else if (*hit_amino != '*' && *query_amino == '*')
  {
    /* Berechnung der Combined-Scores immer fuer die 3 DNA-Basen des
       Triplets Stop-Codon in der Query-Sequnez wird mit -2.0 bestraft;
       Anzahl der beteiligten Hits am Combined-Score um 1 erhoeht */
    for (k = 0; k < 3; k++)
    {
      add_scores(parsestruct_ptr, matrix_row, count_row_fct,
                 current_row_fct, hit_number, position_contig, k,
                 ARGUMENTSSTRUCT(stopcodon_hitseq));
    }
  }
  /* Stop-Codon in der Hit-Sequnez, kein Stop-Codon in der Query-Sequenz */
  else if (*hit_amino == '*' || *query_amino == '*')
  {
    /* Berechnung der Combined-Scores immer fuer die 3 DNA-Basen des
       Triplets Stop-Codon in der Hit-Sequnez wird mit -5.0 bestraft;
       Anzahl der beteiligten Hits am Combined-Score um 1 erhoeht */
    for (k = 0; k < 3; k++)
    {
      add_scores(parsestruct_ptr, matrix_row, count_row_fct,
                 current_row_fct, hit_number, position_contig, k,
                 ARGUMENTSSTRUCT(stopcodon_queryseq));
    }
  }
  /* AS von Hit und Query sind identisch */
  else if (*hit_amino == *query_amino)
  {
    /* Berechnung der Combined-Scores immer fuer die 3 DNA-Basen des
       Triplets */
    for (k = 0; k < 3; k++)
    {
      if (!ARGUMENTSSTRUCT(homology_mode))
      {
        /* DNA-Basen stimmen nicht ueberein - synonymer Austausch; hit_seq
           ist Zeiger auf die erste DNA-Base der Sequenz position ist die
           aktuelle Position in der DNA-Sequenz und k die Position im
           aktuellen Triplet */
        if (tolower(hit_seq[position_hit + k]) !=
            tolower(contig_seq[position_contig + k]))
        {
          /* Bei Uebereinstimmung der AS und keiner Basen-Uebereinstimmung
             wird der Zaheler der Synonymen Basen-Austausche erhoeht, der
             aktuellen Position in der Combined-Score-Matrix der Wert 1
             hinzuaddiert und gleichzeitig die Anzahl der beteiligten Hits
             ebenfalls um 1 erhoeht */
          add_scores(parsestruct_ptr, matrix_row, count_row_fct,
                     current_row_fct, hit_number, position_contig, k,
                     ARGUMENTSSTRUCT(synonomic_value));
          parsestruct_ptr->syn++;
        }
      }
      else
      {
        /* Experimenteller Status - Suche nach Homologien statt
           Orieentierung an Synonymen-Basenaustauschen */
        if (tolower(hit_seq[position_hit + k]) ==
            tolower(contig_seq[position_contig + k]))
        {
          /* Bei Uebereinstimmung der AS und keiner Basen-Uebereinstimmung
             wird der Zaheler der Synonymen Basen-Austausche erhoeht, der
             aktuellen Position in der Combined-Score-Matrix der Wert 1
             hinzuaddiert und gleichzeitig die Anzahl der beteiligten Hits
             ebenfalls um 1 erhoeht */
          add_scores(parsestruct_ptr, matrix_row, count_row_fct,
                     current_row_fct, hit_number, position_contig, k,
                     ARGUMENTSSTRUCT(synonomic_value));
          parsestruct_ptr->syn++;
        }
      }
    }
  }
  /* AS stimmen nicht ueberein - nicht synonymer Base-Austausch */
  else if (*hit_amino != *query_amino)
  {
    /* Berechnung der Combined-Scores immer fuer die 3 DNA-Basen des
       Triplets */
    for (k = 0; k < 3; k++)
    {
      /* DNA-Basen stimmen nicht ueberein - nicht-synonymer Austausch;
         hit_seq ist Zeiger auf die erste DNA-Base der Sequenz position
         ist die aktuelle Position in der DNA-Sequenz und k die Position
         im aktuellen Triplet */
      if (tolower(hit_seq[position_hit + k]) !=
          tolower(contig_seq[position_contig + k]))
      {
        /* Bei Nicht-Uebereinstimmung der AS und keiner
           Basen-Uebereinstimmung wird der Zaheler der Nicht-Synonymen
           Basen-Austausche erhoeht, vom Wert der aktuellen Position in
           der Combined-Score-Matrix der Wert 1 subtrahiert und
           gleichzeitig die Anzahl der beteiligten Hits ebenfalls um 1
           erhoeht */
        add_scores(parsestruct_ptr, matrix_row, count_row_fct,
                   current_row_fct, hit_number, position_contig, k,
                   ARGUMENTSSTRUCT(nonsynonomic_value));
        parsestruct_ptr->non_syn++;
      }
    }
  }
  /* Sequenzende und damit Ende der Bearbeitung des aktuellen Hits
     erreicht */
  if (hit_len - 3 == position_hit)
  {
    /* Falls die Non-Syn Anzahl 0 ist wird sie auf 1 gesetzt, um eine
       Division durch 0 zu vermeiden */
    if (parsestruct_ptr->non_syn == 0.0)
      parsestruct_ptr->non_syn = 1.0;
    /* Falls die Rate von Syn zu Non-Syn Austauschen kleiner als 1 ist,
       werden die Ergebnisse nicht betrachtet, da mit hoher WS keine
       codierende Sequenz vorliegt */
    if (parsestruct_ptr->syn / parsestruct_ptr->non_syn < 1.0)
    {
      parsestruct_ptr->syn = 0.0;
      parsestruct_ptr->non_syn = 0.0;
    }
    /* codierende Sequenz liegt mit hoher WS vor, da syn/non_syn >=1 */
    else
    {
      /* Abspeichern der Hit-Informationen */
      strarray_add_cstr(hit_information->hit_gi,
                        strarray_get(MATRIXSTRUCT(hit_gi_nr), hit_number));
      strarray_add_cstr(hit_information->hit_def,
                        strarray_get(MATRIXSTRUCT(hit_gi_def),
                                     hit_number));
      strarray_add_cstr(hit_information->hit_hsp_nr,
                        strarray_get(MATRIXSTRUCT(hit_num), hit_number));
      strarray_add_cstr(hit_information->hit_from,
                        strarray_get(MATRIXSTRUCT(hit_from), hit_number));
      strarray_add_cstr(hit_information->hit_to,
                        strarray_get(MATRIXSTRUCT(hit_to), hit_number));

      nr_of_strings = strarray_size(hit_information->hit_def) - 1;

      /* Uebertragen der Ergebnisse des aktuellen Hits in die
         Combined-Score Matrix */
      for (j = query_from; j < query_to + 1; j++)
      {
        combinedscore_matrix[current_row_fct][j].matrix_score +=
          matrix_row[j];
        combinedscore_matrix[current_row_fct][j].count += count_row_fct[j];
        array_add(combinedscore_matrix[current_row_fct][j].hit_number,
                  nr_of_strings);
      }
      parsestruct_ptr->syn = 0.0;
      parsestruct_ptr->non_syn = 0.0;
    }
  }
}

static void add_scores(ParseStruct *parsestruct_ptr,
                       double *matrix_row,
                       unsigned long *count_row_fct,
                       short current_row_fct,
                       unsigned long hit_number,
                       unsigned long position,
                       unsigned short k, double score)
{
  /* Wenn der Query-Leserahmen negativ ist, muessen die Combined-Scores
     von Rechts nach Links in die Combined-Score Matrix eingetragen werden
     - also ausgehend vom Query-to Wert */
  if (current_row_fct > 3)
  {
    matrix_row[POSITION(MATRIXSTRUCT(query_to), hit_number, -position, -k)]
      += score;
    count_row_fct[POSITION
                  (MATRIXSTRUCT(query_to), hit_number, -position, -k)] +=
      1;
  }
  else
  {
    matrix_row[POSITION(MATRIXSTRUCT(query_from), hit_number, position, k)]
      += score;
    count_row_fct[POSITION
                  (MATRIXSTRUCT(query_from), hit_number, position, k)] +=
      1;
  }
}
