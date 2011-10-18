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
#include <ctype.h>

#include "mg_compute_gene_prediction.h"

#include "core/unused_api.h"

/* Funktion zur Berechnung der Start-Stop Informationen der genkodierenden
   Bereiche sowie deren Leserahmen
   Parameter: row, column, max-Wert der letzten Spalte der Opt-Path-Matrix,
              Path-Matrix, CombinedScore-Matrix, GtArray des opt. Pfades
              (Leserahmen), Hit-Information, Zeiger auf ParseStruct, Zeiger
              auf RegionStruct, Zeiger auf den Speicherbereich zum
              Abspeichern der beteiligten Hits an einem Ergebnis
   Returnwert: void */
static void gene_prediction(unsigned short,
                            unsigned long,
                            double,
                            PathMatrixEntry **,
                            CombinedScoreMatrixEntry **,
                            GtArray *,
                            HitInformation *,
                            ParseStruct *,
                            RegionStruct **, unsigned long *, GtError *);

/* Funktion zur Vereinigung von genkodierenden Bereichen innerhalb des
   selben Leserahmen
   Parameter: Zeiger auf ParseStruct, Zeiger auf RegionStruct,
              der "reale" Leserahmen
   Returnwert: void */
static void genemergeprocessing(ParseStruct *, RegionStruct **, GtError *);

/* Funktion zur Identifizierung von Frameshifts
   Parameter: Zeiger auf ParseStruct, Zeiger auf RegionStruct,
              der "reale" Leserahmen
   Returnwert: void */
static int frameshiftprocessing(ParseStruct *, RegionStruct **, short,
                                GtError *);

/* Funktion zur Ueberpruefung von Sequenzbereichen auf beinhaltende
   Stop-Codons
   Parameter: Zeiger auf ParseStruct, from-Wert, to-Wert,
              aktueller Leserahmen
   Returnwert: 0 - kein Stop-Codon; 1 - Stop-Codon */
static int check_coding(ParseStruct *,
                        unsigned long, unsigned long, short, GtError *);

/* Funktion zum sortierten Einfuegen neuer Bereichsgrenzen kodierender
   Abschnitte in das real-Frame-Array
   Parameter: Zeiger auf RegionStruct-Struktur, GtArrays mit den From- und
              To-Werten des real-Frames,Arrays mit den From- und To-Werten
              der einzufuegenden Abschnitte, Index des Real-Frames,Index
              des tmp-Frames, real-Frame, Zeilenindex;
   Returnwert: void */
static void merge_array(RegionStruct **,
                        GtArray *,
                        GtArray *,
                        GtArray *,
                        GtArray *,
                        unsigned long,
                        unsigned long, short, unsigned short);

/* Funktion zum sortierten der GtArrays mit den neu einzufuegenden
   Bereichsgrenzen
   Parameter: GtArrays der sortierten From- und To-Werte, GtArrays mit den
              From- und To-Werten der zu sortierenden GtArrays;
   Returnwert: void */
static void sort_realtmp(GtArray *, GtArray *, GtArray *, GtArray *);

int mg_compute_gene_prediction(CombinedScoreMatrixEntry
                               **combinedscore_matrix,
                               PathMatrixEntry **path_matrix,
                               unsigned long contig_len,
                               HitInformation *hit_information,
                               ParseStruct *parsestruct_ptr, GtError * err)
{
  int had_err = 0;

  /* Variablen zur Bestimmung des Maximums der letzten Path-Matrix-Spalte */
  double max_lastcolumn = DBL_MIN,
    max_tmp;

  /* Zaehlvariablen fuer die Zeilen */
  unsigned short row_index,
    row_idx,
    real_frame = 0;

  /* Zaehlvariable fuer die Haeufigkeiten der einzelnen Frames sowie
     Variable fuer den Frame-Score */
  unsigned long *frame_counter,
    real_frame_score;

  /* Struktur zur Speicherung der kodierenden Bereiche (from-to - Werte) */
  RegionStruct **regionmatrix;

  /* GtArray fuer das Speichern der Frames des optimalen Pfades */
  GtArray *frame_path_array;
  frame_path_array = gt_array_new(sizeof (unsigned short));

  gt_error_check(err);

  /* Bestimmen des/der Max-Werte(s) der letzten Spalte der path_matrix */
  for (row_index = 0; row_index < 7; row_index++)
  {
    max_tmp = path_matrix[row_index][contig_len - 1].score;

    if (max_tmp > max_lastcolumn)
    {
      max_lastcolumn = max_tmp;
    }
  }

  /* ausgehend von den Max-Werten startet die Vorhersage codierender
     Bereiche */
  for (row_index = 0; row_index < 7; row_index++)
  {
    /* Matrix-Score der aktuellen Zeile entspricht dem Max-Wert */
    if (path_matrix[row_index][contig_len - 1].score == max_lastcolumn)
    {
      gt_array2dim_calloc(regionmatrix, 7, 1);
      /* frame-counter zaehlt die Auftritts-haeufigkeiten der einzelnen
         Frames im opt. Pfad */
      frame_counter = gt_calloc(7, sizeof (unsigned long));
      real_frame_score = 0;

      for (row_idx = 0; row_idx < 7; row_idx++)
      {
        regionmatrix[row_idx][0].from = gt_array_new(sizeof (unsigned long));
        regionmatrix[row_idx][0].to = gt_array_new(sizeof (unsigned long));
      }
      /* Aufruf der Genvorhersagemethode */
      gene_prediction(row_index,
                      contig_len - 1,
                      max_lastcolumn,
                      path_matrix,
                      combinedscore_matrix,
                      frame_path_array,
                      hit_information,
                      parsestruct_ptr, regionmatrix, frame_counter, err);

      /* Falls es fuer Frame X kodierende Abschnitte gibt, werden diese
         von hinten nach vorne in der regionmatrix abgespeichert - hier
         werden sie nun umgedreht, um sie in aufsteigender Form vorliegen
         zu haben */
      for (row_idx = 0; row_idx < 7; row_idx++)
      {
        if (gt_array_size(regionmatrix[row_idx][0].from) > 0)
        {
          gt_array_reverse(regionmatrix[row_idx][0].from);
          gt_array_reverse(regionmatrix[row_idx][0].to);
        }
      }

      /* XXX bestimmen des vermutlich realen Frames - Verbesserung
         moeglich? XXX */
      for (row_idx = 0; row_idx < 7; row_idx++)
      {
        if ((frame_counter[row_idx] > real_frame_score)
            && (gt_array_size(regionmatrix[row_idx][0].from) > 0))
        {
          real_frame_score = frame_counter[row_idx];
          real_frame = row_idx;
        }
      }

      /* Aufruf der Funktion, in der auf moegliche Frame-Shifts geprueft
         wird */
      had_err =
        frameshiftprocessing(parsestruct_ptr, regionmatrix, real_frame,
                             err);

      if (!had_err)
      {
        /* Aufruf der Funktion, in der auf moegliche Vereinbarkeit
           kodierender Bereiche geprueft wird */
        genemergeprocessing(parsestruct_ptr, regionmatrix, err);

        /* Aufruf der Ausgabefunktion und Ausgabe der Ergebnisse */
        mg_outputwriter(parsestruct_ptr, combinedscore_matrix,
                        hit_information, regionmatrix, 'h', err);

        for (row_idx = 0; row_idx < 7; row_idx++)
        {
          gt_array_delete(regionmatrix[row_idx][0].from);
          gt_array_delete(regionmatrix[row_idx][0].to);
        }
      }
      gt_array2dim_delete(regionmatrix);
      gt_free(frame_counter);
    }
  }
  gt_array_delete(frame_path_array);

  return had_err;
}

static void gene_prediction(unsigned short row,
                            unsigned long column,
                            double max_lastcolumn,
                            PathMatrixEntry **path_matrix,
                            CombinedScoreMatrixEntry
                            **combinedscore_matrix,
                            GtArray * frame_path_array,
                            HitInformation *hit_information,
                            ParseStruct *parsestruct_ptr,
                            RegionStruct **regionmatrix,
                            unsigned long *frame_counter, GtError * err)
{
  unsigned long column_from,
    column_to;

  gt_error_check(err);

  /* die ersten beiden Spalten werden gesondert bearbeitet, da immer zwei
     Spalten zur Beurteilung, ob es sich um kodierende oder
     nicht-kodierende Positionen handelt, benoetigt werden. Zudem werden
     erst so die Berechnungen der Sequenzlaengen kodierender Bereiche bei
     Framewechsel und auch bei Wechsel von kodierendem zu
     nicht-kodierendem Bereich */
  if (column == gt_str_length(MATRIXSTRUCT(query_dna)) - 1)
  {
    GENEPREDSTRUCT(matrixscore_before) = path_matrix[row][column].score;
    GENEPREDSTRUCT(frame_before) = path_matrix[row][column].path_frame;

    /* hinzufuegen des aktuellen Leserahmens - als Zeile - in das
       frame_path-Array */
    gt_array_add(frame_path_array, row);
    /* Zaehlen der Frame-Haeufigkeiten */
    frame_counter[row] += 1.0;
  }
  else if (column == gt_str_length(MATRIXSTRUCT(query_dna)) - 2)
  {
    GENEPREDSTRUCT(matrixscore) = path_matrix[row][column].score;
    GENEPREDSTRUCT(current_frame) = path_matrix[row][column].path_frame;

    /* hinzufuegen des aktuellen Leserahmens - als Zeile - in das
       frame_path-Array */
    gt_array_add(frame_path_array, row);
    /* Zaehlen der Frame-Haeufigkeiten */
    frame_counter[row] += 1.0;
  }
  /* ab hier beginnt die eigentliche Bearbeitung des optimalen Pfades in
     der Combined-Score-Matrix - berechnet wird jeweils die Position
     column+2 */
  else
  {
    /* hinzufuegen des aktuellen Leserahmens - als Zeile - in das
       frame_path-Array */
    gt_array_add(frame_path_array, row);
    /* Zaehlen der Frame-Haeufigkeiten */
    frame_counter[row] += 1.0;

    /* nur wenn der vorherige Matrixscore(column+2) hoeher ist als der
       Aktuelle(column+1) ist, handelt es sich bei Position column+2 um
       eine kodierende Position */
    if (GENEPREDSTRUCT(matrixscore_before) - GENEPREDSTRUCT(matrixscore) >
        0)
    {
      /* Laenge des kodierenden Abschnittes um 1 erhoehen */
      ++GENEPREDSTRUCT(codingcounter);

      /* Wechsel des Leserahmens von Position column+2 zu Position
         column+1; kodierender zu kodierender Abschnitt */
      if (GENEPREDSTRUCT(current_frame) != GENEPREDSTRUCT(frame_before)
          && GENEPREDSTRUCT(noncodingcounter) == 0)
      {
        /* Speichern der Grenzen des kodierenden Bereichs vor dem
           Framewechsel */
        column_from = column + 2;
        column_to = column + GENEPREDSTRUCT(codingcounter) + 1;

        gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].from,
                  column_from);
        gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].to,
                  column_to);

        GENEPREDSTRUCT(codingcounter) = 0;
      }
      /* wenn noncodingcounter > 0 liegt ein Wechsel von nicht-kodierender
         zu kodierender Region vor */
      else if (GENEPREDSTRUCT(noncodingcounter) > 0)
      {
        GENEPREDSTRUCT(noncodingcounter) = 0;
      }
    }
    /* else-Fall: nicht-kodierender Positionen */
    else
    {
      /* Wechsel von kodierendem zu nicht-kodierendem Abschnitt */
      if (GENEPREDSTRUCT(codingcounter) > 0)
      {
        /* Speichern der Grenzen des vorangegangenen kodierenden Bereichs */
        column_from = column + 3;
        column_to = column + GENEPREDSTRUCT(codingcounter) + 2;

        gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].from,
                  column_from);
        gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].to,
                  column_to);

        GENEPREDSTRUCT(codingcounter) = 0;
      }
      GENEPREDSTRUCT(noncodingcounter)++;
    }

    /* Sichern der Vorgaengerinformationen (column+1) und einlesen der
       aktuellen Informationen (column) */
    GENEPREDSTRUCT(matrixscore_before) = GENEPREDSTRUCT(matrixscore);
    GENEPREDSTRUCT(frame_before) = GENEPREDSTRUCT(current_frame);

    GENEPREDSTRUCT(matrixscore) =
      path_matrix[GENEPREDSTRUCT(current_frame)][column].score;
    GENEPREDSTRUCT(current_frame) =
      path_matrix[GENEPREDSTRUCT(current_frame)][column].path_frame;
  }

  /* die ersten beiden Spalten der Combined-Score-Matrix werden in einem
     Fall bearbeitet */
  if (column == 0)
  {
    /* das Function Stop-Flag wird gesetzt */
    GENEPREDSTRUCT(function_stop) = 1;
    /* if-Fall: kodierende Base an Position 1 */

    if (GENEPREDSTRUCT(matrixscore_before) - GENEPREDSTRUCT(matrixscore) >
        0)
    {
      ++GENEPREDSTRUCT(codingcounter);

      /* Wechsel des Leserahmens - kodierender zu kodierender Abschnitt */
      if (GENEPREDSTRUCT(current_frame) != GENEPREDSTRUCT(frame_before)
          && GENEPREDSTRUCT(noncodingcounter == 0))
      {
        column_from = column + 1;
        column_to = column + GENEPREDSTRUCT(codingcounter);

        gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].from,
                  column_from);
        gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].to,
                  column_to);

        GENEPREDSTRUCT(codingcounter) = 0;
      }
      if (GENEPREDSTRUCT(noncodingcounter) > 0)
      {
        GENEPREDSTRUCT(noncodingcounter) = 0;
      }
    }
    /* else-Fall: nicht-kodiernde Base an Position 1 */
    else
    {
      /* der Bereich vor der nicht-kodiernden Position 1 ist kodierend */
      if (GENEPREDSTRUCT(codingcounter) > 0)
      {
        column_from = column + 2;
        column_to = column + GENEPREDSTRUCT(codingcounter) + 1;

        gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].from,
                  column_from);
        gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].to,
                  column_to);

        GENEPREDSTRUCT(codingcounter) = 0;
      }
      GENEPREDSTRUCT(noncodingcounter)++;
    }

    /* Wenn der Wert in der ersten Spalte der Combined-Score-Matrix
       positiv ist, ist die Position 0 kodierend */
    if (GENEPREDSTRUCT(matrixscore) > 0)
    {
      GENEPREDSTRUCT(codingcounter)++;

      /* Falls der noncodingcounter > 0 ist, ist der Bereich vor 0
         nicht-kodierend */
      if (GENEPREDSTRUCT(noncodingcounter) > 0)
      {
        column_from = 0;
        column_to = 0;

        gt_array_add(regionmatrix[GENEPREDSTRUCT(current_frame)][0].from,
                  column_from);
        gt_array_add(regionmatrix[GENEPREDSTRUCT(current_frame)][0].to,
                  column_to);

        GENEPREDSTRUCT(noncodingcounter) = 0;
      }
      else
      {
        /* if-Fall: Wechsel des Leserahmens von kodierendem zu kodierendem
           Abschnitt */
        if (GENEPREDSTRUCT(current_frame) != GENEPREDSTRUCT(frame_before))
        {
          /* Vorgaengerbereich von Position 0 ist kodierend */
          column_from = 1;
          column_to = GENEPREDSTRUCT(codingcounter) - 1;

          gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].from,
                    column_from);
          gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].to,
                    column_to);

          /* Position 0 ist ebenfalls kodierend */
          column_from = 0;
          column_to = 0;

          gt_array_add(regionmatrix[GENEPREDSTRUCT(current_frame)][0].from,
                    column_from);
          gt_array_add(regionmatrix[GENEPREDSTRUCT(current_frame)][0].to,
                    column_to);

          GENEPREDSTRUCT(codingcounter) = 1;
        }
        /* else-Fall: Speichern des kodierenden Bereichs von 0 bis x */
        else
        {
          column_from = 0;
          column_to = GENEPREDSTRUCT(codingcounter) - 1;

          gt_array_add(regionmatrix[GENEPREDSTRUCT(current_frame)][0].from,
                    column_from);
          gt_array_add(regionmatrix[GENEPREDSTRUCT(current_frame)][0].to,
                    column_to);
        }
      }
    }
    /* Position 0 ist nicht-kodierend */
    else
    {
      /* Falls codingcounter > 0, dann ist der vor 0 liegende Bereich
         kodierend */
      if (GENEPREDSTRUCT(codingcounter) > 0)
      {
        column_from = 1;
        column_to = column + GENEPREDSTRUCT(codingcounter);

        gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].from,
                  column_from);
        gt_array_add(regionmatrix[GENEPREDSTRUCT(frame_before)][0].to,
                  column_to);

        GENEPREDSTRUCT(codingcounter) = 0;
      }
    }
    GENEPREDSTRUCT(noncodingcounter) = 0;
    GENEPREDSTRUCT(codingcounter) = 0;
  }
  /* solange das Stop-Flag nicht gesetzt ist wird die CombinedScore-Matrix
     rekursiv spaltenweise abgearbeitet */
  if (!GENEPREDSTRUCT(function_stop))
  {
    if (column > 0)
      gene_prediction(path_matrix[row][column].path_frame,
                      column - 1,
                      max_lastcolumn,
                      path_matrix,
                      combinedscore_matrix,
                      frame_path_array,
                      hit_information,
                      parsestruct_ptr, regionmatrix, frame_counter, err);
    else
      gene_prediction(path_matrix[row][column].path_frame,
                      column,
                      max_lastcolumn,
                      path_matrix,
                      combinedscore_matrix,
                      frame_path_array,
                      hit_information,
                      parsestruct_ptr, regionmatrix, frame_counter, err);
  }
  else
    GENEPREDSTRUCT(function_stop) = 0;
}

static int frameshiftprocessing(ParseStruct *parsestruct_ptr,
                                RegionStruct **regionmatrix,
                                short real_frame, GtError * err)
{
  int had_err = 0;

  unsigned short row_index,
    check_bp = 0;

  unsigned long arraysize,
    arraysize_realframe,
    GT_UNUSED arraysize_real,
    arraysize_tmp,
    array_idx,
    arrayreal_idx = 0,
    from_tmp,
    to_tmp,
    from_real,
    to_real,
    from_min = 0,
    to_min = 0;

  long min_value,
    min_value_tmp = LONG_MAX;

  GtArray *tmp_from;
  GtArray *tmp_to;
  GtArray *real_from;
  GtArray *real_to;
  GtArray *real_fromtmp;
  GtArray *real_totmp;
  GtArray *realfrom;
  GtArray *realto;

  gt_error_check(err);

  tmp_from = gt_array_new(sizeof (unsigned long));
  tmp_to = gt_array_new(sizeof (unsigned long));
  real_from = gt_array_new(sizeof (unsigned long));
  real_to = gt_array_new(sizeof (unsigned long));
  real_fromtmp = gt_array_new(sizeof (unsigned long));
  real_totmp = gt_array_new(sizeof (unsigned long));
  realfrom = gt_array_new(sizeof (unsigned long));
  realto = gt_array_new(sizeof (unsigned long));

  /* Vergleich mit allen anderen Frames - Der reale Frame ist der Frame
     mit der hoechsten Haeufigkeit (frame_counter) */
  for (row_index = 0; row_index < 7; row_index++)
  {
    /* Frameshift nur zwischen realem und von diesem verschiedenen Frames
       moeglich */
    if (row_index != real_frame)
    {
      /* Anzahl der kodierenden Bereiche im zu vergleichendem Leserahmen */
      arraysize = gt_array_size(regionmatrix[row_index][0].from);

      if (arraysize > 0)
      {
        /* alle kodierenden Bereiche des zu vergleichenden Frames werden
           betrachtet */
        for (array_idx = 0; array_idx < arraysize; array_idx++)
        {
          /* Grenzen des zu Betrachtenden Abschnitts - Abschnitte liegen
             in aufsteigender Reihenfolge vor */
          from_tmp =
            *(unsigned long *) gt_array_get(regionmatrix[row_index][0].from,
                                            array_idx);
          to_tmp =
            *(unsigned long *) gt_array_get(regionmatrix[row_index][0].to,
                                            array_idx);

          /* Aufruf der Methode zur Ueberpruefung auf Stop-Codons im
             kodierendem Bereich des zu vergleichendem Frames -
             ueberprueft wird der kodierende Abschnitt bezgl. Stop-Codons
             im realen Frame */
          check_bp = check_coding(parsestruct_ptr,
                                  from_tmp, to_tmp, real_frame, err);

          if (check_bp || !check_bp)
          {
            arraysize_realframe =
              gt_array_size(regionmatrix[real_frame][0].from);

            /* alle kodierenden Bereiche des realen Frames werden
               betrachtet */
            for (arrayreal_idx = 0; arrayreal_idx < arraysize_realframe;
                 arrayreal_idx++)
            {
              /* Grenzen des kodierenden Bereichs des realen Frames */
              from_real =
                *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].
                                             from, arrayreal_idx);
              to_real =
                *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].
                                             to, arrayreal_idx);

              /* Bestimmen der Frames von realem und zu vergleichendem
                 Frame mit minimalem Abstand */
              if (from_real > from_tmp)
              {
                min_value = from_real - from_tmp;
              }
              else
                min_value = from_tmp - from_real;

              /* kodierender Bereich mit minimalem Abstand wird
                 gespeichert */
              if (min_value < min_value_tmp)
              {
                min_value_tmp = min_value;
                from_min = from_real =
                  *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].
                                               from, arrayreal_idx);
                to_min = to_real =
                  *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].
                                               to, arrayreal_idx);
              }
            }

            /* Falls der betrachtete Abschnitt des zu vergleichenden
               Frames innerhalb des max. Abstandes liegt (per GtOption bei
               Programmaufruf einstellbar) und der kodierende Bereich kein
               Stop-Codon enthaelt werden die Grenzen des Abschnittes im
               real_fromtmp- bzw. real_totmp-Array abgespeichert */
            if (!check_bp
                && ((to_tmp + ARGUMENTSSTRUCT(frameshift_span) > from_min)
                    || (to_min + ARGUMENTSSTRUCT(frameshift_span) >
                        from_tmp)))
            {
              gt_array_add(real_fromtmp, from_tmp);
              gt_array_add(real_totmp, to_tmp);
            }
            /* ansonsten im from_tmp- bzw. to_tmp-Array */
            else
            {
              gt_array_add(tmp_from, from_tmp);
              gt_array_add(tmp_to, to_tmp);
            }
          }
          else
          {
            had_err = -1;
          }
        }
      }

      if (!had_err)
      {
        /* die from- to-Grenzen des zu vergleichenden GtArrays werden
           aktualisiert */
        gt_array_delete(regionmatrix[row_index][0].from);
        gt_array_delete(regionmatrix[row_index][0].to);

        regionmatrix[row_index][0].from = gt_array_clone(tmp_from);
        regionmatrix[row_index][0].to = gt_array_clone(tmp_to);
        gt_array_reset(tmp_from);
        gt_array_reset(tmp_to);
      }
    }
  }

  if (!had_err)
  {
    /* Anzahl der kodierenden Bereiche im realen- und im
       real_fromtmp-Array werden bestimmt */
    arraysize_real = gt_array_size(regionmatrix[real_frame][0].from);
    arraysize_tmp = gt_array_size(real_fromtmp);

    /* Falls im tmp-Array Eintraege vorhanden sind, sind diese in das
       GtArray der kodierenden Abschnitte der realen-Frames einzutragen */
    if (arraysize_tmp > 0)
    {
      /* sortieren der neu einzutragenden kodierenden Abschnitte anhand
         des From-Wertes */
      sort_realtmp(realfrom, realto, real_fromtmp, real_totmp);

      /* bisher: rueckwaerts von hinten angegebene Bereichsgrenzen -
         danach: Bereichsgrenzen von vorne nach hinten */
      gt_array_reverse(realfrom);
      gt_array_reverse(realto);

      /* sortiertes einfuegen der neuen Bereichsgrenzen in das GtArray
         bestehender kodierender Bereiche des realen Frames */
      merge_array(regionmatrix,
                  real_from,
                  real_to, realfrom, realto, 0, 0, real_frame, row_index);

      /* aktualisieren der Bereichsgrenzen im GtArray des realen-Frames */
      gt_array_delete(regionmatrix[real_frame][0].from);
      gt_array_delete(regionmatrix[real_frame][0].to);

      regionmatrix[real_frame][0].from = gt_array_clone(real_from);
      regionmatrix[real_frame][0].to = gt_array_clone(real_to);
    }
  }

  gt_array_delete(tmp_from);
  gt_array_delete(tmp_to);
  gt_array_delete(real_from);
  gt_array_delete(real_to);
  gt_array_delete(real_fromtmp);
  gt_array_delete(real_totmp);
  gt_array_delete(realfrom);
  gt_array_delete(realto);

  return had_err;
}

static void genemergeprocessing(ParseStruct *parsestruct_ptr,
                                RegionStruct **regionmatrix, GtError * err)
{
  unsigned short row_index;
  short check_bp;
  unsigned long array_idx,
    arraysize,
    from_tmp,
    to_tmp,
    from_tmp_next,
    to_tmp_real,
    function_stop = 0;
  GtArray *tmp_from;
  GtArray *tmp_to;

  gt_error_check(err);

  tmp_from = gt_array_new(sizeof (unsigned long));
  tmp_to = gt_array_new(sizeof (unsigned long));

  /* Frames werden nacheinander abgearbeitet - zeilenkodiert */
  for (row_index = 0; row_index < 7; row_index++)
  {
    /* bestimmen der Anzahl kodierender Bereiche im aktuell betrachteten
       Frame */
    arraysize = gt_array_size(regionmatrix[row_index][0].from);

    gt_array_reset(tmp_from);
    gt_array_reset(tmp_to);

    if (arraysize > 1)
    {
      /* die kodierenden Bereiche werden nacheinander abgearbeitet */
      for (array_idx = 0; array_idx < arraysize - 1;)
      {
        /* ist gt_array_size des tmp-Arrays > 0 wurden bereits Aenderungen
           vorgenommen und die neuen Bereichsgrenzen-Werte muessen
           verwendet werden */
        if (gt_array_size(tmp_from) > 0)
        {
          from_tmp = *(unsigned long *) gt_array_get_last(tmp_from);
          to_tmp = *(unsigned long *) gt_array_get_last(tmp_to);
        }
        /* ansonsten auslesen der Bereichsgrenzen aus dem "original" */
        else
        {
          from_tmp =
            *(unsigned long *) gt_array_get(regionmatrix[row_index][0].from,
                                            array_idx);
          to_tmp =
            *(unsigned long *) gt_array_get(regionmatrix[row_index][0].to,
                                            array_idx);
        }

        ++array_idx;
        function_stop = 0;

        do
        {
          /* auslesen der Grenzen des nachfolgenden kodierenden Bereiches */
          from_tmp_next =
            *(unsigned long *) gt_array_get(regionmatrix[row_index][0].from,
                                            array_idx);
          to_tmp_real =
            *(unsigned long *) gt_array_get(regionmatrix[row_index][0].to,
                                            array_idx);

          /* zusammenlegen kodierender Bereiche nur, wenn diese innerhalb
             der per GtOption angegebenen Spanbreite liegen */
          if (to_tmp + ARGUMENTSSTRUCT(prediction_span) > from_tmp_next)
          {
            /* es wurden noch keine Eintragungen vorgenommen - der erste
               kodierende Abschnitt wird in das tmp-Array eingetragen */
            if (gt_array_size(tmp_from) == 0)
            {
              gt_array_add(tmp_from, from_tmp);
              gt_array_add(tmp_to, to_tmp);
            }

            /* Ueberpruefen, ob zwischen den kodierenden Abschnitten
               Stop-Codons liegen */
            check_bp = check_coding(parsestruct_ptr,
                                    to_tmp, from_tmp_next, row_index, err);

            /* kein Stop-Codon gefunden */
            if (!check_bp)
            {
              /* befindet sich schon ein Eintrag im temp-Array, muss der
                 letzte Eintrag wieder entfernt werden */
              if (gt_array_size(tmp_from) > 0)
              {
                (void) gt_array_pop(tmp_from);
                (void) gt_array_pop(tmp_to);
              }
              /* Speichern der neuen Bereichsgrenzen */
              gt_array_add(tmp_from, from_tmp);
              gt_array_add(tmp_to, to_tmp_real);

              ++array_idx;
            }
            /* Stop-Codon gefunden - der next-Abschnitt wird ebenfalls in
               das tmp-Array eingetragen; die Bearbeitung kann hier
               abgebrochen werden, da das Stop-Codon auch im Vergleich zu
               den anderen Abschnitten eine Vereinigung verhindert */
            else
            {
              gt_array_add(tmp_from, from_tmp_next);
              gt_array_add(tmp_to, to_tmp_real);

              function_stop = 1;
            }
          }
          /* Abstand der kodierenden Bereiche > als angegebene
             Spannbreite, auch hier kann die Bearbeitung beendet werden */
          else
          {
            gt_array_add(tmp_from, from_tmp_next);
            gt_array_add(tmp_to, to_tmp_real);

            function_stop = 1;
          }
          /* Bearbeitung solange, bis keine weiteren Abschnitte vorliegen
             oder das Stop-Flag gesetzt ist */
        } while (array_idx < arraysize && !function_stop);
      }

      /* aktualisieren der Bereichsgrenzen */
      gt_array_delete(regionmatrix[row_index][0].from);
      gt_array_delete(regionmatrix[row_index][0].to);
      regionmatrix[row_index][0].from = gt_array_clone(tmp_from);
      regionmatrix[row_index][0].to = gt_array_clone(tmp_to);
    }
  }
  gt_array_delete(tmp_from);
  gt_array_delete(tmp_to);
}

static int check_coding(ParseStruct *parsestruct_ptr,
                        unsigned long from,
                        unsigned long to, short current_row, GtError * err)
{
  int had_err = 0;

  unsigned long startpoint = 0,
    endpoint = 0,
    contig_len = 0;

  long diff;

  short current_frame = 0,
    found = 0;

  GtStr *query_seq;

  char *contig_seq_ptr = NULL,
    contig_seq_tri[4] = { '\0', '\0', '\0', '\0' };

  gt_error_check(err);

  /* Start- und Endsequenzpositionen der in AS umzuwandelnden Sequenz */
  startpoint = from;
  endpoint = to;

  contig_len = gt_str_length(MATRIXSTRUCT(query_dna));
  query_seq = gt_str_new_cstr(gt_str_get(MATRIXSTRUCT(query_dna)));
  contig_seq_ptr = gt_str_get(query_seq);

  /* Bestimmung des aktuellen Frames aus der Zeilennummer */
  current_frame = get_current_frame(current_row);

  /* Bestimmen der Laenge der umzuwandelnden Sequenz */
  diff = startpoint - endpoint;
  if (diff < 0)
    diff = (-1) * diff;

  /* nur bei 3 oder mehr DNA-Basen ist eine Umwandlung in eine AS-Sequenz
     sinnvoll bzw. machbar */
  if (!(diff < 3))
  {
    /* if: aktueller Leserahmen ist -1, -2 oder -3 */
    if (current_frame < 0)
    {
      /* fuer die weiteren Berechnungen erfolgt eine Multiplikation mit -1
       */
      current_frame *= -1;
      /* das reverse Komplement der Sequenz wird gebildet - erweiterte
         reverse_complement-Funktion beruecksichtigt das erweiterte
         DNA-Alphabet */
      had_err = mg_reverse_complement(contig_seq_ptr, contig_len, err);

      /* Der Startpunkt muss im Fall negativer Leserahmen neu berechnet
         werden */
      startpoint = contig_len - 1 - to;
      endpoint = contig_len - from;
    }

    if (!had_err)
    {
      /* ist der Startpunkt < 3 ist die Startposition 0, 1 oder 2 */
      if (startpoint < 3)
      {
        startpoint = current_frame - 1;
      }
      /* sonst: Berechnung der Startposition anhand der gegebenen Formel */
      else
      {
        startpoint -= (((startpoint) - current_frame) % 3);
        startpoint -= 1;
      }

      /* abschreiten der Sequenz */
      while ((startpoint <= endpoint - 2) && !found)
      {
        /* DNA-Basen-Triplet einlesen */

        contig_seq_tri[0] = tolower(contig_seq_ptr[startpoint]);
        contig_seq_tri[1] =
          tolower(contig_seq_ptr[startpoint + 1]);
        contig_seq_tri[2] =
          tolower(contig_seq_ptr[startpoint + 2]);

        found = gt_check_stopcodon(contig_seq_tri);

        /* Startwert um 3 Basen weitersetzen */
        startpoint += 3;
      }
    }
    else
    {
      found = -1;
    }
  }
  gt_str_delete(query_seq);

  return found;

}

static void merge_array(RegionStruct **regionmatrix,
                        GtArray * real_from_ar,
                        GtArray * real_to_ar,
                        GtArray * real_fromtmp,
                        GtArray * real_totmp,
                        unsigned long real_index,
                        unsigned long tmp_index,
                        short real_frame, unsigned short row_index)
{
  unsigned long real_from,
    real_to,
    tmp_from,
    tmp_to;

  /* Fall 1: Eintraege im GtArray des realen Frames sind abgearbeitet,
     Eintraege im tmp-Array noch nicht */
  if (!(real_index < gt_array_size(regionmatrix[real_frame][0].from))
      && tmp_index < gt_array_size(real_fromtmp))
  {
    /* Bereichsgrenzen des Realtmp-Arrays - from/to-Werte */
    real_from = *(unsigned long *) gt_array_get(real_fromtmp, tmp_index);
    real_to = *(unsigned long *) gt_array_get(real_totmp, tmp_index);

    /* solange sich noch Eintraege im realtmp-Array befinden */
    while (tmp_index < gt_array_size(real_fromtmp))
    {
      /* hinzufuegen der Bereichsgrenzen */
      gt_array_add(real_from_ar, real_from);
      gt_array_add(real_to_ar, real_to);

      tmp_index++;

      /* nur wenn noch weitere Eintraege vorhanden sind, werden die
         naechsten Bereichsgrenzen ausgelesen */
      if (tmp_index < gt_array_size(real_fromtmp))
      {
        real_from = *(unsigned long *) gt_array_get(real_fromtmp, tmp_index);
        real_to = *(unsigned long *) gt_array_get(real_totmp, tmp_index);
      }
    }
  }
  /* Fall 2: Eintraege im GtArray des realen Frames sind noch nicht
     abgearbeitet, Eintraege im tmp-Array abgearbeitet */
  else if (!(tmp_index < gt_array_size(real_fromtmp))
           && real_index < gt_array_size(regionmatrix[real_frame][0].from))
  {
    /* Bereichsgrenzen des Real-Arrays - from/to-Werte */
    tmp_from =
      *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].from,
                                   real_index);
    tmp_to =
      *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].to,
                                   real_index);

    /* solange sich noch Eintraege im real-Array befinden */
    while (real_index < gt_array_size(regionmatrix[real_frame][0].from))
    {
      /* hinzufuegen der Bereichsgrenzen */
      gt_array_add(real_from_ar, tmp_from);
      gt_array_add(real_to_ar, tmp_to);

      real_index++;

      /* nur wenn noch weitere Eintraege vorhanden sind, werden die
         naechsten Bereichsgrenzen ausgelesen */
      if (real_index < gt_array_size(regionmatrix[real_frame][0].from))
      {
        tmp_from =
          *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].from,
                                       real_index);
        tmp_to =
          *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].to,
                                       real_index);
      }
    }
  }
  /* Fall 3: Es sind noch Eintraege in beiden GtArrays vorhanden */
  else if ((tmp_index < gt_array_size(real_fromtmp))
           && real_index < gt_array_size(regionmatrix[real_frame][0].from))
  {
    /* Bereichsgrenzen des Realtmp- und des real-Arrays - from/to-Werte */
    tmp_from =
      *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].from,
                                   real_index);
    tmp_to =
      *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].to,
                                   real_index);
    real_from = *(unsigned long *) gt_array_get(real_fromtmp, tmp_index);
    real_to = *(unsigned long *) gt_array_get(real_totmp, tmp_index);

    /* Fall: realtmp-Abschnitt liegt vor dem real_frame Abschnitt */
    if (real_from < tmp_from)
    {
      /* solange es noch Eintraege gibt, die vor dem real_frame-Abschnitt
         liegen... */
      while (tmp_index < gt_array_size(real_fromtmp) && real_from < tmp_from)
      {
        /* ... werden die bereichsgrenzen Eingetragen... */
        gt_array_add(real_from_ar, real_from);
        gt_array_add(real_to_ar, real_to);

        tmp_index++;

        /* und der naechste realtmp-Eintrag eingelesen */
        if (tmp_index < gt_array_size(real_fromtmp))
        {
          real_from =
            *(unsigned long *) gt_array_get(real_fromtmp, tmp_index);
          real_to = *(unsigned long *) gt_array_get(real_totmp, tmp_index);
        }
      }
      /* rekursiver Aufruf der Sortierfunktion unter Beruecksichtigung der
         neu Eingetragenen Abschnitte */
      merge_array(regionmatrix,
                  real_from_ar,
                  real_to_ar,
                  real_fromtmp,
                  real_totmp,
                  real_index, tmp_index, real_frame, row_index);
    }
    /* real_frame- vor dem realtmp-Abschnitt */
    else if (tmp_from < real_from)
    {
      /* solange es noch Eintraege gibt, die vor dem realtmp-Abschnitt
         liegen... */
      while (real_index < gt_array_size(regionmatrix[real_frame][0].from)
             && tmp_from < real_from)
      {
        /* ... werden die bereichsgrenzen Eingetragen... */
        gt_array_add(real_from_ar, tmp_from);
        gt_array_add(real_to_ar, tmp_to);

        real_index++;

        /* und der naechste real_frame-Eintrag eingelesen */
        if (real_index < gt_array_size(regionmatrix[real_frame][0].from))
        {
          tmp_from =
            *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].from,
                                         real_index);
          tmp_to =
            *(unsigned long *) gt_array_get(regionmatrix[real_frame][0].to,
                                         real_index);
        }
      }
      /* rekursiver Aufruf der Sortierfunktion unter Beruecksichtigung der
         neu Eingetragenen Abschnitte */
      merge_array(regionmatrix,
                  real_from_ar,
                  real_to_ar,
                  real_fromtmp,
                  real_totmp,
                  real_index, tmp_index, real_frame, row_index);
    }
  }
}

static void sort_realtmp(GtArray * realfrom,
                         GtArray * realto,
                         GtArray * real_fromtmp, GtArray * real_totmp)
{
  unsigned long index_outer,
    index_inner,
    max_value = 0,
    from,
    to,
    from_tmp = 0,
    to_tmp = 0;

  /* durchlaufen des GtArrays mit einer aeusseren und einer inneren Schleife
   */
  for (index_outer = 0; index_outer < gt_array_size(real_fromtmp);
       index_outer++)
  {
    /* die innere Schleife bestimmt den naechsten maximalen Wert im GtArray */
    for (index_inner = 0; index_inner < gt_array_size(real_fromtmp);
         index_inner++)
    {
      from = *(unsigned long *) gt_array_get(real_fromtmp, index_inner);
      to = *(unsigned long *) gt_array_get(real_totmp, index_inner);

      /* Fall: es gibt bereits Eintraege im tmp-Array */
      if (gt_array_size(realfrom) > 0)
      {
        /* der aktuelle from-Wert ist groesser als der tmp-Wert und
           kleiner als der letzte max-Wert */
        if ((from > from_tmp) && (from < max_value))
        {
          from_tmp = from;
          to_tmp = to;
        }
      }
      /* Fall: keine Eintraege im tmp-Array */
      else
      {
        /* bestimmen des max-Wertes */
        if (from > from_tmp)
        {
          from_tmp = from;
          to_tmp = to;
        }
      }
    }
    /* Speichern des max-Wertes */
    gt_array_add(realfrom, from_tmp);
    gt_array_add(realto, to_tmp);

    /* merken des letzten max-Wertes */
    max_value = from_tmp;
    from_tmp = 0;
    to_tmp = 0;
  }
}

short gt_check_stopcodon(char *contig_seq_ptrfct)
{
  unsigned short codon_status = 0;

  /* jeder if-Zweig ueberprueft entsprechend des
     Metagenomethreader-Arguments zur Verwendung alternativer Start-Codons
     auf Start-Codons */
  if (!strcmp(contig_seq_ptrfct, "tga")
      || !strcmp(contig_seq_ptrfct, "taa")
      || !strcmp(contig_seq_ptrfct, "tag")
      || !strcmp(contig_seq_ptrfct, "tar")
      || !strcmp(contig_seq_ptrfct, "uga")
      || !strcmp(contig_seq_ptrfct, "uaa")
      || !strcmp(contig_seq_ptrfct, "uag")
      || !strcmp(contig_seq_ptrfct, "uar"))
  {
    codon_status = 1;
  }

  return codon_status;
}
