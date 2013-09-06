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

#include "mg_computepath.h"
#include "core/mathsupport.h"

/* Funktion zur Berechnung der erlaubten Vorgaenger-Leserahmen
   Parameter: aktueller Leserahmen, Position in der Query-Sequence,
   (Zeiger auf) GtArray der moeglichen Vorgaenger
   Returnwert: void */
static void compute_precursors(short, GtUword, short *);

enum {
  NUM_PRECURSORS = 3,
};

int mg_computepath(CombinedScoreMatrixEntry **combinedscore_matrix,
                   HitInformation *hit_information,
                   GtUword rows,
                   GtUword contig_len,
                   ParseStruct *parsestruct_ptr, GtError * err)
{
  int had_err = 0;

  /* Initialisieren der Matrix fuer die Pfadberechnung */
  PathMatrixEntry **path_matrix;

  /* i: Zaehlvariable fuer die Matrix-Zeilen; k: Zaehlvariable Precursors
     (von 0 bis max 2) maxpath_frame: Speichern des vorherigen Frames von
     dem der max-Wert berechnet wird */
  unsigned short row_index = 0,
    precursor_index = 0,
    precursors_row = 0,
    maxpath_frame = 0;

  /* Position in der Query-DNA */
  GtUword column_index = 0;

  /* Variablen fuer den aktuellen Frame, den vorherigen Frame(speichert
     einen Wert aus precursors[], die Zeile des vorherigen Frames, GtArray
     mit den Precursors-Frames */
  short current_frame = 0,
    precursors_frame = 0,
    precursors[NUM_PRECURSORS];

  /* q ist der Wert, der bei Aus- oder Eintreten in ein Gen auf dem
     Forward- bzw. Reverse-Strang berechnet wird */
  double q = ARGUMENTSSTRUCT(leavegene_value),
    max_new = 1,
    max_old = 1;

  /* Speicherreservierung fuer die Path-Matrix - Groesse entsprechend der
     CombinedScore-Matrix */
  gt_array2dim_calloc(path_matrix, 7, contig_len);

  gt_error_check(err);

  /* fuer die erste Spalte der Path-Matrix wird die erste Spalte der
     CombinedScore-Matrix uebernommen */
  for (row_index = 0; row_index < rows; row_index++)
  {
    path_matrix[row_index][0].score =
      combinedscore_matrix[row_index][0].matrix_score;
    path_matrix[row_index][0].path_frame = row_index;
  }

  /* Spaltenweise Berechnung des opt. Pfades */
  for (column_index = 1; column_index < contig_len; column_index++)
  {
    for (row_index = 0; row_index < rows; row_index++)
    {
      /* Zaehlvariable fuer die Zeile wird umgerechnet in den entsprechenden
         Leserahmen */
      current_frame = get_current_frame(row_index);
      /* Aufruf der Methode zum Berechnen der moeglichen Leserahmen anhand von
         aktuellem Leserahmen und der Query-DNA-Sequenz */
      compute_precursors(current_frame,
                         column_index,
                         precursors);

      /* der max-Wert der moeglichen Vorgaenger wird berechnet */
      for (precursor_index = 0;
           precursor_index < NUM_PRECURSORS
             && (precursors[precursor_index] != UNDEFINED);
           ++precursor_index)
      {
        /* aktueller Vorgaengerleserahmen - es gibt max. 3 moegliche
           Vorgaenger */
        precursors_frame = precursors[precursor_index];
        /* Vorgaengerleserahmen wird umgerechnet in die entsprechende
           Matrix-Zeile */
        precursors_row = get_matrix_row(precursors_frame);

        /* der DP-Algo umfasst 3 moegliche Faelle
           1. Fall: Wechsel vom Reversen- auf den Forward-Strang bzw.
           umgekehrt */
        if ((current_frame < 0 && precursors_frame > 0) ||
            (current_frame > 0 && precursors_frame < 0))
        {
            max_new = path_matrix[precursors_row][column_index-1].score +
                      combinedscore_matrix[row_index][column_index].matrix_score
                      + 2*q;
        }
        /* 2. Fall: Einfacher Wechsel des Leserahmens, also von + zu +
           bzw.- zu - */
        else if (current_frame != 0 && precursors_frame != current_frame)
        {
            max_new = path_matrix[precursors_row][column_index-1].score +
                      combinedscore_matrix[row_index][column_index].matrix_score
                      + q;
        }
        /* 3. Fall: Leserahmen wird beibehalten bzw. Wechsel von kodierend zu
           nicht-kodierend oder umgekehrt */
        else
        {
            max_new = path_matrix[precursors_row][column_index-1].score +
                      combinedscore_matrix[row_index][column_index]
                      .matrix_score;
        }

        /* Bestimmen des Max-Wertes der max. 3 Moeglichkeiten und Speichern der
           Zeile, von der der Max-Wert stammt */
        if (gt_double_compare(max_new, max_old) > 0)
        {
            max_old = max_new;
            maxpath_frame = precursors_row;
        }
      }

      /* Speichern des Max-Wertes und der "Vorgaenger"-Zeile;
         zuruecksetzen der Variablen */
      path_matrix[row_index][column_index].score      = max_old;
      path_matrix[row_index][column_index].path_frame = maxpath_frame;

      max_new = DBL_MIN;
      max_old = DBL_MIN;
      maxpath_frame = 0;
    }
  }

  /* Aufruf der Methode zur Genvorhersage */
  had_err = mg_compute_gene_prediction(combinedscore_matrix,
                                       path_matrix,
                                       contig_len,
                                       hit_information,
                                       parsestruct_ptr, err);

  gt_array2dim_delete(path_matrix);

  return had_err;
}

/* Methode zur Berechnung der moeglichen Vorgaengerleserahmen */
static void compute_precursors(short current_frame,
                               GtUword position,
                               short *precursors_fct)
{
  short j = 0;

  /* Formel zur Bestimmung moeglicher Vorgaengerleserahmen anhand der
     aktuellen Position in der Query-DNA */
  j = (position)%NUM_PRECURSORS + 1;

  /* 3 Faelle zur Bestimmung der Vorgaengermenge */
  if (current_frame == 0)
  {
    precursors_fct[0] = j;
    precursors_fct[1] = 0;
    precursors_fct[2] = (-1)*j;
  }
  else if ((current_frame > 0 ? current_frame : (-1) * current_frame) == j)
  {
    precursors_fct[0] = current_frame;
    precursors_fct[1] = 0;
    precursors_fct[2] = (-1)*current_frame;
  }
  else
  {
    precursors_fct[0] = current_frame;
    precursors_fct[1] = UNDEFINED;
    precursors_fct[2] = UNDEFINED;
  }
}
