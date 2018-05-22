#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "core/ma.h"
#include "core/minmax.h"
#include "affine_dband.h"

typedef int GtAffineScore;

typedef enum
{
  Affine_X = 0, /* unknown, keep it, as it is the default value 0 */
  Affine_R, /* 1 = 01, shift 0 */
  Affine_D, /* 2 = 10, shift 2 */
  Affine_I  /* 3 = 11, shift 4 */
} GtAffineAlignEditOp;

typedef struct
{
  GtAffineScore Rvalue,
                Dvalue,
                Ivalue;
} GtAffineScoreTriple;

/* <GtAffineAlignTracebits> objects describe the backtracing edges
   relating to the last edit operation R,D,I. We store all three values
   between 0 and 3 as bitpairs in one uint8_t-value, using the shift
   values as described above. */

typedef struct
{
  uint8_t trace;
} GtAffineAlignTracebits;

#define GT_AFFINE_SCORE_SUM(V1, V2) ((V1) + (V2))

#define GT_AFFINE_SET_MAX2(MAXVALUE,MAXEDGE,SCORE_A,EDGE_A,SCORE_B,EDGE_B)\
        if (SCORE_B <= SCORE_A)\
        {\
          MAXVALUE = SCORE_A;\
          MAXEDGE = EDGE_A;\
        } else\
        {\
          MAXVALUE = SCORE_B;\
          MAXEDGE = EDGE_B;\
        }

#define GT_AFFINE_EDGE_SET(I,J,VALUE,EDGETYPE)\
        bitmatrix[J][I].trace = ((VALUE) << (2 * ((EDGETYPE)-1)))

#define GT_AFFINE_EDGE_SET_ALL(I,J,RVALUE,DVALUE,IVALUE)\
        bitmatrix[J][I].trace = (((RVALUE) << (2 * (Affine_R-1))) |\
                                 ((DVALUE) << (2 * (Affine_D-1))) |\
                                 ((IVALUE) << (2 * (Affine_I-1))))

#define GT_AFFINE_EDGE_GET(I,J,EDGETYPE)\
        ((bitmatrix[J][I].trace >> (2 * ((EDGETYPE)-1))) & 3)

#define GT_SHOW_AFFINE_DP_COLUMNS(CODELINES) /* Nothing */

GT_SHOW_AFFINE_DP_COLUMNS(
static GtUword show_dband_value(GtUword current_idx,GtUword j,
                                const GtAffineScoreTriple *currentcol,
                                GtUword low_row,GtUword high_row)
{
  GtUword idx;

  printf("column " GT_WU " in range " GT_WU "," GT_WU "\n",j,low_row,high_row);
  for (idx = low_row; idx <= high_row; idx++)
  {
    printf(GT_WU " " GT_WU " %d %d %d\n",
           current_idx++,idx,
           currentcol[idx].Rvalue,
           currentcol[idx].Dvalue,
           currentcol[idx].Ivalue);
  }
  return current_idx;
}
)/*GT_SHOW_AFFINE_DP_COLUMNS*/

#define GT_AFFINE_SCORE_ROW_GET(CC) (ISSPECIAL(CC) ? NULL\
                                                   : scorematrix2D[(int) (CC)])
#define GT_AFFINE_SCORE_ROW_ACCESS(CC)\
        ((score_row == NULL || ISSPECIAL(CC)) \
           ? smallest_score\
           : (GtAffineScore) score_row[(int) (CC)])
#define GT_AFFINE_EQUAL_SYMBOLS(CA,CB)\
        (!ISSPECIAL(CA) && (CA) == (CB))

static GtAffineScore gt_affine_diagonalband_fillDPtab_bits(
                                          GtAffineAlignTracebits **bitmatrix,
                                          GtAffineScoreTriple *currentcol,
                                          int8_t gap_opening, /* > 0 */
                                          int8_t gap_extension, /* > 0 */
                                          const int8_t * const *scorematrix2D,
#ifndef NDEBUG
#else
                                          __attribute__ ((unused))
#endif
                                          GtAffineScore min_align_score,
                                          GtAffineScore smallest_score,
                                          const GtUchar *useq,
                                          GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vlen,
                                          GtWord left_dist,
                                          GtWord right_dist)
{
  GtUword i, j, low_row = 0, high_row;
  const GtAffineScore start_penalty = -(gap_opening+gap_extension);
  GtAffineAlignTracebits *colptr;
  GtUword evaluated_cells = 0;
  GT_SHOW_AFFINE_DP_COLUMNS(
  GtUword current_idx = 0;
  )/*GT_SHOW_AFFINE_DP_COLUMNS*/
#ifndef NDEBUG
  const GtUword band_width = (GtUword) (right_dist - left_dist + 1);
  const GtWord lendiff = (GtWord) vlen - (GtWord) ulen;
#endif

  gt_assert(bitmatrix != NULL && currentcol != NULL &&
            gap_opening > 0 && gap_extension > 0 &&
            left_dist <= MIN(0, lendiff) && left_dist >= -(GtWord) ulen &&
            right_dist >= MAX(0, lendiff) && right_dist <= (GtWord) vlen);
  high_row = (GtUword) -left_dist;

  currentcol->Rvalue = 0;
  currentcol->Dvalue = -gap_opening;
  currentcol->Ivalue = -gap_opening;
  colptr = bitmatrix[0];
  for (i = 1; i <= high_row; i++)
  {
    GT_AFFINE_EDGE_SET(i,0,Affine_D,Affine_D);
    currentcol[i].Rvalue = min_align_score;
    currentcol[i].Dvalue = GT_AFFINE_SCORE_SUM(currentcol[i-1].Dvalue,
                                               -gap_extension);
    currentcol[i].Ivalue = min_align_score;
  }
  evaluated_cells += (1+high_row);

  GT_SHOW_AFFINE_DP_COLUMNS(
  printf("\n");
  )/*GT_SHOW_AFFINE_DP_COLUMNS*/
  for (j = 1; j <= vlen; j++)
  {
    const GtAffineScoreTriple *cptr;
    GtAffineScoreTriple nw;
    const GtUchar cb = vseq[j-1];
    const int8_t *score_row = GT_AFFINE_SCORE_ROW_GET(cb);

    GtAffineScore first_ivalue = min_align_score, score_R, score_D, score_I;
    const GtUword prev_high_row = high_row;

    GT_SHOW_AFFINE_DP_COLUMNS(
    current_idx = show_dband_value(current_idx,j-1,currentcol,low_row,high_row);
    )/*GT_SHOW_AFFINE_DP_COLUMNS*/
    gt_assert(low_row <= high_row &&
              (GtUword) (high_row - low_row + 1) <= band_width);
    colptr += (high_row - low_row + 1);
    bitmatrix[j] = colptr - low_row;
    /* below diagonal band*/
    if (j <= right_dist)
    {
      gt_assert(low_row <= prev_high_row);
      cptr = currentcol + low_row;
      first_ivalue = GT_AFFINE_SCORE_SUM(cptr->Ivalue,-gap_extension);
      GT_AFFINE_EDGE_SET(low_row,j,Affine_I,Affine_I);
    }
    nw = currentcol[low_row];
    currentcol[low_row].Rvalue = min_align_score;
    currentcol[low_row].Dvalue = min_align_score;
    currentcol[low_row].Ivalue = first_ivalue;

    if (high_row < ulen)
    {
      high_row++;
    }

    /* diagonalband */
    for (i = low_row+1; i <= high_row; i++)
    {
      GtAffineScoreTriple currententry;
      const GtUchar ca = useq[i-1];
      GtAffineAlignEditOp rmaxedge = Affine_R, dmaxedge, imaxedge = Affine_X;

      currententry.Ivalue = min_align_score;
      /* compute A_affine(i,j,R) from value in the north west */
      currententry.Rvalue = nw.Rvalue;
      if (currententry.Rvalue < nw.Dvalue)
      {
        currententry.Rvalue = nw.Dvalue;
        rmaxedge = Affine_D;
      }
      if (currententry.Rvalue < nw.Ivalue)
      {
        currententry.Rvalue = nw.Ivalue;
        rmaxedge = Affine_I;
      }
      currententry.Rvalue += GT_AFFINE_SCORE_ROW_ACCESS(ca);

      /* compute A_affine(i,j,D) from value in the north */
      cptr = currentcol + i - 1;
      score_R = GT_AFFINE_SCORE_SUM(cptr->Rvalue,start_penalty);
      score_D = GT_AFFINE_SCORE_SUM(cptr->Dvalue,-gap_extension);
      GT_AFFINE_SET_MAX2(currententry.Dvalue,dmaxedge,score_R,Affine_R,
                                                      score_D,Affine_D);
      /* compute A_affine(i,j,I) from value in the west, if available*/
      if (i <= prev_high_row)
      {
        cptr++;
        score_R = GT_AFFINE_SCORE_SUM(cptr->Rvalue,start_penalty);
        score_I = GT_AFFINE_SCORE_SUM(cptr->Ivalue,-gap_extension);
        GT_AFFINE_SET_MAX2(currententry.Ivalue,imaxedge,score_R,Affine_R,
                                                        score_I,Affine_I);
      }

      nw = currentcol[i];
      currentcol[i] = currententry;
      GT_AFFINE_EDGE_SET_ALL(i,j,rmaxedge,dmaxedge,imaxedge);
    }
    gt_assert(low_row < high_row);
    evaluated_cells += (high_row - low_row);
    if (j > right_dist)
    {
      low_row++;
    }
  }
  GT_SHOW_AFFINE_DP_COLUMNS(
  current_idx = show_dband_value(current_idx,vlen,currentcol,low_row,high_row);
  )/*GT_SHOW_AFFINE_DP_COLUMNS*/
  return currentcol[ulen].Rvalue;
}

static GtAffineScore gt_affine_diagonalband_fillDPtab_scores(
                                           GtAffineScoreTriple **dpmatrix,
                                           GtAffineScoreTriple *currentcol,
                                           int8_t gap_opening, /* > 0 */
                                           int8_t gap_extension, /* > 0 */
                                           const int8_t * const *scorematrix2D,
#ifndef NDEBUG
#else
                                           __attribute__ ((unused))
#endif
                                           GtAffineScore min_align_score,
                                           GtAffineScore smallest_score,
                                           const GtUchar *useq,
                                           GtUword ulen,
                                           const GtUchar *vseq,
                                           GtUword vlen,
                                           GtWord left_dist,
                                           GtWord right_dist)
{
  GtUword i, j, low_row = 0, high_row, evaluated_cells = 0;
  const GtAffineScore start_penalty = -(gap_opening+gap_extension);
  GtAffineScoreTriple *colptr = NULL;
  GT_SHOW_AFFINE_DP_COLUMNS(
  GtUword current_idx = 0;
  )/*GT_SHOW_AFFINE_DP_COLUMNS*/
#ifndef NDEBUG
  const GtUword band_width = (GtUword) (right_dist - left_dist + 1);
  const GtWord lendiff = (GtWord) vlen - (GtWord) ulen;
#endif

  gt_assert(currentcol != NULL &&
            gap_opening > 0 && gap_extension > 0 &&
            left_dist <= MIN(0, lendiff) && left_dist >= -(GtWord) ulen &&
            right_dist >= MAX(0, lendiff) && right_dist <= (GtWord) vlen);
  high_row = (GtUword) -left_dist;

  /* first entry */
  currentcol->Rvalue = 0;
  currentcol->Dvalue = -gap_opening;
  currentcol->Ivalue = -gap_opening;
  /* first column */
  for (i = 1; i <= high_row; i++)
  {
    currentcol[i].Rvalue = min_align_score;
    currentcol[i].Dvalue = GT_AFFINE_SCORE_SUM(currentcol[i-1].Dvalue,
                                               -gap_extension);
    currentcol[i].Ivalue = min_align_score;
  }
  evaluated_cells += (1+high_row);
  if (dpmatrix != NULL)
  {
    GtUword width = high_row - low_row + 1;
    colptr = dpmatrix[0];
    memcpy(colptr,currentcol,width * sizeof *currentcol);
    colptr += width;
  }
  GT_SHOW_AFFINE_DP_COLUMNS(
  printf("\n");
  )/*GT_SHOW_AFFINE_DP_COLUMNS*/
  /* next columns */
  for (j = 1; j <= vlen; j++)
  {
    GtAffineScoreTriple nw;
    const GtUchar cb = vseq[j-1];
    const int8_t *score_row = GT_AFFINE_SCORE_ROW_GET(cb);
    GtAffineScore first_ivalue = min_align_score, score_R, score_D, score_I;
    const GtUword prev_high_row = high_row;

    gt_assert(low_row <= high_row &&
              (GtUword) (high_row - low_row + 1) <= band_width);
    GT_SHOW_AFFINE_DP_COLUMNS(
    current_idx = show_dband_value(current_idx,j-1,currentcol,low_row,high_row);
    )/*GT_SHOW_AFFINE_DP_COLUMNS*/
    if (j <= (GtUword) right_dist)
    {
      gt_assert(low_row <= prev_high_row);
      first_ivalue = GT_AFFINE_SCORE_SUM(currentcol[low_row].Ivalue,
                                         -gap_extension);
    }
    nw = currentcol[low_row];
    currentcol[low_row].Rvalue = min_align_score;
    currentcol[low_row].Dvalue = min_align_score;
    currentcol[low_row].Ivalue = first_ivalue;

    /* do not branch in the inner loop */
    for (i = low_row+1; i <= prev_high_row; i++)
    {
      GtAffineScoreTriple currententry,
                     *cptr = currentcol + i - 1;
      const GtUchar ca = useq[i-1];

      score_R = GT_AFFINE_SCORE_SUM(cptr->Rvalue,start_penalty);
      score_D = GT_AFFINE_SCORE_SUM(cptr->Dvalue,-gap_extension);
      currententry.Dvalue = MAX(score_R,score_D);

      cptr++;
      score_R = GT_AFFINE_SCORE_SUM(cptr->Rvalue,start_penalty);
      score_I = GT_AFFINE_SCORE_SUM(cptr->Ivalue,-gap_extension);
      currententry.Ivalue = MAX(score_R,score_I);

      currententry.Rvalue = MAX3(nw.Rvalue,nw.Dvalue,nw.Ivalue) +
                            GT_AFFINE_SCORE_ROW_ACCESS(ca);
      nw = currentcol[i];
      currentcol[i] = currententry;
    }
    if (high_row < ulen)
    {
      GtAffineScoreTriple currententry,
                     *cptr = currentcol + prev_high_row;
      const GtUchar ca = useq[high_row];

      high_row++;
      score_R = GT_AFFINE_SCORE_SUM(cptr->Rvalue,start_penalty);
      score_D = GT_AFFINE_SCORE_SUM(cptr->Dvalue,-gap_extension);
      currententry.Dvalue = MAX(score_R,score_D);
      currententry.Ivalue = min_align_score;
      currententry.Rvalue = MAX3(nw.Rvalue,nw.Dvalue,nw.Ivalue) +
                            GT_AFFINE_SCORE_ROW_ACCESS(ca);
      currentcol[high_row] = currententry;
    }
    evaluated_cells += (high_row - low_row);
    if (j > (GtUword) right_dist)
    {
      low_row++;
    }
    if (dpmatrix != NULL)
    {
      GtUword width = high_row - low_row + 1;
      dpmatrix[j] = colptr - low_row;
      memcpy(colptr,currentcol + low_row,width * sizeof *currentcol);
      colptr += width;
    }
    gt_assert(low_row < high_row &&
              ((j <= (GtUword) right_dist && low_row == 0) ||
               (j > (GtUword) right_dist &&
                low_row == j - (GtUword) right_dist)) &&
                (high_row == MIN(ulen,j + (GtUword) -left_dist)));
  }
  GT_SHOW_AFFINE_DP_COLUMNS(
  current_idx = show_dband_value(current_idx,vlen,currentcol,low_row,high_row);
  )/*GT_SHOW_AFFINE_DP_COLUMNS*/
  return currentcol[ulen].Rvalue;
}

static void gt_affine_alignment_traceback_bits(GtEoplist *eoplist,
                                               const GtAffineAlignTracebits *
                                                 const *bitmatrix,
                                               const GtUchar *useq,
                                               GtUword ulen,
                                               const GtUchar *vseq,
                                               GtUword vlen)
{
  GtUword i = ulen, j = vlen;
  GtAffineAlignEditOp edge;

  gt_assert(eoplist != NULL && bitmatrix != NULL);
  edge = Affine_R;
  /* backtracing */
  while (i > 0 || j > 0)
  {
    switch (edge)
    {
      case Affine_R:
        gt_assert(i > 0 && j > 0);
        if (GT_AFFINE_EQUAL_SYMBOLS(useq[i-1],vseq[j-1]))
        {
          gt_eoplist_match_add(eoplist,1);
        } else
        {
          gt_eoplist_mismatch_add(eoplist);
        }
        edge = GT_AFFINE_EDGE_GET(i,j,Affine_R);
        i--;
        j--;
        break;
      case Affine_D:
        gt_eoplist_deletion_add(eoplist);
        edge = GT_AFFINE_EDGE_GET(i,j,Affine_D);
        gt_assert(i > 0);
        i--;
        break;
      case Affine_I:
        gt_eoplist_insertion_add(eoplist);
        edge = GT_AFFINE_EDGE_GET(i,j,Affine_I);
        gt_assert(j > 0);
        j--;
        break;
      default:
        gt_assert(false);
    }
  }
}

#ifndef NDEBUG
static GtUword low_row_on_the_fly(GtUword j,GtWord right_dist)
{
  return (j <= (GtUword) right_dist) ? 0
                                     : j - (GtUword) right_dist;
}

static GtUword high_row_on_the_fly(GtUword j,GtUword ulen,GtWord left_dist)
{
  return MIN(ulen,j + (GtUword) -left_dist);
}
#endif

static void gt_affine_alignment_traceback_scores(
                                    GtEoplist *eoplist,
                                    const GtAffineScoreTriple * const *dpmatrix,
                                    const GtUchar *useq,
                                    GtUword ulen,
                                    const GtUchar *vseq,
                                    GtUword vlen,
                                    int8_t gap_opening, /* > 0 */
                                    int8_t gap_extension, /* > 0 */
                                    __attribute__ ((unused)) GtWord left_dist,
                                    __attribute__ ((unused)) GtWord right_dist)
{
  GtUword i = ulen, j = vlen;
  GtAffineAlignEditOp edge = Affine_R;
  const GtAffineScoreTriple *previous;
  GtAffineScore maxvalue;
  const GtAffineScore start_penalty = -(gap_opening+gap_extension);

  gt_assert(eoplist != NULL && dpmatrix != NULL);
  /* backtracing */
  while (i > 0 || j > 0)
  {
    switch (edge)
    {
      case Affine_R:
        gt_assert(i > 0 && j > 0);
        if (GT_AFFINE_EQUAL_SYMBOLS(useq[i-1],vseq[j-1]))
        {
          gt_eoplist_match_add(eoplist,1);
        } else
        {
          gt_eoplist_mismatch_add(eoplist);
        }
        j--;
        i--;
        gt_assert(i >= low_row_on_the_fly(j,right_dist) &&
                  i <= high_row_on_the_fly(j,ulen,left_dist));
        previous = dpmatrix[j] + i;
        maxvalue = MAX3(previous->Rvalue,previous->Dvalue,previous->Ivalue);
        if (maxvalue > previous->Rvalue)
        {
          if (maxvalue == previous->Dvalue)
          {
            edge = Affine_D;
          } else
          {
            if (maxvalue == previous->Ivalue)
            {
              edge = Affine_I;
            }
          }
        }
        break;
      case Affine_D:
        gt_eoplist_deletion_add(eoplist);
        i--;
        gt_assert(i >= low_row_on_the_fly(j,right_dist) &&
                  i <= high_row_on_the_fly(j,ulen,left_dist));
        previous = dpmatrix[j] + i;
        if (previous->Rvalue + start_penalty >=
            previous->Dvalue - gap_extension)
        {
          edge = Affine_R;
        }
        break;
      case Affine_I:
        gt_assert(j > 0);
        gt_eoplist_insertion_add(eoplist);
        j--;
        gt_assert(i >= low_row_on_the_fly(j,right_dist) &&
                  i <= high_row_on_the_fly(j,ulen,left_dist));
        previous = dpmatrix[j] + i;
        if (previous->Rvalue + start_penalty >=
            previous->Ivalue - gap_extension)
        {
          edge = Affine_R;
        }
        break;
      default:
        gt_assert(false);
    }
  }
}

struct GtAffineDPreservoir
{
  GtUword band_width, matrix_cells, max_vlen, max_ulen;
  size_t base_type_size;
  GtAffineScoreTriple *columnspace;
  void *matrix_col_ptr, *matrix_space;
};

GtAffineDPreservoir *gt_affine_diagonalband_new(bool opt_memory,
                                                bool keepcolumns,
                                                GtUword max_ulen,
                                                GtUword max_vlen)
{
  GtAffineDPreservoir *adpr = gt_malloc(sizeof *adpr);

  adpr->band_width = 0;
  adpr->base_type_size = opt_memory ? sizeof (GtAffineAlignTracebits)
                                    : sizeof (GtAffineScoreTriple);
  adpr->max_vlen = max_vlen;
  adpr->max_ulen = max_ulen;
  if (max_ulen > 0)
  {
    adpr->columnspace = gt_malloc((max_ulen+1) * sizeof *adpr->columnspace);
  } else
  {
    adpr->columnspace = NULL;
  }
  if (keepcolumns && max_vlen > 0)
  {
    adpr->matrix_col_ptr = gt_malloc((max_vlen+1) * sizeof (void *));
  } else
  {
    adpr->matrix_col_ptr = NULL;
  }
  adpr->matrix_cells = 0;
  adpr->matrix_space = NULL;
  return adpr;
}

static bool gt_affine_opt_memory(const GtAffineDPreservoir *adpr)
{
  gt_assert(adpr != NULL);
  return adpr->base_type_size == sizeof (GtAffineAlignTracebits) ? true
                                                                 : false;
}

static void *gt_affine_diagonalband_col_ptr(GtAffineDPreservoir *adpr,
                                            GtUword vlen)
{
  gt_assert(adpr != NULL);

  if (vlen > adpr->max_vlen)
  {
    adpr->max_vlen = MAX(vlen,adpr->max_vlen * 1.2 + 128);
    adpr->matrix_col_ptr = gt_realloc(adpr->matrix_col_ptr,
                                      (adpr->max_vlen+1) * sizeof (void *));
  } else
  {
    gt_assert(adpr->matrix_col_ptr != NULL);
  }
  return adpr->matrix_col_ptr;
}

static void *gt_affine_diagonalband_column_space(GtAffineDPreservoir *adpr,
                                                 GtUword ulen)
{
  gt_assert(adpr != NULL);

  if (ulen > adpr->max_ulen)
  {
    adpr->max_ulen = MAX(ulen,adpr->max_ulen * 1.2 + 128);
    adpr->columnspace = gt_realloc(adpr->columnspace,
                                   (adpr->max_ulen+1) *
                                   sizeof *adpr->columnspace);
  } else
  {
    gt_assert(adpr->columnspace != NULL);
  }
  return (void *) adpr->columnspace;
}

static void *gt_affine_diagonalband_matrix_space(GtAffineDPreservoir *adpr,
                                                 GtUword band_width,
                                                 GtUword vlen)
{
  gt_assert(adpr != NULL);
  if (band_width * (vlen+1) >= adpr->matrix_cells)
  {
    adpr->matrix_cells = MAX(band_width * (vlen+1),
                             adpr->matrix_cells * 1.2 + 1024);
    adpr->matrix_space = gt_realloc(adpr->matrix_space,
                                    adpr->matrix_cells * adpr->base_type_size);
  }
  gt_assert(adpr->matrix_space != NULL);
  return (void *) adpr->matrix_space;
}

void gt_affine_diagonalband_delete(GtAffineDPreservoir *adpr)
{
  if (adpr != NULL)
  {
    gt_free(adpr->matrix_col_ptr);
    gt_free(adpr->matrix_space);
    gt_free(adpr->columnspace);
    gt_free(adpr);
  }
}

/* calculate alignment within diagonalband specified by left_dist and
   right_dist.  space and running time is O(bandwidth * vlen) */

static GtAffineScore gt_affine_diagonalband_align(
                                          bool keep_columns,
                                          GtAffineDPreservoir *adpr,
                                          int8_t gap_opening, /* > 0 */
                                          int8_t gap_extension, /* > 0 */
                                          const int8_t * const *scorematrix2D,
                                          GtAffineScore min_align_score,
                                          GtAffineScore smallest_score,
                                          const GtUchar *useq,
                                          GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vlen,
                                          GtWord left_dist,
                                          GtWord right_dist)
{
  GtAffineScore lastcolumnRvalue;
  const GtUword band_width = (GtUword) (right_dist - left_dist + 1);
  GtAffineScoreTriple *columnspace
    = (GtAffineScoreTriple *) gt_affine_diagonalband_column_space(adpr,ulen);

  if (keep_columns)
  {
    if (gt_affine_opt_memory(adpr))
    {
      GtAffineAlignTracebits **bitmatrix
        = (GtAffineAlignTracebits **) gt_affine_diagonalband_col_ptr(adpr,vlen);

      bitmatrix[0]
        = (GtAffineAlignTracebits *)
          gt_affine_diagonalband_matrix_space(adpr,band_width,vlen);
      lastcolumnRvalue
        = gt_affine_diagonalband_fillDPtab_bits(bitmatrix,
                                                columnspace,
                                                gap_opening,
                                                gap_extension,
                                                scorematrix2D,
                                                min_align_score,
                                                smallest_score,
                                                useq,
                                                ulen,
                                                vseq,
                                                vlen,
                                                left_dist,
                                                right_dist);
    } else
    {
      GtAffineScoreTriple **dpmatrix
        = (GtAffineScoreTriple **) gt_affine_diagonalband_col_ptr(adpr,vlen);

      dpmatrix[0] = (GtAffineScoreTriple *)
                    gt_affine_diagonalband_matrix_space(adpr,band_width,vlen);
      lastcolumnRvalue
        = gt_affine_diagonalband_fillDPtab_scores(dpmatrix,
                                                  columnspace,
                                                  gap_opening,
                                                  gap_extension,
                                                  scorematrix2D,
                                                  min_align_score,
                                                  smallest_score,
                                                  useq,
                                                  ulen,
                                                  vseq,
                                                  vlen,
                                                  left_dist,
                                                  right_dist);
    }
  } else
  {
    lastcolumnRvalue
      = gt_affine_diagonalband_fillDPtab_scores(NULL,
                                                columnspace,
                                                gap_opening,
                                                gap_extension,
                                                scorematrix2D,
                                                min_align_score,
                                                smallest_score,
                                                useq,
                                                ulen,
                                                vseq,
                                                vlen,
                                                left_dist,
                                                right_dist);
  }
  return lastcolumnRvalue;
}

static void gt_affine_diagonalband_traceback(GtAffineDPreservoir *adpr,
                                             GtEoplist *eoplist,
                                             const GtUchar *useq,
                                             const GtUword ulen,
                                             const GtUchar *vseq,
                                             const GtUword vlen,
                                             int8_t gap_opening, /* > 0 */
                                             int8_t gap_extension, /* > 0 */
                                             GtWord left_dist,
                                             GtWord right_dist)
{
  gt_eoplist_reset(eoplist);
  if (gt_affine_opt_memory(adpr))
  {
    const GtAffineAlignTracebits *const *bitmatrix
      = (const GtAffineAlignTracebits * const*)
        gt_affine_diagonalband_col_ptr(adpr,vlen);
    gt_affine_alignment_traceback_bits(eoplist,bitmatrix,useq,ulen,vseq,vlen);
  } else
  {
    const GtAffineScoreTriple *const * dpmatrix
      = (const GtAffineScoreTriple * const*)
        gt_affine_diagonalband_col_ptr(adpr,vlen);
    gt_affine_alignment_traceback_scores(eoplist,
                                         dpmatrix,
                                         useq,
                                         ulen,
                                         vseq,
                                         vlen,
                                         gap_opening,
                                         gap_extension,
                                         left_dist,
                                         right_dist);
  }
  gt_eoplist_reverse_end(eoplist,0);
}

GtUword gt_affine_iter_diagonalband_align(GtEoplist *eoplist,
                                       GtAffineDPreservoir *adpr,
                                       int8_t gap_opening, /* > 0 */
                                       int8_t gap_extension, /* > 0 */
                                       const int8_t * const *scorematrix2D,
                                       int8_t smallest_score,
#ifndef NDEBUG
#else
                                       __attribute__ ((unused))
#endif
                                       const GtUchar *useq,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vlen,
                                       bool no_score_run,
                                       GtUword expected_score)
{
  GtUword iteration;
#ifndef NDEBUG
  const GtWord lendiff = (GtWord) vlen - (GtWord) ulen;
#endif
  GtUword band_width = 1 + (GtUword) ((vlen > ulen) ? vlen - ulen
                                                    : ulen - vlen),
          previous_dpscore = GT_UWORD_MAX;
  const GtWord min_align_score = (ulen + vlen) * (GtAffineScore) smallest_score;

  gt_assert(smallest_score < 0 && min_align_score > INT_MIN/2);
  for (iteration = 0; /* Nothing */; iteration++)
  {
    GtAffineScore dpscore;
    GtWord left_dist = - (GtWord) band_width,
           right_dist = (GtWord) band_width;

    gt_assert(left_dist <= MIN(0,lendiff) && right_dist >= MAX(0,lendiff));
    if (left_dist < -(GtWord) ulen)
    {
      left_dist = -(GtWord) ulen;
    }
    if (right_dist > (GtWord) vlen)
    {
      right_dist = (GtWord) vlen;
    }
    dpscore = gt_affine_diagonalband_align(eoplist,
                                           adpr,
                                           gap_opening,
                                           gap_extension,
                                           scorematrix2D,
                                           (GtAffineScore) min_align_score,
                                           smallest_score,
                                           useq,
                                           ulen,
                                           vseq,
                                           vlen,
                                           left_dist,
                                           right_dist);
    if (expected_score == 0 || dpscore >= (GtAffineScore) expected_score ||
        (no_score_run && previous_dpscore == dpscore))
    {
      if (eoplist != NULL)
      {
        gt_affine_diagonalband_traceback(adpr,
                                         eoplist,
                                         useq,
                                         ulen,
                                         vseq,
                                         vlen,
                                         gap_opening, /* > 0 */
                                         gap_extension, /* > 0 */
                                         left_dist,
                                         right_dist);
      }
      return (GtUword) dpscore;
    }
    previous_dpscore = dpscore;
    if (band_width < 4)
    {
      band_width *= 2;
    } else
    {
      if (band_width < 20)
      {
        band_width = (band_width * 3)/2;
      } else
      {
        band_width = (band_width * 5)/4;
      }
    }
    gt_assert (left_dist != -(GtWord) ulen || right_dist != vlen);
  }
  gt_assert(false);
  return 0;
}
