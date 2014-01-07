/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/mathsupport.h"
#include "core/types_api.h"
#include "core/chardef.h"
#include "core/unused_api.h"
#include "core/intbits.h"
#include "core/ma_api.h"
#include "core/readmode.h"
#include "match/sfx-lwcheck.h"

#ifdef _WIN32
#define GEN_ESASUFFIXPTRGET(TAB,IDX) ((GtUword*) TAB)[IDX]
#else
#define GEN_ESASUFFIXPTRGET(TAB,IDX) (unitsize == (size_t) 4\
                                       ? (GtUword) ((unsigned int *) TAB)[IDX] \
                                       : (GtUword) ((GtUword *) TAB)[IDX])
#endif

static GtUword gt_check_for_range_occurrence(const void *suftab,
                                             GT_UNUSED size_t unitsize,
                                             GtUword suffix,
                                             GtUword start,
                                             GtUword end)
{
  GtUword idx;

  for (idx = start; idx <= end; idx++)
  {
    GtUword position = GEN_ESASUFFIXPTRGET(suftab, idx);

    if (suffix == position)
    {
      return idx;
    }
  }
  return ULONG_MAX;
}

typedef struct
{
  GtUword start, end;
  GtUchar firstchar;
} GtRangewithchar;

/* The following funktion implements the linear time algorithm of
   @INPROCEEDINGS{BUR:KAER:2003,
   author = {Burkhardt, S. and K{\"a}rkk{\"a}inen, J.},
   title = {{Fast Lightweight Suffix Array Construction and Checking}},
   booktitle = {{Proceedings of the 14th Annual Symposium on Combinatorial
                 Pattern Matching (CPM)}},
   year = {2003},
   editor = {{Baeza-Yates, R. and Ch{\'a}vez, E. and Crochemore, M.}},
   volume = {2676},
   series = {LNCS},
   pages = {200-210},
   publisher = {Springer-Verlag}
   }
   to check the following suffix-order condition of the sorted suffix array:
   For all characters c, if SA[i, j] contains the suffixes starting
   with charcter c, then SA[i]+1, SA[i+1]+1, \ldots, SA[j]+1 occur
   in SA in this order (but not consecutively in general).
   The running time of the algorithm is independent of the alphabet size.
   The main problem is that it requires random access to the sequence
   which slows it down.
*/

static void gt_suftab_bk_suffixorder(LW_accesschar accesschar,
                                     const void *encseq,
                                     GtReadmode readmode,
                                     GtUword totallength,
                                     unsigned int numofchars,
                                     const void  *suftab,
                                     GT_UNUSED size_t unitsize,
                                     const GtRangewithchar *rangestore,
                                     unsigned int numofranges)
{
  unsigned int rangeidx;
  GtUword idx, *nexttab = gt_calloc((size_t) numofchars, sizeof (*nexttab));

  for (rangeidx = 0; rangeidx < numofranges; rangeidx++)
  {
    nexttab[rangestore[rangeidx].firstchar] = rangestore[rangeidx].start;
  }
  for (idx = 0; idx < totallength; idx++)
  {
    GtUword position = GEN_ESASUFFIXPTRGET(suftab, idx);

    if (position > 0)
    {
      GtUchar cc = accesschar(encseq,position-1,readmode);
      if (ISNOTSPECIAL(cc))
      {
        GtUword checkpos;

        checkpos = GEN_ESASUFFIXPTRGET(suftab, nexttab[(int) cc]) + 1;
        if (checkpos != position)
        {
          fprintf(stderr, "idx="GT_WU", checkpos="GT_WU", position="GT_WU"\n",
                          idx, checkpos, position);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        nexttab[(int) cc]++;
      }
    }
  }
  gt_free(nexttab);
}

/* The following function checks the suffix-order condition described
   above using an O(\sigma n) algorithm where \sigma is the alphabet size
   and n is the length of the suffix array. The algorithm does not
   access the sequence and performs \sigma linear scans of the suffix array.
   This makes it faster than the previous method. It is probably possible
   to perform the check in one linear scan.
*/

static void gt_suftab_sk_suffixorder(GtUword totallength,
                                     unsigned int numofchars,
                                     const void *suftab,
                                     size_t unitsize,
                                     const GtRangewithchar *rangestore,
                                     unsigned int numofranges)
{
  unsigned int rangeidx;
  GtUword numofcomparisons = 0;
  double ratio;

  for (rangeidx = 0; rangeidx < numofranges; rangeidx++)
  {
    GtUword idx, start = 0;

    for (idx = rangestore[rangeidx].start; idx <= rangestore[rangeidx].end;
         idx++)
    {
      GtUword position = GEN_ESASUFFIXPTRGET(suftab, idx);

      if (position + 1 <= totallength)
      {
        GtUword found = gt_check_for_range_occurrence(suftab,
                                                      unitsize,
                                                      position + 1,
                                                      start,
                                                      totallength);
        if (found == ULONG_MAX)
        {
          fprintf(stderr, "Cannot find position+1="GT_WU" in range ["GT_WU", "
                  GT_WU"]\n", position+1, start, totallength);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        numofcomparisons += found - start + 1;
        start = found + 1;
      }
    }
  }
  ratio = (double) numofcomparisons/totallength;
  if (gt_double_compare(ratio, (double) numofchars) > 0)
  {
    fprintf(stderr, "gt_double_compare(%.2f, %u) = %d > 0 not exected\n",
                    ratio, numofchars,
                    gt_double_compare(ratio, (double) numofchars));
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

void gt_suftab_lightweightcheck(LW_accesschar accesschar,
                                LW_charcount charcount,
                                const void *encseq,
                                GtReadmode readmode,
                                GtUword totallength,
                                unsigned int numofchars,
                                const void *suftab,
                                size_t unitsize,
                                GtLogger *logger)
{
  GtUword idx, countbitsset = 0, previouspos = 0,
                firstspecial = totallength, rangestart = 0;
  unsigned int charidx, rangeidx = 0, numofranges;
  GtBitsequence *startposoccurs;
  GtUchar previouscc = 0;
  GtRangewithchar *rangestore;
  const bool skcheck = true;

#ifdef _WIN32
  gt_assert(unitsize == (size_t) 4);
#else
  gt_assert(unitsize == (size_t) 4 || unitsize == (size_t) 8);
#endif
  GT_INITBITTAB(startposoccurs, totallength+1);
  rangestore = gt_malloc(sizeof(*rangestore) * numofchars);
  for (idx = 0; idx < totallength; idx++)
  {
    GtUword position = GEN_ESASUFFIXPTRGET(suftab, idx);
    GtUchar cc;

    if (GT_ISIBITSET(startposoccurs, position))
    {
      fprintf(stderr, "ERROR: suffix with startpos "GT_WU" already occurs\n",
              GEN_ESASUFFIXPTRGET(suftab, idx));
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    GT_SETIBIT(startposoccurs, position);
    countbitsset++;
    cc = accesschar(encseq,position,readmode);
    if (idx > 0)
    {
      if (ISSPECIAL(cc))
      {
        if (firstspecial == totallength)
        {
          firstspecial = idx;
          gt_assert(rangeidx < numofchars);
          rangestore[rangeidx].start = rangestart;
          rangestore[rangeidx].end = idx-1;
          rangestore[rangeidx++].firstchar = previouscc;
        }
        if (ISSPECIAL(previouscc))
        {
          if (previouspos > position)
          {
            fprintf(stderr, "incorrect order: " GT_WU " = " GT_WU "=SPECIAL > "
                    "SPECIAL="GT_WU"  = "GT_WU"\n",
                    idx-1, position, previouspos, idx);
            exit(GT_EXIT_PROGRAMMING_ERROR);
          }
        }
      } else
      {
        if (ISSPECIAL(previouscc))
        {
          fprintf(stderr, "incorrect order: "GT_WU"="GT_WU"=SPECIAL > "
                  "%u="GT_WU"="GT_WU"\n", idx-1, position, (unsigned int) cc,
                  previouspos, idx);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        } else
        {
          if (previouscc > cc)
          {
            fprintf(stderr, "incorrect order: "GT_WU" = "GT_WU"=%u > "
                    "%u="GT_WU"="GT_WU"\n", idx-1, position, (unsigned int)
                    previouscc, (unsigned int) cc, previouspos, idx);
            exit(GT_EXIT_PROGRAMMING_ERROR);
          } else
          {
            if (previouscc < cc)
            {
              gt_assert(rangeidx < numofchars);
              rangestore[rangeidx].start = rangestart;
              rangestore[rangeidx].end = idx-1;
              rangestore[rangeidx++].firstchar = previouscc;
              rangestart = idx;
            }
          }
        }
      }
    } else
    {
      if (ISSPECIAL(cc))
      {
        firstspecial = 0;
      }
    }
    previouscc = cc;
    previouspos = position;
  }
  if (countbitsset != totallength)
  {
    fprintf(stderr, "ERROR: only "GT_WU" of "GT_WU" suffixes occur\n",
            countbitsset, totallength);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_free(startposoccurs);
  if (firstspecial == totallength)
  {
    gt_assert(firstspecial > 0 && rangeidx < numofchars);
    rangestore[rangeidx].start = rangestart;
    rangestore[rangeidx].end = firstspecial-1;
    rangestore[rangeidx++].firstchar = previouscc;
  }
  numofranges = rangeidx;
  for (charidx = 0, rangeidx = 0; charidx < numofchars; charidx++)
  {
    GtUword count;
    GtUchar thisidx = GT_ISDIRCOMPLEMENT(readmode)
                        ? (GtUchar) GT_COMPLEMENTBASE(charidx)
                        : (GtUchar) charidx;

    count = charcount(encseq,thisidx);
    if (count != 0)
    {
      gt_assert(rangeidx < numofranges &&
                rangestore[rangeidx].firstchar == (GtUchar) charidx);
      gt_assert(rangestore[rangeidx].end - rangestore[rangeidx].start + 1
                == count);
      rangeidx++;
    }
  }
  gt_logger_log(logger, "suftab-check, first phase done");
  if (skcheck)
  {
    gt_suftab_sk_suffixorder(totallength,
                             numofchars,
                             suftab,
                             unitsize,
                             rangestore,
                             numofranges);
    gt_logger_log(logger, "suftab-check, second phase (sk-method) done");
  } else
  {
    gt_suftab_bk_suffixorder(accesschar,
                             encseq,
                             readmode,
                             totallength,
                             numofchars,
                             suftab,
                             unitsize,
                             rangestore,
                             numofranges);
    gt_logger_log(logger, "suftab-check, second phase (bk-method) done");
  }
  gt_free(rangestore);
}
