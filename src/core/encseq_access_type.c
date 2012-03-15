/*
  Copyright (c) 2007-2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include <inttypes.h>
#include <string.h>
#include "core/alphabet.h"
#include "core/assert_api.h"
#include "core/bitpackarray.h"
#include "core/chardef.h"
#include "core/encseq.h"
#include "core/encseq_access_type.h"
#include "core/filelengthvalues.h"
#include "core/intbits.h"
#include "core/types_api.h"

typedef struct
{
  GtEncseqAccessType sat;
  char *name;
} GtWrittenAccessType;

static GtWrittenAccessType wpa[] = {
  {GT_ACCESS_TYPE_DIRECTACCESS, "direct"},
  {GT_ACCESS_TYPE_BYTECOMPRESS, "bytecompress"},
  {GT_ACCESS_TYPE_EQUALLENGTH, "eqlen"},
  {GT_ACCESS_TYPE_BITACCESS, "bit"},
  {GT_ACCESS_TYPE_UCHARTABLES, "uchar"},
  {GT_ACCESS_TYPE_USHORTTABLES, "ushort"},
  {GT_ACCESS_TYPE_UINT32TABLES, "uint32"}
};

const char* gt_encseq_access_type_list(void)
{
  return "direct, bytecompress, eqlen, bit, uchar, ushort, uint32";
}

const char* gt_encseq_access_type_str(GtEncseqAccessType at)
{
  gt_assert((int) at < (int) GT_ACCESS_TYPE_UNDEFINED);
  return wpa[at].name;
}

GtEncseqAccessType gt_encseq_access_type_get(const char *str)
{
  size_t i;

  for (i=0; i<sizeof (wpa)/sizeof (wpa[0]); i++)
  {
    gt_assert(i == 0 || wpa[i-1].sat < wpa[i].sat);
    if (strcmp(str,wpa[i].name) == 0)
    {
      return wpa[i].sat;
    }
  }
  return GT_ACCESS_TYPE_UNDEFINED;
}

bool gt_encseq_access_type_isviautables(GtEncseqAccessType sat)
{
  gt_assert(sat != GT_ACCESS_TYPE_UNDEFINED);
  return (sat >= GT_ACCESS_TYPE_UCHARTABLES) ? true : false;
}

#define CHECKANDUPDATE(SAT,IDX)\
        tmp = gt_encseq_determine_size(SAT, totallength,\
                                       numofsequences,\
                                       numofdbfiles,\
                                       lengthofdbfilenames,\
                                       wildcardrangestab[IDX],\
                                       numofchars,\
                                       0,\
                                       lengthofalphadef);\
        if (tmp < cmin)\
        {\
          cmin = tmp;\
          cret = SAT;\
          *specialranges = specialrangestab[IDX];\
          *wildcardranges = wildcardrangestab[IDX];\
        }

static GtEncseqAccessType determinesmallestrep(
                                  unsigned long *specialranges,
                                  unsigned long *wildcardranges,
                                  const Definedunsignedlong *equallength,
                                  unsigned long totallength,
                                  unsigned long numofsequences,
                                  unsigned long numofdbfiles,
                                  unsigned long lengthofdbfilenames,
                                  const unsigned long *specialrangestab,
                                  const unsigned long *wildcardrangestab,
                                  unsigned int numofchars,
                                  unsigned long lengthofalphadef)
{
  GtEncseqAccessType cret;
  uint64_t tmp, cmin;

  cret = GT_ACCESS_TYPE_BITACCESS;
  cmin = gt_encseq_determine_size(cret, totallength,numofsequences,
                                  numofdbfiles, lengthofdbfilenames,
                                  wildcardrangestab[0], numofchars, 0,
                                  lengthofalphadef);
  *specialranges = specialrangestab[0];
  *wildcardranges = wildcardrangestab[0];
  if (equallength != NULL && equallength->defined)
  {
    cret = GT_ACCESS_TYPE_EQUALLENGTH;
  } else
  {
    CHECKANDUPDATE(GT_ACCESS_TYPE_UCHARTABLES, 0);
    CHECKANDUPDATE(GT_ACCESS_TYPE_USHORTTABLES, 1);
    CHECKANDUPDATE(GT_ACCESS_TYPE_UINT32TABLES, 2);
  }
  return cret;
}

int gt_encseq_access_type_determine(unsigned long *specialranges,
                                    unsigned long *wildcardranges,
                                    unsigned long totallength,
                                    unsigned long numofsequences,
                                    unsigned long numofdbfiles,
                                    unsigned long lengthofalphadef,
                                    unsigned long lengthofdbfilenames,
                                    const unsigned long *specialrangestab,
                                    const unsigned long *wildcardrangestab,
                                    const Definedunsignedlong *equallength,
                                    unsigned int numofchars,
                                    const char *str_sat,
                                    GtError *err)
{
  GtEncseqAccessType sat = GT_ACCESS_TYPE_UNDEFINED;
  bool haserr = false;

  *specialranges = specialrangestab[0];
  *wildcardranges = wildcardrangestab[0];
  if (str_sat == NULL)
  {
    if (numofchars == GT_DNAALPHASIZE)
    {
      sat = determinesmallestrep(specialranges,wildcardranges,
                                 equallength,totallength,numofsequences,
                                 numofdbfiles,lengthofdbfilenames,
                                 specialrangestab,wildcardrangestab,numofchars,
                                 lengthofalphadef);
    } else
    {
      sat = GT_ACCESS_TYPE_BYTECOMPRESS;
    }
  } else
  {
    sat = gt_encseq_access_type_get(str_sat);
    if (numofchars == GT_DNAALPHASIZE)
    {
      switch (sat)
      {
        case GT_ACCESS_TYPE_UCHARTABLES:
          *specialranges = specialrangestab[0];
          *wildcardranges = wildcardrangestab[0];
          break;
        case GT_ACCESS_TYPE_USHORTTABLES:
          *specialranges = specialrangestab[1];
          *wildcardranges = wildcardrangestab[1];
           break;
        case GT_ACCESS_TYPE_UINT32TABLES:
          *specialranges = specialrangestab[2];
          *wildcardranges = wildcardrangestab[2];
          break;
        case GT_ACCESS_TYPE_DIRECTACCESS:
        case GT_ACCESS_TYPE_BITACCESS:
          break;
        case GT_ACCESS_TYPE_EQUALLENGTH:
          if (equallength == NULL || !equallength->defined) {
            gt_error_set(err,"illegal argument \"%s\" to option -sat: "
                             "%s is only possible for DNA sequences, if "
                             "all sequences are of equal length and no "
                             "sequence contains a wildcard",str_sat,str_sat);
            haserr = true;
          }
          break;
        case GT_ACCESS_TYPE_BYTECOMPRESS:
          gt_error_set(err,"illegal argument \"%s\" to option -sat: "
                           "cannot use bytecompress on DNA sequences",
                           str_sat);
          haserr = true;
          break;
        default:
          gt_assert(sat == GT_ACCESS_TYPE_UNDEFINED);
          gt_error_set(err,"illegal argument \"%s\" to option -sat: "
                           "must be one of the following keywords: %s",
                           str_sat,gt_encseq_access_type_list());
          haserr = true;
          break;
      }
    } else
    {
      if (sat != GT_ACCESS_TYPE_BYTECOMPRESS &&
          sat != GT_ACCESS_TYPE_DIRECTACCESS)
      {
        gt_error_set(err,"illegal argument \"%s\" to option -sat: "
                        "as the sequence is not DNA, you can choose %s or %s",
                        str_sat,
                        gt_encseq_access_type_str(GT_ACCESS_TYPE_BYTECOMPRESS),
                        gt_encseq_access_type_str(GT_ACCESS_TYPE_DIRECTACCESS));
        haserr = true;
      }
    }
  }
  return haserr ? -1 : (int) sat;
}
