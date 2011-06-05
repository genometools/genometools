/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef ESA_SEQREAD_H
#define ESA_SEQREAD_H
#include <stdbool.h>
#include "core/str.h"
#include "core/error.h"
#include "core/encseq.h"

#include "sarr-def.h"
#include "lcpoverflow.h"

typedef enum
{
  SEQ_scan,
  SEQ_mappedboth,
  SEQ_suftabfrommemory
} Sequentialaccesstype;

#ifdef INLINEDSequentialsuffixarrayreader

typedef struct
{
  Suffixarray *suffixarray;
  unsigned long nonspecials,
         numberofsuffixes,
         nextsuftabindex, /* for SEQ_mappedboth | SEQ_suftabfrommemory */
         nextlcptabindex, /* for SEQ_mappedboth */
         largelcpindex;   /* SEQ_mappedboth */
} Sequentialsuffixarrayreader;

#define NEXTSEQUENTIALLCPTABVALUE(VALUE,SSAR)\
        if ((SSAR)->nextlcptabindex >= (SSAR)->numberofsuffixes)\
        {\
          gt_error_set(err,"missing lcpvalue value");\
          haserr = true;\
          break;\
        } else\
        {\
          GtUchar tmpsmalllcpvalue\
            = (SSAR)->suffixarray->lcptab[(SSAR)->nextlcptabindex++];\
          if (tmpsmalllcpvalue < LCPOVERFLOW)\
          {\
            VALUE = (unsigned long) tmpsmalllcpvalue;\
          } else\
          {\
            VALUE = (SSAR)->suffixarray->llvtab[\
                    (SSAR)->largelcpindex++].value;\
          }\
        }

#define NEXTSEQUENTIALSUFTABVALUE(VALUE,SSAR)\
        VALUE = (SSAR)->suffixarray->suftab[(SSAR)->nextsuftabindex++]

#else

#include "esa-lcpval.h"

struct Sequentialsuffixarrayreader
{
  Suffixarray *suffixarray;
  unsigned long nonspecials,
         numberofsuffixes,
         nextsuftabindex, /* for SEQ_mappedboth | SEQ_suftabfrommemory */
         nextlcptabindex, /* for SEQ_mappedboth */
         largelcpindex;   /* SEQ_mappedboth */
  Sequentialaccesstype seqactype;
  Lcpvalueiterator *lvi;
  const ESASuffixptr *suftab;
  const GtEncseq *encseq;
  GtReadmode readmode;
};

typedef struct Sequentialsuffixarrayreader Sequentialsuffixarrayreader;

/* The following four function can only can be used if
   INLINEDSequentialsuffixarrayreader is not set */

Sequentialsuffixarrayreader *gt_newSequentialsuffixarrayreaderfromRAM(
                                        const GtEncseq *encseq,
                                        GtReadmode readmode);

void gt_updateSequentialsuffixarrayreaderfromRAM(
                    Sequentialsuffixarrayreader *ssar,
                    const ESASuffixptr *suftab,
                    bool firstpage,
                    unsigned long numberofsuffixes);

int gt_nextSequentiallcpvalue(unsigned long *currentlcp,
                           Sequentialsuffixarrayreader *ssar,
                           GtError *err);

int gt_nextSequentialsuftabvalue(unsigned long *currentsuffix,
                                 Sequentialsuffixarrayreader *ssar);

#define NEXTSEQUENTIALLCPTABVALUE(VALUE,SSAR)\
        {\
          GtUchar tmpsmalllcpvalue;\
          if ((SSAR)->seqactype == SEQ_scan)\
          {\
            int retval = readnextGtUcharfromstream(&tmpsmalllcpvalue,\
                                        &(SSAR)->suffixarray->lcptabstream);\
            if (retval > 0)\
            {\
              if (tmpsmalllcpvalue == LCPOVERFLOW)\
              {\
                Largelcpvalue tmpexception;\
                retval = readnextLargelcpvaluefromstream(&tmpexception,\
                                        &(SSAR)->suffixarray->llvtabstream);\
                if (retval == 0)\
                {\
                  gt_error_set(err,"file %s: line %d: unexpected end "\
                                   "of file when reading llvtab",\
                                   __FILE__,__LINE__);\
                  haserr = true;\
                  break;\
                }\
                VALUE = tmpexception.value;\
              } else\
              {\
                VALUE = (unsigned long) tmpsmalllcpvalue;\
              }\
            } else\
            {\
              break;\
            }\
          } else\
          {\
            if ((SSAR)->seqactype == SEQ_mappedboth)\
            {\
              if ((SSAR)->nextlcptabindex < (SSAR)->numberofsuffixes)\
              {\
                tmpsmalllcpvalue\
                  = (SSAR)->suffixarray->lcptab[(SSAR)->nextlcptabindex++];\
                if (tmpsmalllcpvalue == LCPOVERFLOW)\
                {\
                  gt_assert((SSAR)->suffixarray->llvtab[(SSAR)->largelcpindex]\
                             .position == (SSAR)->nextlcptabindex-1);\
                  VALUE = (SSAR)->suffixarray->llvtab\
                                      [(SSAR)->largelcpindex++].value;\
                } else\
                {\
                  VALUE = (unsigned long) tmpsmalllcpvalue;\
                }\
              } else\
              {\
                break;\
              }\
            } else\
            {\
              if ((SSAR)->nextlcptabindex < (SSAR)->numberofsuffixes)\
              {\
                VALUE = gt_nextLcpvalueiterator((SSAR)->lvi,\
                                                true,\
                                                (SSAR)->suftab,\
                                                (SSAR)->numberofsuffixes);\
                (SSAR)->nextlcptabindex++;\
              } else\
              {\
                break;\
              }\
            }\
          }\
        }

#define NEXTSEQUENTIALSUFTABVALUE(VALUE,SSAR)\
        switch ((SSAR)->seqactype)\
        {\
          case SEQ_scan:\
            {\
              GtUlongBufferedfile *buf = &(SSAR)->suffixarray->suftabstream;\
              if (buf->nextread >= buf->nextfree)\
              {\
                buf->nextfree = (unsigned int) fread(buf->bufferedfilespace,\
                                                     sizeof (GtUlong),\
                                                     (size_t) FILEBUFFERSIZE,\
                                                     buf->fp);\
                if (ferror(buf->fp))\
                {\
                  gt_error_set(err,"error when trying to read next GtUlong");\
                  haserr = true;\
                } else\
                {\
                  buf->nextread = 0;\
                  if (buf->nextfree == 0)\
                  {\
                    gt_error_set(err,"Missing value in suftab");\
                    haserr = true;\
                    break;\
                  }\
                }\
              }\
              VALUE = buf->bufferedfilespace[buf->nextread++];\
            }\
            break;\
          case SEQ_mappedboth:\
            VALUE = ESASUFFIXPTRGET((SSAR)->suffixarray->suftab,\
                                    (SSAR)->nextsuftabindex++);\
            break;\
          case SEQ_suftabfrommemory:\
            VALUE = ESASUFFIXPTRGET((SSAR)->suftab,\
                                    (SSAR)->nextsuftabindex++);\
          break;\
        }

#endif /* INLINEDSequentialsuffixarrayreader */

Sequentialsuffixarrayreader *gt_newSequentialsuffixarrayreaderfromfile(
                                        const char *indexname,
                                        unsigned int demand,
                                        Sequentialaccesstype seqactype,
                                        GtError *err);

void gt_freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar);

const GtEncseq *gt_encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar);

GtReadmode gt_readmodeSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar);

unsigned long gt_Sequentialsuffixarrayreader_nonspecials(
                          const Sequentialsuffixarrayreader *ssar);

const ESASuffixptr *gt_suftabSequentialsuffixarrayreader(
                        const Sequentialsuffixarrayreader *ssar);

const Suffixarray *gt_suffixarraySequentialsuffixarrayreader(
              const Sequentialsuffixarrayreader *ssar);

unsigned long gt_Sequentialsuffixarrayreader_totallength(
              const Sequentialsuffixarrayreader *ssar);

#endif
