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

#include "esa-lcpval.h"
#include "sarr-def.h"
#include "lcpoverflow.h"

typedef enum
{
  SEQ_scan,
  SEQ_mappedboth,
  SEQ_suftabfrommemory
} Sequentialaccesstype;

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

#define NEXTSEQUENTIALSUFTABVALUE_SEQ_scan_generic(SUFTABVALUE,SSAR,TYPE)\
        {\
          GtBufferedfile_ ## TYPE *buf\
            = &(SSAR)->suffixarray->suftabstream_ ## TYPE;\
          if (buf->nextread >= buf->nextfree)\
          {\
            buf->nextfree\
              = (unsigned int) fread(buf->bufferedfilespace,\
                                     sizeof (*buf->bufferedfilespace),\
                                     (size_t) FILEBUFFERSIZE,\
                                     buf->fp);\
            if (ferror(buf->fp))\
            {\
              gt_error_set(err,"error when trying to read next %s",#TYPE);\
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
          SUFTABVALUE = (unsigned long)buf->bufferedfilespace[buf->nextread++];\
        }

#ifdef _LP64
#define NEXTSEQUENTIALSUFTABVALUE_SEQ_scan(SUFTABVALUE,SSAR)\
        if ((SSAR)->suffixarray->suftabstream_GtUlong.fp != NULL)\
        {\
          NEXTSEQUENTIALSUFTABVALUE_SEQ_scan_generic(SUFTABVALUE,SSAR,GtUlong);\
        } else\
        {\
          NEXTSEQUENTIALSUFTABVALUE_SEQ_scan_generic(SUFTABVALUE,SSAR,\
                                                     uint32_t);\
        }
#else
#define NEXTSEQUENTIALSUFTABVALUE_SEQ_scan(SUFTABVALUE,SSAR)\
        NEXTSEQUENTIALSUFTABVALUE_SEQ_scan_generic(SUFTABVALUE,SSAR,GtUlong)
#endif

#define NEXTSEQUENTIALSUFTABVALUE(SUFTABVALUE,SSAR)\
        switch ((SSAR)->seqactype)\
        {\
          case SEQ_scan:\
            NEXTSEQUENTIALSUFTABVALUE_SEQ_scan(SUFTABVALUE,SSAR);\
            break;\
          case SEQ_mappedboth:\
            SUFTABVALUE = ESASUFFIXPTRGET((SSAR)->suffixarray->suftab,\
                                          (SSAR)->nextsuftabindex++);\
            break;\
          case SEQ_suftabfrommemory:\
            SUFTABVALUE = ESASUFFIXPTRGET((SSAR)->suftab,\
                                          (SSAR)->nextsuftabindex++);\
          break;\
        }

#define NEXTSEQUENTIALLCPTABVALUE(LCPVALUE,SSAR)\
        {\
          GtUchar tmpsmalllcpvalue;\
          if ((SSAR)->seqactype == SEQ_scan)\
          {\
            int retval = gt_readnextfromstream_GtUchar(&tmpsmalllcpvalue,\
                                        &(SSAR)->suffixarray->lcptabstream);\
            if (retval > 0)\
            {\
              if (tmpsmalllcpvalue < LCPOVERFLOW)\
              {\
                LCPVALUE = (unsigned long) tmpsmalllcpvalue;\
              } else\
              {\
                Largelcpvalue tmpexception;\
                retval = gt_readnextfromstream_Largelcpvalue(&tmpexception,\
                                        &(SSAR)->suffixarray->llvtabstream);\
                if (retval == 0)\
                {\
                  gt_error_set(err,"file %s: line %d: unexpected end "\
                                   "of file when reading llvtab",\
                                   __FILE__,__LINE__);\
                  haserr = true;\
                  break;\
                }\
                LCPVALUE = tmpexception.value;\
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
                if (tmpsmalllcpvalue < LCPOVERFLOW)\
                {\
                  LCPVALUE = (unsigned long) tmpsmalllcpvalue;\
                } else\
                {\
                  gt_assert((SSAR)->suffixarray->llvtab[(SSAR)->largelcpindex]\
                             .position == (SSAR)->nextlcptabindex-1);\
                  LCPVALUE = (SSAR)->suffixarray->llvtab\
                                      [(SSAR)->largelcpindex++].value;\
                }\
              } else\
              {\
                break;\
              }\
            } else\
            {\
              if ((SSAR)->nextlcptabindex < (SSAR)->numberofsuffixes)\
              {\
                LCPVALUE = gt_nextLcpvalueiterator((SSAR)->lvi,\
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

#define NEXTSEQUENTIALLCPTABVALUEWITHLAST(LCPVALUE,LASTSUFTABVALUE,SSAR)\
        {\
          GtUchar tmpsmalllcpvalue;\
          if ((SSAR)->seqactype == SEQ_scan)\
          {\
            int retval = gt_readnextfromstream_GtUchar(&tmpsmalllcpvalue,\
                                        &(SSAR)->suffixarray->lcptabstream);\
            if (retval > 0)\
            {\
              if (tmpsmalllcpvalue < LCPOVERFLOW)\
              {\
                LCPVALUE = (unsigned long) tmpsmalllcpvalue;\
              } else\
              {\
                Largelcpvalue tmpexception;\
                retval = gt_readnextfromstream_Largelcpvalue(&tmpexception,\
                                        &(SSAR)->suffixarray->llvtabstream);\
                if (retval == 0)\
                {\
                  gt_error_set(err,"file %s: line %d: unexpected end "\
                                   "of file when reading llvtab",\
                                   __FILE__,__LINE__);\
                  haserr = true;\
                  break;\
                }\
                LCPVALUE = tmpexception.value;\
              }\
            } else\
            {\
              NEXTSEQUENTIALSUFTABVALUE_SEQ_scan(LASTSUFTABVALUE,SSAR);\
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
                if (tmpsmalllcpvalue < LCPOVERFLOW)\
                {\
                  LCPVALUE = (unsigned long) tmpsmalllcpvalue;\
                } else\
                {\
                  gt_assert((SSAR)->suffixarray->llvtab[(SSAR)->largelcpindex]\
                             .position == (SSAR)->nextlcptabindex-1);\
                  LCPVALUE = (SSAR)->suffixarray->llvtab\
                                      [(SSAR)->largelcpindex++].value;\
                }\
              } else\
              {\
                LASTSUFTABVALUE = ESASUFFIXPTRGET((SSAR)->suffixarray->suftab,\
                                                  (SSAR)->nextsuftabindex++);\
                break;\
              }\
            } else\
            {\
              if ((SSAR)->nextlcptabindex < (SSAR)->numberofsuffixes)\
              {\
                LCPVALUE = gt_nextLcpvalueiterator((SSAR)->lvi,\
                                                   true,\
                                                   (SSAR)->suftab,\
                                                   (SSAR)->numberofsuffixes);\
                (SSAR)->nextlcptabindex++;\
              } else\
              {\
                LASTSUFTABVALUE = ESASUFFIXPTRGET((SSAR)->suftab,\
                                                  (SSAR)->nextsuftabindex++);\
                break;\
              }\
            }\
          }\
        }

Sequentialsuffixarrayreader *gt_newSequentialsuffixarrayreaderfromfile(
                                        const char *indexname,
                                        unsigned int demand,
                                        Sequentialaccesstype seqactype,
                                        GtLogger *logger,
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

unsigned int gt_Sequentialsuffixarrayreader_prefixlength(
              const Sequentialsuffixarrayreader *ssar);

GtBcktab *gt_Sequentialsuffixarrayreader_bcktab(
              const Sequentialsuffixarrayreader *ssar);

#endif
