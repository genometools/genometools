/*
  Copyright (c) 2007-2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2013 Center for Bioinformatics, University of Hamburg

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

struct Sequentialsuffixarrayreader
{
  Suffixarray *suffixarray;
  GtUword nonspecials,
         numberofsuffixes,
         nextsuftabindex, /* for !scanfile */
         nextlcptabindex, /* for !scanfile */
         largelcpindex;   /* for !scanfile */
  const ESASuffixptr *suftab;
  const GtEncseq *encseq;
  bool scanfile;
  void *extrainfo;
  GtReadmode readmode;
};

typedef struct Sequentialsuffixarrayreader Sequentialsuffixarrayreader;

int gt_nextSequentiallcpvalue(GtUword *currentlcp,
                           Sequentialsuffixarrayreader *ssar,
                           GtError *err);

int gt_nextSequentialsuftabvalue(GtUword *currentsuffix,
                                 Sequentialsuffixarrayreader *ssar);

#define SSAR_NEXTSEQUENTIALSUFTABVALUE_SEQ_scan_generic(SUFTABVALUE,SSAR,TYPE)\
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
          SUFTABVALUE = (GtUword)buf->bufferedfilespace[buf->nextread++];\
        }

#if defined (_LP64) || defined (_WIN64)
#define SSAR_NEXTSEQUENTIALSUFTABVALUE_SEQ_scan(SUFTABVALUE,SSAR)\
        if ((SSAR)->suffixarray->suftabstream_GtUword.fp != NULL)\
        {\
          SSAR_NEXTSEQUENTIALSUFTABVALUE_SEQ_scan_generic(SUFTABVALUE,SSAR,\
                                                          GtUword);\
        } else\
        {\
          SSAR_NEXTSEQUENTIALSUFTABVALUE_SEQ_scan_generic(SUFTABVALUE,SSAR,\
                                                          uint32_t);\
        }
#else
#define SSAR_NEXTSEQUENTIALSUFTABVALUE_SEQ_scan(SUFTABVALUE,SSAR)\
        SSAR_NEXTSEQUENTIALSUFTABVALUE_SEQ_scan_generic(SUFTABVALUE,SSAR,\
                                                        GtUword)
#endif

#define SSAR_NEXTSEQUENTIALSUFTABVALUE(SUFTABVALUE,SSAR)\
        if ((SSAR)->scanfile)\
        {\
          SSAR_NEXTSEQUENTIALSUFTABVALUE_SEQ_scan(SUFTABVALUE,SSAR);\
        } else\
        {\
          SUFTABVALUE = ESASUFFIXPTRGET((SSAR)->suffixarray->suftab,\
                                        (SSAR)->nextsuftabindex++);\
        }

#define SSAR_NEXTSEQUENTIALLCPTABVALUE(LCPVALUE,SSAR)\
        {\
          GtUchar tmpsmalllcpvalue;\
          if ((SSAR)->scanfile)\
          {\
            int retval = gt_readnextfromstream_GtUchar(&tmpsmalllcpvalue,\
                                        &(SSAR)->suffixarray->lcptabstream);\
            if (retval > 0)\
            {\
              if (tmpsmalllcpvalue < (GtUchar) LCPOVERFLOW)\
              {\
                LCPVALUE = (GtUword) tmpsmalllcpvalue;\
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
            if ((SSAR)->nextlcptabindex < (SSAR)->numberofsuffixes)\
            {\
              tmpsmalllcpvalue\
                = (SSAR)->suffixarray->lcptab[(SSAR)->nextlcptabindex++];\
              if (tmpsmalllcpvalue < (GtUchar) LCPOVERFLOW)\
              {\
                LCPVALUE = (GtUword) tmpsmalllcpvalue;\
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
          }\
        }

#define SSAR_NEXTSEQUENTIALLCPTABVALUEWITHLAST(LCPVALUE,LASTSUFTABVALUE,SSAR)\
        {\
          GtUchar tmpsmalllcpvalue;\
          if ((SSAR)->scanfile)\
          {\
            int retval = gt_readnextfromstream_GtUchar(&tmpsmalllcpvalue,\
                                        &(SSAR)->suffixarray->lcptabstream);\
            if (retval > 0)\
            {\
              if (tmpsmalllcpvalue < (GtUchar) LCPOVERFLOW)\
              {\
                LCPVALUE = (GtUword) tmpsmalllcpvalue;\
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
              SSAR_NEXTSEQUENTIALSUFTABVALUE_SEQ_scan(LASTSUFTABVALUE,SSAR);\
              break;\
            }\
          } else\
          {\
            if ((SSAR)->nextlcptabindex < (SSAR)->numberofsuffixes)\
            {\
              tmpsmalllcpvalue\
                = (SSAR)->suffixarray->lcptab[(SSAR)->nextlcptabindex++];\
              if (tmpsmalllcpvalue < (GtUchar) LCPOVERFLOW)\
              {\
                LCPVALUE = (GtUword) tmpsmalllcpvalue;\
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
          }\
        }

Sequentialsuffixarrayreader *gt_newSequentialsuffixarrayreaderfromfile(
                                        const char *indexname,
                                        unsigned int demand,
                                        bool scanfile,
                                        GtLogger *logger,
                                        GtError *err);

void gt_freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar);

const GtEncseq *gt_encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar);

GtReadmode gt_readmodeSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar);

GtUword gt_Sequentialsuffixarrayreader_nonspecials(
                          const Sequentialsuffixarrayreader *ssar);

const ESASuffixptr *gt_suftabSequentialsuffixarrayreader(
                        const Sequentialsuffixarrayreader *ssar);

const Suffixarray *gt_suffixarraySequentialsuffixarrayreader(
              const Sequentialsuffixarrayreader *ssar);

GtUword gt_Sequentialsuffixarrayreader_totallength(
              const Sequentialsuffixarrayreader *ssar);

GtUword gt_Sequentialsuffixarrayreader_maxbranchdepth(
              const Sequentialsuffixarrayreader *ssar);

unsigned int gt_Sequentialsuffixarrayreader_prefixlength(
              const Sequentialsuffixarrayreader *ssar);

#endif
