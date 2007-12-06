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
#include "libgtcore/str.h"
#include "libgtcore/error.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "sarr-def.h"

typedef enum
{
  SEQ_mappedboth,
  SEQ_scan,
  SEQ_suftabfrommemory
} Sequentialaccesstype;

#ifdef INLINEDSequentialsuffixarrayreader

typedef struct
{
  Suffixarray *suffixarray;
  Seqpos numberofsuffixes,
         nextsuftabindex, /* for SEQ_mappedboth | SEQ_suftabfrommemory */
         nextlcptabindex, /* for SEQ_mappedboth */
         largelcpindex;   /* SEQ_mappedboth */
} Sequentialsuffixarrayreader;

#define NEXTSEQUENTIALLCPTABVALUE(VALUE,SSAR)\
        if ((SSAR)->nextlcptabindex >= (SSAR)->numberofsuffixes)\
        {\
          break;\
        } else\
        {\
          tmpsmalllcpvalue\
            = (SSAR)->suffixarray->lcptab[(SSAR)->nextlcptabindex++];\
          if (tmpsmalllcpvalue == (Uchar) UCHAR_MAX)\
          {\
            VALUE = (SSAR)->suffixarray->llvtab[\
                    (SSAR)->largelcpindex++].value;\
          } else\
          {\
            VALUE = (Seqpos) tmpsmalllcpvalue;\
          }\
        }

#define NEXTSEQUENTIALSUFTABVALUE(VALUE,SSAR)\
        VALUE = (SSAR)->suffixarray->suftab[(SSAR)->nextsuftabindex++]

#else

typedef struct Sequentialsuffixarrayreader Sequentialsuffixarrayreader;

/* The following only can be used for this case */

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromRAM(
                                        const Encodedsequence *encseq,
                                        Readmode readmode);

/* The following can only be used for this case */

void updateSequentialsuffixarrayreaderfromRAM(
                    Sequentialsuffixarrayreader *ssar,
                    const Seqpos *suftab,
                    bool firstpage,
                    Seqpos numberofsuffixes);

int nextSequentiallcpvalue(Seqpos *currentlcp,
                           Sequentialsuffixarrayreader *ssar,
                           Error *err);

int nextSequentialsuftabvalue(Seqpos *currentsuffix,
                              Sequentialsuffixarrayreader *ssar,
                              Error *err);

#endif

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromfile(
                                        const Str *indexname,
                                        unsigned int demand,
                                        Sequentialaccesstype seqactype,
                                        Error *err);

void freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar);

const Encodedsequence *encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *sarr);

Readmode readmodeSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *sarr);

const Alphabet *alphabetSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *sarr);

unsigned long numofdbsequencesSequentialsuffixarrayreader(
                    const Sequentialsuffixarrayreader *sarr);

unsigned long destablengthSequentialsuffixarrayreader(
              const Sequentialsuffixarrayreader *sarr);

const char *destabSequentialsuffixarrayreader(
              const Sequentialsuffixarrayreader *sarr);

#endif
