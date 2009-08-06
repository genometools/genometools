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

#ifndef SARR_DEF_H
#define SARR_DEF_H

#include <stdio.h>
#include "seqpos-def.h"
#include "intcode-def.h"
#include "encseq-def.h"
#include "bcktab.h"

#define FILEBUFFERSIZE 4096

#define SARR_ESQTAB 1U
#define SARR_SUFTAB (1U << 1)
#define SARR_LCPTAB (1U << 2)
#define SARR_BWTTAB (1U << 3)
#define SARR_DESTAB (1U << 4)
#define SARR_SDSTAB (1U << 5)
#define SARR_BCKTAB (1U << 6)
#define SARR_SSPTAB (1U << 7)

#define SARR_ALLTAB (SARR_ESQTAB |\
                     SARR_SUFTAB |\
                     SARR_LCPTAB |\
                     SARR_BWTTAB |\
                     SARR_DESTAB |\
                     SARR_SSPTAB)

#define DECLAREBufferedfiletype(TYPE)\
        typedef struct\
        {\
          unsigned int nextfree,\
                       nextread;\
          TYPE *bufferedfilespace;\
          FILE *fp;\
        } TYPE ## Bufferedfile

DECLAREBufferedfiletype(Seqpos);

DECLAREBufferedfiletype(GtUchar);

DECLAREBufferedfiletype(Largelcpvalue);

#define DECLAREREADFUNCTION(TYPE)\
        static int readnext ## TYPE ## fromstream(TYPE *val,\
                                                  TYPE ## Bufferedfile *buf)\
        {\
          if (buf->nextread >= buf->nextfree)\
          {\
            buf->nextfree = (unsigned int) fread(buf->bufferedfilespace,\
                                                 sizeof (TYPE),\
                                                 (size_t) FILEBUFFERSIZE,\
                                                 buf->fp);\
            if (ferror(buf->fp))\
            {\
              fprintf(stderr,"error when trying to read next %s",#TYPE);\
              exit(GT_EXIT_PROGRAMMING_ERROR);\
            }\
            buf->nextread = 0;\
            if (buf->nextfree == 0)\
            {\
              return 0;\
            }\
          }\
          *val = buf->bufferedfilespace[buf->nextread++];\
          return 1;\
        }

typedef struct
{
  Encodedsequence *encseq;
  DefinedSeqpos numoflargelcpvalues; /* only in esa-map.c */
  DefinedSeqpos longest; /* for BWT */
  Readmode readmode; /* relevant when reading the encoded sequence */
  /* either with mapped input */
  const Seqpos *suftab;
  const GtUchar *lcptab;
  const Largelcpvalue *llvtab;
  const GtUchar *bwttab;
  unsigned int prefixlength;
  Bcktab *bcktab;
  /* or with streams */
  SeqposBufferedfile suftabstream;
  GtUcharBufferedfile bwttabstream,
                    lcptabstream;
  LargelcpvalueBufferedfile llvtabstream;
} Suffixarray;

#endif
