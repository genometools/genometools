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
#include "alphadef.h"
#include "filelength-def.h"
#include "encseq-def.h"

#define FILEBUFFERSIZE 65536

#define SARR_ESQTAB ((unsigned int) 1)
#define SARR_SUFTAB ((unsigned int) 1 << 1)
#define SARR_LCPTAB ((unsigned int) 1 << 2)
#define SARR_BWTTAB ((unsigned int) 1 << 3)
#define SARR_DESTAB ((unsigned int) 1 << 4)

#define SARR_ALLTAB (SARR_ESQTAB |\
                     SARR_SUFTAB |\
                     SARR_LCPTAB |\
                     SARR_BWTTAB |\
                     SARR_DESTAB)

#define DECLAREBufferedfiletype(TYPE)\
        typedef struct\
        {\
          unsigned int nextfree,\
                       nextread;\
          TYPE bufspace[FILEBUFFERSIZE];\
          FILE *fp;\
        } TYPE ## Bufferedfile

DECLAREBufferedfiletype(Seqpos);

DECLAREBufferedfiletype(Uchar);

DECLAREBufferedfiletype(Largelcpvalue);

#define DECLAREREADFUNCTION(TYPE)\
        static int readnext ## TYPE ## fromstream(TYPE *val,\
                                                  TYPE ## Bufferedfile *buf,\
                                                  Env *env)\
        {\
          if (buf->nextread >= buf->nextfree)\
          {\
            buf->nextfree = (unsigned int) fread(buf->bufspace,\
                                                 sizeof (TYPE),\
                                                 (size_t) FILEBUFFERSIZE,\
                                                 buf->fp);\
            if (ferror(buf->fp))\
            {\
              env_error_set(env,"error when trying to read next %s",#TYPE);\
              return -2;\
            }\
            buf->nextread = 0;\
            if (buf->nextfree == 0)\
            {\
              return 0;\
            }\
          }\
          *val = buf->bufspace[buf->nextread++];\
          return 1;\
        }

typedef struct
{
  unsigned long numofdbsequences; /* XXX: move to encoded sequence */
  StrArray *filenametab;
  Filelengthvalues *filelengthtab;
  unsigned int prefixlength;
  DefinedSeqpos numoflargelcpvalues;
  Encodedsequence *encseq;
  DefinedSeqpos longest;
  Specialcharinfo specialcharinfo; /* XXX: move to encoded sequence */
  Alphabet *alpha;                 /* XXX: move to encoded sequence */
  Readmode readmode; /* relevant when reading the encoded sequence */
  /* either with mapped input */
  const Seqpos *suftab;
  const Uchar *lcptab;
  const Largelcpvalue *llvtab;
  const Uchar *bwttab;
  /* or with streams */
  SeqposBufferedfile suftabstream;
  UcharBufferedfile bwttabstream,
                    lcptabstream;
  LargelcpvalueBufferedfile llvtabstream;
  const char *destab;
  unsigned long destablength;
} Suffixarray;

#endif
