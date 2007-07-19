/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SARR_DEF_H
#define SARR_DEF_H

#include <stdio.h>
#include "types.h"
#include "alphadef.h"
#include "encseq-def.h"

#define FILEBUFFERSIZE 65536

#define DECLAREBufferedfiletype(TYPE)\
        typedef struct\
        {\
          uint32_t nextfree,\
                   nextread;\
          TYPE bufspace[FILEBUFFERSIZE];\
          FILE *fp;\
        } TYPE ## Bufferedfile

DECLAREBufferedfiletype(Seqpos);

DECLAREBufferedfiletype(Uchar);

DECLAREBufferedfiletype(Largelcpvalue);

#define DECLAREREADFUNCTION(TYPE)\
        static int readnext ## TYPE ## fromstream(TYPE *val,\
                                                  TYPE ## Bufferedfile *buf)\
        {\
          if(buf->nextread >= buf->nextfree)\
          {\
            buf->nextfree = (uint32_t) fread(buf->bufspace,\
                                             sizeof(TYPE),\
                                             (size_t) FILEBUFFERSIZE,\
                                             buf->fp);\
            if(ferror(buf->fp))\
            {\
              env_error_set(env,"error when trying to read next %s",#TYPE);\
              return -2;\
            }\
            buf->nextread = 0;\
            if(buf->nextfree == 0)\
            {\
              return 0;\
            }\
          }\
          *val = buf->bufspace[buf->nextread++];\
          return 1;\
        }

typedef struct
{
  uint32_t numofdbsequences;
  StrArray *filenametab;
  Filelengthvalues *filelengthtab;
  uint32_t prefixlength;
  DefinedSeqpos numoflargelcpvalues;
  Encodedsequence *encseq;
  DefinedSeqpos longest;
  Specialcharinfo specialcharinfo;
  Alphabet *alpha;
  bool mappedinput;
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
} Suffixarray;

#endif

