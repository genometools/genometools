/*
  Copyright (c) 2016 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2016 Joerg Winkler <j.winkler@posteo.de>
  Copyright (c) 2016 Center for Bioinformatics, University of Hamburg

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

#ifndef DECLARE_READFUNC_H
#define DECLARE_READFUNC_H

#include <stdio.h>
#include "core/unused_api.h"

#define FILEBUFFERSIZE 4096

#define DECLAREBufferedfiletype(TYPE)\
        typedef struct\
        {\
          unsigned int nextfree,\
                       nextread;\
          TYPE *bufferedfilespace;\
          FILE *fp;\
        } GtBufferedfile_ ## TYPE

#define DECLAREREADFUNCTION(TYPE)\
        GT_UNUSED static int gt_readnextfromstream_ ## TYPE (TYPE *val,\
                                                  GtBufferedfile_ ## TYPE *buf)\
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

#endif
