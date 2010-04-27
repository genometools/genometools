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

#ifndef CHECKENCCHAR_H
#define CHECKENCCHAR_H

#ifdef SKDEBUG
#define GT_CHECKENCCHAR(CC,ENCSEQ,POS,READMODE)\
        {\
          GtUchar cctmp = gt_encseq_get_encoded_char(ENCSEQ,POS, \
                                                            READMODE);\
          if ((CC) != cctmp)\
          {\
            printf("file %s, line %d: pos = %lu:cc = %u != %u = ccreal\n",\
                   __FILE__,__LINE__,\
                   (unsigned long) (POS),\
                   (unsigned int) (CC),\
                   (unsigned int) cctmp);\
            exit(GT_EXIT_PROGRAMMING_ERROR);\
          }\
        }
#else
#define GT_CHECKENCCHAR(CC,ENCSEQ,POS,READMODE)
#endif

#endif
