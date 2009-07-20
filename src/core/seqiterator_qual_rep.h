/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQITERATOR_QUAL_REP_H
#define SEQITERATOR_QUAL_REP_H

#include "core/error.h"
#include "core/filelengthvalues.h"
#include "core/file.h"
#include "core/seqiterator_qual.h"
#include "core/str_array.h"

#define GT_SEQIT_QUAL_INBUFSIZE  8192

struct GtSeqIteratorQualClass {
  size_t        size;
  int           (*next)(GtSeqIteratorQual*,
                        const GtUchar**,
                        const GtUchar**,
                        unsigned long*,
                        char**,
                        GtError*);
  void          (*set_symbolmap)(GtSeqIteratorQual*,
                                 const GtUchar*);
  void           (*set_chardisttab)(GtSeqIteratorQual*,
                                    unsigned long*);
  uint64_t       (*get_lastspeciallength)(const GtSeqIteratorQual*);
  const unsigned long long*
                 (*get_current_counter)(GtSeqIteratorQual*,
                                        unsigned long long);
  void          (*free)(GtSeqIteratorQual*);
};

typedef struct GtSeqIteratorQualMembers GtSeqIteratorQualMembers;

struct GtSeqIteratorQual {
  const GtSeqIteratorQualClass *c_class;
  GtSeqIteratorQualMembers *pvt;
};

struct GtSeqIteratorQualMembers {
  unsigned int filenum;
  uint64_t linenum;
  Filelengthvalues *filelengthtab;
  bool complete,
       use_ungetchar;
  GtStr *sequencebuffer,
        *descbuffer,
        *qualsbuffer;
  GtFile *curfile;
  unsigned long *chardisttab,
                currentfillpos,
                currentinpos,
                curline,
                reference_count;
  uint64_t lastspeciallength;
  unsigned long long maxread,
                     currentread;
  const GtStrArray *filenametab;
  unsigned char ungetchar,
                inbuf[GT_SEQIT_QUAL_INBUFSIZE];
  const GtUchar *symbolmap;
};

GtSeqIteratorQual*
gt_seqiterator_qual_create(const GtSeqIteratorQualClass*);

void*
gt_seqiterator_qual_cast(const GtSeqIteratorQualClass*, GtSeqIteratorQual*);

#endif
