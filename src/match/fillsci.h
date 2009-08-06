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

#ifndef FILLSCI_H
#define FILLSCI_H

#include <stdbool.h>
#include "core/alphabet.h"
#include "core/str_api.h"
#include "core/filelengthvalues.h"
#include "core/str_array.h"
#include "core/arraydef.h"
#include "core/error.h"
#include "core/symboldef.h"
#include "seqpos-def.h"
#include "verbose-def.h"

int fasta2sequencekeyvalues(
        const GtStr *indexname,
        Seqpos *totallength,
        Specialcharinfo *specialcharinfo,
        unsigned int forcetable,
        Seqpos *specialrangestab,
        const GtStrArray *filenametab,
        Filelengthvalues **filelengthtab,
        const GtAlphabet *alpha,
        bool plainformat,
        bool outdestab,
        bool outsdstab,
        bool outkystab,
        unsigned long *characterdistribution,
        bool outssptab,
        ArraySeqpos *sequenceseppos,
        Verboseinfo *verboseinfo,
        GtError *err);

void sequence2specialcharinfo(Specialcharinfo *specialcharinfo,
                              const GtUchar *seq,
                              const Seqpos len,
                              Verboseinfo *verboseinfo);

#endif
