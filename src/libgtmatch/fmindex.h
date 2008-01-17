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

#ifndef FMINDEX_H
#define FMINDEX_H

#include "libgtcore/arraydef.h"
#include "seqpos-def.h"
#include "alphadef.h"
#include "encseq-def.h"
#include "fmi-bwtbound.h"

/*
  TO DO:
  - free index if in the merging process all suffixes are inserted.

  - implement option -check to determine the exact size of the index
    for the different methods.

  - can we guarantee that in the BWT we only access positions smaller
    than
    firstignorespecial = fm->bwtlength - (specialcharacters + 1);
    Then we do not have to store so many special positions.

  - generate special versions of fmqhits for the following cases:
       - Viadirectaccess
       - specialcharacters == 0

    This means that the macros with the corresponding character
    accesses are replaced at compile time.

  - also for single input index generate bwt file if option
    -fmout is used.

  - can we, for a given alphabet, first determine all bounds like in
    binsplitinterval?

  - trenne uniquesub in uniquesub und uniquesub-esa auf.

  - integriere codierte BWTtab in den fmindex

  - Programmierrichtlinien: safe cast
                            single thread f"ahig
                            dynamische Allokation f"ur Dateinamen
                            flexible kombination von Optionen
*/

#define FMASCIIFILESUFFIX ".fma"
#define FMDATAFILESUFFIX  ".fmd"

#define MARKPOSTABLELENGTH(BWTLENGTH,MARKDIST)\
        (1 + ((BWTLENGTH) - 1) / (MARKDIST))

#define TFREQSIZE(MAPSIZE)\
        ((MAPSIZE) + 1)

#define BFREQSIZE(MAPSIZE,NOFBLOCKS)\
        ((MAPSIZE) * (NOFBLOCKS))

#define SUPERBFREQSIZE(MAPSIZE,NOFSUPERBLOCKS)\
        ((MAPSIZE) * (NOFSUPERBLOCKS))

typedef int(*FMprocessqhit)(void *,Seqpos,Seqpos);

typedef struct
{
  Encodedsequence *bwtformatching;
  Uchar *bfreq;            /* bfreq[c][i] = #c in block i */
  Seqpos bwtlength,        /* also totallength + 1 */
         *tfreq,           /* tfreq[c] = #characters < c in text */
         *superbfreq,      /* superbfreq[c][i] = #c in all superblocks */
                           /* which are previous to superblock i */
         *markpostable,    /* sampling of entries from suffix array */
         longestsuffixpos,
         negatebsizeones,
         negatesuperbsizeones,
         markdistminus1;   /* markdist - 1 */
  Specialcharinfo specialcharinfo;
  ArrayPairBwtidx specpos; /* positions of special characters */
  Alphabet *alphabet;
  void *mappedptr; /* NULL or pointer to the mapped space block */
  unsigned int mapsize,      /* copy of alphabet.mapsize, used for searching */
           bsize,            /* size of block */
           bsizehalve,       /* DIV2(fm->bsize) */
           superbsize,       /* size of superblock */
           log2bsize,        /* log_{2}(bsize) */
           log2superbsize,   /* log_{2}(superbsize) */
           log2superbsizeminuslog2bsize, /* log{2}(superbsize)-log{2}(bsize) */
           log2markdist,     /* log_{2}(markdist) */
           suffixlength;     /* len of suffix for which buckets are computed*/
  unsigned long sizeofindex; /* size of the fmindex in bytes */
  Seqpos nofblocks,          /* number of blocks (bwtlength/bsize + 1) */
         nofsuperblocks,     /* number of superblocks (bwtlength/superbsize+2)*/
         markdist,           /* multiple of entry num stored in suffix array */
         numofcodes;         /* number of entries in boundaries */
  Bwtbound *boundarray;      /* corresponding boundaries */
} Fmindex;

typedef struct
{
  Fmindex *fmptr;
  bool storeindexpos;
} Fmindexwithoptions;

#endif
