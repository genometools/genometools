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

#ifndef INL_ENCSEQ_H
#define INL_ENCSEQ_H

typedef struct
{
  Uchar *plainseq;
  Seqpos totallength;
  unsigned long numofdbsequences;

  NEWMAPSPEC(encseq->specialcharinfoptr,Specialcharinfo,1UL);
  numofchars = getencseqAlphabetnumofchars(encseq);
  NEWMAPSPEC(encseq->characterdistribution,Unsignedlong,
             (unsigned long) numofchars);
  switch (encseq->sat)
  bool hasownmemory, mappedfile, hasspecialcharacters;
} Encodedsequence;

typedef struct
{
  bool moveforward, exhausted;
  const Encodedsequence *encseq;
  Seqpos pos,
         lengthofspecialrange;
} Specialrangeiterator;

typedef struct
{
  Readmode readmode;
} Encodedsequencescanstate;

#define getencseqtotallength(ENCSEQ) ((ENCSEQ)->totallength)

#define getencseqnumofdbsequences(ENCSEQ) ((ENCSEQ)->numofdbsequences)

#define MAKECOMPL(CC)\
        (ISSPECIAL(CC) ? (CC) : (Uchar) 3 - (CC))

#define getencodedchar(ENCSEQ,POS,RM)\
        (((RM) == Forwardmode)\
          ? (ENCSEQ)->plainseq[POS]\
          : (((RM) == Reversemode)\
            ? (ENCSEQ)->plainseq[REVERSEPOS((ENCSEQ)->totallength,POS)]\
            : (((RM) == Complementmode) \
              ? MAKECOMPL((ENCSEQ)->plainseq[POS])\
              : (MAKECOMPL((ENCSEQ)->plainseq[\
                           REVERSEPOS((ENCSEQ)->totallength,POS)])\
              )\
            )\
          )\
        )

#define getencodedcharnospecial(ENCSEQ,POS,RM)\
        getencodedchar(ENCSEQ,POS,RM)

#define sequentialgetencodedchar(ENCSEQ,ENCSEQSTATE,POS,READMODE)\
        getencodedchar(ENCSEQ,POS,READMODE)

#endif
