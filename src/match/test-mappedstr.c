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

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/error.h"
#include "core/unused_api.h"
#include "core/encseq.h"
#include "sfx-nextchar.h"
#include "kmer2string.h"
#include "sfx-mappedstr.h"

static GtCodetype qgram2codefillspecial(unsigned int numofchars,
                                      unsigned int kmersize,
                                      const GtEncseq *encseq,
                                      GtReadmode readmode,
                                      unsigned long startpos,
                                      unsigned long totallength)
{
  GtCodetype integercode;
  unsigned long pos;
  bool foundspecial;
  GtUchar cc;

  if (startpos >= totallength)
  {
    integercode = (GtCodetype) (numofchars - 1);
    foundspecial = true;
  } else
  {
    /* for testing */
    cc = gt_encseq_get_encoded_char(encseq,startpos,readmode);
    if (ISSPECIAL(cc))
    {
      integercode = (GtCodetype) (numofchars - 1);
      foundspecial = true;
    } else
    {
      integercode = (GtCodetype) cc;
      foundspecial = false;
    }
  }
  for (pos = startpos + 1; pos < startpos + kmersize; pos++)
  {
    if (foundspecial)
    {
      ADDNEXTCHAR(integercode,numofchars-1,numofchars);
    } else
    {
      if (pos >= totallength)
      {
        ADDNEXTCHAR(integercode,numofchars-1,numofchars);
        foundspecial = true;
      } else
      {
        /* for testing */
        cc = gt_encseq_get_encoded_char(encseq,pos,readmode);
        if (ISSPECIAL(cc))
        {
          ADDNEXTCHAR(integercode,numofchars-1,numofchars);
          foundspecial = true;
        } else
        {
          ADDNEXTCHAR(integercode,cc,numofchars);
        }
      }
    }
  }
  return integercode;
}

GT_DECLAREARRAYSTRUCT(GtCodetype);

static void outkmeroccurrence(void *processinfo,
                              const GtKmercode *kmercode)
{
  GtArrayGtCodetype *codelist = (GtArrayGtCodetype *) processinfo;

  GT_STOREINARRAY(codelist,GtCodetype,1024,kmercode->code);
}

/*
   The function to collect the code from a stream of fasta files
   can only produce the sequence of code in forward mode.
   Hence we compute the corresponding sequence also in GT_READMODE_FORWARD.
   Thus we restrict the call for gt_verifymappedstr to the case where
   the suffix array is in readmode = GT_READMODE_FORWARD.
*/

static void collectkmercode(GtArrayGtCodetype *codelist,
                            const GtEncseq *encseq,
                            unsigned int kmersize,
                            unsigned int numofchars,
                            unsigned long stringtotallength)
{
  unsigned long offset;
  GtCodetype code;

  for (offset=0; offset<=stringtotallength; offset++)
  {
    code = qgram2codefillspecial(numofchars,
                                 kmersize,
                                 encseq,
                                 GT_READMODE_FORWARD,
                                 offset,
                                 stringtotallength);
    GT_STOREINARRAY(codelist,GtCodetype,1024,code);
  }
}

static int comparecodelists(const GtArrayGtCodetype *codeliststream,
                            const GtArrayGtCodetype *codeliststring,
                            unsigned int kmersize,
                            unsigned int numofchars,
                            const char *characters,
                            GtError *err)
{
  unsigned long i;
  char buffer1[64+1], buffer2[64+1];

  gt_error_check(err);
  if (codeliststream->nextfreeGtCodetype != codeliststring->nextfreeGtCodetype)
  {
    gt_error_set(err,"length codeliststream= %lu != %lu =length codeliststring",
                  (unsigned long) codeliststream->nextfreeGtCodetype,
                  (unsigned long) codeliststring->nextfreeGtCodetype);
    return -1;
  }
  for (i=0; i<codeliststream->nextfreeGtCodetype; i++)
  {
    if (codeliststream->spaceGtCodetype[i] !=
          codeliststring->spaceGtCodetype[i])
    {
      gt_fromkmercode2string(buffer1,
                          codeliststream->spaceGtCodetype[i],
                          numofchars,
                          kmersize,
                          characters);
      gt_fromkmercode2string(buffer2,
                          codeliststring->spaceGtCodetype[i],
                          numofchars,
                          kmersize,
                          characters);
      gt_error_set(err,"codeliststream[%lu] = " FormatGtCodetype " != "
                    FormatGtCodetype " = codeliststring[%lu]\n%s != %s",
                    i,
                    codeliststream->spaceGtCodetype[i],
                    codeliststring->spaceGtCodetype[i],
                    i,
                    buffer1,
                    buffer2);
      return -1;
    }
  }
  return 0;
}

static int getfastastreamkmers(const GtStrArray *filenametab,
                               unsigned int numofchars,
                               unsigned int kmersize,
                               const GtUchar *symbolmap,
                               bool plainformat,
                               GtArrayGtCodetype *codeliststream,
                               GtError *err)
{
  GtKmercodeiterator *kmercodeiterator;
  const GtKmercode *kmercodeptr;
  bool haserr = false;

  kmercodeiterator = gt_kmercodeiterator_filetab_new(
                                filenametab,
                                numofchars,
                                kmersize,
                                symbolmap,
                                plainformat,
                                err);
  if (!gt_kmercodeiterator_inputexhausted(kmercodeiterator))
  {
    while (!haserr)
    {
      int retval = gt_kmercodeiterator_filetab_next(&kmercodeptr,
                                                    kmercodeiterator,
                                                    err);
      if (retval < 0)
      {
        haserr = true;
      } else
      {
        if (kmercodeptr != NULL)
        {
          outkmeroccurrence(codeliststream,kmercodeptr);
        } else
        {
          break;
        }
      }
    }
  }
  gt_kmercodeiterator_delete(kmercodeiterator);
  return haserr ? -1 : 0;
}

static int verifycodelists(const GtEncseq *encseq,
                           unsigned int kmersize,
                           unsigned int numofchars,
                           const GtArrayGtCodetype *codeliststream,
                           GtError *err)
{
  bool haserr = false;
  GtArrayGtCodetype codeliststring;
  const GtUchar *characters;
  unsigned long stringtotallength;

  gt_error_check(err);
  stringtotallength = gt_encseq_total_length(encseq);
  characters = gt_alphabet_characters(gt_encseq_alphabet(encseq));
  GT_INITARRAY(&codeliststring,GtCodetype);
  collectkmercode(&codeliststring,
                  encseq,
                  kmersize,
                  numofchars,
                  stringtotallength);
  if (comparecodelists(codeliststream,
                       &codeliststring,
                       kmersize,
                       numofchars,
                       (const char *) characters,
                       err) != 0)
  {
    haserr = true;
  }
  GT_FREEARRAY(&codeliststring,GtCodetype);
  return haserr ? -1 : 0;
}

int gt_verifymappedstr(const GtEncseq *encseq,
                       unsigned int prefixlength,
                       GtError *err)
{
  unsigned int numofchars;
  GtArrayGtCodetype codeliststream;
  bool haserr = false;

  gt_error_check(err);
  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  GT_INITARRAY(&codeliststream,GtCodetype);
  if (getfastastreamkmers(gt_encseq_filenames(encseq),
                          numofchars,
                          prefixlength,
                          gt_alphabet_symbolmap(
                                gt_encseq_alphabet(encseq)),
                          false,
                          &codeliststream,
                          err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (verifycodelists(encseq,
                        prefixlength,
                        numofchars,
                        &codeliststream,
                        err) != 0)
    {
      haserr = true;
    }
  }
  GT_FREEARRAY(&codeliststream,GtCodetype);
  return haserr ? -1 : 0;
}
