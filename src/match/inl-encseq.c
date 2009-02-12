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

#ifdef INLINEDENCSEQ
#include "core/chardef.h"
#include "core/filelengthvalues.h"
#include "core/fa.h"
#include "core/error.h"
#include "core/fastabuffer.h"
#include "core/unused_api.h"
#include "encseq-def.h"
#include "spacedef.h"
#include "opensfxfile.h"

#include "fillsci.pr"

#define TISTABFILESUFFIX ".tis"

int flushencseqfile(const GtStr *indexname,Encodedsequence *encseq,
                    GtError *err)
{
  FILE *fp;
  bool haserr = false;

  gt_error_check(err);
  fp = opensfxfile(indexname,TISTABFILESUFFIX,"wb",err);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (fwrite(encseq->plainseq,sizeof (Uchar),
               (size_t) encseq->totallength, fp)
               != (size_t) encseq->totallength)
    {
      haserr = true;
    }
  }
  gt_xfclose(fp);
  return haserr ? -1 : 0;
}

void freeEncodedsequence(Encodedsequence **encseqptr)
{
  Encodedsequence *encseq = *encseqptr;

  if (encseq == NULL)
  {
    return;
  }
  if (encseq->hasownmemory)
  {
    FREESPACE(encseq->plainseq);
  }
  if (encseq->mappedfile)
  {
    gt_xmunmap((void *) encseq->plainseq);
  }
  FREESPACE(*encseqptr);
}

Encodedsequencescanstate *newEncodedsequencescanstate(void)
{
  Encodedsequencescanstate *esr;

  ALLOCASSIGNSPACE(esr,NULL,Encodedsequencescanstate,1);
  return esr;
}

void initEncodedsequencescanstate(Encodedsequencescanstate *esr,
                                  GT_UNUSED const Encodedsequence *encseq,
                                  Readmode readmode,
                                  Seqpos startpos)
{
  esr->readmode = readmode;
}

void freeEncodedsequencescanstate(Encodedsequencescanstate **esr)
{
  FREESPACE(*esr);
}

static int fillplainseq(Encodedsequence *encseq,GtFastaBuffer *fbs,
                        GtError *err)
{
  Seqpos pos;
  int retval;
  Uchar cc;

  gt_error_check(err);
  ALLOCASSIGNSPACE(encseq->plainseq,NULL,Uchar,encseq->totallength);
  encseq->hasownmemory = true;
  encseq->mappedfile = false;
  encseq->hasspecialcharacters = false;
  for (pos=0; /* Nothing */; pos++)
  {
    retval = gt_fastabuffer_next(fbs,&cc,err);
    if (retval < 0)
    {
      FREESPACE(encseq->plainseq);
      return -1;
    }
    if (retval == 0)
    {
      break;
    }
    if (!encseq->hasspecialcharacters)
    {
      if (ISSPECIAL(cc))
      {
        encseq->hasspecialcharacters = true;
      }
    }
    encseq->plainseq[pos] = cc;
  }
  return 0;
}

/*@null@*/ Encodedsequence *files2encodedsequence(bool withrange,
                                                  const GtStrArray
                                                  *filenametab,
                                                  bool plainformat,
                                                  Seqpos totallength,
                                                  unsigned long numofsequences,
                                                  Seqpos specialranges,
                                                  const SfxAlphabet *alphabet,
                                                  const char *str_sat,
                                                  Verboseinfo *verboseinfo,
                                                  unsigned long
                                                    *characterdistribution,
                                                  GtError *err)
{
  Encodedsequence *encseq;
  GtFastaBuffer *fb = NULL;
  bool haserr = false;

  gt_error_check(err);
  fb = gt_fastabuffer_new(filenametab,
                          plainformat ? NULL : getsymbolmapAlphabet(alphabet),
                          plainformat,
                          NULL,
                          NULL,
                          NULL);
  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  encseq->totallength = totallength;
  if (fillplainseq(encseq,fb,err) != 0)
  {
    freeEncodedsequence(&encseq);
    haserr = true;
  }
  gt_fastabuffer_delete(fb);
  return haserr ? NULL : encseq;
}

/*@null@*/ Encodedsequence *mapencodedsequence(
                                   GT_UNUSED bool withrange,
                                   const GtStr *indexname,
                                   GT_UNUSED Seqpos totallength,
                                   Seqpos specialranges,
                                   GT_UNUSED unsigned int numofchars,
                                   GT_UNUSED Verboseinfo *verboseinfo,
                                   GtError *err)
{
  Encodedsequence *encseq;
  GtStr *tmpfilename;

  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,TISTABFILESUFFIX);
  encseq->plainseq = gt_mmap_read(gt_str_get(tmpfilename),
                                  (size_t *) &encseq->totallength);
  gt_str_delete(tmpfilename);
  encseq->hasownmemory = false;
  encseq->mappedfile = true;
  encseq->hasspecialcharacters = (specialranges > 0) ?  true : false;
  return encseq;
}

Encodedsequence *plain2encodedsequence(
                         GT_UNUSED bool withrange,
                         Specialcharinfo *specialcharinfo,
                         const Uchar *seq1,
                         Seqpos len1,
                         const Uchar *seq2,
                         unsigned long len2,
                         unsigned int numofchars,
                         Verboseinfo *verboseinfo)
{
  Encodedsequence *encseq;
  Uchar *seqptr;
  Seqpos pos, len;

  gt_assert(seq1 != NULL);
  gt_assert(len1 > 0);
  if (seq2 == NULL)
  {
    seqptr = (Uchar *) seq1;
    len = len1;
  } else
  {
    len = len1 + (Seqpos) len2 + 1;
    ALLOCASSIGNSPACE(seqptr,NULL,Uchar,len);
    memcpy(seqptr,seq1,sizeof (Uchar) * len1);
    seqptr[len1] = (Uchar) SEPARATOR;
    memcpy(seqptr + len1 + 1,seq2,sizeof (Uchar) * len2);
  }
  sequence2specialcharinfo(specialcharinfo,seqptr,len,verboseinfo);
  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,1);
  encseq->plainseq = seqptr;
  encseq->mappedfile = false;
  encseq->hasownmemory = (seq2 == NULL) ? false : true;
  encseq->totallength = len;
  encseq->hasspecialcharacters = false;
  for (pos=0; pos < encseq->totallength; pos++)
  {
    if (ISSPECIAL(encseq->plainseq[pos]))
    {
      encseq->hasspecialcharacters = true;
      break;
    }
  }
  return encseq;
}

Specialrangeiterator *newspecialrangeiterator(const Encodedsequence *encseq,
                                              bool moveforward)
{
  Specialrangeiterator *sri;

  ALLOCASSIGNSPACE(sri,NULL,Specialrangeiterator,1);
  sri->moveforward = moveforward;
  sri->exhausted = encseq->hasspecialcharacters ? false : true;
  sri->encseq = encseq;
  sri->lengthofspecialrange = 0;
  if (moveforward)
  {
    sri->pos = 0;
  } else
  {
    sri->pos = encseq->totallength-1;
  }
  gt_assert(sri != NULL);
  return sri;
}

bool hasspecialranges(const Encodedsequence *encseq)
{
  return encseq->hasspecialcharacters;
}

bool hasfastspecialrangeenumerator(GT_UNUSED const Encodedsequence *encseq)
{
  return false;
}

static bool bitanddirectnextspecialrangeiterator(Sequencerange *range,
                                                 Specialrangeiterator *sri)
{
  bool success = false;

  while (!success)
  {
    if (ISSPECIAL(sri->encseq->plainseq[sri->pos]))
    {
      sri->lengthofspecialrange++;
    } else
    {
      if (sri->lengthofspecialrange > 0)
      {
        if (sri->moveforward)
        {
          range->leftpos = sri->pos - sri->lengthofspecialrange;
          range->rightpos = sri->pos;
        } else
        {
          range->leftpos = sri->pos+1;
          range->rightpos = sri->pos+1+sri->lengthofspecialrange;
        }
        success = true;
        sri->lengthofspecialrange = 0;
      }
    }
    if (sri->moveforward)
    {
      if (sri->pos == sri->encseq->totallength - 1)
      {
        if (sri->lengthofspecialrange > 0)
        {
          range->leftpos = sri->encseq->totallength - sri->lengthofspecialrange;
          range->rightpos = sri->encseq->totallength;
          success = true;
        }
        sri->exhausted = true;
        break;
      }
      sri->pos++;
    } else
    {
      if (sri->pos == 0)
      {
        if (sri->lengthofspecialrange > 0)
        {
          range->leftpos = 0;
          range->rightpos = sri->lengthofspecialrange;
          success = true;
        }
        sri->exhausted = true;
        break;
      }
      sri->pos--;
    }
  }
  return success;
}

bool nextspecialrangeiterator(Sequencerange *range,Specialrangeiterator *sri)
{
  if (sri->exhausted)
  {
    return false;
  }
  return bitanddirectnextspecialrangeiterator(range,sri);
}

void freespecialrangeiterator(Specialrangeiterator **sri)
{
  FREESPACE(*sri);
}

const char *encseqaccessname(GT_UNUSED const Encodedsequence *encseq)
{
  return "direct";
}

void encseqextract(Uchar *buffer,
                   const Encodedsequence *encseq,
                   Seqpos frompos,
                   Seqpos topos)
{
  unsigned long idx;
  Seqpos pos;

  gt_assert(frompos < topos && topos < encseq->totallength);
  for (pos=frompos, idx = 0; pos <= topos; pos++, idx++)
  {
    buffer[idx] = getencodedchar(encseq,pos,Forwardmode);
  }
}
#endif /* ifdef INLINEDENCSEQ */
