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

#include <string.h>
#include <inttypes.h>
#include "core/fa.h"
#include "core/alphabet.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/encseq_metadata.h"
#include "core/error.h"
#include "core/str.h"
#include "core/xansi_api.h"
#include "core/ma_api.h"

#include "emimergeesa.h"
#include "esa-fileend.h"
#include "fmindex.h"
#include "esa-map.h"

#include "encseq2offset.pr"
#include "fmi-keyval.pr"
#include "fmi-mapspec.pr"

static int copytheindexfile(const char *destindex,
                            const char *sourceindex,
                            const char *suffix,
                            uint64_t maxlength,
                            GtError *err)
{
  FILE *fpdest = NULL, *fpsource = NULL;
  int cc;
  bool haserr = false;

  gt_error_check(err);
  fpdest = gt_fa_fopen_with_suffix(destindex,suffix,"wb",err);
  if (fpdest == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    fpsource = gt_fa_fopen_with_suffix(sourceindex,suffix,"rb",err);
    if (fpsource == NULL)
    {
      haserr = true;
    }
  }
  printf("# cp %s%s %s%s\n",sourceindex,suffix,destindex,suffix);
  if (!haserr)
  {
    if (maxlength == 0)
    {
      while ((cc = fgetc(fpsource)) != EOF)
      {
        (void) putc(cc,fpdest);
      }
    } else
    {
      uint64_t pos;

      for (pos = 0; pos < maxlength; pos++)
      {
        if ((cc = fgetc(fpsource)) == EOF)
        {
          break;
        }
        (void) putc(cc,fpdest);
      }
    }
  }
  gt_fa_xfclose(fpdest);
  gt_fa_xfclose(fpsource);
  return haserr ? -1 : 0;
}

static void allocatefmtables(Fmindex *fm,
                             const GtSpecialcharinfo *specialcharinfo,
                             bool storeindexpos)
{
  fm->tfreq = gt_malloc(sizeof *fm->tfreq * TFREQSIZE(fm->mapsize));
  fm->superbfreq = gt_malloc(sizeof *fm->superbfreq
                             * SUPERBFREQSIZE(fm->mapsize,fm->nofsuperblocks));
  if (storeindexpos)
  {
    fm->markpostable = gt_malloc(sizeof *fm->markpostable
                                 * MARKPOSTABLELENGTH(fm->bwtlength,
                                                      fm->markdist));
    fm->specpos.nextfreeGtPairBwtidx = 0;
    fm->specpos.allocatedGtPairBwtidx
      = (unsigned long) gt_determinenumberofspecialstostore(specialcharinfo);
    printf("# %lu wildcards in the last %lu characters (%.2f)\n",
           (unsigned long) specialcharinfo->specialcharacters -
                           fm->specpos.allocatedGtPairBwtidx,
           specialcharinfo->specialcharacters,
            (double) (specialcharinfo->specialcharacters -
                      fm->specpos.allocatedGtPairBwtidx)/
                     specialcharinfo->specialcharacters);
    fm->specpos.spaceGtPairBwtidx
      = gt_malloc(sizeof *fm->specpos.spaceGtPairBwtidx
                  * fm->specpos.allocatedGtPairBwtidx);
  } else
  {
    GT_INITARRAY(&fm->specpos,GtPairBwtidx);
    fm->markpostable = NULL;
  }
  fm->bfreq = gt_malloc(sizeof *fm->bfreq
                        * BFREQSIZE(fm->mapsize,fm->nofblocks));
}

static void set0frequencies(Fmindex *fm)
{
  unsigned long i;

  for (i = 0; i < (unsigned long) TFREQSIZE(fm->mapsize); i++)
  {
    fm->tfreq[i] = 0;
  }
  for (i = 0; i < (unsigned long) BFREQSIZE(fm->mapsize,fm->nofblocks); i++)
  {
    fm->bfreq[i] = 0;
  }
  for (i = 0;
       i < (unsigned long) SUPERBFREQSIZE(fm->mapsize,fm->nofsuperblocks);
       i++)
  {
    fm->superbfreq[i] = 0;
  }
}

static void finalizefmfrequencies(Fmindex *fm)
{
  unsigned int j;
  unsigned long i, *freqptr;

  for (j = 2U; j <= fm->mapsize; j++)
  {
    fm->tfreq[j] += fm->tfreq[j - 1];
  }
  freqptr = fm->superbfreq;
  for (j = 0; j < fm->mapsize; j++)
  {
    for (i = (unsigned long) 2; i < fm->nofsuperblocks; i++)
    {
      freqptr[i] += freqptr[i-1];
    }
    freqptr += fm->nofsuperblocks;
  }
}

static void showconstructionmessage(const char *indexname,
                                    unsigned long totallength,
                                    unsigned long fmsize,
                                    unsigned int log2bsize,
                                    unsigned int log2markdist,
                                    unsigned int numofchars)
{
  printf("# construct fmindex \"%s\" for bsize=%lu, superbsize=%lu,",
          indexname,
          (unsigned long) GT_POW2(log2bsize),
          (unsigned long) GT_POW2(log2markdist));
  printf(" len=%lu, alphasize=%u: size ",
          totallength,
          numofchars);
  printf("%lu bytes, space overhead %.2f\n",
          fmsize,
          (double) fmsize/(double) (totallength+1));
}

static int nextesamergedsufbwttabvalues(Definedunsignedlong *longest,
                                       GtUchar *bwtvalue,
                                       unsigned long *suftabvalue,
                                       Emissionmergedesa *emmesa,
                                       const unsigned long *sequenceoffsettable,
                                       unsigned long bwtpos,
                                       GtError *err)
{
  Indexedsuffix indexedsuffix;

  gt_error_check(err);
  if (emmesa->buf.nextaccessidx >= emmesa->buf.nextstoreidx)
  {
    if (emmesa->numofentries == 0)
    {
      return 0;
    }
    if (gt_emissionmergedesa_stepdeleteandinsertothersuffixes(emmesa,err) != 0)
    {
      return -1;
    }
    if (emmesa->buf.nextstoreidx == 0)
    {
      return 0;
    }
    emmesa->buf.nextaccessidx = 0;
  }
  indexedsuffix = emmesa->buf.suftabstore[emmesa->buf.nextaccessidx];
  *suftabvalue = sequenceoffsettable[indexedsuffix.idx] +
                 indexedsuffix.startpos;
  if (indexedsuffix.startpos == 0)
  {
    if (indexedsuffix.idx == 0)
    {
      if (longest->defined)
      {
        gt_error_set(err,"longest is already defined as %lu",
                      longest->valueunsignedlong);
        return -2;
      }
      longest->defined = true;
      longest->valueunsignedlong = bwtpos;
      *bwtvalue = (GtUchar) UNDEFBWTCHAR;
    } else
    {
      *bwtvalue = (GtUchar) SEPARATOR;
    }
  } else
  {
    *bwtvalue
      = gt_encseq_get_encoded_char( /* Random access */
           emmesa->suffixarraytable[indexedsuffix.idx].encseq,
           indexedsuffix.startpos-1,
           emmesa->suffixarraytable[indexedsuffix.idx].readmode);
  }
  emmesa->buf.nextaccessidx++;
  return 1;
}

int gt_sufbwt2fmindex(Fmindex *fmindex,
                   GtSpecialcharinfo *specialcharinfo,
                   unsigned int log2bsize,
                   unsigned int log2markdist,
                   const char *outfmindex,
                   const GtStrArray *indexnametab,
                   bool storeindexpos,
                   GtLogger *logger,
                   GtError *err)
{
  Suffixarray suffixarray;
  Emissionmergedesa emmesa;
  GtUchar cc;
  unsigned long bwtpos,
         totallength = 0,
         suftabvalue = 0,
         *sequenceoffsettable = NULL,
         firstignorespecial = 0,
         nextmark,
         *markptr,
         nextprogress,
         tmpsuftabvalue,
         stepprogress;
  unsigned int numofchars = 0,
               suffixlength = 0,
               numofindexes;
  int retval;
  Definedunsignedlong longest = { false, 0 };
  GtPairBwtidx *pairptr;
  FILE *outbwt = NULL;
  GtStr *tmpfilename = NULL;
  bool haserr = false;

  gt_error_check(err);
  longest.defined = false;
  longest.valueunsignedlong = 0;
  numofindexes = (unsigned int) gt_str_array_size(indexnametab);
  if (numofindexes == 1U)
  {
    const char *indexname = gt_str_array_get(indexnametab,0);

    if (streamsuffixarray(&suffixarray,
                          SARR_BWTTAB | (storeindexpos ? SARR_SUFTAB : 0),
                          indexname,
                          logger,
                          err) != 0)
    {
      haserr = true;
    } else
    {
      totallength = gt_encseq_total_length(suffixarray.encseq);
    }
    if (!haserr && gt_specialcharinfo_read(specialcharinfo,indexname,err) != 0)
    {
      haserr = true;
    }
    if (!haserr)
    {
      numofchars = gt_alphabet_num_of_chars(
                               gt_encseq_alphabet(suffixarray.encseq));
      firstignorespecial = totallength - specialcharinfo->specialcharacters;
      if (gt_alphabet_to_file(gt_encseq_alphabet(suffixarray.encseq),
                              outfmindex, err) != 0)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (copytheindexfile(outfmindex,
                            indexname,
                            BWTTABSUFFIX,
                            (uint64_t) firstignorespecial,
                            err) != 0)
      {
        haserr = true;
      }
    }
  } else
  {
    GtEncseqMetadata *emd = NULL;
    if (gt_emissionmergedesa_init(&emmesa,
                               indexnametab,
                               SARR_ESQTAB | SARR_SUFTAB | SARR_LCPTAB,
                               logger,
                               err) != 0)
    {
      haserr = true;
    }
    if (!haserr)
    {
      const char *indexname = NULL;
      indexname = gt_str_array_get(indexnametab,0);
      emd = gt_encseq_metadata_new(indexname, err);
      if (emd == NULL) {
        haserr = true;
      }
    }
    if (!haserr) {
      suffixlength = 0;
      if (gt_alphabet_to_file(gt_encseq_metadata_alphabet(emd),
                              outfmindex, err) != 0)
      {
        haserr = true;
      }
    }
    gt_encseq_metadata_delete(emd);
    if (!haserr)
    {
      sequenceoffsettable = gt_encseqtable2sequenceoffsets(&totallength,
                                                        specialcharinfo,
                                                        emmesa.suffixarraytable,
                                                        numofindexes);
      if (sequenceoffsettable == NULL)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      longest.defined = false;
      longest.valueunsignedlong = 0;
      outbwt = gt_fa_fopen_with_suffix(outfmindex,BWTTABSUFFIX,"wb",err);
      if (outbwt == NULL)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      numofchars = emmesa.numofchars;
      firstignorespecial = totallength - specialcharinfo->specialcharacters;
    }
  }
  if (!haserr)
  {
    printf("# firstignorespecial=%lu\n",
              firstignorespecial);
    gt_computefmkeyvalues (fmindex,
                        specialcharinfo,
                        totallength+1,
                        log2bsize,
                        log2markdist,
                        numofchars,
                        suffixlength,
                        storeindexpos);
    showconstructionmessage(outfmindex,
                            totallength,
                            fmindex->sizeofindex,
                            log2bsize,
                            log2markdist,
                            numofchars);
    allocatefmtables(fmindex,specialcharinfo,storeindexpos);
    set0frequencies(fmindex);
    if (storeindexpos)
    {
      markptr = fmindex->markpostable;
    } else
    {
      markptr = NULL;
    }
    nextprogress = stepprogress = totallength/78;
    for (bwtpos = 0, nextmark = 0; ; bwtpos++)
    {
      if (numofindexes == 1U)
      {
        if (storeindexpos)
        {
          retval = gt_readnextfromstream_GtUlong(&tmpsuftabvalue,
                                            &suffixarray.suftabstream_GtUlong);
          if (retval == 0)
          {
            break;
          }
          suftabvalue = (unsigned long) tmpsuftabvalue;
        }
        retval = gt_readnextfromstream_GtUchar(&cc,&suffixarray.bwttabstream);
        if (retval == 0)
        {
          break;
        }
      } else
      {
        retval = nextesamergedsufbwttabvalues(&longest,
                                              &cc,
                                              &suftabvalue,
                                              &emmesa,
                                              sequenceoffsettable,
                                              bwtpos,
                                              err);
        if (retval < 0)
        {
          haserr = true;
          break;
        }
        if (retval == 0)
        {
          break;
        }
        gt_xfwrite(&cc, sizeof (GtUchar), (size_t) 1, outbwt);
      }
      if (bwtpos == nextprogress)
      {
        if (bwtpos == stepprogress)
        {
          (void) putchar('#');
        }
        (void) putchar('.');
        (void) fflush(stdout);
        nextprogress += stepprogress;
      }
      if (storeindexpos && bwtpos == nextmark)
      {
        *markptr++ = suftabvalue;
        nextmark += fmindex->markdist;
      }
      if (ISBWTSPECIAL(cc))
      {
        if (storeindexpos && bwtpos < firstignorespecial)
        {
          pairptr = fmindex->specpos.spaceGtPairBwtidx +
                    fmindex->specpos.nextfreeGtPairBwtidx++;
          if (pairptr >= fmindex->specpos.spaceGtPairBwtidx +
                         fmindex->specpos.allocatedGtPairBwtidx)
          {
            gt_error_set(err,"program error: not enough space for specpos");
            haserr = true;
            break;
          }
          pairptr->bwtpos = bwtpos;
          pairptr->suftabvalue = suftabvalue;
        }
      } else
      {
        fmindex->tfreq[cc+1]++;
        fmindex->bfreq[(cc * fmindex->nofblocks) +
                       (bwtpos >> fmindex->log2bsize)]++;
        fmindex->superbfreq[(cc * fmindex->nofsuperblocks) +
                       (bwtpos >> fmindex->log2superbsize) + 1]++;
      }
    }
  }
  if (!haserr)
  {
    if (storeindexpos &&
        fmindex->specpos.allocatedGtPairBwtidx !=
        fmindex->specpos.nextfreeGtPairBwtidx)
    {
      gt_error_set(err,"program error: too much space for specpos: "
                    "allocated = %lu != %lu = used",
                    fmindex->specpos.allocatedGtPairBwtidx,
                    fmindex->specpos.nextfreeGtPairBwtidx);
      haserr = true;
    }
  }
  if (!haserr)
  {
    (void) putchar('\n');
    finalizefmfrequencies(fmindex);
    if (fmindex->suffixlength > 0)
    {
      fmindex->boundarray = gt_malloc(sizeof *fmindex->boundarray
                                      * fmindex->numofcodes);
    }
    if (numofindexes == 1U)
    {
      fmindex->longestsuffixpos = suffixarray.longest.valueunsignedlong;
      gt_freesuffixarray(&suffixarray);
    } else
    {
      if (!longest.defined)
      {
        gt_error_set(err,"longest is not defined after merging");
        haserr = true;
      }
      if (!haserr)
      {
        fmindex->longestsuffixpos = longest.valueunsignedlong;
      }
      gt_fa_xfclose(outbwt);
      gt_emissionmergedesa_wrap(&emmesa);
    }
  }
  gt_free(sequenceoffsettable);
  gt_free(tmpfilename);
  return haserr ? -1 : 0;
}
