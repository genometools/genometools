static GtUword gt_sain_PLAINSEQ_insertSstarsuffixes(GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                                 GtUword *suftab,
                                                 GtLogger *logger)
{
  GtUword position,
                nextcc = GT_UNIQUEINT(sainseq->totallength),
                countSstartype = 0,
                *fillptr = sainseq->bucketfillptr;
  GtSainbuffer *sainbuffer = gt_sainbuffer_new(suftab,fillptr,
                                               sainseq->numofchars,logger);
  bool nextisStype = true;

  gt_sain_endbuckets(sainseq);
  for (position = sainseq->totallength-1; /* Nothing */; position--)
  {
    GtUword currentcc = (GtUword)
plainseq[position];
    bool currentisStype = (currentcc < nextcc ||
                           (currentcc == nextcc && nextisStype)) ? true : false;
    if (!currentisStype && nextisStype)
    {
      countSstartype++;
      if (sainseq->sstarfirstcharcount != NULL)
      {
        sainseq->sstarfirstcharcount[nextcc]++;
      }
      if (sainbuffer != NULL)
      {
        gt_sainbuffer_update(sainbuffer,nextcc,position);
      } else
      {
        suftab[--fillptr[nextcc]] = position;
      }
#undef SAINSHOWSTATE
#ifdef SAINSHOWSTATE
      printf("Sstar.suftab[" GT_LU "]=" GT_LU "\n",fillptr[nextcc],position+1);
#endif
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  gt_sainbuffer_flushall(sainbuffer);
  gt_sainbuffer_delete(sainbuffer);
  gt_assert(GT_MULT2(countSstartype) <= sainseq->totallength);
  return countSstartype;
}

static void gt_sain_PLAINSEQ_fast_induceLtypesuffixes1(GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  gt_assert(sainseq->roundtable != NULL);
  for (suftabptr = suftab, sainseq->currentround = 0;
       suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc;

      if (position >= (GtWord) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtWord) sainseq->totallength;
      }
      currentcc = (GtUword)
plainseq[(GtUword) position];
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword t, leftcontextcc;

          gt_assert(position > 0);
          position--;
          leftcontextcc = (GtUword)
plainseq[(GtUword) position];
          t = (currentcc << 1) | (leftcontextcc < currentcc ? 1UL : 0);
          gt_assert(currentcc > 0 &&
                    sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += (GtWord) sainseq->totallength;
            sainseq->roundtable[t] = sainseq->currentround;
          }
          GT_SAINUPDATEBUCKETPTR(currentcc);
          /* negative => position does not derive L-suffix
             positive => position may derive L-suffix */
          gt_assert(suftabptr < bucketptr);
          *bucketptr++ = (t & 1UL) ? ~position : position;
          *suftabptr = 0;
#ifdef SAINSHOWSTATE
          gt_assert(bucketptr != NULL);
          printf("L-induce: suftab[" GT_LU "]=" GT_LD "\n",
                  (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
        }
      } else
      {
        *suftabptr = 0;
      }
    } else
    {
      if (position < 0)
      {
        *suftabptr = ~position;
      }
    }
  }
}

static void gt_sain_PLAINSEQ_induceLtypesuffixes1(GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  gt_assert(sainseq->roundtable == NULL);
  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc = (GtUword)
plainseq[(GtUword) position];
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword leftcontextcc;

          gt_assert(position > 0);
          position--;
          leftcontextcc = (GtUword)
plainseq[(GtUword) position];
          GT_SAINUPDATEBUCKETPTR(currentcc);
          /* negative => position does not derive L-suffix
             positive => position may derive L-suffix */
          gt_assert(suftabptr < bucketptr);
          *bucketptr++ = (leftcontextcc < currentcc) ? ~position : position;
          *suftabptr = 0;
#ifdef SAINSHOWSTATE
          gt_assert(bucketptr != NULL);
          printf("L-induce: suftab[" GT_LU "]=" GT_LD "\n",
                  (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
        }
      } else
      {
        *suftabptr = 0;
      }
    } else
    {
      if (position < 0)
      {
        *suftabptr = ~position;
      }
    }
  }
}

static void gt_sain_PLAINSEQ_fast_induceStypesuffixes1(GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  gt_assert(sainseq->roundtable != NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtWord) (sainseq->totallength-1));
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes1fromspecialranges(sainseq,
                                           sainseq->seq.encseq,
                                           suftab);
  }
  for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
       suftabptr--)
  {
    if ((position = *suftabptr) > 0)
    {
      if (position >= (GtWord) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtWord) sainseq->totallength;
      }
      if (position > 0)
      {
        GtUword currentcc = (GtUword)
plainseq[(GtUword) position];

        if (currentcc < sainseq->numofchars)
        {
          GtUword t, leftcontextcc;

          position--;
          leftcontextcc = (GtUword)
plainseq[(GtUword) position];
          t = (currentcc << 1) | (leftcontextcc > currentcc ? 1UL : 0);
          gt_assert(sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += sainseq->totallength;
            sainseq->roundtable[t] = sainseq->currentround;
          }
          GT_SAINUPDATEBUCKETPTR(currentcc);
          gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
          *(--bucketptr) = (t & 1UL) ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
          printf("S-induce: suftab[" GT_LU "]=" GT_LD "\n",
                  (GtUword) (bucketptr - suftab),*bucketptr);
#endif
        }
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_PLAINSEQ_induceStypesuffixes1(GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  gt_assert(sainseq->roundtable == NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtWord) (sainseq->totallength-1));
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes1fromspecialranges(sainseq,
                                           sainseq->seq.encseq,
                                           suftab);
  }
  for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
       suftabptr--)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc = (GtUword)
plainseq[(GtUword) position];

      if (currentcc < sainseq->numofchars)
      {
        GtUword leftcontextcc;

        position--;
        leftcontextcc = (GtUword)
plainseq[(GtUword) position];
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (leftcontextcc > currentcc)
                          ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
        printf("S-induce: suftab[" GT_LU "]=" GT_LD "\n",
               (GtUword) (bucketptr - suftab),*bucketptr);
#endif
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_PLAINSEQ_induceLtypesuffixes2(const GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    position = *suftabptr;
    *suftabptr = ~position;
    if (position > 0)
    {
      GtUword currentcc;

      position--;
      currentcc = (GtUword)
plainseq[(GtUword) position];
      if (currentcc < sainseq->numofchars)
      {
        gt_assert(currentcc > 0);
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && suftabptr < bucketptr);
        *bucketptr++ = (position > 0 &&
                        ((GtUword)
plainseq[(GtUword) (position-1)])
                                            < currentcc)
                        ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("L-induce: suftab[" GT_LU "]=" GT_LD "\n",
               (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
      }
    }
  }
}

static void gt_sain_PLAINSEQ_induceStypesuffixes2(const GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  gt_sain_special_singleSinduction2(sainseq,
                                    suftab,
                                    (GtWord) sainseq->totallength,
                                    nonspecialentries);
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes2fromspecialranges(sainseq,
                                           sainseq->seq.encseq,
                                           suftab,
                                           nonspecialentries);
  }
  if (nonspecialentries == 0)
  {
    return;
  }
  for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
       suftabptr--)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc;

      position--;
      currentcc = (GtUword)
plainseq[(GtUword) position];
      if (currentcc < sainseq->numofchars)
      {
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (position == 0 ||
                          ((GtUword)
plainseq[(GtUword) (position-1)])
                                             > currentcc)
                         ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("S-induce: suftab[" GT_LU "]=" GT_LD "\n",
                (GtUword) (bucketptr-suftab),*bucketptr);
#endif
      }
    } else
    {
      *suftabptr = ~position;
    }
  }
}

static void gt_sain_PLAINSEQ_expandorder2original(GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                         GtUword numberofsuffixes,
                                         GtUword *suftab)
{
  GtUword *suftabptr, position,
                writeidx = numberofsuffixes - 1,
                nextcc = GT_UNIQUEINT(sainseq->totallength),
                *sstarsuffixes = suftab + numberofsuffixes;
  GtUword *sstarfirstcharcount = NULL, *bucketsize = NULL;
  bool nextisStype = true;

  if (sainseq->seqtype == GT_SAIN_INTSEQ)
  {
    GtUword charidx;

    gt_assert(sainseq->sstarfirstcharcount == NULL);
    sstarfirstcharcount = sainseq->sstarfirstcharcount
                        = sainseq->bucketfillptr;
    bucketsize = sainseq->bucketsize;
    for (charidx = 0; charidx < sainseq->numofchars; charidx++)
    {
      sstarfirstcharcount[charidx] = 0;
      bucketsize[charidx] = 0;
    }
  }
  for (position = sainseq->totallength-1; /* Nothing */; position--)
  {
    GtUword currentcc = (GtUword)
plainseq[position];
    bool currentisStype = (currentcc < nextcc ||
                           (currentcc == nextcc && nextisStype)) ? true : false;

    if (!currentisStype && nextisStype)
    {
      if (sstarfirstcharcount != NULL)
      {
        sstarfirstcharcount[nextcc]++;
      }
      sstarsuffixes[writeidx--] = position+1;
    }
    if (bucketsize != NULL)
    {
      bucketsize[currentcc]++;
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  for (suftabptr = suftab; suftabptr < suftab + numberofsuffixes; suftabptr++)
  {
    *suftabptr = sstarsuffixes[*suftabptr];
  }
}
static GtUword gt_sain_ENCSEQ_insertSstarsuffixes(GtSainseq *sainseq,
                                                 const GtEncseq *encseq,
                                                 GtUword *suftab,
                                                 GtLogger *logger)
{
  GtUword position,
                nextcc = GT_UNIQUEINT(sainseq->totallength),
                countSstartype = 0,
                *fillptr = sainseq->bucketfillptr;
  GtSainbuffer *sainbuffer = gt_sainbuffer_new(suftab,fillptr,
                                               sainseq->numofchars,logger);
  GtUchar tmpcc;
bool nextisStype = true;

  gt_sain_endbuckets(sainseq);
  for (position = sainseq->totallength-1; /* Nothing */; position--)
  {
    GtUword currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
position,
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;
    bool currentisStype = (currentcc < nextcc ||
                           (currentcc == nextcc && nextisStype)) ? true : false;
    if (!currentisStype && nextisStype)
    {
      countSstartype++;
      if (sainseq->sstarfirstcharcount != NULL)
      {
        sainseq->sstarfirstcharcount[nextcc]++;
      }
      if (sainbuffer != NULL)
      {
        gt_sainbuffer_update(sainbuffer,nextcc,position);
      } else
      {
        suftab[--fillptr[nextcc]] = position;
      }
#undef SAINSHOWSTATE
#ifdef SAINSHOWSTATE
      printf("Sstar.suftab[" GT_LU "]=" GT_LU "\n",fillptr[nextcc],position+1);
#endif
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  gt_sainbuffer_flushall(sainbuffer);
  gt_sainbuffer_delete(sainbuffer);
  gt_assert(GT_MULT2(countSstartype) <= sainseq->totallength);
  return countSstartype;
}

static void gt_sain_ENCSEQ_fast_induceLtypesuffixes1(GtSainseq *sainseq,
                                                 const GtEncseq *encseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtWord position;

  gt_assert(sainseq->roundtable != NULL);
  for (suftabptr = suftab, sainseq->currentround = 0;
       suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc;

      if (position >= (GtWord) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtWord) sainseq->totallength;
      }
      currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) position,
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) position) : (GtUword) tmpcc;
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword t, leftcontextcc;

          gt_assert(position > 0);
          position--;
          leftcontextcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) position,
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) position) : (GtUword) tmpcc;
          t = (currentcc << 1) | (leftcontextcc < currentcc ? 1UL : 0);
          gt_assert(currentcc > 0 &&
                    sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += (GtWord) sainseq->totallength;
            sainseq->roundtable[t] = sainseq->currentround;
          }
          GT_SAINUPDATEBUCKETPTR(currentcc);
          /* negative => position does not derive L-suffix
             positive => position may derive L-suffix */
          gt_assert(suftabptr < bucketptr);
          *bucketptr++ = (t & 1UL) ? ~position : position;
          *suftabptr = 0;
#ifdef SAINSHOWSTATE
          gt_assert(bucketptr != NULL);
          printf("L-induce: suftab[" GT_LU "]=" GT_LD "\n",
                  (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
        }
      } else
      {
        *suftabptr = 0;
      }
    } else
    {
      if (position < 0)
      {
        *suftabptr = ~position;
      }
    }
  }
}

static void gt_sain_ENCSEQ_induceLtypesuffixes1(GtSainseq *sainseq,
                                                 const GtEncseq *encseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtWord position;

  gt_assert(sainseq->roundtable == NULL);
  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) position,
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) position) : (GtUword) tmpcc;
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword leftcontextcc;

          gt_assert(position > 0);
          position--;
          leftcontextcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) position,
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) position) : (GtUword) tmpcc;
          GT_SAINUPDATEBUCKETPTR(currentcc);
          /* negative => position does not derive L-suffix
             positive => position may derive L-suffix */
          gt_assert(suftabptr < bucketptr);
          *bucketptr++ = (leftcontextcc < currentcc) ? ~position : position;
          *suftabptr = 0;
#ifdef SAINSHOWSTATE
          gt_assert(bucketptr != NULL);
          printf("L-induce: suftab[" GT_LU "]=" GT_LD "\n",
                  (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
        }
      } else
      {
        *suftabptr = 0;
      }
    } else
    {
      if (position < 0)
      {
        *suftabptr = ~position;
      }
    }
  }
}

static void gt_sain_ENCSEQ_fast_induceStypesuffixes1(GtSainseq *sainseq,
                                                 const GtEncseq *encseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtWord position;

  gt_assert(sainseq->roundtable != NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtWord) (sainseq->totallength-1));
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes1fromspecialranges(sainseq,
                                           sainseq->seq.encseq,
                                           suftab);
  }
  for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
       suftabptr--)
  {
    if ((position = *suftabptr) > 0)
    {
      if (position >= (GtWord) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtWord) sainseq->totallength;
      }
      if (position > 0)
      {
        GtUword currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) position,
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) position) : (GtUword) tmpcc;

        if (currentcc < sainseq->numofchars)
        {
          GtUword t, leftcontextcc;

          position--;
          leftcontextcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) position,
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) position) : (GtUword) tmpcc;
          t = (currentcc << 1) | (leftcontextcc > currentcc ? 1UL : 0);
          gt_assert(sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += sainseq->totallength;
            sainseq->roundtable[t] = sainseq->currentround;
          }
          GT_SAINUPDATEBUCKETPTR(currentcc);
          gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
          *(--bucketptr) = (t & 1UL) ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
          printf("S-induce: suftab[" GT_LU "]=" GT_LD "\n",
                  (GtUword) (bucketptr - suftab),*bucketptr);
#endif
        }
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_ENCSEQ_induceStypesuffixes1(GtSainseq *sainseq,
                                                 const GtEncseq *encseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtWord position;

  gt_assert(sainseq->roundtable == NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtWord) (sainseq->totallength-1));
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes1fromspecialranges(sainseq,
                                           sainseq->seq.encseq,
                                           suftab);
  }
  for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
       suftabptr--)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) position,
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) position) : (GtUword) tmpcc;

      if (currentcc < sainseq->numofchars)
      {
        GtUword leftcontextcc;

        position--;
        leftcontextcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) position,
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) position) : (GtUword) tmpcc;
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (leftcontextcc > currentcc)
                          ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
        printf("S-induce: suftab[" GT_LU "]=" GT_LD "\n",
               (GtUword) (bucketptr - suftab),*bucketptr);
#endif
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_ENCSEQ_induceLtypesuffixes2(const GtSainseq *sainseq,
                                                 const GtEncseq *encseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtWord position;

  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    position = *suftabptr;
    *suftabptr = ~position;
    if (position > 0)
    {
      GtUword currentcc;

      position--;
      currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) position,
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) position) : (GtUword) tmpcc;
      if (currentcc < sainseq->numofchars)
      {
        gt_assert(currentcc > 0);
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && suftabptr < bucketptr);
        *bucketptr++ = (position > 0 &&
                        (ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position-1),
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) (position-1)) : (GtUword) tmpcc)
                                            < currentcc)
                        ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("L-induce: suftab[" GT_LU "]=" GT_LD "\n",
               (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
      }
    }
  }
}

static void gt_sain_ENCSEQ_induceStypesuffixes2(const GtSainseq *sainseq,
                                                 const GtEncseq *encseq,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtWord position;

  gt_sain_special_singleSinduction2(sainseq,
                                    suftab,
                                    (GtWord) sainseq->totallength,
                                    nonspecialentries);
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes2fromspecialranges(sainseq,
                                           sainseq->seq.encseq,
                                           suftab,
                                           nonspecialentries);
  }
  if (nonspecialentries == 0)
  {
    return;
  }
  for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
       suftabptr--)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc;

      position--;
      currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) position,
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) position) : (GtUword) tmpcc;
      if (currentcc < sainseq->numofchars)
      {
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (position == 0 ||
                          (ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position-1),
sainseq->readmode))
  ? GT_UNIQUEINT((GtUword) (position-1)) : (GtUword) tmpcc)
                                             > currentcc)
                         ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("S-induce: suftab[" GT_LU "]=" GT_LD "\n",
                (GtUword) (bucketptr-suftab),*bucketptr);
#endif
      }
    } else
    {
      *suftabptr = ~position;
    }
  }
}

static void gt_sain_ENCSEQ_expandorder2original(GtSainseq *sainseq,
                                                 const GtEncseq *encseq,
                                         GtUword numberofsuffixes,
                                         GtUword *suftab)
{
  GtUword *suftabptr, position,
                writeidx = numberofsuffixes - 1,
                nextcc = GT_UNIQUEINT(sainseq->totallength),
                *sstarsuffixes = suftab + numberofsuffixes;
  GtUword *sstarfirstcharcount = NULL, *bucketsize = NULL;
  GtUchar tmpcc;
bool nextisStype = true;

  if (sainseq->seqtype == GT_SAIN_INTSEQ)
  {
    GtUword charidx;

    gt_assert(sainseq->sstarfirstcharcount == NULL);
    sstarfirstcharcount = sainseq->sstarfirstcharcount
                        = sainseq->bucketfillptr;
    bucketsize = sainseq->bucketsize;
    for (charidx = 0; charidx < sainseq->numofchars; charidx++)
    {
      sstarfirstcharcount[charidx] = 0;
      bucketsize[charidx] = 0;
    }
  }
  for (position = sainseq->totallength-1; /* Nothing */; position--)
  {
    GtUword currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
position,
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;
    bool currentisStype = (currentcc < nextcc ||
                           (currentcc == nextcc && nextisStype)) ? true : false;

    if (!currentisStype && nextisStype)
    {
      if (sstarfirstcharcount != NULL)
      {
        sstarfirstcharcount[nextcc]++;
      }
      sstarsuffixes[writeidx--] = position+1;
    }
    if (bucketsize != NULL)
    {
      bucketsize[currentcc]++;
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  for (suftabptr = suftab; suftabptr < suftab + numberofsuffixes; suftabptr++)
  {
    *suftabptr = sstarsuffixes[*suftabptr];
  }
}
static GtUword gt_sain_INTSEQ_insertSstarsuffixes(GtSainseq *sainseq,
                                                 const GtUword *array,
                                                 GtUword *suftab,
                                                 GtLogger *logger)
{
  GtUword position,
                nextcc = GT_UNIQUEINT(sainseq->totallength),
                countSstartype = 0,
                *fillptr = sainseq->bucketfillptr;
  GtSainbuffer *sainbuffer = gt_sainbuffer_new(suftab,fillptr,
                                               sainseq->numofchars,logger);
  bool nextisStype = true;

  gt_sain_endbuckets(sainseq);
  for (position = sainseq->totallength-1; /* Nothing */; position--)
  {
    GtUword currentcc = array[position];
    bool currentisStype = (currentcc < nextcc ||
                           (currentcc == nextcc && nextisStype)) ? true : false;
    if (!currentisStype && nextisStype)
    {
      countSstartype++;
      if (sainseq->sstarfirstcharcount != NULL)
      {
        sainseq->sstarfirstcharcount[nextcc]++;
      }
      if (sainbuffer != NULL)
      {
        gt_sainbuffer_update(sainbuffer,nextcc,position);
      } else
      {
        suftab[--fillptr[nextcc]] = position;
      }
#undef SAINSHOWSTATE
#ifdef SAINSHOWSTATE
      printf("Sstar.suftab[" GT_LU "]=" GT_LU "\n",fillptr[nextcc],position+1);
#endif
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  gt_sainbuffer_flushall(sainbuffer);
  gt_sainbuffer_delete(sainbuffer);
  gt_assert(GT_MULT2(countSstartype) <= sainseq->totallength);
  return countSstartype;
}

static void gt_sain_INTSEQ_fast_induceLtypesuffixes1(GtSainseq *sainseq,
                                                 const GtUword *array,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  gt_assert(sainseq->roundtable != NULL);
  for (suftabptr = suftab, sainseq->currentround = 0;
       suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc;

      if (position >= (GtWord) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtWord) sainseq->totallength;
      }
      currentcc = array[(GtUword) position];
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword t, leftcontextcc;

          gt_assert(position > 0);
          position--;
          leftcontextcc = array[(GtUword) position];
          t = (currentcc << 1) | (leftcontextcc < currentcc ? 1UL : 0);
          gt_assert(currentcc > 0 &&
                    sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += (GtWord) sainseq->totallength;
            sainseq->roundtable[t] = sainseq->currentround;
          }
          GT_SAINUPDATEBUCKETPTR(currentcc);
          /* negative => position does not derive L-suffix
             positive => position may derive L-suffix */
          gt_assert(suftabptr < bucketptr);
          *bucketptr++ = (t & 1UL) ? ~position : position;
          *suftabptr = 0;
#ifdef SAINSHOWSTATE
          gt_assert(bucketptr != NULL);
          printf("L-induce: suftab[" GT_LU "]=" GT_LD "\n",
                  (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
        }
      } else
      {
        *suftabptr = 0;
      }
    } else
    {
      if (position < 0)
      {
        *suftabptr = ~position;
      }
    }
  }
}

static void gt_sain_INTSEQ_induceLtypesuffixes1(GtSainseq *sainseq,
                                                 const GtUword *array,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  gt_assert(sainseq->roundtable == NULL);
  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc = array[(GtUword) position];
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword leftcontextcc;

          gt_assert(position > 0);
          position--;
          leftcontextcc = array[(GtUword) position];
          GT_SAINUPDATEBUCKETPTR(currentcc);
          /* negative => position does not derive L-suffix
             positive => position may derive L-suffix */
          gt_assert(suftabptr < bucketptr);
          *bucketptr++ = (leftcontextcc < currentcc) ? ~position : position;
          *suftabptr = 0;
#ifdef SAINSHOWSTATE
          gt_assert(bucketptr != NULL);
          printf("L-induce: suftab[" GT_LU "]=" GT_LD "\n",
                  (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
        }
      } else
      {
        *suftabptr = 0;
      }
    } else
    {
      if (position < 0)
      {
        *suftabptr = ~position;
      }
    }
  }
}

static void gt_sain_INTSEQ_fast_induceStypesuffixes1(GtSainseq *sainseq,
                                                 const GtUword *array,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  gt_assert(sainseq->roundtable != NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtWord) (sainseq->totallength-1));
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes1fromspecialranges(sainseq,
                                           sainseq->seq.encseq,
                                           suftab);
  }
  for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
       suftabptr--)
  {
    if ((position = *suftabptr) > 0)
    {
      if (position >= (GtWord) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtWord) sainseq->totallength;
      }
      if (position > 0)
      {
        GtUword currentcc = array[(GtUword) position];

        if (currentcc < sainseq->numofchars)
        {
          GtUword t, leftcontextcc;

          position--;
          leftcontextcc = array[(GtUword) position];
          t = (currentcc << 1) | (leftcontextcc > currentcc ? 1UL : 0);
          gt_assert(sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += sainseq->totallength;
            sainseq->roundtable[t] = sainseq->currentround;
          }
          GT_SAINUPDATEBUCKETPTR(currentcc);
          gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
          *(--bucketptr) = (t & 1UL) ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
          printf("S-induce: suftab[" GT_LU "]=" GT_LD "\n",
                  (GtUword) (bucketptr - suftab),*bucketptr);
#endif
        }
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_INTSEQ_induceStypesuffixes1(GtSainseq *sainseq,
                                                 const GtUword *array,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  gt_assert(sainseq->roundtable == NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtWord) (sainseq->totallength-1));
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes1fromspecialranges(sainseq,
                                           sainseq->seq.encseq,
                                           suftab);
  }
  for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
       suftabptr--)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc = array[(GtUword) position];

      if (currentcc < sainseq->numofchars)
      {
        GtUword leftcontextcc;

        position--;
        leftcontextcc = array[(GtUword) position];
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (leftcontextcc > currentcc)
                          ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
        printf("S-induce: suftab[" GT_LU "]=" GT_LD "\n",
               (GtUword) (bucketptr - suftab),*bucketptr);
#endif
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_INTSEQ_induceLtypesuffixes2(const GtSainseq *sainseq,
                                                 const GtUword *array,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    position = *suftabptr;
    *suftabptr = ~position;
    if (position > 0)
    {
      GtUword currentcc;

      position--;
      currentcc = array[(GtUword) position];
      if (currentcc < sainseq->numofchars)
      {
        gt_assert(currentcc > 0);
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && suftabptr < bucketptr);
        *bucketptr++ = (position > 0 &&
                        (array[(GtUword) (position-1)])
                                            < currentcc)
                        ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("L-induce: suftab[" GT_LU "]=" GT_LD "\n",
               (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
      }
    }
  }
}

static void gt_sain_INTSEQ_induceStypesuffixes2(const GtSainseq *sainseq,
                                                 const GtUword *array,
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  GtWord position;

  gt_sain_special_singleSinduction2(sainseq,
                                    suftab,
                                    (GtWord) sainseq->totallength,
                                    nonspecialentries);
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes2fromspecialranges(sainseq,
                                           sainseq->seq.encseq,
                                           suftab,
                                           nonspecialentries);
  }
  if (nonspecialentries == 0)
  {
    return;
  }
  for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
       suftabptr--)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc;

      position--;
      currentcc = array[(GtUword) position];
      if (currentcc < sainseq->numofchars)
      {
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (position == 0 ||
                          (array[(GtUword) (position-1)])
                                             > currentcc)
                         ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("S-induce: suftab[" GT_LU "]=" GT_LD "\n",
                (GtUword) (bucketptr-suftab),*bucketptr);
#endif
      }
    } else
    {
      *suftabptr = ~position;
    }
  }
}

static void gt_sain_INTSEQ_expandorder2original(GtSainseq *sainseq,
                                                 const GtUword *array,
                                         GtUword numberofsuffixes,
                                         GtUword *suftab)
{
  GtUword *suftabptr, position,
                writeidx = numberofsuffixes - 1,
                nextcc = GT_UNIQUEINT(sainseq->totallength),
                *sstarsuffixes = suftab + numberofsuffixes;
  GtUword *sstarfirstcharcount = NULL, *bucketsize = NULL;
  bool nextisStype = true;

  if (sainseq->seqtype == GT_SAIN_INTSEQ)
  {
    GtUword charidx;

    gt_assert(sainseq->sstarfirstcharcount == NULL);
    sstarfirstcharcount = sainseq->sstarfirstcharcount
                        = sainseq->bucketfillptr;
    bucketsize = sainseq->bucketsize;
    for (charidx = 0; charidx < sainseq->numofchars; charidx++)
    {
      sstarfirstcharcount[charidx] = 0;
      bucketsize[charidx] = 0;
    }
  }
  for (position = sainseq->totallength-1; /* Nothing */; position--)
  {
    GtUword currentcc = array[position];
    bool currentisStype = (currentcc < nextcc ||
                           (currentcc == nextcc && nextisStype)) ? true : false;

    if (!currentisStype && nextisStype)
    {
      if (sstarfirstcharcount != NULL)
      {
        sstarfirstcharcount[nextcc]++;
      }
      sstarsuffixes[writeidx--] = position+1;
    }
    if (bucketsize != NULL)
    {
      bucketsize[currentcc]++;
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  for (suftabptr = suftab; suftabptr < suftab + numberofsuffixes; suftabptr++)
  {
    *suftabptr = sstarsuffixes[*suftabptr];
  }
}
