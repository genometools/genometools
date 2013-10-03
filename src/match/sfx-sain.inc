static GtUword gt_sain_PLAINSEQ_insertSstarsuffixes(GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                                 GtUsainindextype *suftab,
                                                 GtLogger *logger)
{
  GtUword nextcc = GT_UNIQUEINT(sainseq->totallength),
          countSstartype = 0;
  GtUsainindextype position, *fillptr = sainseq->bucketfillptr;
  GtSainbuffer *sainbuffer = gt_sainbuffer_new(suftab,fillptr,
                                               sainseq->numofchars,logger);
  bool nextisStype = true;

  gt_sain_endbuckets(sainseq);
  for (position = (GtUsainindextype) (sainseq->totallength-1); /* Nothing */;
       position--)
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
      printf("Sstar.suftab[" GT_WU "]=" GT_WU "\n",fillptr[nextcc],position+1);
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
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype position;

  gt_assert(sainseq->roundtable != NULL);
  for (suftabptr = suftab, sainseq->currentround = 0;
       suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc;

      if (position >= (GtSsainindextype) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtSsainindextype) sainseq->totallength;
      }
      currentcc = (GtUword)
plainseq[position];
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword t, leftcontextcc;

          position--;
          leftcontextcc = (GtUword)
plainseq[position];
          t = (currentcc << 1) | (leftcontextcc < currentcc ? 1UL : 0);
          gt_assert(currentcc > 0 &&
                    sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += (GtSsainindextype) sainseq->totallength;
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
          printf("L-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtSsainindextype position;

  gt_assert(sainseq->roundtable == NULL);
  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc = (GtUword)
plainseq[position];
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword leftcontextcc;

          position--;
          leftcontextcc = (GtUword)
plainseq[position];
          GT_SAINUPDATEBUCKETPTR(currentcc);
          /* negative => position does not derive L-suffix
             positive => position may derive L-suffix */
          gt_assert(suftabptr < bucketptr);
          *bucketptr++ = (leftcontextcc < currentcc) ? ~position : position;
          *suftabptr = 0;
#ifdef SAINSHOWSTATE
          gt_assert(bucketptr != NULL);
          printf("L-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtSsainindextype position;

  gt_assert(sainseq->roundtable != NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtSsainindextype)
                                    (sainseq->totallength-1));
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
      if (position >= (GtSsainindextype) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtSsainindextype) sainseq->totallength;
      }
      if (position > 0)
      {
        GtUword currentcc = (GtUword)
plainseq[position];

        if (currentcc < sainseq->numofchars)
        {
          GtUword t, leftcontextcc;

          position--;
          leftcontextcc = (GtUword)
plainseq[position];
          t = (currentcc << 1) | (leftcontextcc > currentcc ? 1UL : 0);
          gt_assert(sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += (GtSsainindextype) sainseq->totallength;
            sainseq->roundtable[t] = sainseq->currentround;
          }
          GT_SAINUPDATEBUCKETPTR(currentcc);
          gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
          *(--bucketptr) = (t & 1UL) ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
          printf("S-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtSsainindextype position;

  gt_assert(sainseq->roundtable == NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtSsainindextype)
                                    (sainseq->totallength-1));
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
plainseq[position];

      if (currentcc < sainseq->numofchars)
      {
        GtUword leftcontextcc;

        position--;
        leftcontextcc = (GtUword)
plainseq[position];
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (leftcontextcc > currentcc)
                          ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
        printf("S-induce: suftab[" GT_WU "]=" GT_WD "\n",
               (GtUword) (bucketptr - suftab),*bucketptr);
#endif
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_PLAINSEQ_induceLtypesuffixes2(const GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtSsainindextype position;

  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    position = *suftabptr;
    *suftabptr = ~position;
    if (position > 0)
    {
      GtUword currentcc;

      position--;
      currentcc = (GtUword)
plainseq[position];
      if (currentcc < sainseq->numofchars)
      {
        gt_assert(currentcc > 0);
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && suftabptr < bucketptr);
        *bucketptr++ = (position > 0 &&
                        ((GtUword)
plainseq[position-1])
                                            < currentcc)
                        ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("L-induce: suftab[" GT_WU "]=" GT_WD "\n",
               (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
      }
    }
  }
}

static void gt_sain_PLAINSEQ_induceStypesuffixes2(const GtSainseq *sainseq,
                                                 const GtUchar *plainseq,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtSsainindextype position;

  gt_sain_special_singleSinduction2(sainseq,
                                    suftab,
                                    (GtSsainindextype) sainseq->totallength,
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
plainseq[position];
      if (currentcc < sainseq->numofchars)
      {
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (position == 0 ||
                          ((GtUword)
plainseq[position-1]) > currentcc)
                         ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("S-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                         GtUsainindextype *suftab)
{
  GtUsainindextype *suftabptr,
                   position,
                   writeidx = (GtUsainindextype) (numberofsuffixes - 1),
                   *sstarsuffixes = suftab + numberofsuffixes,
                   *sstarfirstcharcount = NULL,
                   *bucketsize = NULL;
  GtUword nextcc = GT_UNIQUEINT(sainseq->totallength);
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
  for (position = (GtUsainindextype) (sainseq->totallength-1); /* Nothing */;
       position--)
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
                                                 GtUsainindextype *suftab,
                                                 GtLogger *logger)
{
  GtUword nextcc = GT_UNIQUEINT(sainseq->totallength),
          countSstartype = 0;
  GtUsainindextype position, *fillptr = sainseq->bucketfillptr;
  GtSainbuffer *sainbuffer = gt_sainbuffer_new(suftab,fillptr,
                                               sainseq->numofchars,logger);
  GtUchar tmpcc;
bool nextisStype = true;

  gt_sain_endbuckets(sainseq);
  for (position = (GtUsainindextype) (sainseq->totallength-1); /* Nothing */;
       position--)
  {
    GtUword currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position),
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
      printf("Sstar.suftab[" GT_WU "]=" GT_WU "\n",fillptr[nextcc],position+1);
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
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtUchar tmpcc;
GtSsainindextype position;

  gt_assert(sainseq->roundtable != NULL);
  for (suftabptr = suftab, sainseq->currentround = 0;
       suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc;

      if (position >= (GtSsainindextype) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtSsainindextype) sainseq->totallength;
      }
      currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position),
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword t, leftcontextcc;

          position--;
          leftcontextcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position),
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;
          t = (currentcc << 1) | (leftcontextcc < currentcc ? 1UL : 0);
          gt_assert(currentcc > 0 &&
                    sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += (GtSsainindextype) sainseq->totallength;
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
          printf("L-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtSsainindextype position;

  gt_assert(sainseq->roundtable == NULL);
  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position),
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword leftcontextcc;

          position--;
          leftcontextcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position),
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;
          GT_SAINUPDATEBUCKETPTR(currentcc);
          /* negative => position does not derive L-suffix
             positive => position may derive L-suffix */
          gt_assert(suftabptr < bucketptr);
          *bucketptr++ = (leftcontextcc < currentcc) ? ~position : position;
          *suftabptr = 0;
#ifdef SAINSHOWSTATE
          gt_assert(bucketptr != NULL);
          printf("L-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtSsainindextype position;

  gt_assert(sainseq->roundtable != NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtSsainindextype)
                                    (sainseq->totallength-1));
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
      if (position >= (GtSsainindextype) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtSsainindextype) sainseq->totallength;
      }
      if (position > 0)
      {
        GtUword currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position),
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;

        if (currentcc < sainseq->numofchars)
        {
          GtUword t, leftcontextcc;

          position--;
          leftcontextcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position),
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;
          t = (currentcc << 1) | (leftcontextcc > currentcc ? 1UL : 0);
          gt_assert(sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += (GtSsainindextype) sainseq->totallength;
            sainseq->roundtable[t] = sainseq->currentround;
          }
          GT_SAINUPDATEBUCKETPTR(currentcc);
          gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
          *(--bucketptr) = (t & 1UL) ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
          printf("S-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtSsainindextype position;

  gt_assert(sainseq->roundtable == NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtSsainindextype)
                                    (sainseq->totallength-1));
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
(GtUword) (position),
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;

      if (currentcc < sainseq->numofchars)
      {
        GtUword leftcontextcc;

        position--;
        leftcontextcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position),
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (leftcontextcc > currentcc)
                          ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
        printf("S-induce: suftab[" GT_WU "]=" GT_WD "\n",
               (GtUword) (bucketptr - suftab),*bucketptr);
#endif
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_ENCSEQ_induceLtypesuffixes2(const GtSainseq *sainseq,
                                                 const GtEncseq *encseq,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtSsainindextype position;

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
(GtUword) (position),
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;
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
  ? GT_UNIQUEINT(position-1) : (GtUword) tmpcc)
                                            < currentcc)
                        ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("L-induce: suftab[" GT_WU "]=" GT_WD "\n",
               (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
      }
    }
  }
}

static void gt_sain_ENCSEQ_induceStypesuffixes2(const GtSainseq *sainseq,
                                                 const GtEncseq *encseq,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtUchar tmpcc;
GtSsainindextype position;

  gt_sain_special_singleSinduction2(sainseq,
                                    suftab,
                                    (GtSsainindextype) sainseq->totallength,
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
(GtUword) (position),
sainseq->readmode))
  ? GT_UNIQUEINT(position) : (GtUword) tmpcc;
      if (currentcc < sainseq->numofchars)
      {
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (position == 0 ||
                          (ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position-1),
sainseq->readmode))
  ? GT_UNIQUEINT(position-1) : (GtUword) tmpcc) > currentcc)
                         ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("S-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                         GtUsainindextype *suftab)
{
  GtUsainindextype *suftabptr,
                   position,
                   writeidx = (GtUsainindextype) (numberofsuffixes - 1),
                   *sstarsuffixes = suftab + numberofsuffixes,
                   *sstarfirstcharcount = NULL,
                   *bucketsize = NULL;
  GtUword nextcc = GT_UNIQUEINT(sainseq->totallength);
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
  for (position = (GtUsainindextype) (sainseq->totallength-1); /* Nothing */;
       position--)
  {
    GtUword currentcc = ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(
encseq,
(GtUword) (position),
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
                                                 const GtUsainindextype *array,
                                                 GtUsainindextype *suftab,
                                                 GtLogger *logger)
{
  GtUword nextcc = GT_UNIQUEINT(sainseq->totallength),
          countSstartype = 0;
  GtUsainindextype position, *fillptr = sainseq->bucketfillptr;
  GtSainbuffer *sainbuffer = gt_sainbuffer_new(suftab,fillptr,
                                               sainseq->numofchars,logger);
  bool nextisStype = true;

  gt_sain_endbuckets(sainseq);
  for (position = (GtUsainindextype) (sainseq->totallength-1); /* Nothing */;
       position--)
  {
    GtUword currentcc = (GtUword) array[position];
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
      printf("Sstar.suftab[" GT_WU "]=" GT_WU "\n",fillptr[nextcc],position+1);
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
                                                 const GtUsainindextype *array,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype position;

  gt_assert(sainseq->roundtable != NULL);
  for (suftabptr = suftab, sainseq->currentround = 0;
       suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc;

      if (position >= (GtSsainindextype) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtSsainindextype) sainseq->totallength;
      }
      currentcc = (GtUword) array[position];
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword t, leftcontextcc;

          position--;
          leftcontextcc = (GtUword) array[position];
          t = (currentcc << 1) | (leftcontextcc < currentcc ? 1UL : 0);
          gt_assert(currentcc > 0 &&
                    sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += (GtSsainindextype) sainseq->totallength;
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
          printf("L-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                                 const GtUsainindextype *array,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtSsainindextype position;

  gt_assert(sainseq->roundtable == NULL);
  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc = (GtUword) array[position];
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword leftcontextcc;

          position--;
          leftcontextcc = (GtUword) array[position];
          GT_SAINUPDATEBUCKETPTR(currentcc);
          /* negative => position does not derive L-suffix
             positive => position may derive L-suffix */
          gt_assert(suftabptr < bucketptr);
          *bucketptr++ = (leftcontextcc < currentcc) ? ~position : position;
          *suftabptr = 0;
#ifdef SAINSHOWSTATE
          gt_assert(bucketptr != NULL);
          printf("L-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                                 const GtUsainindextype *array,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtSsainindextype position;

  gt_assert(sainseq->roundtable != NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtSsainindextype)
                                    (sainseq->totallength-1));
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
      if (position >= (GtSsainindextype) sainseq->totallength)
      {
        sainseq->currentround++;
        position -= (GtSsainindextype) sainseq->totallength;
      }
      if (position > 0)
      {
        GtUword currentcc = (GtUword) array[position];

        if (currentcc < sainseq->numofchars)
        {
          GtUword t, leftcontextcc;

          position--;
          leftcontextcc = (GtUword) array[position];
          t = (currentcc << 1) | (leftcontextcc > currentcc ? 1UL : 0);
          gt_assert(sainseq->roundtable[t] <= sainseq->currentround);
          if (sainseq->roundtable[t] < sainseq->currentround)
          {
            position += (GtSsainindextype) sainseq->totallength;
            sainseq->roundtable[t] = sainseq->currentround;
          }
          GT_SAINUPDATEBUCKETPTR(currentcc);
          gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
          *(--bucketptr) = (t & 1UL) ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
          printf("S-induce: suftab[" GT_WU "]=" GT_WD "\n",
                  (GtUword) (bucketptr - suftab),*bucketptr);
#endif
        }
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_INTSEQ_induceStypesuffixes1(GtSainseq *sainseq,
                                                 const GtUsainindextype *array,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtSsainindextype position;

  gt_assert(sainseq->roundtable == NULL);
  gt_sain_special_singleSinduction1(sainseq,
                                    suftab,
                                    (GtSsainindextype)
                                    (sainseq->totallength-1));
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
      GtUword currentcc = (GtUword) array[position];

      if (currentcc < sainseq->numofchars)
      {
        GtUword leftcontextcc;

        position--;
        leftcontextcc = (GtUword) array[position];
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (leftcontextcc > currentcc)
                          ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
        printf("S-induce: suftab[" GT_WU "]=" GT_WD "\n",
               (GtUword) (bucketptr - suftab),*bucketptr);
#endif
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_INTSEQ_induceLtypesuffixes2(const GtSainseq *sainseq,
                                                 const GtUsainindextype *array,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtSsainindextype position;

  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    position = *suftabptr;
    *suftabptr = ~position;
    if (position > 0)
    {
      GtUword currentcc;

      position--;
      currentcc = (GtUword) array[position];
      if (currentcc < sainseq->numofchars)
      {
        gt_assert(currentcc > 0);
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && suftabptr < bucketptr);
        *bucketptr++ = (position > 0 &&
                        ((GtUword) array[position-1])
                                            < currentcc)
                        ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("L-induce: suftab[" GT_WU "]=" GT_WD "\n",
               (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
      }
    }
  }
}

static void gt_sain_INTSEQ_induceStypesuffixes2(const GtSainseq *sainseq,
                                                 const GtUsainindextype *array,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0;
  GtUsainindextype *fillptr = sainseq->bucketfillptr;
  GtSsainindextype *suftabptr, *bucketptr = NULL;
  GtSsainindextype position;

  gt_sain_special_singleSinduction2(sainseq,
                                    suftab,
                                    (GtSsainindextype) sainseq->totallength,
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
      currentcc = (GtUword) array[position];
      if (currentcc < sainseq->numofchars)
      {
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (position == 0 ||
                          ((GtUword) array[position-1]) > currentcc)
                         ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("S-induce: suftab[" GT_WU "]=" GT_WD "\n",
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
                                                 const GtUsainindextype *array,
                                         GtUword numberofsuffixes,
                                         GtUsainindextype *suftab)
{
  GtUsainindextype *suftabptr,
                   position,
                   writeidx = (GtUsainindextype) (numberofsuffixes - 1),
                   *sstarsuffixes = suftab + numberofsuffixes,
                   *sstarfirstcharcount = NULL,
                   *bucketsize = NULL;
  GtUword nextcc = GT_UNIQUEINT(sainseq->totallength);
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
  for (position = (GtUsainindextype) (sainseq->totallength-1); /* Nothing */;
       position--)
  {
    GtUword currentcc = (GtUword) array[position];
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
