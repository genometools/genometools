#!/usr/bin/env ruby

keylist = ["PLAINSEQ","ENCSEQ","INTSEQ"]

def getc_param(key)
  if key == "PLAINSEQ"
    return "const GtUchar *plainseq"
  elsif key == "INTSEQ"
    return "const GtUword *array"
  else
    return "const GtEncseq *encseq"
  end
end

def getc_call(key,posexpr)
  if key == "PLAINSEQ"
    return "(GtUword)\nplainseq[#{posexpr}]"
  elsif key == "INTSEQ"
    return "array[#{posexpr}]"
  else
    return "ISSPECIAL(tmpcc = gt_encseq_get_encoded_char(\n" +
                                  "encseq,\n" +
                                  "#{posexpr},\n" +
                                  "sainseq->readmode))\n" +
               "  ? GT_UNIQUEINT(#{posexpr}) : (GtUword) tmpcc"
  end
end

def declare_tmpcc(key,context)
  if key == "ENCSEQ"
    return "GtUchar tmpcc;\n#{context};"
  else
    return "#{context};"
  end
end

begin
  fo = File.open("src/match/sfx-sain.inc","w")
rescue => err
  STDERR.puts "#{$0}: #{err}"
  exit 1
end

keylist.each do |key|
fo.puts <<CCODE
static GtUword gt_sain_#{key}_insertSstarsuffixes(GtSainseq *sainseq,
                                                 #{getc_param(key)},
                                                 GtUword *suftab,
                                                 GtLogger *logger)
{
  GtUword position,
                nextcc = GT_UNIQUEINT(sainseq->totallength),
                countSstartype = 0,
                *fillptr = sainseq->bucketfillptr;
  GtSainbuffer *sainbuffer = gt_sainbuffer_new(suftab,fillptr,
                                               sainseq->numofchars,logger);
  #{declare_tmpcc(key,"bool nextisStype = true")}

  gt_sain_endbuckets(sainseq);
  for (position = sainseq->totallength-1; /* Nothing */; position--)
  {
    GtUword currentcc = #{getc_call(key,"position")};
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
      printf("Sstar.suftab[" GT_WU "]=" GT_WU "\\n",fillptr[nextcc],position+1);
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

static void gt_sain_#{key}_fast_induceLtypesuffixes1(GtSainseq *sainseq,
                                                 #{getc_param(key)},
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  #{declare_tmpcc(key,"GtWord position")}

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
      currentcc = #{getc_call(key,"(GtUword) position")};
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword t, leftcontextcc;

          gt_assert(position > 0);
          position--;
          leftcontextcc = #{getc_call(key,"(GtUword) position")};
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
          printf("L-induce: suftab[" GT_WU "]=" GT_WD "\\n",
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

static void gt_sain_#{key}_induceLtypesuffixes1(GtSainseq *sainseq,
                                                 #{getc_param(key)},
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  #{declare_tmpcc(key,"GtWord position")}

  gt_assert(sainseq->roundtable == NULL);
  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    if ((position = *suftabptr) > 0)
    {
      GtUword currentcc = #{getc_call(key,"(GtUword) position")};
      if (currentcc < sainseq->numofchars)
      {
        if (position > 0)
        {
          GtUword leftcontextcc;

          gt_assert(position > 0);
          position--;
          leftcontextcc = #{getc_call(key,"(GtUword) position")};
          GT_SAINUPDATEBUCKETPTR(currentcc);
          /* negative => position does not derive L-suffix
             positive => position may derive L-suffix */
          gt_assert(suftabptr < bucketptr);
          *bucketptr++ = (leftcontextcc < currentcc) ? ~position : position;
          *suftabptr = 0;
#ifdef SAINSHOWSTATE
          gt_assert(bucketptr != NULL);
          printf("L-induce: suftab[" GT_WU "]=" GT_WD "\\n",
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

static void gt_sain_#{key}_fast_induceStypesuffixes1(GtSainseq *sainseq,
                                                 #{getc_param(key)},
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  #{declare_tmpcc(key,"GtWord position")}

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
        GtUword currentcc = #{getc_call(key,"(GtUword) position")};

        if (currentcc < sainseq->numofchars)
        {
          GtUword t, leftcontextcc;

          position--;
          leftcontextcc = #{getc_call(key,"(GtUword) position")};
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
          printf("S-induce: suftab[" GT_WU "]=" GT_WD "\\n",
                  (GtUword) (bucketptr - suftab),*bucketptr);
#endif
        }
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_#{key}_induceStypesuffixes1(GtSainseq *sainseq,
                                                 #{getc_param(key)},
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  #{declare_tmpcc(key,"GtWord position")}

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
      GtUword currentcc = #{getc_call(key,"(GtUword) position")};

      if (currentcc < sainseq->numofchars)
      {
        GtUword leftcontextcc;

        position--;
        leftcontextcc = #{getc_call(key,"(GtUword) position")};
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (leftcontextcc > currentcc)
                          ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
        printf("S-induce: suftab[" GT_WU "]=" GT_WD "\\n",
               (GtUword) (bucketptr - suftab),*bucketptr);
#endif
      }
      *suftabptr = 0;
    }
  }
}

static void gt_sain_#{key}_induceLtypesuffixes2(const GtSainseq *sainseq,
                                                 #{getc_param(key)},
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  #{declare_tmpcc(key,"GtWord position")}

  for (suftabptr = suftab; suftabptr < suftab + nonspecialentries; suftabptr++)
  {
    position = *suftabptr;
    *suftabptr = ~position;
    if (position > 0)
    {
      GtUword currentcc;

      position--;
      currentcc = #{getc_call(key,"(GtUword) position")};
      if (currentcc < sainseq->numofchars)
      {
        gt_assert(currentcc > 0);
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && suftabptr < bucketptr);
        *bucketptr++ = (position > 0 &&
                        (#{getc_call(key,"(GtUword) (position-1)")})
                                            < currentcc)
                        ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("L-induce: suftab[" GT_WU "]=" GT_WD "\\n",
               (GtUword) (bucketptr-1-suftab),*(bucketptr-1));
#endif
      }
    }
  }
}

static void gt_sain_#{key}_induceStypesuffixes2(const GtSainseq *sainseq,
                                                 #{getc_param(key)},
                                         GtWord *suftab,
                                         GtUword nonspecialentries)
{
  GtUword lastupdatecc = 0, *fillptr = sainseq->bucketfillptr;
  GtWord *suftabptr, *bucketptr = NULL;
  #{declare_tmpcc(key,"GtWord position")}

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
      currentcc = #{getc_call(key,"(GtUword) position")};
      if (currentcc < sainseq->numofchars)
      {
        GT_SAINUPDATEBUCKETPTR(currentcc);
        gt_assert(bucketptr != NULL && bucketptr - 1 < suftabptr);
        *(--bucketptr) = (position == 0 ||
                          (#{getc_call(key,"(GtUword) (position-1)")})
                                             > currentcc)
                         ? ~position : position;
#ifdef SAINSHOWSTATE
        gt_assert(bucketptr != NULL);
        printf("S-induce: suftab[" GT_WU "]=" GT_WD "\\n",
                (GtUword) (bucketptr-suftab),*bucketptr);
#endif
      }
    } else
    {
      *suftabptr = ~position;
    }
  }
}

static void gt_sain_#{key}_expandorder2original(GtSainseq *sainseq,
                                                 #{getc_param(key)},
                                         GtUword numberofsuffixes,
                                         GtUword *suftab)
{
  GtUword *suftabptr, position,
                writeidx = numberofsuffixes - 1,
                nextcc = GT_UNIQUEINT(sainseq->totallength),
                *sstarsuffixes = suftab + numberofsuffixes;
  GtUword *sstarfirstcharcount = NULL, *bucketsize = NULL;
  #{declare_tmpcc(key,"bool nextisStype = true")}

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
    GtUword currentcc = #{getc_call(key,"position")};
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
CCODE
end
