static GtQuasiSeedExtendSeedPair gt_radixsort_seedpair_bin_get(
                                            const GtRadixbuffer *rbuf,
                                            GtUword binnum)
{
  return rbuf->values.seedpairptr[(binnum << rbuf->log_bufsize) +
                             (GtUword) rbuf->nextidx[binnum]];
}

static void gt_radixsort_seedpair_bin_update(
                                    GtQuasiSeedExtendSeedPair *target,
                                    GtRadixbuffer *rbuf,
                                    GtUword binnum,
                                    GtQuasiSeedExtendSeedPair value)
{
  GtUword binoffset = binnum << rbuf->log_bufsize;

  rbuf->values.seedpairptr[binoffset +
                                      (GtUword) rbuf->nextidx[binnum]] = value;
  if ((GtUword) rbuf->nextidx[binnum] < rbuf->buf_size - 1)
  {
    rbuf->nextidx[binnum]++;
  } else
  {
    GtUword j;
    GtQuasiSeedExtendSeedPair *wtargetptr, *rtargetptr, *rend, *valptr;

    wtargetptr = target + rbuf->endofbin[binnum] - (rbuf->buf_size - 1);
    rtargetptr = wtargetptr + rbuf->buf_size;
    rend = target + rbuf->startofbin[binnum+1];
    valptr = rbuf->values.seedpairptr + binoffset;
    for (j=0; j<rbuf->buf_size; j++)
    {
      *wtargetptr++ = *valptr;
      if (rtargetptr < rend)
      {
        *valptr = *rtargetptr++;
      }
      valptr++;
    }
    rbuf->nextidx[binnum] = 0;
  }
  rbuf->endofbin[binnum]++;
}

static void gt_radixsort_seedpair_cached_shuffle(GtRadixbuffer *rbuf,
                                              GtQuasiSeedExtendSeedPair *source,
                                              GtCountbasetype len,
                                              size_t rightshift)
{
  GtUword binoffset, binnum, bufoffset,
                nextbin, firstnonemptybin = UINT8_MAX+1;
  GtCountbasetype *count, previouscount, current;
  GtQuasiSeedExtendSeedPair *sp, *spend = source + len;

  rbuf->countcached++;
  count = rbuf->startofbin; /* use same memory for count and startofbin */
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    count[binnum] = 0;
    rbuf->nextidx[binnum] = 0;
  }
  for (sp = source; sp < spend; sp++)
  {
    count[GT_RADIX_KEY(UINT8_MAX,rightshift,sp->a)]++;
  }
  for (bufoffset = 0, binoffset = 0, binnum = 0; binnum <= UINT8_MAX;
       bufoffset += rbuf->buf_size, binoffset += count[binnum], binnum++)
  {
    GtUword j;
    const GtUword end = MIN(rbuf->buf_size,(GtUword) count[binnum]);

    if (firstnonemptybin == UINT8_MAX+1 && end > 0)
    {
      firstnonemptybin = binnum;
    }
    for (j=0; j<end; j++)
    {
      rbuf->values.seedpairptr[bufoffset + j] =
        source[binoffset + j];
    }
  }
  previouscount = count[0];
  rbuf->startofbin[0] = rbuf->endofbin[0] = 0;
  nextbin = 0;
  for (binnum = 1UL; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype tmp = rbuf->startofbin[binnum-1] + previouscount;
    previouscount = count[binnum];
    rbuf->startofbin[binnum] = rbuf->endofbin[binnum] = tmp;
  }
  /* to simplify compution of bin end */
  rbuf->startofbin[UINT8_MAX+1] = len;
  for (current = 0, binnum = firstnonemptybin;
       current < len; binnum = nextbin - 1)
  {
    GtQuasiSeedExtendSeedPair currentvalue =
      gt_radixsort_seedpair_bin_get(rbuf,binnum);
    while (true)
    {
      binnum = GT_RADIX_KEY(UINT8_MAX,rightshift,
                            currentvalue.a);
      if (current != rbuf->endofbin[binnum])
      {
        GtQuasiSeedExtendSeedPair tmp = currentvalue;
        currentvalue = gt_radixsort_seedpair_bin_get(rbuf,binnum);
        gt_radixsort_seedpair_bin_update(source,rbuf,binnum,tmp);
      } else
      {
        break;
      }
    }
    gt_radixsort_seedpair_bin_update(source,rbuf,binnum,
                                                currentvalue);
    current++;
    /* skip over empty bins */
    while (nextbin <= UINT8_MAX && current >= rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    /* skip over full bins */
    while (nextbin <= UINT8_MAX &&
           rbuf->endofbin[nextbin-1] == rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    if (current < rbuf->endofbin[nextbin-1])
    {
      current = rbuf->endofbin[nextbin-1];
    }
  }
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    GtUword bufleft = (GtUword) rbuf->nextidx[binnum];

    if (bufleft > 0)
    {
      GtUword j;
      GtQuasiSeedExtendSeedPair *targetptr, *valptr;

      valptr =
        rbuf->values.seedpairptr + (binnum << rbuf->log_bufsize);
      targetptr = source + rbuf->startofbin[binnum+1] - bufleft;
      for (j=0; j < bufleft; j++)
      {
        targetptr[j] = valptr[j];
      }
    }
  }
}

static void gt_radixsort_seedpair_uncached_shuffle(
                                                GtRadixbuffer *rbuf,
                                                GtQuasiSeedExtendSeedPair *source,
                                                GtCountbasetype len,
                                                size_t rightshift)
{
  GtUword binnum, nextbin;
  GtQuasiSeedExtendSeedPair *sp, *spend = source + len;
  GtCountbasetype current, previouscount, *count;

  rbuf->countuncached++;
  count = rbuf->startofbin; /* use same memory for count and startofbin */
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    count[binnum] = 0;
    rbuf->nextidx[binnum] = 0;
  }
  for (sp = source; sp < spend; sp++)
  {
    count[GT_RADIX_KEY(UINT8_MAX,rightshift,sp->a)]++;
  }
  previouscount = count[0];
  rbuf->startofbin[0] = rbuf->endofbin[0] = 0;
  nextbin = 0;
  for (binnum = 1UL; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype tmp = rbuf->startofbin[binnum-1] + previouscount;
    previouscount = count[binnum];
    rbuf->startofbin[binnum] = rbuf->endofbin[binnum] = tmp;
  }
  /* to simplify compution of bin end */
  rbuf->startofbin[UINT8_MAX+1] = len;
  for (current = 0; current < len; /* Nothing */)
  {
    GtQuasiSeedExtendSeedPair currentvalue = source[current];
    GtCountbasetype *binptr;

    while (true)
    {
      binptr = rbuf->endofbin +
               GT_RADIX_KEY(UINT8_MAX,rightshift,
                            currentvalue.a);
      if (current != *binptr)
      {
        GtQuasiSeedExtendSeedPair tmp = currentvalue;
        currentvalue = source[*binptr];
        source[*binptr] = tmp;
        (*binptr)++;
      } else
      {
        break;
      }
    }
    source[current++] = currentvalue;
    (*binptr)++;
    /* skip over empty bins */
    while (nextbin <= UINT8_MAX && current >= rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    /* skip over full bins */
    while (nextbin <= UINT8_MAX &&
           rbuf->endofbin[nextbin-1] == rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    if (current < rbuf->endofbin[nextbin-1])
    {
      current = rbuf->endofbin[nextbin-1];
    }
  }
}

static void gt_radixsort_seedpair_shuffle(GtRadixbuffer *rbuf,
                                       GtQuasiSeedExtendSeedPair *source,
                                       GtCountbasetype len,
                                       size_t rightshift)
{
  gt_assert(rbuf != NULL);
  if ((GtUword) len > rbuf->cachesize)
  {
    gt_radixsort_seedpair_cached_shuffle(rbuf,source,len,rightshift);
  } else
  {
    gt_radixsort_seedpair_uncached_shuffle(rbuf,source,len,
                                                      rightshift);
  }
}

static void
gt_radixsort_seedpair_inplace_insertionsort(GtQuasiSeedExtendSeedPair *a,
                                               GtCountbasetype a_size)
{
  GtQuasiSeedExtendSeedPair *optr, *iptr, *end = a + a_size;

  for (optr = a + 1; optr < end; optr++)
  {
    if (optr->a < (optr-1)->a)
    {
      GtQuasiSeedExtendSeedPair currentElement = *optr;

      *optr = *(optr-1);
      for (iptr = optr-1;
           iptr > a &&
             currentElement.a <
             (iptr-1)->a;
           iptr--)
      {
        *iptr = *(iptr-1);
      }
      *iptr = currentElement;
    }
  }
}

static void gt_radixsort_seedpair_process_bin(
                                     GtStackGtRadixsort_stackelem *stack,
                                     GtRadixbuffer *rbuf,
                                     GtQuasiSeedExtendSeedPair *source,
                                     size_t shift)
{
  GtUword binnum;

  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype width = rbuf->endofbin[binnum] - rbuf->startofbin[binnum];

    if (width >= (GtCountbasetype) 2)
    {
      GtQuasiSeedExtendSeedPair *ptr = source + rbuf->startofbin[binnum];

      if (width == (GtCountbasetype) 2)
      {
        if ((ptr+1)->a < ptr->a)
        {
          GtQuasiSeedExtendSeedPair tmp = *ptr;
          *ptr = *(ptr+1);
          *(ptr+1) = tmp;
        }
      } else
      {
        if (width <= (GtCountbasetype) 32)
        {
          rbuf->countinsertionsort++;
          gt_radixsort_seedpair_inplace_insertionsort(ptr,width);
        } else
        {
          GtRadixsort_stackelem tmpstackelem;

          tmpstackelem.left.seedpairptr = ptr;
          tmpstackelem.len = width;
          tmpstackelem.shift = shift - CHAR_BIT;
          GT_STACK_PUSH(stack,tmpstackelem);
        }
      }
    }
  }
}

static void gt_radixsort_seedpair_sub_inplace(GtRadixbuffer *rbuf,
                                           GtStackGtRadixsort_stackelem *stack)
{
  GtRadixsort_stackelem currentstackelem;

  while (!GT_STACK_ISEMPTY(stack))
  {
    currentstackelem = GT_STACK_POP(stack);
    gt_radixsort_seedpair_shuffle(rbuf,
                         currentstackelem.left.seedpairptr,
                         currentstackelem.len,
                         currentstackelem.shift);
    if (currentstackelem.shift > 0)
    {
      (void) gt_radixsort_seedpair_process_bin(stack,rbuf,
                                   currentstackelem.left.seedpairptr,
                                   currentstackelem.shift);
    }
  }
}
