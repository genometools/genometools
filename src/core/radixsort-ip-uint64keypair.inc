static Gtuint64keyPair gt_radixsort_uint64keypair_bin_get(
                                            const GtRadixbuffer *rbuf,
                                            GtUword binnum)
{
  return rbuf->values.uint64keypairptr[(binnum << rbuf->log_bufsize) +
                             (GtUword) rbuf->nextidx[binnum]];
}

static void gt_radixsort_uint64keypair_bin_update(
                                    Gtuint64keyPair *target,
                                    GtRadixbuffer *rbuf,
                                    GtUword binnum,
                                    Gtuint64keyPair value)
{
  GtUword binoffset = binnum << rbuf->log_bufsize;

  rbuf->values.uint64keypairptr[binoffset +
                                      (GtUword) rbuf->nextidx[binnum]] = value;
  if ((GtUword) rbuf->nextidx[binnum] < rbuf->buf_size - 1)
  {
    rbuf->nextidx[binnum]++;
  } else
  {
    GtUword j;
    Gtuint64keyPair *wtargetptr, *rtargetptr, *rend, *valptr;

    wtargetptr = target + rbuf->endofbin[binnum] - (rbuf->buf_size - 1);
    rtargetptr = wtargetptr + rbuf->buf_size;
    rend = target + rbuf->startofbin[binnum+1];
    valptr = rbuf->values.uint64keypairptr + binoffset;
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

static void gt_radixsort_uint64keypair_cached_shuffle(GtRadixbuffer *rbuf,
                                              Gtuint64keyPair *source,
                                              GtCountbasetype len,
                                              size_t rightshift)
{
  GtUword binoffset, binnum, bufoffset,
                nextbin, firstnonemptybin = UINT8_MAX+1;
  GtCountbasetype *count, previouscount, current;
  Gtuint64keyPair *sp, *spend = source + len;

  rbuf->countcached++;
  count = rbuf->startofbin; /* use same memory for count and startofbin */
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    count[binnum] = 0;
    rbuf->nextidx[binnum] = 0;
  }
  for (sp = source; sp < spend; sp++)
  {
    count[(rightshift > (sizeof (GtUword) - 1) * CHAR_BIT) ?
GT_RADIX_KEY(UINT8_MAX,rightshift - sizeof (GtUword) * CHAR_BIT,
sp->uint64_a) :
GT_RADIX_KEY(UINT8_MAX,rightshift,sp->uint64_b)]++;
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
      rbuf->values.uint64keypairptr[bufoffset + j] =
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
    Gtuint64keyPair currentvalue =
      gt_radixsort_uint64keypair_bin_get(rbuf,binnum);
    while (true)
    {
      binnum = (rightshift > (sizeof (GtUword) - 1) * CHAR_BIT) ?
GT_RADIX_KEY(UINT8_MAX,rightshift - sizeof (GtUword) * CHAR_BIT,
currentvalue.uint64_a) :
GT_RADIX_KEY(UINT8_MAX,rightshift,currentvalue.uint64_b);
      if (current != rbuf->endofbin[binnum])
      {
        Gtuint64keyPair tmp = currentvalue;
        currentvalue = gt_radixsort_uint64keypair_bin_get(rbuf,binnum);
        gt_radixsort_uint64keypair_bin_update(source,rbuf,binnum,tmp);
      } else
      {
        break;
      }
    }
    gt_radixsort_uint64keypair_bin_update(source,rbuf,binnum,
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
      Gtuint64keyPair *targetptr, *valptr;

      valptr =
        rbuf->values.uint64keypairptr + (binnum << rbuf->log_bufsize);
      targetptr = source + rbuf->startofbin[binnum+1] - bufleft;
      for (j=0; j < bufleft; j++)
      {
        targetptr[j] = valptr[j];
      }
    }
  }
}

static void gt_radixsort_uint64keypair_uncached_shuffle(
                       GtRadixbuffer *rbuf,
                       Gtuint64keyPair *source,
                       GtCountbasetype len,
                       size_t rightshift)
{
  GtUword binnum, nextbin;
  Gtuint64keyPair *sp, *spend = source + len;
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
    count[(rightshift > (sizeof (GtUword) - 1) * CHAR_BIT) ?
GT_RADIX_KEY(UINT8_MAX,rightshift - sizeof (GtUword) * CHAR_BIT,
sp->uint64_a) :
GT_RADIX_KEY(UINT8_MAX,rightshift,sp->uint64_b)]++;
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
    Gtuint64keyPair currentvalue = source[current];
    GtCountbasetype *binptr;

    while (true)
    {
      binptr = rbuf->endofbin +
               ((rightshift > (sizeof (GtUword) - 1) * CHAR_BIT) ?
GT_RADIX_KEY(UINT8_MAX,rightshift - sizeof (GtUword) * CHAR_BIT,
currentvalue.uint64_a) :
GT_RADIX_KEY(UINT8_MAX,rightshift,currentvalue.uint64_b));
      if (current != *binptr)
      {
        Gtuint64keyPair tmp = currentvalue;
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

static void gt_radixsort_uint64keypair_shuffle(GtRadixbuffer *rbuf,
                                       Gtuint64keyPair *source,
                                       GtCountbasetype len,
                                       size_t rightshift)
{
  gt_assert(rbuf != NULL);
  if ((GtUword) len > rbuf->cachesize)
  {
    gt_radixsort_uint64keypair_cached_shuffle(rbuf,source,len,rightshift);
  } else
  {
    gt_radixsort_uint64keypair_uncached_shuffle(rbuf,source,len,
                                                      rightshift);
  }
}

static void
gt_radixsort_uint64keypair_inplace_insertionsort(Gtuint64keyPair *a,
                                               GtCountbasetype a_size)
{
  Gtuint64keyPair *optr, *iptr, *end = a + a_size;

  for (optr = a + 1; optr < end; optr++)
  {
    if (gt_radixsort_compare_smaller(optr,(optr-1)))
    {
      Gtuint64keyPair currentElement = *optr;

      *optr = *(optr-1);
      for (iptr = optr-1;
           iptr > a && gt_radixsort_compare_smaller(&currentElement,(iptr-1));
           iptr--)
      {
        *iptr = *(iptr-1);
      }
      *iptr = currentElement;
    }
  }
}

static void gt_radixsort_uint64keypair_process_bin(
                                     GtStackGtRadixsort_stackelem *stack,
                                     GtRadixbuffer *rbuf,
                                     Gtuint64keyPair *source,
                                     size_t shift)
{
  GtUword binnum;

  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype width = rbuf->endofbin[binnum] - rbuf->startofbin[binnum];

    if (width >= (GtCountbasetype) 2)
    {
      Gtuint64keyPair *ptr = source + rbuf->startofbin[binnum];

      if (width == (GtCountbasetype) 2)
      {
        if (gt_radixsort_compare_smaller((ptr+1),ptr))
        {
          Gtuint64keyPair tmp = *ptr;
          *ptr = *(ptr+1);
          *(ptr+1) = tmp;
        }
      } else
      {
        if (width <= (GtCountbasetype) 32)
        {
          rbuf->countinsertionsort++;
          gt_radixsort_uint64keypair_inplace_insertionsort(ptr,width);
        } else
        {
          GtRadixsort_stackelem tmpstackelem;

          tmpstackelem.left.uint64keypairptr = ptr;
          tmpstackelem.len = width;
          tmpstackelem.shift = shift - CHAR_BIT;
          GT_STACK_PUSH(stack,tmpstackelem);
        }
      }
    }
  }
}

static void gt_radixsort_uint64keypair_sub_inplace(GtRadixbuffer *rbuf,
                                           GtStackGtRadixsort_stackelem *stack)
{
  GtRadixsort_stackelem currentstackelem;

  while (!GT_STACK_ISEMPTY(stack))
  {
    currentstackelem = GT_STACK_POP(stack);
    gt_radixsort_uint64keypair_shuffle(rbuf,
                         currentstackelem.left.uint64keypairptr,
                         currentstackelem.len,
                         currentstackelem.shift);
    if (currentstackelem.shift > 0)
    {
      (void) gt_radixsort_uint64keypair_process_bin(stack,rbuf,
                                   currentstackelem.left.uint64keypairptr,
                                   currentstackelem.shift);
    }
  }
}
