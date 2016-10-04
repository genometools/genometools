static GtUword gt_radixsort_ulong_bin_get(
                                            const GtRadixbuffer *rbuf,
                                            GtUword binnum)
{
  return rbuf->values.ulongptr[
                 (binnum << rbuf->log_bufsize) +
                 (GtUword) rbuf->nextidx[binnum]];
}
static void gt_radixsort_ulong_bin_update(
                                    GtUword *source,
                                    GtRadixbuffer *rbuf,
                                    GtUword binnum,
                                    GtUword value)
{
  GtUword binoffset = binnum << rbuf->log_bufsize;

  rbuf->values.ulongptr
[binoffset + (GtUword) rbuf->nextidx[binnum]]=
value;
  if ((GtUword) rbuf->nextidx[binnum] < rbuf->buf_size - 1)
  {
    rbuf->nextidx[binnum]++;
  } else
  {
    GtUword j;
    GtUword *wsourceptr, *rsourceptr, *rend, *valptr;

    wsourceptr = source +
                 (rbuf->endofbin[binnum] - (rbuf->buf_size - 1))
;
    rsourceptr = wsourceptr + rbuf->buf_size;
    rend = source + rbuf->startofbin[binnum+1];
    valptr = rbuf->values.ulongptr +
             binoffset;
    for (j=0; j<rbuf->buf_size; j++)
    {
      *wsourceptr = *valptr;
      wsourceptr++;
      if (rsourceptr < rend)
      {
        *valptr = *rsourceptr;
        rsourceptr++;
      }
      valptr++;
    }
    rbuf->nextidx[binnum] = 0;
  }
  rbuf->endofbin[binnum]++;
}

static void gt_radixsort_ulong_cached_shuffle(GtRadixbuffer *rbuf,
                                              GtUword *source,
                                              GtCountbasetype len,
                                              size_t rightshift)
{
  GtUword binoffset, binnum, bufoffset,
                nextbin, firstnonemptybin = UINT8_MAX+1;
  GtCountbasetype *count, previouscount, currentidx;
  GtUword *sourceptr,
                           *sourceend = source + len;

  rbuf->countcached++;
  count = rbuf->startofbin; /* use same memory for count and startofbin */
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    count[binnum] = 0;
    rbuf->nextidx[binnum] = 0;
  }
  for (sourceptr = source; sourceptr < sourceend; sourceptr++)
  {
    count[GT_RADIX_KEY(UINT8_MAX,rightshift,*sourceptr)]++;
  }
  for (bufoffset = 0, binoffset = 0, binnum = 0; binnum <= UINT8_MAX;
       bufoffset += rbuf->buf_size, binoffset += count[binnum], binnum++)
  {
    const GtUword elems2copy = MIN(rbuf->buf_size,(GtUword) count[binnum]);

    if (elems2copy > 0)
    {
      if (firstnonemptybin == UINT8_MAX+1)
      {
        firstnonemptybin = binnum;
      }
      memcpy(rbuf->values.
             ulongptr + bufoffset,
             source + binoffset,
             (sizeof *source * elems2copy));
    }
  }
  previouscount = count[0];
  rbuf->startofbin[0] = rbuf->endofbin[0] = 0;
  nextbin = 0;
  for (binnum = 1UL; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype temp = rbuf->startofbin[binnum-1] + previouscount;
    previouscount = count[binnum];
    rbuf->startofbin[binnum] = rbuf->endofbin[binnum] = temp;
  }
  /* to simplify compution of bin end */
  rbuf->startofbin[UINT8_MAX+1] = len;
  for (currentidx = 0, binnum = firstnonemptybin;
       currentidx < len; binnum = nextbin - 1)
  {
    GtUword tmpvalue;
    tmpvalue = gt_radixsort_ulong_bin_get(rbuf,binnum);
    while (true)
    {
      binnum = GT_RADIX_KEY(UINT8_MAX,rightshift,tmpvalue);
      if (currentidx != rbuf->endofbin[binnum])
      {
        GtUword tmpswap;
        tmpswap = tmpvalue;
        tmpvalue = gt_radixsort_ulong_bin_get(rbuf,binnum);
        gt_radixsort_ulong_bin_update
                             (source,rbuf,binnum,
                              tmpswap);
      } else
      {
        break;
      }
    }
    gt_radixsort_ulong_bin_update(source,rbuf,binnum,
                                           tmpvalue);
    currentidx++;
    /* skip over empty bins */
    while (nextbin <= UINT8_MAX && currentidx >= rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    /* skip over full bins */
    while (nextbin <= UINT8_MAX &&
           rbuf->endofbin[nextbin-1] == rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    if (currentidx < rbuf->endofbin[nextbin-1])
    {
      currentidx = rbuf->endofbin[nextbin-1];
    }
  }
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    GtUword bufleft = (GtUword) rbuf->nextidx[binnum];

    if (bufleft > 0)
    {
      GtUword *sourceptr, *valptr;

      valptr = rbuf->values.ulongptr +
               (binnum << rbuf->log_bufsize);
      sourceptr = source +
                  (rbuf->startofbin[binnum+1] - bufleft);
      memcpy(sourceptr,valptr,(sizeof *sourceptr * bufleft));
    }
  }
}

static void gt_radixsort_ulong_uncached_shuffle(
                       GtRadixbuffer *rbuf,
                       GtUword *source,
                       GtCountbasetype len,
                       size_t rightshift)
{
  GtUword binnum, nextbin;
  GtCountbasetype currentidx, previouscount, *count;
  GtUword *sourceptr,
                           *sourceend = source + len;

  rbuf->countuncached++;
  count = rbuf->startofbin; /* use same memory for count and startofbin */
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    count[binnum] = 0;
    rbuf->nextidx[binnum] = 0;
  }
  for (sourceptr = source; sourceptr < sourceend; sourceptr++)
  {
    count[GT_RADIX_KEY(UINT8_MAX,rightshift,*sourceptr)]++;
  }
  previouscount = count[0];
  rbuf->startofbin[0] = rbuf->endofbin[0] = 0;
  nextbin = 0;
  for (binnum = 1UL; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype temp = rbuf->startofbin[binnum-1] + previouscount;
    previouscount = count[binnum];
    rbuf->startofbin[binnum] = rbuf->endofbin[binnum] = temp;
  }
  /* to simplify compution of bin end */
  rbuf->startofbin[UINT8_MAX+1] = len;
  for (currentidx = 0; currentidx < len; /* Nothing */)
  {
    GtCountbasetype *binptr;
    GtUword tmpvalue;
    tmpvalue = source[currentidx];

    while (true)
    {
      binptr = rbuf->endofbin +
               (GT_RADIX_KEY(UINT8_MAX,rightshift,tmpvalue));
      binnum = *binptr;
      if (currentidx != binnum)
      {
        GtUword tmpswap;
        tmpswap = tmpvalue;
        tmpvalue = source[binnum];
        source[binnum] = tmpswap;
        (*binptr)++;
      } else
      {
        break;
      }
    }
    source[binnum] = tmpvalue;
    currentidx++;
    (*binptr)++;
    /* skip over empty bins */
    while (nextbin <= UINT8_MAX && currentidx >= rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    /* skip over full bins */
    while (nextbin <= UINT8_MAX &&
           rbuf->endofbin[nextbin-1] == rbuf->startofbin[nextbin])
    {
      nextbin++;
    }
    if (currentidx < rbuf->endofbin[nextbin-1])
    {
      currentidx = rbuf->endofbin[nextbin-1];
    }
  }
}

static void gt_radixsort_ulong_shuffle(GtRadixbuffer *rbuf,
                                       GtUword *source,
                                       GtCountbasetype len,
                                       size_t rightshift)
{
  gt_assert(rbuf != NULL);
  if ((GtUword) len > rbuf->cachesize)
  {
    gt_radixsort_ulong_cached_shuffle(rbuf,source,len,rightshift);
  } else
  {
    gt_radixsort_ulong_uncached_shuffle(rbuf,source,len,
                                                      rightshift);
  }
}

static void
gt_radixsort_ulong_inplace_insertionsort(
                                  GT_UNUSED GtRadixbuffer *rbuf,
                                  GtUword *arr,
                                  GtCountbasetype a_size)
{
  GtUword *optr,
                           *end = arr + a_size;

  for (optr = arr + 1; optr < end;
       optr++)
  {
    GtUword *oprevious = optr - 1;

    if (*optr < *oprevious)
    {
      GtUword *iptr;
      GtUword tmpvalue;
      tmpvalue = *optr;

      *optr = *oprevious;
      for (iptr = oprevious; iptr > arr; iptr -= 1)
      {
        GtUword *iprevious = iptr - 1;
        if (!(tmpvalue < *iprevious))
        {
          break;
        }
        *iptr = *iprevious;
      }
      *iptr = tmpvalue;
    }
  }
}

static void gt_radixsort_ulong_process_bin(
                                     GtStackGtRadixsort_stackelem *stack,
                                     GtRadixbuffer *rbuf,
                                     GtUword *source,
                                     size_t shift)
{
  GtUword binnum;

  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype width = rbuf->endofbin[binnum] - rbuf->startofbin[binnum];

    if (width >= (GtCountbasetype) 2)
    {
      GtUword *ptr
       = source + rbuf->startofbin[binnum];

      if (width == (GtCountbasetype) 2)
      {
        GtUword *nextptr = ptr + 1;
        if (*nextptr < *ptr)
        {
          GtUword tmpswap;
          tmpswap = *ptr;
          *ptr = *nextptr;
          *nextptr = tmpswap;
        }
      } else
      {
        if (width <= (GtCountbasetype) 32)
        {
          rbuf->countinsertionsort++;
          gt_radixsort_ulong_inplace_insertionsort(rbuf,ptr,width);
        } else
        {
          GtRadixsort_stackelem tmpstackelem;

          tmpstackelem.left.ulongptr = ptr;
          tmpstackelem.len = width;
          tmpstackelem.shift = shift - CHAR_BIT;
          GT_STACK_PUSH(stack,tmpstackelem);
        }
      }
    }
  }
}

static void gt_radixsort_ulong_sub_inplace(GtRadixbuffer *rbuf,
                                           GtStackGtRadixsort_stackelem *stack)
{
  GtRadixsort_stackelem currentstackelem;

  while (!GT_STACK_ISEMPTY(stack))
  {
    currentstackelem = GT_STACK_POP(stack);
    gt_radixsort_ulong_shuffle(rbuf,
                         currentstackelem.left.ulongptr,
                         currentstackelem.len,
                         currentstackelem.shift);
    if (currentstackelem.shift > 0)
    {
      (void) gt_radixsort_ulong_process_bin(stack,rbuf,
                                   currentstackelem.left.ulongptr,
                                   currentstackelem.shift);
    }
  }
}
