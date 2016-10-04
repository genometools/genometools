static Gtuint64keyPair gt_radixsort_uint64keypair_bin_get(
                                            const GtRadixbuffer *rbuf,
                                            GtUword binnum)
{
  return rbuf->values.uint64keypairptr[
                 (binnum << rbuf->log_bufsize) +
                 (GtUword) rbuf->nextidx[binnum]];
}
static inline void gt_radixsort_uint64keypair_bin_update(
                                    Gtuint64keyPair *source,
                                    GtRadixbuffer *rbuf,
                                    GtUword binnum,
                                    Gtuint64keyPair value)
{
  GtUword binoffset = binnum << rbuf->log_bufsize;

  rbuf->values.uint64keypairptr
[binoffset + (GtUword) rbuf->nextidx[binnum]]=
value;
  if ((GtUword) rbuf->nextidx[binnum] < rbuf->buf_size - 1)
  {
    rbuf->nextidx[binnum]++;
  } else
  {
    GtUword j;
    Gtuint64keyPair *wsourceptr, *rsourceptr, *rend, *valptr;

    wsourceptr = source +
                 (rbuf->endofbin[binnum] - (rbuf->buf_size - 1))
;
    rsourceptr = wsourceptr + rbuf->buf_size;
    rend = source + rbuf->startofbin[binnum+1];
    valptr = rbuf->values.uint64keypairptr +
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

static void gt_radixsort_uint64keypair_cached_shuffle(GtRadixbuffer *rbuf,
                                              Gtuint64keyPair *source,
                                              GtCountbasetype len,
                                              size_t rightshift)
{
  GtUword binoffset, binnum, bufoffset,
                nextbin, firstnonemptybin = UINT8_MAX+1;
  GtCountbasetype *count, previouscount, currentidx;
  Gtuint64keyPair *sourceptr,
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
    count[(rightshift > (sizeof (GtUword) - 1) * CHAR_BIT) ?
GT_RADIX_KEY(UINT8_MAX,rightshift - sizeof (GtUword) * CHAR_BIT,
sourceptr->uint64_a) :
GT_RADIX_KEY(UINT8_MAX,rightshift,sourceptr->uint64_b)]++;
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
             uint64keypairptr + bufoffset,
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
    Gtuint64keyPair tmpvalue;
    tmpvalue = gt_radixsort_uint64keypair_bin_get(rbuf,binnum);
    while (true)
    {
      binnum = (rightshift > (sizeof (GtUword) - 1) * CHAR_BIT) ?
GT_RADIX_KEY(UINT8_MAX,rightshift - sizeof (GtUword) * CHAR_BIT,
tmpvalue.uint64_a) :
GT_RADIX_KEY(UINT8_MAX,rightshift,tmpvalue.uint64_b);
      if (currentidx != rbuf->endofbin[binnum])
      {
        Gtuint64keyPair tmpswap;
        tmpswap = tmpvalue;
        tmpvalue = gt_radixsort_uint64keypair_bin_get(rbuf,binnum);
        gt_radixsort_uint64keypair_bin_update
                             (source,rbuf,binnum,
                              tmpswap);
      } else
      {
        break;
      }
    }
    gt_radixsort_uint64keypair_bin_update(source,rbuf,binnum,
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
      Gtuint64keyPair *sourceptr, *valptr;

      valptr = rbuf->values.uint64keypairptr +
               (binnum << rbuf->log_bufsize);
      sourceptr = source +
                  (rbuf->startofbin[binnum+1] - bufleft);
      memcpy(sourceptr,valptr,(sizeof *sourceptr * bufleft));
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
  GtCountbasetype currentidx, previouscount, *count;
  Gtuint64keyPair *sourceptr,
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
    count[(rightshift > (sizeof (GtUword) - 1) * CHAR_BIT) ?
GT_RADIX_KEY(UINT8_MAX,rightshift - sizeof (GtUword) * CHAR_BIT,
sourceptr->uint64_a) :
GT_RADIX_KEY(UINT8_MAX,rightshift,sourceptr->uint64_b)]++;
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
    Gtuint64keyPair tmpvalue;
    tmpvalue = source[currentidx];

    while (true)
    {
      binptr = rbuf->endofbin +
               ((rightshift > (sizeof (GtUword) - 1) * CHAR_BIT) ?
GT_RADIX_KEY(UINT8_MAX,rightshift - sizeof (GtUword) * CHAR_BIT,
tmpvalue.uint64_a) :
GT_RADIX_KEY(UINT8_MAX,rightshift,tmpvalue.uint64_b));
      binnum = *binptr;
      if (currentidx != binnum)
      {
        Gtuint64keyPair tmpswap;
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
gt_radixsort_uint64keypair_inplace_insertionsort(
                                  GT_UNUSED GtRadixbuffer *rbuf,
                                  Gtuint64keyPair *arr,
                                  GtCountbasetype a_size)
{
  Gtuint64keyPair *optr,
                           *end = arr + a_size;

  for (optr = arr + 1; optr < end;
       optr++)
  {
    Gtuint64keyPair *oprevious = optr - 1;

    if (gt_radixsort_uint64keypair_smaller(optr,oprevious))
    {
      Gtuint64keyPair *iptr;
      Gtuint64keyPair tmpvalue;
      tmpvalue = *optr;

      *optr = *oprevious;
      for (iptr = oprevious; iptr > arr; iptr -= 1)
      {
        Gtuint64keyPair *iprevious = iptr - 1;
        if (!(gt_radixsort_uint64keypair_smaller(&tmpvalue,iprevious)))
        {
          break;
        }
        *iptr = *iprevious;
      }
      *iptr = tmpvalue;
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
      Gtuint64keyPair *ptr
       = source + rbuf->startofbin[binnum];

      if (width == (GtCountbasetype) 2)
      {
        Gtuint64keyPair *nextptr = ptr + 1;
        if (gt_radixsort_uint64keypair_smaller(nextptr,ptr))
        {
          Gtuint64keyPair tmpswap;
          tmpswap = *ptr;
          *ptr = *nextptr;
          *nextptr = tmpswap;
        }
      } else
      {
        if (width <= (GtCountbasetype) 32)
        {
          rbuf->countinsertionsort++;
          gt_radixsort_uint64keypair_inplace_insertionsort(rbuf,ptr,width);
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
