static uint8_t *gt_radixsort_flba_bin_get(
                                            const GtRadixbuffer *rbuf,
                                            GtUword binnum)
{
  return rbuf->values.flbaptr +
                 ((binnum << rbuf->log_bufsize) +
                  (GtUword) rbuf->nextidx[binnum]) * rbuf->flba_unitsize;
}
static inline void gt_radixsort_flba_bin_update(
                                    uint8_t *source,
                                    GtRadixbuffer *rbuf,
                                    GtUword binnum,
                                    const uint8_t* value)
{
  GtUword binoffset = binnum << rbuf->log_bufsize;

  memcpy(rbuf->values.flbaptr
 + (binoffset + (GtUword) rbuf->nextidx[binnum]) * rbuf->flba_unitsize,
value,rbuf->flba_unitsize);
  if ((GtUword) rbuf->nextidx[binnum] < rbuf->buf_size - 1)
  {
    rbuf->nextidx[binnum]++;
  } else
  {
    GtUword j;
    uint8_t *wsourceptr, *rsourceptr, *rend, *valptr;

    wsourceptr = source +
                 (rbuf->endofbin[binnum] - (rbuf->buf_size - 1))
 * rbuf->flba_unitsize;
    rsourceptr = wsourceptr + rbuf->buf_size * rbuf->flba_unitsize;
    rend = source + rbuf->startofbin[binnum+1] * rbuf->flba_unitsize;
    valptr = rbuf->values.flbaptr +
             binoffset * rbuf->flba_unitsize;
    for (j=0; j<rbuf->buf_size; j++)
    {
      memcpy(wsourceptr,valptr,rbuf->flba_unitsize);
      wsourceptr += rbuf->flba_unitsize;
      if (rsourceptr < rend)
      {
        memcpy(valptr,rsourceptr,rbuf->flba_unitsize);
        rsourceptr += rbuf->flba_unitsize;
      }
      valptr += rbuf->flba_unitsize;
    }
    rbuf->nextidx[binnum] = 0;
  }
  rbuf->endofbin[binnum]++;
}

static void gt_radixsort_flba_cached_shuffle(GtRadixbuffer *rbuf,
                                              uint8_t *source,
                                              GtCountbasetype len,
                                              size_t rightshift)
{
  GtUword binoffset, binnum, bufoffset,
                nextbin, firstnonemptybin = UINT8_MAX+1;
  GtCountbasetype *count, previouscount, currentidx;
  uint8_t *sourceptr,
                           *sourceend = source + len * rbuf->flba_unitsize;

  rbuf->countcached++;
  count = rbuf->startofbin; /* use same memory for count and startofbin */
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    count[binnum] = 0;
    rbuf->nextidx[binnum] = 0;
  }
  for (sourceptr = source; sourceptr < sourceend;
       sourceptr += rbuf->flba_unitsize)
  {
    count[sourceptr[rightshift]]++;
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
             flbaptr + bufoffset * rbuf->flba_unitsize,
             source + binoffset * rbuf->flba_unitsize,
             (sizeof *source * elems2copy) * rbuf->flba_unitsize);
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
    /* no decl. */
    memcpy(rbuf->tmpvalue_ptr,gt_radixsort_flba_bin_get(rbuf,binnum),
rbuf->flba_unitsize);
    while (true)
    {
      binnum = rbuf->tmpvalue_ptr[rightshift];
      if (currentidx != rbuf->endofbin[binnum])
      {
        /* no decl. */
        memcpy(rbuf->tmpswap_ptr,rbuf->tmpvalue_ptr,
rbuf->flba_unitsize);
        memcpy(rbuf->tmpvalue_ptr,gt_radixsort_flba_bin_get(rbuf,binnum),
rbuf->flba_unitsize);
        gt_radixsort_flba_bin_update
                             (source,rbuf,binnum,
                              rbuf->tmpswap_ptr);
      } else
      {
        break;
      }
    }
    gt_radixsort_flba_bin_update(source,rbuf,binnum,
                                           rbuf->tmpvalue_ptr);
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
      uint8_t *sourceptr, *valptr;

      valptr = rbuf->values.flbaptr +
               (binnum << rbuf->log_bufsize) * rbuf->flba_unitsize;
      sourceptr = source +
                  (rbuf->startofbin[binnum+1] - bufleft) * rbuf->flba_unitsize;
      memcpy(sourceptr,valptr,
             (sizeof *sourceptr * bufleft) * rbuf->flba_unitsize);
    }
  }
}

static void gt_radixsort_flba_uncached_shuffle(
                       GtRadixbuffer *rbuf,
                       uint8_t *source,
                       GtCountbasetype len,
                       size_t rightshift)
{
  GtUword binnum, nextbin;
  GtCountbasetype currentidx, previouscount, *count;
  uint8_t *sourceptr,
                           *sourceend = source + len * rbuf->flba_unitsize;

  rbuf->countuncached++;
  count = rbuf->startofbin; /* use same memory for count and startofbin */
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    count[binnum] = 0;
    rbuf->nextidx[binnum] = 0;
  }
  for (sourceptr = source; sourceptr < sourceend;
       sourceptr += rbuf->flba_unitsize)
  {
    count[sourceptr[rightshift]]++;
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
    /* no decl. */
    memcpy(rbuf->tmpvalue_ptr,source + (currentidx) * rbuf->flba_unitsize,
rbuf->flba_unitsize);

    while (true)
    {
      binptr = rbuf->endofbin +
               (rbuf->tmpvalue_ptr[rightshift]);
      binnum = *binptr;
      if (currentidx != binnum)
      {
        /* no decl. */
        memcpy(rbuf->tmpswap_ptr,rbuf->tmpvalue_ptr,
rbuf->flba_unitsize);
        memcpy(rbuf->tmpvalue_ptr,source + (binnum) * rbuf->flba_unitsize,
rbuf->flba_unitsize);
        memcpy(source + (binnum) * rbuf->flba_unitsize,
rbuf->tmpswap_ptr,rbuf->flba_unitsize);
        (*binptr)++;
      } else
      {
        break;
      }
    }
    memcpy(source + (binnum) * rbuf->flba_unitsize,
rbuf->tmpvalue_ptr,rbuf->flba_unitsize);
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

static void gt_radixsort_flba_shuffle(GtRadixbuffer *rbuf,
                                       uint8_t *source,
                                       GtCountbasetype len,
                                       size_t rightshift)
{
  gt_assert(rbuf != NULL);
  if ((GtUword) len > rbuf->cachesize)
  {
    gt_radixsort_flba_cached_shuffle(rbuf,source,len,rightshift);
  } else
  {
    gt_radixsort_flba_uncached_shuffle(rbuf,source,len,
                                                      rightshift);
  }
}

static void
gt_radixsort_flba_inplace_insertionsort(
                                  GT_UNUSED GtRadixbuffer *rbuf,
                                  uint8_t *arr,
                                  GtCountbasetype a_size)
{
  uint8_t *optr,
                           *end = arr + a_size * rbuf->flba_unitsize;

  for (optr = arr + 1 * rbuf->flba_unitsize; optr < end;
       optr += rbuf->flba_unitsize)
  {
    uint8_t *oprevious = optr - 1 * rbuf->flba_unitsize;

    if (memcmp(optr,oprevious,rbuf->flba_unitsize) < 0)
    {
      uint8_t *iptr;
      /* no decl. */
      memcpy(rbuf->tmpvalue_ptr,optr,
rbuf->flba_unitsize);

      memcpy(optr,oprevious,rbuf->flba_unitsize);
      for (iptr = oprevious; iptr > arr; iptr -= 1 * rbuf->flba_unitsize)
      {
        uint8_t *iprevious = iptr - 1 * rbuf->flba_unitsize;
        if (!(memcmp(rbuf->tmpvalue_ptr,iprevious,rbuf->flba_unitsize) < 0))
        {
          break;
        }
        memcpy(iptr,iprevious,rbuf->flba_unitsize);
      }
      memcpy(iptr,rbuf->tmpvalue_ptr,
rbuf->flba_unitsize);
    }
  }
}

static void gt_radixsort_flba_process_bin(
                                     GtStackGtRadixsort_stackelem *stack,
                                     GtRadixbuffer *rbuf,
                                     uint8_t *source,
                                     size_t shift)
{
  GtUword binnum;

  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype width = rbuf->endofbin[binnum] - rbuf->startofbin[binnum];

    if (width >= (GtCountbasetype) 2)
    {
      uint8_t *ptr
       = source + rbuf->startofbin[binnum] * rbuf->flba_unitsize;

      if (width == (GtCountbasetype) 2)
      {
        uint8_t *nextptr = ptr + 1 * rbuf->flba_unitsize;
        if (memcmp(nextptr,ptr,rbuf->flba_unitsize) < 0)
        {
          /* no decl. */
          memcpy(rbuf->tmpswap_ptr,ptr,
rbuf->flba_unitsize);
          memcpy(ptr,nextptr,rbuf->flba_unitsize);
          memcpy(nextptr,rbuf->tmpswap_ptr,
rbuf->flba_unitsize);
        }
      } else
      {
        if (width <= (GtCountbasetype) 32)
        {
          rbuf->countinsertionsort++;
          gt_radixsort_flba_inplace_insertionsort(rbuf,ptr,width);
        } else
        {
          GtRadixsort_stackelem tmpstackelem;

          tmpstackelem.left.flbaptr = ptr;
          tmpstackelem.len = width;
          tmpstackelem.shift = shift+1;
          GT_STACK_PUSH(stack,tmpstackelem);
        }
      }
    }
  }
}

static void gt_radixsort_flba_sub_inplace(GtRadixbuffer *rbuf,
                                           GtStackGtRadixsort_stackelem *stack)
{
  GtRadixsort_stackelem currentstackelem;

  while (!GT_STACK_ISEMPTY(stack))
  {
    currentstackelem = GT_STACK_POP(stack);
    gt_radixsort_flba_shuffle(rbuf,
                         currentstackelem.left.flbaptr,
                         currentstackelem.len,
                         currentstackelem.shift);
    if (currentstackelem.shift < rbuf->sort_bytes-1)
    {
      (void) gt_radixsort_flba_process_bin(stack,rbuf,
                                   currentstackelem.left.flbaptr,
                                   currentstackelem.shift);
    }
  }
}
