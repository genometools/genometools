#!/usr/bin/env ruby

require 'ostruct'
require 'optparse'

def usage(opts,msg)
  STDERR.puts "#{$0}: #{msg}\n#{opts.to_s}"
  exit 1
end

def parseargs(argv)
  options = OpenStruct.new
  options.ulong = true
  opts = OptionParser.new
  opts.on("--ulongpair","generate code sorting pairs of ulongs") do |x|
    options.ulong = false
  end
  rest = opts.parse(argv)
  if not rest.empty?
    usage(opts,"superfluous arguments")
  end
  return options
end

def makekey(options)
  if options.ulong 
    return "ulong"
  else
    return "ulongpair"
  end
end

def maketype(options)
  if options.ulong 
    return "unsigned long"
  else
    return "GtUlongPair"
  end
end

def derefptr(ptr,options)
  if options.ulong 
    return "*#{ptr}"
  else
    return "#{ptr}->a"
  end
end

def derefval(val,options)
  if options.ulong 
    return "#{val}"
  else
    return "#{val}.a"
  end
end

options = parseargs(ARGV)

print <<END_OF_FILE
static #{maketype(options)} gt_radixsort_#{makekey(options)}_bin_get(const GtRadixbuffer *rbuf,
                                            unsigned long binnum)
{
  return rbuf->values.#{makekey(options)}ptr[(binnum << rbuf->log_bufsize) +
                             (unsigned long) rbuf->nextidx[binnum]];
}

static void gt_radixsort_#{makekey(options)}_bin_update(#{maketype(options)} *target,
                                    GtRadixbuffer *rbuf,
                                    unsigned long binnum,
                                    #{maketype(options)} value)
{
  unsigned long binoffset = binnum << rbuf->log_bufsize;

  rbuf->values.#{makekey(options)}ptr[binoffset + (unsigned long) rbuf->nextidx[binnum]]
    = value;
  if ((unsigned long) rbuf->nextidx[binnum] < rbuf->buf_size - 1)
  {
    rbuf->nextidx[binnum]++;
  } else
  {
    unsigned long j;
    #{maketype(options)} *wtargetptr, *rtargetptr, *rend, *valptr;

    wtargetptr = target + rbuf->endofbin[binnum] - (rbuf->buf_size - 1);
    rtargetptr = wtargetptr + rbuf->buf_size;
    rend = target + rbuf->startofbin[binnum+1];
    valptr = rbuf->values.#{makekey(options)}ptr + binoffset;
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

static void gt_radixsort_#{makekey(options)}_cached_shuffle(GtRadixbuffer *rbuf,
                                              #{maketype(options)} *source,
                                              GtCountbasetype len,
                                              size_t rightshift)
{
  unsigned long binoffset, binnum, bufoffset,
                nextbin, firstnonemptybin = UINT8_MAX+1;
  GtCountbasetype *count, previouscount, current;
  #{maketype(options)} *sp, *spend = source + len;

  rbuf->countcached++;
  count = rbuf->startofbin; /* use same memory for count and startofbin */
  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    count[binnum] = 0;
    rbuf->nextidx[binnum] = 0;
  }
  for (sp = source; sp < spend; sp++)
  {
    count[GT_RADIX_KEY(UINT8_MAX,rightshift,#{derefptr("sp",options)})]++;
  }
  for (bufoffset = 0, binoffset = 0, binnum = 0; binnum <= UINT8_MAX;
       bufoffset += rbuf->buf_size, binoffset += count[binnum], binnum++)
  {
    unsigned long j;
    const unsigned long end = MIN(rbuf->buf_size,(unsigned long) count[binnum]);

    if (firstnonemptybin == UINT8_MAX+1 && end > 0)
    {
      firstnonemptybin = binnum;
    }
    for (j=0; j<end; j++)
    {
      rbuf->values.#{makekey(options)}ptr[bufoffset + j] = source[binoffset + j];
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
    #{maketype(options)} currentvalue = gt_radixsort_#{makekey(options)}_bin_get(rbuf,binnum);
    while (true)
    {
      binnum = GT_RADIX_KEY(UINT8_MAX,rightshift,#{derefval("currentvalue",options)});
      if (current != rbuf->endofbin[binnum])
      {
        #{maketype(options)} tmp = currentvalue;
        currentvalue = gt_radixsort_#{makekey(options)}_bin_get(rbuf,binnum);
        gt_radixsort_#{makekey(options)}_bin_update(source,rbuf,binnum,tmp);
      } else
      {
        break;
      }
    }
    gt_radixsort_#{makekey(options)}_bin_update(source,rbuf,binnum,currentvalue);
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
    unsigned long bufleft = (unsigned long) rbuf->nextidx[binnum];

    if (bufleft > 0)
    {
      unsigned long j;
      #{maketype(options)} *targetptr, *valptr;

      valptr = rbuf->values.#{makekey(options)}ptr + (binnum << rbuf->log_bufsize);
      targetptr = source + rbuf->startofbin[binnum+1] - bufleft;
      for (j=0; j < bufleft; j++)
      {
        targetptr[j] = valptr[j];
      }
    }
  }
}

static void gt_radixsort_#{makekey(options)}_uncached_shuffle(GtRadixbuffer *rbuf,
                                                #{maketype(options)} *source,
                                                GtCountbasetype len,
                                                size_t rightshift)
{
  unsigned long binnum, nextbin;
  #{maketype(options)} *sp, *spend = source + len;
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
    count[GT_RADIX_KEY(UINT8_MAX,rightshift,#{derefptr("sp",options)})]++;
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
    #{maketype(options)} currentvalue = source[current];
    GtCountbasetype *binptr;

    while (true)
    {
      binptr = rbuf->endofbin +
               GT_RADIX_KEY(UINT8_MAX,rightshift,#{derefval("currentvalue",options)});
      if (current != *binptr)
      {
        #{maketype(options)} tmp = currentvalue;
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

static void gt_radixsort_#{makekey(options)}_shuffle(GtRadixbuffer *rbuf,
                                       #{maketype(options)} *source,
                                       GtCountbasetype len,
                                       size_t rightshift)
{
  ((unsigned long) len > rbuf->cachesize
     ? gt_radixsort_#{makekey(options)}_cached_shuffle
     : gt_radixsort_#{makekey(options)}_uncached_shuffle) (rbuf,source,len,rightshift);
}

static void
gt_radixsort_#{makekey(options)}_inplace_insertionsort(#{maketype(options)} *a,
                                               GtCountbasetype a_size)
{
  #{maketype(options)} *optr, *iptr, *end = a + a_size;

  for (optr = a + 1; optr < end; optr++)
  {
    if (#{derefptr("optr",options)} < #{derefptr("(optr-1)",options)})
    {
      #{maketype(options)} currentElement = *optr;

      *optr = *(optr-1);
      for (iptr = optr-1; iptr > a && #{derefval("currentElement",options)} < #{derefptr("(iptr-1)",options)}; iptr--)
      {
        *iptr = *(iptr-1);
      }
      *iptr = currentElement;
    }
  }
}

static void gt_radixsort_#{makekey(options)}_process_bin(
                                     GtStackGtRadixsort_stackelem *stack,
                                     GtRadixbuffer *rbuf,
                                     #{maketype(options)} *source,
                                     size_t shift)
{
  unsigned long binnum;

  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype width = rbuf->endofbin[binnum] - rbuf->startofbin[binnum];

    if (width >= (GtCountbasetype) 2)
    {
      #{maketype(options)} *ptr = source + rbuf->startofbin[binnum];

      if (width == (GtCountbasetype) 2)
      {
        if (#{derefptr("ptr",options)} > #{derefptr("(ptr+1)",options)})
        {
          #{maketype(options)} tmp = *ptr;
          *ptr = *(ptr+1);
          *(ptr+1) = tmp;
        }
      } else
      {
        if (width <= (GtCountbasetype) 32)
        {
          rbuf->countinsertionsort++;
          gt_radixsort_#{makekey(options)}_inplace_insertionsort(ptr,width);
        } else
        {
          GtRadixsort_stackelem tmpstackelem;

          tmpstackelem.left.#{makekey(options)}ptr = ptr;
          tmpstackelem.len = width;
          tmpstackelem.shift = shift - CHAR_BIT;
          GT_STACK_PUSH(stack,tmpstackelem);
        }
      }
    }
  }
}

static void gt_radixsort_#{makekey(options)}_sub_inplace(GtRadixbuffer *rbuf,
                                           GtStackGtRadixsort_stackelem *stack)
{
  GtRadixsort_stackelem currentstackelem;

  while (!GT_STACK_ISEMPTY(stack))
  {
    currentstackelem = GT_STACK_POP(stack);
    gt_radixsort_#{makekey(options)}_shuffle(rbuf,currentstackelem.left.#{makekey(options)}ptr,
                         currentstackelem.len,
                         currentstackelem.shift);
    if (currentstackelem.shift > 0)
    {
      (void) gt_radixsort_#{makekey(options)}_process_bin(stack,rbuf,
                                            currentstackelem.left.#{makekey(options)}ptr,
                                            currentstackelem.shift);
    }
  }
}
END_OF_FILE
