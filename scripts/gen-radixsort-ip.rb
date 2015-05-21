#!/usr/bin/env ruby

require 'ostruct'
require 'optparse'

def usage(opts,msg)
  STDERR.puts "#{$0}: #{msg}\n#{opts.to_s}"
  exit 1
end

def parseargs(argv)
  options = OpenStruct.new
  options.typename = :ulong
  ulongpairset = false
  uint64keypairset = false
  opts = OptionParser.new
  opts.on("--ulongpair",
          "generate code sorting key/value pairs of ulongs") do |x|
    ulongpairset = true
    options.typename = :ulongkeyvaluepair
  end
  opts.on("--uint64keypair","generate code sorting pairs of uint64_t keys") do |x|
    uint64keypairset = true
    options.typename = :uint64keypair
  end
  rest = opts.parse(argv)
  if not rest.empty?
    usage(opts,"superfluous arguments")
  end
  if ulongpairset and uint64keypairset
    STDERR.puts "#{$0}: cannot combine option --ulongpair and --uint64keypair"
    exit 1
  end
  return options
end

def makekey(options)
  if options.typename == :ulong
    return "ulong"
  elsif options.typename == :ulongkeyvaluepair
    return "ulongpair"
  else
    return "uint64keypair"
  end
end

def maketype(options)
  if options.typename == :ulong
    return "GtUword"
  elsif options.typename == :ulongkeyvaluepair
    return "GtUwordPair"
  else
    return "Gtuint64keyPair"
  end
end

def derefptr(ptr,options,comp="a")
  if m = ptr.match(/^&(\w+)/)
    var = m[1]
  else
    var = ptr
  end
  if options.typename == :ulong
    if var == ptr
      return "*#{var}"
    else
      return "#{var}"
    end
  elsif options.typename == :ulongkeyvaluepair
    if var == ptr
      return "#{var}->#{comp}"
    else
      return "#{var}.#{comp}"
    end
  else
    if var == ptr
      return "#{var}->uint64_#{comp}"
    else
      return "#{var}.uint64_#{comp}"
    end
  end
end

def radixkey(var,options)
  if options.typename == :ulong or options.typename == :ulongkeyvaluepair
    return "GT_RADIX_KEY(UINT8_MAX,rightshift,#{derefptr(var,options)})"
  else
    return "(rightshift > (sizeof (GtUword) - 1) * CHAR_BIT) ?\n" +
           "GT_RADIX_KEY(UINT8_MAX,rightshift - sizeof (GtUword) * CHAR_BIT,\n" +
           "#{derefptr(var,options)}) :\n" +
           "GT_RADIX_KEY(UINT8_MAX,rightshift,#{derefptr(var,options,"b")})"
  end
end

def compare_smaller(ptr1,ptr2,options)
  if options.typename == :ulong or options.typename == :ulongkeyvaluepair
    return "#{derefptr(ptr1,options)} < #{derefptr(ptr2,options)}"
  else
    return "gt_radixsort_compare_smaller(#{ptr1},#{ptr2})"
  end
end

options = parseargs(ARGV)

print <<END_OF_FILE
static #{maketype(options)} gt_radixsort_#{makekey(options)}_bin_get(
                                            const GtRadixbuffer *rbuf,
                                            GtUword binnum)
{
  return rbuf->values.#{makekey(options)}ptr[(binnum << rbuf->log_bufsize) +
                             (GtUword) rbuf->nextidx[binnum]];
}

static void gt_radixsort_#{makekey(options)}_bin_update(
                                    #{maketype(options)} *target,
                                    GtRadixbuffer *rbuf,
                                    GtUword binnum,
                                    #{maketype(options)} value)
{
  GtUword binoffset = binnum << rbuf->log_bufsize;

  rbuf->values.#{makekey(options)}ptr[binoffset +
                                      (GtUword) rbuf->nextidx[binnum]] = value;
  if ((GtUword) rbuf->nextidx[binnum] < rbuf->buf_size - 1)
  {
    rbuf->nextidx[binnum]++;
  } else
  {
    GtUword j;
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
  GtUword binoffset, binnum, bufoffset,
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
    count[#{radixkey("sp",options)}]++;
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
      rbuf->values.#{makekey(options)}ptr[bufoffset + j] =
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
    #{maketype(options)} currentvalue =
      gt_radixsort_#{makekey(options)}_bin_get(rbuf,binnum);
    while (true)
    {
      binnum = #{radixkey("&currentvalue",options)};
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
    gt_radixsort_#{makekey(options)}_bin_update(source,rbuf,binnum,
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
      #{maketype(options)} *targetptr, *valptr;

      valptr =
        rbuf->values.#{makekey(options)}ptr + (binnum << rbuf->log_bufsize);
      targetptr = source + rbuf->startofbin[binnum+1] - bufleft;
      for (j=0; j < bufleft; j++)
      {
        targetptr[j] = valptr[j];
      }
    }
  }
}

static void gt_radixsort_#{makekey(options)}_uncached_shuffle(
                       GtRadixbuffer *rbuf,
                       #{maketype(options)} *source,
                       GtCountbasetype len,
                       size_t rightshift)
{
  GtUword binnum, nextbin;
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
    count[#{radixkey("sp",options)}]++;
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
               (#{radixkey("&currentvalue",options)});
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
  gt_assert(rbuf != NULL);
  if ((GtUword) len > rbuf->cachesize)
  {
    gt_radixsort_#{makekey(options)}_cached_shuffle(rbuf,source,len,rightshift);
  } else
  {
    gt_radixsort_#{makekey(options)}_uncached_shuffle(rbuf,source,len,
                                                      rightshift);
  }
}

static void
gt_radixsort_#{makekey(options)}_inplace_insertionsort(#{maketype(options)} *a,
                                               GtCountbasetype a_size)
{
  #{maketype(options)} *optr, *iptr, *end = a + a_size;

  for (optr = a + 1; optr < end; optr++)
  {
    if (#{compare_smaller("optr","(optr-1)",options)})
    {
      #{maketype(options)} currentElement = *optr;

      *optr = *(optr-1);
      for (iptr = optr-1;
           iptr > a && #{compare_smaller("&currentElement","(iptr-1)",options)};
           iptr--)
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
  GtUword binnum;

  for (binnum = 0; binnum <= UINT8_MAX; binnum++)
  {
    GtCountbasetype width = rbuf->endofbin[binnum] - rbuf->startofbin[binnum];

    if (width >= (GtCountbasetype) 2)
    {
      #{maketype(options)} *ptr = source + rbuf->startofbin[binnum];

      if (width == (GtCountbasetype) 2)
      {
        if (#{compare_smaller("(ptr+1)","ptr",options)})
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
    gt_radixsort_#{makekey(options)}_shuffle(rbuf,
                         currentstackelem.left.#{makekey(options)}ptr,
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
