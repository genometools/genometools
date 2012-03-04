#!/usr/bin/env ruby

require 'ostruct'
require 'optparse'

def usage(opts,msg)
  STDERR.puts "#{$0}: #{msg}\n#{opts.to_s}"
  exit 1
end

def parseargs(argv)
  options = OpenStruct.new
  options.key = nil
  options.usefile = false
  options.absolute = false
  options.with_process_branching = true
  options.with_process_lcpinterval = true
  options.process_lastvalue = true
  options.withlastfrompreviousbucket = false
  options.gtlcpvaluetypeset = false
  options.nodeclarations = false
  options.additionaluint32bucket = false
  opts = OptionParser.new
  opts.on("-k","--key STRING","use given key as suffix for all symbols") do |x|
    options.key = x
  end
  opts.on("--reader","generate code for the esa-reader") do |x|
    options.usefile = true
  end
  opts.on("--absolute","use absolute start positions to address suffixes") do |x|
    options.absolute = true
  end
  opts.on("--no_process_branchingedge","no processbranchingedge function") do |x|
    options.with_process_branching = false
  end
  opts.on("--no_process_lcpinterval","no processlcpinterval function") do |x|
    options.with_process_lcpinterval = false
  end
  opts.on("--no_process_lastvalue","no proceccing of lastsuftabvalue") do |x|
    options.process_lastvalue = false
  end
  opts.on("--gtlcpvaluetypeset","use GtLcpvaluetype for lcp type") do |x|
    options.gtlcpvaluetypeset = true
  end
  opts.on("--withlastfrompreviousbucket","process last value from previous bucket") do |x|
    options.withlastfrompreviousbucket = true
  end
  opts.on("--additionaluint32bucket","add uint32-bucket argument") do |x|
    options.additionaluint32bucket = true
  end
  opts.on("--no_declarations","do not output declarations") do |x|
    options.nodeclarations = true
  end
  rest = opts.parse(argv)
  if not rest.empty?
    usage(opts,"superfluous arguments")
  end
  if options.key.nil?
    usage(opts,"option --key is mandatory")
  end
  return options
end

def processleafedgeargs(options)
  if options.absolute
    return "unsigned long, /* position */"
  else
    return "unsigned long, /* seqnum */
    unsigned long, /* relpos */"
  end
end

def spacer(width)
  s = ""
  0.upto(width-1) do |idx|
    s += " "
  end
  return s
end

def previoussuffix_expr_get(width,variable,options)
  if options.absolute
    return variable + ","
  else
    return "gt_seqnumrelpos_decode_seqnum(snrp,#{variable}),
           #{" " * width} gt_seqnumrelpos_decode_relpos(snrp,#{variable}),"
  end
end

def previoussuffix_param_get(width,options)
  return previoussuffix_expr_get(width,"previoussuffix",options)
end

def processbranching_call1(key,options)
  if options.with_process_branching
    return "if (processbranchingedge_#{key}(firstedge,
                   TOP_ESA_BOTTOMUP_#{key}.lcp,
                   &TOP_ESA_BOTTOMUP_#{key}.info,
                   lastinterval->lcp,
                   lastinterval->rb - lastinterval->lb + 1,
                   &lastinterval->info,
                   bustate,
                   err) != 0)
        {
          haserr = true;
        }"
  else
    return "/* no call to processbranchingedge_#{key} */"
  end
end

def processbranching_call2(key,options)
  if options.with_process_branching
    return "unsigned long lastintervallcp = lastinterval->lcp,
              lastintervalrb = lastinterval->rb;
        PUSH_ESA_BOTTOMUP_#{key}(lcpvalue,lastintervallb);
        if (processbranchingedge_#{key}(true,
                       TOP_ESA_BOTTOMUP_#{key}.lcp,
                       &TOP_ESA_BOTTOMUP_#{key}.info,
                       lastintervallcp,
                       lastintervalrb - lastintervallb + 1,
                       NULL,
                       bustate,
                       err) != 0)
        {
          haserr = true;
        }"
  else
    return "PUSH_ESA_BOTTOMUP_#{key}(lcpvalue,lastintervallb);"
  end
end

def processlcpinterval_decl(key,options)
  if options.with_process_lcpinterval
    return "static int processlcpinterval_#{key}(unsigned long,
    GtBUinfo_#{key} *,
    GtBUstate_#{key} *,
    GtError *err);"
  else
    return "/* no declaration of processlcpinterval_#{key} */"
  end
end

def processlcpinterval_call1(key,options)
  if options.with_process_lcpinterval
    return "if (processlcpinterval_#{key}(lastinterval->lcp,
                             &lastinterval->info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
      }"
  else
    return "/* no call to processlcpinterval_#{key} */"
  end
end

def processlcpinterval_call2(key,options)
  if options.with_process_lcpinterval
    return "if (processlcpinterval_#{key}(TOP_ESA_BOTTOMUP_#{key}.lcp,
                             &TOP_ESA_BOTTOMUP_#{key}.info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
      }"
  else
    return "/* no call to processlcpinterval_#{key} */"
  end
end

def showidxexpr(options)
  if options.withlastfrompreviousbucket
    return "idx + bustate->idxoffset"
  else
    return "idx"
  end
end

def process_suf_lcp(key,options)
print <<END_OF_FILE
    gt_assert(stack->nextfreeGtBUItvinfo > 0);
    if (lcpvalue <= TOP_ESA_BOTTOMUP_#{key}.lcp)
    {
      if (TOP_ESA_BOTTOMUP_#{key}.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge_#{key}(firstedge,
                          TOP_ESA_BOTTOMUP_#{key}.lcp,
                          &TOP_ESA_BOTTOMUP_#{key}.info,
                          #{previoussuffix_param_get(14,options)}
                          bustate,
                          err) != 0)
      {
        haserr = true;
      }
    }
    gt_assert(lastinterval == NULL);
    while (!haserr && lcpvalue < TOP_ESA_BOTTOMUP_#{key}.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_#{key};
      lastinterval->rb = #{showidxexpr(options)};
      #{processlcpinterval_call1(key,options)}
      if (lcpvalue <= TOP_ESA_BOTTOMUP_#{key}.lcp)
      {
        if (TOP_ESA_BOTTOMUP_#{key}.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
        #{processbranching_call1(key,options)}
        lastinterval = NULL;
      }
    }
    if (!haserr && lcpvalue > TOP_ESA_BOTTOMUP_#{key}.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallb = lastinterval->lb;
        #{processbranching_call2(key,options)}
        lastinterval = NULL;
      } else
      {
        PUSH_ESA_BOTTOMUP_#{key}(lcpvalue,#{showidxexpr(options)});
        if (processleafedge_#{key}(true,
                            TOP_ESA_BOTTOMUP_#{key}.lcp,
                            &TOP_ESA_BOTTOMUP_#{key}.info,
                            #{previoussuffix_param_get(16,options)}
                            bustate,
                            err) != 0)
        {
          haserr = true;
        }
      }
    }
END_OF_FILE
end

def lastsuftabvalue_fromarray(options)
  if not options.usefile
    return "unsigned long lastsuftabvalue = bucketofsuffixes[numberofsuffixes-1];"
  else
    return "/* no assignment to lastsuftabvalue */"
  end
end

def lastsuftabvalue_get(key,options)
  if options.process_lastvalue
print <<END_OF_FILE
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP_#{key}.lcp > 0)
  {
    #{lastsuftabvalue_fromarray(options)}
    if (processleafedge_#{key}(false,
                        TOP_ESA_BOTTOMUP_#{key}.lcp,
                        &TOP_ESA_BOTTOMUP_#{key}.info,
                        #{previoussuffix_expr_get(12,"lastsuftabvalue",options)}
                        bustate,
                        err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_#{key}.rb = #{showidxexpr(options)};
      #{processlcpinterval_call2(key,options)}
    }
  }
END_OF_FILE
  else
print <<END_OF_FILE
  if (!haserr)
  {
    bustate->previousbucketlastsuffix
      = #{accessbucketofsuffixes("numberofsuffixes-1",options)}
    bustate->firstedgefromroot = firstedgefromroot;
  }
END_OF_FILE
  end
end

def seqnumrelpos_include(options)
  if not options.absolute
    return "#include \"seqnumrelpos.h\""
  else
    return "/* no include for seqnumrelpos.h */"
  end
end

def processbranchingedge_decl(key,options)
  if options.with_process_branching
    return "static int processbranchingedge_#{key}(bool firstsucc,
    unsigned long,
    GtBUinfo_#{key} *,
    unsigned long,
    unsigned long,
    GtBUinfo_#{key} *,
    GtBUstate_#{key} *,
    GtError *);"
  else
   return "/* no declaration of processbranchingedge_#{key} */"
  end
end

def return_snrp_decl(options)
  if not options.absolute
    return "const GtSeqnumrelpos *snrp,"
  else
    return "/* no parameter snrp */"
  end
end

def lcptype(options)
  if options.gtlcpvaluetypeset
    return "const GtLcpvaluetype"
  else
    return "const uint16_t"
  end
end

def additionaluint32bucket(options)
  if options.additionaluint32bucket
    return "const uint32_t *bucketofsuffixes_uint32,\n" + (" " * 24) + lcptype(options)
  else
    return "#{lcptype(options)}"
  end
end

def accessbucketofsuffixes(idx,options)
  if options.additionaluint32bucket
    return "bucketofsuffixes != NULL ? bucketofsuffixes[#{idx}]
                                 : (unsigned long)
                                   bucketofsuffixes_uint32[#{idx}];"
  else
    return "bucketofsuffixes[#{idx}];"
  end
end

def formatargv(argv)
  s = "\n  #{argv[0]} #{argv[1]}"
  argv.each_with_index do |arg,i|
    if i > 1
      s += "\n  " + argv[i]
    end
  end
  return s
end

def initfirstinterval(key,options)
  if options.withlastfrompreviousbucket
print <<END_OF_FILE

  if (bustate->previousbucketlastsuffix == ULONG_MAX)
  {
    PUSH_ESA_BOTTOMUP_#{key}(0,0);
    firstedgefromroot = true;
  } else
  {
    firstedgefromroot = bustate->firstedgefromroot;
  }
END_OF_FILE
  else
print <<END_OF_FILE

  PUSH_ESA_BOTTOMUP_#{key}(0,0);
  firstedgefromroot = true;
END_OF_FILE
  end
end

options = parseargs(ARGV)
key = options.key

print <<END_OF_FILE
/*
  Copyright (c) 2011-2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

/*
  THIS FILE IS GENERATED by\n  #{$0}#{formatargv(ARGV)}.
  DO NOT EDIT.
*/

#include <limits.h>
#include "core/ma.h"
#include "esa-seqread.h"
#{seqnumrelpos_include(options)}
END_OF_FILE

if not options.nodeclarations
print <<END_OF_FILE

static void initBUinfo_#{key}(GtBUinfo_#{key} *,GtBUstate_#{key} *);

static void freeBUinfo_#{key}(GtBUinfo_#{key} *,GtBUstate_#{key} *);

static int processleafedge_#{key}(bool,
    unsigned long,
    GtBUinfo_#{key} *,
    #{processleafedgeargs(options)}
    GtBUstate_#{key} *,
    GtError *err);

#{processbranchingedge_decl(key,options)}

#{processlcpinterval_decl(key,options)}

#define TOP_ESA_BOTTOMUP_#{key}\\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo-1]

#define POP_ESA_BOTTOMUP_#{key}\\
        (stack->spaceGtBUItvinfo + (--stack->nextfreeGtBUItvinfo))

#define PUSH_ESA_BOTTOMUP_#{key}(LCP,LB)\\
        if (stack->nextfreeGtBUItvinfo >= stack->allocatedGtBUItvinfo)\\
        {\\
          gt_assert(stack->nextfreeGtBUItvinfo ==\\
                    stack->allocatedGtBUItvinfo);\\
          stack->spaceGtBUItvinfo\\
            = allocateBUstack_#{key}(stack->spaceGtBUItvinfo,\\
                              stack->allocatedGtBUItvinfo,\\
                              stack->allocatedGtBUItvinfo+incrementstacksize,\\
                              bustate);\\
          stack->allocatedGtBUItvinfo += incrementstacksize;\\
        }\\
        gt_assert(stack->spaceGtBUItvinfo != NULL);\\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo].lcp = LCP;\\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo].lb = LB;\\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo++].rb = ULONG_MAX

typedef struct
{
  unsigned long lcp, lb, rb;
  GtBUinfo_#{key} info;
} GtBUItvinfo_#{key};

typedef struct
{
  GtBUItvinfo_#{key} *spaceGtBUItvinfo;
  unsigned long allocatedGtBUItvinfo,
                nextfreeGtBUItvinfo;
} GtArrayGtBUItvinfo_#{key};

GtArrayGtBUItvinfo_#{key} *gt_GtArrayGtBUItvinfo_new_#{key}(void)
{
  GtArrayGtBUItvinfo_#{key} *stack = gt_malloc(sizeof (*stack));

  stack->spaceGtBUItvinfo = NULL;
  stack->allocatedGtBUItvinfo = stack->nextfreeGtBUItvinfo = 0;
  return stack;
}

void gt_GtArrayGtBUItvinfo_delete_#{key}(GtArrayGtBUItvinfo_#{key} *stack,
                                  GtBUstate_#{key} *state)
{
  unsigned long idx;

  for (idx=0; idx<stack->allocatedGtBUItvinfo; idx++)
  {
    freeBUinfo_#{key}(&stack->spaceGtBUItvinfo[idx].info,state);
  }
  gt_free(stack->spaceGtBUItvinfo);
  gt_free(stack);
}

static GtBUItvinfo_#{key} *allocateBUstack_#{key}(GtBUItvinfo_#{key} *ptr,
                                   unsigned long currentallocated,
                                   unsigned long allocated,
                                   GtBUstate_#{key} *state)
{
  unsigned long idx;
  GtBUItvinfo_#{key} *itvinfo;

  itvinfo = gt_realloc(ptr,sizeof (*itvinfo) * allocated);
  gt_assert(allocated > currentallocated);
  for (idx=currentallocated; idx<allocated; idx++)
  {
    initBUinfo_#{key}(&itvinfo[idx].info,state);
  }
  gt_assert(itvinfo != NULL);
  return itvinfo;
}
END_OF_FILE
end

if options.withlastfrompreviousbucket
print <<END_OF_FILE

static int gt_esa_bottomup_RAM_previousfromlast_#{key}(
                        unsigned long previoussuffix,
                        unsigned long lcpvalue,
                        GtArrayGtBUItvinfo_#{key} *stack,
                        GtBUstate_#{key} *bustate,
                        #{return_snrp_decl(options)}
                        GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long idx = 0;
  GtBUItvinfo_#{key} *lastinterval = NULL;
  bool haserr = false, firstedge,
       firstedgefromroot = bustate->firstedgefromroot;

END_OF_FILE

  process_suf_lcp(key,options)

print <<END_OF_FILE
  if (!haserr)
  {
    bustate->firstedgefromroot = firstedgefromroot;
  }
  return haserr ? -1 : 0;
}
END_OF_FILE
end

if options.usefile
print <<END_OF_FILE

static int gt_esa_bottomup_#{key}(Sequentialsuffixarrayreader *ssar,
                    GtBUstate_#{key} *bustate,
                    #{return_snrp_decl(options)}
                    GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix = 0,
                idx,
                numberofsuffixes,
                lastsuftabvalue = 0;
  GtBUItvinfo_#{key} *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot = true;
  GtArrayGtBUItvinfo_#{key} *stack;

  stack = gt_GtArrayGtBUItvinfo_new_#{key}();
  PUSH_ESA_BOTTOMUP_#{key}(0,0);
  numberofsuffixes = gt_Sequentialsuffixarrayreader_nonspecials(ssar);
  for (idx = 0; !haserr && idx < numberofsuffixes; idx++)
  {
    NEXTSEQUENTIALLCPTABVALUEWITHLAST(lcpvalue,lastsuftabvalue,ssar);
    NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
END_OF_FILE
else
print <<END_OF_FILE

static int gt_esa_bottomup_RAM_#{key}(const unsigned long *bucketofsuffixes,
                        #{additionaluint32bucket(options)} *lcptab_bucket,
                        unsigned long numberofsuffixes,
                        GtArrayGtBUItvinfo_#{key} *stack,
                        GtBUstate_#{key} *bustate,
                        #{return_snrp_decl(options)}
                        GtError *err)
{
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix,
                idx;
  GtBUItvinfo_#{key} *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot;
END_OF_FILE

initfirstinterval(key,options)

print <<END_OF_FILE
  gt_assert (numberofsuffixes > 0);
  for (idx = 0; !haserr && idx < numberofsuffixes-1; idx++)
  {
    lcpvalue = (unsigned long) lcptab_bucket[idx+1];
    previoussuffix = #{accessbucketofsuffixes("idx",options)}
END_OF_FILE
end
process_suf_lcp(key,options)
puts "  }"
lastsuftabvalue_get(key,options)
if options.usefile
  puts "  gt_GtArrayGtBUItvinfo_delete_#{key}(stack,bustate);"
elsif options.process_lastvalue
  puts "  stack->nextfreeGtBUItvinfo = 0; /* empty the stack */"
end
puts "  return haserr ? -1 : 0;"
puts "}"
