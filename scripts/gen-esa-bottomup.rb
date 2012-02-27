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
  options.withlastsuftabvalue = true
  options.lcptypeulong = false
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
  opts.on("--nolastvalue","no proceccing of lastsuftabvalue") do |x|
    options.withlastsuftabvalue = false
  end
  opts.on("--lcptypeulong","use unsigned long for lcp type") do |x|
    options.lcptypeulong = true
  end
  rest = opts.parse(argv)
  if not rest.empty?
    usage(opts,"superfluous arguments")
  end
  return options
end

def processleafedgeargs(options)
  if options.absolute
    return "unsigned long, /* position */"
  else
    return "unsigned long, /* seqnum */\n" +
           "unsigned long, /* relpos */"
  end
end

def lcptype(options)
  if options.lcptypeulong
    return "const unsigned long"
  else
    return "const uint16_t"
  end
end

def previoussuffix_param_get(options)
  if options.absolute
    return "previoussuffix,"
  else
    return "previousseqnum,\n" +
           "                          previousrelpos,"
  end
end

def previoussuffix_expr_get(options)
  if options.absolute
    return "lastsuftabvalue,"
  else 
    return "gt_seqnumrelpos_decode_seqnum(snrp,lastsuftabvalue),\n" +
           "                        gt_seqnumrelpos_decode_relpos(snrp,lastsuftabvalue),"
  end
end

def previousseqnumrelpos(options)
  if not options.absolute
    return "previousseqnum = gt_seqnumrelpos_decode_seqnum(snrp,previoussuffix);\n" +
           "    previousrelpos = gt_seqnumrelpos_decode_relpos(snrp,previoussuffix);"
  else
    return "/* Nothing */"
  end
end

def processbranching_call1(key,options)
  if options.with_process_branching
    return "if (processbranchingedge_#{key}(firstedge,\n" +
           "                   TOP_ESA_BOTTOMUP_#{key}.lcp,\n" +
           "                   &TOP_ESA_BOTTOMUP_#{key}.info,\n" +
           "                   lastinterval->lcp,\n" +
           "                   lastinterval->rb - lastinterval->lb + 1,\n" +
           "                   &lastinterval->info,\n" +
           "                   bustate,\n" +
           "                   err) != 0)\n" +
           "        {\n" +
           "          haserr = true;\n" +
           "          break;\n" +
           "        }"
  else
    return "/* Nothing */"
  end
end

def processbranching_call2(key,options)
  if options.with_process_branching
    return "unsigned long lastintervallcp = lastinterval->lcp,\n" +
           "              lastintervalrb = lastinterval->rb;\n" +
           "        PUSH_ESA_BOTTOMUP_#{key}(lcpvalue,lastintervallb);\n" +
           "        if (processbranchingedge_#{key}(true,\n" +
           "                       TOP_ESA_BOTTOMUP_#{key}.lcp,\n" +
           "                       &TOP_ESA_BOTTOMUP_#{key}.info,\n" +
           "                       lastintervallcp,\n" +
           "                       lastintervalrb - lastintervallb + 1,\n" +
           "                       NULL,\n" +
           "                       bustate,\n" +
           "                       err) != 0)\n" +
           "        {\n" +
           "          haserr = true;\n" +
           "          break;\n" +
           "        }"
  else
    return "PUSH_ESA_BOTTOMUP_#{key}(lcpvalue,lastintervallb);"
  end
end

def processlcpinterval_decl(key,options)
  if options.with_process_lcpinterval
    return "static int processlcpinterval_#{key}(unsigned long,\n" +
                       "    GtBUinfo_#{key} *,\n" +
                       "    GtBUstate_#{key} *,\n" +
                       "    GtError *err);"
  else
    return "/* Nothing */"
  end
end

def processlcpinterval_call1(key,options)
  if options.with_process_lcpinterval
    return "if (processlcpinterval_#{key}(lastinterval->lcp,\n" +
           "                             &lastinterval->info,\n" +
           "                             bustate,\n" +
           "                             err) != 0)\n" +
           "      {\n" +
           "        haserr = true;\n" +
           "        break;\n" +
           "      }"
  else
    return "/* Nothing */"
  end
end

def processlcpinterval_call2(key,options)
  if options.with_process_lcpinterval
    return "if (processlcpinterval_#{key}(TOP_ESA_BOTTOMUP_#{key}.lcp,\n" +
           "                             &TOP_ESA_BOTTOMUP_#{key}.info,\n" +
           "                             bustate,\n" +
           "                             err) != 0)\n" +
           "      {\n" +
           "        haserr = true;\n" +
           "      }"
  else
    return "/* Nothing */"
  end
end

def mainloop(key,options)
print <<END_OF_FILE
    #{previousseqnumrelpos(options)}
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
                          #{previoussuffix_param_get(options)}
                          bustate,
                          err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_assert(lastinterval == NULL);
    while (lcpvalue < TOP_ESA_BOTTOMUP_#{key}.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_#{key};
      lastinterval->rb = idx;
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
    if (haserr)
    {
      break;
    }
    if (lcpvalue > TOP_ESA_BOTTOMUP_#{key}.lcp)
    {
      if (lastinterval != NULL)
      {
        unsigned long lastintervallb = lastinterval->lb;
        #{processbranching_call2(key,options)}
        lastinterval = NULL;
      } else
      {
        PUSH_ESA_BOTTOMUP_#{key}(lcpvalue,idx);
        if (processleafedge_#{key}(true,
                            TOP_ESA_BOTTOMUP_#{key}.lcp,
                            &TOP_ESA_BOTTOMUP_#{key}.info,
                            #{previoussuffix_param_get(options)}
                            bustate,
                            err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
END_OF_FILE
end

def lastsuftabvalue_fromarray(options)
  if not options.usefile
    return "lastsuftabvalue = bucketofsuffixes[numberofsuffixes-1];"
  else
    return "/* Nothing */"
  end
end

def lastsuftabvalue_get(key,options)
print <<END_OF_FILE
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP_#{key}.lcp > 0)
  {
    #{lastsuftabvalue_fromarray(options)}
    if (processleafedge_#{key}(false,
                        TOP_ESA_BOTTOMUP_#{key}.lcp,
                        &TOP_ESA_BOTTOMUP_#{key}.info,
                        #{previoussuffix_expr_get(options)}
                        bustate,
                        err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_#{key}.rb = idx;
      #{processlcpinterval_call2(key,options)}
    }
  }
END_OF_FILE
end

def seqnumrelpos_include(options)
  if not options.absolute
    return "#include \"seqnumrelpos.h\""
  else
    return "/* Nothing */"
  end
end

def processbranchingedge_decl(key,options)
  if options.with_process_branching
    return "static int processbranchingedge_#{key}(bool firstsucc,\n" +
           "    unsigned long,\n" +
           "    GtBUinfo_#{key} *,\n" +
           "    unsigned long,\n" +
           "    unsigned long,\n" +
           "    GtBUinfo_#{key} *,\n" +
           "    GtBUstate_#{key} *,\n" +
           "    GtError *);"
  else
   return "/* Nothing*/"
  end
end

def return_snrp_decl(options)
  if not options.absolute
    return "const GtSeqnumrelpos *snrp,"
  else
    return "/* Nothing*/"
  end
end

def return_previous_decl(options)
  if not options.absolute
    return "unsigned long previousseqnum = 0,\n" +
           "              previousrelpos = 0;"
  else
    return "/* Nothing*/"
  end
end

options = parseargs(ARGV)
key=options.key

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

#include <limits.h>
#include "core/ma.h"
#include "esa-seqread.h"
#{seqnumrelpos_include(options)}

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

if options.usefile
print <<END_OF_FILE

int gt_esa_bottomup_#{key}(Sequentialsuffixarrayreader *ssar,
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
  for (idx = 0; idx < numberofsuffixes; idx++)
  {
    NEXTSEQUENTIALLCPTABVALUEWITHLAST(lcpvalue,lastsuftabvalue,ssar);
    NEXTSEQUENTIALSUFTABVALUE(previoussuffix,ssar);
END_OF_FILE
else
print <<END_OF_FILE
int gt_esa_bottomup_RAM_#{key}(const unsigned long *bucketofsuffixes,
                        #{lcptype(options)} *lcptab_bucket,
                        unsigned long numberofsuffixes,
                        GtArrayGtBUItvinfo_#{key} *stack,
                        GtBUstate_#{key} *bustate,
                        #{return_snrp_decl(options)}
                        GtError *err)
{
  #{return_previous_decl(options)}
  const unsigned long incrementstacksize = 32UL;
  unsigned long lcpvalue,
                previoussuffix,
                lastsuftabvalue,
                idx;
  GtBUItvinfo_#{key} *lastinterval = NULL;
  bool haserr = false,
       firstedge,
       firstedgefromroot = true; /* must be part of state */

  gt_assert(numberofsuffixes > 0);
  PUSH_ESA_BOTTOMUP_#{key}(0,0);
  for (idx = 0; idx < numberofsuffixes-1; idx++)
  {
    lcpvalue = (unsigned long) lcptab_bucket[idx+1];
    previoussuffix = bucketofsuffixes[idx];
END_OF_FILE
end
mainloop(key,options)
puts "  }"
if options.withlastsuftabvalue
  lastsuftabvalue_get(key,options)
end
if options.usefile
  puts "  gt_GtArrayGtBUItvinfo_delete_#{key}(stack,bustate);"
else
  puts "  stack->nextfreeGtBUItvinfo = 0; /* empty the stack */"
end
puts "  return haserr ? -1 : 0;"
puts "}"
