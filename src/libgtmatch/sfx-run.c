/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include "spacedef.h"
#include "alphadef.h"
#include "sfx-optdef.h"
#include "encseq-def.h"
#include "measure-time-if.h"
#include "esafileend.h"
#include "intcode-def.h"
#include "filelength-def.h"
#include "chardef.h"

#include "alphabet.pr"
#include "measure-time.pr"
#include "opensfxfile.pr"
#include "sfx-cmpsuf.pr"
#include "sfx-makeprj.pr"
#include "sfx-opt.pr"
#include "sfx-outprj.pr"
#include "sfx-suffixer.pr"
#include "sfx-apfxlen.pr"

typedef struct
{
  FILE *outfpsuftab,
       *outfplcptab,
       *outfpllvtab,
       *outfpbwttab;
  Seqpos lastsuftabentryofpreviouspart, 
         absolutepos,
         numoflargelcpvalues,
         maxbranchdepth;
  const Encodedsequence *encseq;
  DefinedSeqpos longest;
} Outfileinfo;

static void initoutfileinfo(Outfileinfo *outfileinfo)
{
  outfileinfo->outfpsuftab = NULL;
  outfileinfo->outfplcptab = NULL;
  outfileinfo->outfpllvtab = NULL;
  outfileinfo->outfpbwttab = NULL;
  outfileinfo->absolutepos = 0;
  outfileinfo->lastsuftabentryofpreviouspart = 0;
  outfileinfo->numoflargelcpvalues = 0;
  outfileinfo->maxbranchdepth = 0;
  outfileinfo->longest.defined = false;
  outfileinfo->longest.valueseqpos = 0;
}

static int suftab2file(void *info,
                       const Seqpos *suftab,
                       Readmode readmode,
                       Seqpos widthofpart,
                       Env *env)
{
  Outfileinfo *outfileinfo = (Outfileinfo *) info;
  bool haserr = false;
  Seqpos pos;

  env_error_check(env);
  if (outfileinfo->outfpsuftab != NULL)
  {
    if (fwrite(suftab,
              sizeof (*suftab),
              (size_t) widthofpart,
              outfileinfo->outfpsuftab)
              != (size_t) widthofpart)
    {
      env_error_set(env,"cannot write " FormatSeqpos " items of size %u: "
                    "errormsg=\"%s\"",
           PRINTSeqposcast(widthofpart),
           (unsigned int) sizeof (*suftab),
           strerror(errno));
      haserr = true;
    }
  }
  if(!outfileinfo->longest.defined)
  {
    for(pos=0; pos<widthofpart; pos++)
    {
      if(suftab[pos] == 0)
      {
        outfileinfo->longest.defined = true;
        outfileinfo->longest.valueseqpos = outfileinfo->absolutepos + pos;
        break;
      }
    }
  }
  if (!haserr && outfileinfo->outfplcptab != NULL)
  {
    Seqpos lcpvalue;
    Uchar outvalue;
    Largelcpvalue largelcpvalue;
    int cmp;

    outvalue = (Uchar) 0; 
    if (outfileinfo->absolutepos == 0 &&
        fwrite(&outvalue,sizeof(Uchar),(size_t) 1,
               outfileinfo->outfplcptab) != (size_t) 1)
    {
      env_error_set(env,"cannot write 1 item of size %lu: errormsg=\"%s\"",
                         (unsigned long) sizeof(Uchar),
                         strerror(errno));
      haserr = true;
    }
    if(!haserr)
    {
      for(pos=0; pos<widthofpart; pos++)
      {
        if(pos > 0 || outfileinfo->absolutepos > 0)
        {
          cmp = comparetwosuffixes(outfileinfo->encseq,
                                   readmode,
                                   &lcpvalue,
                                   false,
                                   false,
                                   0,
                                   pos > 0 ? suftab[pos-1] 
                                           : outfileinfo->
                                             lastsuftabentryofpreviouspart,
                                   suftab[pos]);
          assert(cmp <= 0);
          if(outfileinfo->maxbranchdepth < lcpvalue)
          {
            outfileinfo->maxbranchdepth = lcpvalue;
          }
          if(lcpvalue >= (Seqpos) UCHAR_MAX)
          {
            outfileinfo->numoflargelcpvalues++;
            largelcpvalue.position = outfileinfo->absolutepos + pos;
            largelcpvalue.value = lcpvalue;
            if (fwrite(&largelcpvalue,sizeof(Largelcpvalue),(size_t) 1,
                       outfileinfo->outfpllvtab) != (size_t) 1)
            {
              env_error_set(env,"cannot write 1 item of size %lu: "
                                "errormsg=\"%s\"",
                                (unsigned long) sizeof(Largelcpvalue),
                                strerror(errno));
              haserr = true;
              break;
            }
            outvalue = (Uchar) UCHAR_MAX; 
          } else
          {
            outvalue = (Uchar) lcpvalue;
          }
          if (!haserr && fwrite(&outvalue,sizeof(Uchar),(size_t) 1,
                                outfileinfo->outfplcptab) != (size_t) 1)
          {
            env_error_set(env,"cannot write 1 item of size %lu: "
                              "errormsg=\"%s\"",
                              (unsigned long) sizeof(Uchar),
                              strerror(errno));
            haserr = true;
            break;
          }
        }
      }
    }
    if(!haserr)
    {
      outfileinfo->lastsuftabentryofpreviouspart = suftab[widthofpart-1];
    }
  }
  if(!haserr)
  {
    outfileinfo->absolutepos += widthofpart;
  }
  if (!haserr && outfileinfo->outfpbwttab != NULL)
  {
    Seqpos startpos;
    Uchar cc;

    for(pos=0; pos < widthofpart; pos++)
    {
      startpos = suftab[pos];
      if(startpos == 0)
      {
        cc = (Uchar) UNDEFBWTCHAR;
      } else
      {
        cc = getencodedchar(outfileinfo->encseq,startpos - 1,readmode);
      }
      if (fwrite(&cc,sizeof(Uchar),(size_t) 1,outfileinfo->outfpbwttab) 
                  != (size_t) 1)
      {
        env_error_set(env,"cannot write 1 item of size %lu: "\
                          "errormsg=\"%s\"",\
                        (unsigned long) sizeof(Uchar),\
                        strerror(errno));\
        haserr = true;
        break;
      }
    }
  }
  return haserr ? -1 : 0;
}

static int outal1file(const Str *indexname,const Alphabet *alpha,Env *env)
{
  FILE *al1fp;
  bool haserr = false;

  env_error_check(env);
  al1fp = opensfxfile(indexname,ALPHABETFILESUFFIX,"wb",env);
  if(al1fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    outputalphabet(al1fp,alpha);
    env_fa_xfclose(al1fp,env);
  }
  return haserr ? -1 : 0;
}

#define INITOUTFILEPTR(PTR,FLAG,SUFFIX)\
        if (!haserr && (FLAG))\
        {\
          outfileinfo.encseq = encseq;\
          PTR = opensfxfile(so->str_indexname,SUFFIX,"wb",env);\
          if((PTR) == NULL)\
          {\
            haserr = true;\
          }\
        }

#define SIZEOFBCKENTRY (2 * sizeof(Seqpos))

static int runsuffixerator(Suffixeratoroptions *so,Env *env)
{
  unsigned char numofchars = 0;
  unsigned long numofsequences;
  Seqpos totallength;
  Alphabet *alpha;
  Specialcharinfo specialcharinfo;
  Filelengthvalues *filelengthtab = NULL;
  Outfileinfo outfileinfo;
  bool haserr = false;
  Encodedsequence *encseq = NULL;
  Measuretime *mtime;

  env_error_check(env);
  inittheclock(&mtime,
               "determining sequence length and number of special symbols",
               env);
  alpha = assigninputalphabet(so->isdna,
                              so->isprotein,
                              so->str_smap,
                              so->filenametab,
                              env);
  if (alpha == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (scanfastasequence(&numofsequences,
                          &totallength,
                          &specialcharinfo,
                          so->filenametab,
                          &filelengthtab,
                          getsymbolmapAlphabet(alpha), 
                          so->isplain,
                          env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    numofchars = (unsigned char) getnumofcharsAlphabet(alpha);
    if (outal1file(so->str_indexname,alpha,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    deliverthetime(stdout,mtime,"computing sequence encoding",env);
    encseq = initencodedseq(true,
                            so->filenametab,
                            so->isplain,
                            NULL,
                            totallength,
                            &specialcharinfo,
                            alpha,
                            str_length(so->str_sat) > 0
                                  ? str_get(so->str_sat)
                                  : NULL,
                            env);
    if (encseq == NULL)
    {
      haserr = true;
    } else
    {
      if (so->outtistab)
      {
        if (flushencseqfile(so->str_indexname,encseq,env) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  initoutfileinfo(&outfileinfo);
  INITOUTFILEPTR(outfileinfo.outfpsuftab,so->outsuftab,SUFTABSUFFIX);
  INITOUTFILEPTR(outfileinfo.outfplcptab,so->outlcptab,LCPTABSUFFIX);
  INITOUTFILEPTR(outfileinfo.outfpllvtab,so->outlcptab,LARGELCPTABSUFFIX);
  INITOUTFILEPTR(outfileinfo.outfpbwttab,so->outbwttab,BWTTABSUFFIX);
  if (!haserr)
  {
    printf("# specialcharacters=" FormatSeqpos "\n",
           PRINTSeqposcast(specialcharinfo.specialcharacters));
    printf("# specialranges=" FormatSeqpos "\n",
           PRINTSeqposcast(specialcharinfo.specialranges));
    if(so->readmode == Complementmode || so->readmode == Reversecomplementmode)
    {
      if(!isdnaalphabet(alpha))
      {
        env_error_set(env,"option %s only can be used for DNA alphabets",
                          so->readmode == Complementmode ? "-cpl" : "rcl");
        haserr = true;
      }
    }
    if (!haserr && (so->outsuftab || so->outbwttab || so->outlcptab))
    {
      if(so->prefixlength == PREFIXLENGTH_AUTOMATIC)
      {
        so->prefixlength = recommendedprefixlength(numofchars,
                                                   totallength,
                                                   SIZEOFBCKENTRY);
        printf("# automatically determined prefixlength = %u\n",
                so->prefixlength);
      } else
      {
        unsigned int maxprefixlen;

        maxprefixlen 
          = whatisthemaximalprefixlength(numofchars,
                                         totallength,
                                         SIZEOFBCKENTRY,
                                         (unsigned int) PREFIXLENBITS);
        if(checkprefixlength(maxprefixlen,so->prefixlength,env) != 0)
        {
          haserr = true;
        } else
        {
          showmaximalprefixlength(maxprefixlen,
                                  recommendedprefixlength(
                                  numofchars,
                                  totallength,
                                  SIZEOFBCKENTRY));
        }
      }
      if(!haserr)
      {
        if (suffixerator(suftab2file,
                         &outfileinfo,
                         specialcharinfo.specialcharacters,
                         specialcharinfo.specialranges,
                         encseq,
                         so->readmode,
                         (uint32_t) numofchars,
                         (uint32_t) so->prefixlength,
                         (uint32_t) so->numofparts,
                         mtime,
                         env) != 0)
        {
          haserr = true;
        }
      }
    } else
    {
      if(so->readmode != Forwardmode)
      {
        env_error_set(env,"readmode = %u only makes sense in combination with "
                          "at least one of the options -suf, -lcp, or -bwt",
                          so->readmode);
        haserr = true;
      }
    }
  }
  env_fa_xfclose(outfileinfo.outfpsuftab,env);
  env_fa_xfclose(outfileinfo.outfplcptab,env);
  env_fa_xfclose(outfileinfo.outfpllvtab,env);
  env_fa_xfclose(outfileinfo.outfpbwttab,env);
  if (!haserr)
  {
    if (outprjfile(so->str_indexname,
                   so->filenametab,
                   so->readmode,
                   filelengthtab,
                   totallength,
                   numofsequences,
                   &specialcharinfo,
                   (uint32_t) so->prefixlength,
                   outfileinfo.numoflargelcpvalues,
                   outfileinfo.maxbranchdepth,
                   &outfileinfo.longest,
                   env) != 0)
    {
      haserr = true;
    }
  }
  FREESPACE(filelengthtab);
  freeAlphabet(&alpha,env);
  freeEncodedsequence(&encseq,env);
  deliverthetime(stdout,mtime,NULL,env);
  return haserr ? -1 : 0;
}

int parseargsandcallsuffixerator(int argc,const char **argv,Env *env)
{
  Suffixeratoroptions so;
  int retval;
  bool haserr = false;

  retval = suffixeratoroptions(&so,argc,argv,env);
  if(retval == 0)
  {
    printf("# sizeof(Seqpos)=%lu\n",
            (unsigned long) (sizeof(Seqpos) * CHAR_BIT));
  }
  if (retval == 0)
  {
    if (runsuffixerator(&so,env) < 0)
    {
      haserr = true;
    }
  } else
  {
    if (retval < 0)
    {
      haserr = true;
    }
  }
  wrapsfxoptions(&so,env);
  return haserr ? -1 : 0;
}
