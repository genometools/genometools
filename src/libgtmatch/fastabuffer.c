/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include "libgtcore/env.h"
#include "libgtcore/strarray.h"
#include "format64.h"
#include "fastabuffer.h"
#include "spacedef.h"
#include "chardef.h"

#define FASTASEPARATOR    '>'
#define NEWLINESYMBOL     '\n'

FastaBuffer* fastabuffer_new(const StrArray *filenametab,
                             const Uchar *symbolmap,
                             bool plainformat,
                             Filelengthvalues **filelengthtab,
                             Sequencedescription *sequencedescription,
                             unsigned long *characterdistribution,
                             Env *env)
{
  FastaBuffer *fb;
  env_error_check(env);
  fb = env_ma_calloc(env, 1, sizeof (FastaBuffer));
  fb->plainformat = plainformat;
  fb->filenum = 0;
  fb->firstoverallseq = true;
  fb->firstseqinfile = true;
  fb->nextfile = true;
  fb->nextread = fb->nextfree = 0;
  fb->filenametab = filenametab;
  fb->symbolmap = symbolmap;
  fb->complete = false;
  fb->lastspeciallength = 0;
  fb->sequencedescription = sequencedescription;
  if (filelengthtab != NULL)
  {
    ALLOCASSIGNSPACE(*filelengthtab,NULL,Filelengthvalues,
                     strarray_size(filenametab));
    fb->filelengthtab = *filelengthtab;
  } else
  {
    fb->filelengthtab = NULL;
  }
  fb->characterdistribution = characterdistribution;
  return fb;
}

static int advancefastabufferstate(FastaBuffer *fb,Env *env)
{
  int currentchar;
  unsigned long currentposition = 0, currentfileadd = 0, currentfileread = 0;
  Uchar charcode;
  char *savebuffer;

  env_error_check(env);
  while (true)
  {
    if (currentposition >= (unsigned long) FILEBUFFERSIZE)
    {
      if (fb->filelengthtab != NULL)
      {
        fb->filelengthtab[fb->filenum].length
          += (uint64_t) currentfileread;
        fb->filelengthtab[fb->filenum].effectivelength
          += (uint64_t) currentfileadd;
      }
      break;
    }
    if (fb->nextfile)
    {
      if (fb->filelengthtab != NULL)
      {
        fb->filelengthtab[fb->filenum].length = 0;
        fb->filelengthtab[fb->filenum].effectivelength = 0;
      }
      fb->nextfile = false;
      fb->indesc = false;
      fb->firstseqinfile = true;
      currentfileadd = 0;
      currentfileread = 0;
      fb->linenum = (uint64_t) 1;
      fb->inputstream = genfile_xopen(strarray_get(fb->filenametab,
                                                  (unsigned long) fb->filenum),
                                       "rb", env);
    } else
    {
      /* XXX: use genfile_xread() */
      currentchar = genfile_getc(fb->inputstream);
      if (currentchar == EOF)
      {
        genfile_xclose(fb->inputstream, env);
        fb->inputstream = NULL;
        if (fb->filelengthtab != NULL)
        {
          fb->filelengthtab[fb->filenum].length += currentfileread;
          fb->filelengthtab[fb->filenum].effectivelength += currentfileadd;
        }
        if ((unsigned long) fb->filenum == strarray_size(fb->filenametab) - 1)
        {
          fb->complete = true;
          break;
        }
        fb->filenum++;
        fb->nextfile = true;
      } else
      {
        currentfileread++;
        if (fb->indesc)
        {
          if (currentchar == NEWLINESYMBOL)
          {
            fb->linenum++;
            fb->indesc = false;
          }
          if (fb->sequencedescription != NULL)
          {
            if (currentchar == NEWLINESYMBOL)
            {
              STOREINARRAY(&fb->sequencedescription->headerbuffer,char,128,
                           '\0');
              ALLOCASSIGNSPACE(savebuffer,NULL,char,
                               fb->sequencedescription->headerbuffer.
                                                         nextfreechar);
              strcpy(savebuffer,
                     fb->sequencedescription->headerbuffer.spacechar);
              queue_add(fb->sequencedescription->descptr,savebuffer,env);
              fb->sequencedescription->headerbuffer.nextfreechar = 0;
            } else
            {
              STOREINARRAY(&fb->sequencedescription->headerbuffer,char,128,
                           currentchar);
            }
          }
        } else
        {
          if (!isspace((int) currentchar))
          {
            if (currentchar == FASTASEPARATOR)
            {
              if (fb->firstoverallseq)
              {
                fb->firstoverallseq = false;
                fb->firstseqinfile = false;
              } else
              {
                if (fb->firstseqinfile)
                {
                  fb->firstseqinfile = false;
                } else
                {
                  currentfileadd++;
                }
                fb->bufspace[currentposition++] = (Uchar) SEPARATOR;
                fb->lastspeciallength++;
              }
              fb->indesc = true;
            } else
            {
              if (fb->symbolmap == NULL)
              {
                fb->bufspace[currentposition++] = (Uchar) currentchar;
              } else
              {
                charcode = fb->symbolmap[(unsigned int) currentchar];
                if (charcode == (Uchar) UNDEFCHAR)
                {
                  env_error_set(env,"illegal character '%c':"
                                    " file \"%s\", line " Formatuint64_t,
                                currentchar,
                                strarray_get(fb->filenametab,
                                             (unsigned long) fb->filenum),
                                PRINTuint64_tcast(fb->linenum));
                  return -1;
                }
                if (ISSPECIAL(charcode))
                {
                  fb->lastspeciallength++;
                } else
                {
                  if (fb->lastspeciallength > 0)
                  {
                    fb->lastspeciallength = 0;
                  }
                  if (fb->characterdistribution != NULL)
                  {
                    fb->characterdistribution[charcode]++;
                  }
                }
                fb->bufspace[currentposition++] = charcode;
              }
              currentfileadd++;
            }
          }
        }
      }
    }
  }
  if (fb->firstoverallseq)
  {
    env_error_set(env,"no sequences in multiple fasta file(s) %s ...",
                  strarray_get(fb->filenametab,0));
    return -2;
  }
  fb->nextfree = currentposition;
  return 0;
}

static int advancePlainbufferstate(FastaBuffer *fb,Env *env)
{
  int currentchar;
  unsigned long currentposition = 0, currentfileread = 0;

  env_error_check(env);
  if (fb->sequencedescription != NULL)
  {
    env_error_set(env,"no headers in plain sequence file");
    return -1;
  }
  while (true)
  {
    if (currentposition >= (unsigned long) FILEBUFFERSIZE)
    {
      if (fb->filelengthtab != NULL)
      {
        fb->filelengthtab[fb->filenum].length
           += (uint64_t) currentfileread;
        fb->filelengthtab[fb->filenum].effectivelength
           += (uint64_t) currentfileread;
      }
      break;
    }
    if (fb->nextfile)
    {
      if (fb->filelengthtab != NULL)
      {
        fb->filelengthtab[fb->filenum].length = 0;
        fb->filelengthtab[fb->filenum].effectivelength = 0;
      }
      fb->nextfile = false;
      fb->firstseqinfile = true;
      currentfileread = 0;
      fb->inputstream = genfile_xopen(strarray_get(fb->filenametab,
                                                  (unsigned long) fb->filenum),
                                       "rb", env);
    } else
    {
      /* XXX: use genfile_xread() */
      currentchar = genfile_getc(fb->inputstream);
      if (currentchar == EOF)
      {
        genfile_xclose(fb->inputstream, env);
        fb->inputstream = NULL;
        if (fb->filelengthtab != NULL)
        {
          fb->filelengthtab[fb->filenum].length
            += (uint64_t) currentfileread;
          fb->filelengthtab[fb->filenum].effectivelength
            += (uint64_t) currentfileread;
        }
        if ((unsigned long) fb->filenum == strarray_size(fb->filenametab) - 1)
        {
          fb->complete = true;
          break;
        }
        fb->filenum++;
        fb->nextfile = true;
      } else
      {
        currentfileread++;
        fb->bufspace[currentposition++] = (Uchar) currentchar;
      }
    }
  }
  if (currentposition == 0)
  {
    env_error_set(env,"no characters in plain file(s) %s ...",
                  strarray_get(fb->filenametab,0));
    return -2;
  }
  fb->nextfree = currentposition;
  return 0;
}

int advanceformatbufferstate(FastaBuffer *fb,Env *env)
{
  env_error_check(env);
  if (fb->plainformat)
  {
    return advancePlainbufferstate(fb,env);
  }
  return advancefastabufferstate(fb,env);
}

void fastabuffer_delete(FastaBuffer *fb, Env *env)
{
  if (!fb) return;
  genfile_xclose(fb->inputstream, env);
  env_ma_free(fb, env);
}
