/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <ctype.h>
#include "libgtcore/env.h"
#include "libgtcore/strarray.h"
#include "format64.h"
#include "fbs-def.h"
#include "spacedef.h"
#include "chardef.h"
#include "gqueue-def.h"

#include "genericqueue.pr"

#define GZIPSUFFIX        ".gz"
#define GZIPSUFFIXLENGTH  (sizeof (GZIPSUFFIX)-1)
#define FASTASEPARATOR    '>'
#define NEWLINESYMBOL     '\n'

static unsigned char checkgzipsuffix(const char *filename)
{
  size_t filenamelength = strlen (filename);

  if (filenamelength < GZIPSUFFIXLENGTH ||
      strcmp (filename + filenamelength - GZIPSUFFIXLENGTH, GZIPSUFFIX) != 0)
  {
    return 0;
  }
  return (unsigned char) 1;
}

static void opengenericstream(Genericstream *inputstream,const char *inputfile)
{
  if (checkgzipsuffix(inputfile))
  {
    inputstream->stream.gzippedstream = gzopen(inputfile,"rb");
    if (inputstream->stream.gzippedstream == NULL)
    {
      fprintf(stderr,"cannot open file \"%s\"\n",inputfile);
      exit(EXIT_FAILURE);
    }
    inputstream->isgzippedstream = true;
  } else
  {
    inputstream->stream.fopenstream = fopen(inputfile,"rb");
    if (inputstream->stream.fopenstream == NULL)
    {
      fprintf(stderr,"cannot open file \"%s\"\n",inputfile);
      exit(EXIT_FAILURE);
    }
    inputstream->isgzippedstream = false;
  }
}

static void closegenericstream(Genericstream *inputstream,const char *inputfile)
{
  int retval;

  if (inputstream->isgzippedstream)
  {
    retval = gzclose(inputstream->stream.gzippedstream);
  } else
  {
    retval = fclose(inputstream->stream.fopenstream);
  }
  if (retval != 0)
  {
    fprintf(stderr,"cannot close file \"%s\"\n",inputfile);
    exit(EXIT_FAILURE);
  }
}

void initformatbufferstate(Fastabufferstate *fbs,
                           const StrArray *filenametab,
                           const Uchar *symbolmap,
                           bool plainformat,
                           Filelengthvalues **filelengthtab,
                           Sequencedescription *sequencedescription,
                           Env *env)
{
  env_error_check(env);
  fbs->plainformat = plainformat;
  fbs->filenum = 0;
  fbs->firstoverallseq = true;
  fbs->firstseqinfile = true;
  fbs->nextfile = true;
  fbs->nextread = fbs->nextfree = 0;
  fbs->filenametab = filenametab;
  fbs->symbolmap = symbolmap;
  fbs->complete = false;
  fbs->lastspeciallength = 0;
  fbs->sequencedescription = sequencedescription;
  if(filelengthtab != NULL)
  {
    ALLOCASSIGNSPACE(*filelengthtab,NULL,Filelengthvalues,
                     strarray_size(filenametab));
    fbs->filelengthtab = *filelengthtab;
  } else
  {
    fbs->filelengthtab = NULL;
  }
}

static int advanceFastabufferstate(Fastabufferstate *fbs,Env *env)
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
      if(fbs->filelengthtab != NULL)
      {
        fbs->filelengthtab[fbs->filenum].length
          += (uint64_t) currentfileread;
        fbs->filelengthtab[fbs->filenum].effectivelength
          += (uint64_t) currentfileadd;
      }
      break;
    }
    if (fbs->nextfile)
    {
      if(fbs->filelengthtab != NULL)
      {
        fbs->filelengthtab[fbs->filenum].length = 0;
        fbs->filelengthtab[fbs->filenum].effectivelength = 0;
      }
      fbs->nextfile = false;
      fbs->indesc = false;
      fbs->firstseqinfile = true;
      currentfileadd = 0;
      currentfileread = 0;
      fbs->linenum = (uint32_t) 1;
      opengenericstream(&fbs->inputstream,
                        strarray_get(fbs->filenametab,
                        (unsigned long) fbs->filenum));
    } else
    {
      if (fbs->inputstream.isgzippedstream)
      {
        currentchar = gzgetc(fbs->inputstream.stream.gzippedstream);
      } else
      {
        currentchar = fgetc(fbs->inputstream.stream.fopenstream);
      }
      if (currentchar == EOF)
      {
        closegenericstream(&fbs->inputstream,
                           strarray_get(fbs->filenametab,
                                        (unsigned long) fbs->filenum));
        if(fbs->filelengthtab != NULL)
        {
          fbs->filelengthtab[fbs->filenum].length += currentfileread;
          fbs->filelengthtab[fbs->filenum].effectivelength += currentfileadd;
        }
        if ((unsigned long) fbs->filenum == strarray_size(fbs->filenametab) - 1)
        {
          fbs->complete = true;
          break;
        }
        fbs->filenum++;
        fbs->nextfile = true;
      } else
      {
        currentfileread++;
        if (fbs->indesc)
        {
          if (currentchar == NEWLINESYMBOL)
          {
            fbs->linenum++;
            fbs->indesc = false;
          }
          if(fbs->sequencedescription != NULL)
          {
            if(currentchar == NEWLINESYMBOL)
            {
              STOREINARRAY(&fbs->sequencedescription->headerbuffer,char,128,
                           '\0');
              ALLOCASSIGNSPACE(savebuffer,NULL,char,
                               fbs->sequencedescription->headerbuffer.
                                                         nextfreechar);
              strcpy(savebuffer,
                     fbs->sequencedescription->headerbuffer.spacechar);
              enqueuegeneric(fbs->sequencedescription->descptr,savebuffer,env);
              fbs->sequencedescription->headerbuffer.nextfreechar = 0;
            } else
            {
              STOREINARRAY(&fbs->sequencedescription->headerbuffer,char,128,
                           currentchar);
            }
          }
        } else
        {
          if (!isspace((int) currentchar))
          {
            if (currentchar == FASTASEPARATOR)
            {
              if (fbs->firstoverallseq)
              {
                fbs->firstoverallseq = false;
                fbs->firstseqinfile = false;
              } else
              {
                if (fbs->firstseqinfile)
                {
                  fbs->firstseqinfile = false;
                } else
                {
                  currentfileadd++;
                }
                fbs->bufspace[currentposition++] = (Uchar) SEPARATOR;
                fbs->lastspeciallength++;
              }
              fbs->indesc = true;
            } else
            {
              if(fbs->symbolmap == NULL)
              {
                fbs->bufspace[currentposition++] = (Uchar) currentchar;
              } else
              {
                charcode = fbs->symbolmap[(uint32_t) currentchar];
                if (charcode == (Uchar) UNDEFCHAR)
                {
                  env_error_set(env,"illegal character '%c':"
                                    " file \"%s\", line " Formatuint64_t,
                                currentchar,
                                strarray_get(fbs->filenametab,
                                             (unsigned long) fbs->filenum),
                                PRINTuint64_tcast(fbs->linenum));
                  return -1;
                }
                if (ISSPECIAL(charcode))
                {
                  fbs->lastspeciallength++;
                } else
                {
                  if (fbs->lastspeciallength > 0)
                  {
                    fbs->lastspeciallength = 0;
                  }
                }
                fbs->bufspace[currentposition++] = charcode;
              }
              currentfileadd++;
            }
          }
        }
      }
    }
  }
  if (fbs->firstoverallseq)
  {
    env_error_set(env,"no sequences in multiple fasta file(s) %s ...",
                  strarray_get(fbs->filenametab,0));
    return -2;
  }
  fbs->nextfree = currentposition;
  return 0;
}

static int advancePlainbufferstate(Fastabufferstate *fbs,Env *env)
{
  int currentchar;
  unsigned long currentposition = 0, currentfileread = 0;

  env_error_check(env);
  if(fbs->sequencedescription != NULL)
  {
    env_error_set(env,"no headers in plain sequence file");
    return -1;
  }
  while (true)
  {
    if (currentposition >= (unsigned long) FILEBUFFERSIZE)
    {
      if(fbs->filelengthtab != NULL)
      {
        fbs->filelengthtab[fbs->filenum].length
           += (uint64_t) currentfileread;
        fbs->filelengthtab[fbs->filenum].effectivelength
           += (uint64_t) currentfileread;
      }
      break;
    }
    if (fbs->nextfile)
    {
      if(fbs->filelengthtab != NULL)
      {
        fbs->filelengthtab[fbs->filenum].length = 0;
        fbs->filelengthtab[fbs->filenum].effectivelength = 0;
      }
      fbs->nextfile = false;
      fbs->firstseqinfile = true;
      currentfileread = 0;
      opengenericstream(&fbs->inputstream,
                        strarray_get(fbs->filenametab,
                        (unsigned long) fbs->filenum));
    } else
    {
      if (fbs->inputstream.isgzippedstream)
      {
        currentchar = gzgetc(fbs->inputstream.stream.gzippedstream);
      } else
      {
        currentchar = fgetc(fbs->inputstream.stream.fopenstream);
      }
      if (currentchar == EOF)
      {
        closegenericstream(&fbs->inputstream,
                           strarray_get(fbs->filenametab,
                                        (unsigned long) fbs->filenum));
        if(fbs->filelengthtab != NULL)
        {
          fbs->filelengthtab[fbs->filenum].length
            += (uint64_t) currentfileread;
          fbs->filelengthtab[fbs->filenum].effectivelength
            += (uint64_t) currentfileread;
        }
        if ((unsigned long) fbs->filenum == strarray_size(fbs->filenametab) - 1)
        {
          fbs->complete = true;
          break;
        }
        fbs->filenum++;
        fbs->nextfile = true;
      } else
      {
        currentfileread++;
        fbs->bufspace[currentposition++] = (Uchar) currentchar;
      }
    }
  }
  if(currentposition == 0)
  {
    env_error_set(env,"no characters in plain file(s) %s ...",
                  strarray_get(fbs->filenametab,0));
    return -2;
  }
  fbs->nextfree = currentposition;
  return 0;
}

int advanceformatbufferstate(Fastabufferstate *fbs,Env *env)
{
  if(fbs->plainformat)
  {
    return advancePlainbufferstate(fbs,env);
  }
  return advanceFastabufferstate(fbs,env);
}
