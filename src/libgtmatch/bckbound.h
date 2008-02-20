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

#ifndef BCKBOUND_H
#define BCKBOUND_H

#include <assert.h>
#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "seqpos-def.h"
#include "intcode-def.h"

typedef struct
{
  Seqpos left;
  unsigned long nonspecialsinbucket,
                specialsinbucket;
} Bucketspecification;

typedef struct
{
  Seqpos offset,
         left,
         right;
} Lcpinterval;

typedef struct
{
  Seqpos totallength,
         *leftborder,
         *countspecialcodes;
  Codetype numofallcodes,
           numofspecialcodes,
           **multimappower,
           *basepower,
           *filltable;
  unsigned long **distpfxidx;
} Bcktab;

int mapbcktab(Bcktab *bcktab,
              const Str *indexname,
              Seqpos totallength,
              unsigned int numofchars,
              unsigned int prefixlength,
              Error *err);

void freebcktab(Bcktab *bcktab,bool mapped);

void initbcktabwithNULL(Bcktab *bcktab);

int allocBcktab(Bcktab *bcktab,
                Seqpos totallength,
                unsigned int numofchars,
                unsigned int prefixlength,
                unsigned int codebits,
                unsigned int maxcodevalue,
                Error *err);

void updatebckspecials(Bcktab *bcktab,
                       Codetype code,
                       unsigned int numofchars,
                       unsigned int prefixindex,
                       unsigned int prefixlength);

Codetype codedownscale(const Bcktab *bcktab,
                       Codetype code,
                       unsigned int prefixindex,
                       unsigned int maxprefixlen);

void addfinalbckspecials(Bcktab *bcktab,unsigned int numofchars,
                         Seqpos specialcharacters);

int bcktab2file(FILE *fp,
                const Bcktab *bcktab,
                Error *err);

unsigned int calcbucketboundsparts(Bucketspecification *bucketspec,
                                   const Bcktab *bcktab,
                                   Codetype code,
                                   Codetype maxcode,
                                   Seqpos totalwidth,
                                   unsigned int rightchar,
                                   unsigned int numofchars);

void calcbucketboundaries(Bucketspecification *bucketspec,
                          const Bcktab *bcktab,
                          Codetype code);

#endif
