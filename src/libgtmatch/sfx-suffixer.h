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

#ifndef SFX_SUFFIXER_H
#define SFX_SUFFIXER_H
#include "libgtcore/env.h"
#include "encseq-def.h"
#include "readmode-def.h"
#include "measure-time-if.h"
#include "seqpos-def.h"

typedef struct Sfxiterator Sfxiterator;

void freeSfxiterator(Sfxiterator **sfxiterator,Env *env);

Sfxiterator *newsfxiterator(Seqpos specialcharacters,
                            Seqpos specialranges,
                            const Encodedsequence *encseq,
                            Readmode readmode,
                            unsigned int numofchars,
                            unsigned int prefixlength,
                            unsigned int numofparts,
                            Measuretime *mtime,
                            Env *env);

const Seqpos *nextSfxiterator(Seqpos *len,Measuretime *mtime,
                              Sfxiterator *sfi,Env *env);

#endif
