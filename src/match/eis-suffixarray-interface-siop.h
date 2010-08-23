/*
  Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#ifndef EIS_SUFFIXARRAY_INTERFACE_SIOP_H
#define EIS_SUFFIXARRAY_INTERFACE_SIOP_H

#include "eis-suffixarray-interface.h"
#include "eis-suffixarray-interface-priv.h"

static inline MRAEnc *
SAINewMRAEnc(const SuffixarrayFileInterface *sai)
{
  return gt_SANewMRAEnc(gt_encseq_alphabet(sai->sa->encseq));
}

static inline const GtEncseq *
SAIGetEncSeq(const SuffixarrayFileInterface *sai)
{
  return sai->sa->encseq;
}

static inline GtReadmode
SAIGetGtReadmode(const SuffixarrayFileInterface *sai)
{
  return sai->sa->readmode;
}

static inline unsigned long
SAIGetLength(const SuffixarrayFileInterface *sai)
{
  return SASSGetLength(constSAI2SASS(sai));
}

#endif
