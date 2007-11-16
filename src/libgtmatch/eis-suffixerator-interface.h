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
/**
 * \file eis-suffixerator-interface.h
 * \brief Methods to call suffixerator functions through one object,
 * but have the same data available to multiple listeners.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
#ifndef EIS_SUFFIXERATOR_INTERFACE_H
#define EIS_SUFFIXERATOR_INTERFACE_H

#include "libgtcore/env.h"
#include "libgtmatch/sfx-suffixer.h"
#include "libgtmatch/sfx-optdef.h"
#include "libgtmatch/eis-mrangealphabet.h"

typedef unsigned listenerID;

typedef struct sfxInterface sfxInterface;

enum sfxDataRequest {
  SFX_REQUEST_NONE = 0,
  SFX_REQUEST_ALL = 7,
  SFX_REQUEST_ANY = SFX_REQUEST_ALL,
  SFX_REQUEST_SUFTAB = 1<<0,
  SFX_REQUEST_LCPTAB = 1<<1,
  SFX_REQUEST_BWTTAB = 1<<2,
};

extern sfxInterface *
newSfxInterface(Suffixeratoroptions *so,
                Verboseinfo *verbosity,
                Env *env);

extern sfxInterface *
newSfxInterfaceWithReaders(Suffixeratoroptions *so,
                           size_t numReaders,
                           enum sfxDataRequest *requests,
                           listenerID *ids,
                           Verboseinfo *verbosity,
                           Env *env);

extern int
deleteSfxInterface(sfxInterface *iface, Env *env);

extern const Alphabet *
getSfxIAlphabet(const sfxInterface *si);

extern const struct seqStats *
getSfxISeqStats(const sfxInterface *si);

extern Seqpos
getSfxILength(const sfxInterface *si);

extern const Uchar *
readSfxIESQRange(sfxInterface *iface, Seqpos start, Seqpos len,
                 Uchar *dest);

/**
 * @return >0 on success, 0 on error
 */
extern int
SfxIRegisterReader(sfxInterface *iface, listenerID *id,
                   enum sfxDataRequest request, Env *env);

extern size_t
readSfxIBWTRange(sfxInterface *iface, listenerID id, size_t len,
                 Uchar *dest, Env *env);

extern size_t
readSfxIBWTRangeSym(sfxInterface *iface, listenerID id, size_t len,
                    Symbol *dest, Env *env);
/**
 * @return actual number of symbols read
 */
extern size_t
readSfxILCPRange(sfxInterface *iface, listenerID id, size_t len,
                 Seqpos *dest, Env *env);

/**
 * @return actual number of symbols read
 */
extern size_t
readSfxISufTabRange(sfxInterface *iface, listenerID id, size_t len,
                    Seqpos *dest, Env *env);

#endif
