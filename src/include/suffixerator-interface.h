/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/
/**
 * \file suffixerator-interface.h
 * \brief Methods to call suffixerator functions through one object,
 * but have the same data available to multiple listeners.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
#ifndef SUFFIXERATOR_INTERFACE_H_INCLUDED
#define SUFFIXERATOR_INTERFACE_H_INCLUDED

#include <libgtcore/env.h>

#include <libgtmatch/sfx-suffixer.h>
#include <libgtmatch/sfx-optdef.h>

#include "mrangealphabet.h"

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
newSfxInterface(Suffixeratoroptions *so, Env *env);

extern sfxInterface *
newSfxInterfaceWithReaders(Suffixeratoroptions *so,
                           size_t numReaders,
                           enum sfxDataRequest *requests,
                           listenerID *ids, Env *env);

extern int
deleteSfxInterface(sfxInterface *iface, Env *env);

extern const Alphabet *
getSfxIAlphabet(const sfxInterface *si);

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

#endif /* SUFFIXERATOR_INTERFACE_H_INCLUDED */
