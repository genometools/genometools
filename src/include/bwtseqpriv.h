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
#ifndef BWTSEQPRIV_H_INCLUDED
#define BWTSEQPRIV_H_INCLUDED
#include "bwtseq.h"
#include "encidxseq.h"

struct BWTSeq
{
  enum seqBaseEncoding type;
  struct encIdxSeq *seqIdx;
  size_t alphabetSize;
  EISHint hint;
  unsigned locateSampleInterval; /*< no sampling if 0*/
  Seqpos *count;
};

#endif /* BWTSEQPRIV_H_INCLUDED */
