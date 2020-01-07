/*
  Copyright (c) 2011 Gordon Gremme <gordon@gremme.org>

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

#ifndef BUFFER_STREAM_API_H
#define BUFFER_STREAM_API_H

#include "extended/node_stream_api.h"

/* The <GtBufferStream> is a <GtNodeStream> that buffers <GtGenomeNode>s. */
typedef struct GtBufferStream GtBufferStream;

const GtNodeStreamClass* gt_buffer_stream_class(void);

/* Create a new <GtBufferStream>, reading from <in_stream>. The stream is
   initially configured to buffer nodes read from the input until
   <gt_buffer_stream_dequeue()> switches it to dequeue mode. In this mode,
   calls to <gt_node_stream_next()> will deliver the stored nodes in FIFO
   order. */
GtNodeStream*            gt_buffer_stream_new(GtNodeStream *in_stream);
/* Switch stream <bs> to dequeue mode. */
void                     gt_buffer_stream_dequeue(GtBufferStream *bs);

#endif
