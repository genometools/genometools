/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef MD5_ENCODER_H
#define MD5_ENCODER_H

typedef struct GtMD5Encoder GtMD5Encoder;

GtMD5Encoder* gt_md5_encoder_new(void);
void          gt_md5_encoder_add_block(GtMD5Encoder *enc, const char *message,
                                       unsigned long len);
void          gt_md5_encoder_finish(GtMD5Encoder *enc, unsigned char *output,
                                    char *outstr);
void          gt_md5_encoder_reset(GtMD5Encoder *enc);
void          gt_md5_encoder_delete(GtMD5Encoder *enc);
#endif
