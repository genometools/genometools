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

#ifndef SEQ_COL_REP_H
#define SEQ_COL_REP_H

#include <stdio.h>
#include "core/seq_col.h"

typedef struct GtSeqColClass GtSeqColClass;

typedef void          (*GtSeqColFreeFunc)(GtSeqCol*);
typedef int           (*GtSeqColGrepDescFunc)(GtSeqCol*, char **seq,
                                              unsigned long start,
                                              unsigned long end,
                                              GtStr *seqid, GtError*);
typedef int           (*GtSeqColGrepDescMD5Func)(GtSeqCol*, const char **md5,
                                                 GtStr *seqid, GtError*);
typedef int           (*GtSeqColGrepDescSeqlenFunc)(GtSeqCol*, unsigned long*,
                                                    GtStr *, GtError*);
typedef int           (*GtSeqColMD5ToSeqFunc)(GtSeqCol*, char **seq,
                                              unsigned long start,
                                              unsigned long end,
                                              GtStr *md5_seqid, GtError *err);
typedef int           (*GtSeqColMD5ToDescFunc)(GtSeqCol*, GtStr *desc,
                                               GtStr *md5_seqid, GtError *err);
typedef int           (*GtSeqColMD5ToSeqlenFunc)(GtSeqCol*, unsigned long*,
                                                 GtStr *md5_seqid,
                                                 GtError *err);
typedef unsigned long (*GtSeqColNumFilesFunc)(const GtSeqCol*);
typedef unsigned long (*GtSeqColNumSeqsFunc)(const GtSeqCol*,
                                             unsigned long filenum);
typedef const char*   (*GtSeqColGetMD5Func)(const GtSeqCol*,
                                            unsigned long filenum,
                                            unsigned long seqnum);
typedef       char*   (*GtSeqColGetSeqFunc)(const GtSeqCol*,
                                            unsigned long filenum,
                                            unsigned long seqnum,
                                            unsigned long start,
                                            unsigned long end);
typedef       char*   (*GtSeqColGetDescFunc)(const GtSeqCol*,
                                             unsigned long filenum,
                                             unsigned long seqnum);
typedef unsigned long (*GtSeqColGetSeqlenFunc)(const GtSeqCol*,
                                               unsigned long filenum,
                                               unsigned long seqnum);

struct GtSeqColClass {
  size_t size;
  GtSeqColFreeFunc free;
  GtSeqColGrepDescFunc grep_desc;
  GtSeqColGrepDescMD5Func grep_desc_md5;
  GtSeqColGrepDescSeqlenFunc grep_desc_seqlen;
  GtSeqColMD5ToSeqFunc md5_to_seq;
  GtSeqColMD5ToDescFunc md5_to_desc;
  GtSeqColMD5ToSeqlenFunc md5_to_seqlen;
  GtSeqColNumFilesFunc num_files;
  GtSeqColNumSeqsFunc num_seqs;
  GtSeqColGetMD5Func get_md5;
  GtSeqColGetSeqFunc get_seq;
  GtSeqColGetDescFunc get_desc;
  GtSeqColGetSeqlenFunc get_seqlen;
};

struct GtSeqCol {
  const GtSeqColClass *c_class;
};

const GtSeqColClass* gt_seq_col_class_new(size_t size,
                                          GtSeqColFreeFunc free,
                                          GtSeqColGrepDescFunc grep_desc,
                                          GtSeqColGrepDescMD5Func grep_desc_md5,
                                          GtSeqColGrepDescSeqlenFunc
                                                               grep_desc_seqlen,
                                          GtSeqColMD5ToSeqFunc md5_to_seq,
                                          GtSeqColMD5ToDescFunc md5_to_desc,
                                          GtSeqColMD5ToSeqlenFunc md5_to_seqlen,
                                          GtSeqColNumFilesFunc num_files,
                                          GtSeqColNumSeqsFunc num_seqs,
                                          GtSeqColGetMD5Func get_md5,
                                          GtSeqColGetSeqFunc get_seq,
                                          GtSeqColGetDescFunc get_desc,
                                          GtSeqColGetSeqlenFunc get_seqlen);
GtSeqCol*      gt_seq_col_create(const GtSeqColClass*);
void*          gt_seq_col_cast(const GtSeqColClass*, const GtSeqCol*);

#endif
