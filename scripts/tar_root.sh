#!/bin/bash
#
# Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2012 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

tar --owner=root -cf $1.tar $1 > /dev/null 2>&1

if [ $? -ne 0 ]
then
  /usr/bin/gnutar --owner=root -cf $1.tar $1 > /dev/null 2>&1
else
  exit 0
fi

if [ $? -ne 0 ]
then
  tar -cf $1.tar $1
else
  exit 0
fi
