#!/bin/sh

if [ -z "$MYEMAIL" ]; then
    echo "Need to set environment variable MYEMAIL"
    exit 1
fi

if [ -z "$MYNAME" ]; then
    echo "Need to set environment variable MYNAME"
    exit 1
fi

YEAR=`date +"%Y"`
if [ -z "$MYAFFILIATION" ]; then
  AFF=""
else
  AFF="Copyright (c) ${YEAR} ${MYAFFILIATION}"
fi

TMPFILE=`mktemp TMP.XXXXXX` || exit 1
for filename in $*
do
  cat > ${TMPFILE} <<EOF
/*
  Copyright (c) ${YEAR} ${MYNAME} <${MYEMAIL}>
  ${AFF}

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
EOF
  cat ${filename} >> ${TMPFILE}
  mv ${TMPFILE} ${filename}
done
