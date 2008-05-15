#! /bin/sh
find src/libgtmatch -type f \( -name '*.[ch]' -o -name '*.pr' \) -print | \
  etags -o TAGS.libgtmatch -
find . -type f \( -name '*.[ch]' -o -name '*.pr' \) -print | \
  etags -
