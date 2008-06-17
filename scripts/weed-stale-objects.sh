#! /bin/sh
#set -x
find obj -type f -name '*.o' -print | \
  while read fname ; do
  fspec="${fname%.o}"
  fspec="${fspec#obj/}"
  if [ \! -f "${fspec}.c" -a \! -f "${fspec}.cpp" \
    -a \! -f "${fspec}.cc" -a \! -f "${fspec}.cxx" ]; then
    echo "removing $fname"
    rm "$fname"
  fi
done