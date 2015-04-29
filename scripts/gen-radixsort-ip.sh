#!/bin/sh

gen-radixsort-ip.rb > src/core/radixsort-ip-ulong.inc
gen-radixsort-ip.rb --ulongpair > src/core/radixsort-ip-ulongpair.inc
gen-radixsort-ip.rb --seedpair > src/core/radixsort-ip-seedpair.inc
