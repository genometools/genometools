#!/bin/sh

cd testsuite
env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb -keywords 'gt_greedyfwdmat'
cd ..
