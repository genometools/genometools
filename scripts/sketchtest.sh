#!/bin/sh
set -e -x

cd testsuite
env -i ./testsuite.rb -keywords 'gt_sketch or gt_ruby or gt_python or gt_scripts'
cd ..
