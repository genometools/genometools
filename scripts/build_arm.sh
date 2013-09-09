#!/bin/sh -ex

# build manuals for distributions first
make cleanup
make cairo=no manuals

# ARMv7_hf
# on Ubuntu: sudo apt-get install gcc-arm-linux-gnueabihf
make clean
make MACHINE=ARMv7_hf                \
     CC=arm-linux-gnueabihf-gcc      \
     AR=arm-linux-gnueabihf-ar       \
     STRIP=arm-linux-gnueabihf-strip \
     cairo=no                        \
     dist $*

# ARMv7_el
# on Ubuntu: sudo apt-get install gcc-arm-linux-gnueabi
make clean
make MACHINE=ARMv7_el              \
     CC=arm-linux-gnueabi-gcc      \
     AR=arm-linux-gnueabi-ar       \
     STRIP=arm-linux-gnueabi-strip \
     cairo=no                      \
     dist $*

# ARMv6_hf
# on Linux:
# cd && mkdir -p arm && cd arm
# git clone --depth 1 https://github.com/raspberrypi/tools.git
make clean
make MACHINE=ARMv6_hf            \
     CC="ccache $HOME/arm/tools/arm-bcm2708/gcc-linaro-arm-linux-gnueabihf-raspbian/bin/arm-linux-gnueabihf-gcc" \
     AR=$HOME/arm/tools/arm-bcm2708/gcc-linaro-arm-linux-gnueabihf-raspbian/bin/arm-linux-gnueabihf-ar           \
     STRIP=$HOME/arm/tools/arm-bcm2708/gcc-linaro-arm-linux-gnueabihf-raspbian/bin/arm-linux-gnueabihf-strip     \
     cairo=no                    \
     dist $*
