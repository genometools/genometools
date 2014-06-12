#!/bin/sh -ex

# create directory for distributions
DISTDIR=distall
mkdir $DISTDIR
TMPFILE=`mktemp` || exit 1


# Linux 32-bit
make cleanup
make amalgamation=yes cairo=no m32=yes MACHINE=i386 \
     CPPFLAGS='-fno-stack-protector -U_FORTIFY_SOURCE -D_GNU_SOURCE' $*
make amalgamation=yes cairo=no m32=yes MACHINE=i386 \
     CPPFLAGS='-fno-stack-protector -U_FORTIFY_SOURCE -D_GNU_SOURCE' manuals
make amalgamation=yes cairo=no m32=yes MACHINE=i386 \
     CPPFLAGS='-fno-stack-protector -U_FORTIFY_SOURCE -D_GNU_SOURCE' dist > $TMPFILE
DISTRIBUTION=`tail -n 1 $TMPFILE`
cp -f $DISTRIBUTION $DISTDIR


# Linux 64-bit (barbone)
make clean
make 64bit=yes wrapmemcpy=yes \
     CPPFLAGS='-fno-stack-protector -U_FORTIFY_SOURCE -D_GNU_SOURCE' \
     amalgamation=yes cairo=no $*
make 64bit=yes wrapmemcpy=yes \
     CPPFLAGS='-fno-stack-protector -U_FORTIFY_SOURCE -D_GNU_SOURCE' \
     amalgamation=yes cairo=no DISTSUFFIX=-barebone dist > $TMPFILE
DISTRIBUTION=`tail -n 1 $TMPFILE`
cp -f $DISTRIBUTION $DISTDIR


# Linux 64-bit (complete)
make clean
make 64bit=yes wrapmemcpy=yes \
     CPPFLAGS='-fno-stack-protector -U_FORTIFY_SOURCE -D_GNU_SOURCE' \
     amalgamation=yes $*
make 64bit=yes wrapmemcpy=yes \
     CPPFLAGS='-fno-stack-protector -U_FORTIFY_SOURCE -D_GNU_SOURCE' \
     amalgamation=yes DISTSUFFIX=-complete dist > $TMPFILE
DISTRIBUTION=`tail -n 1 $TMPFILE`
cp -f $DISTRIBUTION $DISTDIR


# ARMv7_hf
# on Ubuntu: sudo apt-get install gcc-arm-linux-gnueabihf
make clean
make MACHINE=ARMv7                   \
     CC=arm-linux-gnueabihf-gcc      \
     AR=arm-linux-gnueabihf-ar       \
     STRIP=arm-linux-gnueabihf-strip \
     amalgamation=yes                \
     cairo=no                        \
     $*
make MACHINE=ARMv7                   \
     CC=arm-linux-gnueabihf-gcc      \
     AR=arm-linux-gnueabihf-ar       \
     STRIP=arm-linux-gnueabihf-strip \
     amalgamation=yes                \
     cairo=no                        \
     dist > $TMPFILE
DISTRIBUTION=`tail -n 1 $TMPFILE`
cp -f $DISTRIBUTION $DISTDIR


# ARMv6_hf
# on Linux:
# cd && mkdir -p arm && cd arm
# git clone --depth 1 https://github.com/raspberrypi/tools.git
make clean
make MACHINE=ARMv6              \
     CC="ccache $HOME/arm/tools/arm-bcm2708/gcc-linaro-arm-linux-gnueabihf-raspbian/bin/arm-linux-gnueabihf-gcc" \
     AR=$HOME/arm/tools/arm-bcm2708/gcc-linaro-arm-linux-gnueabihf-raspbian/bin/arm-linux-gnueabihf-ar           \
     STRIP=$HOME/arm/tools/arm-bcm2708/gcc-linaro-arm-linux-gnueabihf-raspbian/bin/arm-linux-gnueabihf-strip     \
     amalgamation=yes           \
     cairo=no                   \
     $*
make MACHINE=ARMv6              \
     CC="ccache $HOME/arm/tools/arm-bcm2708/gcc-linaro-arm-linux-gnueabihf-raspbian/bin/arm-linux-gnueabihf-gcc" \
     AR=$HOME/arm/tools/arm-bcm2708/gcc-linaro-arm-linux-gnueabihf-raspbian/bin/arm-linux-gnueabihf-ar           \
     STRIP=$HOME/arm/tools/arm-bcm2708/gcc-linaro-arm-linux-gnueabihf-raspbian/bin/arm-linux-gnueabihf-strip     \
     amalgamation=yes           \
     cairo=no                   \
     dist > $TMPFILE
DISTRIBUTION=`tail -n 1 $TMPFILE`
cp -f $DISTRIBUTION $DISTDIR


# Darwin 32-bit
# https://launchpad.net/~flosoft/+archive/cross-apple
make clean
make SYSTEM=Darwin                     \
     MACHINE=i386                      \
     CC=i686-apple-darwin10-gcc        \
     AR=i686-apple-darwin10-ar         \
     RANLIB=i686-apple-darwin10-ranlib \
     STRIP=i686-apple-darwin10-strip   \
     amalgamation=yes                  \
     cairo=no                          \
     $*
make SYSTEM=Darwin                     \
     MACHINE=i386                      \
     CC=i686-apple-darwin10-gcc        \
     AR=i686-apple-darwin10-ar         \
     RANLIB=i686-apple-darwin10-ranlib \
     STRIP=i686-apple-darwin10-strip   \
     amalgamation=yes                  \
     cairo=no                          \
     dist > $TMPFILE
DISTRIBUTION=`tail -n 1 $TMPFILE`
cp -f $DISTRIBUTION $DISTDIR


# Darwin 64-bit
# https://launchpad.net/~flosoft/+archive/cross-apple
make clean
make SYSTEM=Darwin                     \
     MACHINE=i386                      \
     64bit=yes                         \
     CC=i686-apple-darwin10-gcc        \
     AR=i686-apple-darwin10-ar         \
     RANLIB=i686-apple-darwin10-ranlib \
     STRIP=i686-apple-darwin10-strip   \
     amalgamation=yes                  \
     cairo=no                          \
     $*
make SYSTEM=Darwin                     \
     MACHINE=i386                      \
     64bit=yes                         \
     CC=i686-apple-darwin10-gcc        \
     AR=i686-apple-darwin10-ar         \
     RANLIB=i686-apple-darwin10-ranlib \
     STRIP=i686-apple-darwin10-strip   \
     amalgamation=yes                  \
     cairo=no                          \
     dist > $TMPFILE
DISTRIBUTION=`tail -n 1 $TMPFILE`
cp -f $DISTRIBUTION $DISTDIR
