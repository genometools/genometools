# set prefix for libgenometools
prefix ?= /usr/local

CC:=gcc
GT_CFLAGS:=-Wall -Werror -Wunused-parameter -g -Os -pipe \
           -I$(prefix)/include/genometools \
           -I/usr/include/cairo \
           -I/usr/local/include/cairo \
           -I/sw/include/cairo \
           -I/opt/local/include/cairo
LD:=$(CC)
GT_LDFLAGS:=-L$(prefix)/lib \
            -L/usr/X11R6/lib \
            -L/usr/local/lib \
            -L/opt/local/lib

all: gff3validator

gff3validator: gff3validator.o
	$(LD) $(LDFLAGS) $(GT_LDFLAGS) $< -lm -lgenometools -lcairo -o $@

# generic compilation rule which creates dependency file on the fly
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

# read dependencies
-include $(wildcard *.d)

.PHONY: clean
clean:
	rm -f *.[od]
