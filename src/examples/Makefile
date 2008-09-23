# set prefix for libgenometools
prefix ?= /usr/local

CC:=gcc
CFLAGS:=-Wall -Werror -Wunused-parameter -g -Os -pipe \
        -I$(prefix)/include/genometools \
        -I/usr/include/cairo \
        -I/usr/local/include/cairo \
        -I/opt/local/include/cairo
LD:=$(CC)
LDFLAGS:=-L$(prefix)/lib

all: gff3validator

gff3validator: gff3validator.o
	$(LD) $(LDFLAGS) $< -lm -lgenometools -o $@

# generic compilation rule which creates dependency file on the fly
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

# read dependencies
-include $(wildcard *.d)

.PHONY: clean
clean:
	rm -f *.[od]
