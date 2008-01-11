CWARNS = -Wall -Wextra -pedantic \
        -Waggregate-return \
        -Wcast-align \
        -Wmissing-prototypes \
        -Wnested-externs \
        -Wpointer-arith \
        -Wshadow \
        -Wwrite-strings \
        -Wcast-qual

LUADIR = ..
COPT = -O2 -DNDEBUG
CFLAGS = $(CWARNS) -ansi -I$(LUADIR)/lua -shared -o lpeg.so
CC = gcc

lpeg.so: lpeg.c
	$(CC) $(COPT) $(CFLAGS) lpeg.c

deb:	lpeg.c	
	$(CC) -g $(CFLAGS) lpeg.c; touch deb; rm -f opt

opt:	lpeg.c	
	$(CC) $(COPT) $(CFLAGS) lpeg.c; touch opt; rm -f deb

test: test.lua re.lua lpeg.so
	test.lua
