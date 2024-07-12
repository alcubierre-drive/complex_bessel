CC ?= gcc
CFLAGS +=
LDFLAGS +=
use: use.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ -I.. -L.. -Wl,-rpath=.. -lbessel
