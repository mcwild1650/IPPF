CFLAGS=-Wall -Wfatal-errors

all: ippf

ippf:
	$(CC) $(CFLAGS) config.c driver.c field.c restart.c tools.c -o ippf
