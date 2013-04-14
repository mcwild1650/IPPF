CFLAGS=-Wall -Wfatal-errors
LDFLAGS=-lm
CC=mpicc

all: ippf serial

ippf:
	$(CC) $(CFLAGS) config.c calc.c driver.c field.c restart.c \
	space.c tools.c $(LDFLAGS) -o ippf

serial:
	$(CC) $(CFLAGS) ParabolaFlow_v3_2.c $(LDFLAGS) -o serial

clean:
	rm -rf ippf serial
