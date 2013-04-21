CC=mpicc
CFLAGS=-Wall -Wfatal-errors
LDFLAGS=-lm
HEADERS=config.h calc.h driver.h field.h restart.h space.h tools.h
SRC=config.c calc.c driver.c field.c restart.c space.c tools.c

all: ippf serial

debug:CFLAGS+=-g -DDEBUG
debug: all

ippf: $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o ippf

serial:
	$(CC) $(CFLAGS) ParabolaFlow_v3_2.c $(LDFLAGS) -o serial

clean:
	rm -rf ippf serial
