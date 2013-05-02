CFLAGS=-Wall -Wfatal-errors -DCLOCK_RATE=1.6e9
LDFLAGS=-lm
HEADERS=config.h calc.h driver.h field.h restart.h space.h tools.h grid.h
SRC=config.c calc.c driver.c field.c restart.c space.c tools.c grid.c
CC=mpicc

all: ippf serial

debug:CFLAGS+=-g -DDEBUG
debug: all

ippf: $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o ippf

serial: ParabolaFlow_v3_2.c
	$(CC) $(CFLAGS) ParabolaFlow_v3_2.c $(LDFLAGS) -o serial

clean:
	rm -rf ippf serial

test: ippf serial
	rm -f good.dat out.dat
	mpirun -np 2 ./ippf run testconf.txt
	./serial < testin.txt
	diff good.dat out.dat
