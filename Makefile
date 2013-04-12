CFLAGS=-Wall -Wfatal-errors
LDFLAGS=-lm

all: ippf serial

ippf:
	$(CC) $(CFLAGS) config.c driver.c field.c restart.c tools.c $(LDFLAGS) -o ippf

serial:
	$(CC) $(CFLAGS) ParabolaFlow_v3_2.c $(LDFLAGS) -o serial

clean:
	rm -rf ippf serial
