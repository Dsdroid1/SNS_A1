CC = gcc
CFLAGS = "-lgmp"

prg1: prg1.o helpers_gmp.o
	$(CC) -o $@ $^ $(CFLAGS)

prg2: prg2.o helpers_gmp.o
	$(CC) -o $@ $^ $(CFLAGS)

prg3: prg3.o helpers_gmp.o
	$(CC) -o $@ $^ $(CFLAGS)	

prg4: prg4.o helpers_gmp.o
	$(CC) -o $@ $^ $(CFLAGS)

prg5: prg5.o helpers_gmp.o
	$(CC) -o $@ $^ $(CFLAGS)

prg6: prg6.o helpers_gmp.o
	$(CC) -o $@ $^ $(CFLAGS)

prg7: prg7.o helpers_gmp.o
	$(CC) -o $@ $^ $(CFLAGS)	

prg8: prg8.o helpers_gmp.o
	$(CC) -o $@ $^ $(CFLAGS)

prg9: prg9.o helpers_gmp.o
	$(CC) -o $@ $^ $(CFLAGS)

prg10: prg10.o helpers_gmp.o
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c
	$(CC) -c -o $@ $^ $(CFLAGS)
