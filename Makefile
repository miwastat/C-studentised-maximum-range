#!sh

OBJ=smrng_lq.o smrng_lp.o rng_lp.o nrml_p.o
CC=gcc

# Strip *.exe files in Windows_NT
ifeq ($(OS),Windows_NT)
	EXE=.exe
endif

smrng_tbl: smrng_tbl.o $(OBJ)
	$(CC) smrng_tbl.o $(OBJ) -o smrng_tbl -lm
	strip smrng_tbl$(EXE)

smrng_tbl.o: smrng_tbl.c
	$(CC) -c smrng_tbl.c

smrng_lq_tst: smrng_lq_tst.o $(OBJ)
	$(CC) smrng_lq_tst.o $(OBJ) -o smrng_lq_tst -lm
	strip smrng_lq_tst$(EXE)

smrng_lq_tst.o: smrng_lq_tst.c
	$(CC) -c smrng_lq_tst.c

smrng_lq.o: smrng_lq.c
	$(CC) -c smrng_lq.c

smrng_lp_tst: smrng_lp_tst.o smrng_lp.o rng_lp.o nrml_p.o
	$(CC) smrng_lp_tst.o smrng_lp.o rng_lp.o nrml_p.o -o smrng_lp_tst -lm
	strip smrng_lp_tst$(EXE)

smrng_lp_tst.o: smrng_lp_tst.c
	$(CC) -c smrng_lp_tst.c

smrng_lp.o: smrng_lp.c
	$(CC) -c smrng_lp.c

rng_lp_tst: rng_lp_tst.o rng_lp.o nrml_p.o
	$(CC) rng_lp_tst.o rng_lp.o nrml_p.o -o rng_lp_tst -lm
	strip rng_lp_tst$(EXE)

rng_lp_tst.o: rng_lp_tst.c
	$(CC) -c rng_lp_tst.c

rng_lp.o: rng_lp.c 
	$(CC) -c rng_lp.c

nrml_p.o: nrml_p.c
	$(CC) -c nrml_p.c

