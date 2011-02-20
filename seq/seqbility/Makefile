CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2 -m64
CXXFLAGS=	$(CFLAGS)
DFLAGS=
PROG=		gen_mask apply_mask_l apply_mask_s splitfa
INCLUDES=	
LIBS=		-lm -lz -lpthread -Lbwt_gen -lbwtgen

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

gen_mask:gen_mask.o
		$(CC) $(CFLAGS) gen_mask.o -o $@ -lz

apply_mask_l:apply_mask_l.o
		$(CC) $(CFLAGS) apply_mask_l.o -o $@ -lz

apply_mask_s:apply_mask_s.o
		$(CC) $(CFLAGS) apply_mask_s.o -o $@ -lz

splitfa:splitfa.o
		$(CC) $(CFLAGS) splitfa.o -o $@ -lz

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
