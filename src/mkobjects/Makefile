OPTFLAGS=-Ofast
LDFLAGS=-lm
CFLAGS= $(OPTFLAGS)

all: mkdisk mksphere

mkdisk: mkdisk.o

mksphere: mksphere.o

clean:
	+$(RM) *.o

mostlyclean: clean
	$(RM) *~

INDENTFLAGS=\
	-bad\
	-bap\
	-c0\
	-cd0\
	-cd2\
	-d0\
	-di16\
	-i4\
	-l132\
	-nbbb\
	-nbfda\
	-ncdb\
	-nfca\
	-npsl\
	-nsc\
	-nut\
	-orig\
	-ppi4\
	-sob\
	-ts4\

INDENT=/usr/bin/indent $(INDENTFLAGS)

indent: mksphere.c mkdisk.c
	$(INDENT) $(C-codes) $(header-files)

