CFLAGS=-Wall -std=gnu99 -g -O3
LDFLAGS=-lpthread -lm $(CFLAGS)
TARGETS=find_roots_seq find_roots_para

all: $(TARGETS)

find_roots_para.o : find_roots_para.c

find_roots_seq.o : find_roots_seq.c

find_roots_para : find_roots_para.o

find_roots_seq : find_roots_seq.o

.PHONY : clean

clean :
	-$(RM) *~ *.o
	-find . -name '*_flymake*' -exec rm {} \;
	-$(RM) $(TARGETS)
