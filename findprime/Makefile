EXES=prime_seq prime_para prime_exclusive
OBJS=prime_seq.o prime_para.o prime_exclusive.o
CFLAGS= -Wall -Werror -DNDEBUG
LDFLAGS= -Wall -Werror

all : $(EXES)
	echo $(EXES)

$(OBJS) : %.o : %.c
	$(CC) -c $< $(CFLAGS)

prime_seq : prime_seq.o
	$(CC) -o $@ $^ $(LDFLAGS)

prime_para : prime_para.o
	$(CC) -o $@ $^ -lpthread $(LDFLAGS)

prime_exclusive : prime_exclusive.o
	$(CC) -o $@ $^ -lpthread $(LDFLAGS)

clean :
	-$(RM) *~ $(OBJS) $(EXES)

.PHONY : all clean
