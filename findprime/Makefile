EXES=prime_seq prime_para prime_exclusive
OBJS=prime_seq.o prime_para.o prime_exclusive.o

all : $(EXES)
	echo $(EXES)

$(OBJS) : %.o : %.c
	$(CC) -c $< -Wall

prime_seq : prime_seq.o
	$(CC) -o $@ $^ -Wall

prime_para : prime_para.o
	$(CC) -o $@ $^ -lpthread -Wall

prime_exclusive : prime_exclusive.o
	$(CC) -o $@ $^ -lpthread -Wall

clean :
	-$(RM) *~ $(OBJS) $(EXES)

.PHONY : all clean
