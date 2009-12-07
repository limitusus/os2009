EXES=int_seq int_para
OBJS=int_seq.o int_para.o

all : $(EXES)
	echo $(EXES)

$(OBJS) : %.o : %.c
	$(CC) -c $< -Wall

int_seq : int_seq.o
	$(CC) -o $@ $^ -Wall -lm

int_para : int_para.o
	$(CC) -o $@ $^ -lpthread -Wall -lm

clean :
	-$(RM) *~ $(OBJS) $(EXES)

.PHONY : all clean
