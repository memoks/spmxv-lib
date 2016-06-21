all:	kway.o
	gcc -g3 $(CFLAGS) -o kway kway.o -lpatoh64 -L. -I. -lm

.c .o:
	gcc -g3 $(CFLAGS) -o $@ -c $<
