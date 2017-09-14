METISLIB = ./metis/lib/libmetis.a
INCLUDES = -I./metis/
LIBDIR= -L/usr/lib -L./metis/lib -lmpi -lm
CC=/usr/lib64/lam/bin/mpicc
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(DEBUGFLAGS) $(INCLUDES)

OBJS = ./metis/io.o graph_init.o graph_transfer.o

         
.c.o:
	$(CC) $(CFLAGS) -c $*.c
graph_obj_example: $(OBJS) graph_main.o
	$(CC) -o graph_obj_example $(OBJS) $(METISLIB) graph_main.o $(LIBDIR)

clean:
	rm -f *~ *.o graph_obj_example
