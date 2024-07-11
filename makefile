
RM = rm -f
CC = gcc -Wall $(CNORMALFLAG) -Wno-unused-result
#CFLAGS = -g
CFLAGS = -O3
LIBS = -lm # -pthread
SRCS = proc4opt.c proc3opt.c proc2opt.c stopwatch.c heap.c woeginger.c io4opt.c
OBJS = $(SRCS:%.c=%.o)
HDRS = heap.h io4opt.h proc3opt.h proc4opt.h stopwatch.h tspmain.h
#use this to compute/read costs at beginning and store them in a matrix
DEFINED =

#use this to compute costs on the fly instead of storing them in a matrix. E.g., use for large n
#DEFINED =  -DCOSTS_ONFLY 
EXECUTABLES = KoptLS KoptLSHP # 4opt1move da sistemare...

all : $(EXECUTABLES)

KoptLS : clean ${HDRS} $(SRCS) cmd_convergenceKopt.c
	${CC} ${CFLAGS} ${DEFINED} cmd_convergenceKopt.c ${SRCS} -o $@ $(LIBS)

KoptLSHP : clean ${HDRS} $(SRCS) cmd_convergenceKopt.c
	${CC} ${CFLAGS} ${DEFINED} cmd_convergenceKopt.c ${SRCS} -o $@ $(LIBS)

KoptLSW10 : clean ${HDRS} $(SRCS) cmd_convergenceKopt.c
	${CC} ${CFLAGS} ${DEFINED} cmd_convergenceKopt.c ${SRCS} -o $@ $(LIBS)

KoptLScyg : clean ${HDRS} $(SRCS) cmd_convergenceKopt.c
	${CC} ${CFLAGS} ${DEFINED} cmd_convergenceKopt.c ${SRCS} -o $@ $(LIBS)

clean :
	${RM} ${OBJS}

distclean : clean
	${RM} $(EXECUTABLES)









