include ../Make.inc

PREFIX = ..
COMPAT_DIR = ../compat
INCLUDE_DIR = ..
CFLAGS += -I$(INCLUDE_DIR)
LIBS += -lrt

OBJS =	cmdline.o computeGraph.o globals.o graph-manip.o graph_shared.o graphio.o read-graph-el.o timer.o xmalloc.o xmt-luc-blech.o getUserParameters.o undirected.o genScalData.o

all:	GraphCT_srv

GraphCT_srv:	server.o $(addprefix $(PREFIX)/,$(OBJS)) $(addprefix $(PREFIX)/,$(COMPATLIB))
	$(LD) $(LDFLAGS) -o GraphCT_srv server.o $(OBJS) $(LIBS) $(addprefix $(PREFIX)/,$(COMPATLIB))

server.o:  server.c
	$(CC) $(CFLAGS) -c $<

$(PREFIX)/xmt-luc-blech.o:    $(PREFIX)/xmt-luc-blech.cc
	$(CXX) $(CFLAGS) -c $<

$(PREFIX)/%.o:	$(PREFIX)/%.c
	$(CC) $(CFLAGS) -c $<

$(PREFIX)/$(IOLIB) $(PREFIX)/$(COMPATLIB):
	$(MAKE) -C $(PREFIX)/$(IOLIB_DIR)

clean:
	rm -f *.o GraphCT.pl GraphCT_srv
