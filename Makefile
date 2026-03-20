CXXFLAGS  += -I./include -I./include/star -DHAVE_CONFIG_H -Wconversion
LDLIBS	  += -Lsstar -lsstar -Lnode -lnode -Lstd -lstd -lm
CFLAGS    += -O
CXXFLAGS  += -O

EXE	= starev
#DIRS	= node std sstar dstar
DIRS	= node std sstar dstar rdc

all:	lib $(EXE)

test:
	echo $(CXX)
	echo `which $(CXX)`

Makefile.inc:	Makefile.inc.conf
	@sed s#__BASE__#$(PWD)#g Makefile.inc.conf >Makefile.inc

starev.o: starev.C
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DTOOLBOX -c -o $@ $^

starev: starev.o
	$(CXX) -o $@ $(LDFLAGS) $< $(LDLIBS)

lib:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done

clean:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done
	/bin/rm -f *~ $(EXE) Makefile.inc
