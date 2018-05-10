include flags
SUBDIRS=applications/ baseTypes/
OBJECTFILES=baseTypes/simpleExecutable.o
LIBFILES=symfpu.a
PROGS=test


.PHONY: all subdirs $(SUBDIRS) clean $(PROGS)

all : subdirs $(LIBFILES) $(PROGS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

symfpu.a : $(OBJECTFILES)
	ar rcs $@ $^

clean :
	find . -name '*.o' -exec rm {} \;
	rm -f $(LIBFILES) $(PROGS)

test : applications/test.o $(LIBFILES)
	$(CXX) $(CXXFLAGS) $^ -o $@

cbmcverification : applications/cbmcverification.o $(LIBFILES)
	$(CXX) $(CXXFLAGS) $^ -o $@

generate : applications/generate.o $(LIBFILES)
	$(CXX) $(CXXFLAGS) $^ -o $@

