include ../flags
CXXFLAGS+=-I../../
ALL=test.o

.PHONY : all

all : $(ALL)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

