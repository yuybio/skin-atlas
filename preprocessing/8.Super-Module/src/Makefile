SOURCE = $(wildcard *.cpp)
TARGETS = $(patsubst %.cpp, %, $(SOURCE)) 

CXX = g++
CPPFLAGS = -O3

all:$(TARGETS)

$(TARGETS):%:%.cpp
	$(CXX) $< $(CPPFLAGS) -o $@  
	
.PHONY:clean all  


