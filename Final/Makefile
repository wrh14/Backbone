XX = g++
CFLAGS = -g 
CLIBS = -lpthread

INCLUDE = $(wildcard ../SparseMatrix/*.h)
SOURCES = $(wildcard ./*.cpp ../SparseMatrix/*.cpp)
INCLUDE_DIRS = -I../SparseMatrix/*.h

TARGET = nodeSelect

OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

$(TARGET) : $(OBJECTS)
	$(XX)  -std=c++11 $(CFLAGS) $^ -o $@ $(CLIBS)
$(OBJECTS) : %.o : %.cpp 
	$(XX)  -std=c++11 -c $(CFLAGS) $< -o $@ $(INCLUDE_DIRS)
